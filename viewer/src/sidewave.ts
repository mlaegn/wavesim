/**
 * The side view of one breaking wave: the bed, the water and its surface coloured by how fast
 * the water moves against the wave, played back in slow motion.
 */
import { listSlices, loadSlice, type Slice, speedColour } from "./slice";

function el<T extends HTMLElement>(id: string): T {
  const node = document.getElementById(id);
  if (!node) throw new Error(`missing #${id}`);
  return node as T;
}

const canvas = el<HTMLCanvasElement>("side");
const choose = el<HTMLSelectElement>("slices");
const message = el("message");
const controls = el("controls");
const playButton = el<HTMLButtonElement>("play");
const scrub = el<HTMLInputElement>("scrub");
const clock = el("clock");
const speedSelect = el<HTMLSelectElement>("speed");
const viewSelect = el<HTMLSelectElement>("frame");
const breakButton = el<HTMLButtonElement>("to-break");
const facts = el("facts");
const now = el("now");

let slice: Slice | null = null;
let time = 0;
let playing = false;
let last = 0;

/** World box shown: x and z ranges, and the vertical stretch. */
interface Box {
  x0: number;
  x1: number;
  z0: number;
  z1: number;
  stretch: number;
}

function wave(s: Slice): { depth: number; height: number } {
  const w = s.header.wave;
  const depth = typeof w.depth === "number" ? w.depth : -point(s, 0)[1];
  const height = typeof w.height === "number" ? w.height : 0.3 * depth;
  return { depth, height };
}

/** Point `k` of the bed, which the header check guarantees has at least two. */
function point(s: Slice, k: number): [number, number] {
  return s.header.bed[k] as [number, number];
}

/** The bed's height at `x`, straight between its points. */
function bedAt(s: Slice, x: number): number {
  const bed = s.header.bed;
  let k = 1;
  while (k < bed.length - 1 && point(s, k)[0] < x) k++;
  const [a, b] = [point(s, k - 1), point(s, k)];
  const w = Math.min(1, Math.max(0, (x - a[0]) / (b[0] - a[0])));
  return a[1] + (b[1] - a[1]) * w;
}

/** Where the wave breaks: where its face stood vertical, or failing that the last crest. */
function breakX(s: Slice): number {
  if (s.header.curl) return s.header.curl.x;
  const last = s.header.frame_count - 1;
  const x = s.field(last, 0);
  const z = s.field(last, 1);
  let top = 0;
  for (let i = 1; i < z.length; i++) if ((z[i] as number) > (z[top] as number)) top = i;
  return x[top] as number;
}

function box(s: Slice, kind: string, width: number, height: number): Box {
  const bed = s.header.bed;
  const { depth, height: h } = wave(s);
  if (kind === "whole") {
    const x0 = point(s, 0)[0];
    const x1 = point(s, bed.length - 1)[0];
    const z0 = -depth * 1.05;
    const z1 = h * 1.8;
    const stretch = Math.max(1, ((x1 - x0) / (z1 - z0)) * (height / width));
    return { x0, x1, z0, z1, stretch };
  }
  // A still camera at true proportions on the whole break, from a little behind where the face
  // stands vertical to a little past where the lip lands, the wave running into it.
  const from = breakX(s);
  const to = s.header.landing?.x ?? from;
  const local = -bedAt(s, from);
  const z0 = -Math.max(local * 1.6, 0.5 * h);
  const z1 = 1.9 * h;
  const needed = to - from + 4.5 * h;
  const span = Math.max(((z1 - z0) * width) / height, needed);
  const x0 = from - 3 * h - (span - needed) / 2;
  return { x0, x1: x0 + span, z0, z1, stretch: 1 };
}

function draw(): void {
  const ctx = canvas.getContext("2d");
  if (!ctx) return;
  const dpr = window.devicePixelRatio || 1;
  const width = canvas.clientWidth;
  const height = canvas.clientHeight;
  if (canvas.width !== Math.round(width * dpr) || canvas.height !== Math.round(height * dpr)) {
    canvas.width = Math.round(width * dpr);
    canvas.height = Math.round(height * dpr);
  }
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  ctx.clearRect(0, 0, width, height);
  if (!slice) return;
  const s = slice;
  const b = box(s, viewSelect.value, width, height);
  // Uniform scale in x; z stretched by b.stretch.
  const scale = Math.min(width / (b.x1 - b.x0), height / ((b.z1 - b.z0) * b.stretch));
  const left = (width - (b.x1 - b.x0) * scale) / 2;
  const top = (height - (b.z1 - b.z0) * b.stretch * scale) / 2;
  const px = (x: number): number => left + (x - b.x0) * scale;
  const pz = (z: number): number => top + (b.z1 - z) * b.stretch * scale;

  const surface = s.at(time);
  const bed = s.header.bed;
  const bottom = Math.max(height, pz(b.z0));

  // Water: the surface, down the right wall, back along the bed, up the left wall.
  const water = ctx.createLinearGradient(0, pz(0), 0, pz(-wave(s).depth));
  water.addColorStop(0, "#1f7fae");
  water.addColorStop(1, "#0b3550");
  ctx.beginPath();
  const sx = (i: number): number => px(surface.x[i] as number);
  const sz = (i: number): number => pz(surface.z[i] as number);
  for (let i = 0; i < surface.x.length; i++) {
    if (i === 0) ctx.moveTo(sx(i), sz(i));
    else ctx.lineTo(sx(i), sz(i));
  }
  for (let k = bed.length - 1; k >= 0; k--) ctx.lineTo(px(point(s, k)[0]), pz(point(s, k)[1]));
  ctx.closePath();
  ctx.fillStyle = water;
  ctx.fill();

  // The bed, filled down to the bottom of the view.
  ctx.beginPath();
  const [first, end] = [point(s, 0), point(s, bed.length - 1)];
  ctx.moveTo(px(first[0]), bottom);
  for (const [x, z] of bed) ctx.lineTo(px(x), pz(z));
  ctx.lineTo(px(end[0]), bottom);
  ctx.closePath();
  ctx.fillStyle = "#c9a96e";
  ctx.fill();
  ctx.strokeStyle = "#8a6a35";
  ctx.lineWidth = 1.5;
  ctx.stroke();

  // Past the shallow wall that ends the computed water, the beach carried on at its last slope,
  // faded, so the wall does not read as a cliff.
  if (px(end[0]) < width) {
    const before = point(s, bed.length - 2);
    const slope = (end[1] - before[1]) / (end[0] - before[0]);
    const x = b.x1;
    const z = Math.min(end[1] + slope * (x - end[0]), wave(s).height);
    ctx.beginPath();
    ctx.moveTo(px(end[0]), bottom);
    ctx.lineTo(px(end[0]), pz(end[1]));
    ctx.lineTo(px(x), pz(z));
    ctx.lineTo(px(x), bottom);
    ctx.closePath();
    ctx.fillStyle = "rgba(201, 169, 110, 0.45)";
    ctx.fill();
    ctx.strokeStyle = "rgba(255, 255, 255, 0.55)";
    ctx.setLineDash([3, 4]);
    ctx.beginPath();
    ctx.moveTo(px(end[0]), pz(end[1]));
    ctx.lineTo(px(end[0]), pz(wave(s).height));
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "rgba(255, 255, 255, 0.75)";
    ctx.font = "11px system-ui, sans-serif";
    ctx.fillText("computed water ends", px(end[0]) + 5, pz(wave(s).height));
  }

  // Still-water level.
  ctx.strokeStyle = "rgba(255, 255, 255, 0.25)";
  ctx.setLineDash([4, 6]);
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(px(first[0]), pz(0));
  ctx.lineTo(px(end[0]), pz(0));
  ctx.stroke();
  ctx.setLineDash([]);

  // The surface, each piece coloured by the water's speed against the wave's.
  const c = s.waveSpeed;
  let fastest = 0;
  ctx.lineWidth = 3;
  ctx.lineCap = "round";
  for (let i = 0; i + 1 < surface.x.length; i++) {
    const ratio = ((surface.speed[i] as number) + (surface.speed[i + 1] as number)) / (2 * c);
    fastest = Math.max(fastest, ratio);
    ctx.strokeStyle = speedColour(ratio);
    ctx.beginPath();
    ctx.moveTo(sx(i), sz(i));
    ctx.lineTo(sx(i + 1), sz(i + 1));
    ctx.stroke();
  }

  // Where the face stood vertical and where the lip lands, once they have happened.
  const mark = (x: number, z: number, label: string): void => {
    ctx.fillStyle = "rgba(255, 255, 255, 0.9)";
    ctx.beginPath();
    ctx.arc(px(x), pz(z), 4, 0, 2 * Math.PI);
    ctx.fill();
    ctx.font = "12px system-ui, sans-serif";
    ctx.fillText(label, px(x) + 7, pz(z) - 7);
  };
  const { curl, landing } = s.header;
  if (curl && time >= curl.time) mark(curl.x, 0, `face stood vertical above here, ${curl.time.toFixed(2)} s`);
  if (landing && time >= landing.time) mark(landing.x, landing.z, `lip lands, ${landing.time.toFixed(2)} s`);

  now.textContent =
    `${time.toFixed(2)} s · fastest water ${(fastest * c).toFixed(1)} m/s, ` +
    `${fastest.toFixed(2)} of the wave's speed` +
    (b.stretch > 1.05 ? ` · height stretched ×${b.stretch.toFixed(1)}` : "");
  clock.textContent = `${time.toFixed(2)} / ${s.duration.toFixed(2)} s`;
  scrub.value = String(time);
}

function showFacts(s: Slice): void {
  const { depth, height } = wave(s);
  const w = s.header.wave;
  const rows: [string, string][] = [
    ["Seabed", String(w.case ?? s.header.name)],
    ["Wave", `${height.toFixed(1)} m on ${depth.toFixed(1)} m of water`],
  ];
  const { curl, landing } = s.header;
  if (curl) {
    rows.push(["Face vertical", `${curl.time.toFixed(2)} s, in ${curl.depth.toFixed(2)} m of water, crest ${curl.crest.toFixed(2)} m up`]);
  } else {
    rows.push(["Face vertical", "never: it did not plunge"]);
  }
  if (landing) {
    rows.push(["Lip lands", `${landing.time.toFixed(2)} s, ${landing.throw.toFixed(2)} m further on`]);
    rows.push([
      "Tube",
      `${landing.tube_height.toFixed(2)} m tall, ${landing.tube_width.toFixed(2)} m wide, ${landing.tube_area.toFixed(2)} m² of air`,
    ]);
  }
  if (s.header.ended !== "landed") rows.push(["Ended", s.header.ended]);
  facts.replaceChildren(
    ...rows.flatMap(([k, v]) => {
      const dt = document.createElement("dt");
      dt.textContent = k;
      const dd = document.createElement("dd");
      dd.textContent = v;
      return [dt, dd];
    }),
  );
}

/** Draw once, soon; repeated calls before then draw once. */
let pending = false;
function redraw(): void {
  if (pending) return;
  pending = true;
  requestAnimationFrame(() => {
    pending = false;
    draw();
  });
}

/** Advance and draw while playing; when paused nothing runs until something changes. */
function tick(stamp: number): void {
  if (!slice || !playing) {
    last = 0;
    return;
  }
  const elapsed = last ? (stamp - last) / 1000 : 0;
  last = stamp;
  time += elapsed * Number(speedSelect.value);
  if (time >= slice.duration) {
    time = slice.duration;
    playing = false;
    playButton.textContent = "Play";
  }
  draw();
  if (playing) requestAnimationFrame(tick);
  else last = 0;
}

function play(on: boolean): void {
  // Start the loop only if it is not already running, or it would run twice as fast.
  const start = on && !playing;
  playing = on;
  playButton.textContent = on ? "Pause" : "Play";
  if (start) requestAnimationFrame(tick);
  redraw();
}

async function open(name: string): Promise<void> {
  message.textContent = `Loading ${name}…`;
  message.hidden = false;
  try {
    slice = await loadSlice(`/runs/${encodeURIComponent(name)}`);
  } catch (e) {
    message.textContent = String(e instanceof Error ? e.message : e);
    return;
  }
  message.hidden = true;
  controls.hidden = false;
  scrub.max = String(slice.duration);
  time = 0;
  showFacts(slice);
  play(true);
  const url = new URL(window.location.href);
  url.searchParams.set("slice", name);
  history.replaceState(null, "", url);
}

playButton.addEventListener("click", () => {
  if (!slice) return;
  if (time >= slice.duration) time = 0;
  play(!playing);
});
scrub.addEventListener("input", () => {
  time = Number(scrub.value);
  play(false);
});
breakButton.addEventListener("click", () => {
  if (!slice) return;
  const at = slice.header.curl?.time ?? slice.duration;
  time = Math.max(0, at - 1.5);
  speedSelect.value = "0.1";
  viewSelect.value = "break";
  play(true);
});
choose.addEventListener("change", () => void open(choose.value));
viewSelect.addEventListener("change", redraw);
window.addEventListener("resize", redraw);

void (async () => {
  const names = await listSlices();
  if (!names.length) {
    message.textContent = "No slices yet: run `wavesim slice reef` (or gentle, steep) and reload.";
    return;
  }
  choose.replaceChildren(
    ...names.map((n) => {
      const option = document.createElement("option");
      option.value = n;
      option.textContent = n.replace(/^slice-/, "");
      return option;
    }),
  );
  const wanted = new URL(window.location.href).searchParams.get("slice");
  choose.value = wanted && names.includes(wanted) ? wanted : (names[0] as string);
  await open(choose.value);
})();
