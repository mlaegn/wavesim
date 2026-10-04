import { listRuns, loadRunFromFiles, loadRunFromUrl } from "./load";
import { drawProfile } from "./profile";
import type { Run } from "./run";
import { RunView, type Overlay, type Preset } from "./view";

function el<T extends HTMLElement>(id: string): T {
  const node = document.getElementById(id);
  if (!node) throw new Error(`missing #${id}`);
  return node as T;
}

const canvas = el<HTMLCanvasElement>("view");
const overlay = el("overlay");
const message = el("message");
const progress = el<HTMLProgressElement>("progress");
const controls = el("controls");
const info = el("info");
const facts = el("facts");
const credit = el("credit");
const meta = el("meta");
const runsSelect = el<HTMLSelectElement>("runs");
const playButton = el<HTMLButtonElement>("play");
const scrub = el<HTMLInputElement>("scrub");
const clock = el("clock");
const speedSelect = el<HTMLSelectElement>("speed");
const exag = el<HTMLInputElement>("exag");
const exagOut = el("exag-out");
const foam = el<HTMLInputElement>("foam");
const zones = el<HTMLInputElement>("zones");
const mapKind = el<HTMLSelectElement>("mapkind");
const crop = el<HTMLInputElement>("crop");
const side = el<HTMLInputElement>("side");
const profilePanel = el("profile-panel");
const profileCanvas = el<HTMLCanvasElement>("profile");
const transect = el<HTMLInputElement>("transect");
const transectAt = el("transect-at");
const profileReadout = el("profile-readout");

let view: RunView | null = null;
let time = 0;
let playing = true;

function say(text: string, kind: "info" | "error" = "info"): void {
  overlay.hidden = false;
  overlay.classList.toggle("error", kind === "error");
  message.textContent = text;
}

function hideOverlay(): void {
  overlay.hidden = true;
  progress.hidden = true;
}

const fmt = (v: unknown, unit = "", digits = 1): string | null =>
  typeof v === "number" ? `${Number.isInteger(v) ? v : v.toFixed(digits)}${unit}` : null;

function showFacts(run: Run): void {
  const h = run.header;
  const w = h.waves;
  const rows: [string, string | null][] = [
    ["Wave height", fmt(w.wave_height_m, " m")],
    ["Period", fmt(w.wave_period_s, " s")],
    ["Tide", fmt(w.tide_m, " m")],
    ["Wave-maker depth", fmt(w.maker_depth_m, " m")],
    ["Friction (Manning)", typeof w.manning === "number" && w.manning > 0 ? String(w.manning) : null],
    ["Grid", `${h.nx} × ${h.ny} · ${h.dx} m`],
    ["Frames", `${h.frame_count} (every ${fmt(w.frame_interval_s, " s") ?? "?"})`],
    ["Heading", `${Math.round(h.frame.x_bearing_deg)}° (waves travel)`],
  ];
  facts.replaceChildren();
  for (const [name, value] of rows) {
    if (value === null) continue;
    const dt = document.createElement("dt");
    dt.textContent = name;
    const dd = document.createElement("dd");
    dd.textContent = value;
    facts.append(dt, dd);
  }
  const s = h.source;
  const parts: string[] = [];
  if (s && typeof s.name === "string" && s.name) parts.push(`Bathymetry: ${s.name}`);
  if (s && typeof s.attribution === "string" && s.attribution) parts.push(`Credit: ${s.attribution}`);
  if (s && typeof s.datum === "string" && s.datum) parts.push(`Datum: ${s.datum}`);
  parts.push(`Model: ${typeof w.model === "string" ? w.model : "shallow water"}`);
  credit.textContent = parts.join(" · ");
  info.hidden = false;
}

/** What happened along one line over the whole run: where it got steepest and where it broke. */
function wholeRun(run: Run, row: number): string {
  const { nx, dx } = run.header;
  const box = run.viewBox;
  const first = Math.ceil(box.x0 / dx);
  const last = Math.min(nx - 1, Math.floor(box.x1 / dx));
  const steep = run.steepnessMap();
  let best = { s: 0, i: first };
  for (let i = first; i <= last; i++) {
    const s = steep[row * nx + i] as number;
    if (s > best.s) best = { s, i };
  }
  const out: string[] = [];
  if (best.s > 0) {
    const deg = (Math.atan(best.s) * 180) / Math.PI;
    out.push(`Over the whole run it gets steepest at x = ${((best.i + 0.5) * dx).toFixed(0)} m (slope ${best.s.toFixed(2)}, ${deg.toFixed(0)}°).`);
  }
  const breaks = run.breakMap();
  if (breaks) {
    let from = -1;
    let to = -1;
    let top = { f: 0, i: first };
    for (let i = first; i <= last; i++) {
      const f = breaks[row * nx + i] as number;
      if (f >= 0.05) {
        if (from < 0) from = i;
        to = i;
      }
      if (f > top.f) top = { f, i };
    }
    out.push(
      from < 0
        ? "It never breaks along this line."
        : `It breaks from x = ${((from + 0.5) * dx).toFixed(0)} to ${((to + 0.5) * dx).toFixed(0)} m, most often at x = ${((top.i + 0.5) * dx).toFixed(0)} m (${(100 * top.f).toFixed(0)}% of the time).`,
    );
  }
  return out.join(" ");
}

/** Redraw the side view for the current time and the chosen line along the shore. */
function updateProfile(): void {
  if (!view || profilePanel.hidden) return;
  const run = view.run;
  const row = Number(transect.value);
  const y = (row + 0.5) * run.header.dy;
  transectAt.textContent = `y = ${y.toFixed(0)} m`;
  view.setTransect(y);
  const box = run.viewBox;
  const readout = drawProfile(profileCanvas, run, view.eta, view.breaking, row, box.x0, Math.min(run.shoreline + 40, box.x1));
  const parts: string[] = [];
  if (readout.steepest) {
    const { slope, x, depth } = readout.steepest;
    const degrees = (Math.atan(slope) * 180) / Math.PI;
    parts.push(`Steepest face: slope ${slope.toFixed(2)} (${degrees.toFixed(0)}°) at x = ${x.toFixed(0)} m, ${depth.toFixed(1)} m deep.`);
  }
  if (view.breaking) {
    parts.push(
      readout.breakingFrom
        ? `The model reads breaking from x = ${readout.breakingFrom.x.toFixed(0)} m (${readout.breakingFrom.depth.toFixed(1)} m deep).`
        : "Not breaking along this line right now.",
    );
  }
  parts.push(wholeRun(run, row));
  profileReadout.textContent = parts.filter(Boolean).join(" ");
}

function setRun(run: Run, label: string): void {
  view?.dispose();
  view = new RunView(canvas, run);
  time = run.start;
  scrub.min = String(run.start);
  scrub.max = String(run.duration);
  scrub.value = String(run.start);
  view.setExaggeration(Number(exag.value));
  view.setPreset("oblique");
  view.setFoam(foam.checked);
  view.setZones(zones.checked);
  view.setCrop(crop.checked);
  mapKind.value = "none";
  view.setOverlay("none");
  (mapKind.options[1] as HTMLOptionElement).disabled = !run.hasField("breaking");
  const box = run.viewBox;
  transect.min = String(Math.ceil(box.y0 / run.header.dy));
  transect.max = String(Math.floor(box.y1 / run.header.dy) - 1);
  transect.value = String(Math.round((Number(transect.min) + Number(transect.max)) / 2));
  profilePanel.hidden = !side.checked;
  showFacts(run);
  meta.textContent = `${label} · ${run.header.nx} × ${run.header.ny} cells`;
  controls.hidden = false;
  hideOverlay();
  playing = true;
  playButton.textContent = "Pause";
}

function onProgress(got: number, total: number, what: string): void {
  progress.hidden = false;
  progress.value = got / total;
  message.textContent = `Loading ${what}… ${(got / 1e6).toFixed(0)} of ${(total / 1e6).toFixed(0)} MB`;
}

async function load(task: () => Promise<Run>, label: string): Promise<void> {
  try {
    say("Loading run…");
    setRun(await task(), label);
  } catch (e) {
    progress.hidden = true;
    say(e instanceof Error ? e.message : String(e), "error");
  }
}

function openByName(name: string): Promise<void> {
  const base = /^(https?:)?\/\//.test(name) || name.startsWith("/") ? name : `/runs/${encodeURIComponent(name)}`;
  const url = new URL(location.href);
  url.searchParams.set("run", name);
  history.replaceState(null, "", url);
  return load(() => loadRunFromUrl(base, onProgress), name);
}

// Playback ---------------------------------------------------------------------------

function seek(t: number): void {
  if (!view) return;
  time = Math.min(Math.max(t, view.run.start), view.run.duration);
  view.setTime(time);
  updateProfile();
  scrub.value = String(time);
  clock.textContent = `${time.toFixed(1)} / ${view.run.duration.toFixed(0)} s`;
}

function setPlaying(on: boolean): void {
  playing = on;
  playButton.textContent = on ? "Pause" : "Play";
}

let last = performance.now();
function frame(now: number): void {
  const dt = Math.min((now - last) / 1000, 0.1);
  last = now;
  if (view) {
    if (playing) {
      const { start, duration: end } = view.run;
      const next = time + dt * Number(speedSelect.value);
      seek(next >= end ? start + (next - end) : next);
    }
    view.render();
  }
  requestAnimationFrame(frame);
}

// Controls ---------------------------------------------------------------------------

playButton.addEventListener("click", () => setPlaying(!playing));
scrub.addEventListener("input", () => {
  setPlaying(false);
  seek(Number(scrub.value));
});
exag.addEventListener("input", () => {
  exagOut.textContent = exag.value;
  view?.setExaggeration(Number(exag.value));
});
foam.addEventListener("change", () => view?.setFoam(foam.checked));
zones.addEventListener("change", () => view?.setZones(zones.checked));
mapKind.addEventListener("change", () => view?.setOverlay(mapKind.value as Overlay));
crop.addEventListener("change", () => view?.setCrop(crop.checked));
side.addEventListener("change", () => {
  profilePanel.hidden = !side.checked || !view;
  updateProfile();
});
transect.addEventListener("input", updateProfile);
for (const b of document.querySelectorAll<HTMLButtonElement>("button[data-view]")) {
  b.addEventListener("click", () => view?.setPreset(b.dataset.view as Preset));
}
runsSelect.addEventListener("change", () => void openByName(runsSelect.value));

window.addEventListener("keydown", (e) => {
  if (!view || (e.target instanceof HTMLElement && /^(INPUT|SELECT)$/.test(e.target.tagName) && e.key !== " ")) return;
  const step = view.run.header.times.length > 1 ? (view.run.header.times[1] as number) - (view.run.header.times[0] as number) : 1;
  if (e.key === " ") {
    e.preventDefault();
    setPlaying(!playing);
  } else if (e.key === "ArrowRight") {
    setPlaying(false);
    seek(time + step);
  } else if (e.key === "ArrowLeft") {
    setPlaying(false);
    seek(time - step);
  } else if (e.key === "Home") {
    seek(view.run.start);
  }
});

window.addEventListener("resize", () => {
  view?.resize();
  updateProfile();
});

// Files ------------------------------------------------------------------------------

el<HTMLInputElement>("files").addEventListener("change", (e) => {
  const files = Array.from((e.target as HTMLInputElement).files ?? []);
  if (files.length) void load(() => loadRunFromFiles(files), "local files");
});
window.addEventListener("dragover", (e) => {
  e.preventDefault();
  overlay.hidden = false;
  overlay.classList.add("drag");
});
window.addEventListener("dragleave", (e) => {
  if (e.relatedTarget === null) {
    overlay.classList.remove("drag");
    if (view) hideOverlay();
  }
});
window.addEventListener("drop", (e) => {
  e.preventDefault();
  overlay.classList.remove("drag");
  const files = Array.from(e.dataTransfer?.files ?? []);
  if (files.length) void load(() => loadRunFromFiles(files), "local files");
  else if (view) hideOverlay();
});

// Start ------------------------------------------------------------------------------

declare global {
  interface Window {
    /** For tests and screenshots. */
    wavesim: { seek: (t: number) => void; pause: () => void; view: () => RunView | null };
  }
}
window.wavesim = { seek, pause: () => setPlaying(false), view: () => view };

async function start(): Promise<void> {
  const names = await listRuns();
  runsSelect.replaceChildren(...names.map((n) => new Option(n, n)));
  runsSelect.hidden = names.length === 0;
  const wanted = new URL(location.href).searchParams.get("run");
  const first = wanted ?? names[0];
  if (first) {
    if (names.includes(first)) runsSelect.value = first;
    await openByName(first);
  } else {
    say("No runs found. Run `wavesim run` first, or drop run.json, bed.f32 and frames.f32 here.");
  }
}

requestAnimationFrame(frame);
void start();
