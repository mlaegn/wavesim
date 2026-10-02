import { listRuns, loadRunFromFiles, loadRunFromUrl } from "./load";
import type { Run } from "./run";
import { RunView, type Preset } from "./view";

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

function setRun(run: Run, label: string): void {
  view?.dispose();
  view = new RunView(canvas, run);
  time = 0;
  scrub.max = String(run.duration);
  scrub.value = "0";
  view.setExaggeration(Number(exag.value));
  view.setFoam(foam.checked);
  view.setZones(zones.checked);
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
  time = Math.min(Math.max(t, 0), view.run.duration);
  view.setTime(time);
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
      const end = view.run.duration;
      const next = time + dt * Number(speedSelect.value);
      seek(next >= end ? next - end : next);
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
    seek(0);
  }
});

window.addEventListener("resize", () => view?.resize());

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
