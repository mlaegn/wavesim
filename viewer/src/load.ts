import { parseHeader, Run, RunError, runFromBuffers } from "./run";

export type Progress = (loadedBytes: number, totalBytes: number, what: string) => void;

async function fetchBytes(url: string, what: string, onProgress?: Progress): Promise<ArrayBuffer> {
  const response = await fetch(url);
  if (!response.ok) throw new RunError(`${url}: HTTP ${response.status}`);
  const total = Number(response.headers.get("content-length")) || 0;
  const encoded = response.headers.get("content-encoding");
  if (!response.body || !total || encoded) return await response.arrayBuffer();

  // Stream into one preallocated buffer so a 120 MB run is held in memory once.
  const bytes = new Uint8Array(total);
  const reader = response.body.getReader();
  let got = 0;
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    if (got + value.length > total) throw new RunError(`${url}: more data than the server announced`);
    bytes.set(value, got);
    got += value.length;
    onProgress?.(got, total, what);
  }
  if (got !== total) throw new RunError(`${url}: connection ended after ${got} of ${total} bytes`);
  return bytes.buffer;
}

/** Load a run from a directory URL such as `/runs/pipeline`. */
export async function loadRunFromUrl(base: string, onProgress?: Progress): Promise<Run> {
  const root = base.replace(/\/+$/, "");
  const response = await fetch(`${root}/run.json`);
  if (!response.ok) throw new RunError(`${root}/run.json: HTTP ${response.status}`);
  const raw: unknown = await response.json();
  const header = parseHeader(raw);
  const bed = await fetchBytes(`${root}/${header.bed}`, header.bed, onProgress);
  const frames = await fetchBytes(`${root}/${header.frames}`, header.frames, onProgress);
  const stats = header.stats ? await fetchBytes(`${root}/${header.stats.file}`, header.stats.file, onProgress) : null;
  return runFromBuffers(raw, bed, frames, stats);
}

/** Load a run from files the user picked or dropped: run.json plus the files it names. */
export async function loadRunFromFiles(files: File[]): Promise<Run> {
  const byName = new Map(files.map((f) => [f.name, f]));
  const headerFile = byName.get("run.json");
  if (!headerFile) throw new RunError("include run.json (and the bed and frames files it names)");
  const raw: unknown = JSON.parse(await headerFile.text());
  const header = parseHeader(raw);
  const need = (name: string): File => {
    const f = byName.get(name);
    if (!f) throw new RunError(`run.json names ${name}, which was not included`);
    return f;
  };
  const bed = await need(header.bed).arrayBuffer();
  const frames = await need(header.frames).arrayBuffer();
  // The statistics are optional: without them the maps fall back to the frames.
  const statsFile = header.stats ? byName.get(header.stats.file) : undefined;
  const stats = statsFile ? await statsFile.arrayBuffer() : null;
  return runFromBuffers(raw, bed, frames, stats);
}

/** Names of the runs the dev server can see, or an empty list if there is no listing. */
export async function listRuns(): Promise<string[]> {
  try {
    const response = await fetch("/runs/index.json");
    if (!response.ok) return [];
    const body: unknown = await response.json();
    if (typeof body === "object" && body !== null && "runs" in body && Array.isArray(body.runs)) {
      return body.runs.filter((n): n is string => typeof n === "string");
    }
  } catch {
    // No listing available (a static build, for instance).
  }
  return [];
}
