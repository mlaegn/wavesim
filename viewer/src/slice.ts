/**
 * A wave breaking in a vertical slice, as `wavesim slice` writes it (see docs/data-format.md):
 * `slice.json` and `frames.f32`, the surface as a chain of nodes in every frame over a fixed
 * bed between two walls.
 */

export class SliceError extends Error {}

export interface Curl {
  time: number;
  x: number;
  crest: number;
  depth: number;
}

export interface Landing {
  time: number;
  x: number;
  z: number;
  throw: number;
  tube_area: number;
  tube_width: number;
  tube_height: number;
}

export interface SliceHeader {
  format: "wavesim-slice";
  version: 1;
  name: string;
  bed: [number, number][];
  nodes: number;
  frames: string;
  fields: string[];
  dtype: string;
  frame_count: number;
  times: number[];
  wave: Record<string, unknown>;
  ended: string;
  curl?: Curl;
  landing?: Landing;
}

const isNumber = (v: unknown): v is number => typeof v === "number" && Number.isFinite(v);

function numbers<K extends string>(raw: unknown, keys: readonly K[], what: string): Record<K, number> {
  if (typeof raw !== "object" || raw === null) throw new SliceError(`${what} must be an object`);
  const out = {} as Record<K, number>;
  for (const k of keys) {
    const v = (raw as Record<string, unknown>)[k];
    if (!isNumber(v)) throw new SliceError(`${what}.${k} must be a number`);
    out[k] = v;
  }
  return out;
}

/** Check a parsed `slice.json`, with messages that say what is wrong. */
export function parseSliceHeader(raw: unknown): SliceHeader {
  if (typeof raw !== "object" || raw === null) throw new SliceError("slice.json must be an object");
  const h = raw as Record<string, unknown>;
  if (h.format !== "wavesim-slice" || h.version !== 1) throw new SliceError("not a wavesim-slice version 1");
  const count = (k: string): number => {
    const v = h[k];
    if (!Number.isInteger(v) || (v as number) < 0) throw new SliceError(`${k} must be a whole number`);
    return v as number;
  };
  const nodes = count("nodes");
  const frameCount = count("frame_count");
  if (nodes < 2) throw new SliceError("a slice needs at least two surface nodes");
  const bed = h.bed;
  if (!Array.isArray(bed) || bed.length < 2 || !bed.every((p) => Array.isArray(p) && p.length === 2 && p.every(isNumber))) {
    throw new SliceError("bed must be a list of at least two [x, z] points");
  }
  const times = h.times;
  if (!Array.isArray(times) || times.length !== frameCount || !times.every(isNumber)) {
    throw new SliceError("times must hold one number per frame");
  }
  for (let k = 1; k < times.length; k++) {
    if ((times[k] as number) < (times[k - 1] as number)) throw new SliceError("times must not go backwards");
  }
  const fields = h.fields;
  if (!Array.isArray(fields) || fields.join() !== "x,z,speed") throw new SliceError('fields must be ["x", "z", "speed"]');
  if (h.dtype !== "<f4") throw new SliceError("dtype must be <f4");
  if (typeof h.frames !== "string") throw new SliceError("frames must name a file");
  const curl = h.curl === undefined ? undefined : numbers(h.curl, ["time", "x", "crest", "depth"] as const, "curl");
  const landing =
    h.landing === undefined
      ? undefined
      : numbers(h.landing, ["time", "x", "z", "throw", "tube_area", "tube_width", "tube_height"] as const, "landing");
  return {
    format: "wavesim-slice",
    version: 1,
    name: typeof h.name === "string" ? h.name : "slice",
    bed: bed as [number, number][],
    nodes,
    frames: h.frames,
    fields: fields as string[],
    dtype: "<f4",
    frame_count: frameCount,
    times: times as number[],
    wave: typeof h.wave === "object" && h.wave !== null ? (h.wave as Record<string, unknown>) : {},
    ended: typeof h.ended === "string" ? h.ended : "unknown",
    curl,
    landing,
  };
}

/** The surface at one moment: node positions and the water's speed at each. */
export interface Surface {
  x: Float32Array;
  z: Float32Array;
  speed: Float32Array;
}

export class Slice {
  constructor(
    readonly header: SliceHeader,
    private readonly data: Float32Array,
  ) {}

  get duration(): number {
    const t = this.header.times;
    return t[t.length - 1] ?? 0;
  }

  /** One field of one frame: 0 for x, 1 for z, 2 for speed. */
  field(frame: number, field: number): Float32Array {
    const n = this.header.nodes;
    const start = (frame * 3 + field) * n;
    return this.data.subarray(start, start + n);
  }

  /** The surface at time `t`, straight between the frames either side. */
  at(t: number): Surface {
    const times = this.header.times;
    const time = (k: number): number => times[k] as number;
    const last = times.length - 1;
    let k = 0;
    while (k < last && time(k + 1) <= t) k++;
    const k1 = Math.min(k + 1, last);
    const span = time(k1) - time(k);
    const w = span > 0 ? Math.min(1, Math.max(0, (t - time(k)) / span)) : 0;
    const mix = (field: number): Float32Array => {
      const a = this.field(k, field);
      const b = this.field(k1, field);
      const out = new Float32Array(a.length);
      for (let i = 0; i < a.length; i++) out[i] = (a[i] as number) + ((b[i] as number) - (a[i] as number)) * w;
      return out;
    };
    return { x: mix(0), z: mix(1), speed: mix(2) };
  }

  /** The long-wave speed offshore, `sqrt(g (h + H))`: the speed the wave arrives with. */
  get waveSpeed(): number {
    const w = this.header.wave;
    const g = isNumber(w.g) ? w.g : 9.81;
    const depth = isNumber(w.depth) ? w.depth : -(this.header.bed[0]?.[1] ?? 0);
    const height = isNumber(w.height) ? w.height : 0;
    return Math.sqrt(g * (depth + height));
  }
}

/** A slice from `slice.json` (parsed) and the bytes of its frames file. */
export function sliceFromBuffers(raw: unknown, frames: ArrayBuffer): Slice {
  const header = parseSliceHeader(raw);
  const expected = header.frame_count * header.nodes * 3 * 4;
  if (frames.byteLength !== expected) {
    throw new SliceError(`${header.frames} has ${frames.byteLength} bytes; slice.json promises ${expected}`);
  }
  return new Slice(header, new Float32Array(frames));
}

/** Load a slice from a directory URL such as `/runs/slice-reef`. */
export async function loadSlice(base: string): Promise<Slice> {
  const root = base.replace(/\/+$/, "");
  const response = await fetch(`${root}/slice.json`);
  if (!response.ok) throw new SliceError(`${root}/slice.json: HTTP ${response.status}`);
  const raw: unknown = await response.json();
  const header = parseSliceHeader(raw);
  const frames = await fetch(`${root}/${header.frames}`);
  if (!frames.ok) throw new SliceError(`${root}/${header.frames}: HTTP ${frames.status}`);
  return sliceFromBuffers(raw, await frames.arrayBuffer());
}

/** Names of the slices the dev server can see. */
export async function listSlices(): Promise<string[]> {
  try {
    const response = await fetch("/runs/index.json");
    if (!response.ok) return [];
    const body: unknown = await response.json();
    if (typeof body === "object" && body !== null && "slices" in body && Array.isArray(body.slices)) {
      return body.slices.filter((n): n is string => typeof n === "string");
    }
  } catch {
    // No listing available.
  }
  return [];
}

/**
 * The colour for water moving at `ratio` times the wave's speed: pale blue when it is slow,
 * through white to orange as it nears the wave's speed, red where it outruns the wave, which
 * is where a lip is thrown.
 */
export function speedColour(ratio: number): string {
  const stops: [number, [number, number, number]][] = [
    [0, [160, 210, 235]],
    [0.6, [255, 255, 255]],
    [0.9, [255, 170, 60]],
    [1.2, [215, 25, 20]],
  ];
  const r = Math.max(0, Math.min(1.2, ratio));
  const stop = (k: number): [number, [number, number, number]] => stops[k] as [number, [number, number, number]];
  let k = 0;
  while (k < stops.length - 2 && r > stop(k + 1)[0]) k++;
  const [a, ca] = stop(k);
  const [b, cb] = stop(k + 1);
  const w = (r - a) / (b - a);
  const c = ca.map((v, i) => Math.round(v + ((cb[i] as number) - v) * w));
  return `rgb(${c[0]}, ${c[1]}, ${c[2]})`;
}
