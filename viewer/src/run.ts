/**
 * A wavesim run in memory: the header, the bed and every frame of free-surface
 * elevation. No DOM and no WebGL here, so it can be tested in Node. The byte layout is
 * specified in docs/data-format.md.
 */

export class RunError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "RunError";
  }
}

export interface RunHeader {
  format: string;
  version: number;
  name: string;
  nx: number;
  ny: number;
  dx: number;
  dy: number;
  bed: string;
  frames: string;
  fields: string[];
  dtype: string;
  frame_count: number;
  times: number[];
  frame: {
    crs: string;
    origin_easting: number;
    origin_northing: number;
    x_bearing_deg: number;
  };
  waves: Record<string, unknown>;
  source?: Record<string, unknown> | null;
}

const isRecord = (v: unknown): v is Record<string, unknown> =>
  typeof v === "object" && v !== null && !Array.isArray(v);

function num(o: Record<string, unknown>, key: string): number {
  const v = o[key];
  if (typeof v !== "number" || !Number.isFinite(v)) {
    throw new RunError(`run.json: "${key}" must be a number`);
  }
  return v;
}

function posInt(o: Record<string, unknown>, key: string): number {
  const v = num(o, key);
  if (!Number.isInteger(v) || v <= 0) throw new RunError(`run.json: "${key}" must be a positive integer`);
  return v;
}

function str(o: Record<string, unknown>, key: string): string {
  const v = o[key];
  if (typeof v !== "string") throw new RunError(`run.json: "${key}" must be a string`);
  return v;
}

/** Check a parsed `run.json` and return it typed. Throws a `RunError` saying what is wrong. */
export function parseHeader(raw: unknown): RunHeader {
  if (!isRecord(raw)) throw new RunError("run.json is not an object");
  if (raw.format !== "wavesim-run") {
    throw new RunError(`not a wavesim run: format is ${JSON.stringify(raw.format)}`);
  }
  if (raw.version !== 1) throw new RunError(`unsupported run version ${String(raw.version)} (this viewer reads 1)`);
  if (str(raw, "dtype") !== "f32-le") throw new RunError(`unsupported dtype ${str(raw, "dtype")}`);

  const fields = raw.fields;
  if (!Array.isArray(fields) || !fields.every((f) => typeof f === "string") || !fields.includes("eta")) {
    throw new RunError('run.json: "fields" must be a list that includes "eta"');
  }
  const times = raw.times;
  if (!Array.isArray(times) || !times.every((t) => typeof t === "number" && Number.isFinite(t))) {
    throw new RunError('run.json: "times" must be a list of numbers');
  }
  const frameCount = posInt(raw, "frame_count");
  if (times.length !== frameCount) {
    throw new RunError(`run.json: ${times.length} times but frame_count is ${frameCount}`);
  }
  for (let k = 1; k < times.length; k++) {
    if (!((times[k] as number) > (times[k - 1] as number))) throw new RunError("run.json: times must increase");
  }
  const frame = raw.frame;
  if (!isRecord(frame)) throw new RunError('run.json: "frame" is missing');

  return {
    format: "wavesim-run",
    version: 1,
    name: str(raw, "name"),
    nx: posInt(raw, "nx"),
    ny: posInt(raw, "ny"),
    dx: num(raw, "dx"),
    dy: num(raw, "dy"),
    bed: str(raw, "bed"),
    frames: str(raw, "frames"),
    fields: fields as string[],
    dtype: "f32-le",
    frame_count: frameCount,
    times: times as number[],
    frame: {
      crs: str(frame, "crs"),
      origin_easting: num(frame, "origin_easting"),
      origin_northing: num(frame, "origin_northing"),
      x_bearing_deg: num(frame, "x_bearing_deg"),
    },
    waves: isRecord(raw.waves) ? raw.waves : {},
    source: isRecord(raw.source) ? raw.source : null,
  };
}

/** True when this machine stores floats little-endian, as the files do. */
export function isLittleEndian(): boolean {
  return new Uint8Array(new Uint16Array([1]).buffer)[0] === 1;
}

function floats(buf: ArrayBuffer, expectedCount: number, what: string): Float32Array {
  if (buf.byteLength !== expectedCount * 4) {
    throw new RunError(`${what}: ${buf.byteLength} bytes, expected ${expectedCount * 4} (${expectedCount} values)`);
  }
  return new Float32Array(buf);
}

export class Run {
  readonly header: RunHeader;
  /** Bed elevation in metres above mean sea level, `nx * ny`, y outer, x inner. */
  readonly bed: Float32Array;
  private readonly frames: Float32Array;
  private readonly etaField: number;

  constructor(header: RunHeader, bed: Float32Array, frames: Float32Array) {
    this.header = header;
    this.bed = bed;
    this.frames = frames;
    this.etaField = header.fields.indexOf("eta");
  }

  get cells(): number {
    return this.header.nx * this.header.ny;
  }

  get frameCount(): number {
    return this.header.frame_count;
  }

  /** Simulated seconds covered by the run. */
  get duration(): number {
    return this.header.times[this.frameCount - 1] as number;
  }

  /** Free-surface elevation of stored frame `k` (a view, not a copy). */
  frame(k: number): Float32Array {
    if (!Number.isInteger(k) || k < 0 || k >= this.frameCount) throw new RangeError(`no frame ${k}`);
    const start = (k * this.header.fields.length + this.etaField) * this.cells;
    return this.frames.subarray(start, start + this.cells);
  }

  /**
   * Free-surface elevation at time `t` seconds, written into `out`. Between stored frames
   * it interpolates with a Catmull-Rom cubic, which follows a wave far better than a
   * straight line does: linearly blending frames 2 s apart would shrink a 14 s wave by
   * 10% halfway between them. The result never goes below the bed.
   */
  sampleInto(t: number, out: Float32Array): void {
    if (out.length !== this.cells) throw new RangeError("output has the wrong size");
    const times = this.header.times;
    const last = this.frameCount - 1;
    if (this.frameCount === 1 || t <= (times[0] as number)) return void out.set(this.frame(0));
    if (t >= (times[last] as number)) return void out.set(this.frame(last));

    // Last frame whose time is <= t.
    let lo = 0;
    let hi = last;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if ((times[mid] as number) <= t) lo = mid;
      else hi = mid;
    }
    const u = (t - (times[lo] as number)) / ((times[lo + 1] as number) - (times[lo] as number));
    const p0 = this.frame(Math.max(lo - 1, 0));
    const p1 = this.frame(lo);
    const p2 = this.frame(lo + 1);
    const p3 = this.frame(Math.min(lo + 2, last));
    const u2 = u * u;
    const u3 = u2 * u;
    const c0 = -0.5 * u3 + u2 - 0.5 * u;
    const c1 = 1.5 * u3 - 2.5 * u2 + 1;
    const c2 = -1.5 * u3 + 2 * u2 + 0.5 * u;
    const c3 = 0.5 * u3 - 0.5 * u2;
    const bed = this.bed;
    for (let n = 0; n < out.length; n++) {
      const v = c0 * (p0[n] as number) + c1 * (p1[n] as number) + c2 * (p2[n] as number) + c3 * (p3[n] as number);
      const b = bed[n] as number;
      out[n] = v > b ? v : b;
    }
  }
}

/** Assemble a run from the three files' contents, checking sizes against the header. */
export function runFromBuffers(headerJson: unknown, bedBytes: ArrayBuffer, frameBytes: ArrayBuffer): Run {
  if (!isLittleEndian()) throw new RunError("this device is big-endian; the run files are little-endian");
  const header = parseHeader(headerJson);
  const cells = header.nx * header.ny;
  const bed = floats(bedBytes, cells, header.bed);
  const frames = floats(frameBytes, header.frame_count * header.fields.length * cells, header.frames);
  for (const v of bed) if (!Number.isFinite(v)) throw new RunError(`${header.bed} contains a value that is not finite`);
  return new Run(header, bed, frames);
}
