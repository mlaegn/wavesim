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
  /** Fields recorded at every step over a window of the run, in their own file. */
  stats?: { file: string; fields: string[]; from_s: number; to_s: number } | null;
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
    stats: parseStats(raw.stats),
  };
}

function parseStats(raw: unknown): RunHeader["stats"] {
  if (raw === undefined || raw === null) return null;
  if (!isRecord(raw)) throw new RunError('run.json: "stats" must be an object');
  const fields = raw.fields;
  if (!Array.isArray(fields) || !fields.every((f) => typeof f === "string")) {
    throw new RunError('run.json: "stats.fields" must be a list of names');
  }
  return { file: str(raw, "file"), fields: fields as string[], from_s: num(raw, "from_s"), to_s: num(raw, "to_s") };
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

/** How a wave breaks, from Battjes' surf-similarity number: mellow, heavy, or riding up. */
export type Breaker = "spilling" | "plunging" | "surging";

/**
 * Where a wave breaks on one row: metres along x, still-water depth and its height there, and,
 * for runs that record them, the slope it meets and how it breaks.
 */
export interface BreakPoint {
  x: number;
  depth: number;
  height: number;
  slope?: number;
  surfSimilarity?: number;
  breaker?: Breaker;
}

function breakPoint(p: unknown): BreakPoint | null {
  if (Array.isArray(p) && p.length === 3 && p.every((v) => typeof v === "number")) {
    return { x: p[0] as number, depth: p[1] as number, height: p[2] as number };
  }
  if (!isRecord(p) || typeof p.x_m !== "number" || typeof p.depth_m !== "number" || typeof p.height_m !== "number") {
    return null;
  }
  const breaker = p.breaker === "spilling" || p.breaker === "plunging" || p.breaker === "surging" ? p.breaker : undefined;
  return {
    x: p.x_m,
    depth: p.depth_m,
    height: p.height_m,
    ...(typeof p.slope === "number" ? { slope: p.slope } : {}),
    ...(typeof p.surf_similarity === "number" ? { surfSimilarity: p.surf_similarity } : {}),
    ...(breaker ? { breaker } : {}),
  };
}

export class Run {
  readonly header: RunHeader;
  /** Bed elevation in metres above mean sea level, `nx * ny`, y outer, x inner. */
  readonly bed: Float32Array;
  private readonly frames: Float32Array;
  private readonly stats: Float32Array | null;
  private readonly etaField: number;
  private steepness: Float32Array | undefined;

  constructor(header: RunHeader, bed: Float32Array, frames: Float32Array, stats: Float32Array | null = null) {
    this.header = header;
    this.bed = bed;
    this.frames = frames;
    this.stats = stats;
    this.etaField = header.fields.indexOf("eta");
  }

  get cells(): number {
    return this.header.nx * this.header.ny;
  }

  get frameCount(): number {
    return this.header.frame_count;
  }

  /** Simulated time of the first frame: 0, or later for the fine grid of a nested run. */
  get start(): number {
    return this.header.times[0] as number;
  }

  /** Simulated time of the last frame. */
  get duration(): number {
    return this.header.times[this.frameCount - 1] as number;
  }

  hasField(name: string): boolean {
    return this.header.fields.includes(name);
  }

  private fieldIndex(name: string): number {
    const f = this.header.fields.indexOf(name);
    if (f < 0) throw new RunError(`this run has no "${name}" field`);
    return f;
  }

  private block(f: number, k: number): Float32Array {
    if (!Number.isInteger(k) || k < 0 || k >= this.frameCount) throw new RangeError(`no frame ${k}`);
    const start = (k * this.header.fields.length + f) * this.cells;
    return this.frames.subarray(start, start + this.cells);
  }

  /** Free-surface elevation of stored frame `k` (a view, not a copy). */
  frame(k: number): Float32Array {
    return this.block(this.etaField, k);
  }

  /** Field `name` of stored frame `k` (a view, not a copy). */
  field(name: string, k: number): Float32Array {
    return this.block(this.fieldIndex(name), k);
  }

  /**
   * Free-surface elevation at time `t` seconds, written into `out`. Between stored frames
   * it interpolates with a Catmull-Rom cubic, which follows a wave far better than a
   * straight line does: linearly blending frames 2 s apart would shrink a 14 s wave by
   * 10% halfway between them. The result never goes below the bed.
   */
  sampleInto(t: number, out: Float32Array): void {
    this.sample("eta", t, out);
  }

  /** Any field at time `t`, interpolated the same way. `breaking` stays within 0 to 1. */
  sampleFieldInto(name: string, t: number, out: Float32Array): void {
    this.sample(name, t, out);
  }

  private sample(name: string, t: number, out: Float32Array): void {
    if (out.length !== this.cells) throw new RangeError("output has the wrong size");
    const f = this.fieldIndex(name);
    const times = this.header.times;
    const last = this.frameCount - 1;
    if (this.frameCount === 1 || t <= (times[0] as number)) return void out.set(this.block(f, 0));
    if (t >= (times[last] as number)) return void out.set(this.block(f, last));

    // Last frame whose time is <= t.
    let lo = 0;
    let hi = last;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if ((times[mid] as number) <= t) lo = mid;
      else hi = mid;
    }
    const u = (t - (times[lo] as number)) / ((times[lo + 1] as number) - (times[lo] as number));
    const p0 = this.block(f, Math.max(lo - 1, 0));
    const p1 = this.block(f, lo);
    const p2 = this.block(f, lo + 1);
    const p3 = this.block(f, Math.min(lo + 2, last));
    const u2 = u * u;
    const u3 = u2 * u;
    const c0 = -0.5 * u3 + u2 - 0.5 * u;
    const c1 = 1.5 * u3 - 2.5 * u2 + 1;
    const c2 = -1.5 * u3 + 2 * u2 + 0.5 * u;
    const c3 = 0.5 * u3 - 0.5 * u2;
    const bed = this.bed;
    const isEta = f === this.etaField;
    for (let n = 0; n < out.length; n++) {
      const v = c0 * (p0[n] as number) + c1 * (p1[n] as number) + c2 * (p2[n] as number) + c3 * (p3[n] as number);
      if (isEta) {
        const b = bed[n] as number;
        out[n] = v > b ? v : b;
      } else {
        out[n] = v < 0 ? 0 : v > 1 ? 1 : v;
      }
    }
  }

  /** Statistics field `name`, recorded at every step, or `null` if the run has none. */
  stat(name: string): Float32Array | null {
    const f = this.header.stats?.fields.indexOf(name) ?? -1;
    if (f < 0 || !this.stats) return null;
    return this.stats.subarray(f * this.cells, (f + 1) * this.cells);
  }

  /** The height of each wave, crest to trough, averaged over the waves the run recorded at
   * every step, or `null` for a run that did not record it. */
  waveHeight(): Float32Array | null {
    return this.stat("wave_height");
  }

  /**
   * Where the waves break, row by row along the shore: the point where the wave is tallest
   * (which is where Ting & Kirby's laboratory saw it break), with the still-water depth and
   * the wave's height there; `null` for a row with no break, and an empty list for a run that
   * did not record it.
   */
  get breakLine(): (BreakPoint | null)[] {
    const raw = this.header.waves.break_line;
    return Array.isArray(raw) ? raw.map(breakPoint) : [];
  }

  /**
   * For every cell, the steepest the surface got there during the run: the largest
   * `|grad eta|` over all frames, from central differences. Cells are skipped where they or
   * a neighbour have less than 30 cm of water, because a shoreline is not a wave face.
   * Slope is dimensionless (rise over run); 0.1 is about 6 degrees and 0.4 about 22.
   */
  steepnessMap(): Float32Array {
    if (this.steepness) return this.steepness;
    const { nx, ny, dx, dy } = this.header;
    const map = new Float32Array(this.cells);
    const minDepth = 0.3;
    for (let k = 0; k < this.frameCount; k++) {
      const eta = this.frame(k);
      for (let j = 1; j < ny - 1; j++) {
        for (let i = 1; i < nx - 1; i++) {
          const n = j * nx + i;
          if (
            (eta[n] as number) - (this.bed[n] as number) < minDepth ||
            (eta[n - 1] as number) - (this.bed[n - 1] as number) < minDepth ||
            (eta[n + 1] as number) - (this.bed[n + 1] as number) < minDepth ||
            (eta[n - nx] as number) - (this.bed[n - nx] as number) < minDepth ||
            (eta[n + nx] as number) - (this.bed[n + nx] as number) < minDepth
          ) {
            continue;
          }
          const sx = ((eta[n + 1] as number) - (eta[n - 1] as number)) / (2 * dx);
          const sy = ((eta[n + nx] as number) - (eta[n - nx] as number)) / (2 * dy);
          const s = Math.hypot(sx, sy);
          if (s > (map[n] as number)) map[n] = s;
        }
      }
    }
    return (this.steepness = map);
  }

  /**
   * Where the water meets the land, in metres along x: the median over rows of the first
   * place the bed rises above sea level. The length of the domain if there is no shore.
   */
  get shoreline(): number {
    const { nx, ny, dx } = this.header;
    const crossings: number[] = [];
    for (let j = 0; j < ny; j++) {
      for (let i = 1; i < nx; i++) {
        if ((this.bed[j * nx + i] as number) >= 0 && (this.bed[j * nx + i - 1] as number) < 0) {
          crossings.push((i + 0.5) * dx);
          break;
        }
      }
    }
    if (crossings.length === 0) return nx * dx;
    crossings.sort((a, b) => a - b);
    return crossings[crossings.length >> 1] as number;
  }

  /**
   * The part of the domain that is sea, in metres: past the wave-maker's own bump on the
   * seaward side, and inside the absorbing strips along the two sides. Everything outside
   * it is there for the numerics, not the ocean.
   */
  get viewBox(): { x0: number; x1: number; y0: number; y1: number } {
    const { nx, ny, dx } = this.header;
    const side = this.sideMargin;
    return { x0: this.viewStart, x1: nx * dx, y0: side, y1: ny * this.header.dy - side };
  }

  /**
   * Metres along each side that are there for the numerics: the margin by the walls, which
   * mirror the bed (`side_margin_cells`), or in older runs absorbing strips (`sponge_side_cells`).
   */
  get sideMargin(): number {
    const w = this.header.waves;
    const cells = w.side_margin_cells ?? w.sponge_side_cells;
    return typeof cells === "number" ? cells * this.header.dy : 0;
  }

  /** Where the viewer should start: past the wave-maker's own bump, which is not sea. */
  get viewStart(): number {
    const w = this.header.waves;
    const v = w.near_field_end_m ?? w.maker_x_m;
    return typeof v === "number" ? v : 0;
  }
}

/** Assemble a run from its files' contents, checking sizes against the header. */
export function runFromBuffers(
  headerJson: unknown,
  bedBytes: ArrayBuffer,
  frameBytes: ArrayBuffer,
  statsBytes: ArrayBuffer | null = null,
): Run {
  if (!isLittleEndian()) throw new RunError("this device is big-endian; the run files are little-endian");
  const header = parseHeader(headerJson);
  const cells = header.nx * header.ny;
  const bed = floats(bedBytes, cells, header.bed);
  const frames = floats(frameBytes, header.frame_count * header.fields.length * cells, header.frames);
  for (const v of bed) if (!Number.isFinite(v)) throw new RunError(`${header.bed} contains a value that is not finite`);
  const st = header.stats;
  const stats = st && statsBytes ? floats(statsBytes, st.fields.length * cells, st.file) : null;
  return new Run(header, bed, frames, stats);
}
