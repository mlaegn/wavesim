import fs from "node:fs";
import path from "node:path";
import { describe, expect, it } from "vitest";
import { breakingStrength, parseHeader, Run, RunError, runFromBuffers, type RunHeader } from "../src/run";

const header = (over: Record<string, unknown> = {}) => ({
  format: "wavesim-run",
  version: 1,
  name: "t",
  nx: 3,
  ny: 2,
  dx: 3,
  dy: 3,
  bed: "bed.f32",
  frames: "frames.f32",
  fields: ["eta"],
  dtype: "f32-le",
  layout: "",
  frame_count: 4,
  times: [0, 2, 4, 6],
  frame: { crs: "local", origin_easting: 0, origin_northing: 0, x_bearing_deg: 0 },
  waves: {},
  ...over,
});

const buf = (values: number[]) => new Float32Array(values).buffer as ArrayBuffer;

/** Build a run whose frame k, cell n holds fn(k, n). */
function makeRun(fn: (k: number, n: number) => number, bed = -10, over: Record<string, unknown> = {}) {
  const h = header(over);
  const cells = 6;
  const frames: number[] = [];
  for (let k = 0; k < (h.frame_count as number); k++) for (let n = 0; n < cells; n++) frames.push(fn(k, n));
  return runFromBuffers(h, buf(new Array(cells).fill(bed)), buf(frames));
}

describe("parseHeader", () => {
  it("accepts a well-formed header", () => {
    const h: RunHeader = parseHeader(header());
    expect(h.nx).toBe(3);
    expect(h.times).toEqual([0, 2, 4, 6]);
  });

  it.each([
    [{ format: "something-else" }, /not a wavesim run/],
    [{ version: 2 }, /unsupported run version/],
    [{ dtype: "f64-le" }, /dtype/],
    [{ fields: ["h"] }, /eta/],
    [{ nx: 0 }, /nx/],
    [{ nx: 2.5 }, /nx/],
    [{ times: [0, 2, 4] }, /3 times but frame_count is 4/],
    [{ times: [0, 2, 2, 6] }, /increase/],
    [{ frame: null }, /frame/],
  ])("rejects %j", (over, message) => {
    expect(() => parseHeader(header(over))).toThrowError(message);
  });

  it("rejects things that are not objects", () => {
    expect(() => parseHeader(null)).toThrow(RunError);
    expect(() => parseHeader([])).toThrow(RunError);
  });
});

describe("runFromBuffers", () => {
  it("reads frames from the right offsets", () => {
    const run = makeRun((k, n) => 100 * k + n);
    expect(Array.from(run.frame(0))).toEqual([0, 1, 2, 3, 4, 5]);
    expect(Array.from(run.frame(3))).toEqual([300, 301, 302, 303, 304, 305]);
    expect(run.duration).toBe(6);
    expect(run.frameCount).toBe(4);
  });

  it("finds eta when it is not the first field", () => {
    // Two fields per frame, "eta" second: frame k = [ other (6), eta (6) ].
    const h = header({ fields: ["u", "eta"], frame_count: 2, times: [0, 1] });
    const frames: number[] = [];
    for (let k = 0; k < 2; k++) {
      frames.push(...new Array(6).fill(-1), ...[0, 1, 2, 3, 4, 5].map((n) => 10 * k + n));
    }
    const run = runFromBuffers(h, buf(new Array(6).fill(-10)), buf(frames));
    expect(Array.from(run.frame(1))).toEqual([10, 11, 12, 13, 14, 15]);
  });

  it("rejects files of the wrong size, naming the file", () => {
    const h = header();
    expect(() => runFromBuffers(h, buf([1, 2, 3]), buf(new Array(24).fill(0)))).toThrowError(/bed\.f32: 12 bytes, expected 24/);
    expect(() => runFromBuffers(h, buf(new Array(6).fill(0)), buf(new Array(23).fill(0)))).toThrowError(/frames\.f32/);
  });

  it("rejects a bed with gaps", () => {
    const bed = new Array(6).fill(0);
    bed[2] = NaN;
    expect(() => runFromBuffers(header(), buf(bed), buf(new Array(24).fill(0)))).toThrowError(/not finite/);
  });

  it("refuses frames outside the run", () => {
    const run = makeRun(() => 0);
    expect(() => run.frame(4)).toThrow(RangeError);
    expect(() => run.frame(-1)).toThrow(RangeError);
  });
});

describe("sampleInto", () => {
  const out = () => new Float32Array(6);

  it("returns the stored frames exactly at their times", () => {
    const run = makeRun((k, n) => Math.sin(k + n));
    for (let k = 0; k < 4; k++) {
      const o = out();
      run.sampleInto(2 * k, o);
      expect(Array.from(o)).toEqual(Array.from(run.frame(k)));
    }
  });

  it("holds the first and last frame outside the run", () => {
    const run = makeRun((k) => k);
    const o = out();
    run.sampleInto(-5, o);
    expect(o[0]).toBe(0);
    run.sampleInto(99, o);
    expect(o[0]).toBe(3);
  });

  it("reproduces anything that is linear in time", () => {
    const run = makeRun((k, n) => 0.5 * k + 0.1 * n);
    const o = out();
    run.sampleInto(3.0, o); // 1.5 frames in
    expect(o[0]).toBeCloseTo(0.75, 6);
    expect(o[5]).toBeCloseTo(0.75 + 0.5, 6);
  });

  it("follows a wave between frames much better than a straight line", () => {
    // A 14 s wave sampled every 2 s, as wavesim writes it, with a crest at t = 3: exactly
    // halfway between two frames, the worst place for a straight line.
    const period = 14;
    const f = (t: number) => Math.cos((2 * Math.PI * (t - 3)) / period);
    const h = header({ frame_count: 8, times: [0, 2, 4, 6, 8, 10, 12, 14] });
    const frames: number[] = [];
    for (let k = 0; k < 8; k++) frames.push(...new Array(6).fill(f(2 * k)));
    const run = runFromBuffers(h, buf(new Array(6).fill(-10)), buf(frames));

    let worstCubic = 0;
    let worstLinear = 0;
    const o = out();
    for (let t = 2; t <= 10; t += 0.25) {
      run.sampleInto(t, o);
      worstCubic = Math.max(worstCubic, Math.abs((o[0] as number) - f(t)));
      const k = Math.floor(t / 2);
      const u = t / 2 - k;
      const linear = (1 - u) * f(2 * k) + u * f(2 * k + 2);
      worstLinear = Math.max(worstLinear, Math.abs(linear - f(t)));
    }
    expect(worstLinear).toBeGreaterThan(0.09); // the crest is cut by 10%
    expect(worstCubic).toBeLessThan(0.02);
  });

  it("never puts the water below the bed", () => {
    // Frames 0,1 are dry (surface = bed = 0); frames 2,3 have water. The cubic would
    // undershoot between them without the clamp.
    const run = makeRun((k) => (k < 2 ? 0 : 1), 0);
    const o = out();
    for (let t = 0; t <= 6; t += 0.1) {
      run.sampleInto(t, o);
      expect(o[0]).toBeGreaterThanOrEqual(0);
    }
  });

  it("works with a single frame", () => {
    const h = header({ frame_count: 1, times: [0] });
    const run = runFromBuffers(h, buf(new Array(6).fill(-1)), buf([1, 2, 3, 4, 5, 6]));
    const o = out();
    run.sampleInto(12, o);
    expect(Array.from(o)).toEqual([1, 2, 3, 4, 5, 6]);
  });

  it("rejects an output of the wrong size", () => {
    const run = makeRun(() => 0);
    expect(() => run.sampleInto(0, new Float32Array(5))).toThrow(RangeError);
  });
});

describe("the file written by the Rust solver", () => {
  const dir = path.join(__dirname, "fixtures", "tiny");
  const exists = fs.existsSync(path.join(dir, "run.json"));

  it.skipIf(!exists)("loads and reads back the values the Rust test wrote", () => {
    const bytes = (f: string) => {
      const b = fs.readFileSync(path.join(dir, f));
      return b.buffer.slice(b.byteOffset, b.byteOffset + b.byteLength) as ArrayBuffer;
    };
    const run = runFromBuffers(JSON.parse(fs.readFileSync(path.join(dir, "run.json"), "utf8")), bytes("bed.f32"), bytes("frames.f32"));
    expect(run.header.nx).toBe(6);
    expect(run.header.ny).toBe(4);
    expect(run.header.frame.x_bearing_deg).toBe(130);
    expect(run.frameCount).toBe(4);
    expect(run.header.times).toEqual([0, 2, 4, 6]);
    expect(run.header.fields).toEqual(["eta", "breaking"]);
    expect(run.hasField("breaking")).toBe(true);
    // Cell (i = 3, j = 1): bed = -2 + 0.75 * 3 = 0.25, dry, so the surface equals the bed.
    const n = 1 * 6 + 3;
    expect(run.bed[n]).toBe(0.25);
    for (let k = 0; k < 4; k++) expect(run.frame(k)[n]).toBe(0.25);
    // Cell (i = 1, j = 0): bed = -1.25, wet: eta = 0.125 k + 0.0625 i - 0.25.
    for (let k = 0; k < 4; k++) expect(run.frame(k)[1]).toBeCloseTo(0.125 * k + 0.0625 - 0.25, 7);
    // The second field: wet cells hold 0.25 * ((n + k) % 4), dry cells 0.
    for (let k = 0; k < 4; k++) {
      expect(run.field("breaking", k)[1]).toBe(0.25 * ((1 + k) % 4));
      expect(run.field("breaking", k)[n]).toBe(0);
    }
  });
});

// Keep the type import honest.
export type { Run };

describe("fields and what they are used for", () => {
  /** Frames with two fields, "eta" then "breaking", as the solver writes them. */
  function twoFields(eta: (k: number, n: number) => number, breaking: (k: number, n: number) => number, bed = -10) {
    const h = header({ fields: ["eta", "breaking"] });
    const frames: number[] = [];
    for (let k = 0; k < 4; k++) {
      for (let n = 0; n < 6; n++) frames.push(eta(k, n));
      for (let n = 0; n < 6; n++) frames.push(breaking(k, n));
    }
    return runFromBuffers(h, buf(new Array(6).fill(bed)), buf(frames));
  }

  it("reads each field from its own block", () => {
    const run = twoFields((k, n) => 100 * k + n, (k, n) => 0.1 * k + 0.01 * n);
    expect(Array.from(run.frame(3))).toEqual([300, 301, 302, 303, 304, 305]);
    expect(run.field("breaking", 3)[5]).toBeCloseTo(0.35, 6);
    expect(() => run.field("velocity", 0)).toThrowError(/no "velocity" field/);
    expect(run.hasField("velocity")).toBe(false);
  });

  it("keeps an interpolated breaking field within 0 to 1", () => {
    // 0, 1, 0, 1 would overshoot a cubic between samples; the field is a fraction and must not.
    const run = twoFields(() => 0, (k) => (k % 2 === 0 ? 0 : 1));
    const o = new Float32Array(6);
    for (let t = 0; t <= 6; t += 0.05) {
      run.sampleFieldInto("breaking", t, o);
      expect(o[0]).toBeGreaterThanOrEqual(0);
      expect(o[0]).toBeLessThanOrEqual(1);
    }
    run.sampleFieldInto("breaking", 2, o);
    expect(o[0]).toBe(1);
  });

  it("finds the shoreline and the viewer's start", () => {
    // 3 x 2 cells of 3 m; bed rises through zero between the 2nd and 3rd column.
    const h = header({ waves: { maker_x_m: 3, near_field_end_m: 5.5 } });
    const bed = [-4, -1, 2, -4, -1, 2];
    const run = runFromBuffers(h, buf(bed), buf(new Array(24).fill(0)));
    expect(run.shoreline).toBe(7.5); // centre of the third column
    expect(run.viewStart).toBe(5.5);
    const noNearField = runFromBuffers(header({ waves: { maker_x_m: 3 } }), buf(bed), buf(new Array(24).fill(0)));
    expect(noNearField.viewStart).toBe(3);
    const allWater = runFromBuffers(header(), buf(new Array(6).fill(-5)), buf(new Array(24).fill(0)));
    expect(allWater.shoreline).toBe(9); // no shore: the length of the domain
  });
});

describe("nested and walled runs", () => {
  it("starts at its first frame, which a nested run's fine grid puts after zero", () => {
    const h = header({ frame_count: 2, times: [72, 74] });
    const run = runFromBuffers(h, buf(new Array(6).fill(-5)), buf(new Array(12).fill(0)));
    expect(run.start).toBe(72);
    expect(run.duration).toBe(74);
  });

  it("trims the margin by the walls, or the absorbing strips of older runs", () => {
    const bed = buf(new Array(6).fill(-5));
    const frames = () => buf(new Array(24).fill(0));
    const walled = runFromBuffers(header({ waves: { side_margin_cells: 1 } }), bed, frames());
    expect(walled.sideMargin).toBe(3);
    expect(walled.viewBox.y0).toBe(3);
    const older = runFromBuffers(header({ waves: { sponge_side_cells: 1 } }), bed, frames());
    expect(older.sideMargin).toBe(3);
    expect(runFromBuffers(header(), bed, frames()).sideMargin).toBe(0);
  });
});

describe("statistics recorded at every step", () => {
  it("give the height per wave and the break line", () => {
    const stats = { file: "stats.f32", fields: ["break_fraction", "wave_height"], from_s: 0, to_s: 6 };
    const line = [{ x_m: 4.5, depth_m: 2.0, height_m: 1.6, slope: 0.05, surf_similarity: 0.6, breaker: "plunging" }, null, [1, 2, 3]];
    const h = header({ fields: ["eta", "breaking"], stats, waves: { break_line: line } });
    const bed = buf(new Array(6).fill(-5));
    const frames = buf(new Array(48).fill(0)); // 4 frames x 2 fields x 6 cells
    const height = [0, 0.5, 1.0, 1.5, 2.0, 2.5];
    const run = runFromBuffers(h, bed, frames, buf([...new Array(6).fill(0), ...height]));
    expect(Array.from(run.waveHeight() ?? [])).toEqual(height);
    expect(run.breakLine).toEqual([
      { x: 4.5, depth: 2.0, height: 1.6, slope: 0.05, surfSimilarity: 0.6, breaker: "plunging" },
      null,
      { x: 1, depth: 2, height: 3 }, // the older form, without a type
    ]);
    // Without the statistics there is neither.
    const plain = runFromBuffers(header(), bed, buf(new Array(24).fill(0)));
    expect(plain.waveHeight()).toBeNull();
    expect(plain.breakLine).toEqual([]);
  });

  it("must have the size the header promises", () => {
    const stats = { file: "stats.f32", fields: ["break_fraction"], from_s: 0, to_s: 6 };
    const h = header({ stats });
    expect(() => runFromBuffers(h, buf(new Array(6).fill(-5)), buf(new Array(24).fill(0)), buf([1, 2]))).toThrow(RunError);
  });
});

describe("how strongly the bed forces a break", () => {
  it("runs from 0 at a surf similarity of 0.1 to 1 at 1, logarithmically, with Battjes' 0.4 near 0.6", () => {
    expect(breakingStrength(0.1)).toBe(0);
    expect(breakingStrength(1)).toBe(1);
    expect(breakingStrength(0.4)).toBeCloseTo(Math.log10(4), 12);
    expect(breakingStrength(0.02)).toBe(0);
    expect(breakingStrength(3)).toBe(1);
  });
});

describe("steepness map", () => {
  it("records the steepest the surface got in each cell over the run", () => {
    // 3 x 2 cells of 3 m would be too small for central differences; use a 6 x 5 grid.
    const nx = 6;
    const ny = 5;
    const h = header({ nx, ny, frame_count: 2, times: [0, 2] });
    const cells = nx * ny;
    const frames: number[] = [];
    // Frame 0 slopes 0.1 along x, frame 1 slopes 0.3 along x (cell size 3 m).
    for (const slope of [0.1, 0.3]) for (let n = 0; n < cells; n++) frames.push(slope * 3 * (n % nx));
    const run = runFromBuffers(h, buf(new Array(cells).fill(-50)), buf(frames));
    const map = run.steepnessMap();
    expect(map[2 * nx + 2]).toBeCloseTo(0.3, 6);
    expect(map[1 * nx + 3]).toBeCloseTo(0.3, 6);
    // The edge rows and columns have no central difference.
    expect(map[0]).toBe(0);
    expect(map[2 * nx + 0]).toBe(0);
    expect(map[(ny - 1) * nx + 2]).toBe(0);
    expect(run.steepnessMap()).toBe(map); // computed once
  });

  it("ignores shallow water, where a shoreline is not a wave face", () => {
    const nx = 5;
    const ny = 5;
    const h = header({ nx, ny, frame_count: 1, times: [0] });
    const cells = nx * ny;
    // A steep surface, but the bed comes up to within 10 cm of it everywhere.
    const eta = Array.from({ length: cells }, (_, n) => 1.0 * (n % nx));
    const bed = eta.map((e) => e - 0.1);
    const run = runFromBuffers(h, buf(bed), buf(eta));
    expect(Array.from(run.steepnessMap()).every((s) => s === 0)).toBe(true);
  });
});
