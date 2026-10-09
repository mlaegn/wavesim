import { describe, expect, it } from "vitest";
import { parseSliceHeader, SliceError, sliceFromBuffers, speedColour } from "../src/slice";

/** A two-frame slice of three nodes, as `wavesim slice` writes it. */
function sample(): { header: Record<string, unknown>; frames: ArrayBuffer } {
  const header = {
    format: "wavesim-slice",
    version: 1,
    name: "slice-test",
    bed: [
      [0, -1],
      [5, -1],
      [10, -0.2],
    ],
    nodes: 3,
    frames: "frames.f32",
    fields: ["x", "z", "speed"],
    dtype: "<f4",
    frame_count: 2,
    times: [0, 0.5],
    wave: { depth: 1, height: 0.3, g: 9.81 },
    ended: "landed",
    landing: { time: 0.5, x: 6, z: 0, throw: 0.5, tube_area: 0.1, tube_width: 0.4, tube_height: 0.3 },
  };
  // Field after field within each frame: x, then z, then speed.
  const values = [0, 5, 10, 0, 0.3, 0, 0, 2, 0, 0, 6, 10, 0, 0.5, 0.1, 0, 4, 1];
  return { header, frames: new Float32Array(values).buffer };
}

describe("slice files", () => {
  it("read back field by field, frame by frame", () => {
    const { header, frames } = sample();
    const s = sliceFromBuffers(header, frames);
    expect(Array.from(s.field(1, 0))).toEqual([0, 6, 10]);
    expect(Array.from(s.field(0, 2))).toEqual([0, 2, 0]);
    expect(s.duration).toBe(0.5);
    expect(s.header.landing?.throw).toBe(0.5);
    expect(s.waveSpeed).toBeCloseTo(Math.sqrt(9.81 * 1.3), 12);
  });

  it("give the surface between frames by a straight line, and hold it past the ends", () => {
    const { header, frames } = sample();
    const s = sliceFromBuffers(header, frames);
    const mid = s.at(0.25);
    expect(mid.x[1]).toBeCloseTo(5.5, 6);
    expect(mid.z[1]).toBeCloseTo(0.4, 6);
    expect(mid.speed[1]).toBeCloseTo(3, 6);
    expect(s.at(-1).x[1]).toBe(5);
    expect(s.at(9).x[1]).toBe(6);
  });

  it("are rejected with a reason when they do not fit together", () => {
    const { header, frames } = sample();
    expect(() => sliceFromBuffers(header, frames.slice(0, 8))).toThrow(SliceError);
    expect(() => parseSliceHeader({ ...header, times: [0] })).toThrow(/one number per frame/);
    expect(() => parseSliceHeader({ ...header, times: [0.5, 0] })).toThrow(/backwards/);
    expect(() => parseSliceHeader({ ...header, format: "wavesim-run" })).toThrow(/wavesim-slice/);
    expect(() => parseSliceHeader({ ...header, landing: { time: 1 } })).toThrow(/landing\.x/);
  });
});

describe("the speed colour", () => {
  it("runs pale blue when still, white at 0.6 of the wave's speed, red beyond it", () => {
    expect(speedColour(0)).toBe("rgb(160, 210, 235)");
    expect(speedColour(0.6)).toBe("rgb(255, 255, 255)");
    expect(speedColour(1.2)).toBe("rgb(215, 25, 20)");
    expect(speedColour(5)).toBe("rgb(215, 25, 20)");
    expect(speedColour(-1)).toBe("rgb(160, 210, 235)");
  });
});
