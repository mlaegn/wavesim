import type { Run } from "./run";

export interface Readout {
  /** The steepest face of the surface along the line, where there is real water. */
  steepest: { slope: number; x: number; depth: number } | null;
}

const lerp = (a: number, b: number, t: number) => a + (b - a) * t;

function mix(c0: [number, number, number], c1: [number, number, number], t: number): string {
  const u = Math.min(Math.max(t, 0), 1);
  return `rgb(${Math.round(lerp(c0[0], c1[0], u))}, ${Math.round(lerp(c0[1], c1[1], u))}, ${Math.round(lerp(c0[2], c1[2], u))})`;
}

/**
 * Draw the wave along one line across the break, side-on: the bed, the water, and the
 * surface coloured by how close to breaking the switch says it is, and a white line where the
 * waves break over the run (where they are tallest). `row` is the cell row
 * (the position along the shore). Returns the numbers worth reading off it.
 */
export function drawProfile(
  canvas: HTMLCanvasElement,
  run: Run,
  eta: Float32Array,
  breaking: Float32Array | null,
  row: number,
  from: number,
  to: number,
): Readout {
  const { nx, dx } = run.header;
  const dpr = window.devicePixelRatio || 1;
  const w = canvas.clientWidth;
  const h = canvas.clientHeight;
  if (w === 0 || h === 0) return { steepest: null };
  if (canvas.width !== Math.round(w * dpr) || canvas.height !== Math.round(h * dpr)) {
    canvas.width = Math.round(w * dpr);
    canvas.height = Math.round(h * dpr);
  }
  const g = canvas.getContext("2d");
  if (!g) return { steepest: null };
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.clearRect(0, 0, w, h);

  const base = row * nx;
  const i0 = Math.max(0, Math.floor(from / dx));
  const i1 = Math.min(nx - 1, Math.ceil(to / dx));
  let lowest = 0;
  for (let i = i0; i <= i1; i++) lowest = Math.min(lowest, run.bed[base + i] as number);
  const zMin = Math.min(Math.max(lowest, -10), -3);
  const zMax = 4;
  const m = { left: 34, right: 8, top: 8, bottom: 20 };
  const sx = (w - m.left - m.right) / (i1 - i0 || 1);
  const X = (i: number) => m.left + (i - i0) * sx;
  const Z = (z: number) => m.top + ((zMax - Math.min(Math.max(z, zMin), zMax)) / (zMax - zMin)) * (h - m.top - m.bottom);

  // Sea level and the grid.
  g.font = "10px system-ui, sans-serif";
  g.fillStyle = "rgba(255,255,255,0.55)";
  g.strokeStyle = "rgba(255,255,255,0.12)";
  g.lineWidth = 1;
  for (let z = Math.ceil(zMin / 2) * 2; z <= zMax; z += 2) {
    g.beginPath();
    g.moveTo(m.left, Z(z));
    g.lineTo(w - m.right, Z(z));
    g.stroke();
    g.fillText(`${z}`, 6, Z(z) + 3);
  }
  const x0 = Math.ceil((i0 * dx) / 50) * 50;
  for (let x = x0; x <= i1 * dx; x += 50) {
    g.fillText(`${x}`, X(x / dx) - 10, h - 6);
  }

  // The bed.
  g.beginPath();
  g.moveTo(X(i0), h - m.bottom);
  for (let i = i0; i <= i1; i++) g.lineTo(X(i), Z(run.bed[base + i] as number));
  g.lineTo(X(i1), h - m.bottom);
  g.closePath();
  g.fillStyle = "#a89468";
  g.fill();

  // The water, as one translucent column per cell.
  g.fillStyle = "rgba(60, 140, 200, 0.55)";
  for (let i = i0; i < i1; i++) {
    const e0 = eta[base + i] as number;
    const e1 = eta[base + i + 1] as number;
    const b0 = run.bed[base + i] as number;
    const b1 = run.bed[base + i + 1] as number;
    if (e0 - b0 < 0.01 || e1 - b1 < 0.01) continue;
    g.beginPath();
    g.moveTo(X(i), Z(e0));
    g.lineTo(X(i + 1), Z(e1));
    g.lineTo(X(i + 1), Z(b1));
    g.lineTo(X(i), Z(b0));
    g.closePath();
    g.fill();
  }

  // The surface, coloured by how close to breaking it is.
  let steepest: Readout["steepest"] = null;
  for (let i = i0; i < i1; i++) {
    const e0 = eta[base + i] as number;
    const e1 = eta[base + i + 1] as number;
    const depth = e0 - (run.bed[base + i] as number);
    if (depth < 0.01 || e1 - (run.bed[base + i + 1] as number) < 0.01) continue;
    const b = breaking ? (breaking[base + i] as number) : 0;
    g.strokeStyle = mix([200, 235, 255], [255, 92, 48], (b - 0.3) / 0.6);
    g.lineWidth = 1.5 + 2 * Math.min(Math.max((b - 0.3) / 0.6, 0), 1);
    g.beginPath();
    g.moveTo(X(i), Z(e0));
    g.lineTo(X(i + 1), Z(e1));
    g.stroke();

    if (i > i0 && depth > 0.3) {
      const slope = Math.abs((eta[base + i + 1] as number) - (eta[base + i - 1] as number)) / (2 * dx);
      if (!steepest || slope > steepest.slope) steepest = { slope, x: (i + 0.5) * dx, depth };
    }
  }

  // Mark the steepest point, and where the waves break on this line over the run: where they
  // are tallest, as the run records it.
  if (steepest) {
    const i = steepest.x / dx - 0.5;
    g.fillStyle = "#ffd23f";
    g.beginPath();
    g.arc(X(i), Z(eta[base + Math.round(i)] as number), 4, 0, Math.PI * 2);
    g.fill();
  }
  const point = run.breakLine[row];
  if (point) {
    const i = point.x / dx - 0.5;
    g.strokeStyle = "#ffffff";
    g.lineWidth = 2;
    g.beginPath();
    g.moveTo(X(i), m.top);
    g.lineTo(X(i), h - m.bottom);
    g.stroke();
    g.fillStyle = "#ffffff";
    g.fillText("breaks", X(i) + 4, m.top + 10);
  }
  return { steepest };
}
