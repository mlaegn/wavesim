import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import type { Run } from "./run";

export type Preset = "oblique" | "top" | "beach";
export type Overlay = "none" | "breaks" | "steep";

const WATER_VERTEX = /* glsl */ `
attribute vec2 cell;

uniform highp sampler2D uEta;
uniform highp sampler2D uBed;
uniform highp sampler2D uBreak;
uniform highp sampler2D uMap;
uniform vec2 uSize;
uniform vec2 uSpacing;
uniform float uExag;
uniform float uLevel;

varying vec3 vWorld;
varying float vDepth;
varying float vHeight;
varying float vBreak;
varying float vMap;
varying vec2 vSlope;
varying vec2 vXY;

// Surface elevation of a neighbour, or the centre's own value if the neighbour is dry,
// so the shoreline does not read as a steep slope.
float wetEta(ivec2 p, float fallback) {
  p = clamp(p, ivec2(0), ivec2(uSize) - 1);
  float e = texelFetch(uEta, p, 0).r;
  float b = texelFetch(uBed, p, 0).r;
  return (e - b > 0.02) ? e : fallback;
}

void main() {
  ivec2 c = ivec2(cell + 0.5);
  float eta = texelFetch(uEta, c, 0).r;
  float bed = texelFetch(uBed, c, 0).r;
  vDepth = eta - bed;
  vHeight = eta - uLevel;
  vBreak = texelFetch(uBreak, c, 0).r;
  vMap = texelFetch(uMap, c, 0).r;

  float eL = wetEta(c + ivec2(-1, 0), eta);
  float eR = wetEta(c + ivec2( 1, 0), eta);
  float eD = wetEta(c + ivec2(0, -1), eta);
  float eU = wetEta(c + ivec2(0,  1), eta);
  vSlope = vec2((eR - eL) / (2.0 * uSpacing.x), (eU - eD) / (2.0 * uSpacing.y));

  vec3 p = vec3(position.x, eta * uExag, position.z);
  vWorld = p;
  vXY = vec2(position.x, -position.z);
  gl_Position = projectionMatrix * viewMatrix * vec4(p, 1.0);
}
`;

const WATER_FRAGMENT = /* glsl */ `
uniform float uExag;
uniform vec3 uSunDir;
uniform vec3 uShallow;
uniform vec3 uDeep;
uniform vec3 uFoam;
uniform vec3 uSky;
uniform float uFoamOn;
uniform float uZonesOn;
uniform vec4 uZones;   // offshore zone width, side zone width, domain length, domain width (m)
uniform float uMakerX;
uniform float uLine;   // half-width of the wave-maker line (m)
uniform float uRange;  // surface height that gets the full tint (m)
uniform float uTint;
uniform float uHasBreak; // 1 if the run records where it is breaking
uniform float uMapOn;
uniform float uCropOn;
uniform vec4 uBox;     // the sea: x0, x1, y0, y1 in metres

varying vec3 vWorld;
varying float vDepth;
varying float vHeight;
varying float vBreak;
varying float vMap;
varying vec2 vSlope;
varying vec2 vXY;

void main() {
  if (vDepth < 0.01) discard;
  if (uCropOn > 0.5 && (vXY.x < uBox.x || vXY.y < uBox.z || vXY.y > uBox.w)) discard;

  // Surface normal from the physical slope, with the vertical exaggeration applied.
  vec3 n = normalize(vec3(-uExag * vSlope.x, 1.0, uExag * vSlope.y));
  vec3 v = normalize(cameraPosition - vWorld);
  float ndv = clamp(dot(n, v), 0.0, 1.0);
  float fresnel = 0.02 + 0.98 * pow(1.0 - ndv, 5.0);

  vec3 base = mix(uShallow, uDeep, smoothstep(0.0, 7.0, vDepth));
  float diffuse = clamp(dot(n, uSunDir), 0.0, 1.0);
  vec3 col = base * (0.5 + 0.6 * diffuse);
  vec3 h = normalize(uSunDir + v);
  float spec = pow(clamp(dot(n, h), 0.0, 1.0), 300.0);
  col = mix(col, uSky, fresnel * 0.4) + vec3(spec) * 0.25;

  // Crests lighter, troughs darker, so the wave pattern reads from above.
  float rise = clamp(vHeight / uRange, -1.0, 1.0);
  col *= 1.0 + uTint * rise;
  col = mix(col, uSky, 0.18 * max(rise, 0.0));

  // Foam where the surface is steep (bore fronts) and in the swash right at the shore.
  float steep = length(vSlope);
  // Where the model itself says the wave is breaking, if the run recorded that; otherwise
  // wherever the surface is steep. Plus the swash right at the shore.
  float front = uHasBreak > 0.5 ? smoothstep(0.55, 0.95, vBreak) : smoothstep(0.05, 0.13, steep);
  float foam = uFoamOn * max(front, 0.55 * smoothstep(0.35, 0.0, vDepth));
  col = mix(col, uFoam, foam);

  if (uZonesOn > 0.5) {
    bool inZone = vXY.x < uZones.x || vXY.y < uZones.y || vXY.y > uZones.w - uZones.y;
    if (inZone) col = mix(col, vec3(0.55, 0.6, 0.66), 0.5);
    if (abs(vXY.x - uMakerX) < uLine) col = vec3(1.0, 0.82, 0.25);
  }

  float alpha = mix(0.6, 0.93, smoothstep(0.0, 4.0, vDepth));
  alpha = max(alpha, foam * 0.95);
  alpha *= smoothstep(0.01, 0.15, vDepth);
  if (uMapOn > 0.5) {
    // How much of the run each place spent breaking: yellow for now and then, red for always.
    float m = smoothstep(0.0, 0.2, vMap);
    vec3 heat = mix(vec3(1.0, 0.85, 0.2), vec3(0.9, 0.12, 0.1), smoothstep(0.1, 0.6, vMap));
    col = mix(col, heat, 0.9 * m);
    alpha = max(alpha, 0.97 * m);
  }
  gl_FragColor = vec4(col, alpha);
  #include <colorspace_fragment>
}
`;

/** Terrain colour for a bed elevation in metres above mean sea level. */
const PALETTE: [number, string][] = [
  [-15, "#0f2538"],
  [-8, "#1c4a69"],
  [-3, "#37798c"],
  [-0.6, "#86b5a6"],
  [0, "#dbd1a8"],
  [1.5, "#cdbd8a"],
  [3, "#92a05f"],
  [6, "#628047"],
  [10, "#4b6339"],
];

function terrainColour(z: number, out: THREE.Color): THREE.Color {
  const first = PALETTE[0] as [number, string];
  if (z <= first[0]) return out.set(first[1]);
  for (let s = 1; s < PALETTE.length; s++) {
    const [z1, c1] = PALETTE[s] as [number, string];
    if (z <= z1) {
      const [z0, c0] = PALETTE[s - 1] as [number, string];
      return out.set(c0).lerp(new THREE.Color(c1), (z - z0) / (z1 - z0));
    }
  }
  return out.set((PALETTE[PALETTE.length - 1] as [number, string])[1]);
}

/**
 * The grid the terrain and the water both sit on. The model's `(x, y, up)` becomes
 * three.js `(X, Y, Z) = (x, up, -y)`, which keeps the coordinate system right-handed.
 */
function gridGeometry(run: Run, bedHeights: boolean): THREE.BufferGeometry {
  const { nx, ny, dx, dy } = run.header;
  const position = new Float32Array(nx * ny * 3);
  const cell = new Float32Array(nx * ny * 2);
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const n = j * nx + i;
      position[3 * n] = (i + 0.5) * dx;
      position[3 * n + 1] = bedHeights ? (run.bed[n] as number) : 0;
      position[3 * n + 2] = -(j + 0.5) * dy;
      cell[2 * n] = i;
      cell[2 * n + 1] = j;
    }
  }
  // Two triangles per quad, counter-clockwise seen from above.
  const index = new Uint32Array(6 * (nx - 1) * (ny - 1));
  let q = 0;
  for (let j = 0; j < ny - 1; j++) {
    for (let i = 0; i < nx - 1; i++) {
      const a = j * nx + i;
      const b = a + 1;
      const c = a + nx;
      const d = c + 1;
      index.set([a, b, d, a, d, c], q);
      q += 6;
    }
  }
  const g = new THREE.BufferGeometry();
  g.setAttribute("position", new THREE.BufferAttribute(position, 3));
  g.setAttribute("cell", new THREE.BufferAttribute(cell, 2));
  g.setIndex(new THREE.BufferAttribute(index, 1));
  return g;
}

function floatTexture(data: Float32Array, nx: number, ny: number): THREE.DataTexture {
  const t = new THREE.DataTexture(data, nx, ny, THREE.RedFormat, THREE.FloatType);
  t.minFilter = THREE.NearestFilter;
  t.magFilter = THREE.NearestFilter;
  t.generateMipmaps = false;
  t.needsUpdate = true;
  return t;
}

export class RunView {
  private readonly renderer: THREE.WebGLRenderer;
  private readonly scene = new THREE.Scene();
  private readonly camera = new THREE.PerspectiveCamera(42, 1, 1, 30000);
  private readonly controls: OrbitControls;
  private readonly terrain: THREE.Mesh;
  private readonly water: THREE.Mesh;
  private readonly etaData: Float32Array;
  private readonly etaTexture: THREE.DataTexture;
  private readonly bedTexture: THREE.DataTexture;
  private readonly breakData: Float32Array;
  private readonly breakTexture: THREE.DataTexture;
  private readonly mapData: Float32Array;
  private readonly mapTexture: THREE.DataTexture;
  private readonly cropPlanes: THREE.Plane[];
  private readonly transect: THREE.Line;
  private readonly uniforms: Record<string, THREE.IUniform>;
  private readonly lengthX: number;
  private readonly lengthY: number;

  constructor(
    private readonly canvas: HTMLCanvasElement,
    readonly run: Run,
  ) {
    const { nx, ny, dx, dy } = run.header;
    this.lengthX = nx * dx;
    this.lengthY = ny * dy;

    this.renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.setClearColor(0x000000, 0);
    this.renderer.localClippingEnabled = true;

    // Terrain.
    const terrainGeometry = gridGeometry(run, true);
    const colours = new Float32Array(nx * ny * 3);
    const tmp = new THREE.Color();
    for (let n = 0; n < nx * ny; n++) {
      terrainColour(run.bed[n] as number, tmp);
      colours.set([tmp.r, tmp.g, tmp.b], 3 * n);
    }
    terrainGeometry.setAttribute("color", new THREE.BufferAttribute(colours, 3));
    terrainGeometry.computeVertexNormals();
    // The terrain is cut to the sea box with clipping planes, so the strips that exist only
    // for the numerics do not show. (World X is x; world Z is minus y.)
    const box = run.viewBox;
    this.cropPlanes = [
      new THREE.Plane(new THREE.Vector3(1, 0, 0), -box.x0),
      new THREE.Plane(new THREE.Vector3(0, 0, -1), -box.y0),
      new THREE.Plane(new THREE.Vector3(0, 0, 1), box.y1),
    ];
    this.terrain = new THREE.Mesh(
      terrainGeometry,
      new THREE.MeshStandardMaterial({
        vertexColors: true,
        roughness: 0.95,
        metalness: 0,
        clippingPlanes: this.cropPlanes,
      }),
    );
    this.scene.add(this.terrain);

    const sun = new THREE.DirectionalLight(0xfff2dc, 2.6);
    sun.position.set(-0.45, 0.75, 0.5).multiplyScalar(this.lengthX);
    this.scene.add(sun, new THREE.HemisphereLight(0xcfe6ff, 0x6b5b45, 1.1));

    // Water: the surface is a texture, displaced in the vertex shader.
    this.etaData = new Float32Array(nx * ny);
    this.etaTexture = floatTexture(this.etaData, nx, ny);
    this.bedTexture = floatTexture(run.bed, nx, ny);
    this.breakData = new Float32Array(nx * ny);
    this.breakTexture = floatTexture(this.breakData, nx, ny);
    this.mapData = new Float32Array(nx * ny);
    this.mapTexture = floatTexture(this.mapData, nx, ny);
    const waves = run.header.waves;
    const cellsToMetres = (key: string, spacing: number) => {
      const v = waves[key];
      return typeof v === "number" ? v * spacing : 0;
    };
    const makerX = typeof waves.maker_x_m === "number" ? waves.maker_x_m : -1;
    this.uniforms = {
      uEta: { value: this.etaTexture },
      uBed: { value: this.bedTexture },
      uBreak: { value: this.breakTexture },
      uMap: { value: this.mapTexture },
      uHasBreak: { value: run.hasField("breaking") ? 1 : 0 },
      uMapOn: { value: 0 },
      uCropOn: { value: 1 },
      uBox: { value: new THREE.Vector4(box.x0, box.x1, box.y0, box.y1) },
      uSize: { value: new THREE.Vector2(nx, ny) },
      uSpacing: { value: new THREE.Vector2(dx, dy) },
      uExag: { value: 2 },
      uLevel: { value: typeof waves.tide_m === "number" ? waves.tide_m : 0 },
      uRange: { value: typeof waves.wave_height_m === "number" ? Math.max(waves.wave_height_m / 2, 0.05) : 0.5 },
      uTint: { value: 0.55 },
      uSunDir: { value: sun.position.clone().normalize() },
      uShallow: { value: new THREE.Color("#52c2b8") },
      uDeep: { value: new THREE.Color("#0a3a66") },
      uFoam: { value: new THREE.Color("#f2f8fa") },
      uSky: { value: new THREE.Color("#cfe4f2") },
      uFoamOn: { value: 1 },
      uZonesOn: { value: 0 },
      uZones: {
        value: new THREE.Vector4(
          cellsToMetres("sponge_offshore_cells", dx),
          cellsToMetres("sponge_side_cells", dy),
          this.lengthX,
          this.lengthY,
        ),
      },
      uMakerX: { value: makerX },
      uLine: { value: 1.5 * dx },
    };
    this.water = new THREE.Mesh(
      gridGeometry(run, false),
      new THREE.ShaderMaterial({
        vertexShader: WATER_VERTEX,
        fragmentShader: WATER_FRAGMENT,
        uniforms: this.uniforms,
        transparent: true,
      }),
    );
    this.water.frustumCulled = false; // the shader moves the vertices
    this.scene.add(this.water);

    // A line on the water marking where the side view is taken.
    this.transect = new THREE.Line(
      new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(box.x0, 0, 0), new THREE.Vector3(box.x1, 0, 0)]),
      new THREE.LineBasicMaterial({ color: 0xffd23f, depthTest: false, transparent: true, opacity: 0.9 }),
    );
    this.transect.renderOrder = 10;
    this.transect.frustumCulled = false;
    this.scene.add(this.transect);

    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.08;
    this.controls.maxPolarAngle = Math.PI * 0.495;
    this.setPreset("oblique");
    this.setTime(0);
    this.setExaggeration(2);
    this.resize();
  }

  private uniform(name: string): THREE.IUniform {
    const u = this.uniforms[name];
    if (!u) throw new Error(`no uniform ${name}`);
    return u;
  }

  /** Show the surface at `t` seconds into the run. */
  setTime(t: number): void {
    this.run.sampleInto(t, this.etaData);
    this.etaTexture.needsUpdate = true;
    if (this.run.hasField("breaking")) {
      this.run.sampleFieldInto("breaking", t, this.breakData);
      this.breakTexture.needsUpdate = true;
    }
  }

  /** The free-surface elevation being shown now. */
  get eta(): Float32Array {
    return this.etaData;
  }

  /** How close to breaking each cell is now, or `null` if the run did not record it. */
  get breaking(): Float32Array | null {
    return this.run.hasField("breaking") ? this.breakData : null;
  }

  /** Hide the strips that exist only for the numerics: the wave-maker's bump and the sides. */
  setCrop(on: boolean): void {
    this.uniform("uCropOn").value = on ? 1 : 0;
    const material = this.terrain.material as THREE.MeshStandardMaterial;
    material.clippingPlanes = on ? this.cropPlanes : [];
    material.needsUpdate = true;
  }

  /**
   * Colour the water by what happened there over the whole run: how much of it was spent
   * breaking, or how steep the surface ever got (0.06 and below is clear, 0.4 and above is
   * full red).
   */
  setOverlay(kind: Overlay): void {
    const map =
      kind === "breaks"
        ? this.run.breakMap()
        : kind === "steep"
          ? this.run.steepnessMap().map((s) => Math.min(Math.max((s - 0.06) / 0.34, 0), 1))
          : null;
    this.mapData.set(map ?? new Float32Array(this.mapData.length));
    this.mapTexture.needsUpdate = true;
    this.uniform("uMapOn").value = map ? 1 : 0;
  }

  /** Mark the line the side view is taken along, at `y` metres. */
  setTransect(y: number): void {
    this.transect.position.z = -y;
    this.transect.position.y = 2 * (this.uniform("uExag").value as number);
  }

  /** Stretch heights (bed and water together) to make small waves visible. */
  setExaggeration(factor: number): void {
    this.terrain.scale.y = factor;
    this.uniform("uExag").value = factor;
    this.transect.position.y = 2 * factor;
  }

  setFoam(on: boolean): void {
    this.uniform("uFoamOn").value = on ? 1 : 0;
  }

  /** Mark the wave-maker line and the absorbing zones, which are numerical, not sea. */
  setZones(on: boolean): void {
    this.uniform("uZonesOn").value = on ? 1 : 0;
  }

  setPreset(preset: Preset): void {
    const { x0, x1, y0, y1 } = this.run.viewBox;
    const shore = Math.min(this.run.shoreline, x1);
    const sea = Math.max(shore - x0, 60); // how far the sea reaches from the wave-maker to the shore
    const width = y1 - y0;
    const midZ = -(y0 + y1) / 2;
    const e = this.uniform("uExag").value as number; // heights are stretched by this
    const target = new THREE.Vector3();
    const position = new THREE.Vector3();
    switch (preset) {
      case "oblique": // from the south-west corner, along the break
        target.set(x0 + 0.75 * sea, 0, -(y0 + 0.3 * width));
        position.set(x0 - 0.1 * sea, (0.28 * sea + 12) * Math.max(e / 3, 1), -y0 + 0.45 * sea);
        break;
      case "top": {
        // The shore runs up the screen, so the alongshore width is the vertical extent.
        const reach = Math.max(width, (sea + 80) / this.camera.aspect);
        const height = (0.55 * reach) / Math.tan(THREE.MathUtils.degToRad(this.camera.fov / 2));
        target.set(x0 + 0.5 * (sea + 80), 0, midZ);
        position.set(target.x, height, midZ + 1);
        break;
      }
      case "beach": // standing at the waterline, looking out along the lineup
        target.set(x0 + 0.3 * sea, 0, -(y0 + 0.62 * width));
        position.set(shore - 0.04 * sea, 5 * e, -(y0 + 0.18 * width));
        break;
    }
    this.camera.position.copy(position);
    this.controls.target.copy(target);
    this.controls.update();
  }

  resize(): void {
    const w = this.canvas.clientWidth;
    const h = this.canvas.clientHeight;
    if (w === 0 || h === 0) return;
    this.renderer.setSize(w, h, false);
    this.camera.aspect = w / h;
    this.camera.updateProjectionMatrix();
  }

  render(): void {
    this.controls.update();
    this.renderer.render(this.scene, this.camera);
  }

  dispose(): void {
    this.controls.dispose();
    this.terrain.geometry.dispose();
    this.water.geometry.dispose();
    (this.terrain.material as THREE.Material).dispose();
    (this.water.material as THREE.Material).dispose();
    this.etaTexture.dispose();
    this.bedTexture.dispose();
    this.breakTexture.dispose();
    this.mapTexture.dispose();
    this.transect.geometry.dispose();
    (this.transect.material as THREE.Material).dispose();
    this.renderer.dispose();
  }
}
