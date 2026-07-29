/* ============================================================
   Sol Atlas
   The solar system, live and to scale.

   Data:
     - NASA/JPL Small-Body Database Query API  (orbital elements)
       https://ssd-api.jpl.nasa.gov/doc/sbdb_query.html
     - NASA/JPL CNEOS Close-Approach Data API
       https://ssd-api.jpl.nasa.gov/doc/cad.html

   Rendering: hand-rolled WebGL1, zero dependencies.
   Positions are propagated from Keplerian elements in the
   heliocentric ecliptic J2000 frame, solved per frame.
   ============================================================ */
"use strict";

/* ============================================================
   1. Orbital mechanics
   ============================================================ */
const TWO_PI = Math.PI * 2;
const DEG = Math.PI / 180;
const J2000 = 2451545.0;            // JD of epoch J2000.0
const GAUSS_K = 0.01720209895;      // Gaussian gravitational constant, rad/day

function jdNow() {
  return Date.now() / 86400000 + 2440587.5;
}

/* Solve Kepler's equation  E - e·sinE = M  (radians) via Newton. */
function solveKepler(M, e) {
  M = M % TWO_PI;
  if (M > Math.PI) M -= TWO_PI;
  else if (M < -Math.PI) M += TWO_PI;
  let E = e < 0.8 ? M : Math.PI * (M < 0 ? -1 : 1);
  for (let j = 0; j < 12; j++) {
    const s = Math.sin(E);
    const c = Math.cos(E);
    const d = (E - e * s - M) / (1 - e * c);
    E -= d;
    if (Math.abs(d) < 1e-9) break;
  }
  return E;
}

/* Perifocal→ecliptic basis vectors P (periapsis dir) and Q, from
   argument of periapsis w, longitude of node om, inclination i (rad). */
function perifocalBasis(w, om, i, out) {
  const cw = Math.cos(w), sw = Math.sin(w);
  const co = Math.cos(om), so = Math.sin(om);
  const ci = Math.cos(i), si = Math.sin(i);
  out[0] = co * cw - so * sw * ci;   // Px
  out[1] = so * cw + co * sw * ci;   // Py
  out[2] = sw * si;                  // Pz
  out[3] = -co * sw - so * cw * ci;  // Qx
  out[4] = -so * sw + co * cw * ci;  // Qy
  out[5] = cw * si;                  // Qz
  return out;
}

/* Position on the ellipse from eccentric anomaly. */
function ellipsePoint(a, e, b, P, E, out) {
  const xp = a * (Math.cos(E) - e);
  const yp = b * Math.sin(E);
  out[0] = xp * P[0] + yp * P[3];
  out[1] = xp * P[1] + yp * P[4];
  out[2] = xp * P[2] + yp * P[5];
  return out;
}

/* ---- Major planets ----
   Keplerian elements + per-century rates, valid 1800–2050 AD.
   Source: JPL "Approximate Positions of the Planets" (Table 1),
   https://ssd.jpl.nasa.gov/planets/approx_pos.html
   Order: a [au], e, I [deg], L [deg], long.peri [deg], long.node [deg] */
const PLANETS = [
  { name: "Mercury", color: [0.72, 0.69, 0.65], size: 0.035,
    el: [0.38709927, 0.20563593, 7.00497902, 252.25032350, 77.45779628, 48.33076593],
    rate: [0.00000037, 0.00001906, -0.00594749, 149472.67411175, 0.16047689, -0.12534081] },
  { name: "Venus", color: [0.96, 0.85, 0.63], size: 0.055,
    el: [0.72333566, 0.00677672, 3.39467605, 181.97909950, 131.60246718, 76.67984255],
    rate: [0.00000390, -0.00004107, -0.00078890, 58517.81538729, 0.00268329, -0.27769418] },
  { name: "Earth", color: [0.44, 0.71, 1.0], size: 0.058,
    el: [1.00000261, 0.01671123, -0.00001531, 100.46457166, 102.93768193, 0.0],
    rate: [0.00000562, -0.00004392, -0.01294668, 35999.37244981, 0.32327364, 0.0] },
  { name: "Mars", color: [1.0, 0.54, 0.36], size: 0.045,
    el: [1.52371034, 0.09339410, 1.84969142, -4.55343205, -23.94362959, 49.55953891],
    rate: [0.00001847, 0.00007882, -0.00813131, 19140.30268499, 0.44441088, -0.29257343] },
  { name: "Jupiter", color: [0.91, 0.72, 0.54], size: 0.16,
    el: [5.20288700, 0.04838624, 1.30439695, 34.39644051, 14.72847983, 100.47390909],
    rate: [-0.00011607, -0.00013253, -0.00183714, 3034.74612775, 0.21252668, 0.20469106] },
  { name: "Saturn", color: [0.95, 0.84, 0.56], size: 0.14,
    el: [9.53667594, 0.05386179, 2.48599187, 49.95424423, 92.59887831, 113.66242448],
    rate: [-0.00125060, -0.00050991, 0.00193609, 1222.49362201, -0.41897216, -0.28867794] },
  { name: "Uranus", color: [0.61, 0.85, 0.88], size: 0.10,
    el: [19.18916464, 0.04725744, 0.77263783, 313.23810451, 170.95427630, 74.01692503],
    rate: [-0.00196176, -0.00004397, -0.00242939, 428.48202785, 0.40805281, 0.04240589] },
  { name: "Neptune", color: [0.50, 0.60, 1.0], size: 0.10,
    el: [30.06992276, 0.00859048, 1.77004347, -55.12002969, 44.96476227, 131.78422574],
    rate: [0.00026291, 0.00005105, 0.00035372, 218.45945325, -0.32241464, -0.00508664] },
];

/* Planet Keplerian elements {a,e,i,om,w,M} (rad where angular) at JD. */
function planetElements(p, jd, out) {
  const T = (jd - J2000) / 36525;
  const a = p.el[0] + p.rate[0] * T;
  const e = p.el[1] + p.rate[1] * T;
  const i = (p.el[2] + p.rate[2] * T) * DEG;
  const L = (p.el[3] + p.rate[3] * T) * DEG;
  const lp = (p.el[4] + p.rate[4] * T) * DEG;
  const om = (p.el[5] + p.rate[5] * T) * DEG;
  out.a = a; out.e = e; out.i = i; out.om = om;
  out.w = lp - om;
  out.M = L - lp;
  return out;
}

/* Heliocentric ecliptic position of a planet at JD (au). */
function planetPosition(p, jd, out) {
  const el = planetElements(p, jd, { a: 0, e: 0, i: 0, om: 0, w: 0, M: 0 });
  const P = perifocalBasis(el.w, el.om, el.i, new Float64Array(6));
  const E = solveKepler(el.M, el.e);
  return ellipsePoint(el.a, el.e, el.a * Math.sqrt(1 - el.e * el.e), P, E, out);
}

/* Node test hook (no effect in the browser). */
if (typeof module !== "undefined" && module.exports) {
  module.exports = { jdNow, solveKepler, perifocalBasis, ellipsePoint, planetElements, planetPosition, PLANETS, J2000, DEG, GAUSS_K };
}

/* ============================================================
   Everything below runs only in the browser.
   ============================================================ */
if (typeof document !== "undefined") (() => {

/* ============================================================
   2. Configuration / populations
   ============================================================ */
// The JPL SSD APIs don't send CORS headers, so the browser can't query them
// directly. A GitHub Actions workflow (.github/workflows/data-refresh.yml)
// pulls the same queries server-side and commits this snapshot.
const DATA_URL = "data/asteroids.json";

const GROUPS = [
  { key: "NEO", label: "Near-Earth", color: [1.0, 0.42, 0.42], css: "#ff6b6b",
    desc: "Asteroids whose orbits come within 1.3 au of the Sun — the Atira, Aten, Apollo and Amor classes. Some pass close to Earth." },
  { key: "MCA", label: "Mars-crossers", color: [1.0, 0.66, 0.30], css: "#ffa94d",
    desc: "Asteroids whose orbits cross the orbit of Mars." },
  { key: "MBA", label: "Main belt", color: [1.0, 0.83, 0.47], css: "#ffd479",
    desc: "The classical asteroid belt between Mars and Jupiter, roughly 2.0–3.3 au from the Sun. Home to most known asteroids." },
  { key: "TJN", label: "Jupiter Trojans", color: [0.37, 0.92, 0.83], css: "#5eead4",
    desc: "Asteroids sharing Jupiter's orbit, held in two camps 60° ahead of and behind the planet (the L4 / L5 Lagrange points)." },
  { key: "CEN", label: "Centaurs", color: [0.75, 0.52, 0.99], css: "#c084fc",
    desc: "Icy bodies on unstable orbits between Jupiter and Neptune — part asteroid, part comet." },
  { key: "TNO", label: "Trans-Neptunian", color: [0.49, 0.65, 1.0], css: "#7da7ff",
    desc: "Everything orbiting beyond Neptune: the Kuiper Belt and the scattered disk." },
  { key: "DWF", label: "Dwarf planets", color: [0.83, 0.76, 1.0], css: "#d4c2ff",
    desc: "Worlds massive enough for gravity to pull them round, but which never cleared their orbits — Ceres, Pluto, Eris, Haumea, Makemake and the leading candidates." },
  { key: "COM", label: "Comets", color: [0.55, 0.93, 1.0], css: "#8ce9ff",
    desc: "Periodic comets on closed orbits — icy nuclei that grow comas and tails when they near the Sun." },
  { key: "LPC", label: "Long-period comets", color: [0.42, 0.78, 0.92], css: "#6bc7eb",
    desc: "Near-parabolic comets (e ≥ 0.995) falling in from thousands of au — the observational evidence for the Oort cloud. Select one and follow its orbit outward." },
  { key: "ISO", label: "Interstellar", color: [1.0, 0.42, 0.85], css: "#ff6bd9",
    desc: "Visitors from other star systems on unbound hyperbolic paths — they pass through once and never return." },
  { key: "OTH", label: "Other", color: [0.58, 0.64, 0.72], css: "#94a3b8",
    desc: "Objects that don't fit any of the classes above." },
];
// key → group index, so we don't hard-code positions when buckets change.
const GI = Object.fromEntries(GROUPS.map((g, i) => [g.key, i]));
const CLASS_TO_GROUP = {
  IEO: GI.NEO, ATE: GI.NEO, APO: GI.NEO, AMO: GI.NEO,
  MCA: GI.MCA,
  IMB: GI.MBA, MBA: GI.MBA, OMB: GI.MBA,
  TJN: GI.TJN,
  CEN: GI.CEN,
  TNO: GI.TNO,
};
const CLASS_NAMES = {
  IEO: "Atira (interior-Earth orbit)", ATE: "Aten near-Earth asteroid",
  APO: "Apollo near-Earth asteroid", AMO: "Amor near-Earth asteroid",
  MCA: "Mars-crossing asteroid", IMB: "Inner main-belt asteroid",
  MBA: "Main-belt asteroid", OMB: "Outer main-belt asteroid",
  TJN: "Jupiter trojan", CEN: "Centaur", TNO: "Trans-Neptunian object",
  AST: "Asteroid", PAA: "Parabolic asteroid", HYA: "Hyperbolic asteroid",
  COM: "Periodic comet", LPC: "Long-period comet", ISO: "Interstellar object",
};

// Planet facts for the info card (diameter km, sidereal rotation hours —
// negative = retrograde). Moon counts come from the loaded moons.json.
const PLANET_FACTS = [
  { d: 4879, rot: 1407.6 },    // Mercury
  { d: 12104, rot: -5832.5 },  // Venus
  { d: 12742, rot: 23.93 },    // Earth
  { d: 6779, rot: 24.62 },     // Mars
  { d: 139820, rot: 9.93 },    // Jupiter
  { d: 116460, rot: 10.66 },   // Saturn
  { d: 50724, rot: -17.24 },   // Uranus
  { d: 49244, rot: 16.11 },    // Neptune
];
// planet / moon name → [img/bodies file, credit]
const PLANET_IMG = [
  ["mercury.jpg", "NASA/JHUAPL · MESSENGER"], ["venus.jpg", "NASA/JPL · Mariner 10"],
  ["earth.jpg", "NASA · Apollo 17"], ["mars.png", "ESA · Rosetta/OSIRIS"],
  ["jupiter.png", "NASA/ESA · Hubble"], ["saturn.jpg", "NASA/JPL/SSI · Cassini"],
  ["uranus.png", "NASA/JPL · Voyager 2"], ["neptune.png", "NASA/JPL · Voyager 2"],
];
const MOON_IMG = {
  moon: ["luna.jpg", "NASA/GSFC · LRO"],
  phobos: ["phobos.jpg", "NASA/JPL/UA · MRO"], deimos: ["deimos.jpg", "NASA/JPL/UA · MRO"],
  io: ["io.jpg", "NASA/JPL · Galileo"], europa: ["europa-moon.png", "NASA/JPL · Galileo"],
  ganymede: ["ganymede.png", "NASA/JPL · Juno"], callisto: ["callisto.png", "NASA/JPL · Galileo"],
  amalthea: ["amalthea.png", "NASA/JPL · Galileo"],
  mimas: ["mimas.jpg", "NASA/JPL/SSI · Cassini"], enceladus: ["enceladus.jpg", "NASA/JPL/SSI · Cassini"],
  tethys: ["tethys.png", "NASA/JPL/SSI · Cassini"], dione: ["dione.jpg", "NASA/JPL/SSI · Cassini"],
  rhea: ["rhea.jpg", "NASA/JPL/SSI · Cassini"], titan: ["titan.jpg", "NASA/JPL/SSI · Cassini"],
  hyperion: ["hyperion.jpg", "NASA/JPL/SSI · Cassini"], iapetus: ["iapetus.jpg", "NASA/JPL/SSI · Cassini"],
  phoebe: ["phoebe-moon.jpg", "NASA/JPL/SSI · Cassini"],
  miranda: ["miranda.png", "NASA/JPL · Voyager 2"], ariel: ["ariel.jpg", "NASA/JPL · Voyager 2"],
  umbriel: ["umbriel.jpg", "NASA/JPL · Voyager 2"], titania: ["titania.png", "NASA/JPL · Voyager 2"],
  oberon: ["oberon.jpg", "NASA/JPL · Voyager 2"],
  triton: ["triton.jpg", "NASA/JPL · Voyager 2"], proteus: ["proteus.jpg", "NASA/JPL · Voyager 2"],
  charon: ["charon.png", "NASA/JHUAPL/SwRI · New Horizons"],
};

// IAU dwarf planets + leading candidates. SBDB often lacks their diameter,
// so we carry known values for display and to scale the preview. Keyed by
// primary designation. These are already in the data — this just labels them.
const DWARFS = {
  "1": { name: "Ceres", diam: 939 },
  "134340": { name: "Pluto", diam: 2377 },
  "136199": { name: "Eris", diam: 2326 },
  "136108": { name: "Haumea", diam: 1560 },
  "136472": { name: "Makemake", diam: 1430 },
  "225088": { name: "Gonggong", diam: 1230 },
  "50000": { name: "Quaoar", diam: 1090 },
  "90377": { name: "Sedna", diam: 1000 },
  "90482": { name: "Orcus", diam: 910 },
  "120347": { name: "Salacia", diam: 850 },
};

// The ~21 small bodies that have genuinely been resolved (visited or radar).
// Everything else uses the procedural preview — we only show a real photo
// where one exists. credits live in img/bodies/manifest.json.
const BODY_IMAGES = {
  ceres: { file: "img/bodies/ceres.jpg", credit: "NASA/JPL-Caltech · Dawn" },
  vesta: { file: "img/bodies/vesta.jpg", credit: "NASA/JPL-Caltech · Dawn" },
  pluto: { file: "img/bodies/pluto.png", credit: "NASA/JHUAPL/SwRI · New Horizons" },
  bennu: { file: "img/bodies/bennu.png", credit: "NASA/Goddard · OSIRIS-REx" },
  ryugu: { file: "img/bodies/ryugu.jpg", credit: "JAXA · Hayabusa2" },
  eros: { file: "img/bodies/eros.jpg", credit: "NASA/JPL · NEAR" },
  itokawa: { file: "img/bodies/itokawa.jpg", credit: "JAXA · Hayabusa" },
  ida: { file: "img/bodies/ida.jpg", credit: "NASA/JPL · Galileo" },
  gaspra: { file: "img/bodies/gaspra.jpg", credit: "NASA/JPL · Galileo" },
  mathilde: { file: "img/bodies/mathilde.jpg", credit: "NASA/JPL · NEAR" },
  lutetia: { file: "img/bodies/lutetia.jpg", credit: "ESA · Rosetta/OSIRIS" },
  steins: { file: "img/bodies/steins.jpg", credit: "ESA · Rosetta/OSIRIS" },
  dinkinesh: { file: "img/bodies/dinkinesh.png", credit: "NASA/Goddard/SwRI · Lucy" },
  dimorphos: { file: "img/bodies/dimorphos.png", credit: "NASA/JHUAPL · DART" },
  arrokoth: { file: "img/bodies/arrokoth.png", credit: "NASA/JHUAPL/SwRI · New Horizons" },
  halley: { file: "img/bodies/halley.jpg", credit: "NASA/NSSDC · 1986" },
  churyumov: { file: "img/bodies/churyumov.jpg", credit: "ESA · Rosetta/NavCam" },
  tempel: { file: "img/bodies/tempel.jpg", credit: "NASA/JPL/UMD · Deep Impact" },
  wild: { file: "img/bodies/wild.jpg", credit: "NASA/JPL · Stardust" },
  hartley: { file: "img/bodies/hartley.jpg", credit: "NASA/JPL/UMD · EPOXI" },
  borrelly: { file: "img/bodies/borrelly.jpg", credit: "NASA/JPL · Deep Space 1" },
};
// primary designation → BODY_IMAGES key. Matched by designation only —
// names collide (e.g. asteroid 2688 Halley ≠ comet 1P/Halley).
const REAL_IMG = {
  "1": "ceres", "4": "vesta", "134340": "pluto", "101955": "bennu",
  "162173": "ryugu", "433": "eros", "25143": "itokawa", "243": "ida",
  "951": "gaspra", "253": "mathilde", "21": "lutetia", "2867": "steins",
  "152830": "dinkinesh", "65803": "dimorphos", "486958": "arrokoth",
  "1P": "halley", "67P": "churyumov", "9P": "tempel", "81P": "wild",
  "103P": "hartley", "19P": "borrelly",
};
const sentryMap = new Map();   // designation → risk record (filled at load)

// Per-object element record stride inside group.el:
// [a, e, b, M0, n, epochD, Px, Py, Pz, Qx, Qy, Qz]
const STRIDE = 12;

/* ============================================================
   3. DOM handles & app state
   ============================================================ */
const $ = (id) => document.getElementById(id);
const canvas = $("scene");
const sunHalo = $("sun-halo");
const labelsEl = $("labels");

const state = {
  simJD: jdNow(),
  playing: true,
  dps: 1 / 86400,          // open at the real present, ticking at REAL rate;
                           // faster rates are a user choice (chips below)
  dir: 1,
  cam: {
    yaw: -1.1, pitch: 0.55, dist: 7.8, tYaw: -1.1, tPitch: 0.55, tDist: 7.8,
    target: new Float64Array(3),   // current look-at point
    fFrom: new Float64Array(3),    // look-at when focus last changed
    fBlend: 1,                     // 0→1 transition into the new focus
  },
  topDown: false,
  savedPitch: 0.55,
  selected: null,          // { group, index }
  selPlanet: null,         // selected planet index
  selMoon: null,           // selected moon index
  selSun: false,           // the Sun is selected
  selSat: null,            // selected satellite index
  selCraft: null,          // selected spacecraft index
  showSun: true,           // Sun population (point, halo, picking)
  focus: null,             // camera target: null=Sun | {planet:i} | {small:[gi,k]}
  camEaseRate: 9,          // lowered during the launch dolly-in, restored after
  showPlanets: true,       // planets population (points, orbits, labels)
  showPHA: false,          // highlight potentially-hazardous asteroids
  trueScale: false,        // render bodies at real angular size (the honest emptiness)
  totalLoaded: 0,
  shownCount: 0,
  needFullUpdate: true,
  updateSlices: 1,
  sliceCursor: 0,
};

const groups = GROUPS.map((g) => ({
  ...g,
  count: 0,
  el: new Float32Array(0),       // packed elements, STRIDE per object
  pos: new Float32Array(0),      // xyz per object, recomputed
  sizes: new Float32Array(0),
  meta: [],                      // { name, cls, a, e, i, H, diam }
  visible: true,
  posBuf: null,
  sizeBuf: null,
  dirty: true,
}));

const seen = new Set();          // pdes dedupe across queries
const phaList = [];              // { gi, k } for potentially-hazardous asteroids

// Planetary moons (loaded from data/moons.json, propagated around parents).
const moons = {
  count: 0, visible: true,
  el: new Float32Array(0),       // same STRIDE layout as asteroid groups
  parentIdx: [],                 // planetState index, or -1 (small-body parent)
  parentSmall: [],               // [gi,k] for small-body parents (Pluto), else null
  pos: new Float32Array(0),
  sizes: new Float32Array(0),
  meta: [],                      // { name, parentName, a, e, i, n, radius, img }
  posBuf: null, sizeBuf: null, dirty: true,
};
let plutoRef = null;             // [gi,k] of Pluto in the asteroid groups

// Active satellites (CelesTrak TLEs, propagated with vendored SGP4). Earth's
// artificial moons — only loaded and drawn when the camera reaches Earth.
const EARTH_IDX = 2;             // planetState index of Earth
const OBLIQ = 23.4392911 * DEG;  // ecliptic obliquity (equatorial → ecliptic)
const AU_KM = 149597870.7;
const sats = {
  loaded: false, loading: false, visible: true,    // payloads show on arrival at Earth — the marquee moment
  debrisVisible: false,          // tracked debris, separate toggle + colour, off by default
  payloadCount: 0,               // split index: [0,payloadCount)=payloads, [payloadCount,count)=debris
  count: 0, recs: [], names: [],
  pos: new Float32Array(0), rel: new Float32Array(0), sizes: new Float32Array(0),
  posBuf: null, sizeBuf: null, dirty: false,
  sliceCursor: 0, fullUpdate: true,
  epochJD: 0,                    // snapshot reference epoch (TLEs valid near here)
};

// In-flight deep-space probes (Horizons trajectory polylines, interpolated).
const craft = [];                // { name, info, jd, verts, n, pos, inRange, trajBuf, posBuf, sizeBuf, labelEl }
const CRAFT_CSS = "#8ef0b0";     // spacecraft population colour (green)
let craftVisible = false;        // off by default, like the PHA highlight
let satCountKnown = 0;           // active payloads (from the tiny count file or full load)
let debrisCountKnown = 0;        // tracked debris
// header splits man-made into "functional" (payloads + in-flight spacecraft) and "debris"
function updateManMade() {
  $("stat-craft").textContent = (satCountKnown + craft.length).toLocaleString("en-US");
  $("stat-debris").textContent = debrisCountKnown.toLocaleString("en-US");
}

/* ============================================================
   4. WebGL renderer
   ============================================================ */
const gl = canvas.getContext("webgl", { antialias: true, alpha: false, powerPreference: "high-performance" });
if (!gl) {
  $("loader-status").textContent = "";
  $("loader-error-msg").textContent = "Your browser does not support WebGL, which this visualization needs.";
  $("loader-error").hidden = false;
  $("retry-btn").hidden = true;
  return;
}

function compile(type, src) {
  const s = gl.createShader(type);
  gl.shaderSource(s, src);
  gl.compileShader(s);
  if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) {
    throw new Error("Shader: " + gl.getShaderInfoLog(s));
  }
  return s;
}
function program(vs, fs) {
  const p = gl.createProgram();
  gl.attachShader(p, compile(gl.VERTEX_SHADER, vs));
  gl.attachShader(p, compile(gl.FRAGMENT_SHADER, fs));
  gl.linkProgram(p);
  if (!gl.getProgramParameter(p, gl.LINK_STATUS)) {
    throw new Error("Program: " + gl.getProgramInfoLog(p));
  }
  return p;
}

const POINT_VS = `
attribute vec3 aPos;
attribute float aSize;
attribute float aPhase;
uniform mat4 uVP;
uniform float uPixScale, uTime, uMinPx, uMaxPx, uCull;
varying float vTw;
void main() {
  gl_Position = uVP * vec4(aPos, 1.0);
  float px = uPixScale * aSize / max(gl_Position.w, 1e-4);
  if (px < uCull) { gl_Position = vec4(2.0, 2.0, 2.0, 1.0); }  // true-scale: cull sub-pixel bodies (drivers floor PointSize at 1px)
  gl_PointSize = clamp(px, uMinPx, uMaxPx);
  vTw = aPhase > 0.0 ? 0.72 + 0.28 * sin(uTime * 1.8 + aPhase * 7.0) : 1.0;
}`;
const POINT_FS = `
precision mediump float;
uniform vec3 uColor;
uniform float uAlpha;
varying float vTw;
void main() {
  vec2 d = gl_PointCoord - 0.5;
  float r2 = dot(d, d) * 4.0;
  float a = exp(-r2 * 3.2) * smoothstep(1.0, 0.55, r2) * uAlpha * vTw;
  gl_FragColor = vec4(uColor * a, a);
}`;
const LINE_VS = `
attribute vec3 aPos;
uniform mat4 uVP;
void main() { gl_Position = uVP * vec4(aPos, 1.0); }`;
const LINE_FS = `
precision mediump float;
uniform vec3 uColor;
uniform float uAlpha;
void main() { gl_FragColor = vec4(uColor * uAlpha, uAlpha); }`;

const ptProg = program(POINT_VS, POINT_FS);
const lnProg = program(LINE_VS, LINE_FS);
const PT = {
  aPos: gl.getAttribLocation(ptProg, "aPos"),
  aSize: gl.getAttribLocation(ptProg, "aSize"),
  aPhase: gl.getAttribLocation(ptProg, "aPhase"),
  uVP: gl.getUniformLocation(ptProg, "uVP"),
  uPixScale: gl.getUniformLocation(ptProg, "uPixScale"),
  uTime: gl.getUniformLocation(ptProg, "uTime"),
  uMinPx: gl.getUniformLocation(ptProg, "uMinPx"),
  uMaxPx: gl.getUniformLocation(ptProg, "uMaxPx"),
  uCull: gl.getUniformLocation(ptProg, "uCull"),
  uColor: gl.getUniformLocation(ptProg, "uColor"),
  uAlpha: gl.getUniformLocation(ptProg, "uAlpha"),
};
const LN = {
  aPos: gl.getAttribLocation(lnProg, "aPos"),
  uVP: gl.getUniformLocation(lnProg, "uVP"),
  uColor: gl.getUniformLocation(lnProg, "uColor"),
  uAlpha: gl.getUniformLocation(lnProg, "uAlpha"),
};

// Textured, Lambert-lit sphere for the focused planet (Sun at world origin).
const SPHERE_VS = `
attribute vec3 aPos; attribute vec3 aNorm; attribute vec2 aUV;
uniform mat4 uMVP, uModel;
varying vec3 vN; varying vec3 vW; varying vec2 vUV;
void main() {
  gl_Position = uMVP * vec4(aPos, 1.0);
  vN = mat3(uModel) * aNorm;            // uniform scale → fine as a normal basis
  vW = (uModel * vec4(aPos, 1.0)).xyz;  // world position; Sun is at the origin
  vUV = aUV;
}`;
const SPHERE_FS = `
precision mediump float;
uniform sampler2D uTex; uniform float uAmbient, uAlpha;
varying vec3 vN; varying vec3 vW; varying vec2 vUV;
void main() {
  vec3 N = normalize(vN);
  vec3 L = normalize(-vW);              // toward the Sun (origin)
  float d = max(dot(N, L), 0.0);
  vec3 c = texture2D(uTex, vUV).rgb * (uAmbient + (1.0 - uAmbient) * d);
  gl_FragColor = vec4(c, uAlpha);
}`;
// Procedural ring disc (no texture): radius set per-planet via uR0/uR1.
const RING_VS = `
attribute vec3 aRing;                   // x=cosθ, y=sinθ, z=t (0 inner → 1 outer)
uniform mat4 uMVP; uniform float uR0, uR1;
varying float vR;
void main() {
  float rad = mix(uR0, uR1, aRing.z);
  gl_Position = uMVP * vec4(aRing.x * rad, aRing.y * rad, 0.0, 1.0);
  vR = aRing.z;
}`;
const RING_FS = `
precision mediump float;
uniform vec3 uColor; uniform float uAlpha, uCassini;
varying float vR;
void main() {
  float edge = smoothstep(0.0, 0.05, vR) * smoothstep(1.0, 0.93, vR);
  float gap = 1.0 - uCassini * (1.0 - smoothstep(0.04, 0.10, abs(vR - 0.55)));
  gl_FragColor = vec4(uColor, uAlpha * edge * gap);
}`;
const spProg = program(SPHERE_VS, SPHERE_FS);
const rnProg = program(RING_VS, RING_FS);
const SP = {
  aPos: gl.getAttribLocation(spProg, "aPos"), aNorm: gl.getAttribLocation(spProg, "aNorm"), aUV: gl.getAttribLocation(spProg, "aUV"),
  uMVP: gl.getUniformLocation(spProg, "uMVP"), uModel: gl.getUniformLocation(spProg, "uModel"),
  uTex: gl.getUniformLocation(spProg, "uTex"), uAmbient: gl.getUniformLocation(spProg, "uAmbient"), uAlpha: gl.getUniformLocation(spProg, "uAlpha"),
};
const RN = {
  aRing: gl.getAttribLocation(rnProg, "aRing"),
  uMVP: gl.getUniformLocation(rnProg, "uMVP"), uR0: gl.getUniformLocation(rnProg, "uR0"), uR1: gl.getUniformLocation(rnProg, "uR1"),
  uColor: gl.getUniformLocation(rnProg, "uColor"), uAlpha: gl.getUniformLocation(rnProg, "uAlpha"), uCassini: gl.getUniformLocation(rnProg, "uCassini"),
};

gl.disable(gl.DEPTH_TEST);
gl.enable(gl.BLEND);
gl.blendFunc(gl.ONE, gl.ONE);   // additive — glowing space dust
gl.clearColor(0, 0, 0, 1);

/* ---- tiny mat4 helpers (column-major) ---- */
function mat4Perspective(out, fovY, aspect, near, far) {
  const f = 1 / Math.tan(fovY / 2);
  out.fill(0);
  out[0] = f / aspect; out[5] = f;
  out[10] = (far + near) / (near - far); out[11] = -1;
  out[14] = (2 * far * near) / (near - far);
  return out;
}
function mat4LookAt(out, eye, target, up) {
  let zx = eye[0] - target[0], zy = eye[1] - target[1], zz = eye[2] - target[2];
  let l = Math.hypot(zx, zy, zz); zx /= l; zy /= l; zz /= l;
  let xx = up[1] * zz - up[2] * zy, xy = up[2] * zx - up[0] * zz, xz = up[0] * zy - up[1] * zx;
  l = Math.hypot(xx, xy, xz) || 1; xx /= l; xy /= l; xz /= l;
  const yx = zy * xz - zz * xy, yy = zz * xx - zx * xz, yz = zx * xy - zy * xx;
  out[0] = xx; out[1] = yx; out[2] = zx; out[3] = 0;
  out[4] = xy; out[5] = yy; out[6] = zy; out[7] = 0;
  out[8] = xz; out[9] = yz; out[10] = zz; out[11] = 0;
  out[12] = -(xx * eye[0] + xy * eye[1] + xz * eye[2]);
  out[13] = -(yx * eye[0] + yy * eye[1] + yz * eye[2]);
  out[14] = -(zx * eye[0] + zy * eye[1] + zz * eye[2]);
  out[15] = 1;
  return out;
}
function mat4Mul(out, a, b) {
  for (let c = 0; c < 4; c++) {
    const b0 = b[c * 4], b1 = b[c * 4 + 1], b2 = b[c * 4 + 2], b3 = b[c * 4 + 3];
    out[c * 4] = a[0] * b0 + a[4] * b1 + a[8] * b2 + a[12] * b3;
    out[c * 4 + 1] = a[1] * b0 + a[5] * b1 + a[9] * b2 + a[13] * b3;
    out[c * 4 + 2] = a[2] * b0 + a[6] * b1 + a[10] * b2 + a[14] * b3;
    out[c * 4 + 3] = a[3] * b0 + a[7] * b1 + a[11] * b2 + a[15] * b3;
  }
  return out;
}

const proj = new Float32Array(16);
const view = new Float32Array(16);
const vp = new Float32Array(16);
const viewSky = new Float32Array(16);   // rotation-only view for the star dome
const vpSky = new Float32Array(16);
const FOV = 55 * DEG;
const MAX_DIST = 150000;                // au — past the Oort cloud's outer edge
// On the way out, once you cross into interstellar scale, settle the camera into
// a composed view so the far signposts (interstellar space, Alpha Centauri) are
// already framed — no dragging around to find them.
const FAR_ORIENT_DIST = 55000;          // au — where those labels start to fade in
const FAR_YAW = 2.45, FAR_PITCH = 0.55; // frames interstellar space + Alpha Centauri

let dpr = 1, cssW = 0, cssH = 0;
function resize() {
  dpr = Math.min(window.devicePixelRatio || 1, 2);
  cssW = canvas.clientWidth || window.innerWidth;
  cssH = canvas.clientHeight || window.innerHeight;
  canvas.width = Math.round(cssW * dpr);
  canvas.height = Math.round(cssH * dpr);
  gl.viewport(0, 0, canvas.width, canvas.height);
  mat4Perspective(proj, FOV, canvas.width / canvas.height, 0.01, 4000);
}
window.addEventListener("resize", resize);
resize();

/* ---- static buffers: stars, sun, planets, orbit lines ---- */
function makeBuffer(data, usage) {
  const b = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, b);
  gl.bufferData(gl.ARRAY_BUFFER, data, usage || gl.STATIC_DRAW);
  return b;
}
function makeIndexBuffer(arr) {
  const e = gl.createBuffer();
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, e);
  gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(arr), gl.STATIC_DRAW);
  return e;
}

/* ---- focused-planet sphere + rings (lazy textures, depth-bracketed draw) ---- */
function buildUVSphere(LAT, LON) {
  const pos = [], nrm = [], uv = [], idx = [];
  for (let y = 0; y <= LAT; y++) {
    const v = y / LAT, th = v * Math.PI, st = Math.sin(th), ct = Math.cos(th);
    for (let x = 0; x <= LON; x++) {
      const u = x / LON, ph = u * TWO_PI, nx = Math.cos(ph) * st, ny = Math.sin(ph) * st, nz = ct;
      pos.push(nx, ny, nz); nrm.push(nx, ny, nz); uv.push(u, v);
    }
  }
  const row = LON + 1;
  for (let y = 0; y < LAT; y++) for (let x = 0; x < LON; x++) {
    const a = y * row + x, b = a + row;
    idx.push(a, b, a + 1, a + 1, b, b + 1);
  }
  return { posBuf: makeBuffer(new Float32Array(pos)), normBuf: makeBuffer(new Float32Array(nrm)),
           uvBuf: makeBuffer(new Float32Array(uv)), idxBuf: makeIndexBuffer(idx), count: idx.length };
}
function buildRingMesh(SEG) {
  const v = [], idx = [];
  for (let s = 0; s <= SEG; s++) { const a = s / SEG * TWO_PI, ca = Math.cos(a), sa = Math.sin(a); v.push(ca, sa, 0, ca, sa, 1); }
  for (let s = 0; s < SEG; s++) { const i = s * 2; idx.push(i, i + 1, i + 2, i + 1, i + 3, i + 2); }
  return { buf: makeBuffer(new Float32Array(v)), idxBuf: makeIndexBuffer(idx), count: idx.length };
}
const SPHERE_MESH = buildUVSphere(48, 96);
const RING_MESH = buildRingMesh(160);

const planetTex = {};
function loadTex(file) {
  if (planetTex[file]) return planetTex[file];
  const t = gl.createTexture();
  gl.bindTexture(gl.TEXTURE_2D, t);
  gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, 1, 1, 0, gl.RGBA, gl.UNSIGNED_BYTE, new Uint8Array([90, 92, 100, 255]));
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
  gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  const img = new Image();
  img.onload = () => {
    gl.bindTexture(gl.TEXTURE_2D, t);
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);   // equirectangular maps are top-down
    gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, img);
    gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, false);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);     // NPOT-safe: no mipmap
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
    gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
  };
  img.src = "img/bodies/" + file;
  return (planetTex[file] = t);
}

// axial obliquity to orbit (deg) — fixed-axis tilt about world X (no node data → approximation)
const PLANET_TILT = [0.034, 177.4, 23.44, 25.19, 3.13, 26.73, 97.77, 28.32];
// rings, radii in planet-radii: Saturn full system; the others a faint band
const RING_DEF = {
  4: { r0: 1.4, r1: 1.8, col: [0.85, 0.80, 0.70], a: 0.05, cas: 0 },     // Jupiter
  5: { r0: 1.18, r1: 2.27, col: [0.86, 0.80, 0.66], a: 0.6, cas: 0.5 },  // Saturn
  6: { r0: 1.6, r1: 2.0, col: [0.70, 0.85, 0.88], a: 0.08, cas: 0 },     // Uranus
  7: { r0: 1.6, r1: 2.0, col: [0.60, 0.70, 1.0], a: 0.07, cas: 0 },      // Neptune
};
const _planetProj = new Float32Array(16), _planetVP = new Float32Array(16);
const _sphereM = new Float32Array(16), _sphereMVP = new Float32Array(16);
const _ringM = new Float32Array(16), _ringMVP = new Float32Array(16);
function smoothstep(a, b, x) { x = clamp((x - a) / (b - a), 0, 1); return x * x * (3 - 2 * x); }
// model matrix: translate(p) · Rx(tilt) · scale(s), column-major
function tiltM(out, p, ct, st, s) {
  out[0] = s; out[1] = 0; out[2] = 0; out[3] = 0;
  out[4] = 0; out[5] = ct * s; out[6] = st * s; out[7] = 0;
  out[8] = 0; out[9] = -st * s; out[10] = ct * s; out[11] = 0;
  out[12] = p[0]; out[13] = p[1]; out[14] = p[2]; out[15] = 1;
  return out;
}
function bindSphere() {
  gl.bindBuffer(gl.ARRAY_BUFFER, SPHERE_MESH.posBuf); gl.enableVertexAttribArray(SP.aPos); gl.vertexAttribPointer(SP.aPos, 3, gl.FLOAT, false, 0, 0);
  gl.bindBuffer(gl.ARRAY_BUFFER, SPHERE_MESH.normBuf); gl.enableVertexAttribArray(SP.aNorm); gl.vertexAttribPointer(SP.aNorm, 3, gl.FLOAT, false, 0, 0);
  gl.bindBuffer(gl.ARRAY_BUFFER, SPHERE_MESH.uvBuf); gl.enableVertexAttribArray(SP.aUV); gl.vertexAttribPointer(SP.aUV, 2, gl.FLOAT, false, 0, 0);
  gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, SPHERE_MESH.idxBuf);
}
// the whole pipeline runs depth-off + additive; bracket the solid sphere/rings
// with depth-test + alpha blend, then restore exactly (or the next frame breaks).
function drawFocusedPlanet(fi, sphereA) {
  const ps = planetState[fi];
  const rAU = (PLANET_FACTS[fi].d * 0.5) / AU_KM;
  const th = PLANET_TILT[fi] * DEG, ct = Math.cos(th), st = Math.sin(th);
  // dedicated tight projection so the tiny sphere/rings get usable depth precision
  mat4Perspective(_planetProj, FOV, canvas.width / canvas.height, Math.max(state.cam.dist * 0.1, 1e-6), state.cam.dist * 4 + rAU * 8);
  mat4Mul(_planetVP, _planetProj, view);
  tiltM(_sphereM, ps.pos, ct, st, rAU);
  mat4Mul(_sphereMVP, _planetVP, _sphereM);
  gl.enable(gl.DEPTH_TEST); gl.depthFunc(gl.LEQUAL); gl.clear(gl.DEPTH_BUFFER_BIT);
  gl.blendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA);
  gl.useProgram(spProg);
  gl.uniformMatrix4fv(SP.uMVP, false, _sphereMVP);
  gl.uniformMatrix4fv(SP.uModel, false, _sphereM);
  gl.uniform1f(SP.uAmbient, 0.13);   // night side dim-but-visible so the full disc reads
  gl.uniform1f(SP.uAlpha, sphereA);
  gl.activeTexture(gl.TEXTURE0); gl.bindTexture(gl.TEXTURE_2D, loadTex(PLANET_IMG[fi][0])); gl.uniform1i(SP.uTex, 0);
  bindSphere();
  gl.drawElements(gl.TRIANGLES, SPHERE_MESH.count, gl.UNSIGNED_SHORT, 0);
  const rd = RING_DEF[fi];
  if (rd) {
    tiltM(_ringM, ps.pos, ct, st, 1);
    mat4Mul(_ringMVP, _planetVP, _ringM);
    gl.depthMask(false);                 // translucent rings: test against the sphere, don't self-sort
    gl.useProgram(rnProg);
    gl.uniformMatrix4fv(RN.uMVP, false, _ringMVP);
    gl.uniform1f(RN.uR0, rAU * rd.r0); gl.uniform1f(RN.uR1, rAU * rd.r1);
    gl.uniform3fv(RN.uColor, rd.col); gl.uniform1f(RN.uAlpha, rd.a * sphereA); gl.uniform1f(RN.uCassini, rd.cas);
    gl.bindBuffer(gl.ARRAY_BUFFER, RING_MESH.buf); gl.enableVertexAttribArray(RN.aRing); gl.vertexAttribPointer(RN.aRing, 3, gl.FLOAT, false, 0, 0);
    gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, RING_MESH.idxBuf);
    gl.drawElements(gl.TRIANGLES, RING_MESH.count, gl.UNSIGNED_SHORT, 0);
    gl.depthMask(true);
  }
  gl.disable(gl.DEPTH_TEST);
  gl.blendFunc(gl.ONE, gl.ONE);          // restore the global additive state
  gl.useProgram(ptProg); gl.uniformMatrix4fv(PT.uVP, false, vp);   // hand back to the point pipeline
}

function makeStars(count, banded) {
  const pos = new Float32Array(count * 3);
  const size = new Float32Array(count);
  const phase = new Float32Array(count);
  const R = 1300;
  // tilted plane for a "milky way" band
  const bn = [0.48, 0.2, 0.85];
  const bl = Math.hypot(bn[0], bn[1], bn[2]);
  bn[0] /= bl; bn[1] /= bl; bn[2] /= bl;
  for (let k = 0; k < count; k++) {
    let x, y, z;
    do {
      x = Math.random() * 2 - 1; y = Math.random() * 2 - 1; z = Math.random() * 2 - 1;
    } while (x * x + y * y + z * z > 1 || x * x + y * y + z * z < 1e-4);
    const l = Math.hypot(x, y, z);
    x /= l; y /= l; z /= l;
    if (banded) {
      // squash toward the band plane
      const d = x * bn[0] + y * bn[1] + z * bn[2];
      const f = d * (1 - Math.pow(Math.random(), 2.2) * 0.92);
      x -= bn[0] * (d - f); y -= bn[1] * (d - f); z -= bn[2] * (d - f);
      const l2 = Math.hypot(x, y, z);
      x /= l2; y /= l2; z /= l2;
    }
    pos[k * 3] = x * R; pos[k * 3 + 1] = y * R; pos[k * 3 + 2] = z * R;
    size[k] = (banded ? 1.1 : 1.6) * (0.5 + Math.pow(Math.random(), 3) * 2.6);
    phase[k] = 0.05 + Math.random();
  }
  return {
    count,
    posBuf: makeBuffer(pos),
    sizeBuf: makeBuffer(size),
    phaseBuf: makeBuffer(phase),
  };
}
const starsCool = makeStars(1900, false);
const starsBand = makeStars(2400, true);
const starsWarm = makeStars(450, false);

const sunBuf = makeBuffer(new Float32Array([0, 0, 0]));
const oneBuf = makeBuffer(new Float32Array([0]));        // phase=0
const sunSizeBuf = makeBuffer(new Float32Array([0.45]));

/* ---- Oort cloud — inferred representation, NOT data ----
   No Oort-cloud object has ever been observed; the shell is a statistical
   sketch at its true inferred scale: a flattened inner (Hills) cloud from
   ~2,000–20,000 au and an isotropic outer shell to ~100,000 au, density
   falling as ~r^-3.5. Static points (orbital periods are millions of
   years), seeded so every visitor sees the same cloud. */
const OORT = { count: 90000, visible: true };
(function buildOort() {
  const rnd = mulberry32(20260612);
  const N = OORT.count;
  const pos = new Float32Array(N * 3);
  const size = new Float32Array(N);
  // shell number-density ∝ r^-3.5 → inverse-CDF sample on r^-0.5
  const sampleR = (r0, r1, u) => {
    const a = 1 / Math.sqrt(r0), b = 1 / Math.sqrt(r1);
    const s = a - u * (a - b);
    return 1 / (s * s);
  };
  for (let k = 0; k < N; k++) {
    const inner = k < N * 0.45;               // Hills cloud share
    // inner cloud: centrally condensed (r^-3.5); outer shell: spread evenly
    // along r so the vast sphere stays legible at full zoom-out
    const r = inner ? sampleR(2000, 20000, rnd()) : 20000 + rnd() * 80000;
    const zScale = inner ? 0.38 : 1;          // inner cloud hugs the ecliptic
    let x, y, z, l2;
    do {
      x = rnd() * 2 - 1; y = rnd() * 2 - 1; z = rnd() * 2 - 1;
      l2 = x * x + y * y + z * z;
    } while (l2 > 1 || l2 < 1e-4);
    z *= zScale;
    const l = Math.hypot(x, y, z);
    pos[k * 3] = (x / l) * r; pos[k * 3 + 1] = (y / l) * r; pos[k * 3 + 2] = (z / l) * r;
    size[k] = 60 + rnd() * 160;               // representational, not physical
  }
  OORT.posBuf = makeBuffer(pos);
  OORT.sizeBuf = makeBuffer(size);
})();

const planetState = PLANETS.map((p) => ({
  def: p,
  pos: new Float64Array(3),
  sizeBuf: makeBuffer(new Float32Array([p.size])),
  posBuf: makeBuffer(new Float32Array(3), gl.DYNAMIC_DRAW),
  orbitBuf: null,
  orbitCount: 0,
  labelEl: null,
}));

function buildPlanetOrbits() {
  const jd = state.simJD;
  const SEG = 180;
  const tmp = new Float64Array(3);
  for (const ps of planetState) {
    const el = planetElements(ps.def, jd, { a: 0, e: 0, i: 0, om: 0, w: 0, M: 0 });
    const P = perifocalBasis(el.w, el.om, el.i, new Float64Array(6));
    const b = el.a * Math.sqrt(1 - el.e * el.e);
    const verts = new Float32Array(SEG * 3);
    for (let s = 0; s < SEG; s++) {
      ellipsePoint(el.a, el.e, b, P, (s / SEG) * TWO_PI, tmp);
      verts[s * 3] = tmp[0]; verts[s * 3 + 1] = tmp[1]; verts[s * 3 + 2] = tmp[2];
    }
    if (ps.orbitBuf) gl.deleteBuffer(ps.orbitBuf);
    ps.orbitBuf = makeBuffer(verts);
    ps.orbitCount = SEG;
  }
}
buildPlanetOrbits();

/* selected-orbit line (dynamic) — closed ellipse, or open hyperbola for ISOs */
let selOrbitBuf = gl.createBuffer();
let selOrbitCount = 0;
let selOrbitStrip = false;
function buildSelectedOrbit(rec) {
  const SEG = 256;
  const verts = new Float32Array(SEG * 3);
  const P = [rec.Px, rec.Py, rec.Pz, rec.Qx, rec.Qy, rec.Qz];
  if (rec.e > 1) {
    selOrbitStrip = true;
    const aAbs = Math.abs(rec.a);
    const Fmax = Math.min(5, Math.acosh((80 / aAbs + 1) / rec.e) || 4);
    for (let s = 0; s < SEG; s++) {
      const F = -Fmax + 2 * Fmax * (s / (SEG - 1));
      const xp = aAbs * (rec.e - Math.cosh(F));
      const yp = rec.b * Math.sinh(F);
      verts[s * 3] = xp * P[0] + yp * P[3];
      verts[s * 3 + 1] = xp * P[1] + yp * P[4];
      verts[s * 3 + 2] = xp * P[2] + yp * P[5];
    }
  } else {
    selOrbitStrip = false;
    const tmp = new Float64Array(3);
    for (let s = 0; s < SEG; s++) {
      ellipsePoint(rec.a, rec.e, rec.b, P, (s / SEG) * TWO_PI, tmp);
      verts[s * 3] = tmp[0]; verts[s * 3 + 1] = tmp[1]; verts[s * 3 + 2] = tmp[2];
    }
  }
  gl.bindBuffer(gl.ARRAY_BUFFER, selOrbitBuf);
  gl.bufferData(gl.ARRAY_BUFFER, verts, gl.DYNAMIC_DRAW);
  selOrbitCount = SEG;
}

/* selected moon's orbit, rebuilt as its parent planet moves */
let moonOrbitBuf = gl.createBuffer();
let moonOrbitCount = 0;
const _moVerts = new Float32Array(128 * 3);
function buildMoonOrbit(m) {
  if (!moonParentPos(m, _mpp)) { moonOrbitCount = 0; return; }
  const o = m * STRIDE, el = moons.el;
  const a = el[o], e = el[o + 1], b = el[o + 2];
  const SEG = 128;
  for (let s = 0; s < SEG; s++) {
    const E = (s / SEG) * TWO_PI;
    const xp = a * (Math.cos(E) - e);
    const yp = b * Math.sin(E);
    _moVerts[s * 3] = _mpp[0] + xp * el[o + 6] + yp * el[o + 9];
    _moVerts[s * 3 + 1] = _mpp[1] + xp * el[o + 7] + yp * el[o + 10];
    _moVerts[s * 3 + 2] = _mpp[2] + xp * el[o + 8] + yp * el[o + 11];
  }
  gl.bindBuffer(gl.ARRAY_BUFFER, moonOrbitBuf);
  gl.bufferData(gl.ARRAY_BUFFER, _moVerts, gl.DYNAMIC_DRAW);
  moonOrbitCount = SEG;
}

/* ---- draw helpers ---- */
function bindPointAttrs(posBuf, sizeBuf, phaseBuf, singlePhase) {
  gl.bindBuffer(gl.ARRAY_BUFFER, posBuf);
  gl.enableVertexAttribArray(PT.aPos);
  gl.vertexAttribPointer(PT.aPos, 3, gl.FLOAT, false, 0, 0);
  gl.bindBuffer(gl.ARRAY_BUFFER, sizeBuf);
  gl.enableVertexAttribArray(PT.aSize);
  gl.vertexAttribPointer(PT.aSize, 1, gl.FLOAT, false, 0, 0);
  if (phaseBuf) {
    gl.bindBuffer(gl.ARRAY_BUFFER, phaseBuf);
    gl.enableVertexAttribArray(PT.aPhase);
    gl.vertexAttribPointer(PT.aPhase, 1, gl.FLOAT, false, 0, 0);
  } else {
    gl.disableVertexAttribArray(PT.aPhase);
    gl.vertexAttrib1f(PT.aPhase, singlePhase || 0);
  }
}

/* ============================================================
   5. Position propagation (CPU, chunk-rotated for mobile)
   ============================================================ */
function updateGroupPositions(g, from, to, t) {
  const el = g.el, pos = g.pos;
  const IT = g.kepIters || 8;
  for (let k = from; k < to; k++) {
    const o = k * STRIDE;
    const a = el[o], e = el[o + 1], b = el[o + 2];
    let M = el[o + 3] + el[o + 4] * (t - el[o + 5]);
    M = M % TWO_PI;
    if (M > Math.PI) M -= TWO_PI; else if (M < -Math.PI) M += TWO_PI;
    let E = e < 0.8 ? M : Math.PI * (M < 0 ? -1 : 1);
    for (let j = 0; j < IT; j++) {
      const s = Math.sin(E), c = Math.cos(E);
      const d = (E - e * s - M) / (1 - e * c);
      E -= d;
      if (d < 1e-7 && d > -1e-7) break;
    }
    const xp = a * (Math.cos(E) - e);
    const yp = b * Math.sin(E);
    pos[k * 3] = xp * el[o + 6] + yp * el[o + 9];
    pos[k * 3 + 1] = xp * el[o + 7] + yp * el[o + 10];
    pos[k * 3 + 2] = xp * el[o + 8] + yp * el[o + 11];
  }
}

// Hyperbolic version for interstellar objects (e > 1, a < 0). Solves the
// hyperbolic Kepler equation  M = e·sinh(F) − F  by Newton iteration.
function updateHyperbolic(g, from, to, t) {
  const el = g.el, pos = g.pos;
  for (let k = from; k < to; k++) {
    const o = k * STRIDE;
    const aAbs = -el[o], e = el[o + 1], b = el[o + 2];
    const M = el[o + 3] + el[o + 4] * (t - el[o + 5]);
    let F = Math.asinh(M / e) || M;
    for (let j = 0; j < 24; j++) {
      const d = (e * Math.sinh(F) - F - M) / (e * Math.cosh(F) - 1);
      F -= d;
      if (d < 1e-9 && d > -1e-9) break;
    }
    const xp = aAbs * (e - Math.cosh(F));   // = q at F = 0
    const yp = b * Math.sinh(F);
    pos[k * 3] = xp * el[o + 6] + yp * el[o + 9];
    pos[k * 3 + 1] = xp * el[o + 7] + yp * el[o + 10];
    pos[k * 3 + 2] = xp * el[o + 8] + yp * el[o + 11];
  }
}

/* World position of the camera's focus body (Sun when unfocused). */
const _sunOrigin = new Float64Array(3);
const _focusScratch = new Float64Array(3);
function focusPosition() {
  const f = state.focus;
  if (!f) return _sunOrigin;
  if (f.planet != null) return planetState[f.planet].pos;
  const [gi, k] = f.small;
  const p = groups[gi].pos;
  _focusScratch[0] = p[k * 3]; _focusScratch[1] = p[k * 3 + 1]; _focusScratch[2] = p[k * 3 + 2];
  return _focusScratch;
}

/* Moon position of moon m's parent body, written into out. */
function moonParentPos(m, out) {
  const pi = moons.parentIdx[m];
  if (pi >= 0) {
    const p = planetState[pi].pos;
    out[0] = p[0]; out[1] = p[1]; out[2] = p[2];
    return true;
  }
  const ref = moons.parentSmall[m];
  if (!ref) return false;
  const p = groups[ref[0]].pos;
  out[0] = p[ref[1] * 3]; out[1] = p[ref[1] * 3 + 1]; out[2] = p[ref[1] * 3 + 2];
  return true;
}

/* Kepler around the parent planet, then offset by the parent's position. */
const _mpp = new Float64Array(3);
function updateMoonPositions(jd) {
  if (!moons.count) return;
  const el = moons.el, pos = moons.pos;
  const t = jd - J2000;
  for (let k = 0; k < moons.count; k++) {
    if (!moonParentPos(k, _mpp)) { pos[k * 3] = 1e9; continue; }
    const o = k * STRIDE;
    const a = el[o], e = el[o + 1], b = el[o + 2];
    let M = el[o + 3] + el[o + 4] * (t - el[o + 5]);
    M = M % TWO_PI;
    if (M > Math.PI) M -= TWO_PI; else if (M < -Math.PI) M += TWO_PI;
    let E = e < 0.8 ? M : Math.PI * (M < 0 ? -1 : 1);
    for (let j = 0; j < 10; j++) {
      const s = Math.sin(E), cE = Math.cos(E);
      const d = (E - e * s - M) / (1 - e * cE);
      E -= d;
      if (d < 1e-8 && d > -1e-8) break;
    }
    const xp = a * (Math.cos(E) - e);
    const yp = b * Math.sin(E);
    pos[k * 3] = _mpp[0] + xp * el[o + 6] + yp * el[o + 9];
    pos[k * 3 + 1] = _mpp[1] + xp * el[o + 7] + yp * el[o + 10];
    pos[k * 3 + 2] = _mpp[2] + xp * el[o + 8] + yp * el[o + 11];
  }
  moons.dirty = true;
}

function propagate() {
  const t = state.simJD - J2000;
  const t0 = performance.now();
  const slices = state.needFullUpdate ? 1 : state.updateSlices;
  const cursor = state.sliceCursor;
  for (const g of groups) {
    if (!g.count) continue;
    const upd = g.hyperbolic ? updateHyperbolic : updateGroupPositions;
    if (slices === 1 || g.hyperbolic) {
      upd(g, 0, g.count, t);
    } else {
      const span = Math.ceil(g.count / slices);
      const from = Math.min(cursor * span, g.count);
      const to = Math.min(from + span, g.count);
      upd(g, from, to, t);
    }
    g.dirty = true;
  }
  state.sliceCursor = (cursor + 1) % Math.max(state.updateSlices, 1);
  state.needFullUpdate = false;
  // adapt: if a full pass is slow, rotate through the population in slices
  const took = performance.now() - t0;
  if (slices === 1 && took > 7) state.updateSlices = Math.min(4, Math.ceil(took / 5));
}

/* ============================================================
   6. Render loop
   ============================================================ */
const eye = new Float64Array(3);
let lastFrame = performance.now();
let simAtLastPropagate = NaN;
const _propTmp = new Float64Array(3), _planetUp = new Float32Array(3);   // per-frame planet scratch (no realloc)

// True-scale: real physical radii (au) fed into the point-size attribute. The
// shader already projects size·focalPx/distance, so real radii make most bodies
// shrink to sub-pixel and vanish — the honest emptiness of space. Built lazily.
const PLANET_RADIUS_KM = [2439.7, 6051.8, 6371.0, 3389.5, 69911, 58232, 25362, 24622];
const SUN_RADIUS_KM = 695700;
let sunTrueSizeBuf = null, trueSizesBuilt = false;
function ensureTrueSizes() {
  for (let i = 0; i < planetState.length; i++)
    planetState[i].trueSizeBuf = makeBuffer(new Float32Array([PLANET_RADIUS_KM[i] / AU_KM]));
  sunTrueSizeBuf = makeBuffer(new Float32Array([SUN_RADIUS_KM / AU_KM]));
  for (const g of groups) {
    if (!g.count) continue;
    const a = new Float32Array(g.count);
    for (let k = 0; k < g.count; k++) {
      const m = g.meta[k];
      const alb = isFinite(m.albedo) && m.albedo > 0 ? m.albedo : 0.14;
      const dkm = isFinite(m.diam) ? m.diam
        : isFinite(m.H) ? 1329 / Math.sqrt(alb) * Math.pow(10, -m.H / 5) : 1;
      a[k] = (dkm * 0.5) / AU_KM;
    }
    g.trueSizeBuf = makeBuffer(a);
  }
  if (moons.count) {
    const a = new Float32Array(moons.count);
    for (let k = 0; k < moons.count; k++) a[k] = (moons.meta[k].radius || 1) / AU_KM;
    moons.trueSizeBuf = makeBuffer(a);
  }
  trueSizesBuilt = true;
}

// Comet tails — the Sun is at the origin, so the anti-sunward direction is just
// pos/|pos|. Near-perihelion comets (r < TAIL_R0) get an additive streak whose
// length grows ~1/r; one batched LINES draw covers both comet populations.
const TAIL_R0 = 3.5, TAIL_LEN = 0.25, TAIL_MAX = 0.35;   // au — short accents, not spears
let tailVerts = new Float32Array(0), tailBuf = null;
function drawCometTails() {
  let cap = 0;
  for (const gi of [GI.COM, GI.LPC]) { const g = groups[gi]; if (g.count && g.visible) cap += g.count; }
  if (!cap) return;
  if (tailVerts.length < cap * 6) tailVerts = new Float32Array(cap * 6);
  let o = 0;
  for (const gi of [GI.COM, GI.LPC]) {
    const g = groups[gi]; if (!g.count || !g.visible) continue;
    for (let k = 0; k < g.count; k++) {
      const x = g.pos[k * 3], y = g.pos[k * 3 + 1], z = g.pos[k * 3 + 2];
      const r = Math.hypot(x, y, z);
      if (!(r > 0 && r < TAIL_R0)) continue;
      const inv = 1 / r;
      let L = TAIL_LEN * (inv - 1 / TAIL_R0); if (L > TAIL_MAX) L = TAIL_MAX;
      tailVerts[o++] = x; tailVerts[o++] = y; tailVerts[o++] = z;                       // nucleus
      tailVerts[o++] = x + x * inv * L; tailVerts[o++] = y + y * inv * L; tailVerts[o++] = z + z * inv * L; // tip (anti-sunward)
    }
  }
  if (!o) return;
  gl.useProgram(lnProg);
  gl.uniformMatrix4fv(LN.uVP, false, vp);
  gl.enableVertexAttribArray(LN.aPos);
  if (!tailBuf) tailBuf = gl.createBuffer();
  gl.bindBuffer(gl.ARRAY_BUFFER, tailBuf);
  gl.bufferData(gl.ARRAY_BUFFER, tailVerts.subarray(0, o), gl.DYNAMIC_DRAW);
  gl.vertexAttribPointer(LN.aPos, 3, gl.FLOAT, false, 0, 0);
  const c = groups[GI.COM].color;
  gl.uniform3f(LN.uColor, c[0], c[1], c[2]);
  gl.uniform1f(LN.uAlpha, 0.4);
  gl.drawArrays(gl.LINES, 0, o / 3);
}

function render(now) {
  const dt = Math.min((now - lastFrame) / 1000, 0.25);
  lastFrame = now;

  // time machine
  if (state.playing) state.simJD += state.dir * state.dps * dt;
  updateBrandOrrery(dt);
  // deep-link: consume lazy sat restores, then mirror the live view into the URL (debounced)
  if (pendingSatSel != null && sats.loaded) { selectSatellite(pendingSatSel); pendingSatSel = null; }
  if (pendingSatLayers && sats.loaded) { sats.visible = pendingSatLayers.wantSat; sats.debrisVisible = pendingSatLayers.wantDeb; pendingSatLayers = null; renderLegend(); }
  if (loadedOnce && !isApplying && now - lastHashCheck > 200) {   // ~5 Hz; the write is already debounced
    lastHashCheck = now;
    const h = encodeState(); if (h !== lastHash) { lastHash = h; scheduleHashWrite(h); }
  }

  // propagate asteroid + planet + moon positions when time moved — BEFORE the
  // camera reads focusPosition(), or a focused planet is aimed at one frame late
  if (state.simJD !== simAtLastPropagate || state.needFullUpdate) {
    propagate();
    simAtLastPropagate = state.simJD;
    for (const ps of planetState) {
      planetPosition(ps.def, state.simJD, _propTmp);
      ps.pos.set(_propTmp);
      _planetUp.set(_propTmp);   // narrow Float64 → Float32 for the attribute upload
      gl.bindBuffer(gl.ARRAY_BUFFER, ps.posBuf);
      gl.bufferData(gl.ARRAY_BUFFER, _planetUp, gl.DYNAMIC_DRAW);
    }
    updateMoonPositions(state.simJD);
    if (state.selMoon != null) buildMoonOrbit(state.selMoon);
    for (const c of craft) {
      craftPositionAt(c, state.simJD);
      gl.bindBuffer(gl.ARRAY_BUFFER, c.posBuf);
      gl.bufferData(gl.ARRAY_BUFFER, c.pos, gl.DYNAMIC_DRAW);
    }
  }

  // satellites: live, only at Earth, while the TLE snapshot is fresh
  const earthFocused = state.focus && state.focus.planet === EARTH_IDX;
  if (sats.loaded && (sats.visible || sats.debrisVisible) && earthFocused && satEpochFade() > 0.01) {
    updateSatellites();   // live: advance at real-time every frame, not on the sim clock
  }

  // camera easing (orientation, distance, and focus target)
  const c = state.cam;
  // hand control back at normal speed once the launch dolly has settled
  if (state.camEaseRate !== 9 && Math.abs(c.dist - c.tDist) < Math.max(c.tDist * 0.04, 1e-5)) {
    state.camEaseRate = 9;
  }
  const ease = 1 - Math.exp(-dt * state.camEaseRate);
  c.yaw += (c.tYaw - c.yaw) * ease;
  c.pitch += (c.tPitch - c.pitch) * ease;
  // during the tour, glide distance in LOG space so huge zoom changes sweep
  // uniformly through every scale (a continuous zoom) instead of snapping the near part
  if (tourActive && c.dist > 0 && c.tDist > 0) {
    c.dist = Math.exp(Math.log(c.dist) + (Math.log(c.tDist) - Math.log(c.dist)) * ease);
  } else {
    c.dist += (c.tDist - c.dist) * ease;
  }
  // blend into the focus body once, then track it rigidly — easing the target
  // per-frame would trail a moving planet by more than a zoomed-in view width
  const fp = focusPosition();
  c.fBlend = Math.min(1, c.fBlend + dt * 2.2);
  const bl = c.fBlend * c.fBlend * (3 - 2 * c.fBlend);
  c.target[0] = c.fFrom[0] + (fp[0] - c.fFrom[0]) * bl;
  c.target[1] = c.fFrom[1] + (fp[1] - c.fFrom[1]) * bl;
  c.target[2] = c.fFrom[2] + (fp[2] - c.fFrom[2]) * bl;

  eye[0] = c.target[0] + c.dist * Math.cos(c.pitch) * Math.cos(c.yaw);
  eye[1] = c.target[1] + c.dist * Math.cos(c.pitch) * Math.sin(c.yaw);
  eye[2] = c.target[2] + c.dist * Math.sin(c.pitch);
  // near plane follows zoom so a focused moon system isn't clipped away;
  // far plane covers the Oort shell and long-period comet orbits
  mat4Perspective(proj, FOV, canvas.width / canvas.height, Math.min(0.01, c.dist * 0.2), 300000);
  mat4LookAt(view, eye, c.target, [0, 0, 1]);
  mat4Mul(vp, proj, view);
  // star dome: same orientation, no translation — a sky at any zoom level
  viewSky.set(view);
  viewSky[12] = viewSky[13] = viewSky[14] = 0;
  mat4Mul(vpSky, proj, viewSky);

  gl.clear(gl.COLOR_BUFFER_BIT);
  const pixScale = canvas.height / (2 * Math.tan(FOV / 2));
  const timeS = now / 1000;

  /* ---- orbit lines ---- */
  gl.useProgram(lnProg);
  gl.uniformMatrix4fv(LN.uVP, false, vp);
  gl.enableVertexAttribArray(LN.aPos);
  if (state.showPlanets) {
    for (const ps of planetState) {
      gl.bindBuffer(gl.ARRAY_BUFFER, ps.orbitBuf);
      gl.vertexAttribPointer(LN.aPos, 3, gl.FLOAT, false, 0, 0);
      gl.uniform3f(LN.uColor, ps.def.color[0], ps.def.color[1], ps.def.color[2]);
      gl.uniform1f(LN.uAlpha, 0.14);
      gl.drawArrays(gl.LINE_LOOP, 0, ps.orbitCount);
    }
  }
  if (state.selected && selOrbitCount) {
    gl.bindBuffer(gl.ARRAY_BUFFER, selOrbitBuf);
    gl.vertexAttribPointer(LN.aPos, 3, gl.FLOAT, false, 0, 0);
    const gc = groups[state.selected.group].color;
    gl.uniform3f(LN.uColor, gc[0], gc[1], gc[2]);
    gl.uniform1f(LN.uAlpha, 0.55);
    gl.drawArrays(selOrbitStrip ? gl.LINE_STRIP : gl.LINE_LOOP, 0, selOrbitCount);
  }
  // selected moon's orbit around its parent
  if (state.selMoon != null && moonOrbitCount) {
    gl.bindBuffer(gl.ARRAY_BUFFER, moonOrbitBuf);
    gl.vertexAttribPointer(LN.aPos, 3, gl.FLOAT, false, 0, 0);
    gl.uniform3f(LN.uColor, 0.75, 0.82, 0.95);
    gl.uniform1f(LN.uAlpha, 0.5);
    gl.drawArrays(gl.LINE_LOOP, 0, moonOrbitCount);
  }
  // spacecraft flight paths
  if (craftVisible) {
    for (let i = 0; i < craft.length; i++) {
      const c = craft[i];
      gl.bindBuffer(gl.ARRAY_BUFFER, c.trajBuf);
      gl.vertexAttribPointer(LN.aPos, 3, gl.FLOAT, false, 0, 0);
      const sel = state.selCraft === i;
      gl.uniform3f(LN.uColor, 0.56, 0.94, 0.69);
      gl.uniform1f(LN.uAlpha, sel ? 0.7 : 0.22);
      gl.drawArrays(gl.LINE_STRIP, 0, c.n);
    }
  }

  /* ---- points ---- */
  gl.useProgram(ptProg);
  gl.uniformMatrix4fv(PT.uVP, false, vpSky);   // stars ride the rotation-only dome
  gl.uniform1f(PT.uPixScale, pixScale);
  gl.uniform1f(PT.uTime, timeS);
  gl.uniform1f(PT.uCull, 0);          // default: never cull (true-scale blocks set it per-draw)

  // stars (sizes are in px, not world units: fake it with min=max clamp window)
  gl.uniform1f(PT.uMinPx, 0.6 * dpr);
  gl.uniform1f(PT.uMaxPx, 3.4 * dpr);
  gl.uniform1f(PT.uAlpha, 0.9);
  gl.uniform3f(PT.uColor, 0.78, 0.85, 1.0);
  bindPointAttrs(starsCool.posBuf, starsCool.sizeBuf, starsCool.phaseBuf);
  gl.drawArrays(gl.POINTS, 0, starsCool.count);
  gl.uniform3f(PT.uColor, 0.62, 0.74, 0.98);
  gl.uniform1f(PT.uAlpha, 0.55);
  bindPointAttrs(starsBand.posBuf, starsBand.sizeBuf, starsBand.phaseBuf);
  gl.drawArrays(gl.POINTS, 0, starsBand.count);
  gl.uniform3f(PT.uColor, 1.0, 0.86, 0.66);
  gl.uniform1f(PT.uAlpha, 0.8);
  bindPointAttrs(starsWarm.posBuf, starsWarm.sizeBuf, starsWarm.phaseBuf);
  gl.drawArrays(gl.POINTS, 0, starsWarm.count);
  gl.uniformMatrix4fv(PT.uVP, false, vp);      // back to world space

  // true-scale: real angular sizes — tiny bodies vanish (built lazily on first use)
  const ts = state.trueScale;
  if (ts && !trueSizesBuilt) ensureTrueSizes();

  // asteroids — recede politely while the camera is parked at a planet, so
  // the focused world isn't buried in out-of-focus foreground confetti
  const focusFade = state.focus ? 0.22 : 1;
  gl.uniform1f(PT.uMinPx, ts ? 0 : 1.25 * dpr);
  gl.uniform1f(PT.uMaxPx, ts ? 40 * dpr : (state.focus ? 3.5 : 9) * dpr);
  gl.uniform1f(PT.uCull, ts ? 1 : 0);   // true-scale: vanish unless zoomed onto a real body
  for (const g of groups) {
    if (!g.count || !g.visible) continue;
    if (g.dirty) {
      gl.bindBuffer(gl.ARRAY_BUFFER, g.posBuf);
      gl.bufferData(gl.ARRAY_BUFFER, g.pos, gl.DYNAMIC_DRAW);
      g.dirty = false;
    }
    gl.uniform3f(PT.uColor, g.color[0], g.color[1], g.color[2]);
    gl.uniform1f(PT.uAlpha, 0.85 * focusFade);
    bindPointAttrs(g.posBuf, ts ? g.trueSizeBuf : g.sizeBuf, null, 0);
    gl.drawArrays(gl.POINTS, 0, g.count);
  }

  // Oort cloud representation — invisible inside the planetary system,
  // resolves as the camera makes the journey out
  if (OORT.visible) {
    const oortA = Math.min((state.cam.dist - 150) / 2500, 0.5);
    if (oortA > 0.004) {
      gl.uniform1f(PT.uMinPx, 0.7 * dpr);
      gl.uniform1f(PT.uMaxPx, 3 * dpr);
      gl.uniform3f(PT.uColor, 0.62, 0.7, 0.85);
      gl.uniform1f(PT.uAlpha, oortA);
      bindPointAttrs(OORT.posBuf, OORT.sizeBuf, null, 0);
      gl.drawArrays(gl.POINTS, 0, OORT.count);
    }
  }

  // PHA highlight overlay (orange, on top of the normal points)
  if (state.showPHA && phaList.length) {
    if (phaPos.length !== phaList.length * 3) phaPos = new Float32Array(phaList.length * 3);
    for (let i = 0; i < phaList.length; i++) {
      const g = groups[phaList[i].gi], k = phaList[i].k;
      phaPos[i * 3] = g.pos[k * 3]; phaPos[i * 3 + 1] = g.pos[k * 3 + 1]; phaPos[i * 3 + 2] = g.pos[k * 3 + 2];
    }
    if (!phaPosBuf) phaPosBuf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, phaPosBuf);
    gl.bufferData(gl.ARRAY_BUFFER, phaPos, gl.DYNAMIC_DRAW);
    gl.enableVertexAttribArray(PT.aPos);
    gl.vertexAttribPointer(PT.aPos, 3, gl.FLOAT, false, 0, 0);
    gl.disableVertexAttribArray(PT.aSize); gl.vertexAttrib1f(PT.aSize, 0.05);
    gl.disableVertexAttribArray(PT.aPhase); gl.vertexAttrib1f(PT.aPhase, 0);
    gl.uniform1f(PT.uMinPx, 3.2 * dpr);
    gl.uniform1f(PT.uMaxPx, 13 * dpr);
    gl.uniform3f(PT.uColor, 1.0, 0.5, 0.16);
    gl.uniform1f(PT.uAlpha, 0.92);
    gl.drawArrays(gl.POINTS, 0, phaList.length);
  }

  // moons
  if (moons.count && moons.visible) {
    if (moons.dirty) {
      gl.bindBuffer(gl.ARRAY_BUFFER, moons.posBuf);
      gl.bufferData(gl.ARRAY_BUFFER, moons.pos, gl.DYNAMIC_DRAW);
      moons.dirty = false;
    }
    gl.uniform1f(PT.uMinPx, ts ? 0 : 1.4 * dpr);
    gl.uniform1f(PT.uMaxPx, ts ? 40 * dpr : 10 * dpr);
    gl.uniform1f(PT.uCull, ts ? 1 : 0);
    gl.uniform3f(PT.uColor, 0.78, 0.84, 0.94);
    gl.uniform1f(PT.uAlpha, 0.95);
    bindPointAttrs(moons.posBuf, ts ? moons.trueSizeBuf : moons.sizeBuf, null, 0);
    gl.drawArrays(gl.POINTS, 0, moons.count);
  }

  // focused planet → a lit sphere fades in over the sprite as you approach
  const fpi = (state.focus && state.focus.planet != null) ? state.focus.planet : -1;
  let sphereA = 0;
  if (fpi >= 0) {
    const rAU = (PLANET_FACTS[fpi].d * 0.5) / AU_KM;
    sphereA = smoothstep(2, 6, pixScale * rAU / Math.max(state.cam.dist, 1e-6));
  }

  // planets — a soft halo pass first so the eight anchors read above the
  // 50k-point cloud at overview zoom, then the bright core sprite
  if (state.showPlanets) {
    gl.uniform1f(PT.uCull, 0);
    if (!ts) {
      gl.uniform1f(PT.uMinPx, 7 * dpr);
      gl.uniform1f(PT.uMaxPx, 40 * dpr);
      for (let i = 0; i < planetState.length; i++) {
        const ps = planetState[i];
        gl.uniform3f(PT.uColor, ps.def.color[0], ps.def.color[1], ps.def.color[2]);
        gl.uniform1f(PT.uAlpha, 0.22 * (i === fpi ? 1 - sphereA : 1.0));
        bindPointAttrs(ps.posBuf, ps.sizeBuf, null, 0);
        gl.drawArrays(gl.POINTS, 0, 1);
      }
    }
    gl.uniform1f(PT.uMinPx, ts ? 1.5 * dpr : 3.6 * dpr);   // keep planets as faint points in true-scale
    gl.uniform1f(PT.uMaxPx, ts ? 200 * dpr : 26 * dpr);
    for (let i = 0; i < planetState.length; i++) {
      const ps = planetState[i];
      gl.uniform3f(PT.uColor, ps.def.color[0], ps.def.color[1], ps.def.color[2]);
      gl.uniform1f(PT.uAlpha, i === fpi ? 1 - sphereA : 1.0);   // fade the sprite under the lit sphere
      bindPointAttrs(ps.posBuf, ts ? ps.trueSizeBuf : ps.sizeBuf, null, 0);
      gl.drawArrays(gl.POINTS, 0, 1);
    }
  }
  if (fpi >= 0 && sphereA > 0.004) drawFocusedPlanet(fpi, sphereA);

  // satellites + debris (Earth's artificial moons) — only in the Earth-focus
  // view, drawn AFTER the globe so the swarm reads against it; when the lit
  // sphere is up we reuse its depth buffer (and its projection, so depths
  // compare correctly) to occlude satellites passing behind the planet.
  if (sats.loaded && (sats.visible || sats.debrisVisible) && earthFocused) {
    const sf = satEpochFade();
    if (sf > 0.01) {
      if (sats.dirty) {
        gl.bindBuffer(gl.ARRAY_BUFFER, sats.posBuf);
        gl.bufferData(gl.ARRAY_BUFFER, sats.pos, gl.DYNAMIC_DRAW);
        sats.dirty = false;
      }
      const occlude = fpi === EARTH_IDX && sphereA > 0.004;
      if (occlude) {
        gl.enable(gl.DEPTH_TEST); gl.depthFunc(gl.LEQUAL); gl.depthMask(false);
        gl.uniformMatrix4fv(PT.uVP, false, _planetVP);
      }
      gl.uniform1f(PT.uMinPx, 0.9 * dpr);
      gl.uniform1f(PT.uMaxPx, 3.5 * dpr);
      bindPointAttrs(sats.posBuf, sats.sizeBuf, null, 0);
      const pc = sats.payloadCount;
      if (sats.visible && pc > 0) {
        gl.uniform3f(PT.uColor, 0.85, 0.95, 1.0);
        gl.uniform1f(PT.uAlpha, 0.5 * sf);   // dense LEO must not white out the globe
        gl.drawArrays(gl.POINTS, 0, pc);
      }
      if (sats.debrisVisible && sats.count > pc) {
        gl.uniform3f(PT.uColor, 1.0, 0.47, 0.29);   // debris — orange
        gl.uniform1f(PT.uAlpha, 0.55 * sf);
        gl.drawArrays(gl.POINTS, pc, sats.count - pc);
      }
      if (occlude) {
        gl.disable(gl.DEPTH_TEST); gl.depthMask(true);
        gl.uniformMatrix4fv(PT.uVP, false, vp);
      }
    }
  }

  // sun core
  if (state.showSun) {
    gl.uniform1f(PT.uMinPx, ts ? 3 * dpr : 10 * dpr);      // Sun stays a visible anchor
    gl.uniform1f(PT.uMaxPx, ts ? 1024 * dpr : 58 * dpr);
    gl.uniform1f(PT.uCull, 0);
    gl.uniform3f(PT.uColor, 1.0, 0.93, 0.78);
    gl.uniform1f(PT.uAlpha, 1.0);
    bindPointAttrs(sunBuf, ts ? sunTrueSizeBuf : sunSizeBuf, null, 0);
    gl.drawArrays(gl.POINTS, 0, 1);
  }

  // spacecraft markers
  if (craftVisible) {
    gl.uniform1f(PT.uMinPx, 3.5 * dpr);
    gl.uniform1f(PT.uMaxPx, 10 * dpr);
    gl.uniform3f(PT.uColor, 0.56, 0.94, 0.69);
    gl.uniform1f(PT.uAlpha, 1.0);
    for (const c of craft) {
      if (!c.inRange) continue;
      bindPointAttrs(c.posBuf, c.sizeBuf, null, 0);
      gl.drawArrays(gl.POINTS, 0, 1);
    }
  }

  // comet tails (anti-sunward) — drawn last so they layer over everything
  // (additive); parked at a planet they recede with the rest of the field
  if (!state.focus && (groups[GI.COM].visible || groups[GI.LPC].visible)) drawCometTails();

  updateOverlays(pixScale);
  updateHUD();
  requestAnimationFrame(render);
}

/* project world → CSS px; returns null when behind the camera */
const _pv = new Float64Array(4);
function project(x, y, z) {
  _pv[0] = vp[0] * x + vp[4] * y + vp[8] * z + vp[12];
  _pv[1] = vp[1] * x + vp[5] * y + vp[9] * z + vp[13];
  _pv[3] = vp[3] * x + vp[7] * y + vp[11] * z + vp[15];
  if (_pv[3] <= 1e-6) return null;
  return [
    (_pv[0] / _pv[3] * 0.5 + 0.5) * cssW,
    (1 - (_pv[1] / _pv[3] * 0.5 + 0.5)) * cssH,
    _pv[3],
  ];
}

/* ---- HTML overlays: sun halo, planet labels, selection marker ---- */
let selMarkerEl = null;
let selMarkerLastRect = null;
// Region labels: anchored at world-space points, fading with camera distance.
// "→ Alpha Centauri" sits along the star's TRUE J2000 direction (RA 14h39.6m,
// Dec −60.8° rotated into the ecliptic frame) — you have to face it to see it.
const REGION_LABELS = [
  { text: "Asteroid Belt", cls: "belt-label", pos: [1.95, -1.95, 0.1],
    in0: 3.5, in1: 6, out0: 45, out1: 110, show: () => groups[GI.MBA].visible, onClick: () => selectRegion("asteroid") },
  { text: "Kuiper Belt", cls: "kuiper-label", pos: [31, -31, 2],
    in0: 15, in1: 28, out0: 250, out1: 800, show: () => groups[GI.TNO].visible, onClick: () => selectRegion("kuiper") },
  { text: "Oort cloud · inferred", cls: "oort-label", pos: [26000, 8000, 9000],
    in0: 1500, in1: 4000, show: () => OORT.visible, onClick: () => selectRegion("oort") },
  { text: "interstellar space", cls: "oort-label", pos: [110000, 60000, 30000],
    in0: 55000, in1: 105000, show: () => true, onClick: () => selectRegion("interstellar") },
  { text: "→ Alpha Centauri · 276,000 au", cls: "alpha-label", pos: [-44100, -74800, -79900],
    in0: 55000, in1: 105000, show: () => true, onClick: () => selectRegion("alphacen") },
];
const dwarfList = [];   // { gi, k, name, el } for persistent dwarf-planet labels
const moonLabelList = []; // { k, name, el } for major moons (radius ≥ 200 km)
let sunLabelEl = null;
function buildDwarfList() {
  dwarfList.length = 0;
  for (let gi = 0; gi < groups.length; gi++) {
    const g = groups[gi];
    for (let k = 0; k < g.count; k++) {
      if (g.meta[k].dwarf) {
        dwarfList.push({ gi, k, name: g.meta[k].dwarf, el: null });
        if (g.meta[k].dwarf === "Pluto") plutoRef = [gi, k];
      }
    }
  }
}
function buildPhaList() {
  phaList.length = 0;
  for (let gi = 0; gi < groups.length; gi++) {
    const g = groups[gi];
    for (let k = 0; k < g.count; k++) if (g.meta[k].pha) phaList.push({ gi, k });
  }
}
let phaPos = new Float32Array(0), phaPosBuf = null;
/* screen-space label collision: first placed wins; later labels that overlap
   an occupied rect stay hidden this frame instead of stacking */
const placedRects = [];
function rectFree(x, y, text) {
  const w = text.length * 6.8 + 10;
  const r = [x - w / 2, y - 26, x + w / 2, y - 4];
  for (const q of placedRects) {
    if (r[0] < q[2] && r[2] > q[0] && r[1] < q[3] && r[3] > q[1]) return false;
  }
  placedRects.push(r);
  return true;
}

function updateOverlays(pixScale) {
  placedRects.length = 0;
  // the selection marker always wins: reserve its footprint (from last frame)
  if (selMarkerLastRect) placedRects.push(selMarkerLastRect);
  // sun halo (fades out once the Sun is a sub-pixel speck)
  const sp = state.showSun ? project(0, 0, 0) : null;
  if (sp) {
    const px = (pixScale / dpr) * 0.55 / Math.max(sp[2], 0.05);
    const sc = Math.min(Math.max(px / 160, 0.22), 2.6);
    sunHalo.style.opacity = String(Math.min(px / 24, 1));
    sunHalo.style.transform = `translate(${sp[0] - cssW / 2}px, ${sp[1] - cssH / 2}px) scale(${sc})`;
  } else {
    sunHalo.style.opacity = "0";
  }
  // Sun label (persistent, like the planets)
  if (!sunLabelEl) {
    sunLabelEl = document.createElement("div");
    sunLabelEl.className = "pl-label";
    sunLabelEl.textContent = "Sun";
    labelsEl.appendChild(sunLabelEl);
  }
  // fade out as the camera pulls past the planets — otherwise it just sits as an
  // unreadable smudge on the collapsed inner-system blob (like the planet labels do)
  const sunFade = clamp(1 - (state.cam.dist - 8) / 12, 0, 1);   // full ≤8 au, gone by ~20 au
  const sshow = !state.selSun && sp && sunFade > 0.02 && sp[0] > -40 && sp[0] < cssW + 40 && sp[1] > -20 && sp[1] < cssH + 20 &&
    rectFree(sp[0], sp[1], "Sun");
  sunLabelEl.style.opacity = sshow ? sunFade.toFixed(2) : "0";
  if (sshow) sunLabelEl.style.transform = `translate(${sp[0]}px, ${sp[1]}px) translate(-50%,-150%)`;
  // planet labels (the selected planet's ambient label yields to the marker)
  for (let i = 0; i < planetState.length; i++) {
    const ps = planetState[i];
    if (!ps.labelEl) {
      ps.labelEl = document.createElement("div");
      ps.labelEl.className = "pl-label";
      ps.labelEl.textContent = ps.def.name;
      labelsEl.appendChild(ps.labelEl);
    }
    const p = project(ps.pos[0], ps.pos[1], ps.pos[2]);
    if (!p) { ps.labelEl.style.opacity = "0"; continue; }
    const px = (pixScale / dpr) * ps.def.size / p[2];
    // parked at a planet, background planet labels stand down
    const bgMuted = state.focus && !(state.focus.planet === i);
    const show = state.showPlanets && state.selPlanet !== i && !bgMuted && px > 1.6 &&
      p[0] > -40 && p[0] < cssW + 40 && p[1] > -20 && p[1] < cssH + 20 &&
      rectFree(p[0], p[1], ps.def.name);
    ps.labelEl.style.opacity = show ? "1" : "0";
    if (show) ps.labelEl.style.transform = `translate(${p[0]}px, ${p[1]}px) translate(-50%,-150%)`;
  }
  // dwarf-planet labels — fade out as the camera leaves the planetary
  // region (they'd pile up at the center of the Oort shot), back on return
  const dwarfFade = clamp(1 - (state.cam.dist - 250) / 550, 0, 1);
  for (const d of dwarfList) {
    const g = groups[d.gi];
    if (!d.el) {
      d.el = document.createElement("div");
      d.el.className = "pl-label dwarf-label";
      d.el.textContent = d.name;
      labelsEl.appendChild(d.el);
    }
    const isSel = state.selected && state.selected.group === d.gi && state.selected.index === d.k;
    const isFocusTarget = state.focus && state.focus.small && state.focus.small[0] === d.gi && state.focus.small[1] === d.k;
    if (!g.visible || dwarfFade === 0 || isSel || (state.focus && !isFocusTarget)) { d.el.style.opacity = "0"; continue; }
    const p = project(g.pos[d.k * 3], g.pos[d.k * 3 + 1], g.pos[d.k * 3 + 2]);
    const show = p && p[0] > -40 && p[0] < cssW + 40 && p[1] > -20 && p[1] < cssH + 20 &&
      rectFree(p[0], p[1], d.name);
    d.el.style.opacity = show ? (0.85 * dwarfFade).toFixed(2) : "0";
    if (show) d.el.style.transform = `translate(${p[0]}px, ${p[1]}px) translate(-50%,-150%)`;
  }
  // major-moon labels — only when the moon is visually separated from its parent
  for (const ml of moonLabelList) {
    if (!ml.el) {
      ml.el = document.createElement("div");
      ml.el.className = "pl-label moon-label";
      ml.el.textContent = ml.name;
      labelsEl.appendChild(ml.el);
    }
    if (!moons.visible) { ml.el.style.opacity = "0"; continue; }
    const k = ml.k;
    const p = project(moons.pos[k * 3], moons.pos[k * 3 + 1], moons.pos[k * 3 + 2]);
    if (!p) { ml.el.style.opacity = "0"; continue; }
    let sep = 1e9;
    if (moonParentPos(k, _mpp)) {
      const pp = project(_mpp[0], _mpp[1], _mpp[2]);
      if (pp) sep = Math.hypot(p[0] - pp[0], p[1] - pp[1]);
    }
    const show = sep > 30 && state.selMoon !== ml.k && p[0] > -40 && p[0] < cssW + 40 && p[1] > -20 && p[1] < cssH + 20 &&
      rectFree(p[0], p[1], ml.name);
    ml.el.style.opacity = show ? "0.85" : "0";
    if (show) ml.el.style.transform = `translate(${p[0]}px, ${p[1]}px) translate(-50%,-150%)`;
  }
  // spacecraft labels (always-on while the craft is in its data window)
  for (const c of craft) {
    if (!c.labelEl) {
      c.labelEl = document.createElement("div");
      c.labelEl.className = "pl-label craft-label";
      c.labelEl.textContent = c.name;
      labelsEl.appendChild(c.labelEl);
    }
    if (!craftVisible || !c.inRange || (state.selCraft != null && craft[state.selCraft] === c)) { c.labelEl.style.opacity = "0"; continue; }
    const p = project(c.pos[0], c.pos[1], c.pos[2]);
    const show = p && p[0] > -60 && p[0] < cssW + 60 && p[1] > -20 && p[1] < cssH + 20 &&
      rectFree(p[0], p[1], c.name);
    c.labelEl.style.opacity = show ? "0.9" : "0";
    if (show) c.labelEl.style.transform = `translate(${p[0]}px, ${p[1]}px) translate(-50%,-150%)`;
  }
  // region labels (Kuiper Belt, Oort cloud, interstellar space, Alpha Cen)
  for (const rl of REGION_LABELS) {
    if (!rl.el) {
      rl.el = document.createElement("div");
      rl.el.className = "pl-label " + rl.cls + (rl.onClick ? " region-clickable" : "");
      rl.el.textContent = rl.text;
      if (rl.onClick) rl.el.addEventListener("click", rl.onClick);
      labelsEl.appendChild(rl.el);
    }
    let fade = rl.show() ? clamp((state.cam.dist - rl.in0) / (rl.in1 - rl.in0), 0, 1) : 0;
    if (fade > 0 && rl.out0) fade *= clamp(1 - (state.cam.dist - rl.out0) / (rl.out1 - rl.out0), 0, 1);
    let shown = false;
    if (fade > 0.01) {
      const p = project(rl.pos[0], rl.pos[1], rl.pos[2]);
      if (p && p[0] > 40 && p[0] < cssW - 40 && p[1] > 30 && p[1] < cssH - 30) {
        rl.el.style.opacity = (0.8 * fade).toFixed(2);
        rl.el.style.transform = `translate(${p[0]}px, ${p[1]}px) translate(-50%,-150%)`;
        shown = true;
      } else rl.el.style.opacity = "0";
    } else rl.el.style.opacity = "0";
    // clickable region labels only intercept clicks while actually visible,
    // so an opacity-0 label never swallows a click in empty space
    if (rl.onClick) rl.el.style.pointerEvents = shown ? "auto" : "none";
  }

  // selection marker (asteroid, planet, or moon) + live distance readout
  let selPos = null, selName = "";
  if (state.selected) {
    const g = groups[state.selected.group], k = state.selected.index;
    selPos = [g.pos[k * 3], g.pos[k * 3 + 1], g.pos[k * 3 + 2]];
    selName = g.meta[k].name;
  } else if (state.selPlanet != null) {
    const ps = planetState[state.selPlanet];
    selPos = [ps.pos[0], ps.pos[1], ps.pos[2]];
    selName = ps.def.name;
  } else if (state.selMoon != null) {
    const k = state.selMoon;
    selPos = [moons.pos[k * 3], moons.pos[k * 3 + 1], moons.pos[k * 3 + 2]];
    selName = moons.meta[k].name;
  } else if (state.selSun) {
    selPos = [0, 0, 0];
    selName = "Sun";
  } else if (state.selSat != null && sats.pos[state.selSat * 3] < 1e8) {
    const k = state.selSat;
    selPos = [sats.pos[k * 3], sats.pos[k * 3 + 1], sats.pos[k * 3 + 2]];
    selName = sats.names[k];
  } else if (state.selCraft != null && craft[state.selCraft].inRange) {
    const c = craft[state.selCraft];
    selPos = [c.pos[0], c.pos[1], c.pos[2]];
    selName = c.name;
  }
  if (selPos && state.selSat == null) $("info-r").textContent = Math.hypot(selPos[0], selPos[1], selPos[2]).toFixed(3) + " au";
  // live distance + speed for a selected spacecraft
  if (state.selCraft != null) {
    const c = craft[state.selCraft];
    if (c.inRange) {
      const f = craftFacts(c);
      $("craft-r").textContent = f.r.toFixed(2) + " au";
      $("craft-v").textContent = isFinite(f.v) ? f.v.toFixed(1) + " km/s" : "—";
    } else {
      $("craft-r").textContent = "outside data range";
      $("craft-v").textContent = "—";
    }
  }
  updateScaleBar();
  updateZoomBar();
  if (selPos) {
    if (!selMarkerEl) {
      selMarkerEl = document.createElement("div");
      selMarkerEl.className = "pl-label sel-marker";
      labelsEl.appendChild(selMarkerEl);
    }
    const p = project(selPos[0], selPos[1], selPos[2]);
    selMarkerEl.textContent = selName;
    if (p) {
      selMarkerEl.style.opacity = "1";
      selMarkerEl.style.transform = `translate(${p[0]}px, ${p[1]}px) translate(-50%,-150%)`;
      const w = selName.length * 6.8 + 10;
      selMarkerLastRect = [p[0] - w / 2, p[1] - 26, p[0] + w / 2, p[1] - 4];
    } else { selMarkerEl.style.opacity = "0"; selMarkerLastRect = null; }
  } else {
    if (selMarkerEl) selMarkerEl.style.opacity = "0";
    selMarkerLastRect = null;
  }
}

/* ---- zoom position indicator: log scale over the standard zoom range
   (0.18 au → Oort edge); focused moon-system dives peg the thumb left ---- */
const ZB_LO = Math.log10(0.18), ZB_HI = Math.log10(MAX_DIST);
const zbFrac = (d) => (Math.log10(d) - ZB_LO) / (ZB_HI - ZB_LO);
(function initZoomBar() {
  $("zb-earth").style.left = (zbFrac(1) * 100).toFixed(1) + "%";
  $("zb-neptune").style.left = (zbFrac(30.07) * 100).toFixed(1) + "%";
  $("zb-kuiper").style.left = (zbFrac(44) * 100).toFixed(1) + "%";
  const l = zbFrac(2000) * 100, r = zbFrac(100000) * 100;
  $("zb-oort").style.left = l.toFixed(1) + "%";
  $("zb-oort").style.width = (r - l).toFixed(1) + "%";
})();
const zoomThumbEl = $("zoom-thumb");
function updateZoomBar() {
  const f = clamp(zbFrac(state.cam.dist), 0, 1);
  zoomThumbEl.style.left = (f * 100).toFixed(2) + "%";
  zoomBarEl.setAttribute("aria-valuenow", Math.round(f * 100));
}

/* ---- click / tap / drag the zoom bar to jump to that zoom level (no wheel needed) ---- */
const zoomBarEl = $("zoom-bar");
let zbDragging = false;
function zoomBarTo(clientX) {
  const r = zoomBarEl.getBoundingClientRect();
  if (r.width <= 0) return;
  const f = clamp((clientX - r.left) / r.width, 0, 1);
  const prev = state.cam.tDist;
  state.camEaseRate = 9;
  state.cam.tDist = clamp(Math.pow(10, ZB_LO + f * (ZB_HI - ZB_LO)), minDist(), MAX_DIST);
  farOrientIfCrossed(prev);
}
zoomBarEl.addEventListener("pointerdown", (ev) => {
  zbDragging = true;
  try { zoomBarEl.setPointerCapture(ev.pointerId); } catch (_) { /* best-effort */ }
  zoomBarTo(ev.clientX);
  dismissHint();
  ev.preventDefault();
});
zoomBarEl.addEventListener("pointermove", (ev) => { if (zbDragging) zoomBarTo(ev.clientX); });
function zbEnd(ev) { zbDragging = false; try { zoomBarEl.releasePointerCapture(ev.pointerId); } catch (_) { /* noop */ } }
zoomBarEl.addEventListener("pointerup", zbEnd);
zoomBarEl.addEventListener("pointercancel", zbEnd);
// keyboard operability for the slider (focused only); stopPropagation so arrows
// don't also reach the global cycleSpeed handler
zoomBarEl.addEventListener("keydown", (ev) => {
  const cur = clamp(zbFrac(state.cam.dist), 0, 1);
  let f = cur;
  if (ev.key === "ArrowRight" || ev.key === "ArrowUp") f = clamp(cur + 0.04, 0, 1);
  else if (ev.key === "ArrowLeft" || ev.key === "ArrowDown") f = clamp(cur - 0.04, 0, 1);
  else if (ev.key === "Home") f = 0;
  else if (ev.key === "End") f = 1;
  else return;
  ev.preventDefault(); ev.stopPropagation();
  const prev = state.cam.tDist;
  state.camEaseRate = 9;
  state.cam.tDist = clamp(Math.pow(10, ZB_LO + f * (ZB_HI - ZB_LO)), minDist(), MAX_DIST);
  farOrientIfCrossed(prev);
  dismissHint();
});

/* ---- brand-mark mini-orrery, spin synced to the time machine's speed & direction ---- */
const bmOrbs = ["bm-o1", "bm-o2", "bm-o3"].map((c) => document.querySelector(".brand-mark ." + c));
const bmPeriod = [14, 9, 5.5];   // seconds per revolution at the 1× multiplier (inner = faster)
const bmAng = [0, 0, 0];
const bmReduced = window.matchMedia && window.matchMedia("(prefers-reduced-motion: reduce)").matches;
function updateBrandOrrery(dt) {
  if (bmReduced) return;
  // compress the huge dps range (1 → 1826 d/s) into a lively but legible spin
  const mult = 0.3 + Math.log10(state.dps + 1) * 0.7;
  const sgn = state.playing ? state.dir : 0;   // freeze when paused, reverse when rewinding
  for (let i = 0; i < 3; i++) {
    if (!bmOrbs[i]) continue;
    bmAng[i] = (bmAng[i] + sgn * (360 / bmPeriod[i]) * mult * dt) % 360;
    bmOrbs[i].style.transform = "rotate(" + bmAng[i].toFixed(2) + "deg)";
  }
}

/* ---- scale bar: a nice 1-2-5 length at the focus-plane depth ---- */
let lastScaleTxt = "";
function updateScaleBar() {
  const auPerPx = (2 * state.cam.dist * Math.tan(FOV / 2)) / cssH;
  const raw = auPerPx * 130;                 // target ≈130 px
  let unit = "au", scale = 1;
  if (raw < 0.1) { unit = "km"; scale = 149597870.7; }
  const rawU = raw * scale;
  const p10 = Math.pow(10, Math.floor(Math.log10(rawU)));
  const mant = rawU / p10;
  const nice = (mant < 1.5 ? 1 : mant < 3.5 ? 2 : mant < 7.5 ? 5 : 10) * p10;
  $("scale-line").style.width = (nice / scale / auPerPx).toFixed(1) + "px";
  const txt = nice.toLocaleString("en-US") + " " + unit;
  if (txt !== lastScaleTxt) { $("scale-text").textContent = txt; lastScaleTxt = txt; }
}

/* ---- HUD: animated counter + clock ---- */
const MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
let lastDateStr = "";
function updateHUD() {
  if (state.shownCount !== state.totalLoaded) {
    const d = state.totalLoaded - state.shownCount;
    state.shownCount += Math.abs(d) < 4 ? d : Math.round(d * 0.08);
    $("stat-count").textContent = state.shownCount.toLocaleString("en-US");
  }
  const date = new Date((state.simJD - 2440587.5) * 86400000);
  if (!isFinite(date.getTime())) return;
  const s = `${String(date.getUTCDate()).padStart(2, "0")} ${MONTHS[date.getUTCMonth()]} ${date.getUTCFullYear()} · ${String(date.getUTCHours()).padStart(2, "0")}:${String(date.getUTCMinutes()).padStart(2, "0")} UTC`;
  if (s !== lastDateStr) { $("sim-date").textContent = s; lastDateStr = s; }
}

/* ============================================================
   7. Data loading — NASA/JPL SBDB
   ============================================================ */
function fetchJSON(url) {
  return fetch(url, { headers: { Accept: "application/json" } }).then((r) => {
    if (!r.ok) throw new Error("HTTP " + r.status);
    return r.json();
  });
}

// Build the per-object record kept for the info panel & preview, tagging
// dwarf planets, Sentry-listed objects, and any real imagery we have.
function makeMeta(r, kind) {
  const m = {
    name: r.name, pdes: r.pdes, cls: r.cls, a: r.a, e: r.e, i: r.inc,
    H: r.H, diam: r.diam, albedo: r.albedo, rot: r.rot,
    spec: r.spec || "", pha: !!r.pha, moid: r.moid,
    q: r.a * (1 - r.e), kind: kind || "a", dwarf: null, sentry: null, img: null,
    ucode: r.ucode != null ? r.ucode : null, arc: r.arc, lastObs: r.lastObs || null,
  };
  const dw = DWARFS[r.pdes];
  if (dw) {
    m.dwarf = dw.name;
    if (!isFinite(m.diam)) m.diam = dw.diam;
  }
  m.sentry = sentryMap.get(r.pdes) || sentryMap.get((r.name || "").trim()) || null;
  m.img = REAL_IMG[r.pdes] || null;
  return m;
}

function ingest(resp) {
  if (!resp || !Array.isArray(resp.fields) || !Array.isArray(resp.data)) return 0;
  const ix = {};
  resp.fields.forEach((f, k) => (ix[f] = k));
  const need = ["pdes", "a", "e", "i", "om", "w", "ma", "epoch"];
  for (const f of need) if (!(f in ix)) return 0;
  const num = (v) => (v == null || v === "" ? NaN : +v);

  // bucket rows per group first, then append in bulk
  const buckets = groups.map(() => []);
  for (const row of resp.data) {
    const pdes = row[ix.pdes];
    if (pdes == null || seen.has(pdes)) continue;
    const a = num(row[ix.a]), e = num(row[ix.e]), inc = num(row[ix.i]);
    const om = num(row[ix.om]), w = num(row[ix.w]), ma = num(row[ix.ma]), ep = num(row[ix.epoch]);
    if (!(a > 0) || !(e >= 0) || e >= 0.995 || !isFinite(inc) || !isFinite(om) || !isFinite(w) || !isFinite(ma) || !isFinite(ep)) continue;
    seen.add(pdes);
    const cls = (row[ix.class] || "").trim();
    // dwarf planets get their own population regardless of dynamical class
    const gi = DWARFS[pdes] ? GI.DWF
      : cls in CLASS_TO_GROUP ? CLASS_TO_GROUP[cls] : GI.OTH;
    const H = num(row[ix.H]);
    const diam = num(row[ix.diameter]);
    const name = (row[ix.name] || "").trim() || pdes;
    buckets[gi].push({
      pdes, name, cls, a, e, inc, om, w, ma, ep,
      H: isFinite(H) ? H : NaN, diam: isFinite(diam) ? diam : NaN,
      albedo: num(row[ix.albedo]), rot: num(row[ix.rot_per]),
      spec: (row[ix.spec_T] || row[ix.spec_B] || "").trim(),
      pha: (row[ix.pha] || "") === "Y", moid: num(row[ix.moid]),
      ucode: row[ix.condition_code], arc: num(row[ix.data_arc]), lastObs: row[ix.last_obs] || null,
    });
  }

  let added = 0;
  const Pb = new Float64Array(6);
  buckets.forEach((rows, gi) => {
    if (!rows.length) return;
    const g = groups[gi];
    const n0 = g.count, n1 = n0 + rows.length;
    const el = new Float32Array(n1 * STRIDE); el.set(g.el);
    const pos = new Float32Array(n1 * 3); pos.set(g.pos);
    const sizes = new Float32Array(n1); sizes.set(g.sizes);
    rows.forEach((r, j) => {
      const k = n0 + j, o = k * STRIDE;
      perifocalBasis(r.w * DEG, r.om * DEG, r.inc * DEG, Pb);
      el[o] = r.a;
      el[o + 1] = r.e;
      el[o + 2] = r.a * Math.sqrt(1 - r.e * r.e);
      el[o + 3] = r.ma * DEG;
      el[o + 4] = GAUSS_K / Math.pow(r.a, 1.5);   // rad/day
      el[o + 5] = r.ep - J2000;
      el[o + 6] = Pb[0]; el[o + 7] = Pb[1]; el[o + 8] = Pb[2];
      el[o + 9] = Pb[3]; el[o + 10] = Pb[4]; el[o + 11] = Pb[5];
      const H = isFinite(r.H) ? r.H : 16;
      const m = makeMeta(r, gi === GI.COM ? "c" : "a");
      sizes[k] = Math.min(0.006 * Math.pow(1.32, Math.max(17 - H, 0)), 0.13);
      if (m.dwarf) sizes[k] = Math.max(sizes[k], 0.085);   // make the named dwarfs visible
      g.meta.push(m);
    });
    g.el = el; g.pos = pos; g.sizes = sizes; g.count = n1;
    if (!g.posBuf) g.posBuf = gl.createBuffer();
    if (g.sizeBuf) gl.deleteBuffer(g.sizeBuf);
    g.sizeBuf = makeBuffer(sizes);
    g.dirty = true;
    added += rows.length;
  });
  if (added) {
    state.totalLoaded += added;
    state.needFullUpdate = true;
    renderLegend();
  }
  return added;
}

// Comets are given as perihelion distance q + time of perihelion passage tp.
// Convert to the same {a, M, n, epoch} record the propagator already uses:
// a = q/(1-e); mean anomaly is 0 at perihelion, so reference the epoch to tp.
// e < 0.995 → periodic comets; 0.995 ≤ e < 0.9999 → long-period comets, the
// Oort-cloud evidence layer. Truly parabolic/hyperbolic records are skipped.
function ingestComets(resp) {
  if (!resp || !Array.isArray(resp.fields) || !Array.isArray(resp.data)) return 0;
  const ix = {};
  resp.fields.forEach((f, k) => (ix[f] = k));
  for (const f of ["pdes", "e", "i", "om", "w", "q", "tp"]) if (!(f in ix)) return 0;
  const num = (v) => (v == null || v === "" ? NaN : +v);

  const buckets = { [GI.COM]: [], [GI.LPC]: [] };
  for (const row of resp.data) {
    const pdes = row[ix.pdes];
    if (pdes == null || seen.has(pdes)) continue;
    const e = num(row[ix.e]), q = num(row[ix.q]), tp = num(row[ix.tp]);
    const inc = num(row[ix.i]), om = num(row[ix.om]), w = num(row[ix.w]);
    if (!(q > 0) || !(e >= 0) || e >= 0.9999 ||
        !isFinite(inc) || !isFinite(om) || !isFinite(w) || !isFinite(tp)) continue;
    const a = q / (1 - e);
    if (a > 80000) continue;                     // beyond even the Oort shell
    seen.add(pdes);
    const lpc = e >= 0.995;
    const diam = num(row[ix.diameter]);
    buckets[lpc ? GI.LPC : GI.COM].push({
      pdes, name: (row[ix.name] || "").trim() || pdes, cls: lpc ? "LPC" : "COM",
      a, e, inc, om, w, tp,
      diam: isFinite(diam) ? diam : NaN, rot: num(row[ix.rot_per]),
      albedo: NaN, spec: "", pha: false, moid: NaN,
      ucode: row[ix.condition_code], arc: num(row[ix.data_arc]), lastObs: row[ix.last_obs] || null,
    });
  }

  let added = 0;
  const Pb = new Float64Array(6);
  for (const giStr of Object.keys(buckets)) {
    const gi = +giStr, rows = buckets[gi];
    if (!rows.length) continue;
    const g = groups[gi];
    const n0 = g.count, n1 = n0 + rows.length;
    const el = new Float32Array(n1 * STRIDE); el.set(g.el);
    const pos = new Float32Array(n1 * 3); pos.set(g.pos);
    const sizes = new Float32Array(n1); sizes.set(g.sizes);
    rows.forEach((r, j) => {
      const k = n0 + j, o = k * STRIDE;
      perifocalBasis(r.w * DEG, r.om * DEG, r.inc * DEG, Pb);
      el[o] = r.a;
      el[o + 1] = r.e;
      el[o + 2] = r.a * Math.sqrt(1 - r.e * r.e);
      el[o + 3] = 0;                               // M = 0 at perihelion
      el[o + 4] = GAUSS_K / Math.pow(r.a, 1.5);    // rad/day
      el[o + 5] = r.tp - J2000;                    // reference epoch = perihelion passage
      el[o + 6] = Pb[0]; el[o + 7] = Pb[1]; el[o + 8] = Pb[2];
      el[o + 9] = Pb[3]; el[o + 10] = Pb[4]; el[o + 11] = Pb[5];
      sizes[k] = 0.05;
      g.meta.push(makeMeta(r, "c"));
    });
    g.el = el; g.pos = pos; g.sizes = sizes; g.count = n1;
    if (gi === GI.LPC) g.kepIters = 28;            // near-parabolic Newton is slow
    if (!g.posBuf) g.posBuf = gl.createBuffer();
    if (g.sizeBuf) gl.deleteBuffer(g.sizeBuf);
    g.sizeBuf = makeBuffer(sizes);
    g.dirty = true;
    added += rows.length;
  }
  if (added) {
    state.totalLoaded += added;
    state.needFullUpdate = true;
    renderLegend();
  }
  return added;
}

// Interstellar visitors: e > 1, so the orbit is an open hyperbola and the body
// never returns. Same q/tp → record conversion as comets, but a is negative;
// the group is flagged so the propagator uses hyperbolic Kepler.
function ingestInterstellar(resp) {
  if (!resp || !Array.isArray(resp.fields) || !Array.isArray(resp.data)) return 0;
  const ix = {};
  resp.fields.forEach((f, k) => (ix[f] = k));
  for (const f of ["pdes", "e", "i", "om", "w", "q", "tp"]) if (!(f in ix)) return 0;
  const num = (v) => (v == null || v === "" ? NaN : +v);

  const rows = [];
  for (const row of resp.data) {
    const pdes = row[ix.pdes];
    if (pdes == null || seen.has(pdes)) continue;
    const e = num(row[ix.e]), q = num(row[ix.q]), tp = num(row[ix.tp]);
    const inc = num(row[ix.i]), om = num(row[ix.om]), w = num(row[ix.w]);
    if (!(q > 0) || !(e > 1) || !isFinite(inc) || !isFinite(om) || !isFinite(w) || !isFinite(tp)) continue;
    seen.add(pdes);
    const diam = num(row[ix.diameter]);
    rows.push({
      pdes, name: (row[ix.name] || "").trim() || pdes, cls: "ISO",
      a: q / (1 - e), e, inc, om, w, tp,
      diam: isFinite(diam) ? diam : NaN, rot: num(row[ix.rot_per]),
      albedo: NaN, spec: "", pha: false, moid: NaN,
      ucode: row[ix.condition_code], arc: num(row[ix.data_arc]), lastObs: row[ix.last_obs] || null,
    });
  }
  if (!rows.length) return 0;

  const g = groups[GI.ISO];
  g.hyperbolic = true;
  const n0 = g.count, n1 = n0 + rows.length;
  const el = new Float32Array(n1 * STRIDE); el.set(g.el);
  const pos = new Float32Array(n1 * 3); pos.set(g.pos);
  const sizes = new Float32Array(n1); sizes.set(g.sizes);
  const Pb = new Float64Array(6);
  rows.forEach((r, j) => {
    const k = n0 + j, o = k * STRIDE;
    const aAbs = Math.abs(r.a);
    perifocalBasis(r.w * DEG, r.om * DEG, r.inc * DEG, Pb);
    el[o] = r.a;                                  // negative for hyperbola
    el[o + 1] = r.e;
    el[o + 2] = aAbs * Math.sqrt(r.e * r.e - 1);  // b
    el[o + 3] = 0;                                // M = 0 at perihelion
    el[o + 4] = GAUSS_K / Math.pow(aAbs, 1.5);    // rad/day
    el[o + 5] = r.tp - J2000;                     // reference epoch = perihelion passage
    el[o + 6] = Pb[0]; el[o + 7] = Pb[1]; el[o + 8] = Pb[2];
    el[o + 9] = Pb[3]; el[o + 10] = Pb[4]; el[o + 11] = Pb[5];
    sizes[k] = 0.07;
    g.meta.push(makeMeta(r, "i"));
  });
  g.el = el; g.pos = pos; g.sizes = sizes; g.count = n1;
  if (!g.posBuf) g.posBuf = gl.createBuffer();
  if (g.sizeBuf) gl.deleteBuffer(g.sizeBuf);
  g.sizeBuf = makeBuffer(sizes);
  g.dirty = true;
  state.totalLoaded += rows.length;
  state.needFullUpdate = true;
  renderLegend();
  return rows.length;
}

// Planetary moons from data/moons.json: Keplerian elements relative to the
// parent planet (ecliptic J2000, from Horizons), packed into the same STRIDE
// layout the asteroid propagator uses. Parent digit 3..8 → planetState index;
// 9 (Pluto) resolves to its entry in the asteroid groups.
const PARENT_NAMES = { 3: "Earth", 4: "Mars", 5: "Jupiter", 6: "Saturn", 7: "Uranus", 8: "Neptune", 9: "Pluto" };
function buildMoons(data) {
  if (!data || !Array.isArray(data.moons) || !data.moons.length) return 0;
  const list = data.moons;
  const n = list.length;
  moons.el = new Float32Array(n * STRIDE);
  moons.pos = new Float32Array(n * 3);
  moons.sizes = new Float32Array(n);
  moons.parentIdx = new Int16Array(n);
  moons.parentSmall = new Array(n).fill(null);
  moons.meta = [];
  const Pb = new Float64Array(6);
  list.forEach((r, k) => {
    const o = k * STRIDE;
    perifocalBasis(r.w * DEG, r.om * DEG, r.i * DEG, Pb);
    moons.el[o] = r.a;
    moons.el[o + 1] = r.e;
    moons.el[o + 2] = r.a * Math.sqrt(1 - r.e * r.e);
    moons.el[o + 3] = r.ma * DEG;
    moons.el[o + 4] = r.n * DEG;            // Horizons n is deg/day
    moons.el[o + 5] = r.epoch - J2000;
    moons.el[o + 6] = Pb[0]; moons.el[o + 7] = Pb[1]; moons.el[o + 8] = Pb[2];
    moons.el[o + 9] = Pb[3]; moons.el[o + 10] = Pb[4]; moons.el[o + 11] = Pb[5];
    moons.parentIdx[k] = r.parent >= 3 && r.parent <= 8 ? r.parent - 1 : -1;
    if (r.parent === 9 && plutoRef) moons.parentSmall[k] = plutoRef;
    moons.sizes[k] = 0.0015 + Math.min((r.radius || 8) / 2700, 1) * 0.02;
    const img = MOON_IMG[r.name.toLowerCase()] || null;
    moons.meta.push({
      name: r.name, parentName: PARENT_NAMES[r.parent] || "?",
      a: r.a, e: r.e, i: r.i, n: r.n, radius: r.radius, img,
    });
    if (r.radius >= 200) moonLabelList.push({ k, name: r.name, el: null });
  });
  moons.count = n;
  if (!moons.posBuf) moons.posBuf = gl.createBuffer();
  if (moons.sizeBuf) gl.deleteBuffer(moons.sizeBuf);
  moons.sizeBuf = makeBuffer(moons.sizes);
  moons.dirty = true;
  return n;
}

/* Lazy-load the active-satellite catalogue on first arrival at Earth. */
function loadSatellites() {
  if (sats.loaded || sats.loading || typeof satellite === "undefined") return;
  sats.loading = true;
  fetchJSON("data/satellites.json").then((data) => {
    const list = (data && data.sats) || [];
    for (const [name, l1, l2] of list) {
      const rec = satellite.twoline2satrec(l1, l2);
      if (rec.error) continue;
      sats.recs.push(rec);
      sats.names.push(name);
    }
    sats.count = sats.recs.length;
    sats.payloadCount = data.payloadCount != null ? Math.min(data.payloadCount, sats.count) : sats.count;
    sats.pos = new Float32Array(sats.count * 3);
    sats.rel = new Float32Array(sats.count * 3);   // Earth-relative offsets (SGP4); pos = ep + rel each frame
    sats.sizes = new Float32Array(sats.count).fill(0.012);
    sats.posBuf = gl.createBuffer();
    sats.sizeBuf = makeBuffer(sats.sizes);
    sats.epochJD = data.generated ? new Date(data.generated).getTime() / 86400000 + 2440587.5 : jdNow();
    sats.loaded = true;
    sats.loading = false;
    sats.fullUpdate = true;
    satCountKnown = sats.payloadCount;
    debrisCountKnown = sats.count - sats.payloadCount;
    updateManMade();
    renderLegend();
    // if we're already at Earth, surface the count badge now that data's in
    if (state.selPlanet === EARTH_IDX && !$("info-badges").querySelector(".badge-sat")) {
      $("info-badges").insertAdjacentHTML("beforeend",
        `<span class="badge badge-sat">🛰 ${sats.payloadCount.toLocaleString("en-US")} satellites</span>`);
    }
  }).catch((err) => { console.warn("satellite load failed:", err); sats.loading = false; });
}

// Satellites are shown LIVE (real-time), not on the sim clock — you can't watch a
// 90-minute LEO orbit at 1 month/second. So the fade only guards data staleness:
// full while the daily TLE snapshot is fresh, gone if it ever falls >2 yr behind.
function satEpochFade() {
  const d = Math.abs(jdNow() - sats.epochJD);       // days since the TLE snapshot
  return clamp(1 - (d - 365) / 365, 0, 1);
}

// SGP4 → TEME km → ecliptic au → + Earth's heliocentric position. Propagated to
// WALL-CLOCK now (live), so the rings drift smoothly instead of aliasing into
// staggered dashes under the fast time machine. Sliced so the ~15k-object pass
// doesn't stall the frame; at real-time the per-slice lag is sub-pixel.
const SAT_SLICES = 4;
const _satDate = new Date();
function updateSatellites() {
  if (!sats.count) return;
  _satDate.setTime(Date.now());
  const ce = Math.cos(OBLIQ), se = Math.sin(OBLIQ);
  // SGP4 the EARTH-RELATIVE offsets (TEME→ecliptic au) — expensive, so sliced.
  const slices = sats.fullUpdate ? 1 : SAT_SLICES;
  const span = Math.ceil(sats.count / slices);
  const from = sats.fullUpdate ? 0 : sats.sliceCursor * span;
  const to = Math.min(from + span, sats.count);
  for (let k = from; k < to; k++) {
    const pv = satellite.propagate(sats.recs[k], _satDate);
    const p = pv && pv.position;
    if (!p || !isFinite(p.x)) { sats.rel[k * 3] = 1e9; continue; }
    const xe = p.x / AU_KM, ye = p.y / AU_KM, ze = p.z / AU_KM;
    sats.rel[k * 3] = xe;
    sats.rel[k * 3 + 1] = ye * ce + ze * se;
    sats.rel[k * 3 + 2] = -ye * se + ze * ce;
  }
  sats.sliceCursor = (sats.sliceCursor + 1) % SAT_SLICES;
  sats.fullUpdate = false;
  // Recombine ALL satellites with the CURRENT Earth position every frame (cheap
  // additions, no SGP4). Without this, sliced updates bake a stale Earth position
  // into 3/4 of the points, tearing the swarm into ghost clusters as Earth moves.
  const ep = planetState[EARTH_IDX].pos;
  for (let k = 0; k < sats.count; k++) {
    if (sats.rel[k * 3] > 1e8) { sats.pos[k * 3] = 1e9; continue; }
    sats.pos[k * 3] = ep[0] + sats.rel[k * 3];
    sats.pos[k * 3 + 1] = ep[1] + sats.rel[k * 3 + 1];
    sats.pos[k * 3 + 2] = ep[2] + sats.rel[k * 3 + 2];
  }
  sats.dirty = true;
}

// Spacecraft trajectories from Horizons → flat arrays for interpolation.
function ingestCraft(data) {
  craft.length = 0;
  const list = (data && data.craft) || [];
  for (const c of list) {
    const raw = c.points || [];
    const n = raw.length;
    if (!n) continue;
    const jd = new Float64Array(n);
    const verts = new Float32Array(n * 3);
    for (let i = 0; i < n; i++) {
      jd[i] = raw[i][0];
      verts[i * 3] = raw[i][1]; verts[i * 3 + 1] = raw[i][2]; verts[i * 3 + 2] = raw[i][3];
    }
    craft.push({
      name: c.name, info: c, jd, verts, n,
      pos: new Float32Array(3), inRange: false,
      trajBuf: makeBuffer(verts),
      posBuf: makeBuffer(new Float32Array(3), gl.DYNAMIC_DRAW),
      sizeBuf: makeBuffer(new Float32Array([0.06])),
      labelEl: null,
    });
  }
  return craft.length;
}
// linear-interpolate a craft's heliocentric position at sim JD
function craftPositionAt(c, jd) {
  if (jd < c.jd[0] || jd > c.jd[c.n - 1]) { c.inRange = false; return; }
  c.inRange = true;
  let lo = 0, hi = c.n - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (c.jd[mid] <= jd) lo = mid; else hi = mid; }
  const t0 = c.jd[lo], t1 = c.jd[hi];
  const f = t1 > t0 ? (jd - t0) / (t1 - t0) : 0;
  for (let d = 0; d < 3; d++) c.pos[d] = c.verts[lo * 3 + d] + (c.verts[hi * 3 + d] - c.verts[lo * 3 + d]) * f;
}
function craftFacts(c) {
  const r = Math.hypot(c.pos[0], c.pos[1], c.pos[2]);
  // speed from the bracketing samples (au/day → km/s)
  let lo = 0, hi = c.n - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (c.jd[mid] <= state.simJD) lo = mid; else hi = mid; }
  const dt = c.jd[hi] - c.jd[lo];
  let v = NaN;
  if (dt > 0) {
    const dx = c.verts[hi * 3] - c.verts[lo * 3], dy = c.verts[hi * 3 + 1] - c.verts[lo * 3 + 1], dz = c.verts[hi * 3 + 2] - c.verts[lo * 3 + 2];
    v = Math.hypot(dx, dy, dz) / dt * 1731.46;   // 1 au/day = 1731.46 km/s
  }
  return { r, v };
}

// derive orbit regime + altitudes from a satrec for the info card
function satFacts(k) {
  const rec = sats.recs[k];
  const a = Math.pow(398600.4418 / (rec.no * rec.no / 3600), 1 / 3);  // semi-major (km), no in rad/min
  const peri = a * (1 - rec.ecco) - 6371;
  const apo = a * (1 + rec.ecco) - 6371;
  const periodMin = (2 * Math.PI) / rec.no;
  const meanAlt = (peri + apo) / 2;
  const regime = meanAlt < 2000 ? "Low Earth orbit"
    : meanAlt < 34000 ? "Medium Earth orbit"
    : meanAlt < 37000 ? "Geostationary / GEO"
    : "High / elliptical orbit";
  return { regime, peri, apo, periodMin, inc: rec.inclo / DEG, norad: rec.satnum };
}

let loadedOnce = false;
async function loadAsteroids() {
  const fill = $("loader-fill");
  const status = $("loader-status");
  $("loader-error").hidden = true;
  status.textContent = "downloading NASA/JPL data snapshot…";
  fill.style.width = "15%";
  try {
    const [snap, moonData, craftData] = await Promise.all([
      fetchJSON(DATA_URL),
      fetchJSON("data/moons.json").catch(() => null),       // moons are optional
      fetchJSON("data/spacecraft.json").catch(() => null),  // spacecraft optional
    ]);
    fill.style.width = "60%";
    status.textContent = "propagating orbits…";
    buildSentryMap(snap.sentry);           // before ingest so objects get flagged
    let added = 0;
    for (const q of snap.queries || []) added += ingest(q);
    added += ingestComets(snap.comets);
    added += ingestInterstellar(snap.interstellar);
    if (!added) throw new Error("snapshot contained no usable records");
    buildDwarfList();
    buildPhaList();
    buildMoons(moonData);
    ingestCraft(craftData);
    // the headline "natural objects" counts the small bodies plus moons,
    // the eight planets and the Sun (satellites/spacecraft are tallied separately)
    state.totalLoaded += moons.count + planetState.length + 1;
    updateManMade();
    renderLegend();
    fill.style.width = "100%";
    const total = +snap.totalKnown;
    if (total > 0) $("stat-known").textContent = total.toLocaleString("en-US");
    if (snap.generated) {
      const d = new Date(snap.generated);
      if (isFinite(d.getTime())) {
        $("snap-date").textContent =
          ` · snapshot ${String(d.getUTCDate()).padStart(2, "0")} ${MONTHS[d.getUTCMonth()]} ${d.getUTCFullYear()}`;
      }
    }
    renderCloseApproaches(snap.cad);
    renderSentry(snap.sentry);
    renderEvents(snap);
    loaderReady();
    applyState(location.hash);        // restore a shared/deep-linked view, if any
    lastHash = encodeState();         // baseline so the writer doesn't immediately rewrite a clean URL
  } catch (err) {
    console.warn("snapshot load failed:", err);
    status.textContent = "";
    $("loader-error-msg").textContent =
      "Could not load the asteroid data snapshot (data/asteroids.json). Check your connection and try again.";
    $("loader-error").hidden = false;
  }
}

/* Title screen holds until the user taps/presses, then the dissolve hands
   off to a camera dolly: start far above the system and dive to the
   default view in one continuous shot. */
function loaderReady() {
  loadedOnce = true;
  const loader = $("loader");
  loader.classList.add("ready");
  $("loader-status").textContent = "ready";
  $("loader-begin").hidden = false;
  function launch() {
    loader.removeEventListener("click", launch);
    window.removeEventListener("keydown", onKey);
    loader.classList.add("done");
    setTimeout(loadSatellites, 5000);   // pre-warm the live-satellite layer for search & Earth
    document.body.classList.add("launched");
    if (skipDolly) {                 // arrived via a deep link — keep the restored view, no dive
      state.camEaseRate = 9;
      setTimeout(dismissHint, 6000);
      return;
    }
    const c = state.cam;
    c.dist = 205;
    c.pitch = 1.28;
    c.yaw = c.tYaw + 0.6;
    state.camEaseRate = 1.5;        // slow exponential ease ≈ 2.5 s dive
    setTimeout(dismissHint, 9000);  // start the hint clock at launch, not load
  }
  function onKey(ev) {
    if (["Shift", "Control", "Alt", "Meta"].includes(ev.key)) return;
    launch();
  }
  // "click" (not pointerdown): the release must land on the loader too, or it
  // falls through to the canvas as an orphan pointerup and picks an object
  loader.addEventListener("click", launch);
  window.addEventListener("keydown", onKey);
}

/* the mobile bottom sheet / FAB rail react to the info panel's visibility */
new MutationObserver(() => {
  document.body.classList.toggle("sheet-open", !$("panel-info").hidden);
}).observe($("panel-info"), { attributes: true, attributeFilter: ["hidden"] });

/* ---- close approaches (JPL CNEOS, from the snapshot) ---- */
function renderCloseApproaches(resp) {
  const list = $("cad-list");
  try {
    const ix = {};
    (resp.fields || []).forEach((f, k) => (ix[f] = k));
    if (!("des" in ix) || !("jd" in ix)) throw new Error("no cad data");
    const now = jdNow();
    const rows = (resp.data || []).filter((r) => +r[ix.jd] >= now - 0.5).slice(0, 30);
    if (!rows.length) {
      list.innerHTML = '<li class="cad-empty">no approaches within 10 LD in the next 60 days</li>';
      return;
    }
    list.innerHTML = "";
    for (const r of rows) {
      const des = r[ix.des] || "?";
      const cd = (r[ix.cd] || "").replace(/^(\d{4})-(\w{3})-(\d{2})\s+(\d{2}:\d{2})$/, "$3 $2 $1 · $4 UTC");
      const ld = (+r[ix.dist] || 0) / 0.002569;
      const v = +r[ix.v_rel];
      const cls = ld < 1 ? "cad-near" : ld < 4 ? "cad-mid" : "cad-far";
      const li = document.createElement("li");
      li.innerHTML =
        `<span class="cad-name"></span>` +
        `<span class="cad-dist ${cls}">${ld.toFixed(2)} LD<small>${isFinite(v) ? v.toFixed(1) + " km/s" : ""}</small></span>` +
        `<span class="cad-when">${cd}</span>`;
      li.querySelector(".cad-name").textContent = des;
      list.appendChild(li);
    }
  } catch (err) {
    list.innerHTML = '<li class="cad-empty">close-approach feed unavailable</li>';
  }
}

/* ============================================================
   "Upcoming" events feed — all computed client-side
   ============================================================ */
const EV_DOT = { cad: "#ff6b6b", com: "#8ce9ff", opp: "#9ec5ff", met: "#fde68a" };
const SHOWERS = [
  { m: 0, d: 3, name: "Quadrantids", zhr: 110 }, { m: 3, d: 22, name: "Lyrids", zhr: 18 },
  { m: 4, d: 6, name: "Eta Aquariids", zhr: 50 }, { m: 7, d: 12, name: "Perseids", zhr: 100 },
  { m: 9, d: 21, name: "Orionids", zhr: 20 }, { m: 10, d: 17, name: "Leonids", zhr: 15 },
  { m: 11, d: 14, name: "Geminids", zhr: 150 }, { m: 11, d: 22, name: "Ursids", zhr: 10 },
];
function jumpToJD(jd) { state.simJD = jd; state.needFullUpdate = true; buildPlanetOrbits(); }
function selectByName(name) {
  const q = (name || "").toLowerCase().replace(/[()]/g, "").trim();
  if (!q) return false;
  for (let gi = 0; gi < groups.length; gi++) {
    const g = groups[gi];
    for (let k = 0; k < g.count; k++) {
      const nm = g.meta[k].name.toLowerCase();
      if (nm === q || nm.includes(q)) { if (!g.visible) { g.visible = true; renderLegend(); } selectObject(gi, k); return true; }
    }
  }
  return false;
}
function evCloseApproaches(snap) {
  const out = [], c = snap.cad; if (!c || !c.fields) return out;
  const ix = {}; c.fields.forEach((f, k) => (ix[f] = k)); const now = jdNow();
  for (const r of (c.data || [])) {
    const jd = +r[ix.jd]; if (!(jd >= now - 0.5)) continue;
    const des = r[ix.des] || "?", ld = (+r[ix.dist] || 0) / 0.002569;
    out.push({ jd, kind: "cad", label: des + " — Earth flyby", sub: ld.toFixed(2) + " LD",
      act: () => { if (!selectByName(des)) selectPlanet(EARTH_IDX); jumpToJD(jd); } });
  }
  return out;
}
function evCometPerihelia(now, end) {
  const out = [];
  for (const gi of [GI.COM, GI.LPC]) {
    const g = groups[gi];
    for (let k = 0; k < g.count; k++) {
      const o = k * STRIDE, tp = g.el[o + 5] + J2000, n = g.el[o + 4];
      if (!(n > 0)) continue;
      let peri = tp;
      if (gi === GI.COM) { const period = TWO_PI / n; peri = tp + Math.ceil((now - tp) / period) * period; }
      if (peri < now - 1 || peri > end) continue;
      const m = g.meta[k], ggi = gi, kk = k;
      // "238P/Read", not a bare "Read" that scans as a UI verb
      const evName = m.pdes && !m.name.includes(String(m.pdes)) && /^\d+[PDI]/.test(String(m.pdes))
        ? m.pdes + "/" + m.name : m.name;
      out.push({ jd: peri, kind: "com", label: evName + " — perihelion", sub: "q " + m.q.toFixed(2) + " au",
        act: () => { if (!groups[ggi].visible) { groups[ggi].visible = true; renderLegend(); } selectObject(ggi, kk); jumpToJD(peri); } });
    }
  }
  out.sort((a, b) => a.jd - b.jd);
  return out.slice(0, 14);   // soonest few, so comets don't flood the list
}
function dLon(p, jd) {
  const E = [0, 0, 0], M = [0, 0, 0];
  planetPosition(PLANETS[EARTH_IDX], jd, E); planetPosition(p, jd, M);
  let d = Math.atan2(M[1], M[0]) - Math.atan2(E[1], E[0]);
  while (d > Math.PI) d -= TWO_PI; while (d < -Math.PI) d += TWO_PI;
  return d;
}
function evOppositions(now, end) {
  const out = [];
  for (let i = 3; i < PLANETS.length; i++) {            // Mars … Neptune (superior planets)
    const p = PLANETS[i]; let prev = dLon(p, now), found = false;
    for (let jd = now + 2; jd <= end && !found; jd += 2) {
      const cur = dLon(p, jd);
      if ((prev <= 0 && cur > 0) || (prev >= 0 && cur < 0)) {   // heliocentric lons coincide → opposition
        let lo = jd - 2, hi = jd;
        for (let b = 0; b < 8; b++) { const mid = (lo + hi) / 2; if (Math.sign(dLon(p, mid)) === Math.sign(prev)) lo = mid; else hi = mid; }
        const oj = (lo + hi) / 2, ii = i;
        out.push({ jd: oj, kind: "opp", label: p.name + " at opposition", sub: "closest & brightest", act: () => { selectPlanet(ii); jumpToJD(oj); } });
        found = true;
      }
      prev = cur;
    }
  }
  return out;
}
function evShowers(now, end) {
  const out = [], y0 = new Date((now - 2440587.5) * 86400000).getUTCFullYear();
  for (const s of SHOWERS) {
    for (const y of [y0, y0 + 1]) {
      const jd = new Date(Date.UTC(y, s.m, s.d, 4)).getTime() / 86400000 + 2440587.5;
      if (jd < now - 1 || jd > end) continue;
      out.push({ jd, kind: "met", label: s.name + " — meteor peak", sub: "~" + s.zhr + "/hr", act: () => jumpToJD(jd) });
      break;
    }
  }
  return out;
}
function renderEvents(snap) {
  const list = $("events-list"); if (!list) return;
  const now = jdNow(), end = now + 800;
  let evs = [...evCloseApproaches(snap), ...evCometPerihelia(now, end), ...evOppositions(now, end), ...evShowers(now, end)];
  evs = evs.filter((e) => e.jd >= now - 1 && e.jd <= end).sort((a, b) => a.jd - b.jd);
  if (!evs.length) { list.innerHTML = '<li class="cad-empty">no events in the window</li>'; return; }
  list.innerHTML = "";
  for (const e of evs) {
    const d = new Date((e.jd - 2440587.5) * 86400000);
    const when = `${String(d.getUTCDate()).padStart(2, "0")} ${MONTHS[d.getUTCMonth()]} ${d.getUTCFullYear()}`;
    const li = document.createElement("li");
    li.dataset.jd = e.jd;
    li.innerHTML = `<span class="ev-dot" style="background:${EV_DOT[e.kind]};box-shadow:0 0 7px ${EV_DOT[e.kind]}"></span>` +
      `<span class="ev-name"></span><span class="ev-sub">${e.sub}</span><span class="ev-date">${when}</span>`;
    li.querySelector(".ev-name").textContent = e.label;
    li.addEventListener("click", e.act);
    list.appendChild(li);
  }
}
/* events already behind the simulated clock gray out as "past" */
setInterval(() => {
  const now = state.simJD;
  document.querySelectorAll("#events-list li[data-jd]").forEach((li) => {
    li.classList.toggle("past", +li.dataset.jd < now - 0.5);
  });
}, 2500);

/* ---- Sentry impact-risk watchlist (JPL CNEOS, from the snapshot) ---- */
function buildSentryMap(sentry) {
  sentryMap.clear();
  for (const o of (sentry && sentry.objects) || []) {
    if (o.des != null) sentryMap.set(String(o.des), o);
  }
}
function shortNum(n) {
  if (n >= 1e6) return (n / 1e6).toFixed(n < 1e7 ? 1 : 0) + "M";
  if (n >= 1e3) return Math.round(n / 1e3) + "k";
  return String(n);
}
/* fly to a Sentry object if it happens to be in the rendered set */
function sentrySelect(des) {
  for (let gi = 0; gi < groups.length; gi++) {
    const g = groups[gi];
    for (let k = 0; k < g.count; k++) {
      const sr = g.meta[k].sentry;
      if (sr && String(sr.des) === String(des)) {
        if (!g.visible) { g.visible = true; renderLegend(); }
        selectObject(gi, k);
        if (isFinite(g.meta[k].a)) state.cam.tDist = clamp(g.meta[k].a * 2.4, 1.5, 240);
        return true;
      }
    }
  }
  return false;
}
function renderSentry(sentry) {
  const list = $("sentry-list");
  const objs = (sentry && sentry.objects) || [];
  if (!objs.length) { list.innerHTML = '<li class="cad-empty">sentry risk table unavailable</li>'; return; }
  $("sentry-count").textContent = objs.length.toLocaleString("en-US");
  list.innerHTML = "";
  for (const o of objs.slice(0, 50)) {
    const name = (o.name || o.des || "?").replace(/^\((.*)\)$/, "$1");
    const torino = o.ts != null && o.ts > 0;
    const odds = o.ip > 0 ? "1 in " + shortNum(Math.round(1 / o.ip)) : "—";
    const diam = isFinite(o.diam) ? formatKm(o.diam) : "—";
    const cls = torino ? "sn-hi" : (o.ps != null && o.ps > -2 ? "sn-mid" : "sn-lo");
    const li = document.createElement("li");
    li.innerHTML =
      `<span class="sn-name"></span>` +
      `<span class="sn-risk ${cls}">${torino ? "Torino " + o.ts : "Palermo " + (o.ps != null ? o.ps.toFixed(1) : "?")}</span>` +
      `<span class="sn-meta">${diam} · ${odds}${o.range ? " · " + o.range : ""}</span>`;
    li.querySelector(".sn-name").textContent = name;
    li.addEventListener("click", () => {
      if (!sentrySelect(o.des)) {
        li.classList.add("sn-flash");
        setTimeout(() => li.classList.remove("sn-flash"), 600);
      }
    });
    list.appendChild(li);
  }
}

/* ============================================================
   8. Legend / filters
   ============================================================ */
/* One legend row: dot, label, count, an ⓘ that expands a description,
   and the row-click visibility toggle. */
function legendRow(label, count, css, desc, onToggle, isOn) {
  const li = document.createElement("li");
  li.tabIndex = 0;
  li.setAttribute("role", "button");
  li.setAttribute("aria-pressed", String(!!isOn));
  if (!isOn) li.classList.add("off");
  li.innerHTML =
    `<span class="dot" style="background:${css};box-shadow:0 0 8px ${css}"></span>` +
    `<span class="lname">${label}</span>` +
    `<span class="lcount">${count.toLocaleString("en-US")}</span>` +
    `<button class="linfo" aria-label="About ${label}" title="What is this?" tabindex="-1">ⓘ</button>` +
    `<p class="ldesc" hidden></p>`;
  li.querySelector(".ldesc").textContent = desc;
  li.addEventListener("click", (ev) => {
    if (ev.target.closest(".linfo")) {
      const d = li.querySelector(".ldesc");
      d.hidden = !d.hidden;
      return;
    }
    if (ev.target.closest(".ldesc")) return;   // reading isn't toggling
    onToggle(li);
    li.setAttribute("aria-pressed", String(!li.classList.contains("off")));
  });
  li.addEventListener("keydown", (ev) => {
    if (ev.key === "Enter" || ev.key === " ") { ev.preventDefault(); li.click(); }
  });
  return li;
}

function renderLegend() {
  const ul = $("legend-list");
  ul.innerHTML = "";
  ul.appendChild(legendRow("Sun", 1, "#ffd479",
    "Our star — a G2V main-sequence star holding 99.86% of the solar system's mass. Click it for details.",
    (li) => {
      state.showSun = !state.showSun;
      li.classList.toggle("off", !state.showSun);
      if (!state.showSun && state.selSun) clearSelection();
    }, state.showSun));
  ul.appendChild(legendRow("Planets", PLANETS.length, "#9ec5ff",
    "The eight major planets, placed by the JPL approximate ephemeris (date-accurate 1800–2050). Click one and the camera flies to it.",
    (li) => {
      state.showPlanets = !state.showPlanets;
      li.classList.toggle("off", !state.showPlanets);
      if (!state.showPlanets && state.selPlanet != null) clearSelection();
    }, state.showPlanets));
  // Moons sit right under Planets (they're natural bodies, not an overlay)
  const moonsRow = legendRow("Moons", moons.count, "#c7d2e8",
    "Every planetary satellite in JPL Horizons, orbiting its parent planet. Click a planet and zoom in to explore its moon system.",
    (el) => { moons.visible = !moons.visible; if (!moons.visible && state.selMoon != null) clearSelection(); el.classList.toggle("off", !moons.visible); },
    moons.visible);
  moonsRow.classList.add("legend-toggle");
  ul.appendChild(moonsRow);
  groups.forEach((g, gi) => {
    if (!g.count) return;
    ul.appendChild(legendRow(g.label, g.count, g.css, g.desc, (li) => {
      g.visible = !g.visible;
      li.classList.toggle("off", !g.visible);
      if (!g.visible && state.selected && state.selected.group === gi) clearSelection();
    }, g.visible));
  });

  ul.appendChild(legendRow("Oort cloud", "inferred", "#9aa7c4",
    "A statistical representation, not data — no Oort-cloud object has ever been directly observed. Its existence and extent (≈2,000–100,000 au) are inferred from the orbits of long-period comets. Zoom all the way out to make the journey — and once you're past the shell, look around: there's a signpost to our nearest neighbouring star.",
    (li) => { OORT.visible = !OORT.visible; li.classList.toggle("off", !OORT.visible); }, OORT.visible));

  // overlay toggles (not populations)
  const overlays = [
    { label: "True scale", always: true, css: "#7dd3fc", count: "—",
      desc: "Render every body at its real angular size instead of an exaggerated dot. The Sun and planets shrink to faint points; asteroids and moons vanish entirely — the honest emptiness of space. (They're still searchable and clickable.)",
      get on() { return state.trueScale; },
      toggle() { state.trueScale = !state.trueScale; } },
    { label: "Highlight PHAs", count: phaList.length, css: "#ff7e2a",
      desc: "Potentially Hazardous Asteroids — larger than ~140 m with orbits passing within 0.05 au (~19 lunar distances) of Earth's orbit. Lights them up in orange.",
      get on() { return state.showPHA; },
      toggle() { state.showPHA = !state.showPHA; } },
    { label: "Satellites", always: true, css: "#bfe3ff",
      count: sats.loaded ? sats.payloadCount : "—",
      desc: "Active artificial satellites — Earth's man-made moons, from CelesTrak (live TLEs, propagated with SGP4). Click Earth to fly down and see the LEO swarm, the GPS/MEO ring and the geostationary belt. Shown live: positions update in real time, independent of the time machine.",
      get on() { return sats.visible; },
      toggle() { if (!sats.loaded) { loadSatellites(); sats.visible = true; return; } sats.visible = !sats.visible; if (!sats.visible && state.selSat != null && state.selSat < sats.payloadCount) clearSelection(); } },
    { label: "Debris", always: true, css: "#ff7849",
      count: sats.loaded ? (sats.count - sats.payloadCount) : "—",
      desc: "Tracked orbital debris — fragments from satellite breakups, catalogued by CelesTrak: the Fengyun-1C anti-satellite test (2007), the Iridium-33 / Cosmos-2251 collision (2009) and the Cosmos-1408 ASAT test (2021). The most significant clouds of man-made 'space junk' still circling Earth, shown live alongside the working satellites.",
      get on() { return sats.debrisVisible; },
      toggle() { if (!sats.loaded) { loadSatellites(); sats.debrisVisible = true; return; } sats.debrisVisible = !sats.debrisVisible; if (!sats.debrisVisible && state.selSat != null && state.selSat >= sats.payloadCount) clearSelection(); } },
    { label: "Spacecraft", count: craft.length, css: CRAFT_CSS,
      desc: "Humanity's in-flight deep-space probes — Voyager 1 & 2, the Pioneers, New Horizons, Parker Solar Probe, Lucy and Psyche — with their real flight paths from JPL Horizons. The Voyagers are out past the planets, heading toward interstellar space.",
      get on() { return craftVisible; },
      toggle() { craftVisible = !craftVisible; if (!craftVisible && state.selCraft != null) clearSelection(); } },
  ];
  for (const o of overlays) {
    if (!o.always && !o.count) continue;
    const li = legendRow(o.label, o.count, o.css, o.desc, (el) => {
      o.toggle();
      el.classList.toggle("off", !o.on);
      $("fab-legend").classList.toggle("overlay-active", overlays.some((x) => x.on));
    }, o.on);
    li.classList.add("legend-toggle");
    ul.appendChild(li);
  }
  // signal on the closed legend FAB that a non-default overlay (true-scale, PHA, …) is changing the view
  $("fab-legend").classList.toggle("overlay-active", overlays.some((o) => o.on));
}

/* ============================================================
   9. Selection, info card & search
   ============================================================ */
function setRow(rowId, ddId, value) {
  const row = $(rowId);
  if (value == null) { row.hidden = true; }
  else { row.hidden = false; $(ddId).textContent = value; }
}

// tint the info card to match the selected body's legend colour
const PLANET_CSS = "#9ec5ff", SUN_CSS = "#ffd479", MOON_CSS = "#c7d2e8", SAT_CSS = "#bfe3ff", DEBRIS_CSS = "#ff7849";
function setCardAccent(css) { $("panel-info").style.setProperty("--sel", css); }

// ---- deep-link / shareable-URL state (defined here so select* can set currentRegion) ----
let currentRegion = null;      // which REGIONS key is showing (regions set no state.sel* flag)
let isApplying = false, lastHash = "", hashTimer = 0, lastHashCheck = 0;
let pendingSatSel = null, pendingSatLayers = null, skipDolly = false;

function selectObject(gi, k) {
  if (tourActive && !tourStepping) cancelTour();
  const g = groups[gi];
  const m = g.meta[k];
  state.selected = { group: gi, index: k };
  state.selPlanet = null;
  state.selMoon = null;
  state.selSun = false;
  state.selSat = null;
  state.selCraft = null;
  satLayout(false);
  // selecting a small body with its own moons (Pluto) focuses the camera on it
  if (plutoRef && gi === plutoRef[0] && k === plutoRef[1]) focusOn({ small: plutoRef });
  else flyToSmallBody(g, k);   // otherwise, always fly to frame the object
  const o = k * STRIDE, el = g.el;
  buildSelectedOrbit({
    a: el[o], e: el[o + 1], b: el[o + 2],
    Px: el[o + 6], Py: el[o + 7], Pz: el[o + 8],
    Qx: el[o + 9], Qy: el[o + 10], Qz: el[o + 11],
  });
  const hyper = m.e >= 1;
  setCardAccent(g.css);
  $("info-name").textContent = m.name;
  $("info-class").textContent =
    (m.dwarf ? "Dwarf planet · " : "") + (CLASS_NAMES[m.cls] || m.cls || "asteroid");

  // hazard / status badges
  const badges = $("info-badges");
  badges.innerHTML = "";
  if (m.dwarf) badges.insertAdjacentHTML("beforeend", `<span class="badge badge-dwarf">● Dwarf planet</span>`);
  if (m.sentry) {
    const s = m.sentry;
    const txt = s.ts != null && s.ts > 0 ? "Torino " + s.ts
      : "Palermo " + (s.ps != null ? s.ps.toFixed(1) : "?");
    const odds = s.ip > 0 ? " · 1 in " + shortNum(Math.round(1 / s.ip)) : "";
    badges.insertAdjacentHTML("beforeend", `<span class="badge badge-sentry">☢ Sentry risk · ${txt}${odds}</span>`);
  } else if (m.pha) {
    badges.insertAdjacentHTML("beforeend", `<span class="badge badge-pha">⚠ Potentially hazardous</span>`);
  }

  // orbital
  $("info-a").textContent = hyper ? "—" : m.a.toFixed(3) + " au";
  $("info-e").textContent = m.e.toFixed(3);
  $("info-i").textContent = m.i.toFixed(1) + "°";
  $("info-per").textContent = hyper ? "unbound · escaping" : formatPeriod(Math.pow(m.a, 1.5));

  // diameter: measured if we have it, else estimated from absolute magnitude
  let dTxt = "—", dLabel = "diameter";
  if (isFinite(m.diam)) dTxt = formatKm(m.diam);
  else if (isFinite(m.H)) {
    const alb = isFinite(m.albedo) && m.albedo > 0 ? m.albedo : 0.14;
    dTxt = "~" + formatKm(1329 / Math.sqrt(alb) * Math.pow(10, -m.H / 5));
    dLabel = "est. diameter";
  }
  $("info-d").textContent = dTxt;
  $("info-d-label").textContent = dLabel;

  // optional rows — hidden when the datum is missing
  const showQ = gi === GI.NEO || gi === GI.MCA || gi === GI.ISO || m.kind === "c";
  setRow("info-row-q", "info-q", showQ ? m.q.toFixed(3) + " au" : null);
  setRow("info-row-albedo", "info-albedo", isFinite(m.albedo) ? m.albedo.toFixed(2) : null);
  setRow("info-row-rot", "info-rot", isFinite(m.rot) && m.rot > 0 ? formatHours(m.rot) : null);
  setRow("info-row-spec", "info-spec", m.spec || null);
  setRow("info-row-moid", "info-moid", isFinite(m.moid)
    ? (m.moid < 0.001 ? m.moid.toFixed(4) + " au · " + (m.moid * 389.17).toFixed(2) + " LD" : m.moid.toFixed(3) + " au")
    : null);
  // honesty: how stale / well-constrained the orbit solution is
  setRow("info-row-lastobs", "info-lastobs", m.lastObs);
  setRow("info-row-arc", "info-arc", isFinite(m.arc) ? formatArc(m.arc) : null);
  const U = m.ucode == null ? NaN : parseInt(m.ucode, 10);
  const shortArc = (isFinite(U) && U >= 5) || (isFinite(m.arc) && m.arc < 30);
  const cav = $("info-caveat");
  if (shortArc) {
    const why = isFinite(m.arc) && m.arc < 30
      ? (m.arc < 1.5 ? "single-night arc" : Math.round(m.arc) + "-day arc")
      : "uncertainty code U" + U;
    cav.textContent = "Short observation arc (" + why + ") — this orbit is poorly constrained; the plotted path is approximate.";
    cav.hidden = false;
  } else cav.hidden = true;

  const link = $("info-link");
  link.href = "https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/?sstr=" + encodeURIComponent(m.name);
  link.hidden = false;

  renderPreview(m);
  $("panel-info").hidden = false;
}
/* Sun-centred fly-to for a small body: swing the camera to its side of the
   sky and frame its current solar distance so the object and its orbit
   context land together (selection must always arrive somewhere). */
function flyToSmallBody(g, k) {
  if (state.focus) { state.focus = null; retarget(); }
  const x = g.pos[k * 3], y = g.pos[k * 3 + 1], z = g.pos[k * 3 + 2];
  const r = Math.hypot(x, y, z);
  if (!(r > 0.05)) return;
  const az = Math.atan2(y, x);
  let dYaw = az - state.cam.tYaw;
  dYaw = ((dYaw + Math.PI) % TWO_PI + TWO_PI) % TWO_PI - Math.PI;
  state.cam.tYaw += dYaw;
  state.cam.tPitch = clamp(0.42 + Math.atan2(z, Math.hypot(x, y)) * 0.5, -1.3, 1.3);
  state.cam.tDist = clamp(r * 2.3, 0.6, MAX_DIST);
  state.topDown = false;
  $("fab-view").classList.remove("on");
}

/* Camera focus: re-target the orbit camera onto a body and pick a zoom that
   frames its moon system (major moons ≥ 400 km set the scale). */
function retarget() {
  state.cam.fFrom.set(state.cam.target);
  state.cam.fBlend = 0;
}
function focusOn(f) {
  state.focus = f;
  retarget();
  if (f.planet === EARTH_IDX) sats.fullUpdate = true;   // force a full SGP4 pass on arrival
  if (f.planet != null) {
    // arrive AT the planet, Eyes-style: the lit disc (and rings/satellite
    // shells) fills ~26% of the viewport; the moon system is a scroll away
    const rAU = (PLANET_FACTS[f.planet].d * 0.5) / AU_KM;
    state.cam.tDist = clamp(rAU * 14.8, 1e-4, 0.5);
    // approach from sunward so the disc arrives lit, not as a night-side crescent
    const pp = planetState[f.planet].pos;
    const sunYaw = Math.atan2(-pp[1], -pp[0]);
    let dYaw = sunYaw - state.cam.tYaw;
    dYaw = ((dYaw + Math.PI) % TWO_PI + TWO_PI) % TWO_PI - Math.PI;
    state.cam.tYaw += dYaw;
    state.cam.tPitch = 0.22;
    return;
  }
  let aMax = 0, aAny = 0;
  for (let k = 0; k < moons.count; k++) {
    const isChild = moons.parentSmall[k] && f.small && moons.parentSmall[k][0] === f.small[0] && moons.parentSmall[k][1] === f.small[1];
    if (!isChild) continue;
    aAny = Math.max(aAny, moons.meta[k].a);
    if (moons.meta[k].radius >= 400) aMax = Math.max(aMax, moons.meta[k].a);
  }
  const span = aMax || Math.min(aAny, 0.05) || 0.008;
  state.cam.tDist = clamp(span * 2.6, 2e-4, 0.5);
}
function unfocus() {
  state.focus = null;
  retarget();
  state.cam.tDist = Math.max(state.cam.tDist, 2.5);
}

function selectPlanet(i) {
  if (tourActive && !tourStepping) cancelTour();
  const ps = planetState[i];
  const facts = PLANET_FACTS[i];
  state.selPlanet = i;
  state.selected = null;
  state.selMoon = null;
  state.selSun = false;
  state.selSat = null;
  state.selCraft = null;
  satLayout(false);
  selOrbitCount = 0;
  focusOn({ planet: i });
  if (i === EARTH_IDX) loadSatellites();   // Earth's artificial moons, on arrival
  let nMoons = 0;
  for (let k = 0; k < moons.count; k++) if (moons.parentIdx[k] === i) nMoons++;
  setCardAccent(PLANET_CSS);
  $("info-name").textContent = ps.def.name;
  $("info-class").textContent = "Planet";
  const badges = $("info-badges");
  badges.innerHTML = "";
  if (nMoons) badges.insertAdjacentHTML("beforeend", `<span class="badge badge-dwarf">● ${nMoons} moon${nMoons > 1 ? "s" : ""} mapped</span>`);
  if (i === EARTH_IDX && sats.loaded) badges.insertAdjacentHTML("beforeend", `<span class="badge badge-sat">🛰 ${sats.payloadCount.toLocaleString("en-US")} satellites</span>`);
  const el = planetElements(ps.def, state.simJD, { a: 0, e: 0, i: 0, om: 0, w: 0, M: 0 });
  $("info-a").textContent = el.a.toFixed(3) + " au";
  $("info-e").textContent = el.e.toFixed(4);
  const iDeg = el.i / DEG;
  $("info-i").textContent = (Math.abs(iDeg) < 0.005 ? 0 : iDeg).toFixed(2) + "°";
  $("info-per").textContent = formatPeriod(Math.pow(el.a, 1.5));
  $("info-d-label").textContent = "diameter";
  $("info-d").textContent = facts.d.toLocaleString("en-US") + " km";
  setRow("info-row-q", "info-q", null);
  setRow("info-row-albedo", "info-albedo", null);
  setRow("info-row-rot", "info-rot", formatRotation(facts.rot));
  setRow("info-row-spec", "info-spec", null);
  setRow("info-row-moid", "info-moid", null);
  setRow("info-row-lastobs", "info-lastobs", null);
  setRow("info-row-arc", "info-arc", null);
  $("info-link").hidden = true;
  showPhotoPreview(PLANET_IMG[i][0], PLANET_IMG[i][1]);
  $("panel-info").hidden = false;
}

function selectSun() {
  if (tourActive && !tourStepping) cancelTour();
  state.selSun = true;
  state.selected = null;
  state.selPlanet = null;
  state.selMoon = null;
  state.selSat = null;
  state.selCraft = null;
  satLayout(false);
  selOrbitCount = 0;
  moonOrbitCount = 0;
  setCardAccent(SUN_CSS);
  $("info-name").textContent = "Sun";
  $("info-class").textContent = "G2V main-sequence star";
  $("info-badges").innerHTML =
    `<span class="badge badge-dwarf">● 99.86% of the system's mass</span>`;
  $("info-a").textContent = "—";
  $("info-e").textContent = "—";
  $("info-i").textContent = "—";
  $("info-per").textContent = "≈230 Myr (galactic)";
  $("info-d-label").textContent = "diameter";
  $("info-d").textContent = "1,392,700 km";
  setRow("info-row-q", "info-q", null);
  setRow("info-row-albedo", "info-albedo", null);
  setRow("info-row-rot", "info-rot", "25.05 d (equator)");
  setRow("info-row-spec", "info-spec", "G2V");
  setRow("info-row-moid", "info-moid", null);
  setRow("info-row-lastobs", "info-lastobs", null);
  setRow("info-row-arc", "info-arc", null);
  $("info-link").hidden = true;
  showPhotoPreview("sun.jpg", "NASA/SDO · extreme UV");
  $("panel-info").hidden = false;
}

/* Regions — belts & clouds. Not bodies, so they draw no orbit: clicking one and
   seeing no orbit ring (every other target shows one) is itself the point —
   these are populations/shells, not a single object on a path. Each is framed
   by how well we actually know it: belts are directly observed, the Oort cloud
   is inferred. */
const REGIONS = {
  asteroid: {
    name: "Asteroid Belt", cls: "Observed · rock & metal between Mars and Jupiter", accent: "#ffd479",
    badge: `<span class="badge" style="color:#ffd479;background:rgba(255,212,121,0.12)">● the system's densest population</span>`,
    desc: "The solar system's densest band of rock and metal — millions of bodies orbiting between Mars and Jupiter, roughly 2.1 to 3.3 au out, left over from a planet that Jupiter's gravity never let form. They run from the dwarf planet Ceres down to dust, yet their combined mass is only a few percent of the Moon's. It's a population, not a single object — every member rides its own orbit, which is why selecting it traces none.",
    facts: [["inner edge", "≈ 2.1 au"], ["outer edge", "≈ 3.3 au"], ["largest", "Ceres · 940 km"], ["total mass", "~3% of the Moon"]],
    src: "directly observed · over a million catalogued; this map plots the brightest tens of thousands",
    preview: { name: "Belt asteroid", kind: "a", spec: "S", palette: [150, 128, 98], diam: 120, albedo: 0.12, rot: NaN, img: null },
    cap: "representation · a stony main-belt asteroid (scaled to size & class)",
  },
  kuiper: {
    name: "Kuiper Belt", cls: "Observed · icy worlds beyond Neptune", accent: "#7da7ff",
    badge: `<span class="badge" style="color:#7da7ff;background:rgba(125,167,255,0.12)">● home of Pluto &amp; the dwarf planets</span>`,
    desc: "A broad ring of icy worlds beyond Neptune, from about 30 to 50 au — the source of the short-period comets and home to Pluto and most of the dwarf planets. Its bodies are ancient ice and rock, often reddened by irradiated organics (tholins). Thousands are now catalogued and one, Arrokoth, has been seen up close — unlike the Oort cloud, this belt is directly observed. It too is a population: countless separate orbits, not one path, so clicking it traces no single ring.",
    facts: [["inner edge", "≈ 30 au"], ["outer edge", "≈ 50 au"], ["largest", "Pluto · 2,377 km"], ["known members", "~3,000+"]],
    src: "directly observed · New Horizons flew through it — Pluto 2015, Arrokoth 2019",
    preview: { name: "Kuiper body", kind: "c", noComa: true, palette: [150, 116, 104], diam: 80, albedo: 0.1, rot: NaN, img: null },
    cap: "representation · a reddish icy Kuiper body (scaled to size & class)",
  },
  oort: {
    name: "Oort Cloud", cls: "Inferred · leftover icy planetesimals", accent: "#9aa7c4",
    badge: `<span class="badge" style="color:#9aa7c4;background:rgba(154,167,196,0.12)">◌ inferred · never directly observed</span>`,
    desc: "A vast spherical shell of leftover icy planetesimals — the Sun's construction debris, flung to the edge of the system by the young giant planets some 4.6 billion years ago. Each is a few kilometres of frozen volatiles — water, methane, ammonia and carbon-monoxide ice — laced with silicate and organic dust. No Oort-cloud body has ever been directly seen; its existence and extent are inferred from the long-period comets that fall in from these distances. Unlike everything else on this map it has no orbit to trace — it's a diffuse shell, not a body on a path.",
    facts: [["inner edge", "≈ 2,000 au"], ["outer edge", "≈ 100,000 au · 1.5 ly"], ["est. population", "~10¹²–10¹³ bodies"], ["typical size", "1–20 km across"]],
    src: "inferred from long-period comet orbits · never directly observed",
    preview: { name: "Oort planetesimal", kind: "c", noComa: true, palette: [120, 140, 165], diam: 8, albedo: 0.06, rot: NaN, img: null },
    cap: "representation · a dormant icy planetesimal — none has ever been imaged",
  },
  interstellar: {
    name: "Interstellar Space", cls: "Observed · the medium between the stars", accent: "#9aa7c4",
    badge: `<span class="badge" style="color:#9aa7c4;background:rgba(154,167,196,0.12)">● both Voyagers are out here</span>`,
    desc: "Where the Sun's domain finally gives out. The solar wind blows a vast bubble — the heliosphere — through the thin gas drifting between the stars; its edge, the heliopause, lies around 120 au. Past it the Sun is just the brightest star in a cold near-vacuum of hydrogen and helium laced with dust — the Local Interstellar Cloud the whole solar system is currently drifting through. Voyager 1 crossed into it in 2012 and Voyager 2 in 2018, the only craft ever to leave. The Sun's gravity reaches far beyond even here, though — the Oort cloud is still bound to it.",
    facts: [["heliopause", "≈ 120 au"], ["first crossed", "Voyager 1 · 2012"], ["medium", "Local Interstellar Cloud"], ["to nearest star", "4.24 ly"]],
    src: "heliopause distances from the Voyager interstellar mission (NASA/JPL)",
  },
  alphacen: {
    name: "Alpha Centauri", cls: "Observed · the nearest star system", accent: "#7dd3fc",
    badge: `<span class="badge" style="color:#7dd3fc;background:rgba(125,211,252,0.13)">● our nearest stellar neighbour</span>`,
    desc: "The closest star system to our own — a gravitationally bound trio. Alpha Centauri A is a near-twin of the Sun (G2V); B is a smaller orange dwarf (K1V); and faint red dwarf Proxima Centauri, a touch closer at 4.24 light-years, is the single nearest star to the Sun. Proxima holds at least one known world — the roughly Earth-mass Proxima b, in its habitable zone. This label points along the system's true direction in the sky; at Voyager's speed the trip would take about 75,000 years.",
    facts: [["distance", "4.37 ly · 276,000 au"], ["nearest star", "Proxima · 4.24 ly"], ["system", "triple · G2V·K1V·M5.5V"], ["nearest planet", "Proxima b"]],
    src: "distances & spectral types from Gaia / Hipparcos parallax",
  },
};

// A region (belt/cloud) — not a body, so it deliberately draws no orbit ring.
function selectRegion(key) {
  const R = REGIONS[key];
  if (!R) return;
  state.selected = null;
  state.selPlanet = null;
  state.selMoon = null;
  state.selSun = false;
  state.selSat = null;
  state.selCraft = null;
  selOrbitCount = 0;
  moonOrbitCount = 0;
  currentRegion = key;
  setCardAccent(R.accent);
  $("info-name").textContent = R.name;
  $("info-class").textContent = R.cls;
  $("info-badges").innerHTML = R.badge || "";
  $("info-link").hidden = true;
  panelMode("region");
  $("region-desc").textContent = R.desc;
  $("region-grid").innerHTML = R.facts.map(([k, v]) => `<div><dt>${k}</dt><dd>${v}</dd></div>`).join("");
  $("region-src").textContent = R.src;
  // empty space / a distant star has nothing honest to depict — skip the preview
  if (R.preview) {
    $("info-preview").hidden = false;
    renderPreview(R.preview);
    $("preview-cap").textContent = R.cap;
  } else {
    stopPreview();
    $("info-preview").hidden = true;
  }
  $("panel-info").hidden = false;
}

// satellites use a dedicated panel layout; toggle it vs the body grid/preview
// swap the card body between the default grid+preview, the satellite block,
// the spacecraft block, and the region block. mode: "body" | "sat" | "craft" | "region"
function panelMode(mode) {
  $("info-grid").hidden = mode !== "body";
  $("info-preview").hidden = !(mode === "body" || mode === "region");
  $("info-sat").hidden = mode !== "sat";
  $("info-craft").hidden = mode !== "craft";
  $("info-region").hidden = mode !== "region";
  $("info-caveat").hidden = true;   // only selectObject re-shows it (short-arc orbits)
}
const satLayout = (on) => panelMode(on ? "sat" : "body");

function selectSatellite(k) {
  if (tourActive && !tourStepping) cancelTour();
  state.selSat = k;
  state.selected = null;
  state.selPlanet = null;
  state.selMoon = null;
  state.selSun = false;
  state.selCraft = null;
  selOrbitCount = 0;
  moonOrbitCount = 0;
  stopPreview();
  const f = satFacts(k);
  const isDebris = k >= sats.payloadCount;
  setCardAccent(isDebris ? DEBRIS_CSS : SAT_CSS);
  $("info-name").textContent = sats.names[k];
  $("info-class").textContent = (isDebris ? "Orbital debris · NORAD " : "Artificial satellite · NORAD ") + f.norad;
  $("info-badges").innerHTML = isDebris
    ? `<span class="badge" style="color:#ff7849;background:rgba(255,120,73,0.13)">⚠ debris · ${f.regime}</span>`
    : `<span class="badge badge-sat">🛰 ${f.regime}</span>`;
  $("info-link").hidden = true;
  panelMode("sat");
  $("sat-regime").textContent = f.regime;
  $("sat-alt").textContent = f.peri < -50 ? "—"
    : Math.round(Math.max(f.peri, 0)).toLocaleString("en-US") + " × " + Math.round(f.apo).toLocaleString("en-US") + " km";
  $("sat-period").textContent = f.periodMin < 1440 ? f.periodMin.toFixed(0) + " min" : (f.periodMin / 1440).toFixed(2) + " days";
  $("sat-inc").textContent = f.inc.toFixed(1) + "°";
  $("panel-info").hidden = false;
}

function selectCraft(i) {
  if (tourActive && !tourStepping) cancelTour();
  const c = craft[i];
  if (!c) return;
  state.selCraft = i;
  state.selected = null;
  state.selPlanet = null;
  state.selMoon = null;
  state.selSun = false;
  state.selSat = null;
  selOrbitCount = 0;
  moonOrbitCount = 0;
  stopPreview();
  const info = c.info;
  setCardAccent(CRAFT_CSS);
  $("info-name").textContent = c.name;
  $("info-class").textContent = "Spacecraft" + (info.op ? " · " + info.op : "");
  $("info-badges").innerHTML = `<span class="badge badge-craft">🚀 In-flight probe</span>`;
  $("info-link").hidden = true;
  panelMode("craft");
  // most craft cards have no photo; Starman (the Roadster) carries a real SpaceX shot
  if (info.img) { $("info-preview").hidden = false; showPhotoPreview(info.img[0], info.img[1]); }
  $("craft-note").textContent = info.note || "";
  $("craft-launch").textContent = info.launch || "—";
  $("craft-status").textContent = info.status === "silent" ? "silent · still coasting" : "active";
  const f = craftFacts(c);
  $("craft-r").textContent = c.inRange ? f.r.toFixed(2) + " au" : "outside data range";
  $("craft-v").textContent = c.inRange && isFinite(f.v) ? f.v.toFixed(1) + " km/s" : "—";
  $("panel-info").hidden = false;
}

function selectMoon(k) {
  if (tourActive && !tourStepping) cancelTour();
  const m = moons.meta[k];
  state.selMoon = k;
  state.selected = null;
  state.selPlanet = null;
  state.selSun = false;
  state.selSat = null;
  state.selCraft = null;
  satLayout(false);
  selOrbitCount = 0;
  buildMoonOrbit(k);
  setCardAccent(MOON_CSS);
  $("info-name").textContent = m.name;
  $("info-class").textContent = "Moon of " + m.parentName;
  $("info-badges").innerHTML = "";
  $("info-a").textContent = Math.round(m.a * 149597870.7).toLocaleString("en-US") + " km";
  $("info-e").textContent = m.e.toFixed(4);
  $("info-i").textContent = m.i.toFixed(1) + "°";
  $("info-per").textContent = formatPeriod(360 / m.n / 365.25);
  $("info-d-label").textContent = "diameter";
  $("info-d").textContent = m.radius ? formatKm(m.radius * 2) : "—";
  setRow("info-row-q", "info-q", null);
  setRow("info-row-albedo", "info-albedo", null);
  setRow("info-row-rot", "info-rot", null);
  setRow("info-row-spec", "info-spec", null);
  setRow("info-row-moid", "info-moid", null);
  setRow("info-row-lastobs", "info-lastobs", null);
  setRow("info-row-arc", "info-arc", null);
  $("info-link").hidden = true;
  if (m.img) showPhotoPreview(m.img[0], m.img[1]);
  else renderPreview({ name: m.name, kind: "a", spec: "", dwarf: true, diam: m.radius ? m.radius * 2 : NaN, albedo: NaN, rot: NaN, img: null });
  $("panel-info").hidden = false;
}

/* show a known photo in the preview frame (planets & major moons) */
function showPhotoPreview(file, credit) {
  stopPreview();
  $("preview-img").src = "img/bodies/" + file;
  $("preview-img").hidden = false;
  $("preview-canvas").hidden = true;
  $("preview-cap").textContent = "📷 real image · " + credit;
  $("preview-cap").classList.add("real");
}
function formatRotation(h) {
  const abs = Math.abs(h);
  const txt = abs < 48 ? abs.toFixed(2) + " h" : (abs / 24).toFixed(1) + " d";
  return h < 0 ? txt + " (retrograde)" : txt;
}

function clearSelection() {
  state.selected = null;
  state.selPlanet = null;
  state.selMoon = null;
  state.selSun = false;
  state.selSat = null;
  state.selCraft = null;
  currentRegion = null;
  satLayout(false);
  selOrbitCount = 0;
  moonOrbitCount = 0;
  stopPreview();
  $("panel-info").hidden = true;
}

/* ============================================================
   Deep-link / shareable URLs — encode the live view in the hash
   ============================================================ */
const GROUP_KEY = GROUPS.map((g) => g.key);     // group index → layer key

function selectionKey() {
  if (state.selSun) return "sun";
  if (state.selPlanet != null) return "p" + state.selPlanet;
  if (state.selMoon != null) return "m" + state.selMoon;
  if (state.selSat != null) return "sat" + state.selSat;
  if (state.selCraft != null) return "c" + state.selCraft;
  if (state.selected) return "o" + encodeURIComponent(groups[state.selected.group].meta[state.selected.index].name);
  if (currentRegion) return "r" + currentRegion;
  return "";
}
function encodeState() {
  const c = state.cam, p = ["v=1"];
  const sel = selectionKey(); if (sel) p.push("sel=" + sel);
  p.push("cam=" + c.tYaw.toFixed(3) + "," + c.tPitch.toFixed(3) + "," + (+c.tDist.toPrecision(5)));
  if (state.focus) p.push("foc=" + (state.focus.planet != null
    ? "p" + state.focus.planet : "s" + state.focus.small[0] + "." + state.focus.small[1]));
  if (!state.playing) p.push("jd=" + state.simJD.toFixed(4));   // capture a frozen moment; live links stay live
  p.push("t=" + (state.dps * state.dir));
  if (!state.playing) p.push("pl=0");
  const off = [];
  if (!state.showSun) off.push("sun");
  if (!state.showPlanets) off.push("pla");
  if (!moons.visible) off.push("moon");
  if (!OORT.visible) off.push("oort");
  groups.forEach((g, gi) => { if (g.count && !g.visible) off.push(GROUP_KEY[gi]); });
  if (off.length) p.push("off=" + off.join(","));
  const on = [];
  if (state.showPHA) on.push("pha");
  if (sats.visible) on.push("sat");
  if (sats.debrisVisible) on.push("deb");
  if (craftVisible) on.push("craft");
  if (state.trueScale) on.push("ts");
  if (on.length) p.push("on=" + on.join(","));
  return "#" + p.join("&");
}
function resolveSelection(tok) {
  if (tok === "sun") return selectSun();
  if (tok.slice(0, 3) === "sat") { const k = +tok.slice(3); focusOn({ planet: EARTH_IDX });
    if (sats.loaded) selectSatellite(k); else { loadSatellites(); pendingSatSel = k; } return; }
  if (tok[0] === "p") return selectPlanet(+tok.slice(1));
  if (tok[0] === "m") { const k = +tok.slice(1), pi = moons.parentIdx[k];
    if (pi != null && pi >= 0) focusOn({ planet: pi }); return selectMoon(k); }
  if (tok[0] === "c") { craftVisible = true; return selectCraft(+tok.slice(1)); }
  if (tok[0] === "r") return selectRegion(tok.slice(1));
  if (tok[0] === "o") { const name = decodeURIComponent(tok.slice(1));
    for (let gi = 0; gi < groups.length; gi++) { const g = groups[gi];
      for (let k = 0; k < g.count; k++) if (g.meta[k].name === name) return selectObject(gi, k); } }
}
function applyState(hash) {
  if (!hash || hash.length < 2) return;
  try {
    isApplying = true;
    const map = {};
    hash.replace(/^#/, "").split("&").forEach((kv) => { const i = kv.indexOf("="); if (i > 0) map[kv.slice(0, i)] = kv.slice(i + 1); });
    if (map.v !== "1") return;
    const off = (map.off || "").split(",").filter(Boolean);
    const on = (map.on || "").split(",").filter(Boolean);
    state.showSun = !off.includes("sun");
    state.showPlanets = !off.includes("pla");
    moons.visible = !off.includes("moon");
    OORT.visible = !off.includes("oort");
    groups.forEach((g, gi) => { g.visible = !off.includes(GROUP_KEY[gi]); });
    state.showPHA = on.includes("pha");
    craftVisible = on.includes("craft");
    state.trueScale = on.includes("ts");
    if (on.includes("sat") || on.includes("deb")) { loadSatellites(); pendingSatLayers = { wantSat: on.includes("sat"), wantDeb: on.includes("deb") }; }
    renderLegend();
    if (map.t != null && isFinite(+map.t)) {
      const t = +map.t; state.dps = Math.abs(t) || 30; state.dir = t < 0 ? -1 : 1;
      for (const x of $("speed-chips").children) x.classList.toggle("on", +x.dataset.dps === state.dps);
    }
    state.playing = map.pl !== "0"; syncPlayBtn();
    state.simJD = (map.jd != null && isFinite(+map.jd)) ? +map.jd : jdNow();
    state.needFullUpdate = true; buildPlanetOrbits();
    if (map.sel) resolveSelection(map.sel);
    if (map.foc && !state.focus) {
      if (map.foc[0] === "p") focusOn({ planet: +map.foc.slice(1) });
      else if (map.foc[0] === "s") { const a = map.foc.slice(1).split("."); focusOn({ small: [+a[0], +a[1]] }); }
    }
    if (map.cam) {
      const a = map.cam.split(",").map(Number);
      if (a.length === 3 && a.every(isFinite)) {
        const c = state.cam; c.tYaw = c.yaw = a[0]; c.tPitch = c.pitch = a[1]; c.tDist = c.dist = a[2];
        c.fBlend = 1; state.camEaseRate = 9; skipDolly = true;
      }
    }
  } catch (e) { console.warn("deep-link restore failed:", e); }
  finally { isApplying = false; }
}
function scheduleHashWrite(h) {
  clearTimeout(hashTimer);
  hashTimer = setTimeout(() => { try { history.replaceState(null, "", h); } catch (_) { /* ignore */ } }, 250);
}
(function wireShare() {
  const b = $("btn-share");
  if (b) b.addEventListener("click", async () => {
    const url = location.origin + location.pathname + encodeState();
    try { await navigator.clipboard.writeText(url); } catch (_) { prompt("Copy this link to share the view:", url); }
    b.classList.add("copied"); b.textContent = "Copied";
    setTimeout(() => { b.classList.remove("copied"); b.textContent = "Share"; }, 1600);
  });
  window.addEventListener("hashchange", () => { if (!isApplying) { applyState(location.hash); lastHash = encodeState(); } });
})();

function formatArc(d) {
  if (d < 1.5) return "single night";
  if (d < 90) return Math.round(d) + " days";
  const yr = d / 365.25;
  return yr < 10 ? yr.toFixed(1) + " yr" : Math.round(yr) + " yr";
}
function formatPeriod(yr) {
  const d = yr * 365.25;
  if (d < 10) return d.toFixed(2) + " days";
  if (yr >= 1e6) return (yr / 1e6).toFixed(1) + " Myr";
  if (yr >= 1e4) return Math.round(yr / 1e3) + " kyr";
  if (yr >= 1e3) return Math.round(yr).toLocaleString("en-US") + " yr";
  return yr < 1.5 ? d.toFixed(1) + " days" : yr.toFixed(yr < 10 ? 2 : 1) + " yr";
}
function formatKm(km) {
  return km < 1 ? Math.round(km * 1000) + " m" : km < 10 ? km.toFixed(1) + " km" : Math.round(km) + " km";
}
function formatHours(h) {
  if (h < 1) return Math.round(h * 60) + " min";
  if (h < 48) return h.toFixed(h < 10 ? 2 : 1) + " h";
  return (h / 24).toFixed(1) + " d";
}

/* ============================================================
   9b. Object preview — real photo where one exists, else a
   procedural representation driven by the object's real
   diameter, albedo, spectral class and rotation.
   ============================================================ */
let previewRAF = 0;
function stopPreview() { if (previewRAF) { cancelAnimationFrame(previewRAF); previewRAF = 0; } }

function renderPreview(m) {
  stopPreview();
  const cv = $("preview-canvas");
  const img = $("preview-img");
  const cap = $("preview-cap");
  const key = m.img;
  if (key && BODY_IMAGES[key]) {
    img.src = BODY_IMAGES[key].file;
    img.hidden = false; cv.hidden = true;
    cap.textContent = "📷 real image · " + BODY_IMAGES[key].credit;
    cap.classList.add("real");
    return;
  }
  img.hidden = true; img.removeAttribute("src"); cv.hidden = false;
  cap.classList.remove("real");
  cap.textContent = "representation · scaled to size, albedo & class — not a photograph";
  drawProcedural(cv, m);
}

function hashSeed(str) {
  let h = 2166136261;
  for (let i = 0; i < str.length; i++) { h ^= str.charCodeAt(i); h = Math.imul(h, 16777619); }
  return h >>> 0;
}
function mulberry32(a) {
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
const _rgb = (c, a) => `rgba(${c[0] | 0},${c[1] | 0},${c[2] | 0},${a})`;
const _cl = (c) => c.map((v) => (v < 0 ? 0 : v > 255 ? 255 : v));
function bodyPalette(m) {
  if (m.palette) return m.palette;                            // explicit override (region previews)
  const s = (m.spec || "").toUpperCase();
  if (m.kind === "c") return [120, 140, 165];                 // icy comet nucleus
  if (/^[CBFGPDT]/.test(s)) return [86, 74, 64];              // carbonaceous / primitive — dark
  if (/^[SQAVR]/.test(s)) return [165, 130, 96];              // stony / silicaceous — tan
  if (/^[MXE]/.test(s)) return [150, 150, 162];               // metallic — grey
  if (m.dwarf) return [184, 174, 162];                        // icy dwarf — pale
  return [126, 118, 110];                                     // unknown — neutral
}

function drawProcedural(cv, m) {
  const dprL = Math.min(window.devicePixelRatio || 1, 2);
  const W = cv.clientWidth || 300, H = cv.clientHeight || 170;
  cv.width = W * dprL; cv.height = H * dprL;
  const ctx = cv.getContext("2d");
  ctx.setTransform(dprL, 0, 0, dprL, 0, 0);
  const cx = W / 2, cy = H / 2;
  const rnd = mulberry32(hashSeed(m.name));
  const base = bodyPalette(m);
  const isComet = m.kind === "c";
  const dk = isFinite(m.diam) ? m.diam : 5;
  const round = isComet ? 0.65 : Math.min(0.85, 0.25 + Math.log10(Math.max(dk, 1)) / 4);
  const R = Math.min(W, H) * 0.30;
  const albk = isFinite(m.albedo) && m.albedo > 0 ? Math.min(1.4, 0.5 + m.albedo * 2.2) : 0.85;
  const spin = isFinite(m.rot) && m.rot > 0 ? Math.min(0.6, 2.0 / m.rot) : 0.22;

  const harm = [];
  for (let i = 0; i < 4; i++) harm.push({ k: i + 2, amp: (1 - round) * (0.10 + rnd() * 0.13) / (i + 1), ph: rnd() * TWO_PI });
  const silR = (a) => { let r = R; for (const h of harm) r += R * h.amp * Math.sin(h.k * a + h.ph); return r; };

  const craters = [];
  const nC = isComet ? 4 : 8 + Math.floor(rnd() * 8);
  for (let i = 0; i < nC; i++) craters.push({ lon: rnd() * TWO_PI, lat: (rnd() - 0.5) * 1.4, r: (0.06 + rnd() * 0.13) * R, d: 0.45 + rnd() * 0.3 });
  const speck = [];
  for (let i = 0; i < 55; i++) speck.push({ lon: rnd() * TWO_PI, lat: (rnd() - 0.5) * Math.PI * 0.9, r: 0.5 + rnd() * 1.5, t: rnd() * 0.5 - 0.25 });

  const project = (lon, lat, r) => {
    const cl = Math.cos(lat), front = Math.cos(lon) * cl;
    if (front <= 0.03) return null;
    return { x: Math.sin(lon) * cl * r, y: -Math.sin(lat) * r * 0.96, f: front };
  };
  const tracePath = () => {
    ctx.beginPath();
    for (let s = 0; s <= 64; s++) {
      const a = (s / 64) * TWO_PI, r = silR(a);
      const x = cx + Math.cos(a) * r, y = cy + Math.sin(a) * r * 0.96;
      s ? ctx.lineTo(x, y) : ctx.moveTo(x, y);
    }
    ctx.closePath();
  };

  let t0 = null;
  function frame(ts) {
    if (t0 == null) t0 = ts;
    const phase = ((ts - t0) / 1000) * spin;
    ctx.clearRect(0, 0, W, H);
    if (isComet && !m.noComa) drawComa(ctx, cx, cy, R, phase, albk);

    ctx.save();
    tracePath();
    ctx.clip();
    ctx.fillStyle = _rgb(base, 1);
    ctx.fillRect(cx - R * 1.7, cy - R * 1.7, R * 3.4, R * 3.4);
    for (const sp of speck) {
      const p = project(sp.lon + phase, sp.lat, R);
      if (!p) continue;
      ctx.fillStyle = _rgb(_cl(base.map((v) => v * (1 + sp.t))), 0.5);
      ctx.beginPath(); ctx.arc(cx + p.x, cy + p.y, sp.r * p.f, 0, TWO_PI); ctx.fill();
    }
    for (const cr of craters) {
      const p = project(cr.lon + phase, cr.lat, R);
      if (!p) continue;
      const rr = cr.r * p.f; if (rr < 1) continue;
      const g = ctx.createRadialGradient(cx + p.x - rr * 0.3, cy + p.y - rr * 0.3, rr * 0.1, cx + p.x, cy + p.y, rr);
      g.addColorStop(0, _rgb(_cl(base.map((v) => v * cr.d)), 0.9));
      g.addColorStop(0.7, _rgb(_cl(base.map((v) => v * (cr.d + 0.25))), 0.45));
      g.addColorStop(1, _rgb(base, 0));
      ctx.fillStyle = g;
      ctx.beginPath(); ctx.arc(cx + p.x, cy + p.y, rr, 0, TWO_PI); ctx.fill();
    }
    // light from upper-left, shadow opposite
    const lg = ctx.createRadialGradient(cx - 0.55 * R, cy - 0.5 * R, R * 0.1, cx, cy, R * 1.5);
    lg.addColorStop(0, `rgba(255,250,236,${0.45 * albk})`);
    lg.addColorStop(0.45, "rgba(255,250,236,0)");
    lg.addColorStop(1, "rgba(0,0,8,0.74)");
    ctx.fillStyle = lg;
    ctx.fillRect(cx - R * 1.8, cy - R * 1.8, R * 3.6, R * 3.6);
    ctx.restore();

    ctx.save();
    tracePath();
    ctx.lineWidth = 1.4;
    ctx.strokeStyle = "rgba(150,180,235,0.16)";
    ctx.stroke();
    ctx.restore();

    previewRAF = requestAnimationFrame(frame);
  }
  previewRAF = requestAnimationFrame(frame);
}
function drawComa(ctx, cx, cy, R, phase, albk) {
  const flick = 0.85 + 0.15 * Math.sin(phase * 4);
  const tg = ctx.createLinearGradient(cx, cy, cx + R * 5, cy - R * 1.0);
  tg.addColorStop(0, `rgba(150,232,255,${0.32 * flick})`);
  tg.addColorStop(1, "rgba(120,200,255,0)");
  ctx.fillStyle = tg;
  ctx.beginPath();
  ctx.moveTo(cx, cy - R * 0.6); ctx.lineTo(cx + R * 5, cy - R * 1.5);
  ctx.lineTo(cx + R * 5, cy + R * 0.3); ctx.lineTo(cx, cy + R * 0.5);
  ctx.closePath(); ctx.fill();
  const cg = ctx.createRadialGradient(cx, cy, R * 0.2, cx, cy, R * 2.0);
  cg.addColorStop(0, `rgba(205,246,255,${0.5 * flick})`);
  cg.addColorStop(1, "rgba(150,220,255,0)");
  ctx.fillStyle = cg;
  ctx.beginPath(); ctx.arc(cx, cy, R * 2.0, 0, TWO_PI); ctx.fill();
}

function pickAt(x, y) {
  let best = null, bestD = 22 * 22;
  for (let gi = 0; gi < groups.length; gi++) {
    const g = groups[gi];
    if (!g.visible) continue;
    for (let k = 0; k < g.count; k++) {
      const p = project(g.pos[k * 3], g.pos[k * 3 + 1], g.pos[k * 3 + 2]);
      if (!p) continue;
      const dx = p[0] - x, dy = p[1] - y;
      const d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; best = [gi, k]; }
    }
  }
  // the Sun
  let bestSun = false;
  if (state.showSun) {
    const p = project(0, 0, 0);
    if (p) {
      const dx = p[0] - x, dy = p[1] - y, d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; bestSun = true; best = null; }
    }
  }
  // planets
  let bestPlanet = -1;
  if (state.showPlanets) {
    for (let i = 0; i < planetState.length; i++) {
      const ps = planetState[i];
      const p = project(ps.pos[0], ps.pos[1], ps.pos[2]);
      if (!p) continue;
      const dx = p[0] - x, dy = p[1] - y, d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; bestPlanet = i; best = null; }
    }
  }
  // moons — only pickable once visually separated from their parent,
  // so clicking an unzoomed planet never grabs an invisible moon
  let bestMoon = -1;
  if (moons.visible) {
    for (let k = 0; k < moons.count; k++) {
      const p = project(moons.pos[k * 3], moons.pos[k * 3 + 1], moons.pos[k * 3 + 2]);
      if (!p) continue;
      const dx = p[0] - x, dy = p[1] - y, d = dx * dx + dy * dy;
      if (d >= bestD) continue;
      if (moonParentPos(k, _mpp)) {
        const pp = project(_mpp[0], _mpp[1], _mpp[2]);
        if (pp && Math.hypot(p[0] - pp[0], p[1] - pp[1]) < 10) continue;
      }
      bestD = d; bestMoon = k; bestPlanet = -1; best = null;
    }
  }
  // satellites + debris — only in the Earth-focus view, where they're the whole point
  let bestSat = -1;
  if (sats.loaded && (sats.visible || sats.debrisVisible) && state.focus && state.focus.planet === EARTH_IDX && satEpochFade() > 0.01) {
    for (let k = 0; k < sats.count; k++) {
      if (sats.pos[k * 3] > 1e8) continue;
      if (!(k < sats.payloadCount ? sats.visible : sats.debrisVisible)) continue;  // skip hidden type
      const p = project(sats.pos[k * 3], sats.pos[k * 3 + 1], sats.pos[k * 3 + 2]);
      if (!p) continue;
      const dx = p[0] - x, dy = p[1] - y, d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; bestSat = k; bestMoon = -1; bestPlanet = -1; best = null; }
    }
  }
  // spacecraft (in-flight probes)
  let bestCraft = -1;
  if (craftVisible) {
    for (let i = 0; i < craft.length; i++) {
      const c = craft[i];
      if (!c.inRange) continue;
      const p = project(c.pos[0], c.pos[1], c.pos[2]);
      if (!p) continue;
      const dx = p[0] - x, dy = p[1] - y, d = dx * dx + dy * dy;
      if (d < bestD) { bestD = d; bestCraft = i; bestSat = -1; bestMoon = -1; bestPlanet = -1; best = null; }
    }
  }
  if (bestCraft >= 0) selectCraft(bestCraft);
  else if (bestSat >= 0) selectSatellite(bestSat);
  else if (bestMoon >= 0) selectMoon(bestMoon);
  else if (bestPlanet >= 0) selectPlanet(bestPlanet);
  else if (bestSun) selectSun();
  else if (best) selectObject(best[0], best[1]);
  else clearSelection();
}

/* ---- search ---- */
const searchEl = $("search");
const resultsEl = $("search-results");
/* rank: exact name > starts-with > word-start > substring, so "iss" finds
   the ISS before it finds Lar-iss-a; ties keep category order */
function matchScore(name, q) {
  const n = name.toLowerCase();
  if (n === q) return 0;
  if (n.startsWith(q)) return 1;
  const w = n.indexOf(q);
  if (w < 0) return Infinity;
  const prev = n[w - 1];
  return prev === " " || prev === "(" || prev === "/" || prev === "-" ? 2 : 3;
}
function runSearch(q) {
  q = q.trim().toLowerCase();
  resultsEl.innerHTML = "";
  if (q.length < 2) { resultsEl.hidden = true; return; }
  if (!sats.loaded) loadSatellites();   // so live satellites are findable
  const hits = [];   // { name, cls, select, s, ord }
  let ord = 0;
  const add = (name, matchName, cls, select) => {
    const s = matchScore(matchName, q);
    if (s < Infinity) hits.push({ name, cls, select, s, ord: ord++ });
  };
  add("Sun", "sun", "STAR", selectSun);
  for (let i = 0; i < PLANETS.length; i++) {
    const ii = i;
    add(PLANETS[i].name, PLANETS[i].name, "PLANET", () => selectPlanet(ii));
  }
  let taken = 0;
  for (let k = 0; k < moons.count && taken < 14; k++) {
    if (matchScore(moons.meta[k].name, q) < Infinity) {
      const kk = k;
      add(moons.meta[k].name + " · " + moons.meta[k].parentName, moons.meta[k].name, "MOON", () => selectMoon(kk));
      taken++;
    }
  }
  for (let i = 0; i < craft.length; i++) {
    const ii = i;
    add(craft[i].name, craft[i].name, "CRAFT", () => selectCraft(ii));
  }
  if (sats.loaded) {
    taken = 0;
    for (let k = 0; k < sats.count && taken < 12; k++) {
      if (matchScore(sats.names[k], q) < Infinity) {
        const kk = k;
        add(sats.names[k], sats.names[k], "SAT", () => { selectPlanet(EARTH_IDX); selectSatellite(kk); });
        taken++;
      }
    }
  }
  taken = 0;
  outer:
  for (let gi = 0; gi < groups.length; gi++) {
    const g = groups[gi];
    for (let k = 0; k < g.count; k++) {
      if (matchScore(g.meta[k].name, q) < Infinity) {
        const ggi = gi, kk = k;
        add(g.meta[k].name, g.meta[k].name, g.meta[k].cls, () => selectObject(ggi, kk));
        if (++taken >= 30) break outer;
      }
    }
  }
  hits.sort((a, b) => a.s - b.s || a.ord - b.ord);
  if (!hits.length) { resultsEl.hidden = true; return; }
  for (const h of hits.slice(0, 8)) {
    const li = document.createElement("li");
    const nameSpan = document.createElement("span");
    nameSpan.textContent = h.name;
    const clsSpan = document.createElement("span");
    clsSpan.className = "sr-class";
    clsSpan.textContent = h.cls;
    li.append(nameSpan, clsSpan);
    li.addEventListener("pointerdown", (ev) => {
      ev.preventDefault();
      h.select();
      resultsEl.hidden = true;
      searchEl.blur();
    });
    resultsEl.appendChild(li);
  }
  resultsEl.hidden = false;
}
searchEl.addEventListener("input", () => runSearch(searchEl.value));
searchEl.addEventListener("blur", () => setTimeout(() => (resultsEl.hidden = true), 150));
searchEl.addEventListener("keydown", (ev) => {
  if (ev.key === "Enter") {
    const first = resultsEl.querySelector("li");
    if (first) first.dispatchEvent(new PointerEvent("pointerdown"));
  }
  if (ev.key === "Escape") searchEl.blur();
});
window.addEventListener("keydown", (ev) => {
  if (ev.key === "/" && document.activeElement !== searchEl &&
      !/^(INPUT|TEXTAREA)$/.test(document.activeElement.tagName)) {
    ev.preventDefault();
    searchEl.focus();
  }
});

/* ============================================================
   10. Camera input — mouse / touch / wheel
   ============================================================ */
const pointers = new Map();
let dragDist = 0, pinchD0 = 0, distAtPinch = 0;
function dismissHint() { $("hint").classList.add("gone"); }

canvas.addEventListener("pointerdown", (ev) => {
  if (tourActive) cancelTour();            // any scene interaction ends the tour
  try { canvas.setPointerCapture(ev.pointerId); } catch (_) { /* capture is best-effort */ }
  pointers.set(ev.pointerId, { x: ev.clientX, y: ev.clientY });
  dragDist = 0;
  state.camEaseRate = 9;            // grabbing the camera cancels the launch dolly
  if (pointers.size === 2) {
    const [p1, p2] = [...pointers.values()];
    pinchD0 = Math.hypot(p1.x - p2.x, p1.y - p2.y);
    distAtPinch = state.cam.tDist;
  }
});
canvas.addEventListener("pointermove", (ev) => {
  const p = pointers.get(ev.pointerId);
  if (!p) return;
  const dx = ev.clientX - p.x, dy = ev.clientY - p.y;
  dragDist += Math.abs(dx) + Math.abs(dy);
  p.x = ev.clientX; p.y = ev.clientY;
  if (pointers.size === 1) {
    state.cam.tYaw -= dx * 0.0052;
    state.cam.tPitch = clamp(state.cam.tPitch + dy * 0.0045, -1.45, 1.45);
    if (state.topDown && dragDist > 12) setTopDown(false, true);
    dismissHint();
  } else if (pointers.size === 2) {
    const [p1, p2] = [...pointers.values()];
    const d = Math.hypot(p1.x - p2.x, p1.y - p2.y);
    if (pinchD0 > 0 && d > 0) {
      const prev = state.cam.tDist;
      state.cam.tDist = clamp(distAtPinch * (pinchD0 / d), minDist(), MAX_DIST);
      farOrientIfCrossed(prev);
    }
    dismissHint();
  }
});
function endPointer(ev) {
  if (!pointers.has(ev.pointerId)) return;   // orphan release (started elsewhere)
  const wasDrag = dragDist > 7 || pointers.size > 1;
  pointers.delete(ev.pointerId);
  if (!wasDrag && ev.type === "pointerup" && ev.target === canvas) {
    pickAt(ev.clientX, ev.clientY);
  }
}
canvas.addEventListener("pointerup", endPointer);
canvas.addEventListener("pointercancel", endPointer);
canvas.addEventListener("wheel", (ev) => {
  ev.preventDefault();
  state.camEaseRate = 9;
  // bigger strides once past the planets, or the Oort journey takes 80 notches
  const accel = state.cam.tDist > 120 ? 2.2 : 1;
  const prev = state.cam.tDist;
  state.cam.tDist = clamp(state.cam.tDist * Math.exp(ev.deltaY * 0.0011 * accel), minDist(), MAX_DIST);
  farOrientIfCrossed(prev);
  dismissHint();
}, { passive: false });
function resetToSun() {
  state.focus = null;
  retarget();
  Object.assign(state.cam, { tYaw: -1.1, tPitch: 0.55, tDist: 7.8 });
  setTopDown(false, true);
}
canvas.addEventListener("dblclick", resetToSun);
function clamp(v, lo, hi) { return v < lo ? lo : v > hi ? hi : v; }
function minDist() { return state.focus ? 5e-5 : 0.18; }
// snap to the composed interstellar view the first time you cross out past the
// Oort cloud; a later drag overrides it, and zooming back in re-arms the snap
function farOrientIfCrossed(prev) {
  if (prev <= FAR_ORIENT_DIST && state.cam.tDist > FAR_ORIENT_DIST) {
    state.cam.tYaw = FAR_YAW;
    state.cam.tPitch = FAR_PITCH;
  }
}

/* ============================================================
   11. UI wiring — time machine, panels, view
   ============================================================ */
const btnPlay = $("btn-play");
function syncPlayBtn() {
  btnPlay.textContent = state.playing ? "❚❚" : "▶";
  btnPlay.classList.toggle("on", state.playing);
}
btnPlay.addEventListener("click", () => { state.playing = !state.playing; syncPlayBtn(); });
$("btn-rev").addEventListener("click", (ev) => {
  state.dir *= -1;
  ev.currentTarget.classList.toggle("on", state.dir < 0);
});
$("btn-now").addEventListener("click", () => {
  state.simJD = jdNow();
  state.needFullUpdate = true;
  buildPlanetOrbits();
});
$("speed-chips").addEventListener("click", (ev) => {
  const b = ev.target.closest("button[data-dps]");
  if (!b) return;
  state.dps = b.dataset.dps === "real" ? 1 / 86400 : +b.dataset.dps;
  if (b.dataset.dps === "real") { state.simJD = jdNow(); state.needFullUpdate = true; }  // "live" snaps back to the actual present
  if (!state.playing) { state.playing = true; syncPlayBtn(); }
  for (const x of $("speed-chips").children) x.classList.toggle("on", x === b);
});
syncPlayBtn();

/* ---- guided tour: set camera targets + selection per step; the render-loop
   easing does the flying. A generation counter makes stale timers no-ops. ---- */
const TOUR = [
  { label: "The Sun · 99.86% of the system's mass", dwell: 5500, go() { resetToSun(); selectSun(); } },
  { label: "The asteroid belt · millions of worlds", dwell: 5500,
    go() { state.focus = null; retarget(); Object.assign(state.cam, { tYaw: -1.1, tPitch: 0.45, tDist: 8 }); selectRegion("asteroid"); } },
  { label: "Jupiter · king of the planets", dwell: 6500, go() { selectPlanet(4); } },
  { label: "Saturn · the rings", dwell: 7000, go() { selectPlanet(5); state.cam.tDist = 0.004; } },
  { label: "Pluto · the Kuiper frontier", dwell: 6500,
    go() { if (plutoRef) selectObject(plutoRef[0], plutoRef[1]);
      else { state.focus = null; retarget(); Object.assign(state.cam, { tYaw: -1.1, tPitch: 0.5, tDist: 42 }); selectRegion("kuiper"); } } },
  { label: "Voyager 1 · humanity's farthest machine", dwell: 6500,
    go() { craftVisible = true; renderLegend(); state.focus = null; retarget();
      Object.assign(state.cam, { tYaw: FAR_YAW, tPitch: 0.5, tDist: 200 }); if (craft[0]) selectCraft(0); } },
  { label: "The Oort cloud · the long fall outward", dwell: 9500,
    go() { state.focus = null; retarget(); Object.assign(state.cam, { tYaw: -1.1, tPitch: 0.5, tDist: 3000 }); selectRegion("oort"); } },
  { label: "Alpha Centauri · our nearest neighbour", dwell: 8000,
    go() { state.focus = null; retarget(); Object.assign(state.cam, { tYaw: FAR_YAW, tPitch: FAR_PITCH, tDist: 80000 }); selectRegion("alphacen"); } },
];
let tourActive = false, tourGen = 0, tourTimer = null, tourIdx = 0, tourStepping = false;
const tourHud = $("tour-hud");
function runStep(gen, i) {
  if (gen !== tourGen) return;                 // a cancelled/restarted tour
  if (i >= TOUR.length) { cancelTour(); return; }
  tourIdx = i; const s = TOUR[i];
  tourStepping = true; s.go(); tourStepping = false;   // the tour's own selects don't self-cancel
  state.camEaseRate = 1.6;   // slow, opening-style glide between steps instead of a snap
  if (tourHud) tourHud.textContent = s.label + "  ·  " + (i + 1) + " / " + TOUR.length + "  ·  tap or Esc to exit";
  tourTimer = setTimeout(() => runStep(gen, i + 1), s.dwell);
}
function startTour() {
  if (tourActive) return;
  tourActive = true; tourGen++; tourIdx = 0;
  if (!state.playing) { state.playing = true; syncPlayBtn(); }
  if (tourHud) tourHud.hidden = false;
  $("btn-tour").classList.add("on");
  runStep(tourGen, 0);
}
function cancelTour() {
  if (!tourActive) return;
  tourActive = false; tourGen++;
  clearTimeout(tourTimer); tourTimer = null;
  if (tourHud) tourHud.hidden = true;
  $("btn-tour").classList.remove("on");
}
$("btn-tour").addEventListener("click", () => { tourActive ? cancelTour() : startTour(); });
if (tourHud) tourHud.addEventListener("click", (e) => { e.stopPropagation(); clearTimeout(tourTimer); runStep(tourGen, tourIdx + 1); });

/* ---- keyboard shortcuts ---- */
function cycleSpeed(dir) {
  const chips = [...$("speed-chips").children];
  let idx = chips.findIndex((c) => c.classList.contains("on"));
  idx = clamp((idx < 0 ? 2 : idx) + dir, 0, chips.length - 1);
  const b = chips[idx];
  state.dps = +b.dataset.dps;
  if (!state.playing) { state.playing = true; syncPlayBtn(); }
  chips.forEach((x) => x.classList.toggle("on", x === b));
}
function hasSelection() {
  return state.selected || state.selPlanet != null || state.selSun || state.selMoon != null ||
         state.selSat != null || state.selCraft != null || currentRegion;
}
window.addEventListener("keydown", (ev) => {
  if (ev.ctrlKey || ev.metaKey || ev.altKey) return;
  const t = document.activeElement;
  if (t === searchEl || /^(INPUT|TEXTAREA)$/.test(t.tagName)) return;
  switch (ev.key) {
    case "Escape":
      if (tourActive) { cancelTour(); ev.preventDefault(); }
      else if (hasSelection()) { clearSelection(); ev.preventDefault(); }
      return;
    case " ": ev.preventDefault(); state.playing = !state.playing; syncPlayBtn(); return;
    case "ArrowRight": case ".": case "]": ev.preventDefault(); cycleSpeed(1); return;
    case "ArrowLeft": case ",": case "[": ev.preventDefault(); cycleSpeed(-1); return;
    case "n": case "N": state.simJD = jdNow(); state.needFullUpdate = true; buildPlanetOrbits(); return;
  }
});

/* panels: close buttons + mobile FAB toggles */
document.querySelectorAll(".panel-close").forEach((b) =>
  b.addEventListener("click", () => {
    const p = $(b.dataset.close);
    if (p.id === "panel-info") { clearSelection(); return; }
    p.classList.remove("open");
    p.classList.add("closed");
  })
);
const FAB_PANELS = [["panel-legend", "fab-legend"], ["panel-cad", "fab-cad"], ["panel-sentry", "fab-sentry"], ["panel-events", "fab-events"]];
function togglePanel(id, fab) {
  const p = $(id);
  const isOpen = p.classList.contains("open") || (!p.classList.contains("closed") && getComputedStyle(p).display !== "none");
  if (!isOpen) {
    // opening one closes the other FAB-owned panels so they never overlap
    for (const [pid, fid] of FAB_PANELS) {
      if (pid === id) continue;
      $(pid).classList.add("closed"); $(pid).classList.remove("open");
      $(fid).classList.remove("on");
    }
  }
  if (isOpen) { p.classList.remove("open"); p.classList.add("closed"); }
  else { p.classList.add("open"); p.classList.remove("closed"); }
  if (fab) fab.classList.toggle("on", !isOpen);
}
$("fab-legend").addEventListener("click", (ev) => togglePanel("panel-legend", ev.currentTarget));
// the legend ships open on desktop so the color code is always explained
if (window.matchMedia && window.matchMedia("(min-width: 861px)").matches) $("fab-legend").classList.add("on");
$("fab-cad").addEventListener("click", (ev) => togglePanel("panel-cad", ev.currentTarget));
$("fab-sentry").addEventListener("click", (ev) => togglePanel("panel-sentry", ev.currentTarget));
$("fab-events").addEventListener("click", (ev) => togglePanel("panel-events", ev.currentTarget));

function setTopDown(on, skipCamera) {
  state.topDown = on;
  $("fab-view").classList.toggle("on", on);
  if (skipCamera) return;
  if (on) {
    state.savedPitch = state.cam.tPitch;
    state.cam.tPitch = 1.45;
  } else {
    state.cam.tPitch = state.savedPitch;
  }
}
$("fab-sun").addEventListener("click", resetToSun);
$("fab-view").addEventListener("click", () => setTopDown(!state.topDown));
$("retry-btn").addEventListener("click", loadAsteroids);

/* keep the planet-orbit lines honest if the user time-travels far */
setInterval(() => {
  buildPlanetOrbits();
}, 30000);

/* ============================================================
   12. Lift-off
   ============================================================ */
window.addEventListener("error", (ev) => {
  // surface fatal init errors instead of an eternal spinner
  const loader = $("loader");
  if (loader && !loader.classList.contains("done") && !loadedOnce) {
    $("loader-status").textContent = "";
    $("loader-error-msg").textContent = "Something went wrong while starting up: " + (ev.message || "unknown error");
    $("loader-error").hidden = false;
  }
});
requestAnimationFrame((t) => { lastFrame = t; requestAnimationFrame(render); });
loadAsteroids();

// preload just the satellite COUNT (tiny) so the man-made figure is correct
// from the first frame, without the 2.5 MB catalogue
fetchJSON("data/sat-count.json")
  .then((d) => {
    if (!d || d.count == null) return;
    satCountKnown = d.payloads != null ? d.payloads : d.count;   // functional payloads
    debrisCountKnown = d.debris || 0;
    updateManMade();
  })
  .catch(() => {});

})();
