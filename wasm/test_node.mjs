/**
 * HPOP WASM smoke tests
 * Tests both simple JSON and VCM-schema JSON input paths
 */
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';

const __dirname = dirname(fileURLToPath(import.meta.url));

const HPOP_WASM = (await import(join(__dirname, 'node', 'hpop_wasm.js'))).default;
const Module = await HPOP_WASM();

console.log('✅ HPOP WASM loaded');
console.log('   Version:', Module.getVersion());
console.log('   Force models:', Module.getForceModels());

// Test 1: Simple JSON propagation
console.log('\n📡 Test 1: Simple JSON propagation...');
const t0 = performance.now();
const result1 = Module.propagate(JSON.stringify({
    epoch: "2024-01-01T00:00:00.000Z",
    x: 6778.0, y: 0, z: 0, vx: 0, vy: 7.669, vz: 0,
    duration_days: 0.1, step_seconds: 60.0,
    mass_kg: 260, drag_area: 22, cd: 2.2, srp_area: 22, cr: 1.3,
    f107: 150.0
}));
console.log(`   Done in ${(performance.now() - t0).toFixed(1)} ms, ${result1.length} bytes`);

// Test 2: VCM-schema JSON propagation
console.log('\n📡 Test 2: VCM-schema JSON propagation...');
const t1 = performance.now();
const result2 = Module.propagateVCM(JSON.stringify({
    OBJECT_NAME: "ISS (ZARYA)",
    OBJECT_ID: "1998-067A",
    epoch: "2024-01-01T00:00:00.000Z",
    X: -2517.232, Y: -663.786, Z: 6691.339,
    X_DOT: -1.768, Y_DOT: 7.315, Z_DOT: 1.439,
    MASS: 420000.0,
    DRAG_AREA: 1600.0,
    DRAG_COEFF: 2.2,
    SOLAR_RAD_AREA: 2500.0,
    SOLAR_RAD_COEFF: 1.3,
    NORAD_CAT_ID: 25544,
    duration_days: 0.1,
    step_seconds: 60.0,
    f107: 150.0
}));
console.log(`   Done in ${(performance.now() - t1).toFixed(1)} ms, ${result2.length} bytes`);

// Test 3: VCM with Keplerian elements
console.log('\n📡 Test 3: Keplerian elements...');
const t2 = performance.now();
const result3 = Module.propagateVCM(JSON.stringify({
    OBJECT_NAME: "STARLINK-1234",
    epoch: "2024-06-21T12:00:00.000Z",
    SEMI_MAJOR_AXIS: 6921.0,
    ECCENTRICITY: 0.0001,
    INCLINATION: 53.0,
    RA_OF_ASC_NODE: 0.0,
    ARG_OF_PERICENTER: 0.0,
    ANOMALY: 0.0,
    MASS: 260.0,
    DRAG_AREA: 22.0,
    DRAG_COEFF: 2.2,
    duration_days: 0.1,
    step_seconds: 60.0
}));
console.log(`   Done in ${(performance.now() - t2).toFixed(1)} ms, ${result3.length} bytes`);

console.log('\n✅ All smoke tests passed');
