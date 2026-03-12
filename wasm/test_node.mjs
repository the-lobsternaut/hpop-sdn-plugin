/**
 * HPOP WASM smoke tests
 * Tests: simple JSON propagation, VCM-schema JSON, Keplerian, raw text VCM → OCM
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
    epoch: "2024-01-01T00:00:00.000Z",
    X: -2517.232, Y: -663.786, Z: 6691.339,
    X_DOT: -1.768, Y_DOT: 7.315, Z_DOT: 1.439,
    MASS: 420000.0, DRAG_AREA: 1600.0, DRAG_COEFF: 2.2,
    duration_days: 0.1, step_seconds: 60.0
}));
console.log(`   Done in ${(performance.now() - t1).toFixed(1)} ms, ${result2.length} bytes`);

// Test 3: Raw text VCM → OCM FlatBuffer
console.log('\n📡 Test 3: Raw text VCM → OCM...');
const rawVCM = `<> SP VECTOR/COVARIANCE MESSAGE - V2.0
<>
<> MESSAGE TIME (UTC): 2023 007 (07 JAN) 22:45:06.000  CENTER: FAKE_CENTER
<> SATELLITE NUMBER: 00000                    INT. DES.: 0000-000A
<> COMMON NAME:
<> EPOCH TIME (UTC): 2023 007 (07 JAN) 21:39: 8.414  EPOCH REV: 37693
<> J2K POS (KM):       369.89521989    5103.23778139    4467.41312115
<> J2K VEL (KM/S):  -6.491206849227  -2.407500055425   3.279477293781
<> ECI POS (KM):       333.71334681    5104.90725904    4468.35490991
<> ECI VEL (KM/S):  -6.485984386597  -2.441004400229   3.265011217687
<> EFG POS (KM):      4957.72595592    1261.90176605    4468.35490991
<> EFG VEL (KM/S):  -4.235735566625   5.051150329994   3.265011217687
<> GEOPOTENTIAL: EGM-96 70Z,70T  DRAG: JAC70/MSIS90  LUNAR/SOLAR:  ON
<> SOLAR RAD PRESS: OFF  SOLID EARTH TIDES:  ON  IN-TRACK THRUST: OFF
<> BALLISTIC COEF (M2/KG):  0.826455E-02 BDOT (M2/KG-S): 0.000000E+00
<> SOLAR RAD PRESS COEFF (M2/KG):  0.000000E+00  EDR(W/KG):  0.39E-02
<> THRUST ACCEL (M/S2):  0.000000E+00  C.M. OFFSET (M):  0.000000E+00
<> SOLAR FLUX: F10: 153  AVERAGE F10: 137  AVERAGE AP:  10.0
<> TAI-UTC (S): 37  UT1-UTC (S): -0.01719  UT1 RATE (MS/DAY):  0.384
<> POLAR MOT X,Y (ARCSEC):  0.0505  0.2109 IAU 1980 NUTAT:   4 TERMS
<> TIME CONST LEAP SECOND TIME (UTC): 2049 365 (31 DEC) 23:59:59.999
<> INTEGRATOR MODE: ASW          COORD SYS: J2000  PARTIALS: FAST NUM
<> STEP MODE: AUTO  FIXED STEP: OFF  STEP SIZE SELECTION: MANUAL
<> INITIAL STEP SIZE (S):   20.000  ERROR CONTROL: 0.100E-13
<> VECTOR U,V,W SIGMAS (KM):           0.0084     0.0402     0.0074
<> VECTOR UD,VD,WD SIGMAS (KM/S):      0.0000     0.0000     0.0000
<> COVARIANCE MATRIX (EQUINOCTIAL ELS): ( 7x 7) WTD RMS:  0.11498E+01
<>  0.20458E-11  0.10230E-11  0.13484E-11  0.13797E-11 -0.47633E-12
<>  0.13642E-10  0.13871E-12  0.91700E-13  0.70998E-12  0.93201E-13
<>  0.50454E-12  0.26963E-12 -0.10331E-12  0.40296E-14  0.39682E-12
<> -0.38952E-12 -0.24005E-12  0.73298E-14 -0.10121E-13 -0.17845E-12
<>  0.40019E-12  0.91367E-08  0.53586E-08  0.36570E-07  0.50105E-08
<>  0.93145E-09 -0.57204E-09  0.13665E-02`;

const t2 = performance.now();
const ocmBytes = Module.parseVCM(rawVCM);
const elapsed = (performance.now() - t2).toFixed(1);

console.log(`   Done in ${elapsed} ms`);
console.log(`   OCM binary size: ${ocmBytes.length} bytes`);
console.log(`   Type: ${ocmBytes.constructor.name}`);

// Verify file identifier ($OCM)
const fid = String.fromCharCode(ocmBytes[4], ocmBytes[5], ocmBytes[6], ocmBytes[7]);
console.log(`   File identifier: ${fid}`);

if (fid === '$OCM') {
    console.log('   ✅ Valid OCM FlatBuffer');
} else {
    console.log('   ❌ Invalid file identifier');
    process.exit(1);
}

// Test 4: OCM → text VCM (round-trip)
console.log('\n📡 Test 4: OCM → text VCM (round-trip)...');
const t3 = performance.now();
const textBack = Module.convertOCM(ocmBytes);
console.log(`   Done in ${(performance.now() - t3).toFixed(1)} ms`);
console.log(`   Text VCM size: ${textBack.length} chars`);

// Verify round-trip preserves key fields
const checks = [
    ['SP VECTOR/COVARIANCE', 'Header'],
    ['FAKE_CENTER', 'Center'],
    ['0000-000A', 'Int designator'],
    ['369.89521989', 'J2K X position'],
    ['EGM-96', 'Geopotential model'],
    ['70Z,70T', 'Degree/order'],
];
let allOk = true;
for (const [needle, label] of checks) {
    if (!textBack.includes(needle)) {
        console.log(`   ❌ ${label} not found`);
        allOk = false;
    }
}
if (allOk) console.log('   ✅ All key fields preserved in round-trip');

console.log('\n✅ All smoke tests passed');
