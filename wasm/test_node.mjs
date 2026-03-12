/**
 * Quick smoke test: load HPOP WASM in Node.js and run a propagation
 */
import { createRequire } from 'module';
import { readFileSync } from 'fs';
import { fileURLToPath } from 'url';
import { dirname, join } from 'path';

const __dirname = dirname(fileURLToPath(import.meta.url));

// Load the WASM module
const HPOP_WASM = (await import(join(__dirname, 'node', 'hpop_wasm.js'))).default;
const Module = await HPOP_WASM();

console.log('✅ HPOP WASM loaded');
console.log('   Version:', Module.getVersion());
console.log('   Force models:', Module.getForceModels());

// Run a propagation
const request = JSON.stringify({
    epoch: "2024-01-01T00:00:00.000Z",
    x: 6778.0, y: 0, z: 0,
    vx: 0, vy: 7.669, vz: 0,
    duration_days: 0.1,
    step_seconds: 60.0,
    mass_kg: 260, drag_area: 22, cd: 2.2, srp_area: 22, cr: 1.3,
    forces: ["twobody", "j2", "drag", "srp", "sun", "moon"],
    f107: 150.0,
    model: "harris-priester"
});

console.log('\n📡 Propagating ISS-like orbit for 0.1 days...');
const t0 = performance.now();
const oemBinary = Module.propagate(request);
const dt = performance.now() - t0;

console.log(`   Done in ${dt.toFixed(1)} ms`);
console.log(`   OEM binary size: ${oemBinary.length} bytes`);
console.log(`   File identifier: ${oemBinary.substring(4, 8)}`);

console.log('\n✅ All smoke tests passed');
