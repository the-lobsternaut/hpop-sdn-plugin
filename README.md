# HPOP SDN Plugin

High-Precision Orbit Propagator for the [Space Data Network](https://github.com/the-lobsternaut/space-data-network). Parses 18 SPCS text VCMs, propagates SP orbits with a full force model stack, outputs OEM/OCM FlatBuffers.

**270 KB WebAssembly** — runs in Node.js and browsers.

## Force Models

| Model | Implementation | Reference |
|---|---|---|
| Two-body gravity | μ = 398600.4418 km³/s² | — |
| Geopotential | EGM-96, configurable degree/order | Montenbruck & Gill Ch. 3.2 |
| Atmospheric drag | Harris-Priester + DCA hook | M&G Table 3.3 |
| Solar radiation pressure | Cannonball, conical shadow | M&G Section 3.4.2 |
| Third-body Sun | Low-precision ephemeris, Battin formulation | Vallado, AstroLib |
| Third-body Moon | Low-precision ephemeris, Battin formulation | Vallado, AstroLib |
| Integrator | RKF-78 adaptive step | Fehlberg 1968 |

## Data Flow

```
┌─────────────┐     parse()      ┌──────────┐     propagate     ┌──────────┐
│ Text VCM    │ ──────────────→  │ OCM $OCM │                   │          │
│ (18 SPCS)   │                  │ FlatBuf  │                   │          │
└─────────────┘                  └──────────┘                   │          │
                                                                │ OEM $OEM │
┌─────────────┐     parse()      ┌──────────┐   propagate()     │ FlatBuf  │
│ VCM $VCM    │ ──────────────→  │ State +  │ ──────────────→  │          │
│ FlatBuffer  │                  │ Config   │                   │          │
└─────────────┘                  └──────────┘                   └──────────┘

┌─────────────┐   propagate()    ┌──────────┐
│ JSON        │ ──────────────→  │ OEM $OEM │
│ (simple)    │                  │ FlatBuf  │
└─────────────┘                  └──────────┘

┌──────────┐     convert()       ┌─────────────┐
│ OCM $OCM │ ──────────────→    │ Text VCM    │
│ FlatBuf  │                    │ (18 SPCS)   │
└──────────┘                    └─────────────┘
```

## Usage

### Node.js

```js
import HPOP from 'hpop-sdn-plugin';

const Module = await HPOP();

// 1. Parse raw text VCM → OCM FlatBuffer
const ocmBinary = Module.parseVCM(vcmText);  // Uint8Array

// 2. Convert OCM back to text VCM
const textVcm = Module.convertOCM(ocmBinary);

// 3. Propagate from JSON → OEM FlatBuffer
const oem = Module.propagate(JSON.stringify({
    epoch: "2024-01-01T00:00:00Z",
    x: 6778, y: 0, z: 0, vx: 0, vy: 7.669, vz: 0,
    duration_days: 1.0, step_seconds: 60,
    mass_kg: 260, drag_area: 22, cd: 2.2
}));

// 4. Propagate from VCM-schema JSON
const oem2 = Module.propagateVCM(JSON.stringify({
    OBJECT_NAME: "ISS",
    epoch: "2024-01-01T00:00:00Z",
    X: -2517, Y: -663, Z: 6691,
    X_DOT: -1.768, Y_DOT: 7.315, Z_DOT: 1.439,
    MASS: 420000, DRAG_AREA: 1600, DRAG_COEFF: 2.2,
    duration_days: 1.0
}));
```

### SDN Plugin C ABI

```c
// parse(): auto-detects input format
//   text VCM (<> prefix) → $OCM FlatBuffer
//   $VCM FlatBuffer → propagate → $OEM FlatBuffer
//   JSON → propagate → $OEM FlatBuffer
int32_t parse(const char* input, size_t input_len,
              uint8_t* output, size_t output_len);

// convert(): auto-detects FlatBuffer type
//   $OCM → text VCM
//   $OEM → CCSDS OEM text
int32_t convert(const uint8_t* input, size_t input_len,
                uint8_t* output, size_t output_len);
```

## Building

```bash
# Build WASM (auto-installs emsdk)
./build.sh

# Run native tests
cd src/cpp/build && cmake .. && make -j4
for t in test_*; do ./$t; done

# Run WASM smoke tests
node wasm/test_node.mjs
```

## VCM Format Support

Parses the full 18 SPCS V2.0 fixed-format VCM:

- State vectors: J2K, ECI, EFG coordinate systems
- Geopotential: model name, zonal degree, tesseral order
- Drag: JACCHIA 70, JB2008, NRLMSISE-00, HASDM, Harris-Priester
- Perturbations: Lunar/Solar, SRP, Solid Earth Tides, In-track Thrust
- Physical: Ballistic coeff (Cd·A/m), SRP coeff, EDR, thrust accel
- Space weather: F10.7, average F10.7, average Ap
- Time corrections: TAI-UTC, UT1-UTC, UT1 rate, polar motion
- Integrator: mode, coordinate system, partials type, step config
- Covariance: equinoctial elements, dynamic matrix size (6×6 to 12×12+)

## Tests

| Suite | Tests | Coverage |
|---|---|---|
| test_twobody | 13 | Kepler orbits, periods, circular/elliptical/hyperbolic |
| test_geopotential | 13 | J2 secular drift, EGM-96 validation |
| test_drag | 23 | Harris-Priester, ISS decay rate, altitude-dependent |
| test_srp | 17 | Cannonball SRP, conical shadow, eclipse transitions |
| test_thirdbody | 17 | Sun/Moon positions, Battin formulation, perturbation magnitude |
| test_wasm_api | 20 | Time conversion, propagation, OEM FlatBuffer generation |
| test_vcm_input | 38 | VCM FlatBuffer/JSON parsing, force model config |
| test_vcm_parser | 79 | Text VCM parsing, round-trip, covariance, edge cases |
| **Total** | **220** | |

## License

Apache-2.0
