# HPOP SDN Plugin — High Precision Orbit Propagator

## Overview

Numerical orbit propagator with pluggable force models, designed for SP-compatible orbit determination. Compiles to native and WASM. Integrates into the OD plugin as an alternative to SGP4.

## Architecture

```
hpop-sdn-plugin/
  src/cpp/
    include/hpop/
      propagator.h        # Cowell integrator (RKF-45, RKF-78, Dormand-Prince 853)
      gravity.h           # Spherical harmonics geopotential (EGM-96)
      atmosphere.h        # Atmosphere model interface
      jb2008.h            # JB2008 atmosphere (HASDM baseline)
      nrlmsise00.h        # NRLMSISE-00 atmosphere
      srp.h               # Solar radiation pressure
      thirdbody.h         # Sun/Moon third-body perturbations (DE430)
      frames.h            # TEME ↔ J2000 ↔ ECEF transforms
      state.h             # State vector types
      force_model.h       # Force model interface
    src/
      propagator.cpp      # RKF-45/78 + Dormand-Prince 853 integrators
      gravity.cpp         # EGM-96 spherical harmonics acceleration
      atmosphere.cpp      # Atmosphere model dispatcher
      jb2008.cpp          # JB2008 implementation
      nrlmsise00.cpp      # NRLMSISE-00 implementation
      srp.cpp             # Solar radiation pressure (dual-cone)
      thirdbody.cpp       # Sun/Moon ephemeris + perturbation
      frames.cpp          # Reference frame conversions
      wasm_api.cpp        # C API for WASM exports
    deps/
      vallado-astpert/    # Vallado's astPert (reference, AGPL-3.0)
      vallado-astrolib/   # Vallado's AstroLib
      vallado-mathtime/   # Vallado's MathTimeLib
      vallado-msis/       # NRLMSISE-00 from Vallado
      egm96/              # EGM-96 coefficients file
    tests/
      test_twobody.cpp         # Two-body Kepler orbit (analytical solution)
      test_j2.cpp              # J2-only propagation
      test_geopotential.cpp    # Full geopotential (vs Tudat reference)
      test_drag.cpp            # Atmospheric drag (vs OreKit reference)
      test_srp.cpp             # SRP (vs OreKit reference)
      test_thirdbody.cpp       # Sun/Moon perturbations
      test_full_perturbation.cpp  # All forces combined
      test_integrators.cpp     # RKF-45/78/DP853 accuracy & step control
      data/
        egm96_to70.txt         # EGM-96 coefficients (70×70)
        de430_excerpt.bin      # DE430 Sun/Moon ephemeris (test window)
        tudat_reference/       # Reference trajectories from Tudat
        orekit_reference/      # Reference trajectories from OreKit
  dist/
    hpop_wasm.js
    hpop_wasm.wasm
  plugin-manifest.json
  README.md
```

## Test Strategy

### Level 1: Unit Tests (analytical reference)
- **Two-body**: Propagate circular/elliptical orbits, compare to Kepler equation (< 1e-12 km)
- **J2-only**: Known J2 secular rates for RAAN and argument of perigee
- **Integrator accuracy**: RKF-45 vs RK4 on known ODE (Tudat's Matlab reference data)

### Level 2: Component Tests (cross-validated)
- **Geopotential**: EGM-96 36×36 acceleration at known points (vs Tudat SphericalHarmonicsGravityModel)
- **Drag**: JB2008/NRLMSISE-00 density at known altitudes/conditions
- **SRP**: Radiation pressure acceleration at known Sun angles
- **Third-body**: Sun/Moon acceleration at known epochs (vs DE430)

### Level 3: Integration Tests (full propagation)
- **Tudat reference orbit**: a=7500km, e=0.1, i=85.3°, with point mass + J2 + SRP + third-body
  - RK4 step=10s, propagate 1 day, compare Cartesian state (< 1e-6 km)
- **OreKit reference orbit**: LEO with HolmesFeatherstone 36×36 + DTM2000 + SRP + Sun/Moon
  - DormandPrince853, propagate 1 day, compare (< 1e-4 km)
- **Starlink-like orbit**: a=6921km (550km alt), i=53°, full perturbation suite
  - Compare to MEME ephemeris truth data (same test data as SGP4 plugin)

### Level 4: SP Compatibility
- Output VCM with correct model metadata
- Output OCM with state time-series
- Verify VCM fields match Space-Track format expectations

## Force Models — Implementation Order

1. **Two-body** (point mass gravity) — baseline, analytical validation
2. **Geopotential** (EGM-96, configurable degree/order up to 70×70)
3. **Atmospheric drag** (JB2008 first, then NRLMSISE-00)
4. **Solar radiation pressure** (dual-cone model)
5. **Third-body** (Sun + Moon from DE430 or analytical)
6. **Solid Earth tides** (optional, for SP completeness)
7. **General relativity** (optional, for SP completeness)

## Integrators

- **RKF-45** — Primary workhorse, adaptive step (Vallado has this)
- **RKF-78** — Higher order, better for long propagations (Vallado has this)
- **Dormand-Prince 853** — OreKit's default, 8th order (implement from Butcher tableau)

## SP Compatibility Requirements

Space Force SP uses HASDM with DCA (Dynamic Calibration of the Atmosphere):
- HASDM wraps JB2008 with real-time density corrections from ~75 calibration satellites
- HASDM itself is ITAR-controlled and not publicly available
- SP compatibility = same VCM format + same force model categories + comparable fidelity
- We implement JB2008 (the baseline HASDM wraps) and flag it in VCM metadata

VCM output must include:
- `atmosphericModel`: JB2008 or NRLMSISE_00
- `geopotentialModel`: EGM96
- `propagatorType`: COWELL
- Full 6×6 covariance matrix (21 lower-triangular elements)
- Equinoctial + Keplerian + Cartesian state

## Integration with OD Plugin

The OD plugin (`od-sdn-plugin`) needs a propagator interface:

```cpp
// In od-sdn-plugin
class Propagator {
public:
    virtual bool init(const OrbitalState& state) = 0;
    virtual bool propagate(double dt_sec, double r[3], double v[3]) = 0;
    virtual bool propagate_batch(const std::vector<double>& times,
                                 std::vector<PropState>& out) = 0;
};

class SGP4Propagator : public Propagator { ... };  // existing
class HPOPPropagator : public Propagator { ... };  // new, links hpop-sdn-plugin
```

The OD batch least squares / EKF code then works with either propagator.

## Dependencies

- Eigen3 (linear algebra)
- Vallado astPert/AstroLib/MathTimeLib (reference implementations, AGPL-3.0)
- EGM-96 coefficients (public domain, from NGA)
- DE430 Sun/Moon ephemeris (public domain, from JPL)
- Emscripten (WASM compilation)

## References

- Vallado — *Fundamentals of Astrodynamics and Applications* 5th ed.
- Montenbruck & Gill — *Satellite Orbits: Models, Methods, Applications*
- Tudat — `tudat-team/tudat` (BSD-3, test vectors)
- OreKit — `CS-SI/Orekit` (Apache-2.0, test vectors)
- EGM-96 — https://earth-info.nga.mil/GandG/wgs84/gravitymod/egm96/
- DE430 — https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/
- JB2008 — Bowman et al. 2008, https://sol.spacenvironment.net/jb2008/
- NRLMSISE-00 — Picone et al. 2002
