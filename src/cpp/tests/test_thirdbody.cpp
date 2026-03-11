/**
 * HPOP Third-Body Perturbation Tests
 *
 * Test strategy:
 *   1. Sun position: validate against known almanac values at multiple epochs
 *   2. Moon position: validate distance and direction at known epochs
 *   3. Sun perturbation magnitude and direction
 *   4. Moon perturbation magnitude (should be ~2× Sun for LEO)
 *   5. Combined Sun+Moon effect on GEO orbit
 *   6. Lunar resonance: GEO inclination growth
 *   7. Battin vs standard formulation equivalence
 *
 * Reference values:
 *   - Sun perturbation on LEO: ~5e-7 m/s² (Montenbruck & Gill)
 *   - Moon perturbation on LEO: ~1e-6 m/s² (about 2× Sun)
 *   - GEO inclination drift: ~0.75-0.95°/year from Sun+Moon
 */

#include "hpop/propagator.h"
#include "hpop/thirdbody.h"
#include <cstdio>
#include <cmath>
#include <memory>

using namespace hpop;

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", msg); tests_failed++; } \
    else { printf("  PASS: %s\n", msg); tests_passed++; } \
} while(0)

#define CHECK_TOL(val, expected, tol, msg) do { \
    double _v = (val), _e = (expected), _t = (tol); \
    if (std::abs(_v - _e) > _t) { \
        printf("  FAIL: %s (got %.6e, expected %.6e, diff %.3e)\n", msg, _v, _e, std::abs(_v-_e)); \
        tests_failed++; \
    } else { \
        printf("  PASS: %s (diff %.3e < %.3e)\n", msg, std::abs(_v-_e), _t); \
        tests_passed++; \
    } \
} while(0)

// ── Test 1: Sun Position ──

void test_sun_position() {
    printf("\n=== Test 1: Sun Position (ECI) ===\n");

    double r[3];

    // J2000 epoch (Jan 1, 2000 12:00 TT)
    sun_position_eci(2451545.0, r);
    double dist = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double ra = std::atan2(r[1], r[0]) * RAD2DEG;
    if (ra < 0) ra += 360;
    double dec = std::asin(r[2] / dist) * RAD2DEG;

    printf("  J2000: dist=%.0f km (%.4f AU), RA=%.1f°, Dec=%.1f°\n",
           dist, dist/AU_KM, ra, dec);

    // J2000 Sun: RA≈281°, Dec≈-23°, dist≈0.983 AU
    CHECK_TOL(dist / AU_KM, 0.983, 0.005, "Sun distance ≈ 0.983 AU at J2000");
    CHECK_TOL(ra, 281.3, 1.0, "Sun RA ≈ 281° at J2000");
    CHECK_TOL(dec, -23.0, 1.0, "Sun Dec ≈ -23° at J2000");

    // Summer solstice 2000 (Jun 21): Sun RA≈90°, Dec≈+23.4°
    sun_position_eci(2451716.5, r);
    dist = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    ra = std::atan2(r[1], r[0]) * RAD2DEG;
    if (ra < 0) ra += 360;
    dec = std::asin(r[2] / dist) * RAD2DEG;
    printf("  Jun 21: dist=%.4f AU, RA=%.1f°, Dec=%.1f°\n", dist/AU_KM, ra, dec);
    CHECK_TOL(dec, 23.4, 1.0, "Sun Dec ≈ +23.4° at summer solstice");
}

// ── Test 2: Moon Position ──

void test_moon_position() {
    printf("\n=== Test 2: Moon Position (ECI) ===\n");

    double r[3];

    // J2000 epoch
    moon_position_eci(2451545.0, r);
    double dist = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double ra = std::atan2(r[1], r[0]) * RAD2DEG;
    if (ra < 0) ra += 360;
    double dec = std::asin(r[2] / dist) * RAD2DEG;

    printf("  J2000: dist=%.0f km, RA=%.1f°, Dec=%.1f°\n", dist, ra, dec);

    // Moon distance should be 356,500-406,700 km (perigee to apogee)
    CHECK(dist > 350000 && dist < 410000, "Moon distance in valid range (350-410 Mm)");

    // Check at different epoch (1 month later)
    moon_position_eci(2451545.0 + 29.53, r);
    double dist2 = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    printf("  +1 synodic month: dist=%.0f km\n", dist2);
    CHECK(dist2 > 350000 && dist2 < 410000, "Moon distance valid after 1 month");

    // Verify Moon moves (~13°/day in ecliptic longitude)
    moon_position_eci(2451545.0, r);
    double r1[3];
    moon_position_eci(2451546.0, r1);
    double cos_angle = (r[0]*r1[0] + r[1]*r1[1] + r[2]*r1[2])
                     / (std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
                       * std::sqrt(r1[0]*r1[0]+r1[1]*r1[1]+r1[2]*r1[2]));
    double daily_motion = std::acos(std::max(-1.0, std::min(1.0, cos_angle))) * RAD2DEG;
    printf("  Daily motion: %.1f°/day (expected ~13°)\n", daily_motion);
    CHECK_TOL(daily_motion, 13.0, 3.0, "Moon daily motion ≈ 13°/day");
}

// ── Test 3: Sun Perturbation Magnitude ──

void test_sun_perturbation_mag() {
    printf("\n=== Test 3: Sun Third-Body Perturbation ===\n");

    ThirdBodyForce sun_force(ThirdBody::SUN);

    // LEO satellite at 400 km
    double r[3] = {RE_KM + 400.0, 0, 0};
    double v[3] = {0, std::sqrt(MU_EARTH / (RE_KM + 400.0)), 0};

    double a[3];
    sun_force.acceleration(r, v, 2451545.0, nullptr, a);
    double a_mag = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

    printf("  |a_sun| = %.3e km/s² = %.3e m/s²\n", a_mag, a_mag * 1000);

    // Expected: ~5e-7 m/s² = ~5e-10 km/s² for LEO
    CHECK(a_mag > 1e-11 && a_mag < 1e-8, "Sun perturbation in expected range (1e-11 to 1e-8 km/s²)");

    // Compare to two-body gravity
    double g = MU_EARTH / ((RE_KM + 400.0) * (RE_KM + 400.0));
    printf("  a_sun/g = %.3e\n", a_mag / g);
    CHECK(a_mag / g < 1e-6, "Sun perturbation << gravity for LEO");
}

// ── Test 4: Moon Perturbation Magnitude ──

void test_moon_perturbation_mag() {
    printf("\n=== Test 4: Moon Third-Body Perturbation ===\n");

    ThirdBodyForce moon_force(ThirdBody::MOON);

    double r[3] = {RE_KM + 400.0, 0, 0};
    double v[3] = {0, std::sqrt(MU_EARTH / (RE_KM + 400.0)), 0};

    double a[3];
    moon_force.acceleration(r, v, 2451545.0, nullptr, a);
    double a_mag_moon = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);

    printf("  |a_moon| = %.3e km/s² = %.3e m/s²\n", a_mag_moon, a_mag_moon * 1000);

    // Moon perturbation should be roughly 2× Sun for LEO
    ThirdBodyForce sun_force(ThirdBody::SUN);
    double a_sun[3];
    sun_force.acceleration(r, v, 2451545.0, nullptr, a_sun);
    double a_mag_sun = std::sqrt(a_sun[0]*a_sun[0] + a_sun[1]*a_sun[1] + a_sun[2]*a_sun[2]);

    double ratio = a_mag_moon / a_mag_sun;
    printf("  Moon/Sun ratio: %.2f (expected ~2.2)\n", ratio);
    CHECK(ratio > 1.0 && ratio < 5.0, "Moon perturbation 1-5× Sun for LEO");
}

// ── Test 5: GEO Orbit with Sun+Moon ──

void test_geo_perturbation() {
    printf("\n=== Test 5: GEO Orbit — Sun+Moon Perturbations ===\n");

    // GEO orbit: a ≈ 42164 km, e ≈ 0, i ≈ 0
    KeplerianElements kep = {42164.0, 0.0001, 0.1 * DEG2RAD, 0.0, 0.0, 0.0};
    State s0 = keplerian_to_cartesian(kep);

    double epoch_jd0 = 2451545.0;
    double epoch_jdf = epoch_jd0 + 30.0;  // 30 days

    // Two-body only
    Propagator prop_2b;
    prop_2b.add_force(std::make_unique<TwoBodyForce>());
    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF78;
    cfg.step_size = 60.0;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-12;
    prop_2b.set_config(cfg);
    State sf_2b = prop_2b.propagate(s0, epoch_jd0, epoch_jdf);

    // Two-body + Sun + Moon
    Propagator prop_3b;
    prop_3b.add_force(std::make_unique<TwoBodyForce>());
    prop_3b.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::SUN));
    prop_3b.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::MOON));
    prop_3b.set_config(cfg);
    State sf_3b = prop_3b.propagate(s0, epoch_jd0, epoch_jdf);

    auto kep_2b = cartesian_to_keplerian(sf_2b);
    auto kep_3b = cartesian_to_keplerian(sf_3b);

    printf("  2-body: a=%.3f km, e=%.6f, i=%.4f°\n",
           kep_2b.a, kep_2b.e, kep_2b.i * RAD2DEG);
    printf("  3-body: a=%.3f km, e=%.6f, i=%.4f°\n",
           kep_3b.a, kep_3b.e, kep_3b.i * RAD2DEG);

    double di = (kep_3b.i - kep_2b.i) * RAD2DEG;
    double de = kep_3b.e - kep_2b.e;
    printf("  Δi = %.4f° (30 days)\n", di);
    printf("  Δe = %.6f (30 days)\n", de);

    // GEO inclination grows ~0.75-0.95°/year from Sun+Moon
    // In 30 days: ~0.06-0.08°
    CHECK(std::abs(di) > 0.01, "Measurable inclination change in 30 days");
    CHECK(std::abs(di) < 0.5, "Inclination change < 0.5° in 30 days (sanity)");

    // Position difference between 2-body and 3-body
    double dr[3] = {sf_3b.r[0]-sf_2b.r[0], sf_3b.r[1]-sf_2b.r[1], sf_3b.r[2]-sf_2b.r[2]};
    double dr_mag = std::sqrt(dr[0]*dr[0] + dr[1]*dr[1] + dr[2]*dr[2]);
    printf("  |Δr| = %.1f km after 30 days\n", dr_mag);
    CHECK(dr_mag > 1.0, "Third-body causes > 1 km position difference in 30 days");
    CHECK(dr_mag < 10000.0, "Position difference < 10000 km (sanity)");
}

// ── Test 6: LEO Propagation with Full Force Model ──

void test_leo_full_force() {
    printf("\n=== Test 6: LEO Full Force Model (2-body + J2 + Sun + Moon) ===\n");

    // ISS-like orbit
    KeplerianElements kep = {6778.0, 0.001, 51.6 * DEG2RAD, 0.0, 0.0, 0.0};
    State s0 = keplerian_to_cartesian(kep);

    double epoch_jd0 = 2451545.0;
    double epoch_jdf = epoch_jd0 + 1.0;  // 1 day

    Propagator prop;
    prop.add_force(std::make_unique<TwoBodyForce>());
    prop.add_force(std::make_unique<J2Force>());
    prop.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::SUN));
    prop.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::MOON));

    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF78;
    cfg.step_size = 30.0;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-12;
    prop.set_config(cfg);

    State sf = prop.propagate(s0, epoch_jd0, epoch_jdf);
    auto kep_f = cartesian_to_keplerian(sf);

    printf("  Initial: a=%.3f km, e=%.6f, i=%.4f°\n", kep.a, kep.e, kep.i * RAD2DEG);
    printf("  Final:   a=%.3f km, e=%.6f, i=%.4f°\n", kep_f.a, kep_f.e, kep_f.i * RAD2DEG);

    // J2 causes RAAN precession (~7°/day for ISS) and argp rotation
    // Osculating SMA oscillates due to J2 short-period variations (~few km)
    // but there is no secular SMA change from conservative forces
    CHECK(std::abs(kep_f.a - kep.a) < 10.0, "Osculating SMA change < 10 km (J2 short-period oscillation)");

    // RAAN should precess (J2 dominant effect)
    double draan = (kep_f.raan - kep.raan) * RAD2DEG;
    printf("  ΔRAAN = %.2f°/day\n", draan);
    // ISS RAAN precession: ~-7°/day
    CHECK(std::abs(draan) > 1.0, "RAAN precession > 1°/day from J2");
}

// ── Test 7: Battin vs Standard Formulation ──

void test_battin_equivalence() {
    printf("\n=== Test 7: Battin vs Standard Formulation ===\n");

    // Standard third-body: a = μ₃ · (r₃sat/|r₃sat|³ - r₃/|r₃|³)
    double r_sat[3] = {RE_KM + 400.0, 0, 0};
    double r_sun[3];
    sun_position_eci(2451545.0, r_sun);

    // Standard formulation
    double r3sat[3] = {r_sun[0]-r_sat[0], r_sun[1]-r_sat[1], r_sun[2]-r_sat[2]};
    double r3sat_mag = std::sqrt(r3sat[0]*r3sat[0] + r3sat[1]*r3sat[1] + r3sat[2]*r3sat[2]);
    double rsun_mag = std::sqrt(r_sun[0]*r_sun[0] + r_sun[1]*r_sun[1] + r_sun[2]*r_sun[2]);

    double a_std[3];
    for (int i = 0; i < 3; i++) {
        a_std[i] = MU_SUN * (r3sat[i]/(r3sat_mag*r3sat_mag*r3sat_mag) - r_sun[i]/(rsun_mag*rsun_mag*rsun_mag));
    }

    // Battin formulation (via ThirdBodyForce)
    ThirdBodyForce sun_force(ThirdBody::SUN);
    double v[3] = {0, 0, 0};  // velocity doesn't matter for third-body
    double a_battin[3];
    sun_force.acceleration(r_sat, v, 2451545.0, nullptr, a_battin);

    // They should agree to high precision
    double mag_std = std::sqrt(a_std[0]*a_std[0] + a_std[1]*a_std[1] + a_std[2]*a_std[2]);
    double mag_bat = std::sqrt(a_battin[0]*a_battin[0] + a_battin[1]*a_battin[1] + a_battin[2]*a_battin[2]);

    printf("  Standard:  |a| = %.10e km/s²\n", mag_std);
    printf("  Battin:    |a| = %.10e km/s²\n", mag_bat);
    printf("  Rel diff: %.3e\n", std::abs(mag_bat - mag_std) / mag_std);

    // Component-wise comparison
    for (int i = 0; i < 3; i++) {
        double rel_diff = std::abs(a_battin[i] - a_std[i]) / (std::abs(a_std[i]) + 1e-30);
        printf("  Component %d: std=%.10e, bat=%.10e, rel_diff=%.3e\n",
               i, a_std[i], a_battin[i], rel_diff);
    }

    CHECK_TOL(mag_bat, mag_std, mag_std * 1e-6, "Battin matches standard to 1 ppm");
}

int main() {
    printf("============================================================\n");
    printf("HPOP Third-Body Perturbation Tests\n");
    printf("============================================================\n");

    test_sun_position();
    test_moon_position();
    test_sun_perturbation_mag();
    test_moon_perturbation_mag();
    test_geo_perturbation();
    test_leo_full_force();
    test_battin_equivalence();

    printf("\n============================================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("============================================================\n");

    return tests_failed > 0 ? 1 : 0;
}
