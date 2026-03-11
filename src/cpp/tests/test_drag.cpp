/**
 * HPOP Atmospheric Drag Tests
 *
 * Test strategy:
 *   1. Exponential atmosphere density sanity check
 *   2. Harris-Priester density values match M&G table
 *   3. Harris-Priester diurnal variation (day > night)
 *   4. Drag acceleration direction (opposing velocity)
 *   5. Drag-induced orbit decay (circular LEO should lose altitude)
 *   6. Drag magnitude scaling with Cd·A/m (ballistic coefficient)
 *   7. Sun position accuracy vs known almanac values
 *
 * Reference test values from:
 *   - Montenbruck & Gill (2000), Table 3.3
 *   - Tudat exponential atmosphere test
 *   - Known orbital decay rates for ISS-like objects
 */

#include "hpop/propagator.h"
#include "hpop/atmosphere.h"
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

// ── Test 1: Exponential Atmosphere ──

void test_exponential() {
    printf("\n=== Test 1: Exponential Atmosphere ===\n");

    // ISS altitude range: ~400 km
    ExponentialAtmosphere atm(3.725e-12, 400.0, 58.515);

    double rho_400 = atm.density(400.0, 0, 0, 2451545.0);
    double rho_500 = atm.density(500.0, 0, 0, 2451545.0);
    double rho_300 = atm.density(300.0, 0, 0, 2451545.0);

    printf("  rho(400 km) = %.3e kg/m³\n", rho_400);
    printf("  rho(500 km) = %.3e kg/m³\n", rho_500);
    printf("  rho(300 km) = %.3e kg/m³\n", rho_300);

    CHECK_TOL(rho_400, 3.725e-12, 1e-20, "ρ(400 km) = reference density");
    CHECK(rho_500 < rho_400, "ρ decreases with altitude");
    CHECK(rho_300 > rho_400, "ρ increases below reference");

    // Check exponential scale
    double ratio = rho_400 / rho_500;
    double expected_ratio = std::exp(100.0 / 58.515);
    CHECK_TOL(ratio, expected_ratio, 1e-6, "Exponential ratio correct");
}

// ── Test 2: Harris-Priester Density Table ──

void test_harris_priester_table() {
    printf("\n=== Test 2: Harris-Priester Density Table ===\n");

    // F10.7 = 150 (moderate solar activity)
    HarrisPriester hp(150.0, 2);

    // At 300 km, nighttime ρ_min ≈ 1.708e-11, daytime ρ_max ≈ 3.526e-11
    // With F10.7 scaling: multiply by (1 + 0.5*((150-70)/(250-70))) ≈ 1.222
    double rho_300 = hp.density(300.0, 0, 0, 2451545.0);
    printf("  rho(300 km) = %.3e kg/m³\n", rho_300);

    // Should be between scaled min and max
    double min_scaled = 1.708e-11 * 1.222;
    double max_scaled = 3.526e-11 * 1.222;
    CHECK(rho_300 >= min_scaled * 0.9 && rho_300 <= max_scaled * 1.1,
          "ρ(300 km) within expected range");

    // At 500 km
    double rho_500 = hp.density(500.0, 0, 0, 2451545.0);
    printf("  rho(500 km) = %.3e kg/m³\n", rho_500);
    CHECK(rho_500 < rho_300, "ρ(500 km) < ρ(300 km)");

    // At 200 km
    double rho_200 = hp.density(200.0, 0, 0, 2451545.0);
    printf("  rho(200 km) = %.3e kg/m³\n", rho_200);
    CHECK(rho_200 > rho_300, "ρ(200 km) > ρ(300 km)");

    // Order of magnitude checks
    CHECK(rho_200 > 1e-11 && rho_200 < 1e-9, "ρ(200 km) ~ 1e-10 kg/m³");
    CHECK(rho_500 > 1e-14 && rho_500 < 1e-11, "ρ(500 km) ~ 1e-12 kg/m³");
}

// ── Test 3: Diurnal Variation ──

void test_diurnal_variation() {
    printf("\n=== Test 3: Harris-Priester Diurnal Variation ===\n");

    HarrisPriester hp(150.0, 2);

    // Get sun position at J2000
    double sun_ra, sun_dec;
    HarrisPriester::sun_position(2451545.0, sun_ra, sun_dec);
    printf("  Sun at J2000: RA=%.1f°, Dec=%.1f°\n", sun_ra * RAD2DEG, sun_dec * RAD2DEG);

    // Subsolar point (where bulge should be strongest)
    double rho_day = hp.density(400.0, sun_dec, sun_ra + 30.0 * DEG2RAD, 2451545.0);

    // Anti-solar point (nighttime minimum)
    double rho_night = hp.density(400.0, -sun_dec, sun_ra + M_PI + 30.0 * DEG2RAD, 2451545.0);

    printf("  rho_day(400 km) = %.3e kg/m³\n", rho_day);
    printf("  rho_night(400 km) = %.3e kg/m³\n", rho_night);
    printf("  Day/night ratio: %.2f\n", rho_day / rho_night);

    CHECK(rho_day > rho_night, "Daytime density > nighttime density");
    CHECK(rho_day / rho_night > 1.5, "Day/night ratio > 1.5 at 400 km");
    CHECK(rho_day / rho_night < 10.0, "Day/night ratio < 10 (physical limit)");
}

// ── Test 4: Drag Acceleration Direction ──

void test_drag_direction() {
    printf("\n=== Test 4: Drag Acceleration Direction ===\n");

    auto atm = std::make_unique<ExponentialAtmosphere>(3.725e-12, 400.0, 58.515);
    DragForce drag(std::move(atm), 2.2);

    // ISS-like orbit at 400 km, circular
    double r[3] = {6778.0, 0.0, 0.0};  // 400 km altitude
    double v_circ = std::sqrt(MU_EARTH / 6778.0);
    double v[3] = {0.0, v_circ, 0.0};

    SpacecraftProperties sc;
    sc.mass_kg = 420000.0;  // ISS mass
    sc.drag_area = 1600.0;    // ISS cross-section area

    double a_drag[3];
    drag.acceleration(r, v, 2451545.0, &sc, a_drag);

    printf("  v = [%.3f, %.3f, %.3f] km/s\n", v[0], v[1], v[2]);
    printf("  a_drag = [%.3e, %.3e, %.3e] km/s²\n", a_drag[0], a_drag[1], a_drag[2]);

    // Drag should oppose velocity
    double dot = a_drag[0]*v[0] + a_drag[1]*v[1] + a_drag[2]*v[2];
    printf("  a_drag · v = %.3e (should be negative)\n", dot);
    CHECK(dot < 0, "Drag opposes velocity (a·v < 0)");

    // Drag magnitude should be tiny compared to gravity
    double a_mag = std::sqrt(a_drag[0]*a_drag[0] + a_drag[1]*a_drag[1] + a_drag[2]*a_drag[2]);
    double g = MU_EARTH / (6778.0 * 6778.0);
    printf("  |a_drag| = %.3e km/s², g = %.3e km/s²\n", a_mag, g);
    printf("  a_drag/g = %.3e\n", a_mag / g);
    CHECK(a_mag / g < 1e-6, "Drag/gravity ratio < 1e-6 for ISS");
    CHECK(a_mag > 0, "Drag magnitude is positive");
}

// ── Test 5: Drag-Induced Orbit Decay ──

void test_orbit_decay() {
    printf("\n=== Test 5: Drag-Induced Orbit Decay ===\n");

    // ISS-like orbit: 400 km circular
    KeplerianElements kep = {6778.0, 0.001, 51.6 * DEG2RAD, 0.0, 0.0, 0.0};
    State s0 = keplerian_to_cartesian(kep);

    double epoch_jd0 = 2451545.0;
    double epoch_jdf = epoch_jd0 + 1.0;  // 1 day

    // Propagator with two-body + drag
    Propagator prop;
    prop.add_force(std::make_unique<TwoBodyForce>());

    auto atm = std::make_unique<HarrisPriester>(150.0, 2);
    SpacecraftProperties sc;
    sc.mass_kg = 420000.0;
    sc.drag_area = 1600.0;
    prop.add_force(std::make_unique<DragForce>(std::move(atm), 2.2));
    prop.set_spacecraft(sc);

    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF78;
    cfg.step_size = 60.0;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-12;
    prop.set_config(cfg);

    State sf = prop.propagate(s0, epoch_jd0, epoch_jdf);

    auto kep_f = cartesian_to_keplerian(sf);
    double da = kep_f.a - kep.a;  // Change in semi-major axis

    printf("  Initial a = %.3f km, Final a = %.3f km\n", kep.a, kep_f.a);
    printf("  Δa = %.4f km/day\n", da);

    // ISS loses ~50-200 m/day in SMA depending on solar activity
    // With our moderate F10.7=150, expect ~50-500 m/day
    CHECK(da < 0, "Semi-major axis decreases (orbit decays)");
    CHECK(std::abs(da) > 0.001, "Decay rate > 1 m/day");
    CHECK(std::abs(da) < 10.0, "Decay rate < 10 km/day (sanity)");

    printf("  Decay rate: %.1f m/day\n", std::abs(da) * 1000.0);
}

// ── Test 6: Ballistic Coefficient Scaling ──

void test_bc_scaling() {
    printf("\n=== Test 6: Ballistic Coefficient Scaling ===\n");

    auto atm1 = std::make_unique<ExponentialAtmosphere>(3.725e-12, 400.0, 58.515);
    auto atm2 = std::make_unique<ExponentialAtmosphere>(3.725e-12, 400.0, 58.515);
    DragForce drag1(std::move(atm1), 2.2);
    DragForce drag2(std::move(atm2), 2.2);

    double r[3] = {6778.0, 0.0, 0.0};
    double v[3] = {0.0, std::sqrt(MU_EARTH / 6778.0), 0.0};

    SpacecraftProperties sc1, sc2;
    sc1.mass_kg = 1000.0; sc1.drag_area = 10.0;  // Cd·A/m = 2.2*10/1000 = 0.022
    sc2.mass_kg = 1000.0; sc2.drag_area = 20.0;  // Cd·A/m = 2.2*20/1000 = 0.044

    double a1[3], a2[3];
    drag1.acceleration(r, v, 2451545.0, &sc1, a1);
    drag2.acceleration(r, v, 2451545.0, &sc2, a2);

    double mag1 = std::sqrt(a1[0]*a1[0] + a1[1]*a1[1] + a1[2]*a1[2]);
    double mag2 = std::sqrt(a2[0]*a2[0] + a2[1]*a2[1] + a2[2]*a2[2]);

    printf("  A=10 m²: |a_drag| = %.3e km/s²\n", mag1);
    printf("  A=20 m²: |a_drag| = %.3e km/s²\n", mag2);
    printf("  Ratio: %.4f (expected 2.0)\n", mag2 / mag1);

    CHECK_TOL(mag2 / mag1, 2.0, 0.001, "Doubling area doubles drag");
}

// ── Test 7: Sun Position ──

void test_sun_position() {
    printf("\n=== Test 7: Sun Position ===\n");

    double ra, dec;

    // J2000 epoch: Jan 1, 2000 12:00 TT
    // Sun RA ≈ 281.3° (18h 45m), Dec ≈ -23.0°
    HarrisPriester::sun_position(2451545.0, ra, dec);
    printf("  J2000: Sun RA=%.1f°, Dec=%.1f°\n", ra * RAD2DEG, dec * RAD2DEG);
    CHECK_TOL(ra * RAD2DEG, 281.3, 2.0, "J2000 Sun RA within 2° of expected");
    CHECK_TOL(dec * RAD2DEG, -23.0, 2.0, "J2000 Sun Dec within 2° of expected");

    // Vernal equinox: ~March 20, Sun RA ≈ 0°, Dec ≈ 0°
    // March 20, 2000 ≈ JD 2451623.5
    HarrisPriester::sun_position(2451623.5, ra, dec);
    printf("  Mar equinox: Sun RA=%.1f°, Dec=%.1f°\n", ra * RAD2DEG, dec * RAD2DEG);
    CHECK(std::abs(dec * RAD2DEG) < 3.0, "Equinox: Sun Dec near 0°");

    // Summer solstice: ~June 21, Sun Dec ≈ +23.4°
    // June 21, 2000 ≈ JD 2451716.5
    HarrisPriester::sun_position(2451716.5, ra, dec);
    printf("  Jun solstice: Sun RA=%.1f°, Dec=%.1f°\n", ra * RAD2DEG, dec * RAD2DEG);
    CHECK_TOL(dec * RAD2DEG, 23.4, 2.0, "Solstice: Sun Dec near +23.4°");
}

int main() {
    printf("============================================================\n");
    printf("HPOP Atmospheric Drag Tests\n");
    printf("============================================================\n");

    test_exponential();
    test_harris_priester_table();
    test_diurnal_variation();
    test_drag_direction();
    test_orbit_decay();
    test_bc_scaling();
    test_sun_position();

    printf("\n============================================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("============================================================\n");

    return tests_failed > 0 ? 1 : 0;
}
