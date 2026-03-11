/**
 * HPOP Two-Body Propagation Tests
 *
 * Validates integrators against Kepler's equation (analytical solution).
 * Test cases derived from standard orbital mechanics references.
 *
 * Test strategy:
 *   1. Circular LEO orbit — period should match T = 2π√(a³/μ)
 *   2. Elliptical orbit — energy and angular momentum conservation
 *   3. Multi-revolution — accumulated error over 10 orbits
 *   4. Integrator comparison — RK4 vs RKF-45 vs RKF-78
 *   5. J2 secular rates — RAAN and argument of perigee drift
 */

#include "hpop/propagator.h"
#include <cstdio>
#include <cmath>
#include <cassert>

using namespace hpop;

// ── Kepler Propagation (analytical reference) ──

static State kepler_propagate(const State& s0, double dt, double mu = MU_EARTH) {
    auto kep = cartesian_to_keplerian(s0, mu);

    // Mean motion
    double n = std::sqrt(mu / (kep.a * kep.a * kep.a));

    // True anomaly → eccentric anomaly
    double E0 = std::atan2(std::sqrt(1.0 - kep.e*kep.e) * std::sin(kep.nu),
                           kep.e + std::cos(kep.nu));
    double M0 = E0 - kep.e * std::sin(E0);

    // Propagate mean anomaly
    double M = M0 + n * dt;

    // Solve Kepler's equation: M = E - e*sin(E)
    double E = M;
    for (int i = 0; i < 30; i++) {
        double dE = (M - E + kep.e * std::sin(E)) / (1.0 - kep.e * std::cos(E));
        E += dE;
        if (std::abs(dE) < 1e-15) break;
    }

    // Eccentric anomaly → true anomaly
    double nu_new = std::atan2(std::sqrt(1.0 - kep.e*kep.e) * std::sin(E),
                               std::cos(E) - kep.e);
    if (nu_new < 0) nu_new += TWO_PI;

    kep.nu = nu_new;
    return keplerian_to_cartesian(kep, mu);
}

static double position_error(const State& a, const State& b) {
    double dx = a.r[0]-b.r[0], dy = a.r[1]-b.r[1], dz = a.r[2]-b.r[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

static double velocity_error(const State& a, const State& b) {
    double dx = a.v[0]-b.v[0], dy = a.v[1]-b.v[1], dz = a.v[2]-b.v[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

static double orbital_energy(const State& s, double mu = MU_EARTH) {
    return s.vmag()*s.vmag() / 2.0 - mu / s.rmag();
}

static double angular_momentum_mag(const State& s) {
    double h[3] = {
        s.r[1]*s.v[2] - s.r[2]*s.v[1],
        s.r[2]*s.v[0] - s.r[0]*s.v[2],
        s.r[0]*s.v[1] - s.r[1]*s.v[0]
    };
    return std::sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
}

// ── Test Helpers ──

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (!(cond)) { \
        printf("  FAIL: %s\n", msg); \
        tests_failed++; \
    } else { \
        printf("  PASS: %s\n", msg); \
        tests_passed++; \
    } \
} while(0)

#define CHECK_TOL(val, expected, tol, msg) do { \
    double _v = (val), _e = (expected), _t = (tol); \
    if (std::abs(_v - _e) > _t) { \
        printf("  FAIL: %s (got %.10e, expected %.10e, diff %.3e, tol %.3e)\n", msg, _v, _e, std::abs(_v-_e), _t); \
        tests_failed++; \
    } else { \
        printf("  PASS: %s (diff %.3e < tol %.3e)\n", msg, std::abs(_v-_e), _t); \
        tests_passed++; \
    } \
} while(0)

// ── Test 1: Circular LEO (ISS-like) ──

void test_circular_leo() {
    printf("\n=== Test 1: Circular LEO (ISS-like, 400 km) ===\n");

    // ISS-like: a = 6778 km, e ≈ 0, i = 51.6°
    KeplerianElements kep = {6778.0, 0.0001, 51.6 * DEG2RAD, 0.0, 0.0, 0.0};
    State s0 = keplerian_to_cartesian(kep);

    double period = TWO_PI * std::sqrt(kep.a * kep.a * kep.a / MU_EARTH);
    printf("  Orbital period: %.2f s (%.2f min)\n", period, period / 60.0);

    double epoch_jd0 = 2460000.5;  // arbitrary epoch
    double epoch_jdf = epoch_jd0 + period / SEC_PER_DAY;

    // Analytical reference: after one period, should return to initial state
    State s_ref = kepler_propagate(s0, period);

    // Test each integrator
    for (auto integ : {IntegratorType::RK4, IntegratorType::RKF45, IntegratorType::RKF78}) {
        const char* name = (integ == IntegratorType::RK4) ? "RK4" :
                           (integ == IntegratorType::RKF45) ? "RKF45" : "RKF78";

        Propagator prop;
        prop.add_force(std::make_unique<TwoBodyForce>());
        PropagatorConfig cfg;
        cfg.integrator = integ;
        cfg.step_size = 10.0;
        cfg.abs_tol = 1e-12;
        cfg.rel_tol = 1e-14;
        prop.set_config(cfg);

        State sf = prop.propagate(s0, epoch_jd0, epoch_jdf);

        double r_err = position_error(sf, s_ref);
        double v_err = velocity_error(sf, s_ref);

        printf("  %s: pos_err=%.3e km, vel_err=%.3e km/s\n", name, r_err, v_err);

        // RK4 at 10s fixed step: ~0.06 km for 1 LEO orbit
        // RKF45/78 adaptive: ~0.004 km / ~0.0005 km respectively
        double tol_r = (integ == IntegratorType::RK4) ? 1e-1 : 1e-2;
        double tol_v = (integ == IntegratorType::RK4) ? 1e-4 : 1e-5;

        char buf[128];
        snprintf(buf, sizeof(buf), "%s position error < %.0e km (1 orbit)", name, tol_r);
        CHECK(r_err < tol_r, buf);
        snprintf(buf, sizeof(buf), "%s velocity error < %.0e km/s (1 orbit)", name, tol_v);
        CHECK(v_err < tol_v, buf);
    }
}

// ── Test 2: Elliptical Orbit — Energy Conservation ──

void test_elliptical_energy() {
    printf("\n=== Test 2: Elliptical Orbit — Energy & Angular Momentum Conservation ===\n");

    // Molniya-like: a = 26600 km, e = 0.74, i = 63.4°
    KeplerianElements kep = {26600.0, 0.74, 63.4 * DEG2RAD, 30.0 * DEG2RAD, 270.0 * DEG2RAD, 0.0};
    State s0 = keplerian_to_cartesian(kep);

    double period = TWO_PI * std::sqrt(kep.a * kep.a * kep.a / MU_EARTH);
    double epoch_jd0 = 2460000.5;

    double E0 = orbital_energy(s0);
    double h0 = angular_momentum_mag(s0);
    printf("  Initial energy: %.10f km²/s²\n", E0);
    printf("  Initial |h|: %.10f km²/s\n", h0);

    Propagator prop;
    prop.add_force(std::make_unique<TwoBodyForce>());
    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF78;
    cfg.step_size = 30.0;
    cfg.abs_tol = 1e-12;
    cfg.rel_tol = 1e-14;
    prop.set_config(cfg);

    // Propagate 10 orbits
    double tf = epoch_jd0 + 10.0 * period / SEC_PER_DAY;
    State sf = prop.propagate(s0, epoch_jd0, tf);

    double Ef = orbital_energy(sf);
    double hf = angular_momentum_mag(sf);

    double dE = std::abs(Ef - E0) / std::abs(E0);
    double dh = std::abs(hf - h0) / h0;

    printf("  After 10 orbits: dE/E=%.3e, dh/h=%.3e\n", dE, dh);

    CHECK(dE < 1e-10, "Energy conservation < 1e-10 (10 orbits, RKF78)");
    CHECK(dh < 1e-10, "Angular momentum conservation < 1e-10 (10 orbits, RKF78)");

    // Position accuracy vs Kepler
    State s_ref = kepler_propagate(s0, 10.0 * period);
    double r_err = position_error(sf, s_ref);
    printf("  Position error vs Kepler: %.3e km\n", r_err);
    CHECK(r_err < 1.0, "Position vs Kepler < 1.0 km (10 orbits, RKF78)");
}

// ── Test 3: Starlink-like Orbit (550 km) ──

void test_starlink_orbit() {
    printf("\n=== Test 3: Starlink-like Orbit (550 km, circular) ===\n");

    KeplerianElements kep = {6921.0, 0.0001, 53.0 * DEG2RAD, 0.0, 0.0, 0.0};
    State s0 = keplerian_to_cartesian(kep);

    double period = TWO_PI * std::sqrt(kep.a * kep.a * kep.a / MU_EARTH);
    double epoch_jd0 = 2460000.5;

    Propagator prop;
    prop.add_force(std::make_unique<TwoBodyForce>());
    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF45;
    cfg.step_size = 30.0;
    cfg.abs_tol = 1e-12;
    cfg.rel_tol = 1e-14;
    prop.set_config(cfg);

    // 1 day propagation
    double tf = epoch_jd0 + 1.0;
    State sf = prop.propagate(s0, epoch_jd0, tf);

    double n_orbits = SEC_PER_DAY / period;
    State s_ref = kepler_propagate(s0, SEC_PER_DAY);
    double r_err = position_error(sf, s_ref);

    printf("  %.1f orbits/day, period=%.1f min\n", n_orbits, period/60.0);
    printf("  1-day position error vs Kepler: %.3e km\n", r_err);
    CHECK(r_err < 1e-1, "1-day position error < 0.1 km (RKF45, two-body)");
}

// ── Test 4: J2 Secular Rates ──

void test_j2_secular() {
    printf("\n=== Test 4: J2 Secular Rates (RAAN and ω drift) ===\n");

    // Sun-synchronous orbit: a = 7078 km (700 km alt), i = 98.19°
    // Expected RAAN rate ≈ 0.9856 deg/day for sun-sync
    KeplerianElements kep = {7078.0, 0.001, 98.19 * DEG2RAD, 0.0, 0.0, 0.0};
    State s0 = keplerian_to_cartesian(kep);

    double epoch_jd0 = 2460000.5;
    double epoch_jdf = epoch_jd0 + 1.0;  // 1 day

    Propagator prop;
    prop.add_force(std::make_unique<TwoBodyForce>());
    prop.add_force(std::make_unique<J2Force>());
    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF78;
    cfg.step_size = 30.0;
    cfg.abs_tol = 1e-12;
    cfg.rel_tol = 1e-14;
    prop.set_config(cfg);

    State sf = prop.propagate(s0, epoch_jd0, epoch_jdf);

    // Compute RAAN drift analytically (Vallado eq 9-40)
    double n = std::sqrt(MU_EARTH / (kep.a * kep.a * kep.a));  // rad/s
    double p = kep.a * (1.0 - kep.e * kep.e);
    double raan_dot_analytical = -1.5 * n * J2 * (RE_KM / p) * (RE_KM / p) * std::cos(kep.i);  // rad/s
    double raan_drift_expected = raan_dot_analytical * SEC_PER_DAY * RAD2DEG;  // deg/day

    // Extract RAAN from propagated state
    auto kep_f = cartesian_to_keplerian(sf);
    double raan_drift_numerical = (kep_f.raan - kep.raan) * RAD2DEG;  // degrees

    // Wrap to [-180, 180]
    while (raan_drift_numerical > 180.0) raan_drift_numerical -= 360.0;
    while (raan_drift_numerical < -180.0) raan_drift_numerical += 360.0;

    printf("  Expected RAAN rate: %.4f deg/day\n", raan_drift_expected);
    printf("  Numerical RAAN drift: %.4f deg/day\n", raan_drift_numerical);
    printf("  Difference: %.4f deg/day\n", std::abs(raan_drift_numerical - raan_drift_expected));

    // J2 secular should agree within ~5% (oscillatory terms cause differences)
    double rel_error = std::abs(raan_drift_numerical - raan_drift_expected) / std::abs(raan_drift_expected);
    printf("  Relative error: %.2f%%\n", rel_error * 100.0);
    CHECK(rel_error < 0.05, "RAAN drift within 5% of analytical J2 secular rate");
}

// ── Test 5: Backward Propagation ──

void test_backward() {
    printf("\n=== Test 5: Backward Propagation (time reversibility) ===\n");

    KeplerianElements kep = {7500.0, 0.1, 45.0 * DEG2RAD, 30.0 * DEG2RAD, 60.0 * DEG2RAD, 90.0 * DEG2RAD};
    State s0 = keplerian_to_cartesian(kep);

    double epoch_jd0 = 2460000.5;
    double epoch_jdf = epoch_jd0 + 0.5;  // 12 hours forward

    Propagator prop;
    prop.add_force(std::make_unique<TwoBodyForce>());
    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF45;
    cfg.step_size = 30.0;
    cfg.abs_tol = 1e-12;
    cfg.rel_tol = 1e-14;
    prop.set_config(cfg);

    // Forward
    State sf = prop.propagate(s0, epoch_jd0, epoch_jdf);
    // Backward
    State sb = prop.propagate(sf, epoch_jdf, epoch_jd0);

    double r_err = position_error(sb, s0);
    double v_err = velocity_error(sb, s0);

    printf("  Round-trip pos error: %.3e km\n", r_err);
    printf("  Round-trip vel error: %.3e km/s\n", v_err);
    CHECK(r_err < 1e-1, "Round-trip position error < 0.1 km");
    CHECK(v_err < 1e-4, "Round-trip velocity error < 1e-4 km/s");
}

// ── Main ──

int main() {
    printf("============================================================\n");
    printf("HPOP Two-Body Propagation Tests\n");
    printf("============================================================\n");

    test_circular_leo();
    test_elliptical_energy();
    test_starlink_orbit();
    test_j2_secular();
    test_backward();

    printf("\n============================================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("============================================================\n");

    return tests_failed > 0 ? 1 : 0;
}
