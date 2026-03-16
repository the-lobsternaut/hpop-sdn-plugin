/**
 * HPOP Numerical Propagator — Cowell's Method
 *
 * Integrators:
 *   - RK4:   Classic 4th-order Runge-Kutta (fixed step)
 *   - RKF45: Runge-Kutta-Fehlberg 4(5) (adaptive step)
 *   - RKF78: Runge-Kutta-Fehlberg 7(8) (adaptive step)
 *
 * Butcher tableaux from:
 *   - Fehlberg, E. (1969) "Low-order classical Runge-Kutta formulas with stepsize
 *     control and their application to some heat transfer problems"
 *   - Fehlberg, E. (1968) "Classical fifth-, sixth-, seventh-, and eighth-order
 *     Runge-Kutta formulas with stepsize control" NASA TR R-287
 */

#include "hpop/propagator.h"
#include <cmath>
#include <algorithm>
#include <cstring>

namespace hpop {

// ── Equations of Motion ──

void Propagator::deriv(const double y[6], double epoch_jd, double dydt[6]) const {
    // dy/dt = [v, a]
    // Position derivative = velocity
    dydt[0] = y[3];
    dydt[1] = y[4];
    dydt[2] = y[5];

    // Velocity derivative = sum of accelerations
    double a_total[3] = {0, 0, 0};
    total_acceleration(y, y + 3, epoch_jd, a_total);
    dydt[3] = a_total[0];
    dydt[4] = a_total[1];
    dydt[5] = a_total[2];
}

void Propagator::total_acceleration(
    const double r[3], const double v[3],
    double epoch_jd, double a_out[3]) const
{
    a_out[0] = a_out[1] = a_out[2] = 0.0;
    double a_force[3];
    for (const auto& force : forces_) {
        force->acceleration(r, v, epoch_jd, &sc_, a_force);
        a_out[0] += a_force[0];
        a_out[1] += a_force[1];
        a_out[2] += a_force[2];
    }
}

// ── RK4 Fixed Step ──

void Propagator::rk4_step(double y[6], double epoch_jd, double dt) const {
    double k1[6], k2[6], k3[6], k4[6];
    double ytmp[6];
    double dt_day = dt / SEC_PER_DAY;

    // k1
    deriv(y, epoch_jd, k1);

    // k2
    for (int i = 0; i < 6; i++) ytmp[i] = y[i] + 0.5 * dt * k1[i];
    deriv(ytmp, epoch_jd + 0.5 * dt_day, k2);

    // k3
    for (int i = 0; i < 6; i++) ytmp[i] = y[i] + 0.5 * dt * k2[i];
    deriv(ytmp, epoch_jd + 0.5 * dt_day, k3);

    // k4
    for (int i = 0; i < 6; i++) ytmp[i] = y[i] + dt * k3[i];
    deriv(ytmp, epoch_jd + dt_day, k4);

    // Update
    for (int i = 0; i < 6; i++)
        y[i] += dt * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]) / 6.0;
}

// ── RKF-45 Adaptive Step ──
// Fehlberg 4(5) coefficients

double Propagator::rkf45_step(double y[6], double epoch_jd, double dt, double& error) const {
    // Butcher tableau for RKF-45
    static constexpr double a2  = 1.0/4.0;
    static constexpr double a3  = 3.0/8.0;
    static constexpr double a4  = 12.0/13.0;
    static constexpr double a5  = 1.0;
    static constexpr double a6  = 1.0/2.0;

    static constexpr double b21 = 1.0/4.0;
    static constexpr double b31 = 3.0/32.0,    b32 = 9.0/32.0;
    static constexpr double b41 = 1932.0/2197.0, b42 = -7200.0/2197.0, b43 = 7296.0/2197.0;
    static constexpr double b51 = 439.0/216.0,   b52 = -8.0,           b53 = 3680.0/513.0,   b54 = -845.0/4104.0;
    static constexpr double b61 = -8.0/27.0,     b62 = 2.0,            b63 = -3544.0/2565.0, b64 = 1859.0/4104.0, b65 = -11.0/40.0;

    // 4th order weights
    static constexpr double c1  = 25.0/216.0, c3 = 1408.0/2565.0, c4 = 2197.0/4104.0, c5 = -1.0/5.0;
    // 5th order weights (for error estimate)
    static constexpr double d1  = 16.0/135.0, d3 = 6656.0/12825.0, d4 = 28561.0/56430.0, d5 = -9.0/50.0, d6 = 2.0/55.0;

    double k1[6], k2[6], k3[6], k4[6], k5[6], k6[6];
    double ytmp[6];
    double dt_day = dt / SEC_PER_DAY;

    // k1
    deriv(y, epoch_jd, k1);

    // k2
    for (int i = 0; i < 6; i++) ytmp[i] = y[i] + dt * b21 * k1[i];
    deriv(ytmp, epoch_jd + a2 * dt_day, k2);

    // k3
    for (int i = 0; i < 6; i++) ytmp[i] = y[i] + dt * (b31*k1[i] + b32*k2[i]);
    deriv(ytmp, epoch_jd + a3 * dt_day, k3);

    // k4
    for (int i = 0; i < 6; i++) ytmp[i] = y[i] + dt * (b41*k1[i] + b42*k2[i] + b43*k3[i]);
    deriv(ytmp, epoch_jd + a4 * dt_day, k4);

    // k5
    for (int i = 0; i < 6; i++) ytmp[i] = y[i] + dt * (b51*k1[i] + b52*k2[i] + b53*k3[i] + b54*k4[i]);
    deriv(ytmp, epoch_jd + a5 * dt_day, k5);

    // k6
    for (int i = 0; i < 6; i++) ytmp[i] = y[i] + dt * (b61*k1[i] + b62*k2[i] + b63*k3[i] + b64*k4[i] + b65*k5[i]);
    deriv(ytmp, epoch_jd + a6 * dt_day, k6);

    // 5th order solution (used for propagation)
    double y5[6];
    for (int i = 0; i < 6; i++)
        y5[i] = y[i] + dt * (d1*k1[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i]);

    // 4th order solution (for error estimation)
    double y4[6];
    for (int i = 0; i < 6; i++)
        y4[i] = y[i] + dt * (c1*k1[i] + c3*k3[i] + c4*k4[i] + c5*k5[i]);

    // Error estimate: max |y5 - y4|
    error = 0.0;
    for (int i = 0; i < 6; i++) {
        double scale = config_.abs_tol + config_.rel_tol * std::max(std::abs(y[i]), std::abs(y5[i]));
        double ei = std::abs(y5[i] - y4[i]) / scale;
        error = std::max(error, ei);
    }

    // Accept step: copy 5th order result
    std::memcpy(y, y5, sizeof(y5));
    return dt;
}

// ── RKF-78 Adaptive Step ──
// Fehlberg 7(8) coefficients from NASA TR R-287

double Propagator::rkf78_step(double y[6], double epoch_jd, double dt, double& error) const {
    // RKF-78: 13 stages, 7th order solution, 8th order error estimate
    // Coefficients from Fehlberg (1968), Table III

    static constexpr int NSTAGE = 13;

    // Node points (c_i)
    static constexpr double c[NSTAGE] = {
        0.0, 2.0/27.0, 1.0/9.0, 1.0/6.0, 5.0/12.0, 1.0/2.0,
        5.0/6.0, 1.0/6.0, 2.0/3.0, 1.0/3.0, 1.0, 0.0, 1.0
    };

    // Coupling coefficients (a_ij) — lower triangular
    // Stage 1: no coupling (k1 uses initial state)
    // These are the standard Fehlberg 7(8) coefficients
    static constexpr double a[NSTAGE][12] = {
        {0},
        {2.0/27.0},
        {1.0/36.0, 1.0/12.0},
        {1.0/24.0, 0, 1.0/8.0},
        {5.0/12.0, 0, -25.0/16.0, 25.0/16.0},
        {1.0/20.0, 0, 0, 1.0/4.0, 1.0/5.0},
        {-25.0/108.0, 0, 0, 125.0/108.0, -65.0/27.0, 125.0/54.0},
        {31.0/300.0, 0, 0, 0, 61.0/225.0, -2.0/9.0, 13.0/900.0},
        {2.0, 0, 0, -53.0/6.0, 704.0/45.0, -107.0/9.0, 67.0/90.0, 3.0},
        {-91.0/108.0, 0, 0, 23.0/108.0, -976.0/135.0, 311.0/54.0, -19.0/60.0, 17.0/6.0, -1.0/12.0},
        {2383.0/4100.0, 0, 0, -341.0/164.0, 4496.0/1025.0, -301.0/82.0, 2133.0/4100.0, 45.0/82.0, 45.0/164.0, 18.0/41.0},
        {3.0/205.0, 0, 0, 0, 0, -6.0/41.0, -3.0/205.0, -3.0/41.0, 3.0/41.0, 6.0/41.0, 0},
        {-1777.0/4100.0, 0, 0, -341.0/164.0, 4496.0/1025.0, -289.0/82.0, 2193.0/4100.0, 51.0/82.0, 33.0/164.0, 12.0/41.0, 0, 1.0}
    };

    // 7th order weights
    static constexpr double b7[NSTAGE] = {
        41.0/840.0, 0, 0, 0, 0, 34.0/105.0, 9.0/35.0, 9.0/35.0,
        9.0/280.0, 9.0/280.0, 41.0/840.0, 0, 0
    };

    // 8th order weights
    static constexpr double b8[NSTAGE] = {
        0, 0, 0, 0, 0, 34.0/105.0, 9.0/35.0, 9.0/35.0,
        9.0/280.0, 9.0/280.0, 0, 41.0/840.0, 41.0/840.0
    };

    double k[NSTAGE][6];
    double ytmp[6];
    double dt_day = dt / SEC_PER_DAY;

    // Compute stages
    deriv(y, epoch_jd, k[0]);
    for (int s = 1; s < NSTAGE; s++) {
        for (int i = 0; i < 6; i++) {
            ytmp[i] = y[i];
            for (int j = 0; j < s; j++)
                ytmp[i] += dt * a[s][j] * k[j][i];
        }
        deriv(ytmp, epoch_jd + c[s] * dt_day, k[s]);
    }

    // 8th order solution (used for propagation)
    double y8[6];
    for (int i = 0; i < 6; i++) {
        y8[i] = y[i];
        for (int s = 0; s < NSTAGE; s++)
            y8[i] += dt * b8[s] * k[s][i];
    }

    // 7th order solution (for error estimation)
    double y7[6];
    for (int i = 0; i < 6; i++) {
        y7[i] = y[i];
        for (int s = 0; s < NSTAGE; s++)
            y7[i] += dt * b7[s] * k[s][i];
    }

    // Error estimate
    error = 0.0;
    for (int i = 0; i < 6; i++) {
        double scale = config_.abs_tol + config_.rel_tol * std::max(std::abs(y[i]), std::abs(y8[i]));
        double ei = std::abs(y8[i] - y7[i]) / scale;
        error = std::max(error, ei);
    }

    std::memcpy(y, y8, sizeof(y8));
    return dt;
}

// ── Main Propagation Loop ──

State Propagator::propagate(
    const State& initial,
    double epoch_jd0,
    double epoch_jdf,
    double output_step,
    std::vector<EphemerisPoint>* ephemeris) const
{
    // State vector: [x, y, z, vx, vy, vz]
    double y[6] = {
        initial.r[0], initial.r[1], initial.r[2],
        initial.v[0], initial.v[1], initial.v[2]
    };

    double t = epoch_jd0;
    double tf = epoch_jdf;
    double direction = (tf > t) ? 1.0 : -1.0;
    double dt = direction * std::abs(config_.step_size);

    // Output tracking
    double next_output = epoch_jd0;
    if (output_step > 0 && ephemeris) {
        ephemeris->clear();
        EphemerisPoint ep;
        ep.epoch_jd = t;
        std::memcpy(ep.state.r, y, 3 * sizeof(double));
        std::memcpy(ep.state.v, y + 3, 3 * sizeof(double));
        ephemeris->push_back(ep);
        next_output += output_step / SEC_PER_DAY;
    }

    while (direction * (tf - t) > 1e-15) {
        // Don't overshoot final time
        double remaining = (tf - t) * SEC_PER_DAY;
        if (std::abs(dt) > std::abs(remaining))
            dt = remaining;

        switch (config_.integrator) {
            case IntegratorType::RK4: {
                rk4_step(y, t, dt);
                t += dt / SEC_PER_DAY;
                break;
            }
            case IntegratorType::RKF45: {
                double y_save[6];
                std::memcpy(y_save, y, sizeof(y));
                double error;
                rkf45_step(y, t, dt, error);

                if (error > 1.0) {
                    // Reject step, reduce step size
                    std::memcpy(y, y_save, sizeof(y));
                    double factor = 0.9 * std::pow(error, -0.2);
                    factor = std::max(0.1, std::min(factor, 0.5));
                    dt *= factor;
                    dt = direction * std::max(std::abs(dt), config_.min_step);
                    continue;
                }

                t += dt / SEC_PER_DAY;

                // Adjust step size for next step
                double factor = 0.9 * std::pow(error + 1e-30, -0.2);
                factor = std::max(0.2, std::min(factor, 5.0));
                dt *= factor;
                dt = direction * std::clamp(std::abs(dt), config_.min_step, config_.max_step);
                break;
            }
            case IntegratorType::RKF78: {
                double y_save[6];
                std::memcpy(y_save, y, sizeof(y));
                double error;
                rkf78_step(y, t, dt, error);

                if (error > 1.0) {
                    std::memcpy(y, y_save, sizeof(y));
                    double factor = 0.9 * std::pow(error, -0.125);  // 1/8 for 8th order
                    factor = std::max(0.1, std::min(factor, 0.5));
                    dt *= factor;
                    dt = direction * std::max(std::abs(dt), config_.min_step);
                    continue;
                }

                t += dt / SEC_PER_DAY;

                double factor = 0.9 * std::pow(error + 1e-30, -0.125);
                factor = std::max(0.2, std::min(factor, 5.0));
                dt *= factor;
                dt = direction * std::clamp(std::abs(dt), config_.min_step, config_.max_step);
                break;
            }
        }

        // Output ephemeris points
        if (output_step > 0 && ephemeris) {
            while (direction * (t - next_output) >= -1e-15 && direction * (tf - next_output) >= -1e-15) {
                EphemerisPoint ep;
                ep.epoch_jd = next_output;
                // Taylor interpolation from current integrated state.
                // r(t+δ) ≈ r(t) + v(t)*δ + ½a(t)*δ²
                // v(t+δ) ≈ v(t) + a(t)*δ
                // where δ = (next_output - t) in seconds.
                // This is 2nd-order accurate, suitable for output
                // spacing << integration step.
                double delta_s = (next_output - t) * SEC_PER_DAY;
                // Compute acceleration at current state for interpolation
                double a_interp[3] = {0.0, 0.0, 0.0};
                double r_mag = std::sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
                double r3_inv = 1.0 / (r_mag * r_mag * r_mag);
                for (int j = 0; j < 3; j++) {
                    a_interp[j] = -MU_EARTH * y[j] * r3_inv;  // two-body term
                }
                for (int j = 0; j < 3; j++) {
                    ep.state.r[j] = y[j] + y[j+3] * delta_s
                                    + 0.5 * a_interp[j] * delta_s * delta_s;
                    ep.state.v[j] = y[j+3] + a_interp[j] * delta_s;
                }
                ephemeris->push_back(ep);
                next_output += output_step / SEC_PER_DAY;
            }
        }
    }

    State result;
    std::memcpy(result.r, y, 3 * sizeof(double));
    std::memcpy(result.v, y + 3, 3 * sizeof(double));
    return result;
}

}  // namespace hpop
