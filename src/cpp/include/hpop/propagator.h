#pragma once
/**
 * HPOP Numerical Propagator
 *
 * Cowell's method: direct numerical integration of the equations of motion.
 * Integrators: RKF-45 (adaptive), RKF-78 (adaptive), RK4 (fixed step).
 *
 * Reference: Vallado, "Fundamentals of Astrodynamics and Applications" 5th ed.
 *            Montenbruck & Gill, "Satellite Orbits" Ch. 4
 */

#include "state.h"
#include "force_model.h"
#include <vector>
#include <memory>
#include <functional>

namespace hpop {

// ── Integrator Type ──

enum class IntegratorType {
    RK4,           // Classic 4th order Runge-Kutta (fixed step)
    RKF45,         // Runge-Kutta-Fehlberg 4(5) (adaptive step)
    RKF78,         // Runge-Kutta-Fehlberg 7(8) (adaptive step)
};

// ── Propagator Configuration ──

struct PropagatorConfig {
    IntegratorType integrator = IntegratorType::RKF45;
    double step_size     = 60.0;     // initial/fixed step size (seconds)
    double min_step      = 0.01;     // minimum step size (seconds)
    double max_step      = 300.0;    // maximum step size (seconds)
    double abs_tol       = 1e-10;    // absolute error tolerance (km or km/s)
    double rel_tol       = 1e-12;    // relative error tolerance
};

// ── Ephemeris Point ──

struct EphemerisPoint {
    double epoch_jd;   // Julian Date (TDB)
    State  state;      // position + velocity
};

// ── Propagator ──

class Propagator {
public:
    Propagator() = default;

    /// Add a force model (ownership transferred)
    void add_force(std::unique_ptr<ForceModel> force) {
        forces_.push_back(std::move(force));
    }

    /// Configure integrator settings
    void set_config(const PropagatorConfig& config) { config_ = config; }

    /// Set spacecraft properties
    void set_spacecraft(const SpacecraftProperties& sc) { sc_ = sc; }

    /// Propagate from initial state to target epoch
    /// Returns final state. Optionally records ephemeris at given step.
    State propagate(
        const State& initial,
        double epoch_jd0,         // initial epoch (JD)
        double epoch_jdf,         // final epoch (JD)
        double output_step = 0.0, // output ephemeris step (seconds), 0 = no output
        std::vector<EphemerisPoint>* ephemeris = nullptr
    ) const;

    /// Compute total acceleration from all force models
    void total_acceleration(
        const double r[3], const double v[3],
        double epoch_jd,
        double a_out[3]
    ) const;

private:
    std::vector<std::unique_ptr<ForceModel>> forces_;
    PropagatorConfig config_;
    SpacecraftProperties sc_;

    // ── Integrator Implementations ──

    /// RK4 fixed-step integration
    void rk4_step(
        double y[6], double epoch_jd, double dt
    ) const;

    /// RKF-45 adaptive step integration
    /// Returns actual step taken and error estimate
    double rkf45_step(
        double y[6], double epoch_jd, double dt,
        double& error
    ) const;

    /// RKF-78 adaptive step integration
    double rkf78_step(
        double y[6], double epoch_jd, double dt,
        double& error
    ) const;

    /// Equations of motion: dy/dt = f(t, y)
    void deriv(
        const double y[6], double epoch_jd,
        double dydt[6]
    ) const;
};

}  // namespace hpop
