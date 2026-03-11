#pragma once
/**
 * HPOP Force Model Interface
 *
 * All force models compute acceleration in J2000 ECI frame (km/s²).
 * Force models are composable — the propagator sums all accelerations.
 */

#include "state.h"

namespace hpop {

// ── Force Model Interface ──

class ForceModel {
public:
    virtual ~ForceModel() = default;

    /// Compute acceleration (km/s²) given position, velocity, and epoch (JD)
    virtual void acceleration(
        const double r[3],       // position (km), J2000 ECI
        const double v[3],       // velocity (km/s), J2000 ECI
        double epoch_jd,         // epoch (Julian Date, TDB)
        const SpacecraftProperties* sc,  // spacecraft properties (may be null)
        double a_out[3]          // output acceleration (km/s²)
    ) const = 0;

    /// Human-readable name for logging
    virtual const char* name() const = 0;
};

// ── Two-Body (Point Mass Gravity) ──

class TwoBodyForce : public ForceModel {
public:
    explicit TwoBodyForce(double mu = MU_EARTH) : mu_(mu) {}

    void acceleration(
        const double r[3], const double v[3],
        double epoch_jd, const SpacecraftProperties* sc,
        double a_out[3]) const override
    {
        double rmag = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double r3 = rmag * rmag * rmag;
        double coeff = -mu_ / r3;
        a_out[0] = coeff * r[0];
        a_out[1] = coeff * r[1];
        a_out[2] = coeff * r[2];
    }

    const char* name() const override { return "TwoBody"; }

private:
    double mu_;
};

// ── J2 Perturbation (zonal harmonic) ──

class J2Force : public ForceModel {
public:
    J2Force(double mu = MU_EARTH, double re = RE_KM, double j2 = J2)
        : mu_(mu), re_(re), j2_(j2) {}

    void acceleration(
        const double r[3], const double v[3],
        double epoch_jd, const SpacecraftProperties* sc,
        double a_out[3]) const override
    {
        double rmag = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
        double r2 = rmag * rmag;
        double r5 = r2 * r2 * rmag;
        double re2 = re_ * re_;

        double z2_r2 = (r[2] * r[2]) / r2;
        double coeff = -1.5 * j2_ * mu_ * re2 / r5;

        // J2 acceleration (Vallado eq 8-22)
        a_out[0] = coeff * r[0] * (1.0 - 5.0 * z2_r2);
        a_out[1] = coeff * r[1] * (1.0 - 5.0 * z2_r2);
        a_out[2] = coeff * r[2] * (3.0 - 5.0 * z2_r2);
    }

    const char* name() const override { return "J2"; }

private:
    double mu_, re_, j2_;
};

}  // namespace hpop
