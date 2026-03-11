#pragma once
/**
 * HPOP Third-Body Perturbations (Sun and Moon)
 *
 * Computes the gravitational acceleration from a third body using
 * Battin's formulation to avoid numerical cancellation:
 *
 *   a = μ₃ · (r₃_sat · q - r_sat / |r₃|³)
 *
 * Where q is computed to maintain precision when |r_sat| << |r₃|
 * (which is always the case for Sun perturbations, and nearly always for Moon).
 *
 * Ephemerides:
 *   - Sun: Low-precision analytical (Vallado, Algorithm 29)
 *   - Moon: Low-precision analytical (Vallado, Algorithm 31)
 *   Both accurate to ~0.01° in RA/Dec, sufficient for perturbation computation.
 *
 * References:
 *   - Vallado (2013), Section 8.6.4, Eq. 8-34
 *   - Battin, "An Introduction to the Mathematics and Methods of Astrodynamics" (1999)
 *   - Montenbruck & Gill (2000), Section 3.3
 */

#include "hpop/force_model.h"
#include <cmath>

namespace hpop {

// ── Lunar/Solar Ephemeris ──

/**
 * Compute Sun position in ECI (km) from Julian Date.
 * Low-precision algorithm from Vallado (2013), Algorithm 29.
 * Accuracy: ~0.01° in RA/Dec, ~0.01 AU in distance.
 */
void sun_position_eci(double jd, double r_sun[3]);

/**
 * Compute Moon position in ECI (km) from Julian Date.
 * Low-precision algorithm from Vallado (2013), Algorithm 31.
 * Accuracy: ~0.3° in ecliptic longitude, ~0.2° in latitude.
 */
void moon_position_eci(double jd, double r_moon[3]);

// ── Third-Body Force Model ──

enum class ThirdBody {
    SUN,
    MOON,
};

class ThirdBodyForce : public ForceModel {
public:
    explicit ThirdBodyForce(ThirdBody body) : body_(body) {}

    void acceleration(const double r[3], const double v[3],
                      double epoch_jd, const SpacecraftProperties* sc,
                      double a_out[3]) const override;

    const char* name() const override {
        return body_ == ThirdBody::SUN ? "ThirdBody_Sun" : "ThirdBody_Moon";
    }

private:
    ThirdBody body_;

    /// Battin's third-body acceleration formulation
    /// Avoids numerical cancellation for |r_sat| << |r_3body|
    static void battin_acceleration(
        const double r_sat[3],       // satellite position (km, ECI)
        const double r_3body[3],     // third body position (km, ECI)
        double mu_3body,             // gravitational parameter (km³/s²)
        double a_out[3]              // output acceleration (km/s²)
    );
};

// ── Convenience: Combined Sun+Moon Force ──
// Just add both ThirdBodyForce(SUN) and ThirdBodyForce(MOON) to the propagator

}  // namespace hpop
