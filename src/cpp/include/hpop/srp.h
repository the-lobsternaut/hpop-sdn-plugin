#pragma once
/**
 * HPOP Solar Radiation Pressure (SRP) Model
 *
 * Cannonball model: a_srp = -ν · P_sun · Cr · (A/m) · (AU/|r_sat_sun|)² · r̂_sat_sun
 *
 * Where:
 *   ν          = shadow function (0 = umbra, 0-1 = penumbra, 1 = sunlit)
 *   P_sun      = solar radiation pressure at 1 AU (4.56e-6 N/m²)
 *   Cr         = reflectivity coefficient (1.0 = absorb, 2.0 = perfect specular)
 *   A          = SRP cross-section area (m²)
 *   m          = spacecraft mass (kg)
 *   r_sat_sun  = vector from satellite to Sun (km)
 *   AU         = 1 astronomical unit (km)
 *
 * Shadow models:
 *   - None (always sunlit) — simplest, for testing
 *   - Cylindrical — hard on/off boundary, Earth casts a cylinder shadow
 *   - Conical (dual-cone) — penumbra + umbra, physically accurate
 *
 * The conical model computes the apparent overlap of Earth's disk
 * with the Sun's disk as seen from the satellite position.
 *
 * References:
 *   - Vallado (2013), Section 8.6.3, p. 526
 *   - Montenbruck & Gill (2000), Section 3.4, Algorithm 3.3
 *   - Montenbruck & Gill (2000), Section 3.4.2 (shadow function)
 */

#include "hpop/force_model.h"
#include "hpop/atmosphere.h"  // for HarrisPriester::sun_position
#include <cmath>
#include <memory>

namespace hpop {

// ── Shadow Model ──

enum class ShadowType {
    NONE,         // No shadow (always sunlit, ν = 1)
    CYLINDRICAL,  // Hard shadow cylinder behind Earth
    CONICAL,      // Dual-cone: penumbra + umbra (most accurate)
};

/**
 * Compute shadow function ν ∈ [0, 1].
 *
 * @param r_sat     Satellite position in ECI (km)
 * @param r_sun     Sun position in ECI (km)
 * @param type      Shadow model to use
 * @return ν: 0 = full shadow (umbra), 1 = full sunlight, (0,1) = penumbra
 */
double shadow_function(const double r_sat[3], const double r_sun[3], ShadowType type);

// ── SRP Force Model ──

class SRPForce : public ForceModel {
public:
    /**
     * @param shadow  Shadow model type (default: CONICAL for best accuracy)
     */
    explicit SRPForce(ShadowType shadow = ShadowType::CONICAL)
        : shadow_type_(shadow) {}

    void acceleration(const double r[3], const double v[3],
                      double epoch_jd, const SpacecraftProperties* sc,
                      double a_out[3]) const override;

    const char* name() const override { return "SRP"; }

    void set_shadow_type(ShadowType t) { shadow_type_ = t; }

    /// Compute Sun position in ECI (km) from Julian Date
    static void sun_position_eci(double jd, double r_sun[3]);

private:
    ShadowType shadow_type_;
};

}  // namespace hpop
