/**
 * HPOP Solar Radiation Pressure — Implementation
 *
 * Shadow function algorithm from Montenbruck & Gill (2000), Section 3.4.2.
 * The conical shadow model computes the apparent angular radii of the Sun
 * and Earth as seen from the satellite, then determines the fractional
 * illumination based on disk overlap geometry.
 *
 * SRP acceleration from Vallado (2013), Section 8.6.3.
 *
 * Unit conventions:
 *   - Positions in km (ECI J2000)
 *   - Accelerations in km/s²
 *   - P_sun in N/m² → converted to km/s² via area/mass
 */

#include "hpop/srp.h"

namespace hpop {

// ── Physical Constants ──

static constexpr double R_SUN_KM  = 696000.0;     // Solar radius (km)
static constexpr double R_EARTH_KM = RE_KM;        // Earth radius (km) — from state.h

// P_SUN is already defined in state.h as 4.56e-6 N/m² at 1 AU
// Units: a = P * Cr * A / m → [N/m² · m² / kg] = m/s² → /1000 for km/s²

// ── Shadow Function ──

static inline double vec_mag(const double v[3]) {
    return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

static inline double vec_dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

double shadow_function(const double r_sat[3], const double r_sun[3], ShadowType type) {
    if (type == ShadowType::NONE) return 1.0;

    // Vector from satellite to Sun
    double s[3] = {r_sun[0] - r_sat[0], r_sun[1] - r_sat[1], r_sun[2] - r_sat[2]};
    double s_mag = vec_mag(s);
    double r_mag = vec_mag(r_sat);

    if (type == ShadowType::CYLINDRICAL) {
        // Cylindrical shadow: satellite is in shadow if it's behind Earth
        // relative to the Sun direction.
        //
        // Project satellite position onto Sun direction.
        // If the projection is negative (satellite behind Earth from Sun's POV)
        // AND the perpendicular distance to the Sun line is less than Earth's radius,
        // then satellite is in shadow.
        double s_hat[3] = {s[0]/s_mag, s[1]/s_mag, s[2]/s_mag};

        // Component of r_sat along Sun direction (from satellite toward Sun)
        // Negative means satellite is on the anti-Sun side
        double proj = -(r_sat[0]*s_hat[0] + r_sat[1]*s_hat[1] + r_sat[2]*s_hat[2]);

        if (proj > 0) {
            // Satellite is between Earth and Sun or on Sun side — sunlit
            // Actually we need to check: is the angle between r_sat and r_sun > 90°?
            // If dot(r_sat, r_sun) > 0, satellite is on the Sun side of Earth
            double dot_rs = vec_dot(r_sat, r_sun);
            if (dot_rs > 0) return 1.0;  // Sun side
        }

        // Perpendicular distance from satellite to Earth-Sun line
        // d² = |r_sat|² - (r_sat · ŝ_sun)²
        double sun_hat[3] = {r_sun[0]/vec_mag(r_sun), r_sun[1]/vec_mag(r_sun), r_sun[2]/vec_mag(r_sun)};
        double r_dot_s = vec_dot(r_sat, sun_hat);
        double d_sq = r_mag * r_mag - r_dot_s * r_dot_s;

        if (r_dot_s < 0 && d_sq < R_EARTH_KM * R_EARTH_KM) {
            return 0.0;  // In cylindrical shadow
        }
        return 1.0;
    }

    // ── Conical Shadow (Montenbruck & Gill, Section 3.4.2) ──
    //
    // Compute apparent angular radii of Sun and Earth as seen from satellite,
    // then determine fractional illumination from disk overlap.

    double r_sun_mag = vec_mag(r_sun);

    // Apparent angular radius of the Sun as seen from the satellite
    double a_sun = std::asin(R_SUN_KM / s_mag);

    // Apparent angular radius of Earth as seen from the satellite
    double a_earth = std::asin(R_EARTH_KM / r_mag);

    // Angle between satellite-to-Sun and satellite-to-Earth-center directions
    // Earth center is at origin, so satellite-to-Earth is -r_sat
    // satellite-to-Sun is s = r_sun - r_sat
    double neg_r[3] = {-r_sat[0], -r_sat[1], -r_sat[2]};
    double cos_theta = vec_dot(neg_r, s) / (r_mag * s_mag);
    // Clamp for numerical safety
    cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
    double theta = std::acos(cos_theta);

    // Check illumination conditions
    // Case 1: No eclipse — Earth disk doesn't overlap Sun disk
    if (theta >= a_sun + a_earth) {
        return 1.0;  // Full sunlight
    }

    // Case 2: Total eclipse (umbra) — Earth disk completely covers Sun disk
    if (theta <= a_earth - a_sun && a_earth > a_sun) {
        return 0.0;  // Full shadow
    }

    // Case 3: Annular eclipse — Sun disk larger than Earth disk, but Earth centered in Sun
    // This is rare for satellites, but handle it:
    if (theta <= a_sun - a_earth && a_sun > a_earth) {
        // Earth covers a fraction of the Sun disk
        double ratio = (a_earth * a_earth) / (a_sun * a_sun);
        return 1.0 - ratio;
    }

    // Case 4: Partial eclipse (penumbra) — disk overlap
    // Use the area of intersection of two circles on the sky
    // A_overlap / A_sun gives the fraction blocked
    //
    // Area of intersection of two circles with radii a_sun, a_earth
    // separated by angle theta:
    //
    // A = a_sun² · arccos((theta² + a_sun² - a_earth²) / (2·theta·a_sun))
    //   + a_earth² · arccos((theta² + a_earth² - a_sun²) / (2·theta·a_earth))
    //   - 0.5·sqrt((-theta+a_sun+a_earth)(theta+a_sun-a_earth)(theta-a_sun+a_earth)(theta+a_sun+a_earth))

    double x = (theta*theta + a_sun*a_sun - a_earth*a_earth) / (2.0*theta*a_sun);
    double y = (theta*theta + a_earth*a_earth - a_sun*a_sun) / (2.0*theta*a_earth);

    // Clamp for numerical safety
    x = std::max(-1.0, std::min(1.0, x));
    y = std::max(-1.0, std::min(1.0, y));

    double term1 = a_sun * a_sun * std::acos(x);
    double term2 = a_earth * a_earth * std::acos(y);

    double s1 = (-theta + a_sun + a_earth);
    double s2 = ( theta + a_sun - a_earth);
    double s3 = ( theta - a_sun + a_earth);
    double s4 = ( theta + a_sun + a_earth);
    double term3 = 0.5 * std::sqrt(std::max(0.0, s1 * s2 * s3 * s4));

    double A_overlap = term1 + term2 - term3;
    double A_sun = M_PI * a_sun * a_sun;

    double nu = 1.0 - A_overlap / A_sun;
    return std::max(0.0, std::min(1.0, nu));
}

// ── Sun Position (ECI, km) ──

void SRPForce::sun_position_eci(double jd, double r_sun[3]) {
    // Low-precision Sun position in ECI (km)
    // Reuses the same algorithm as HarrisPriester but outputs Cartesian coordinates
    double T = (jd - 2451545.0) / 36525.0;

    // Mean anomaly
    double M = (357.5256 + 35999.049 * T) * DEG2RAD;

    // Ecliptic longitude
    double lambda = (280.460 + 36000.770 * T) * DEG2RAD
                    + 1.915 * DEG2RAD * std::sin(M)
                    + 0.020 * DEG2RAD * std::sin(2.0 * M);

    // Obliquity
    double epsilon = (23.4393 - 0.0130 * T) * DEG2RAD;

    // Distance (AU) — simplified
    double r_au = 1.00014 - 0.01671 * std::cos(M) - 0.00014 * std::cos(2.0 * M);
    double r_km = r_au * AU_KM;

    // ECI coordinates
    r_sun[0] = r_km * std::cos(lambda);
    r_sun[1] = r_km * std::sin(lambda) * std::cos(epsilon);
    r_sun[2] = r_km * std::sin(lambda) * std::sin(epsilon);
}

// ── SRP Acceleration ──

void SRPForce::acceleration(
    const double r[3], const double v[3],
    double epoch_jd, const SpacecraftProperties* sc,
    double a_out[3]) const
{
    double area = sc ? sc->srp_area : 22.0;
    double mass = sc ? sc->mass_kg : 260.0;
    double cr   = sc ? sc->cr : 1.3;

    // Sun position in ECI (km)
    double r_sun[3];
    sun_position_eci(epoch_jd, r_sun);

    // Shadow function
    double nu = shadow_function(r, r_sun, shadow_type_);
    if (nu < 1e-10) {
        a_out[0] = a_out[1] = a_out[2] = 0.0;
        return;
    }

    // Vector from satellite to Sun
    double r_sat_sun[3] = {
        r_sun[0] - r[0],
        r_sun[1] - r[1],
        r_sun[2] - r[2]
    };
    double dist = vec_mag(r_sat_sun);

    // SRP acceleration (cannonball model):
    // a = -ν · P_sun · Cr · (A/m) · (AU/dist)² · (r_sat_sun / dist)
    //
    // P_sun is in N/m², A in m², m in kg → gives m/s²
    // Divide by 1000 to get km/s²
    // The (AU/dist) factor accounts for the inverse-square law
    // But since dist ≈ AU (satellite is near Earth), this is ≈ 1

    double au_factor = (AU_KM / dist) * (AU_KM / dist);
    double coeff = -nu * P_SUN * cr * (area / mass) * au_factor / dist / 1000.0;

    a_out[0] = coeff * r_sat_sun[0];
    a_out[1] = coeff * r_sat_sun[1];
    a_out[2] = coeff * r_sat_sun[2];
}

}  // namespace hpop
