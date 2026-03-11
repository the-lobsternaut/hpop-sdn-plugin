/**
 * HPOP Third-Body Perturbations — Implementation
 *
 * Sun and Moon ephemerides from Vallado (2013), Algorithms 29 and 31.
 * Third-body acceleration uses Battin's formulation for numerical precision.
 *
 * The Battin q-factor avoids catastrophic cancellation in the standard
 * formulation when |r_sat| << |r_3body| (always true for Sun, usually for Moon).
 *
 * Standard formulation:  a = μ₃ · (r₃sat/|r₃sat|³ - r₃/|r₃|³)
 * Battin formulation:    a = μ₃ · (r₃sat · q - r_sat/|r₃|³)
 *   where q = f(r_sat, r₃sat, r₃) avoids the near-cancellation
 */

#include "hpop/thirdbody.h"

namespace hpop {

// ── Sun Position (Vallado Algorithm 29) ──

void sun_position_eci(double jd, double r_sun[3]) {
    double T = (jd - 2451545.0) / 36525.0;

    // Mean longitude (deg)
    double meanlong = std::fmod(280.460 + 36000.77 * T, 360.0);
    if (meanlong < 0) meanlong += 360.0;

    // Mean anomaly (rad)
    double M = std::fmod(357.5277233 + 35999.05034 * T, 360.0) * DEG2RAD;
    if (M < 0) M += TWO_PI;

    // Ecliptic longitude (deg → rad)
    double eclplong = (meanlong + 1.914666471 * std::sin(M)
                       + 0.019994643 * std::sin(2.0 * M)) * DEG2RAD;

    // Obliquity (deg → rad)
    double obliquity = (23.439291 - 0.0130042 * T) * DEG2RAD;

    // Distance in AU
    double r_au = 1.000140612 - 0.016708617 * std::cos(M)
                  - 0.000139589 * std::cos(2.0 * M);

    // Convert to km
    double r_km = r_au * AU_KM;

    // ECI components
    r_sun[0] = r_km * std::cos(eclplong);
    r_sun[1] = r_km * std::cos(obliquity) * std::sin(eclplong);
    r_sun[2] = r_km * std::sin(obliquity) * std::sin(eclplong);
}

// ── Moon Position (Vallado Algorithm 31) ──

void moon_position_eci(double jd, double r_moon[3]) {
    double T = (jd - 2451545.0) / 36525.0;

    // Ecliptic longitude (deg)
    double eclplong = 218.32 + 481267.8813 * T
        + 6.29 * std::sin((134.9 + 477198.85 * T) * DEG2RAD)
        - 1.27 * std::sin((259.2 - 413335.38 * T) * DEG2RAD)
        + 0.66 * std::sin((235.7 + 890534.23 * T) * DEG2RAD)
        + 0.21 * std::sin((269.9 + 954397.70 * T) * DEG2RAD)
        - 0.19 * std::sin((357.5 + 35999.05  * T) * DEG2RAD)
        - 0.11 * std::sin((186.6 + 966404.05 * T) * DEG2RAD);

    // Ecliptic latitude (deg)
    double eclplat = 5.13 * std::sin((93.3  + 483202.03 * T) * DEG2RAD)
        + 0.28 * std::sin((228.2 + 960400.87 * T) * DEG2RAD)
        - 0.28 * std::sin((318.3 + 6003.18   * T) * DEG2RAD)
        - 0.17 * std::sin((217.6 - 407332.20 * T) * DEG2RAD);

    // Horizontal parallax (deg)
    double hzparal = 0.9508
        + 0.0518 * std::cos((134.9 + 477198.85 * T) * DEG2RAD)
        + 0.0095 * std::cos((259.2 - 413335.38 * T) * DEG2RAD)
        + 0.0078 * std::cos((235.7 + 890534.23 * T) * DEG2RAD)
        + 0.0028 * std::cos((269.9 + 954397.70 * T) * DEG2RAD);

    // Convert to radians
    eclplong *= DEG2RAD;
    eclplat  *= DEG2RAD;
    hzparal  *= DEG2RAD;

    // Obliquity (rad)
    double obliquity = (23.439291 - 0.0130042 * T) * DEG2RAD;

    // Direction cosines
    double l = std::cos(eclplat) * std::cos(eclplong);
    double m = std::cos(obliquity) * std::cos(eclplat) * std::sin(eclplong)
             - std::sin(obliquity) * std::sin(eclplat);
    double n = std::sin(obliquity) * std::cos(eclplat) * std::sin(eclplong)
             + std::cos(obliquity) * std::sin(eclplat);

    // Distance (km) — from horizontal parallax
    // sin(parallax) = R_earth / distance
    double dist_km = RE_KM / std::sin(hzparal);

    // ECI components
    r_moon[0] = dist_km * l;
    r_moon[1] = dist_km * m;
    r_moon[2] = dist_km * n;
}

// ── Battin's Third-Body Acceleration ──
// Reference: Vallado (2013), Eq. 8-34
// Avoids catastrophic cancellation when |r_sat| << |r_3body|

void ThirdBodyForce::battin_acceleration(
    const double r_sat[3],
    const double r_3body[3],
    double mu_3body,
    double a_out[3])
{
    // r₃_sat = r_3body - r_sat (vector from satellite to third body)
    double r3sat[3] = {
        r_3body[0] - r_sat[0],
        r_3body[1] - r_sat[1],
        r_3body[2] - r_sat[2]
    };

    double rmag_sq = r_sat[0]*r_sat[0] + r_sat[1]*r_sat[1] + r_sat[2]*r_sat[2];
    double rmag = std::sqrt(rmag_sq);

    double r3mag = std::sqrt(r_3body[0]*r_3body[0] + r_3body[1]*r_3body[1] + r_3body[2]*r_3body[2]);
    double r3satmag = std::sqrt(r3sat[0]*r3sat[0] + r3sat[1]*r3sat[1] + r3sat[2]*r3sat[2]);

    // Battin's q factor
    // dot(r_sat, r3sat) — note: r3sat = r_3body - r_sat
    double dot_r_r3sat = r_sat[0]*r3sat[0] + r_sat[1]*r3sat[1] + r_sat[2]*r3sat[2];

    double q = (rmag_sq + 2.0 * dot_r_r3sat)
             * (r3mag * r3mag + r3mag * r3satmag + r3satmag * r3satmag)
             / (r3mag * r3mag * r3mag * r3satmag * r3satmag * r3satmag * (r3mag + r3satmag));

    // a = μ₃ · (r3sat · q - r_sat / |r₃|³)
    double inv_r3_cubed = 1.0 / (r3mag * r3mag * r3mag);

    a_out[0] = mu_3body * (r3sat[0] * q - r_sat[0] * inv_r3_cubed);
    a_out[1] = mu_3body * (r3sat[1] * q - r_sat[1] * inv_r3_cubed);
    a_out[2] = mu_3body * (r3sat[2] * q - r_sat[2] * inv_r3_cubed);
}

// ── ThirdBodyForce::acceleration ──

void ThirdBodyForce::acceleration(
    const double r[3], const double v[3],
    double epoch_jd, const SpacecraftProperties* sc,
    double a_out[3]) const
{
    double r_body[3];
    double mu;

    if (body_ == ThirdBody::SUN) {
        sun_position_eci(epoch_jd, r_body);
        mu = MU_SUN;
    } else {
        moon_position_eci(epoch_jd, r_body);
        mu = MU_MOON;
    }

    battin_acceleration(r, r_body, mu, a_out);
}

}  // namespace hpop
