/**
 * HPOP State Vector — Implementation
 *
 * The State, KeplerianElements, and SpacecraftProperties structs are POD types
 * defined in state.h. The conversion functions cartesian_to_keplerian() and
 * keplerian_to_cartesian() are implemented inline in state.h for performance.
 *
 * This translation unit provides additional state utility functions that are
 * not performance-critical and benefit from being in a separate TU to reduce
 * header bloat.
 *
 * References:
 *   - Vallado, "Fundamentals of Astrodynamics and Applications" (2013)
 *     Section 2.4: Classical Orbital Elements
 *     Algorithm 9: RV2COE (Cartesian to Keplerian)
 *     Algorithm 10: COE2RV (Keplerian to Cartesian)
 *   - Bate, Mueller, White, "Fundamentals of Astrodynamics" (1971)
 *     Section 2.5: Orbital Elements from Position and Velocity
 */

#include "hpop/state.h"
#include <cmath>
#include <algorithm>

namespace hpop {

// ── Geodetic Conversions ──

/**
 * Convert ECI position to geodetic coordinates (lat, lon, alt).
 * Uses iterative Bowring method for geodetic latitude.
 *
 * Reference: Bowring, B. R. (1976) "Transformation from spatial to
 * geographical coordinates" Survey Review 23(181), pp. 323-327
 *
 * @param r       ECI position vector (km)
 * @param jd      Julian Date (for GMST computation)
 * @param lat_rad Output geodetic latitude (rad)
 * @param lon_rad Output geodetic longitude (rad)
 * @param alt_km  Output altitude above WGS-84 ellipsoid (km)
 */
void eci_to_geodetic(const double r[3], double jd,
                     double& lat_rad, double& lon_rad, double& alt_km) {
    // WGS-84 ellipsoid parameters
    static constexpr double a = 6378.137;          // equatorial radius (km)
    static constexpr double f = 1.0 / 298.257223563;  // flattening
    static constexpr double b = a * (1.0 - f);     // polar radius (km)
    static constexpr double e2 = 1.0 - (b * b) / (a * a);  // eccentricity squared

    // Greenwich Mean Sidereal Time (GMST) from Julian Date
    // Vallado (2013), Algorithm 15
    double T = (jd - 2451545.0) / 36525.0;
    double gmst = 67310.54841 + (876600.0 * 3600.0 + 8640184.812866) * T
                  + 0.093104 * T * T - 6.2e-6 * T * T * T;
    gmst = std::fmod(gmst, 86400.0);
    if (gmst < 0) gmst += 86400.0;
    double theta_gmst = gmst / 86400.0 * TWO_PI;  // radians

    // Longitude
    lon_rad = std::atan2(r[1], r[0]) - theta_gmst;
    // Normalize to [-π, π]
    while (lon_rad > M_PI) lon_rad -= TWO_PI;
    while (lon_rad < -M_PI) lon_rad += TWO_PI;

    // Distance from Z axis
    double rho = std::sqrt(r[0] * r[0] + r[1] * r[1]);

    // Iterative Bowring method for geodetic latitude
    // Initial approximation using geocentric latitude
    double phi = std::atan2(r[2], rho);
    for (int i = 0; i < 5; ++i) {  // converges in 2-3 iterations
        double sinphi = std::sin(phi);
        double N = a / std::sqrt(1.0 - e2 * sinphi * sinphi);
        phi = std::atan2(r[2] + e2 * N * sinphi, rho);
    }
    lat_rad = phi;

    // Altitude
    double sinphi = std::sin(lat_rad);
    double N = a / std::sqrt(1.0 - e2 * sinphi * sinphi);
    if (std::abs(lat_rad) < M_PI / 4.0) {
        alt_km = rho / std::cos(lat_rad) - N;
    } else {
        alt_km = r[2] / std::sin(lat_rad) - N * (1.0 - e2);
    }
}

/**
 * Compute orbital period from semi-major axis.
 *
 * @param a   Semi-major axis (km)
 * @param mu  Gravitational parameter (km³/s²)
 * @return    Orbital period (seconds)
 */
double orbital_period(double a, double mu) {
    return TWO_PI * std::sqrt(a * a * a / mu);
}

/**
 * Compute specific orbital energy (vis-viva).
 *
 * @param r   Position magnitude (km)
 * @param v   Velocity magnitude (km/s)
 * @param mu  Gravitational parameter (km³/s²)
 * @return    Specific orbital energy (km²/s²)
 */
double specific_energy(double r, double v, double mu) {
    return v * v / 2.0 - mu / r;
}

/**
 * Convert true anomaly to eccentric anomaly.
 *
 * @param nu  True anomaly (rad)
 * @param e   Eccentricity
 * @return    Eccentric anomaly (rad)
 */
double true_to_eccentric_anomaly(double nu, double e) {
    double E = std::atan2(std::sqrt(1.0 - e * e) * std::sin(nu),
                          e + std::cos(nu));
    if (E < 0) E += TWO_PI;
    return E;
}

/**
 * Convert eccentric anomaly to mean anomaly (Kepler's equation).
 *
 * @param E   Eccentric anomaly (rad)
 * @param e   Eccentricity
 * @return    Mean anomaly (rad)
 */
double eccentric_to_mean_anomaly(double E, double e) {
    return E - e * std::sin(E);
}

/**
 * Solve Kepler's equation M = E - e sin(E) for E.
 * Uses Newton-Raphson iteration.
 *
 * Reference: Vallado (2013), Algorithm 2
 *
 * @param M   Mean anomaly (rad)
 * @param e   Eccentricity
 * @return    Eccentric anomaly (rad)
 */
double mean_to_eccentric_anomaly(double M, double e) {
    // Initial guess
    double E = (e < 0.8) ? M : M_PI;

    // Newton-Raphson iteration
    for (int i = 0; i < 50; ++i) {
        double dE = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
        E -= dE;
        if (std::abs(dE) < 1e-14) break;
    }
    return E;
}

/**
 * Convert mean anomaly to true anomaly via eccentric anomaly.
 *
 * @param M   Mean anomaly (rad)
 * @param e   Eccentricity
 * @return    True anomaly (rad)
 */
double mean_to_true_anomaly(double M, double e) {
    double E = mean_to_eccentric_anomaly(M, e);
    double nu = std::atan2(std::sqrt(1.0 - e * e) * std::sin(E),
                           std::cos(E) - e);
    if (nu < 0) nu += TWO_PI;
    return nu;
}

}  // namespace hpop
