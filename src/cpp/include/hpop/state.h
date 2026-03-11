#pragma once
/**
 * HPOP State Vector Types
 *
 * Orbital state representations and conversions.
 * All internal computation in J2000 ECI frame, km and km/s.
 */

#include <array>
#include <cmath>

namespace hpop {

// ── Constants (WGS-84 / EGM-96) ──

static constexpr double MU_EARTH = 398600.4418;           // km³/s² (WGS-84)
static constexpr double RE_KM    = 6378.137;              // Earth equatorial radius (km, WGS-84)
static constexpr double J2       = 1.08262982e-3;         // J2 zonal harmonic
static constexpr double AU_KM    = 149597870.7;           // Astronomical unit (km)
static constexpr double C_LIGHT  = 299792.458;            // Speed of light (km/s)
static constexpr double P_SUN    = 4.56e-6;               // Solar radiation pressure at 1 AU (N/m² = kPa*1e3)
static constexpr double MU_SUN   = 1.32712440018e11;      // km³/s²
static constexpr double MU_MOON  = 4902.800066;           // km³/s²

static constexpr double DEG2RAD  = M_PI / 180.0;
static constexpr double RAD2DEG  = 180.0 / M_PI;
static constexpr double TWO_PI   = 2.0 * M_PI;
static constexpr double SEC_PER_DAY = 86400.0;

// ── State Vector (position + velocity, 6 elements) ──

struct State {
    double r[3];  // position (km), J2000 ECI
    double v[3];  // velocity (km/s), J2000 ECI

    double rmag() const {
        return std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    }
    double vmag() const {
        return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
    }
};

// ── Keplerian Elements ──

struct KeplerianElements {
    double a;      // semi-major axis (km)
    double e;      // eccentricity
    double i;      // inclination (rad)
    double raan;   // right ascension of ascending node (rad)
    double argp;   // argument of perigee (rad)
    double nu;     // true anomaly (rad)
};

// ── Spacecraft Properties ──

struct SpacecraftProperties {
    double mass_kg    = 260.0;   // wet mass (kg) — Starlink default
    double drag_area  = 22.0;    // drag cross-section (m²)
    double cd         = 2.2;     // drag coefficient
    double srp_area   = 22.0;    // SRP cross-section (m²)
    double cr         = 1.3;     // reflectivity coefficient (1.0 = absorb, 2.0 = reflect)
};

// ── Conversions ──

inline KeplerianElements cartesian_to_keplerian(const State& s, double mu = MU_EARTH) {
    KeplerianElements kep{};

    double r = s.rmag();
    double v = s.vmag();

    // Angular momentum
    double h[3] = {
        s.r[1]*s.v[2] - s.r[2]*s.v[1],
        s.r[2]*s.v[0] - s.r[0]*s.v[2],
        s.r[0]*s.v[1] - s.r[1]*s.v[0]
    };
    double hmag = std::sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);

    // Node vector
    double n[3] = {-h[1], h[0], 0.0};
    double nmag = std::sqrt(n[0]*n[0] + n[1]*n[1]);

    // Eccentricity vector
    double rdotv = s.r[0]*s.v[0] + s.r[1]*s.v[1] + s.r[2]*s.v[2];
    double evec[3];
    for (int j = 0; j < 3; j++)
        evec[j] = ((v*v - mu/r) * s.r[j] - rdotv * s.v[j]) / mu;
    kep.e = std::sqrt(evec[0]*evec[0] + evec[1]*evec[1] + evec[2]*evec[2]);

    // Semi-major axis
    double energy = v*v / 2.0 - mu / r;
    kep.a = -mu / (2.0 * energy);

    // Inclination
    kep.i = std::acos(std::clamp(h[2] / hmag, -1.0, 1.0));

    // RAAN
    if (nmag > 1e-12) {
        kep.raan = std::acos(std::clamp(n[0] / nmag, -1.0, 1.0));
        if (n[1] < 0) kep.raan = TWO_PI - kep.raan;
    }

    // Argument of perigee
    if (nmag > 1e-12 && kep.e > 1e-12) {
        double ndote = (n[0]*evec[0] + n[1]*evec[1] + n[2]*evec[2]) / (nmag * kep.e);
        kep.argp = std::acos(std::clamp(ndote, -1.0, 1.0));
        if (evec[2] < 0) kep.argp = TWO_PI - kep.argp;
    }

    // True anomaly
    if (kep.e > 1e-12) {
        double edotr = (evec[0]*s.r[0] + evec[1]*s.r[1] + evec[2]*s.r[2]) / (kep.e * r);
        kep.nu = std::acos(std::clamp(edotr, -1.0, 1.0));
        if (rdotv < 0) kep.nu = TWO_PI - kep.nu;
    }

    return kep;
}

inline State keplerian_to_cartesian(const KeplerianElements& kep, double mu = MU_EARTH) {
    State s{};

    double p = kep.a * (1.0 - kep.e * kep.e);
    double r = p / (1.0 + kep.e * std::cos(kep.nu));

    // Position in perifocal frame
    double rp[3] = {r * std::cos(kep.nu), r * std::sin(kep.nu), 0.0};

    // Velocity in perifocal frame
    double sqrt_mu_p = std::sqrt(mu / p);
    double vp[3] = {
        -sqrt_mu_p * std::sin(kep.nu),
         sqrt_mu_p * (kep.e + std::cos(kep.nu)),
         0.0
    };

    // Rotation matrix: perifocal → ECI
    double cO = std::cos(kep.raan), sO = std::sin(kep.raan);
    double cw = std::cos(kep.argp), sw = std::sin(kep.argp);
    double ci = std::cos(kep.i),    si = std::sin(kep.i);

    double R[3][3] = {
        {cO*cw - sO*sw*ci, -cO*sw - sO*cw*ci,  sO*si},
        {sO*cw + cO*sw*ci, -sO*sw + cO*cw*ci, -cO*si},
        {sw*si,             cw*si,              ci    }
    };

    for (int j = 0; j < 3; j++) {
        s.r[j] = R[j][0]*rp[0] + R[j][1]*rp[1] + R[j][2]*rp[2];
        s.v[j] = R[j][0]*vp[0] + R[j][1]*vp[1] + R[j][2]*vp[2];
    }

    return s;
}

}  // namespace hpop
