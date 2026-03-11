/**
 * HPOP Spherical Harmonics Geopotential
 *
 * Computes gravitational acceleration from normalized spherical harmonic coefficients.
 *
 * The gravitational potential is:
 *   U(r,φ,λ) = (GM/r) Σ_{n=0}^{N} Σ_{m=0}^{n} (R_e/r)^n P_nm(sinφ) [C_nm cos(mλ) + S_nm sin(mλ)]
 *
 * We compute the gradient ∇U in Cartesian body-fixed coordinates using the
 * Cunningham (1970) / Montenbruck & Gill formulation with V/W auxiliary functions.
 *
 * The V/W recursion produces:
 *   V[n][m] = (R_e^n / r^{n+1}) * P_nm(sinφ) * cos(mλ)
 *   W[n][m] = (R_e^n / r^{n+1}) * P_nm(sinφ) * sin(mλ)
 *
 * The acceleration components are (M&G eq 3.33):
 *   a_x = (GM/R_e) Σ_n Σ_m [...]
 *
 * Reference: Montenbruck & Gill (2000), "Satellite Orbits" Ch. 3.2
 */

#include "hpop/gravity.h"
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>

namespace hpop {

bool GravityField::load_gfc(const std::string& filepath, int max_degree, int max_order) {
    std::ifstream file(filepath);
    if (!file.is_open()) return false;

    std::string line;
    int file_max_degree = 0;

    while (std::getline(file, line)) {
        if (line.find("earth_gravity_constant") != std::string::npos) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> gm_;
            gm_ /= 1e9;  // m³/s² → km³/s²
        }
        if (line.substr(0, 6) == "radius") {
            std::istringstream iss(line);
            std::string key;
            double val;
            iss >> key >> val;
            radius_ = val / 1e3;  // m → km
        }
        if (line.find("max_degree") != std::string::npos) {
            std::istringstream iss(line);
            std::string key;
            iss >> key >> file_max_degree;
        }
        if (line.find("end_of_head") != std::string::npos) break;
    }

    max_degree_ = (max_degree > 0) ? std::min(max_degree, file_max_degree) : file_max_degree;
    max_order_ = (max_order > 0) ? std::min(max_order, max_degree_) : max_degree_;

    int N = max_degree_ + 1;
    C_.assign(N, std::vector<double>(N, 0.0));
    S_.assign(N, std::vector<double>(N, 0.0));
    C_[0][0] = 1.0;

    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string key;
        int n, m;
        double c_val, s_val;
        iss >> key >> n >> m >> c_val >> s_val;
        if (key != "gfc") continue;
        if (n > max_degree_ || m > max_order_) continue;
        C_[n][m] = c_val;
        S_[n][m] = s_val;
    }

    return true;
}

void GravityField::set_coefficients(int max_degree, int max_order,
                                     double gm, double radius,
                                     const std::vector<std::vector<double>>& C,
                                     const std::vector<std::vector<double>>& S) {
    max_degree_ = max_degree;
    max_order_ = max_order;
    gm_ = gm;
    radius_ = radius;
    C_ = C;
    S_ = S;
}

// ── Gravitational Acceleration in ECEF ──
//
// Direct computation via potential gradient in spherical coordinates,
// then transform to Cartesian. This is simpler and more debuggable
// than the V/W recursion for a first implementation.
//
// Reference: Vallado, "Fundamentals of Astrodynamics" 5th ed, Section 8.7

void GravityField::accel_ecef(const double r_ecef[3], double a_ecef[3]) const {
    double x = r_ecef[0], y = r_ecef[1], z = r_ecef[2];
    double r2 = x*x + y*y + z*z;
    double r = std::sqrt(r2);
    double r_xy = std::sqrt(x*x + y*y);

    // Geocentric latitude and longitude
    double sinphi = z / r;
    double cosphi = r_xy / r;
    double lambda = std::atan2(y, x);

    int N = max_degree_;

    // Compute FULLY NORMALIZED associated Legendre functions P̄[n][m](sinphi)
    // Using the stable recursion from Holmes & Featherstone (2002)
    // The ICGEM GFC coefficients are already fully normalized.
    //
    // Normalization: P̄_nm = sqrt((2n+1)(2-δ_{0m})(n-m)!/(n+m)!) * P_nm
    //
    // Recursion (Holmes & Featherstone 2002, eq 11-13):
    //   Sectoral:  P̄_mm = sqrt((2m+1)/(2m)) * cosphi * P̄_{m-1,m-1}
    //   Sub-diag:  P̄_{m+1,m} = sqrt(2m+3) * sinphi * P̄_mm
    //   General:   P̄_nm = a_nm * sinphi * P̄_{n-1,m} - b_nm * P̄_{n-2,m}
    //     where a_nm = sqrt((4n²-1)/(n²-m²))
    //           b_nm = sqrt(((2n+1)(n-1-m)(n-1+m))/((2n-3)(n²-m²)))

    int sz = N + 2;
    std::vector<std::vector<double>> P(sz, std::vector<double>(sz, 0.0));

    P[0][0] = 1.0;  // P̄_00 = 1

    for (int m = 0; m < sz; m++) {
        // Sectoral: P̄_mm
        if (m > 0) {
            P[m][m] = std::sqrt((2.0*m + 1.0) / (2.0*m)) * cosphi * P[m-1][m-1];
        }
        // Sub-diagonal: P̄_{m+1,m}
        if (m + 1 < sz) {
            P[m+1][m] = std::sqrt(2.0*m + 3.0) * sinphi * P[m][m];
        }
        // General recursion: P̄_nm for n = m+2 .. sz-1
        for (int n = m + 2; n < sz; n++) {
            double nn = static_cast<double>(n);
            double mm = static_cast<double>(m);
            double a = std::sqrt((4.0*nn*nn - 1.0) / (nn*nn - mm*mm));
            double b = std::sqrt(((2.0*nn + 1.0) * (nn - 1.0 - mm) * (nn - 1.0 + mm))
                                / ((2.0*nn - 3.0) * (nn*nn - mm*mm)));
            P[n][m] = a * sinphi * P[n-1][m] - b * P[n-2][m];
        }
    }

    // Compute dP̄/dφ using the recursion:
    //   dP̄_nm/dφ = sqrt(n²-m²) * sqrt((2n+1)/(2n-1)) * P̄_{n-1,m} - n*sinφ/cosφ * P̄_nm  [m≠n]
    // But simpler: use the identity
    //   dP̄_nm/dφ = n*sinφ*P̄_nm/(cos²φ) - sqrt((n+m)(n-m+1)) * ... [complex]
    // Simplest stable approach for dP̄:
    //   dP̄_nm/dφ = (1/cosφ) * [n*sinφ*P̄_nm - sqrt((2n+1)/(2n-1)) * sqrt((n+m)(n-m)) * P̄_{n-1,m}]  ... still messy
    //
    // Better: use the derivative relation for normalized functions:
    //   cosφ * dP̄_nm/dφ = -m*tanφ*P̄_nm + sqrt((n-m)(n+m+1)) * P̄_{n,m+1}   (for m < n)
    //                    = -n*tanφ*P̄_nn                                       (for m = n)
    // So: dP̄/dφ = [-m*(sinφ/cosφ)*P̄_nm + sqrt((n-m)(n+m+1)) * P̄_{n,m+1}]  / cosφ ... nope, too complex
    //
    // Simplest correct approach: use numerical differentiation or the direct formula:
    //   dP̄_nm/dφ = 0.5 * [sqrt((n-m)(n+m+1)) * P̄_{n,m+1} - sqrt((n+m)(n-m+1)) * P̄_{n,m-1}]  (for m>0)
    //   dP̄_n0/dφ = -sqrt(n(n+1)/2) * P̄_{n,1}
    // This is the standard formula (e.g., Heiskanen & Moritz)

    std::vector<std::vector<double>> dP(sz, std::vector<double>(sz, 0.0));
    for (int n = 1; n < sz; n++) {
        // m = 0
        dP[n][0] = -std::sqrt(n * (n + 1.0) / 2.0) * P[n][1];
        // m > 0
        for (int m = 1; m <= n && m < sz - 1; m++) {
            double nn = n, mm = m;
            double term_plus = std::sqrt((nn - mm) * (nn + mm + 1.0)) * ((m + 1 < sz) ? P[n][m+1] : 0.0);
            double term_minus = std::sqrt((nn + mm) * (nn - mm + 1.0)) * P[n][m-1];
            dP[n][m] = 0.5 * (term_plus - term_minus);
        }
        // m = n: only the minus term
        if (n < sz) {
            double nn = n, mm = n;
            dP[n][n] = -0.5 * std::sqrt((nn + mm) * (nn - mm + 1.0)) * P[n][n-1];
        }
    }

    // Acceleration in spherical coordinates (r, phi, lambda)
    // dU/dr, dU/dphi, dU/dlambda
    double dU_dr = 0.0;
    double dU_dphi = 0.0;
    double dU_dlambda = 0.0;

    double Re_over_r = radius_ / r;
    double Re_over_r_n = 1.0;  // (R_e/r)^n

    for (int n = 0; n <= N; n++) {
        Re_over_r_n *= Re_over_r;  // now (R_e/r)^(n+1)... but we want (R_e/r)^n
        // Actually: after this line, Re_over_r_n = (R_e/r)^{n} for n>=1, need to fix
    }

    // Redo: compute properly
    Re_over_r_n = 1.0;
    for (int n = 0; n <= N; n++) {
        for (int m = 0; m <= std::min(n, max_order_); m++) {
            double Cnm = C_[n][m];
            double Snm = S_[n][m];
            double cos_ml = std::cos(m * lambda);
            double sin_ml = std::sin(m * lambda);
            double CS = Cnm * cos_ml + Snm * sin_ml;
            double CS_d = -m * Cnm * sin_ml + m * Snm * cos_ml;

            dU_dr     += -(n + 1) * Re_over_r_n * P[n][m] * CS;
            dU_dphi   += Re_over_r_n * dP[n][m] * CS;
            dU_dlambda += Re_over_r_n * P[n][m] * CS_d;
        }
        Re_over_r_n *= Re_over_r;
    }

    // Scale by GM/r²
    double GM_r2 = gm_ / r2;
    dU_dr *= GM_r2;
    dU_dphi *= GM_r2;
    dU_dlambda *= GM_r2;

    // Convert spherical acceleration to Cartesian ECEF
    // a_r = dU/dr (radial)
    // a_phi = (1/r) * dU/dphi (latitude)
    // a_lambda = (1/(r*cos(phi))) * dU/dlambda (longitude)

    double a_r = dU_dr;
    double a_phi = dU_dphi / r;
    double a_lambda = (cosphi > 1e-15) ? dU_dlambda / (r * cosphi) : 0.0;

    // Spherical → Cartesian (ECEF)
    double cos_l = std::cos(lambda), sin_l = std::sin(lambda);

    // Unit vectors in spherical coords (at position r,phi,lambda):
    // r̂ = [cosphi*cos_l, cosphi*sin_l, sinphi]
    // φ̂ = [-sinphi*cos_l, -sinphi*sin_l, cosphi]
    // λ̂ = [-sin_l, cos_l, 0]

    a_ecef[0] = a_r * cosphi * cos_l + a_phi * (-sinphi * cos_l) + a_lambda * (-sin_l);
    a_ecef[1] = a_r * cosphi * sin_l + a_phi * (-sinphi * sin_l) + a_lambda * cos_l;
    a_ecef[2] = a_r * sinphi          + a_phi * cosphi;
}

// ── ForceModel Interface ──

void GravityField::acceleration(
    const double r[3], const double v[3],
    double epoch_jd, const SpacecraftProperties* sc,
    double a_out[3]) const
{
    // GMST rotation: ECI → ECEF
    double T = epoch_jd - 2451545.0;
    double gmst_deg = std::fmod(280.46061837 + 360.98564736629 * T, 360.0);
    if (gmst_deg < 0) gmst_deg += 360.0;
    double gmst = gmst_deg * DEG2RAD;

    double cg = std::cos(gmst), sg = std::sin(gmst);

    double r_ecef[3] = {
         cg * r[0] + sg * r[1],
        -sg * r[0] + cg * r[1],
         r[2]
    };

    double a_ecef[3];
    accel_ecef(r_ecef, a_ecef);

    // ECEF → ECI
    a_out[0] = cg * a_ecef[0] - sg * a_ecef[1];
    a_out[1] = sg * a_ecef[0] + cg * a_ecef[1];
    a_out[2] = a_ecef[2];
}

}  // namespace hpop
