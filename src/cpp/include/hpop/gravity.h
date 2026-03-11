#pragma once
/**
 * HPOP Spherical Harmonics Geopotential Model
 *
 * Implements gravitational acceleration from a spherical harmonics expansion.
 * Uses fully-normalized coefficients (geodesy normalization).
 *
 * Algorithm: Montenbruck & Gill, "Satellite Orbits" Ch. 3.2
 *   - Legendre functions via stable column recursion (Holmes & Featherstone 2002)
 *   - Acceleration via gradient of the potential in body-fixed frame, rotated to ECI
 *
 * Reference: Holmes, S.A. and Featherstone, W.E. (2002) "A unified approach to the
 *   Clenshaw summation and the recursive computation of very high degree and order
 *   normalised associated Legendre functions"
 */

#include "force_model.h"
#include <vector>
#include <string>

namespace hpop {

class GravityField : public ForceModel {
public:
    /// Load coefficients from ICGEM GFC format file
    /// max_degree/max_order: truncate to this degree/order (0 = use file max)
    bool load_gfc(const std::string& filepath, int max_degree = 0, int max_order = 0);

    /// Set coefficients directly (for testing)
    void set_coefficients(int max_degree, int max_order,
                          double gm, double radius,
                          const std::vector<std::vector<double>>& C,
                          const std::vector<std::vector<double>>& S);

    /// ForceModel interface
    void acceleration(
        const double r[3], const double v[3],
        double epoch_jd, const SpacecraftProperties* sc,
        double a_out[3]) const override;

    const char* name() const override { return "Geopotential"; }

    int degree() const { return max_degree_; }
    int order() const { return max_order_; }
    double gm() const { return gm_; }
    double radius() const { return radius_; }

private:
    int max_degree_ = 0;
    int max_order_ = 0;
    double gm_ = MU_EARTH * 1e9;   // m³/s² (ICGEM uses SI)
    double radius_ = RE_KM * 1e3;  // m
    // Note: internally we work in km for compatibility with the rest of HPOP

    // Fully-normalized coefficients C[n][m] and S[n][m]
    // Stored as (max_degree+1) x (max_degree+1)
    std::vector<std::vector<double>> C_;
    std::vector<std::vector<double>> S_;

    /// Compute acceleration in ECEF frame, then rotate to ECI
    /// For now, simplified: assume ECI ≈ ECEF (no Earth rotation)
    /// TODO: Add proper GMST rotation for epoch_jd
    void accel_ecef(
        const double r_ecef[3],
        double a_ecef[3]) const;
};

}  // namespace hpop
