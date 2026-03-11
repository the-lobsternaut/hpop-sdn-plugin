#pragma once
/**
 * HPOP Atmosphere Models
 *
 * Provides atmospheric density for drag computation.
 *
 * Models (in order of increasing fidelity):
 *   1. ExponentialAtmosphere — ρ = ρ₀ exp(-h/H), constant scale height
 *   2. HarrisPriester — altitude-dependent density with diurnal bulge,
 *      parameterized by F10.7 solar flux (Montenbruck & Gill, 2000)
 *   3. [Future] Jacchia 1971 — the HASDM baseline model
 *   4. [Future] JB2008 — Jacchia-Bowman 2008 (modern HASDM)
 *
 * HASDM (High Accuracy Satellite Drag Model) uses a modified Jacchia model
 * as its baseline, with DCA (Dynamic Calibration of the Atmosphere) applying
 * real-time corrections from ~75 calibration satellites.
 *
 * References:
 *   - Montenbruck & Gill, "Satellite Orbits" (2000), Section 3.5
 *   - Harris & Priester, "Time-dependent structure of the upper atmosphere" (1962)
 *   - Jacchia, "New static models of the thermosphere" (1971), SAO Special Report 332
 *   - Bowman et al., "A New Empirical Thermospheric Density Model JB2008" (2008), AIAA 2008-6438
 *   - Storz et al., "High Accuracy Satellite Drag Model (HASDM)" (2005), AIAA 2005-6386
 */

#include "hpop/state.h"
#include "hpop/force_model.h"
#include <cmath>
#include <vector>

namespace hpop {

// ── Abstract Atmosphere Interface ──

class AtmosphereModel {
public:
    virtual ~AtmosphereModel() = default;

    /**
     * Get atmospheric density at a given position.
     *
     * @param alt_km      Geodetic altitude (km)
     * @param lat_rad     Geocentric latitude (rad)
     * @param lon_rad     Geocentric longitude (rad)
     * @param epoch_jd    Julian date of epoch
     * @return density in kg/m³
     */
    virtual double density(double alt_km, double lat_rad, double lon_rad,
                           double epoch_jd) const = 0;
};

// ── Exponential Atmosphere ──
// ρ(h) = ρ₀ exp(-(h - h₀) / H)
// Simple model for testing and low-fidelity propagation.

class ExponentialAtmosphere : public AtmosphereModel {
public:
    // Default: Earth sea-level values
    // ρ₀ = 1.225 kg/m³, H = 8.5 km (tropospheric scale height)
    // For LEO drag, use ρ₀ = 3.614e-13 kg/m³ at h₀ = 700 km, H = 88.667 km
    ExponentialAtmosphere(double rho0 = 1.225, double h0 = 0.0, double H = 8.5)
        : rho0_(rho0), h0_(h0), H_(H) {}

    double density(double alt_km, double lat_rad, double lon_rad,
                   double epoch_jd) const override {
        return rho0_ * std::exp(-(alt_km - h0_) / H_);
    }

private:
    double rho0_;  // Reference density (kg/m³)
    double h0_;    // Reference altitude (km)
    double H_;     // Scale height (km)
};

// ── Harris-Priester Atmosphere ──
//
// Altitude-dependent density with diurnal variation (day/night asymmetry).
// Parameterized by F10.7 solar radio flux for solar cycle dependence.
//
// The model interpolates between minimum (nighttime) and maximum (daytime)
// density profiles, with the bulge centered on the apex of the diurnal
// variation (displaced from subsolar point by ~30° in RA).
//
// Reference: Montenbruck & Gill (2000), Table 3.3 and Algorithm 3.4
//
// Valid range: 100-1000 km altitude

class HarrisPriester : public AtmosphereModel {
public:
    /**
     * @param f107    F10.7 solar radio flux (SFU, typically 70-250)
     * @param n_pratt Pratt exponent for diurnal density variation (2-6, default 2)
     *                Higher values = sharper day/night contrast
     */
    HarrisPriester(double f107 = 150.0, int n_pratt = 2);

    double density(double alt_km, double lat_rad, double lon_rad,
                   double epoch_jd) const override;

    void set_f107(double f107) { f107_ = f107; update_tables(); }
    void set_n_pratt(int n) { n_pratt_ = n; }

    // Sun position needed for diurnal bulge calculation
    // Approximate sun RA/Dec from epoch JD
    static void sun_position(double jd, double& ra_rad, double& dec_rad);

private:
    double f107_;
    int n_pratt_;

    // Density tables (altitude, rho_min, rho_max) in kg/m³
    struct DensityEntry {
        double alt_km;
        double rho_min;  // nighttime (kg/m³)
        double rho_max;  // daytime (kg/m³)
    };

    std::vector<DensityEntry> table_;

    void update_tables();
};

// ── Drag Force Model ──
//
// a_drag = -½ · (Cd·A/m) · ρ · |v_rel| · v̂_rel
//
// Where v_rel is the satellite velocity relative to the co-rotating atmosphere:
//   v_rel = v_eci - ω_earth × r_eci
//
// This inherits from ForceModel so it can be added to the propagator.

class DragForce : public ForceModel {
public:
    /**
     * @param atm       Atmosphere model (ownership transferred)
     * @param cd        Drag coefficient (typically 2.0-2.5 for LEO)
     */
    DragForce(std::unique_ptr<AtmosphereModel> atm, double cd = 2.2)
        : atm_(std::move(atm)), cd_(cd) {}

    void acceleration(const double r[3], const double v[3],
                      double epoch_jd, const SpacecraftProperties* sc,
                      double a_out[3]) const override;

    const char* name() const override { return "Drag"; }

    void set_cd(double cd) { cd_ = cd; }

    // DCA (Dynamic Calibration of the Atmosphere) multiplier
    // HASDM applies this as a density correction factor
    void set_dca_factor(double f) { dca_factor_ = f; }

private:
    std::unique_ptr<AtmosphereModel> atm_;
    double cd_;
    double dca_factor_ = 1.0;  // HASDM density correction
};

}  // namespace hpop
