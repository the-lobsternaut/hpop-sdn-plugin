/**
 * HPOP Atmosphere Models Implementation
 *
 * Harris-Priester density table from Montenbruck & Gill (2000), Table 3.3
 * These are for moderate solar activity (F10.7 ≈ 150 SFU). The model
 * scales linearly between low (F10.7=70) and high (F10.7=250) activity.
 *
 * Drag force computation follows Vallado (2013), Section 8.6.2
 */

#include "hpop/atmosphere.h"
#include <algorithm>

namespace hpop {

// ── Harris-Priester Density Table ──
// Source: Montenbruck & Gill (2000), Table 3.3
// Format: {altitude(km), rho_min(kg/m³), rho_max(kg/m³)}
// rho_min = nighttime, rho_max = daytime

static const struct { double h; double rho_min; double rho_max; } HP_TABLE[] = {
    { 100, 4.974e-07, 4.974e-07},
    { 120, 2.490e-08, 2.490e-08},
    { 130, 8.377e-09, 8.710e-09},
    { 140, 3.899e-09, 4.059e-09},
    { 150, 2.122e-09, 2.215e-09},
    { 160, 1.263e-09, 1.344e-09},
    { 170, 8.008e-10, 8.758e-10},
    { 180, 5.283e-10, 6.010e-10},
    { 190, 3.617e-10, 4.297e-10},
    { 200, 2.557e-10, 3.162e-10},
    { 210, 1.839e-10, 2.396e-10},
    { 220, 1.341e-10, 1.853e-10},
    { 230, 9.949e-11, 1.455e-10},
    { 240, 7.488e-11, 1.157e-10},
    { 250, 5.709e-11, 9.308e-11},
    { 260, 4.403e-11, 7.555e-11},
    { 270, 3.430e-11, 6.182e-11},
    { 280, 2.697e-11, 5.095e-11},
    { 290, 2.139e-11, 4.226e-11},
    { 300, 1.708e-11, 3.526e-11},
    { 320, 1.099e-11, 2.511e-11},
    { 340, 7.214e-12, 1.819e-11},
    { 360, 4.824e-12, 1.337e-11},
    { 380, 3.274e-12, 9.955e-12},
    { 400, 2.249e-12, 7.492e-12},
    { 420, 1.558e-12, 5.684e-12},
    { 440, 1.091e-12, 4.355e-12},
    { 460, 7.701e-13, 3.362e-12},
    { 480, 5.474e-13, 2.612e-12},
    { 500, 3.916e-13, 2.042e-12},
    { 520, 2.819e-13, 1.605e-12},
    { 540, 2.042e-13, 1.267e-12},
    { 560, 1.488e-13, 1.005e-12},
    { 580, 1.092e-13, 7.997e-13},
    { 600, 8.070e-14, 6.390e-13},
    { 620, 6.012e-14, 5.123e-13},
    { 640, 4.519e-14, 4.121e-13},
    { 660, 3.430e-14, 3.325e-13},
    { 680, 2.632e-14, 2.691e-13},
    { 700, 2.043e-14, 2.185e-13},
    { 720, 1.607e-14, 1.779e-13},
    { 740, 1.281e-14, 1.452e-13},
    { 760, 1.036e-14, 1.190e-13},
    { 780, 8.496e-15, 9.776e-14},
    { 800, 7.069e-15, 8.059e-14},
    { 840, 4.680e-15, 5.741e-14},
    { 880, 3.200e-15, 4.210e-14},
    { 920, 2.210e-15, 3.130e-14},
    { 960, 1.560e-15, 2.360e-14},
    {1000, 1.150e-15, 1.810e-14},
};

static constexpr int HP_TABLE_SIZE = sizeof(HP_TABLE) / sizeof(HP_TABLE[0]);

// ── Harris-Priester Implementation ──

HarrisPriester::HarrisPriester(double f107, int n_pratt)
    : f107_(f107), n_pratt_(n_pratt)
{
    update_tables();
}

void HarrisPriester::update_tables() {
    table_.resize(HP_TABLE_SIZE);

    // Scale factor for solar activity
    // F10.7 ranges from ~70 (solar min) to ~250 (solar max)
    // The M&G table is for moderate activity; we interpolate
    // between min and max based on F10.7
    double f_scale = std::max(0.0, std::min(1.0, (f107_ - 70.0) / (250.0 - 70.0)));

    for (int i = 0; i < HP_TABLE_SIZE; i++) {
        table_[i].alt_km  = HP_TABLE[i].h;
        // Scale density based on solar activity
        // At solar min, use rho_min for both day and night
        // At solar max, use full rho_max for daytime
        table_[i].rho_min = HP_TABLE[i].rho_min * (1.0 + 0.5 * f_scale);
        table_[i].rho_max = HP_TABLE[i].rho_max * (1.0 + 0.5 * f_scale);
    }
}

void HarrisPriester::sun_position(double jd, double& ra_rad, double& dec_rad) {
    // Low-precision Sun position (Montenbruck & Gill, Section 3.3.2)
    // Accurate to ~1° — sufficient for diurnal bulge calculation
    double T = (jd - 2451545.0) / 36525.0;  // Julian centuries from J2000

    // Mean anomaly (deg)
    double M = 357.5256 + 35999.049 * T;
    M = std::fmod(M, 360.0);
    if (M < 0) M += 360.0;
    double M_rad = M * DEG2RAD;

    // Ecliptic longitude (deg)
    double lambda = 280.460 + 36000.770 * T + 1.915 * std::sin(M_rad) + 0.020 * std::sin(2.0 * M_rad);
    lambda = std::fmod(lambda, 360.0);
    if (lambda < 0) lambda += 360.0;
    double lambda_rad = lambda * DEG2RAD;

    // Obliquity of ecliptic (deg)
    double epsilon = 23.4393 - 0.0130 * T;
    double eps_rad = epsilon * DEG2RAD;

    // Right ascension and declination
    ra_rad = std::atan2(std::cos(eps_rad) * std::sin(lambda_rad), std::cos(lambda_rad));
    if (ra_rad < 0) ra_rad += 2.0 * M_PI;
    dec_rad = std::asin(std::sin(eps_rad) * std::sin(lambda_rad));
}

double HarrisPriester::density(double alt_km, double lat_rad, double lon_rad,
                                double epoch_jd) const {
    // Clamp altitude
    if (alt_km < table_.front().alt_km) return table_.front().rho_max;
    if (alt_km > table_.back().alt_km) {
        // Extrapolate exponentially above 1000 km
        double h_scale = (table_.back().alt_km - table_[HP_TABLE_SIZE-2].alt_km);
        double rho_ratio = table_.back().rho_min / table_[HP_TABLE_SIZE-2].rho_min;
        double H = -h_scale / std::log(rho_ratio);
        return table_.back().rho_min * std::exp(-(alt_km - table_.back().alt_km) / H);
    }

    // Find bracketing entries
    int i = 0;
    while (i < HP_TABLE_SIZE - 2 && table_[i+1].alt_km < alt_km) i++;

    // Log-linear interpolation in altitude
    double h0 = table_[i].alt_km;
    double h1 = table_[i+1].alt_km;
    double dh = h1 - h0;

    // Scale heights for min and max
    double H_min = dh / std::log(table_[i].rho_min / table_[i+1].rho_min);
    double H_max = dh / std::log(table_[i].rho_max / table_[i+1].rho_max);

    double rho_min = table_[i].rho_min * std::exp(-(alt_km - h0) / H_min);
    double rho_max = table_[i].rho_max * std::exp(-(alt_km - h0) / H_max);

    // Diurnal variation using Pratt's cosine power law
    // Compute the cosine of the angle between the satellite position
    // and the diurnal density bulge apex
    double sun_ra, sun_dec;
    sun_position(epoch_jd, sun_ra, sun_dec);

    // Bulge apex is displaced ~30° in RA from the subsolar point (lag angle)
    double bulge_ra = sun_ra + 30.0 * DEG2RAD;
    double bulge_dec = sun_dec;

    // Satellite direction (from lat/lon)
    double sat_x = std::cos(lat_rad) * std::cos(lon_rad);
    double sat_y = std::cos(lat_rad) * std::sin(lon_rad);
    double sat_z = std::sin(lat_rad);

    // Bulge apex direction
    double bulge_x = std::cos(bulge_dec) * std::cos(bulge_ra);
    double bulge_y = std::cos(bulge_dec) * std::sin(bulge_ra);
    double bulge_z = std::sin(bulge_dec);

    // Cosine of angle between satellite and bulge
    double cos_psi = sat_x * bulge_x + sat_y * bulge_y + sat_z * bulge_z;

    // Pratt's formula: diurnal variation factor
    // cos_psi_half² = (1 + cos_psi) / 2  (half-angle identity)
    double cos_psi_half_sq = 0.5 * (1.0 + cos_psi);
    double diurnal = std::pow(std::max(0.0, cos_psi_half_sq), n_pratt_ / 2.0);

    // Interpolate between nighttime (min) and daytime (max)
    return rho_min + (rho_max - rho_min) * diurnal;
}

// ── Drag Force Implementation ──

void DragForce::acceleration(
    const double r[3], const double v[3],
    double epoch_jd, const SpacecraftProperties* sc,
    double a_out[3]) const
{
    // Get spacecraft properties
    double area = sc ? sc->drag_area : 22.0;    // default 22 m² (Starlink)
    double mass = sc ? sc->mass_kg : 260.0;     // default 260 kg (Starlink)
    double cd = sc ? sc->cd : cd_;

    // Compute geodetic altitude, latitude, longitude
    double rmag = std::sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);
    double alt_km = rmag - RE_KM;  // Spherical approximation
    double lat_rad = std::asin(r[2] / rmag);
    double lon_rad = std::atan2(r[1], r[0]);

    // Relative velocity: satellite velocity minus atmosphere co-rotation
    // v_atm = ω_earth × r_eci
    // ω_earth = [0, 0, 7.2921159e-5] rad/s = [0, 0, 7.2921159e-5 * 1e-3] rad/ms... no
    // ω_earth = 7.2921159e-5 rad/s
    // In km/s: v_atm = ω × r where ω is in rad/s and r in km → v in km/s
    static constexpr double OMEGA_EARTH = 7.2921159e-5;  // rad/s

    double v_rel[3] = {
        v[0] + OMEGA_EARTH * r[1],   // v_x - (-ω*y) = v_x + ω*y
        v[1] - OMEGA_EARTH * r[0],   // v_y - (ω*x) = v_y - ω*x
        v[2]
    };

    double vmag = std::sqrt(v_rel[0]*v_rel[0] + v_rel[1]*v_rel[1] + v_rel[2]*v_rel[2]);

    // Get atmospheric density
    double rho = atm_->density(alt_km, lat_rad, lon_rad, epoch_jd);

    // Apply DCA correction factor (HASDM)
    rho *= dca_factor_;

    // Drag acceleration: a = -½ · (Cd·A/m) · ρ · |v_rel| · v_rel
    // ρ is in kg/m³, area in m², mass in kg
    // v_rel is in km/s, so multiply by 1000 to get m/s for the magnitude squared
    // Result needs to be in km/s²
    //
    // a (m/s²) = -½ · (Cd·A/m) · ρ · |v_rel_m/s|² · v̂_rel
    //          = -½ · (Cd·A/m) · ρ · (vmag*1000)² · v̂_rel
    // a (km/s²) = a (m/s²) / 1000
    //           = -½ · (Cd·A/m) · ρ · vmag² · 1e6 / 1e3 · v̂_rel
    //           = -½ · (Cd·A/m) · ρ · vmag² · 1e3 · v̂_rel

    double factor = -0.5 * (cd * area / mass) * rho * vmag * 1e3;  // 1e3 from unit conversion

    a_out[0] = factor * v_rel[0];
    a_out[1] = factor * v_rel[1];
    a_out[2] = factor * v_rel[2];
}

}  // namespace hpop
