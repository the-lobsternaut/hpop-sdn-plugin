/**
 * OCM FlatBuffer → Raw Text VCM Writer — Implementation
 *
 * Reconstructs the 18 SPCS fixed-format ASCII VCM from OCM FlatBuffer data.
 * Maps OCM fields back to VCM line format.
 */

#include "hpop/vcm_writer.h"
#include "sds/ocm_generated.h"

#include <cstdio>
#include <cstring>
#include <sstream>
#include <cmath>

namespace hpop {

// ── Helpers ──

// ISO 8601 → VCM date format: "yyyy ddd (dd mmm) hh:mm:ss.sss"
static std::string iso_to_vcm_date(const char* iso) {
    if (!iso || std::strlen(iso) < 10) return "";

    int year, month, day, hour = 0, min = 0;
    double sec = 0.0;
    std::sscanf(iso, "%d-%d-%dT%d:%d:%lf", &year, &month, &day, &hour, &min, &sec);

    // Month names
    static const char* months[] = {
        "JAN","FEB","MAR","APR","MAY","JUN",
        "JUL","AUG","SEP","OCT","NOV","DEC"
    };

    // Day of year
    int month_days[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    bool leap = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    if (leap) month_days[1] = 29;

    int doy = day;
    for (int m = 0; m < month - 1; m++) doy += month_days[m];

    int isec = static_cast<int>(sec);
    int msec = static_cast<int>((sec - isec) * 1000 + 0.5);
    if (msec >= 1000) { msec -= 1000; isec++; }

    char buf[64];
    std::snprintf(buf, sizeof(buf), "%04d %03d (%02d %s) %02d:%02d:%02d.%03d",
                  year, doy, day, months[month - 1], hour, min, isec, msec);
    return std::string(buf);
}

// ATM enum → drag model string
static std::string atm_to_drag_string(AtmosphericModelFamily model, int year) {
    switch (model) {
        case AtmosphericModelFamily_JXX:
            if (year == 1970) return "JACCHIA 70";
            if (year == 1971) return "JACCHIA 71";
            return "JACCHIA " + std::to_string(year);
        case AtmosphericModelFamily_JR71: return "JR71";
        case AtmosphericModelFamily_JB08: return "JB2008";
        case AtmosphericModelFamily_JAC_HASDM: return "JAC70/HASDM";
        case AtmosphericModelFamily_NRLMSIS00E: return "NRLMSISE-00";
        case AtmosphericModelFamily_MSISE_90: return "MSIS90";
        case AtmosphericModelFamily_MSIS_86: return "MSIS86";
        case AtmosphericModelFamily_DTM_XX: return "DTM-" + std::to_string(year);
        case AtmosphericModelFamily_HP: return "HARRIS-PRIESTER";
        default: return "UNKNOWN";
    }
}

// Find user-defined parameter value
static std::string find_param(const OCM* ocm, const char* name) {
    auto* params = ocm->USER_DEFINED_PARAMETERS();
    if (!params) return "";
    for (unsigned i = 0; i < params->size(); i++) {
        auto* p = params->Get(i);
        if (p->PARAM_NAME() && std::strcmp(p->PARAM_NAME()->c_str(), name) == 0) {
            return p->PARAM_VALUE() ? p->PARAM_VALUE()->str() : "";
        }
    }
    return "";
}

static double find_param_double(const OCM* ocm, const char* name) {
    auto s = find_param(ocm, name);
    if (s.empty()) return 0.0;
    return std::strtod(s.c_str(), nullptr);
}

// ── Main Writer ──

std::string ocm_to_text_vcm(const uint8_t* data, size_t size) {
    if (size < 8) return "";

    auto* ocm = flatbuffers::GetRoot<OCM>(data);
    if (!ocm) return "";

    auto* hdr = ocm->HEADER();
    auto* meta = ocm->METADATA();
    auto* state = ocm->STATE_DATA();
    auto* phys = ocm->PHYSICAL_PROPERTIES();
    auto* pert = ocm->PERTURBATIONS();
    auto* cov = ocm->COVARIANCE_DATA();

    std::ostringstream out;

    // Line 1: Header
    out << "<> SP VECTOR/COVARIANCE MESSAGE - V2.0\n";

    // Line 2: Real/Test indicator
    out << "<>\n";

    // Line 3: Message time + center
    std::string msg_time;
    if (hdr && hdr->CREATION_DATE())
        msg_time = iso_to_vcm_date(hdr->CREATION_DATE()->c_str());
    std::string center = (hdr && hdr->ORIGINATOR()) ? hdr->ORIGINATOR()->str() : "";
    char line[256];
    std::snprintf(line, sizeof(line),
        "<> MESSAGE TIME (UTC): %-35s CENTER: %s",
        msg_time.c_str(), center.c_str());
    out << line << "\n";

    // Line 4: Satellite number + Int Des
    std::string sat_num = find_param(ocm, "NORAD_CAT_ID");
    std::string obj_des = (meta && meta->OBJECT_DESIGNATOR()) ? meta->OBJECT_DESIGNATOR()->str() : "00000";
    std::string intl_des = (meta && meta->INTERNATIONAL_DESIGNATOR()) ? meta->INTERNATIONAL_DESIGNATOR()->str() : "";
    if (sat_num.empty()) sat_num = obj_des;
    std::snprintf(line, sizeof(line),
        "<> SATELLITE NUMBER: %-25s INT. DES.: %s",
        sat_num.c_str(), intl_des.c_str());
    out << line << "\n";

    // Line 5: Common name
    std::string name = (meta && meta->OBJECT_NAME()) ? meta->OBJECT_NAME()->str() : "";
    out << "<> COMMON NAME: " << name << "\n";

    // Line 6: Epoch time + rev
    std::string epoch;
    if (meta && meta->EPOCH_TZERO())
        epoch = iso_to_vcm_date(meta->EPOCH_TZERO()->c_str());
    std::string rev = find_param(ocm, "EPOCH_REV");
    std::snprintf(line, sizeof(line),
        "<> EPOCH TIME (UTC): %-35s EPOCH REV: %s",
        epoch.c_str(), rev.c_str());
    out << line << "\n";

    // Lines 7-8: J2K POS/VEL
    if (state && state->size() >= 6) {
        std::snprintf(line, sizeof(line),
            "<> J2K POS (KM):   %16.8f %16.8f %16.8f",
            state->Get(0), state->Get(1), state->Get(2));
        out << line << "\n";
        std::snprintf(line, sizeof(line),
            "<> J2K VEL (KM/S): %16.12f %16.12f %16.12f",
            state->Get(3), state->Get(4), state->Get(5));
        out << line << "\n";
    }

    // Lines 9-10: ECI POS/VEL (from user params)
    std::string eci_pos = find_param(ocm, "ECI_POS_KM");
    std::string eci_vel = find_param(ocm, "ECI_VEL_KMS");
    if (!eci_pos.empty()) {
        double ex, ey, ez;
        std::sscanf(eci_pos.c_str(), "%lf %lf %lf", &ex, &ey, &ez);
        std::snprintf(line, sizeof(line),
            "<> ECI POS (KM):   %16.8f %16.8f %16.8f", ex, ey, ez);
        out << line << "\n";
    }
    if (!eci_vel.empty()) {
        double evx, evy, evz;
        std::sscanf(eci_vel.c_str(), "%lf %lf %lf", &evx, &evy, &evz);
        std::snprintf(line, sizeof(line),
            "<> ECI VEL (KM/S): %16.12f %16.12f %16.12f", evx, evy, evz);
        out << line << "\n";
    }

    // Lines 11-12: EFG POS/VEL (from user params)
    std::string efg_pos = find_param(ocm, "EFG_POS_KM");
    std::string efg_vel = find_param(ocm, "EFG_VEL_KMS");
    if (!efg_pos.empty()) {
        double fx, fy, fz;
        std::sscanf(efg_pos.c_str(), "%lf %lf %lf", &fx, &fy, &fz);
        std::snprintf(line, sizeof(line),
            "<> EFG POS (KM):   %16.8f %16.8f %16.8f", fx, fy, fz);
        out << line << "\n";
    }
    if (!efg_vel.empty()) {
        double fvx, fvy, fvz;
        std::sscanf(efg_vel.c_str(), "%lf %lf %lf", &fvx, &fvy, &fvz);
        std::snprintf(line, sizeof(line),
            "<> EFG VEL (KM/S): %16.12f %16.12f %16.12f", fvx, fvy, fvz);
        out << line << "\n";
    }

    // Line 13: Geopotential, drag, lunar/solar
    std::string geo_model = (pert && pert->GRAVITY_MODEL()) ? pert->GRAVITY_MODEL()->str() : "";
    int geo_deg = pert ? pert->GRAVITY_DEGREE() : 0;
    int geo_ord = pert ? pert->GRAVITY_ORDER() : 0;

    std::string drag = "UNKNOWN";
    if (pert && pert->ATMOSPHERIC_MODEL()) {
        drag = atm_to_drag_string(
            pert->ATMOSPHERIC_MODEL()->MODEL(),
            pert->ATMOSPHERIC_MODEL()->YEAR());
    }

    bool ls_on = false;
    if (pert && pert->N_BODY_PERTURBATIONS() && pert->N_BODY_PERTURBATIONS()->size() > 0)
        ls_on = true;

    std::snprintf(line, sizeof(line),
        "<> GEOPOTENTIAL: %s %dZ,%dT  DRAG: %-12s LUNAR/SOLAR: %s",
        geo_model.c_str(), geo_deg, geo_ord, drag.c_str(), ls_on ? " ON" : "OFF");
    out << line << "\n";

    // Line 14: SRP, solid tides, thrust
    std::string srp_status = (pert && pert->SOLAR_RAD_PRESSURE()) ? pert->SOLAR_RAD_PRESSURE()->str() : "OFF";
    std::string tides = find_param(ocm, "SOLID_EARTH_TIDES");
    if (tides.empty()) {
        // Check if solid tides model is set
        tides = (pert && pert->SOLID_TIDES_MODEL()) ? "ON" : "OFF";
    }
    std::snprintf(line, sizeof(line),
        "<> SOLAR RAD PRESS: %s  SOLID EARTH TIDES: %s  IN-TRACK THRUST: OFF",
        srp_status.c_str(), tides.c_str());
    out << line << "\n";

    // Line 15: Ballistic coeff, BDOT
    double bc = phys ? phys->DRAG_COEFF_NOM() : 0;
    double bdot = find_param_double(ocm, "BDOT_M2_KG_S");
    std::snprintf(line, sizeof(line),
        "<> BALLISTIC COEF (M2/KG): %13.6E BDOT (M2/KG-S): %12.6E", bc, bdot);
    out << line << "\n";

    // Line 16: SRP coeff, EDR
    double srpc = phys ? phys->SOLAR_RAD_COEFF() : 0;
    double edr = find_param_double(ocm, "EDR_W_KG");
    std::snprintf(line, sizeof(line),
        "<> SOLAR RAD PRESS COEFF (M2/KG): %13.6E  EDR(W/KG): %9.2E", srpc, edr);
    out << line << "\n";

    // Line 17: Thrust accel, CM offset
    double ta = find_param_double(ocm, "THRUST_ACCEL_MS2");
    double cm = find_param_double(ocm, "CM_OFFSET_M");
    std::snprintf(line, sizeof(line),
        "<> THRUST ACCEL (M/S2): %13.6E  C.M. OFFSET (M): %13.6E", ta, cm);
    out << line << "\n";

    // Line 18: Solar flux
    double f10 = pert ? pert->FIXED_F10P7() : 0;
    double avg_f10 = pert ? pert->FIXED_F10P7_MEAN() : 0;
    double avg_ap = pert ? pert->FIXED_GEOMAG_KP() : 0;
    std::snprintf(line, sizeof(line),
        "<> SOLAR FLUX: F10: %3.0f  AVERAGE F10: %3.0f  AVERAGE AP: %5.1f",
        f10, avg_f10, avg_ap);
    out << line << "\n";

    // Line 19: TAI-UTC, UT1-UTC, UT1 rate
    double tai = meta ? meta->TAIMUTC_AT_TZERO() : 0;
    double ut1 = meta ? meta->UT1MUTC_AT_TZERO() : 0;
    double ut1_rate = find_param_double(ocm, "UT1_RATE_MS_DAY");
    std::snprintf(line, sizeof(line),
        "<> TAI-UTC (S): %2.0f  UT1-UTC (S): %8.5f  UT1 RATE (MS/DAY): %6.3f",
        tai, ut1, ut1_rate);
    out << line << "\n";

    // Line 20: Polar motion, nutation
    std::string polar = find_param(ocm, "POLAR_MOTION_XY_ARCSEC");
    double px = 0, py = 0;
    if (!polar.empty()) std::sscanf(polar.c_str(), "%lf %lf", &px, &py);
    std::snprintf(line, sizeof(line),
        "<> POLAR MOT X,Y (ARCSEC): %7.4f %7.4f IAU 1980 NUTAT:   4 TERMS",
        px, py);
    out << line << "\n";

    // Line 21: Leap second time
    std::string ls_time = "2049 365 (31 DEC) 23:59:59.999";
    if (meta && meta->NEXT_LEAP_EPOCH()) {
        std::string t = iso_to_vcm_date(meta->NEXT_LEAP_EPOCH()->c_str());
        if (!t.empty()) ls_time = t;
    }
    out << "<> TIME CONST LEAP SECOND TIME (UTC): " << ls_time << "\n";

    // Line 22: Integrator mode
    std::string int_mode = find_param(ocm, "INTEGRATOR_MODE");
    std::string coord = find_param(ocm, "COORD_SYS");
    std::string partials = find_param(ocm, "PARTIALS");
    std::snprintf(line, sizeof(line),
        "<> INTEGRATOR MODE: %-12s COORD SYS: %-6s PARTIALS: %s",
        int_mode.c_str(), coord.c_str(), partials.c_str());
    out << line << "\n";

    // Line 23: Step mode
    std::string step_mode = find_param(ocm, "STEP_MODE");
    std::snprintf(line, sizeof(line),
        "<> STEP MODE: %-5s FIXED STEP: OFF  STEP SIZE SELECTION: AUTO",
        step_mode.c_str());
    out << line << "\n";

    // Line 24: Step size, error control
    double step = find_param_double(ocm, "INIT_STEP_SIZE_S");
    double err = find_param_double(ocm, "ERROR_CONTROL");
    std::snprintf(line, sizeof(line),
        "<> INITIAL STEP SIZE (S): %8.3f  ERROR CONTROL: %9.3E", step, err);
    out << line << "\n";

    // Line 25: UVW position sigmas
    std::string uvw_pos = find_param(ocm, "UVW_POS_SIGMA_KM");
    double u1 = 0, v1 = 0, w1 = 0;
    if (!uvw_pos.empty()) std::sscanf(uvw_pos.c_str(), "%lf %lf %lf", &u1, &v1, &w1);
    std::snprintf(line, sizeof(line),
        "<> VECTOR U,V,W SIGMAS (KM):        %10.4f %10.4f %10.4f", u1, v1, w1);
    out << line << "\n";

    // Line 26: UVW velocity sigmas
    std::string uvw_vel = find_param(ocm, "UVW_VEL_SIGMA_KMS");
    double u2 = 0, v2 = 0, w2 = 0;
    if (!uvw_vel.empty()) std::sscanf(uvw_vel.c_str(), "%lf %lf %lf", &u2, &v2, &w2);
    std::snprintf(line, sizeof(line),
        "<> VECTOR UD,VD,WD SIGMAS (KM/S):   %10.4f %10.4f %10.4f", u2, v2, w2);
    out << line << "\n";

    // Lines 27+: Covariance matrix
    if (cov && cov->size() > 0) {
        // Determine matrix size from element count: n(n+1)/2 = size
        int n = 0;
        for (int k = 1; k * (k + 1) / 2 <= static_cast<int>(cov->size()); k++) {
            if (k * (k + 1) / 2 == static_cast<int>(cov->size())) n = k;
        }

        std::string rms = find_param(ocm, "WTD_RMS");
        std::snprintf(line, sizeof(line),
            "<> COVARIANCE MATRIX (EQUINOCTIAL ELS): (%2dx%2d) WTD RMS: %13s",
            n, n, rms.c_str());
        out << line << "\n";

        // Write 5 values per line
        for (size_t i = 0; i < cov->size(); i++) {
            if (i % 5 == 0) out << "<> ";
            std::snprintf(line, sizeof(line), " %13.5E", cov->Get(i));
            out << line;
            if ((i + 1) % 5 == 0 || i + 1 == cov->size()) out << "\n";
        }
    } else {
        out << "<> COVARIANCE MATRIX (EQUINOCTIAL ELS): ( 0x 0) WTD RMS:  0.00000E+00\n";
    }

    return out.str();
}

}  // namespace hpop
