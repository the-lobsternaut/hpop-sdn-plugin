/**
 * Raw Text VCM → OCM FlatBuffer Parser — Implementation
 *
 * Parses 18 SPCS fixed-format ASCII VCM (V2.0) line by line.
 * Maps to spacedatastandards.org OCM schema.
 *
 * Line identification uses prefix matching after stripping "<> ".
 * Positional parsing uses fixed column offsets per the VCM format spec.
 */

#include "hpop/vcm_parser.h"
#include "sds/ocm_generated.h"

#include <sstream>
#include <algorithm>
#include <cstring>
#include <cstdlib>
#include <cmath>

namespace hpop {

// ── Helpers ──

static std::string trim(const std::string& s) {
    auto start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    auto end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

static std::string strip_prefix(const std::string& line) {
    // Remove "<> " prefix
    if (line.size() >= 3 && line[0] == '<' && line[1] == '>') {
        size_t i = 2;
        if (i < line.size() && line[i] == ' ') i++;
        return line.substr(i);
    }
    return line;
}

static bool starts_with(const std::string& s, const char* prefix) {
    return s.compare(0, std::strlen(prefix), prefix) == 0;
}

static double parse_double(const std::string& s) {
    std::string trimmed = trim(s);
    if (trimmed.empty()) return 0.0;
    return std::strtod(trimmed.c_str(), nullptr);
}

static int parse_int(const std::string& s) {
    std::string trimmed = trim(s);
    if (trimmed.empty()) return 0;
    return std::atoi(trimmed.c_str());
}

// Parse VCM date format: "yyyy ddd (dd mmm) hh:mm:ss.sss"
// Returns ISO 8601 UTC string
static std::string parse_vcm_date(const std::string& s) {
    // Find year and day-of-year
    int year = 0, doy = 0, hour = 0, min = 0;
    double sec = 0.0;

    // "2023 007 (07 JAN) 21:39: 8.414"
    auto trimmed = trim(s);
    if (trimmed.size() < 25) return "";

    year = parse_int(trimmed.substr(0, 4));
    doy = parse_int(trimmed.substr(5, 3));

    // Find time after the closing paren
    auto paren_close = trimmed.find(')');
    if (paren_close != std::string::npos && paren_close + 2 < trimmed.size()) {
        std::string time_part = trim(trimmed.substr(paren_close + 1));
        // "hh:mm:ss.sss" or "hh:mm: s.sss" (space-padded seconds)
        std::sscanf(time_part.c_str(), "%d:%d:%lf", &hour, &min, &sec);
    }

    // Convert DOY to month/day
    int month_days[] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    bool leap = (year % 4 == 0 && (year % 100 != 0 || year % 400 == 0));
    if (leap) month_days[1] = 29;

    int month = 1, day = doy;
    for (int m = 0; m < 12 && day > month_days[m]; m++) {
        day -= month_days[m];
        month++;
    }

    int isec = static_cast<int>(sec);
    int msec = static_cast<int>((sec - isec) * 1000 + 0.5);

    // Handle carry
    if (msec >= 1000) { msec -= 1000; isec++; }
    if (isec >= 60) { isec -= 60; min++; }
    if (min >= 60) { min -= 60; hour++; }

    char buf[64];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%02d.%03dZ",
                  year, month, day, hour, min, isec, msec);
    return std::string(buf);
}

// Parse scientific notation: " 0.826455E-02" or "-0.826455E-02" or "+0.826455E+02"
static double parse_sci(const std::string& s) {
    return parse_double(s);
}

// Extract value after a label, up to the next label or end
static std::string extract_after(const std::string& line, const char* label) {
    auto pos = line.find(label);
    if (pos == std::string::npos) return "";
    pos += std::strlen(label);
    // Skip colon and spaces
    while (pos < line.size() && (line[pos] == ':' || line[pos] == ' ')) pos++;
    // Find end: next alphabetic label or end of line
    auto end = line.size();
    return trim(line.substr(pos, end - pos));
}

// Extract value between two labels
static std::string extract_between(const std::string& line, const char* after, const char* before) {
    auto pos1 = line.find(after);
    if (pos1 == std::string::npos) return "";
    pos1 += std::strlen(after);
    while (pos1 < line.size() && (line[pos1] == ':' || line[pos1] == ' ')) pos1++;

    auto pos2 = line.find(before, pos1);
    if (pos2 == std::string::npos) return trim(line.substr(pos1));
    return trim(line.substr(pos1, pos2 - pos1));
}

// Map VCM drag model string to ATM enum
static AtmosphericModelFamily map_atm_model(const std::string& drag) {
    if (drag.find("JAC70") != std::string::npos || drag.find("JACCHIA 70") != std::string::npos)
        return AtmosphericModelFamily_JXX;
    if (drag.find("MSIS") != std::string::npos || drag.find("NRLMSISE") != std::string::npos)
        return AtmosphericModelFamily_NRLMSIS00E;
    if (drag.find("JB2008") != std::string::npos || drag.find("JB08") != std::string::npos)
        return AtmosphericModelFamily_JB08;
    if (drag.find("JR71") != std::string::npos || drag.find("JACCHIA-ROBERTS") != std::string::npos)
        return AtmosphericModelFamily_JR71;
    if (drag.find("HASDM") != std::string::npos)
        return AtmosphericModelFamily_JAC_HASDM;
    if (drag.find("DTM") != std::string::npos)
        return AtmosphericModelFamily_DTM_XX;
    if (drag.find("HP") != std::string::npos || drag.find("HARRIS") != std::string::npos)
        return AtmosphericModelFamily_HP;
    return AtmosphericModelFamily_JXX;  // default
}

// Map drag model to year
static int map_atm_year(const std::string& drag) {
    if (drag.find("70") != std::string::npos) return 1970;
    if (drag.find("71") != std::string::npos) return 1971;
    if (drag.find("90") != std::string::npos) return 1990;
    if (drag.find("2000") != std::string::npos || drag.find("00") != std::string::npos) return 2000;
    if (drag.find("2008") != std::string::npos || drag.find("08") != std::string::npos) return 2008;
    return 0;
}

// ── Main Parser ──

VCMParseResult parse_text_vcm(const std::string& vcm_text) {
    VCMParseResult result;

    // Split into lines
    std::vector<std::string> lines;
    std::istringstream stream(vcm_text);
    std::string line;
    while (std::getline(stream, line)) {
        // Remove trailing CR
        if (!line.empty() && line.back() == '\r') line.pop_back();
        lines.push_back(line);
    }

    if (lines.empty()) {
        result.error = "Empty VCM text";
        return result;
    }

    // Parsed fields
    std::string msg_time_utc, center;
    uint32_t sat_num = 0;
    std::string intl_des, common_name, epoch_utc;
    uint32_t epoch_rev = 0;

    // J2K state vector (preferred)
    double j2k_pos[3] = {0}, j2k_vel[3] = {0};
    // ECI state vector
    double eci_pos[3] = {0}, eci_vel[3] = {0};
    // EFG state vector
    double efg_pos[3] = {0}, efg_vel[3] = {0};

    // Perturbation model info
    std::string geopotential_model;
    int geo_degree = 0, geo_order = 0;
    std::string drag_model;
    bool lunar_solar_on = false;
    bool srp_on = false;
    bool solid_tides_on = false;
    bool thrust_on = false;

    // Physical properties
    double ballistic_coeff = 0, bdot = 0;
    double srp_coeff = 0, edr = 0;
    double thrust_accel = 0, cm_offset = 0;

    // Space weather
    double f10 = 0, avg_f10 = 0, avg_ap = 0;

    // Time corrections
    int tai_utc = 0;
    double ut1_utc = 0, ut1_rate = 0;
    double polar_x = 0, polar_y = 0;
    int nutation_terms = 0;
    std::string leap_second_time;

    // Integrator
    std::string integrator_mode, coord_sys, partials;
    std::string step_mode, fixed_step, step_selection;
    double init_step = 0, error_control = 0;

    // Sigmas
    double uvw_pos[3] = {0}, uvw_vel[3] = {0};

    // Covariance
    int cov_rows = 0, cov_cols = 0;
    double wtd_rms = 0;
    std::vector<double> covariance;

    // Parse line by line
    bool in_covariance = false;

    for (const auto& raw_line : lines) {
        std::string l = strip_prefix(raw_line);
        if (l.empty()) continue;

        if (starts_with(l, "SP VECTOR/COVARIANCE")) {
            // Header line
            continue;
        }

        if (starts_with(l, "MESSAGE TIME")) {
            auto time_str = extract_between(l, "MESSAGE TIME (UTC)", "CENTER");
            msg_time_utc = parse_vcm_date(time_str);
            center = extract_after(l, "CENTER");
            continue;
        }

        if (starts_with(l, "SATELLITE NUMBER")) {
            auto sat_str = extract_between(l, "SATELLITE NUMBER", "INT. DES.");
            sat_num = static_cast<uint32_t>(parse_int(sat_str));
            intl_des = extract_after(l, "INT. DES.");
            continue;
        }

        if (starts_with(l, "COMMON NAME")) {
            common_name = trim(l.substr(12));
            if (common_name.size() > 0 && common_name[0] == ':')
                common_name = trim(common_name.substr(1));
            continue;
        }

        if (starts_with(l, "EPOCH TIME")) {
            auto time_str = extract_between(l, "EPOCH TIME (UTC)", "EPOCH REV");
            epoch_utc = parse_vcm_date(time_str);
            auto rev_str = extract_after(l, "EPOCH REV");
            epoch_rev = static_cast<uint32_t>(parse_int(rev_str));
            continue;
        }

        if (starts_with(l, "J2K POS")) {
            // Fixed columns: positions start after "J2K POS (KM):"
            auto data = extract_after(l, "J2K POS (KM)");
            std::sscanf(data.c_str(), "%lf %lf %lf", &j2k_pos[0], &j2k_pos[1], &j2k_pos[2]);
            continue;
        }

        if (starts_with(l, "J2K VEL")) {
            auto data = extract_after(l, "J2K VEL (KM/S)");
            std::sscanf(data.c_str(), "%lf %lf %lf", &j2k_vel[0], &j2k_vel[1], &j2k_vel[2]);
            continue;
        }

        if (starts_with(l, "ECI POS")) {
            auto data = extract_after(l, "ECI POS (KM)");
            std::sscanf(data.c_str(), "%lf %lf %lf", &eci_pos[0], &eci_pos[1], &eci_pos[2]);
            continue;
        }

        if (starts_with(l, "ECI VEL")) {
            auto data = extract_after(l, "ECI VEL (KM/S)");
            std::sscanf(data.c_str(), "%lf %lf %lf", &eci_vel[0], &eci_vel[1], &eci_vel[2]);
            continue;
        }

        if (starts_with(l, "EFG POS")) {
            auto data = extract_after(l, "EFG POS (KM)");
            std::sscanf(data.c_str(), "%lf %lf %lf", &efg_pos[0], &efg_pos[1], &efg_pos[2]);
            continue;
        }

        if (starts_with(l, "EFG VEL")) {
            auto data = extract_after(l, "EFG VEL (KM/S)");
            std::sscanf(data.c_str(), "%lf %lf %lf", &efg_vel[0], &efg_vel[1], &efg_vel[2]);
            continue;
        }

        if (starts_with(l, "GEOPOTENTIAL")) {
            // "GEOPOTENTIAL: EGM-96 70Z,70T  DRAG: JAC70/MSIS90  LUNAR/SOLAR:  ON"
            auto geo_str = extract_between(l, "GEOPOTENTIAL", "DRAG");
            // Parse: "EGM-96 70Z,70T"
            if (!geo_str.empty()) {
                // Parse "EGM-96 70Z,70T" — find the "NNZ,NNT" pattern
                // Look for "Z," which marks the zonal/tesseral spec
                auto z_pos = geo_str.find("Z,");
                if (z_pos != std::string::npos) {
                    // Walk back to find the start of the degree number
                    auto deg_start = z_pos;
                    while (deg_start > 0 && std::isdigit(geo_str[deg_start - 1])) deg_start--;
                    geopotential_model = trim(geo_str.substr(0, deg_start));
                    std::sscanf(geo_str.c_str() + deg_start, "%dZ,%dT", &geo_degree, &geo_order);
                } else {
                    geopotential_model = geo_str;
                }
            }

            drag_model = extract_between(l, "DRAG", "LUNAR/SOLAR");
            auto ls = extract_after(l, "LUNAR/SOLAR");
            lunar_solar_on = (ls.find("ON") != std::string::npos);
            continue;
        }

        if (starts_with(l, "SOLAR RAD PRESS:") || starts_with(l, "SOLAR RAD PRESS ")) {
            if (l.find("COEFF") != std::string::npos) {
                // "SOLAR RAD PRESS COEFF (M2/KG):  0.000000E+00  EDR(W/KG):  0.39E-02"
                auto srp_str = extract_between(l, "SOLAR RAD PRESS COEFF (M2/KG)", "EDR");
                srp_coeff = parse_sci(srp_str);
                auto edr_str = extract_after(l, "EDR(W/KG)");
                edr = parse_sci(edr_str);
            } else {
                // "SOLAR RAD PRESS: OFF  SOLID EARTH TIDES:  ON  IN-TRACK THRUST: OFF"
                auto srp_str = extract_between(l, "SOLAR RAD PRESS", "SOLID EARTH");
                srp_on = (srp_str.find("ON") != std::string::npos);
                auto tides_str = extract_between(l, "SOLID EARTH TIDES", "IN-TRACK");
                solid_tides_on = (tides_str.find("ON") != std::string::npos);
                auto thrust_str = extract_after(l, "IN-TRACK THRUST");
                thrust_on = (thrust_str.find("ON") != std::string::npos);
            }
            continue;
        }

        if (starts_with(l, "BALLISTIC COEF")) {
            auto bc_str = extract_between(l, "BALLISTIC COEF (M2/KG)", "BDOT");
            ballistic_coeff = parse_sci(bc_str);
            auto bdot_str = extract_after(l, "BDOT (M2/KG-S)");
            bdot = parse_sci(bdot_str);
            continue;
        }

        if (starts_with(l, "THRUST ACCEL")) {
            auto ta_str = extract_between(l, "THRUST ACCEL (M/S2)", "C.M. OFFSET");
            thrust_accel = parse_sci(ta_str);
            auto cm_str = extract_after(l, "C.M. OFFSET (M)");
            cm_offset = parse_sci(cm_str);
            continue;
        }

        if (starts_with(l, "SOLAR FLUX")) {
            auto f10_str = extract_between(l, "F10", "AVERAGE F10");
            f10 = parse_double(f10_str);
            auto avg_f10_str = extract_between(l, "AVERAGE F10", "AVERAGE AP");
            avg_f10 = parse_double(avg_f10_str);
            auto ap_str = extract_after(l, "AVERAGE AP");
            avg_ap = parse_double(ap_str);
            continue;
        }

        if (starts_with(l, "TAI-UTC")) {
            auto tai_str = extract_between(l, "TAI-UTC (S)", "UT1-UTC");
            tai_utc = parse_int(tai_str);
            auto ut1_str = extract_between(l, "UT1-UTC (S)", "UT1 RATE");
            ut1_utc = parse_double(ut1_str);
            auto rate_str = extract_after(l, "UT1 RATE (MS/DAY)");
            ut1_rate = parse_double(rate_str);
            continue;
        }

        if (starts_with(l, "POLAR MOT")) {
            auto xy_str = extract_between(l, "POLAR MOT X,Y (ARCSEC)", "IAU");
            std::sscanf(xy_str.c_str(), "%lf %lf", &polar_x, &polar_y);
            auto nut_str = extract_after(l, "IAU 1980 NUTAT");
            nutation_terms = parse_int(nut_str);
            continue;
        }

        if (starts_with(l, "TIME CONST LEAP")) {
            auto ls_str = extract_after(l, "TIME CONST LEAP SECOND TIME (UTC)");
            leap_second_time = parse_vcm_date(ls_str);
            continue;
        }

        if (starts_with(l, "INTEGRATOR MODE")) {
            integrator_mode = extract_between(l, "INTEGRATOR MODE", "COORD SYS");
            coord_sys = extract_between(l, "COORD SYS", "PARTIALS");
            partials = extract_after(l, "PARTIALS");
            continue;
        }

        if (starts_with(l, "STEP MODE")) {
            step_mode = extract_between(l, "STEP MODE", "FIXED STEP");
            fixed_step = extract_between(l, "FIXED STEP", "STEP SIZE");
            step_selection = extract_after(l, "STEP SIZE SELECTION");
            continue;
        }

        if (starts_with(l, "INITIAL STEP SIZE")) {
            auto step_str = extract_between(l, "INITIAL STEP SIZE (S)", "ERROR CONTROL");
            init_step = parse_double(step_str);
            auto err_str = extract_after(l, "ERROR CONTROL");
            error_control = parse_sci(err_str);
            continue;
        }

        if (starts_with(l, "VECTOR U,V,W SIGMAS (KM):")) {
            auto data = extract_after(l, "VECTOR U,V,W SIGMAS (KM)");
            std::sscanf(data.c_str(), "%lf %lf %lf", &uvw_pos[0], &uvw_pos[1], &uvw_pos[2]);
            continue;
        }

        if (starts_with(l, "VECTOR UD,VD,WD SIGMAS")) {
            auto data = extract_after(l, "VECTOR UD,VD,WD SIGMAS (KM/S)");
            std::sscanf(data.c_str(), "%lf %lf %lf", &uvw_vel[0], &uvw_vel[1], &uvw_vel[2]);
            continue;
        }

        if (starts_with(l, "COVARIANCE MATRIX")) {
            // "COVARIANCE MATRIX (EQUINOCTIAL ELS): ( 7x 7) WTD RMS:  0.11498E+01"
            auto size_str = extract_between(l, "(EQUINOCTIAL ELS)", "WTD RMS");
            // Parse "( 7x 7)"
            auto paren = size_str.find('(');
            if (paren != std::string::npos) {
                std::sscanf(size_str.c_str() + paren, "(%dx%d)", &cov_rows, &cov_cols);
                if (cov_rows == 0) std::sscanf(size_str.c_str() + paren, "( %dx %d)", &cov_rows, &cov_cols);
            }
            auto rms_str = extract_after(l, "WTD RMS");
            wtd_rms = parse_sci(rms_str);
            in_covariance = true;
            continue;
        }

        // Covariance data lines (after COVARIANCE MATRIX header)
        if (in_covariance) {
            // Parse space-separated scientific notation values
            std::istringstream iss(l);
            double val;
            while (iss >> val) {
                covariance.push_back(val);
            }
            continue;
        }
    }

    // ── Build OCM FlatBuffer ──

    flatbuffers::FlatBufferBuilder fbb(4096);

    // Header
    auto header = CreateHeaderDirect(fbb,
        "2.0",          // CCSDS_OCM_VERS
        nullptr,        // COMMENT
        nullptr,        // CLASSIFICATION
        msg_time_utc.c_str(),  // CREATION_DATE
        center.c_str(), // ORIGINATOR
        nullptr         // MESSAGE_ID
    );

    // Metadata
    auto metadata = CreateMetadataDirect(fbb,
        nullptr,                    // COMMENT
        common_name.c_str(),        // OBJECT_NAME
        intl_des.c_str(),           // INTERNATIONAL_DESIGNATOR
        nullptr,                    // CATALOG_NAME
        std::to_string(sat_num).c_str(),  // OBJECT_DESIGNATOR
        nullptr,                    // ALTERNATE_NAMES
        nullptr, nullptr, nullptr, nullptr, nullptr,  // originator contact
        nullptr, nullptr, nullptr, nullptr, nullptr, nullptr,  // tech contact
        nullptr, nullptr,           // prev/next message IDs
        nullptr, nullptr, nullptr, nullptr, nullptr,  // links
        nullptr, nullptr, nullptr, nullptr,  // operator, owner, country, constellation
        nullptr,                    // OBJECT_TYPE
        "UTC",                      // TIME_SYSTEM
        epoch_utc.c_str(),          // EPOCH_TZERO
        nullptr,                    // OPS_STATUS
        nullptr,                    // ORBIT_CATEGORY
        nullptr,                    // OCM_DATA_ELEMENTS
        0.0, 0.0,                   // SCLK
        nullptr, nullptr,           // prev/next message epochs
        epoch_utc.c_str(),          // START_TIME
        epoch_utc.c_str(),          // STOP_TIME
        0.0,                        // TIME_SPAN
        static_cast<double>(tai_utc),  // TAIMUTC_AT_TZERO
        leap_second_time.empty() ? nullptr : leap_second_time.c_str(),
        0.0,                        // NEXT_LEAP_TAIMUTC
        ut1_utc                     // UT1MUTC_AT_TZERO
    );

    // State data (J2K, compact array: X,Y,Z,X_DOT,Y_DOT,Z_DOT)
    std::vector<double> state_data = {
        j2k_pos[0], j2k_pos[1], j2k_pos[2],
        j2k_vel[0], j2k_vel[1], j2k_vel[2]
    };

    // Physical properties
    auto phys = CreatePhysicalPropertiesDirect(fbb,
        nullptr,    // COMMENT
        0.0,        // WET_MASS
        0.0,        // DRY_MASS
        nullptr,    // MASS_UNITS
        0, 0, 0, 0, // OEB quaternion
        0, 0, 0,    // OEB dimensions
        0, 0, 0,    // areas along OEB
        nullptr,    // AREA_UNITS
        0.0,        // DRAG_CONST_AREA (not directly in VCM text)
        ballistic_coeff,  // DRAG_COEFF_NOM (ballistic coefficient M2/KG)
        0.0,        // DRAG_UNCERTAINTY
        0.0,        // SRP_CONST_AREA
        srp_coeff,  // SOLAR_RAD_COEFF
        0.0         // SRP_UNCERTAINTY
    );

    // Perturbations
    auto atm_model = CreateATM(fbb,
        map_atm_model(drag_model),
        map_atm_year(drag_model)
    );

    // N-body perturbations list
    std::vector<flatbuffers::Offset<flatbuffers::String>> nbody_vec;
    if (lunar_solar_on) {
        nbody_vec.push_back(fbb.CreateString("MOON"));
        nbody_vec.push_back(fbb.CreateString("SUN"));
    }

    auto perturbations = CreatePerturbationsDirect(fbb,
        nullptr,                    // COMMENT
        atm_model,                  // ATMOSPHERIC_MODEL
        geopotential_model.c_str(), // GRAVITY_MODEL
        geo_degree,                 // GRAVITY_DEGREE
        geo_order,                  // GRAVITY_ORDER
        398600.4418,                // GM
        lunar_solar_on ? &nbody_vec : nullptr,  // N_BODY_PERTURBATIONS
        nullptr,                    // OCEAN_TIDES_MODEL
        solid_tides_on ? "ON" : nullptr,  // SOLID_TIDES_MODEL
        nullptr,                    // ATMOSPHERIC_TIDES_MODEL
        nullptr,                    // GEOPOTENTIAL_MODEL (duplicates GRAVITY_MODEL)
        srp_on ? "ON" : "OFF",     // SOLAR_RAD_PRESSURE
        nullptr,                    // ALBEDO
        nullptr,                    // THERMAL
        nullptr,                    // RELATIVITY
        nullptr,                    // ATMOSPHERIC_DRAG
        avg_ap,                     // FIXED_GEOMAG_KP
        f10,                        // FIXED_F10P7
        avg_f10                     // FIXED_F10P7_MEAN
    );

    // User-defined parameters for fields not in OCM schema
    std::vector<flatbuffers::Offset<UserDefinedParameters>> user_params;

    // ECI state
    if (eci_pos[0] != 0 || eci_pos[1] != 0 || eci_pos[2] != 0) {
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%.8f %.8f %.8f", eci_pos[0], eci_pos[1], eci_pos[2]);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "ECI_POS_KM", buf));
        std::snprintf(buf, sizeof(buf), "%.14f %.14f %.14f", eci_vel[0], eci_vel[1], eci_vel[2]);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "ECI_VEL_KMS", buf));
    }

    // EFG state
    if (efg_pos[0] != 0 || efg_pos[1] != 0 || efg_pos[2] != 0) {
        char buf[256];
        std::snprintf(buf, sizeof(buf), "%.8f %.8f %.8f", efg_pos[0], efg_pos[1], efg_pos[2]);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "EFG_POS_KM", buf));
        std::snprintf(buf, sizeof(buf), "%.14f %.14f %.14f", efg_vel[0], efg_vel[1], efg_vel[2]);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "EFG_VEL_KMS", buf));
    }

    // UVW sigmas
    {
        char buf[128];
        std::snprintf(buf, sizeof(buf), "%.4f %.4f %.4f", uvw_pos[0], uvw_pos[1], uvw_pos[2]);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "UVW_POS_SIGMA_KM", buf));
        std::snprintf(buf, sizeof(buf), "%.4f %.4f %.4f", uvw_vel[0], uvw_vel[1], uvw_vel[2]);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "UVW_VEL_SIGMA_KMS", buf));
    }

    // BDOT, EDR, thrust, CM offset
    {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.6e", bdot);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "BDOT_M2_KG_S", buf));
        std::snprintf(buf, sizeof(buf), "%.2e", edr);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "EDR_W_KG", buf));
        std::snprintf(buf, sizeof(buf), "%.6e", thrust_accel);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "THRUST_ACCEL_MS2", buf));
        std::snprintf(buf, sizeof(buf), "%.6e", cm_offset);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "CM_OFFSET_M", buf));
    }

    // Integrator settings
    if (!integrator_mode.empty()) {
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "INTEGRATOR_MODE", integrator_mode.c_str()));
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "COORD_SYS", coord_sys.c_str()));
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "PARTIALS", partials.c_str()));
    }
    if (!step_mode.empty()) {
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "STEP_MODE", step_mode.c_str()));
    }
    {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.3f", init_step);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "INIT_STEP_SIZE_S", buf));
        std::snprintf(buf, sizeof(buf), "%.3e", error_control);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "ERROR_CONTROL", buf));
    }

    // Epoch rev
    {
        user_params.push_back(CreateUserDefinedParametersDirect(fbb,
            "EPOCH_REV", std::to_string(epoch_rev).c_str()));
    }

    // Weighted RMS
    if (wtd_rms != 0) {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.5e", wtd_rms);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "WTD_RMS", buf));
    }

    // Polar motion
    {
        char buf[64];
        std::snprintf(buf, sizeof(buf), "%.4f %.4f", polar_x, polar_y);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "POLAR_MOTION_XY_ARCSEC", buf));
        std::snprintf(buf, sizeof(buf), "%.3f", ut1_rate);
        user_params.push_back(CreateUserDefinedParametersDirect(fbb, "UT1_RATE_MS_DAY", buf));
    }

    // NORAD CAT ID
    user_params.push_back(CreateUserDefinedParametersDirect(fbb,
        "NORAD_CAT_ID", std::to_string(sat_num).c_str()));

    // Build OCM
    auto ocm = CreateOCM(fbb,
        header,
        metadata,
        fbb.CreateString("ESTIMATED"),  // TRAJ_TYPE
        0.0,                             // STATE_STEP_SIZE (single point)
        6,                               // STATE_VECTOR_SIZE
        fbb.CreateVector(state_data),    // STATE_DATA
        covariance.empty() ? 0 : fbb.CreateVector(covariance),  // COVARIANCE_DATA
        phys,                            // PHYSICAL_PROPERTIES
        0,                               // MANEUVER_DATA
        perturbations,                   // PERTURBATIONS
        0,                               // ORBIT_DETERMINATION
        fbb.CreateVector(user_params)    // USER_DEFINED_PARAMETERS
    );

    fbb.Finish(ocm, "$OCM");

    // Copy to result
    result.ocm_buffer.assign(
        fbb.GetBufferPointer(),
        fbb.GetBufferPointer() + fbb.GetSize()
    );
    result.valid = true;
    result.object_name = common_name;
    result.intl_designator = intl_des;
    result.norad_id = sat_num;
    result.epoch_utc = epoch_utc;

    return result;
}

}  // namespace hpop
