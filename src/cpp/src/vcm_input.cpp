/**
 * HPOP VCM Input Parser — Implementation
 *
 * Parses VCM FlatBuffer binary and JSON into PropagationInput.
 * Maps VCM atmospheric model enums to HPOP atmosphere model strings.
 * Maps VCM propagator config to HPOP force model selection.
 */

#include "hpop/vcm_input.h"
#include "sds/vcm_generated.h"

#include <cstring>
#include <cstdlib>
#include <cmath>

namespace hpop {

// ── Time conversion (shared with wasm_api.cpp) ──

static double iso8601_to_jd_internal(const std::string& iso) {
    int year, month, day, hour = 0, min = 0;
    double sec = 0.0;
    if (std::sscanf(iso.c_str(), "%d-%d-%dT%d:%d:%lf",
                    &year, &month, &day, &hour, &min, &sec) < 3) {
        return 0.0;
    }
    if (month <= 2) { year--; month += 12; }
    int A = year / 100;
    int B = 2 - A + A / 4;
    return std::floor(365.25 * (year + 4716))
         + std::floor(30.6001 * (month + 1))
         + day + B - 1524.5
         + (hour + min / 60.0 + sec / 3600.0) / 24.0;
}

// ── Simple JSON helpers ──

static std::string json_str(const std::string& json, const std::string& key) {
    std::string search = "\"" + key + "\"";
    auto pos = json.find(search);
    if (pos == std::string::npos) return "";
    pos = json.find(':', pos + search.size());
    if (pos == std::string::npos) return "";
    pos = json.find('"', pos + 1);
    if (pos == std::string::npos) return "";
    auto end = json.find('"', pos + 1);
    if (end == std::string::npos) return "";
    return json.substr(pos + 1, end - pos - 1);
}

static double json_num(const std::string& json, const std::string& key, double def = 0.0) {
    std::string search = "\"" + key + "\"";
    auto pos = json.find(search);
    if (pos == std::string::npos) return def;
    pos = json.find(':', pos + search.size());
    if (pos == std::string::npos) return def;
    pos++;
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t')) pos++;
    return std::strtod(json.c_str() + pos, nullptr);
}

static bool json_has(const std::string& json, const std::string& key) {
    return json.find("\"" + key + "\"") != std::string::npos;
}

// ── VCM FlatBuffer Parser ──

PropagationInput parse_vcm_flatbuffer(const uint8_t* data, size_t size) {
    PropagationInput input;

    if (size < 8) {
        input.error = "Buffer too small for FlatBuffer";
        return input;
    }

    // Verify file identifier
    if (!flatbuffers::BufferHasIdentifier(data, "VCM\0") &&
        !flatbuffers::BufferHasIdentifier(data, "$VCM")) {
        // Try without identifier
    }

    auto* vcm = flatbuffers::GetRoot<VCM>(data);
    if (!vcm) {
        input.error = "Failed to parse VCM FlatBuffer";
        return input;
    }

    // Metadata
    if (vcm->OBJECT_NAME()) input.object_name = vcm->OBJECT_NAME()->str();
    if (vcm->OBJECT_ID()) input.object_id = vcm->OBJECT_ID()->str();
    if (vcm->ORIGINATOR()) input.originator = vcm->ORIGINATOR()->str();
    if (vcm->REF_FRAME()) input.ref_frame = vcm->REF_FRAME()->str();
    input.norad_cat_id = vcm->NORAD_CAT_ID();

    // State vector
    auto* sv = vcm->STATE_VECTOR();
    if (sv) {
        input.initial_state.r[0] = sv->X();
        input.initial_state.r[1] = sv->Y();
        input.initial_state.r[2] = sv->Z();
        input.initial_state.v[0] = sv->X_DOT();
        input.initial_state.v[1] = sv->Y_DOT();
        input.initial_state.v[2] = sv->Z_DOT();

        if (sv->EPOCH()) {
            input.epoch_iso = sv->EPOCH()->str();
            input.epoch_jd = iso8601_to_jd_internal(input.epoch_iso);
        }
    }

    // Keplerian elements
    auto* kep = vcm->KEPLERIAN_ELEMENTS();
    if (kep) {
        input.has_keplerian = true;
        input.keplerian.a = kep->SEMI_MAJOR_AXIS();
        input.keplerian.e = kep->ECCENTRICITY();
        input.keplerian.i = kep->INCLINATION() * DEG2RAD;
        input.keplerian.raan = kep->RA_OF_ASC_NODE() * DEG2RAD;
        input.keplerian.argp = kep->ARG_OF_PERICENTER() * DEG2RAD;
        input.keplerian.nu = kep->ANOMALY() * DEG2RAD;

        // If no state vector, convert from Keplerian
        if (!sv || (sv->X() == 0 && sv->Y() == 0 && sv->Z() == 0)) {
            // Handle mean anomaly → true anomaly conversion if needed
            if (kep->ANOMALY_TYPE() == anomalyType_MEAN_ANOMALY) {
                // Simple Newton iteration: M → E → ν
                double M = kep->ANOMALY() * DEG2RAD;
                double e = kep->ECCENTRICITY();
                double E = M;
                for (int i = 0; i < 20; i++) {
                    double dE = (E - e * std::sin(E) - M) / (1.0 - e * std::cos(E));
                    E -= dE;
                    if (std::abs(dE) < 1e-14) break;
                }
                double nu = 2.0 * std::atan2(
                    std::sqrt(1.0 + e) * std::sin(E / 2.0),
                    std::sqrt(1.0 - e) * std::cos(E / 2.0)
                );
                input.keplerian.nu = nu;
            }
            input.initial_state = keplerian_to_cartesian(input.keplerian);
        }
    }

    // Spacecraft properties
    if (vcm->MASS() > 0) input.spacecraft.mass_kg = vcm->MASS();
    if (vcm->DRAG_AREA() > 0) input.spacecraft.drag_area = vcm->DRAG_AREA();
    if (vcm->DRAG_COEFF() > 0) input.spacecraft.cd = vcm->DRAG_COEFF();
    if (vcm->SOLAR_RAD_AREA() > 0) input.spacecraft.srp_area = vcm->SOLAR_RAD_AREA();
    if (vcm->SOLAR_RAD_COEFF() > 0) input.spacecraft.cr = vcm->SOLAR_RAD_COEFF();

    // SRP on/off
    input.use_srp = (vcm->SRP() == perturbationStatus_ON);

    // Atmospheric model data
    auto* atm = vcm->ATMOSPHERIC_MODEL_DATA();
    if (atm) {
        switch (atm->ATMOSPHERIC_MODEL()) {
            case atmosphericModel_JACCHIA_70:
                input.atmosphere_model = "jacchia-70";
                break;
            case atmosphericModel_JB2008:
                input.atmosphere_model = "jb2008";
                break;
            case atmosphericModel_HASDM:
                input.atmosphere_model = "hasdm";
                break;
            case atmosphericModel_NRLMSISE_00:
                input.atmosphere_model = "nrlmsise-00";
                break;
            case atmosphericModel_NONE:
                input.use_drag = false;
                break;
            default:
                input.atmosphere_model = "harris-priester";
                break;
        }

        // Lunar/solar perturbations
        if (atm->LUNAR_SOLAR_PERTURBATION() == perturbationStatus_OFF) {
            input.use_sun = false;
            input.use_moon = false;
        }

        // SRP from atmospheric model data
        if (atm->SOLAR_RADIATION_PRESSURE() == perturbationStatus_OFF) {
            input.use_srp = false;
        }
    }

    // Propagator settings
    auto* prop = vcm->PROPAGATOR_SETTINGS();
    if (prop) {
        if (prop->TIME_STEP() > 0) input.step_seconds = prop->TIME_STEP();

        // Zonal harmonics → geopotential config
        auto* zonals = prop->ZONAL_HARMONIC_TERMS();
        if (zonals && zonals->size() > 0) {
            input.use_j2 = false;  // Use full geopotential instead
            input.use_geopotential = true;
            // Find max degree
            int max_degree = 0;
            for (size_t i = 0; i < zonals->size(); i++) {
                int d = static_cast<int>(zonals->Get(i));
                if (d > max_degree) max_degree = d;
            }
            input.geopotential_degree = max_degree + 1;  // J2=degree 2, etc.
            input.geopotential_order = 0;  // Zonal only
        }
    }

    input.valid = true;
    return input;
}

// ── JSON Parser ──

PropagationInput parse_vcm_json(const std::string& json) {
    PropagationInput input;

    // Epoch
    input.epoch_iso = json_str(json, "epoch");
    if (input.epoch_iso.empty()) {
        // Try VCM-style nested EPOCH in STATE_VECTOR
        input.epoch_iso = json_str(json, "EPOCH");
    }
    if (!input.epoch_iso.empty()) {
        input.epoch_jd = iso8601_to_jd_internal(input.epoch_iso);
    }

    // State vector (flat or nested)
    if (json_has(json, "x") || json_has(json, "X")) {
        input.initial_state.r[0] = json_has(json, "x") ? json_num(json, "x") : json_num(json, "X");
        input.initial_state.r[1] = json_has(json, "y") ? json_num(json, "y") : json_num(json, "Y");
        input.initial_state.r[2] = json_has(json, "z") ? json_num(json, "z") : json_num(json, "Z");
        input.initial_state.v[0] = json_has(json, "vx") ? json_num(json, "vx") : json_num(json, "X_DOT");
        input.initial_state.v[1] = json_has(json, "vy") ? json_num(json, "vy") : json_num(json, "Y_DOT");
        input.initial_state.v[2] = json_has(json, "vz") ? json_num(json, "vz") : json_num(json, "Z_DOT");
    }

    // Keplerian elements
    if (json_has(json, "SEMI_MAJOR_AXIS") || json_has(json, "semi_major_axis")) {
        input.has_keplerian = true;
        input.keplerian.a = json_has(json, "SEMI_MAJOR_AXIS")
            ? json_num(json, "SEMI_MAJOR_AXIS") : json_num(json, "semi_major_axis");
        input.keplerian.e = json_has(json, "ECCENTRICITY")
            ? json_num(json, "ECCENTRICITY") : json_num(json, "eccentricity");
        input.keplerian.i = json_num(json, "INCLINATION", json_num(json, "inclination")) * DEG2RAD;
        input.keplerian.raan = json_num(json, "RA_OF_ASC_NODE", json_num(json, "raan")) * DEG2RAD;
        input.keplerian.argp = json_num(json, "ARG_OF_PERICENTER", json_num(json, "argp")) * DEG2RAD;
        input.keplerian.nu = json_num(json, "ANOMALY", json_num(json, "anomaly")) * DEG2RAD;

        // Convert to Cartesian if no state vector provided
        if (input.initial_state.r[0] == 0 && input.initial_state.r[1] == 0) {
            input.initial_state = keplerian_to_cartesian(input.keplerian);
        }
    }

    // Spacecraft properties
    input.spacecraft.mass_kg = json_num(json, "mass_kg", json_num(json, "MASS", 260.0));
    input.spacecraft.drag_area = json_num(json, "drag_area", json_num(json, "DRAG_AREA", 22.0));
    input.spacecraft.cd = json_num(json, "cd", json_num(json, "DRAG_COEFF", 2.2));
    input.spacecraft.srp_area = json_num(json, "srp_area", json_num(json, "SOLAR_RAD_AREA", 22.0));
    input.spacecraft.cr = json_num(json, "cr", json_num(json, "SOLAR_RAD_COEFF", 1.3));

    // Propagation settings
    input.duration_days = json_num(json, "duration_days", 1.0);
    input.step_seconds = json_num(json, "step_seconds", json_num(json, "TIME_STEP", 60.0));

    // Atmosphere
    input.f107 = json_num(json, "f107", 150.0);
    std::string atm = json_str(json, "atmosphere_model");
    if (atm.empty()) atm = json_str(json, "model");
    if (!atm.empty()) input.atmosphere_model = atm;

    // Metadata
    input.object_name = json_str(json, "OBJECT_NAME");
    if (input.object_name.empty()) input.object_name = json_str(json, "object_name");
    input.object_id = json_str(json, "OBJECT_ID");
    input.norad_cat_id = static_cast<uint32_t>(json_num(json, "NORAD_CAT_ID", 0));

    input.valid = (input.epoch_jd > 0 &&
                   (input.initial_state.rmag() > 0 || input.has_keplerian));

    if (!input.valid) {
        input.error = "Missing epoch or state vector/keplerian elements";
    }

    return input;
}

// ── Auto-detect and parse ──

PropagationInput parse_vcm(const uint8_t* data, size_t size) {
    // Try FlatBuffer first (check for file identifier)
    if (size >= 8) {
        // FlatBuffer file identifiers are at offset 4-7
        const char* id = reinterpret_cast<const char*>(data + 4);
        if (std::strncmp(id, "$VCM", 4) == 0 || std::strncmp(id, "VCM\0", 4) == 0) {
            return parse_vcm_flatbuffer(data, size);
        }
    }

    // Try JSON
    std::string text(reinterpret_cast<const char*>(data), size);
    // Check if it looks like JSON
    for (size_t i = 0; i < text.size(); i++) {
        char c = text[i];
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') continue;
        if (c == '{') return parse_vcm_json(text);
        break;
    }

    PropagationInput input;
    input.error = "Unrecognized input format (expected VCM FlatBuffer or JSON)";
    return input;
}

}  // namespace hpop
