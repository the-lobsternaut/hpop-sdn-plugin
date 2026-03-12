/**
 * HPOP WASM API
 *
 * SDN Plugin interface: accepts initial state (JSON or FlatBuffers OEM),
 * propagates the orbit using the full HPOP force model, and outputs
 * OEM FlatBuffers aligned binary.
 *
 * Exports:
 *   C ABI:   parse(), convert(), malloc(), free()
 *   Embind:  propagate(), propagateJSON(), getForceModels(), getVersion()
 *
 * Input JSON format:
 * {
 *   "epoch": "2024-01-01T00:00:00.000Z",
 *   "state": { "x": 6778.0, "y": 0, "z": 0, "vx": 0, "vy": 7.669, "vz": 0 },
 *   "duration_days": 1.0,
 *   "step_seconds": 60.0,
 *   "spacecraft": { "mass_kg": 260, "drag_area": 22, "cd": 2.2, "srp_area": 22, "cr": 1.3 },
 *   "forces": ["twobody", "j2", "drag", "srp", "sun", "moon"],
 *   "atmosphere": { "model": "harris-priester", "f107": 150.0 },
 *   "output_format": "oem"
 * }
 */

#include "hpop/propagator.h"
#include "hpop/atmosphere.h"
#include "hpop/srp.h"
#include "hpop/thirdbody.h"

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#include <emscripten/bind.h>
using namespace emscripten;
#endif

// OEM FlatBuffers generated header
#include "main_generated.h"

#include <string>
#include <vector>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <cmath>

using namespace hpop;

// ============================================================================
// Time Utilities
// ============================================================================

// Parse ISO 8601 date to Julian Date (TDB)
// Handles: "2024-01-01T00:00:00.000Z" format
static double iso8601_to_jd(const std::string& iso) {
    int year, month, day, hour, min;
    double sec;
    if (std::sscanf(iso.c_str(), "%d-%d-%dT%d:%d:%lf",
                    &year, &month, &day, &hour, &min, &sec) < 5) {
        // Try date-only
        if (std::sscanf(iso.c_str(), "%d-%d-%d", &year, &month, &day) < 3) {
            return 0.0;
        }
        hour = min = 0; sec = 0.0;
    }

    // Algorithm from Meeus, "Astronomical Algorithms" 2nd ed, p. 61
    if (month <= 2) { year--; month += 12; }
    int A = year / 100;
    int B = 2 - A + A / 4;
    double JD = std::floor(365.25 * (year + 4716))
              + std::floor(30.6001 * (month + 1))
              + day + B - 1524.5
              + (hour + min / 60.0 + sec / 3600.0) / 24.0;
    return JD;
}

// Convert Julian Date to ISO 8601 string
static std::string jd_to_iso8601(double jd) {
    // Algorithm from Meeus
    double jd_plus = jd + 0.5;
    int Z = static_cast<int>(jd_plus);
    double F = jd_plus - Z;

    int A;
    if (Z < 2299161) { A = Z; }
    else {
        int alpha = static_cast<int>((Z - 1867216.25) / 36524.25);
        A = Z + 1 + alpha - alpha / 4;
    }

    int B = A + 1524;
    int C = static_cast<int>((B - 122.1) / 365.25);
    int D = static_cast<int>(365.25 * C);
    int E = static_cast<int>((B - D) / 30.6001);

    int day = B - D - static_cast<int>(30.6001 * E);
    int month = (E < 14) ? E - 1 : E - 13;
    int year = (month > 2) ? C - 4716 : C - 4715;

    double total_seconds = F * 86400.0;
    int hour = static_cast<int>(total_seconds / 3600.0);
    total_seconds -= hour * 3600.0;
    int min = static_cast<int>(total_seconds / 60.0);
    double sec = total_seconds - min * 60.0;
    // Handle rounding to 60.000
    if (sec >= 59.9995) { sec = 0.0; min++; }
    if (min >= 60) { min -= 60; hour++; }
    if (hour >= 24) { hour -= 24; day++; }

    char buf[64];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%06.3fZ",
                  year, month, day, hour, min, sec);
    return std::string(buf);
}

// ============================================================================
// Simple JSON Parser (minimal, no external dependency)
// ============================================================================

// Extract a string value for a key from JSON
static std::string json_get_string(const std::string& json, const std::string& key) {
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

// Extract a numeric value for a key from JSON
static double json_get_number(const std::string& json, const std::string& key, double def = 0.0) {
    std::string search = "\"" + key + "\"";
    auto pos = json.find(search);
    if (pos == std::string::npos) return def;
    pos = json.find(':', pos + search.size());
    if (pos == std::string::npos) return def;
    // Skip whitespace
    pos++;
    while (pos < json.size() && (json[pos] == ' ' || json[pos] == '\t')) pos++;
    return std::strtod(json.c_str() + pos, nullptr);
}

// Check if a string array contains a value
static bool json_array_contains(const std::string& json, const std::string& key, const std::string& value) {
    std::string search = "\"" + key + "\"";
    auto pos = json.find(search);
    if (pos == std::string::npos) return true;  // Default: include all forces
    auto arr_start = json.find('[', pos);
    if (arr_start == std::string::npos) return true;
    auto arr_end = json.find(']', arr_start);
    if (arr_end == std::string::npos) return true;
    std::string arr = json.substr(arr_start, arr_end - arr_start + 1);
    return arr.find("\"" + value + "\"") != std::string::npos;
}

// ============================================================================
// Propagation Engine
// ============================================================================

struct PropagationRequest {
    State initial_state;
    double epoch_jd0;
    double duration_days;
    double step_seconds;
    SpacecraftProperties spacecraft;
    bool use_twobody = true;
    bool use_j2 = true;
    bool use_drag = true;
    bool use_srp = true;
    bool use_sun = true;
    bool use_moon = true;
    double f107 = 150.0;
    std::string atm_model = "harris-priester";
};

struct PropagationResult {
    std::vector<EphemerisPoint> ephemeris;
    double epoch_jd0;
    double epoch_jdf;
};

static PropagationResult run_propagation(const PropagationRequest& req) {
    Propagator prop;

    // Build force model stack
    if (req.use_twobody) {
        prop.add_force(std::make_unique<TwoBodyForce>());
    }
    if (req.use_j2) {
        prop.add_force(std::make_unique<J2Force>());
    }
    if (req.use_drag) {
        std::unique_ptr<AtmosphereModel> atm;
        if (req.atm_model == "exponential") {
            atm = std::make_unique<ExponentialAtmosphere>(3.725e-12, 400.0, 58.515);
        } else {
            atm = std::make_unique<HarrisPriester>(req.f107, 2);
        }
        prop.add_force(std::make_unique<DragForce>(std::move(atm), req.spacecraft.cd));
    }
    if (req.use_srp) {
        prop.add_force(std::make_unique<SRPForce>(ShadowType::CONICAL));
    }
    if (req.use_sun) {
        prop.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::SUN));
    }
    if (req.use_moon) {
        prop.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::MOON));
    }

    prop.set_spacecraft(req.spacecraft);

    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF78;
    cfg.step_size = 30.0;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-12;
    prop.set_config(cfg);

    double epoch_jdf = req.epoch_jd0 + req.duration_days;
    std::vector<EphemerisPoint> ephemeris;

    prop.propagate(req.initial_state, req.epoch_jd0, epoch_jdf,
                   req.step_seconds, &ephemeris);

    return {std::move(ephemeris), req.epoch_jd0, epoch_jdf};
}

// ============================================================================
// FlatBuffers OEM Output
// ============================================================================

static std::string build_oem_flatbuffer(const PropagationResult& result,
                                         const std::string& object_name = "HPOP_SAT") {
    flatbuffers::FlatBufferBuilder fbb(64 * 1024);  // 64KB initial

    // Build ephemeris data lines
    std::vector<flatbuffers::Offset<ephemerisDataLine>> lines;
    lines.reserve(result.ephemeris.size());

    for (const auto& pt : result.ephemeris) {
        auto epoch_str = fbb.CreateString(jd_to_iso8601(pt.epoch_jd));
        auto line = CreateephemerisDataLine(fbb,
            epoch_str,
            pt.state.r[0], pt.state.r[1], pt.state.r[2],
            pt.state.v[0], pt.state.v[1], pt.state.v[2]
        );
        lines.push_back(line);
    }

    auto lines_vec = fbb.CreateVector(lines);

    // Build ephemeris data block
    auto center_name = fbb.CreateString("EARTH");
    auto start_time = fbb.CreateString(jd_to_iso8601(result.epoch_jd0));
    auto stop_time = fbb.CreateString(jd_to_iso8601(result.epoch_jdf));
    auto interp = fbb.CreateString("HERMITE");

    auto block = CreateephemerisDataBlock(fbb,
        0,                          // COMMENT
        0,                          // OBJECT (CAT)
        center_name,                // CENTER_NAME
        0,                          // REFERENCE_FRAME (RFM)
        0,                          // REFERENCE_FRAME_EPOCH
        0,                          // COV_REFERENCE_FRAME
        timeSystem_UTC,             // TIME_SYSTEM
        start_time,                 // START_TIME
        0,                          // USEABLE_START_TIME
        0,                          // USEABLE_STOP_TIME
        stop_time,                  // STOP_TIME
        interp,                     // INTERPOLATION
        7,                          // INTERPOLATION_DEGREE
        0.0,                        // STEP_SIZE (0 = non-uniform, use lines)
        6,                          // STATE_VECTOR_SIZE
        0,                          // EPHEMERIS_DATA (compact)
        lines_vec,                  // EPHEMERIS_DATA_LINES
        0                           // COVARIANCE_MATRIX_LINES
    );

    std::vector<flatbuffers::Offset<ephemerisDataBlock>> blocks = {block};
    auto blocks_vec = fbb.CreateVector(blocks);

    // Build OEM
    auto creation_date = fbb.CreateString(jd_to_iso8601(result.epoch_jd0));
    auto originator = fbb.CreateString("HPOP-SDN-PLUGIN");

    auto oem = CreateOEM(fbb,
        0,                // CLASSIFICATION
        2.0,              // CCSDS_OEM_VERS
        creation_date,
        originator,
        blocks_vec
    );

    fbb.Finish(oem, "$OEM");

    const uint8_t* buf = fbb.GetBufferPointer();
    size_t size = fbb.GetSize();
    return std::string(reinterpret_cast<const char*>(buf), size);
}

// ============================================================================
// JSON API
// ============================================================================

static PropagationRequest parse_json_request(const std::string& json) {
    PropagationRequest req;

    // Parse epoch
    std::string epoch_str = json_get_string(json, "epoch");
    if (!epoch_str.empty()) {
        req.epoch_jd0 = iso8601_to_jd(epoch_str);
    } else {
        req.epoch_jd0 = 2451545.0;  // J2000 default
    }

    // Parse state
    req.initial_state.r[0] = json_get_number(json, "x", 6778.0);
    req.initial_state.r[1] = json_get_number(json, "y", 0.0);
    req.initial_state.r[2] = json_get_number(json, "z", 0.0);
    req.initial_state.v[0] = json_get_number(json, "vx", 0.0);
    req.initial_state.v[1] = json_get_number(json, "vy", 7.669);
    req.initial_state.v[2] = json_get_number(json, "vz", 0.0);

    // Parse duration and step
    req.duration_days = json_get_number(json, "duration_days", 1.0);
    req.step_seconds = json_get_number(json, "step_seconds", 60.0);

    // Parse spacecraft properties
    req.spacecraft.mass_kg = json_get_number(json, "mass_kg", 260.0);
    req.spacecraft.drag_area = json_get_number(json, "drag_area", 22.0);
    req.spacecraft.cd = json_get_number(json, "cd", 2.2);
    req.spacecraft.srp_area = json_get_number(json, "srp_area", 22.0);
    req.spacecraft.cr = json_get_number(json, "cr", 1.3);

    // Parse force selection
    req.use_twobody = json_array_contains(json, "forces", "twobody");
    req.use_j2 = json_array_contains(json, "forces", "j2");
    req.use_drag = json_array_contains(json, "forces", "drag");
    req.use_srp = json_array_contains(json, "forces", "srp");
    req.use_sun = json_array_contains(json, "forces", "sun");
    req.use_moon = json_array_contains(json, "forces", "moon");

    // Atmosphere settings
    req.f107 = json_get_number(json, "f107", 150.0);
    req.atm_model = json_get_string(json, "model");
    if (req.atm_model.empty()) req.atm_model = "harris-priester";

    return req;
}

// Propagate from JSON input, return OEM FlatBuffer binary
static std::string propagate_json(const std::string& json) {
    auto req = parse_json_request(json);
    auto result = run_propagation(req);
    return build_oem_flatbuffer(result);
}

// Get version info
static std::string get_version() {
    return "HPOP-SDN-Plugin 1.0.0";
}

// Get available force models
static std::string get_force_models() {
    return R"(["twobody","j2","drag","srp","sun","moon"])";
}

// ============================================================================
// SDN Plugin ABI: extern "C" exports
// ============================================================================

#ifdef __EMSCRIPTEN__

extern "C" {

/**
 * Parse JSON propagation request → OEM FlatBuffer binary
 */
EMSCRIPTEN_KEEPALIVE
int32_t parse(const char* input, size_t input_len,
              uint8_t* output, size_t output_len) {
    try {
        std::string json(input, input_len);
        std::string fb = propagate_json(json);

        if (fb.size() > output_len) return -2;
        std::memcpy(output, fb.data(), fb.size());
        return static_cast<int32_t>(fb.size());
    } catch (const std::exception&) {
        return -1;
    }
}

/**
 * Convert OEM FlatBuffer → CCSDS OEM text
 */
EMSCRIPTEN_KEEPALIVE
int32_t convert(const uint8_t* input, size_t input_len,
                uint8_t* output, size_t output_len) {
    try {
        flatbuffers::Verifier verifier(input, input_len);
        auto* oem = flatbuffers::GetRoot<OEM>(input);
        if (!oem) return -3;

        std::ostringstream oss;
        oss << std::fixed;
        oss << "CCSDS_OEM_VERS = " << oem->CCSDS_OEM_VERS() << "\n";
        if (oem->CREATION_DATE()) oss << "CREATION_DATE = " << oem->CREATION_DATE()->str() << "\n";
        if (oem->ORIGINATOR()) oss << "ORIGINATOR = " << oem->ORIGINATOR()->str() << "\n\n";

        auto* blocks = oem->EPHEMERIS_DATA_BLOCK();
        if (blocks) {
            for (size_t b = 0; b < blocks->size(); b++) {
                auto* block = blocks->Get(b);
                if (!block) continue;

                oss << "META_START\n";
                if (block->CENTER_NAME()) oss << "CENTER_NAME = " << block->CENTER_NAME()->str() << "\n";
                if (block->START_TIME()) oss << "START_TIME = " << block->START_TIME()->str() << "\n";
                if (block->STOP_TIME()) oss << "STOP_TIME = " << block->STOP_TIME()->str() << "\n";
                oss << "META_STOP\n\n";

                auto* lines = block->EPHEMERIS_DATA_LINES();
                if (lines) {
                    for (size_t i = 0; i < lines->size(); i++) {
                        auto* line = lines->Get(i);
                        if (!line || !line->EPOCH()) continue;
                        oss << line->EPOCH()->str()
                            << std::setprecision(6)
                            << " " << line->X()
                            << " " << line->Y()
                            << " " << line->Z()
                            << std::setprecision(9)
                            << " " << line->X_DOT()
                            << " " << line->Y_DOT()
                            << " " << line->Z_DOT()
                            << "\n";
                    }
                }
                oss << "\n";
            }
        }

        std::string text = oss.str();
        if (text.size() > output_len) return -2;
        std::memcpy(output, text.data(), text.size());
        return static_cast<int32_t>(text.size());
    } catch (const std::exception&) {
        return -1;
    }
}

} // extern "C"

// ============================================================================
// Embind API
// ============================================================================

EMSCRIPTEN_BINDINGS(hpop) {
    function("propagate", &propagate_json);
    function("getVersion", &get_version);
    function("getForceModels", &get_force_models);
}

#endif // __EMSCRIPTEN__
