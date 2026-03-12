#pragma once
/**
 * HPOP VCM Input Parser
 *
 * Parses Space Force Vector Covariance Messages (VCM) and extracts
 * everything needed to configure and run an HPOP propagation:
 *
 *   - Initial state vector (position + velocity)
 *   - Spacecraft properties (mass, drag area, Cd, SRP area, Cr)
 *   - Atmospheric model selection
 *   - Propagator configuration (force models, step size)
 *   - Keplerian elements (as alternative to state vector)
 *
 * Input formats:
 *   1. VCM FlatBuffer binary ($VCM file identifier)
 *   2. JSON (same fields as VCM schema)
 *
 * The parser populates a PropagationInput struct that maps directly
 * to the HPOP propagator API.
 *
 * Reference: CCSDS 502.0-B-3 (Orbit Data Messages)
 * Schema: spacedatastandards.org/VCM
 */

#include "hpop/state.h"
#include <string>
#include <vector>

namespace hpop {

// ── Propagation Input (parsed from VCM) ──

struct PropagationInput {
    // Initial state
    State initial_state;
    double epoch_jd = 0.0;
    std::string epoch_iso;

    // Spacecraft properties
    SpacecraftProperties spacecraft;

    // Force model configuration
    bool use_twobody = true;
    bool use_j2 = true;
    bool use_geopotential = false;
    int geopotential_degree = 4;
    int geopotential_order = 4;
    bool use_drag = true;
    bool use_srp = true;
    bool use_sun = true;
    bool use_moon = true;

    // Atmosphere
    std::string atmosphere_model = "harris-priester";  // or "jacchia-70", "jb2008"
    double f107 = 150.0;

    // Propagation settings
    double step_seconds = 60.0;
    double duration_days = 1.0;

    // Metadata
    std::string object_name;
    std::string object_id;
    uint32_t norad_cat_id = 0;
    std::string ref_frame = "TEME";
    std::string originator;

    // Keplerian elements (alternative representation)
    bool has_keplerian = false;
    KeplerianElements keplerian;

    // Validity
    bool valid = false;
    std::string error;
};

/**
 * Parse VCM from FlatBuffer binary.
 * Detects $VCM file identifier.
 */
PropagationInput parse_vcm_flatbuffer(const uint8_t* data, size_t size);

/**
 * Parse VCM from JSON string.
 * Accepts both VCM-schema JSON and simplified propagation JSON.
 */
PropagationInput parse_vcm_json(const std::string& json);

/**
 * Auto-detect format and parse.
 * Tries FlatBuffer first (checks file identifier), falls back to JSON.
 */
PropagationInput parse_vcm(const uint8_t* data, size_t size);

}  // namespace hpop
