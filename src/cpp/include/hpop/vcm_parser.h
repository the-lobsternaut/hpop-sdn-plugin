#pragma once
/**
 * Raw Text VCM → OCM FlatBuffer Parser
 *
 * Parses the 18 SPCS fixed-format ASCII VCM (V2.0) into a
 * spacedatastandards.org OCM ($OCM) FlatBuffer binary.
 *
 * Input: Raw text VCM (JANAP-128 format with <> line prefixes)
 * Output: OCM FlatBuffer binary ($OCM file identifier)
 *
 * Fields mapped:
 *   VCM line → OCM field
 *   ─────────────────────────────────────────────
 *   SATELLITE NUMBER       → Metadata.OBJECT_DESIGNATOR
 *   INT. DES.              → Metadata.INTERNATIONAL_DESIGNATOR
 *   COMMON NAME            → Metadata.OBJECT_NAME
 *   EPOCH TIME             → Metadata.START_TIME / STOP_TIME
 *   J2K POS/VEL            → STATE_DATA (compact array)
 *   GEOPOTENTIAL           → Perturbations.GRAVITY_MODEL/DEGREE/ORDER
 *   DRAG                   → Perturbations.ATMOSPHERIC_MODEL (ATM)
 *   LUNAR/SOLAR            → Perturbations.N_BODY_PERTURBATIONS
 *   SOLAR RAD PRESS        → Perturbations.SOLAR_RAD_PRESSURE
 *   BALLISTIC COEF         → PhysicalProperties.DRAG_COEFF_NOM
 *   SOLAR RAD PRESS COEFF  → PhysicalProperties.SOLAR_RAD_COEFF
 *   F10/AP                 → Perturbations.FIXED_F10P7 / FIXED_GEOMAG_KP
 *   VECTOR U,V,W SIGMAS    → (stored in UserDefinedParameters)
 *   COVARIANCE             → COVARIANCE_DATA (lower triangular)
 *
 * Reference: CCSDS 502.0-B-3, spacedatastandards.org/VCM
 */

#include <string>
#include <vector>
#include <cstdint>

namespace hpop {

struct VCMParseResult {
    std::vector<uint8_t> ocm_buffer;  // OCM FlatBuffer binary
    bool valid = false;
    std::string error;
    
    // Parsed metadata (for inspection)
    std::string object_name;
    std::string intl_designator;
    uint32_t norad_id = 0;
    std::string epoch_utc;
};

/**
 * Parse raw text VCM → OCM FlatBuffer binary.
 * @param vcm_text  Raw VCM text (JANAP-128 format, <> prefixed lines)
 * @return VCMParseResult with OCM FlatBuffer binary
 */
VCMParseResult parse_text_vcm(const std::string& vcm_text);

}  // namespace hpop
