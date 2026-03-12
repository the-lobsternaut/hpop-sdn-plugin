#pragma once
/**
 * OCM FlatBuffer → Raw Text VCM Writer
 *
 * Converts spacedatastandards.org OCM ($OCM) FlatBuffer binary back to
 * 18 SPCS fixed-format ASCII VCM (V2.0).
 *
 * This is the inverse of vcm_parser.h.
 */

#include <string>
#include <vector>
#include <cstdint>

namespace hpop {

/**
 * Convert OCM FlatBuffer binary → raw text VCM string.
 * @param data  OCM FlatBuffer binary
 * @param size  Size of the buffer
 * @return Raw text VCM string (empty on error)
 */
std::string ocm_to_text_vcm(const uint8_t* data, size_t size);

}  // namespace hpop
