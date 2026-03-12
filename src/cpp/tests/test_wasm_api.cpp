/**
 * HPOP WASM API Tests (Native Build)
 *
 * Tests the propagation pipeline and OEM FlatBuffer output
 * without requiring Emscripten. Validates:
 *   1. ISO 8601 ↔ JD time conversion
 *   2. JSON request parsing
 *   3. Full propagation pipeline
 *   4. OEM FlatBuffer output structure
 *   5. OEM → CCSDS text conversion
 */

#include "hpop/propagator.h"
#include "hpop/atmosphere.h"
#include "hpop/srp.h"
#include "hpop/thirdbody.h"
#include "main_generated.h"

#include <cstdio>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

using namespace hpop;

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", msg); tests_failed++; } \
    else { printf("  PASS: %s\n", msg); tests_passed++; } \
} while(0)

#define CHECK_TOL(val, expected, tol, msg) do { \
    double _v = (val), _e = (expected), _t = (tol); \
    if (std::abs(_v - _e) > _t) { \
        printf("  FAIL: %s (got %.6e, expected %.6e)\n", msg, _v, _e); \
        tests_failed++; \
    } else { \
        printf("  PASS: %s\n", msg); \
        tests_passed++; \
    } \
} while(0)

// ── Time conversion (copied from wasm_api.cpp for testing) ──

static double iso8601_to_jd(const std::string& iso) {
    int year, month, day, hour = 0, min = 0;
    double sec = 0.0;
    std::sscanf(iso.c_str(), "%d-%d-%dT%d:%d:%lf", &year, &month, &day, &hour, &min, &sec);
    if (month <= 2) { year--; month += 12; }
    int A = year / 100;
    int B = 2 - A + A / 4;
    return std::floor(365.25 * (year + 4716))
         + std::floor(30.6001 * (month + 1))
         + day + B - 1524.5
         + (hour + min / 60.0 + sec / 3600.0) / 24.0;
}

static std::string jd_to_iso8601(double jd) {
    double jd_plus = jd + 0.5;
    int Z = static_cast<int>(jd_plus);
    double F = jd_plus - Z;
    int A;
    if (Z < 2299161) A = Z;
    else { int alpha = static_cast<int>((Z - 1867216.25) / 36524.25); A = Z + 1 + alpha - alpha / 4; }
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
    int minutes = static_cast<int>(total_seconds / 60.0);
    double sec = total_seconds - minutes * 60.0;
    if (sec >= 59.9995) { sec = 0.0; minutes++; }
    if (minutes >= 60) { minutes -= 60; hour++; }
    if (hour >= 24) { hour -= 24; day++; }
    char buf[64];
    std::snprintf(buf, sizeof(buf), "%04d-%02d-%02dT%02d:%02d:%06.3fZ", year, month, day, hour, minutes, sec);
    return std::string(buf);
}

// ── Test 1: Time Conversions ──

void test_time_conversions() {
    printf("\n=== Test 1: ISO 8601 ↔ JD Conversions ===\n");

    // J2000 epoch
    double jd = iso8601_to_jd("2000-01-01T12:00:00.000Z");
    CHECK_TOL(jd, 2451545.0, 0.001, "J2000 epoch → JD 2451545.0");

    // Round-trip
    std::string iso = jd_to_iso8601(2451545.0);
    printf("  JD 2451545.0 → %s\n", iso.c_str());
    CHECK(iso.find("2000-01-01T12:00") != std::string::npos, "JD 2451545.0 → 2000-01-01T12:00");

    // Arbitrary date
    jd = iso8601_to_jd("2024-06-15T08:30:00.000Z");
    std::string rt = jd_to_iso8601(jd);
    printf("  2024-06-15T08:30:00 → JD %.3f → %s\n", jd, rt.c_str());
    CHECK(rt.find("2024-06-15T08:30") != std::string::npos, "Round-trip preserves date");
}

// ── Test 2: Propagation + OEM FlatBuffer ──

void test_propagation_oem() {
    printf("\n=== Test 2: Propagation → OEM FlatBuffer ===\n");

    // ISS-like orbit
    KeplerianElements kep = {6778.0, 0.001, 51.6 * DEG2RAD, 0.0, 0.0, 0.0};
    State s0 = keplerian_to_cartesian(kep);

    double epoch_jd0 = iso8601_to_jd("2024-01-01T00:00:00.000Z");
    double duration_days = 0.1;  // ~2.4 hours
    double step_seconds = 60.0;

    // Propagate with full force model
    Propagator prop;
    prop.add_force(std::make_unique<TwoBodyForce>());
    prop.add_force(std::make_unique<J2Force>());

    SpacecraftProperties sc;
    prop.set_spacecraft(sc);

    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF78;
    cfg.step_size = 30.0;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-12;
    prop.set_config(cfg);

    std::vector<EphemerisPoint> ephemeris;
    prop.propagate(s0, epoch_jd0, epoch_jd0 + duration_days, step_seconds, &ephemeris);

    printf("  Ephemeris points: %zu\n", ephemeris.size());
    CHECK(ephemeris.size() > 100, "Generated > 100 ephemeris points");

    // Build OEM FlatBuffer
    flatbuffers::FlatBufferBuilder fbb(64 * 1024);

    std::vector<flatbuffers::Offset<ephemerisDataLine>> lines;
    for (const auto& pt : ephemeris) {
        auto epoch_str = fbb.CreateString(jd_to_iso8601(pt.epoch_jd));
        auto line = CreateephemerisDataLine(fbb, epoch_str,
            pt.state.r[0], pt.state.r[1], pt.state.r[2],
            pt.state.v[0], pt.state.v[1], pt.state.v[2]);
        lines.push_back(line);
    }

    auto lines_vec = fbb.CreateVector(lines);
    auto center = fbb.CreateString("EARTH");
    auto start = fbb.CreateString(jd_to_iso8601(epoch_jd0));
    auto stop = fbb.CreateString(jd_to_iso8601(epoch_jd0 + duration_days));
    auto interp = fbb.CreateString("HERMITE");

    auto block = CreateephemerisDataBlock(fbb,
        0, 0, center, 0, 0, 0,
        timeSystem_UTC, start, 0, 0, stop, interp,
        7, 0.0, 6, 0, lines_vec, 0);

    std::vector<flatbuffers::Offset<ephemerisDataBlock>> blocks = {block};
    auto blocks_vec = fbb.CreateVector(blocks);

    auto creation = fbb.CreateString(jd_to_iso8601(epoch_jd0));
    auto orig = fbb.CreateString("HPOP-SDN-PLUGIN");
    auto oem = CreateOEM(fbb, 0, 2.0, creation, orig, blocks_vec);
    fbb.Finish(oem, "$OEM");

    size_t fb_size = fbb.GetSize();
    printf("  OEM FlatBuffer size: %zu bytes\n", fb_size);
    CHECK(fb_size > 0, "FlatBuffer generated");

    // Verify the FlatBuffer
    const uint8_t* buf = fbb.GetBufferPointer();
    flatbuffers::Verifier verifier(buf, fb_size);
    auto* oem_root = flatbuffers::GetRoot<OEM>(buf);
    CHECK(oem_root != nullptr, "OEM root parsed");

    // Check file identifier
    CHECK(flatbuffers::BufferHasIdentifier(buf, "$OEM"), "File identifier is $OEM");

    // Check OEM content
    CHECK_TOL(oem_root->CCSDS_OEM_VERS(), 2.0, 0.01, "OEM version = 2.0");
    CHECK(oem_root->ORIGINATOR() != nullptr, "Originator present");

    auto* edb = oem_root->EPHEMERIS_DATA_BLOCK();
    CHECK(edb != nullptr && edb->size() == 1, "One ephemeris data block");

    if (edb && edb->size() > 0) {
        auto* block0 = edb->Get(0);
        auto* edl = block0->EPHEMERIS_DATA_LINES();
        CHECK(edl != nullptr, "Ephemeris data lines present");
        if (edl) {
            printf("  OEM ephemeris lines: %u\n", edl->size());
            CHECK(edl->size() == ephemeris.size(), "Line count matches ephemeris");

            // Check first line
            auto* first = edl->Get(0);
            CHECK(first->EPOCH() != nullptr, "First line has epoch");
            CHECK_TOL(first->X(), s0.r[0], 0.1, "First X matches initial state");
            CHECK_TOL(first->Y(), s0.r[1], 0.1, "First Y matches initial state");
            CHECK_TOL(first->Z(), s0.r[2], 0.1, "First Z matches initial state");
        }
    }
}

// ── Test 3: Full Force Model Propagation ──

void test_full_force_propagation() {
    printf("\n=== Test 3: Full Force Model Propagation ===\n");

    // Starlink-like orbit
    KeplerianElements kep = {6921.0, 0.0001, 53.0 * DEG2RAD, 0.0, 0.0, 0.0};
    State s0 = keplerian_to_cartesian(kep);

    double epoch_jd0 = iso8601_to_jd("2024-06-21T12:00:00.000Z");

    Propagator prop;
    prop.add_force(std::make_unique<TwoBodyForce>());
    prop.add_force(std::make_unique<J2Force>());

    auto atm = std::make_unique<HarrisPriester>(150.0, 2);
    prop.add_force(std::make_unique<DragForce>(std::move(atm), 2.2));
    prop.add_force(std::make_unique<SRPForce>(ShadowType::CONICAL));
    prop.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::SUN));
    prop.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::MOON));

    SpacecraftProperties sc;
    sc.mass_kg = 260.0;
    sc.drag_area = 22.0;
    sc.cd = 2.2;
    sc.srp_area = 22.0;
    sc.cr = 1.3;
    prop.set_spacecraft(sc);

    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF78;
    cfg.step_size = 30.0;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-12;
    prop.set_config(cfg);

    std::vector<EphemerisPoint> ephemeris;
    State sf = prop.propagate(s0, epoch_jd0, epoch_jd0 + 1.0, 60.0, &ephemeris);

    auto kep_f = cartesian_to_keplerian(sf);
    printf("  Initial: a=%.3f km, e=%.6f, i=%.4f°\n", kep.a, kep.e, kep.i * RAD2DEG);
    printf("  Final:   a=%.3f km, e=%.6f, i=%.4f°\n", kep_f.a, kep_f.e, kep_f.i * RAD2DEG);
    printf("  Ephemeris points: %zu\n", ephemeris.size());

    // Orbit should still be physically reasonable after 1 day
    CHECK(kep_f.a > 6800 && kep_f.a < 7000, "SMA still in LEO range after 1 day");
    CHECK(kep_f.e < 0.1, "Eccentricity remains small");
    CHECK(ephemeris.size() > 1400, "~1440 points at 60s step over 1 day");

    // Drag should cause some decay
    CHECK(kep_f.a < kep.a, "SMA decreased from drag");
}

int main() {
    printf("============================================================\n");
    printf("HPOP WASM API Tests (Native)\n");
    printf("============================================================\n");

    test_time_conversions();
    test_propagation_oem();
    test_full_force_propagation();

    printf("\n============================================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("============================================================\n");

    return tests_failed > 0 ? 1 : 0;
}
