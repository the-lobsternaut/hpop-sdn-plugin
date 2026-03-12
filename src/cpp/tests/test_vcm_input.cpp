/**
 * HPOP VCM Input Parser Tests
 *
 * Test strategy:
 *   1. Parse VCM FlatBuffer with state vector
 *   2. Parse VCM FlatBuffer with Keplerian elements (mean anomaly → true anomaly)
 *   3. Parse VCM with spacecraft properties (mass, drag, SRP)
 *   4. Parse VCM atmospheric model enums
 *   5. Parse JSON with simplified format
 *   6. Parse JSON with VCM-schema fields
 *   7. Auto-detect format (FlatBuffer vs JSON)
 *   8. Full pipeline: VCM → parse → propagate → OEM
 */

#include "hpop/vcm_input.h"
#include "hpop/propagator.h"
#include "hpop/atmosphere.h"
#include "hpop/srp.h"
#include "hpop/thirdbody.h"
#include "sds/vcm_generated.h"

#include <cstdio>
#include <cmath>
#include <memory>

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

// ── Helper: Build a VCM FlatBuffer ──

static std::string build_test_vcm(
    double x, double y, double z,
    double vx, double vy, double vz,
    const char* epoch = "2024-01-01T00:00:00.000Z",
    double mass = 260.0,
    double drag_area = 22.0,
    double drag_coeff = 2.2,
    double srp_area = 22.0,
    double srp_coeff = 1.3,
    atmosphericModel atm_model = atmosphericModel_HASDM,
    bool srp_on = true)
{
    flatbuffers::FlatBufferBuilder fbb(4096);

    // State vector
    auto sv = CreateVCMStateVectorDirect(fbb, epoch, x, y, z, vx, vy, vz);

    // Atmospheric model data
    auto atm = CreateVCMAtmosphericModelData(fbb,
        atm_model,
        geopotentialModel_EGM96,
        perturbationStatus_ON,   // lunar-solar
        lunarPerturbationModel_DE430,
        solarPerturbationModel_DE430,
        srp_on ? perturbationStatus_ON : perturbationStatus_OFF,
        solarRadiationPressureModel_SPHERICAL_MODEL,
        resonanceModel_NONE
    );

    // Object name
    auto obj_name = fbb.CreateString("ISS (ZARYA)");
    auto obj_id = fbb.CreateString("1998-067A");
    auto originator_str = fbb.CreateString("18 SPCS");
    auto center = fbb.CreateString("EARTH");
    auto frame = fbb.CreateString("TEME");
    auto time_sys = fbb.CreateString("UTC");

    auto vcm = CreateVCM(fbb,
        2.0,              // CCSDS_OMM_VERS
        0,                // CREATION_DATE
        originator_str,
        obj_name,
        obj_id,
        center,
        frame,
        time_sys,
        sv,               // STATE_VECTOR
        0,                // KEPLERIAN_ELEMENTS
        0,                // EQUINOCTIAL_ELEMENTS
        398600.4418,      // GM
        atm,              // ATMOSPHERIC_MODEL_DATA
        0,                // PROPAGATOR_SETTINGS
        0,                // UVW_SIGMAS
        mass,
        srp_area,
        srp_coeff,
        drag_area,
        drag_coeff,
        srp_on ? perturbationStatus_ON : perturbationStatus_OFF
    );

    fbb.Finish(vcm, "$VCM");
    const uint8_t* buf = fbb.GetBufferPointer();
    return std::string(reinterpret_cast<const char*>(buf), fbb.GetSize());
}

// ── Test 1: Parse VCM FlatBuffer with State Vector ──

void test_vcm_flatbuffer_state() {
    printf("\n=== Test 1: Parse VCM FlatBuffer — State Vector ===\n");

    auto vcm_buf = build_test_vcm(6778.0, 0, 0, 0, 7.669, 0);
    auto input = parse_vcm_flatbuffer(
        reinterpret_cast<const uint8_t*>(vcm_buf.data()), vcm_buf.size());

    CHECK(input.valid, "VCM parsed successfully");
    CHECK_TOL(input.initial_state.r[0], 6778.0, 0.001, "X = 6778 km");
    CHECK_TOL(input.initial_state.v[1], 7.669, 0.001, "VY = 7.669 km/s");
    CHECK(input.epoch_jd > 2460000, "Epoch JD is reasonable");
    printf("  Object: %s (%s)\n", input.object_name.c_str(), input.object_id.c_str());
    CHECK(input.object_name == "ISS (ZARYA)", "Object name parsed");
    CHECK(input.object_id == "1998-067A", "Object ID parsed");
}

// ── Test 2: Parse VCM with Keplerian Elements ──

void test_vcm_keplerian() {
    printf("\n=== Test 2: Parse VCM FlatBuffer — Keplerian Elements ===\n");

    flatbuffers::FlatBufferBuilder fbb(4096);

    // ISS orbit: a=6778, e=0.001, i=51.6°
    auto kep = CreatekeplerianElements(fbb,
        6778.0,    // SMA (km)
        0.001,     // ecc
        51.6,      // inc (deg)
        0.0,       // RAAN (deg)
        0.0,       // AoP (deg)
        anomalyType_TRUE_ANOMALY,
        0.0        // ν (deg)
    );

    auto sv = CreateVCMStateVectorDirect(fbb, "2024-01-01T00:00:00.000Z", 0, 0, 0, 0, 0, 0);
    auto obj_name = fbb.CreateString("TEST_SAT");
    auto obj_id = fbb.CreateString("2024-001A");
    auto center = fbb.CreateString("EARTH");
    auto frame = fbb.CreateString("TEME");
    auto time_sys = fbb.CreateString("UTC");

    auto vcm = CreateVCM(fbb,
        2.0, 0, 0, obj_name, obj_id, center, frame, time_sys,
        sv, kep, 0, 398600.4418);
    fbb.Finish(vcm, "$VCM");

    auto input = parse_vcm_flatbuffer(fbb.GetBufferPointer(), fbb.GetSize());

    CHECK(input.valid, "VCM with Keplerian parsed");
    CHECK(input.has_keplerian, "Keplerian elements present");
    CHECK_TOL(input.keplerian.a, 6778.0, 0.01, "SMA = 6778 km");

    // State should have been converted from Keplerian
    double rmag = input.initial_state.rmag();
    printf("  Converted r = [%.1f, %.1f, %.1f] km, |r| = %.1f km\n",
           input.initial_state.r[0], input.initial_state.r[1], input.initial_state.r[2], rmag);
    CHECK_TOL(rmag, 6778.0 * (1 - 0.001), 10.0, "Position magnitude ≈ a(1-e)");
}

// ── Test 3: Spacecraft Properties ──

void test_vcm_spacecraft() {
    printf("\n=== Test 3: Parse VCM — Spacecraft Properties ===\n");

    auto vcm_buf = build_test_vcm(6778.0, 0, 0, 0, 7.669, 0,
        "2024-01-01T00:00:00.000Z", 420000.0, 1600.0, 2.5, 2500.0, 1.8);
    auto input = parse_vcm_flatbuffer(
        reinterpret_cast<const uint8_t*>(vcm_buf.data()), vcm_buf.size());

    CHECK_TOL(input.spacecraft.mass_kg, 420000.0, 1.0, "Mass = 420000 kg (ISS)");
    CHECK_TOL(input.spacecraft.drag_area, 1600.0, 1.0, "Drag area = 1600 m²");
    CHECK_TOL(input.spacecraft.cd, 2.5, 0.01, "Cd = 2.5");
    CHECK_TOL(input.spacecraft.srp_area, 2500.0, 1.0, "SRP area = 2500 m²");
    CHECK_TOL(input.spacecraft.cr, 1.8, 0.01, "Cr = 1.8");
}

// ── Test 4: Atmospheric Model Mapping ──

void test_vcm_atmosphere() {
    printf("\n=== Test 4: Parse VCM — Atmospheric Model Mapping ===\n");

    // HASDM
    auto buf1 = build_test_vcm(6778.0, 0, 0, 0, 7.669, 0,
        "2024-01-01T00:00:00Z", 260, 22, 2.2, 22, 1.3, atmosphericModel_HASDM);
    auto in1 = parse_vcm_flatbuffer(reinterpret_cast<const uint8_t*>(buf1.data()), buf1.size());
    printf("  HASDM → %s\n", in1.atmosphere_model.c_str());
    CHECK(in1.atmosphere_model == "hasdm", "HASDM mapped correctly");
    CHECK(in1.use_drag, "Drag enabled for HASDM");

    // Jacchia 70
    auto buf2 = build_test_vcm(6778.0, 0, 0, 0, 7.669, 0,
        "2024-01-01T00:00:00Z", 260, 22, 2.2, 22, 1.3, atmosphericModel_JACCHIA_70);
    auto in2 = parse_vcm_flatbuffer(reinterpret_cast<const uint8_t*>(buf2.data()), buf2.size());
    CHECK(in2.atmosphere_model == "jacchia-70", "Jacchia 70 mapped correctly");

    // NONE → drag disabled
    auto buf3 = build_test_vcm(6778.0, 0, 0, 0, 7.669, 0,
        "2024-01-01T00:00:00Z", 260, 22, 2.2, 22, 1.3, atmosphericModel_NONE);
    auto in3 = parse_vcm_flatbuffer(reinterpret_cast<const uint8_t*>(buf3.data()), buf3.size());
    CHECK(!in3.use_drag, "Drag disabled when atmosphere = NONE");

    // SRP off
    auto buf4 = build_test_vcm(6778.0, 0, 0, 0, 7.669, 0,
        "2024-01-01T00:00:00Z", 260, 22, 2.2, 22, 1.3, atmosphericModel_HASDM, false);
    auto in4 = parse_vcm_flatbuffer(reinterpret_cast<const uint8_t*>(buf4.data()), buf4.size());
    CHECK(!in4.use_srp, "SRP disabled when SRP = OFF");
}

// ── Test 5: JSON Parsing (Simplified) ──

void test_json_simple() {
    printf("\n=== Test 5: Parse JSON — Simplified Format ===\n");

    std::string json = R"({
        "epoch": "2024-06-21T12:00:00.000Z",
        "x": 6921.0, "y": 0, "z": 0,
        "vx": 0, "vy": 7.59, "vz": 0,
        "mass_kg": 260, "drag_area": 22, "cd": 2.2,
        "srp_area": 22, "cr": 1.3,
        "duration_days": 0.5,
        "step_seconds": 30.0,
        "f107": 200.0
    })";

    auto input = parse_vcm_json(json);

    CHECK(input.valid, "JSON parsed successfully");
    CHECK_TOL(input.initial_state.r[0], 6921.0, 0.01, "X = 6921 km");
    CHECK_TOL(input.spacecraft.mass_kg, 260.0, 0.01, "Mass = 260 kg");
    CHECK_TOL(input.duration_days, 0.5, 0.001, "Duration = 0.5 days");
    CHECK_TOL(input.step_seconds, 30.0, 0.01, "Step = 30 s");
    CHECK_TOL(input.f107, 200.0, 0.01, "F10.7 = 200");
}

// ── Test 6: JSON Parsing (VCM Schema) ──

void test_json_vcm_schema() {
    printf("\n=== Test 6: Parse JSON — VCM Schema Fields ===\n");

    std::string json = R"({
        "OBJECT_NAME": "STARLINK-1234",
        "OBJECT_ID": "2020-001A",
        "EPOCH": "2024-01-15T06:00:00.000Z",
        "X": 6778.0, "Y": 100.0, "Z": 200.0,
        "X_DOT": -0.1, "Y_DOT": 7.6, "Z_DOT": 0.5,
        "MASS": 260.0,
        "DRAG_AREA": 22.0,
        "DRAG_COEFF": 2.2,
        "SOLAR_RAD_AREA": 22.0,
        "SOLAR_RAD_COEFF": 1.3,
        "NORAD_CAT_ID": 44235,
        "duration_days": 1.0,
        "TIME_STEP": 120.0
    })";

    auto input = parse_vcm_json(json);

    CHECK(input.valid, "VCM-schema JSON parsed");
    CHECK(input.object_name == "STARLINK-1234", "Object name from VCM field");
    CHECK(input.norad_cat_id == 44235, "NORAD ID parsed");
    CHECK_TOL(input.step_seconds, 120.0, 0.01, "Step from TIME_STEP = 120 s");
    CHECK_TOL(input.initial_state.r[1], 100.0, 0.01, "Y from VCM-style field");
}

// ── Test 7: Auto-detect Format ──

void test_auto_detect() {
    printf("\n=== Test 7: Auto-detect Format ===\n");

    // FlatBuffer
    auto vcm_buf = build_test_vcm(6778.0, 0, 0, 0, 7.669, 0);
    auto input_fb = parse_vcm(
        reinterpret_cast<const uint8_t*>(vcm_buf.data()), vcm_buf.size());
    CHECK(input_fb.valid, "FlatBuffer auto-detected");
    CHECK_TOL(input_fb.initial_state.r[0], 6778.0, 0.01, "FlatBuffer state correct");

    // JSON
    std::string json = R"({"epoch":"2024-01-01T00:00:00Z","x":6778,"y":0,"z":0,"vx":0,"vy":7.669,"vz":0})";
    auto input_json = parse_vcm(
        reinterpret_cast<const uint8_t*>(json.data()), json.size());
    CHECK(input_json.valid, "JSON auto-detected");
    CHECK_TOL(input_json.initial_state.r[0], 6778.0, 0.01, "JSON state correct");
}

// ── Test 8: Full Pipeline (VCM → Parse → Propagate) ──

void test_full_pipeline() {
    printf("\n=== Test 8: Full Pipeline — VCM → Propagate ===\n");

    auto vcm_buf = build_test_vcm(
        -2517.232, -663.786, 6691.339,    // ISS-like position
        -1.768, 7.315, 1.439,             // ISS-like velocity
        "2024-01-01T00:00:00.000Z",
        420000.0, 1600.0, 2.2, 2500.0, 1.3);

    auto input = parse_vcm(
        reinterpret_cast<const uint8_t*>(vcm_buf.data()), vcm_buf.size());
    CHECK(input.valid, "VCM parsed for pipeline");

    // Build propagator from parsed input
    Propagator prop;
    if (input.use_twobody) prop.add_force(std::make_unique<TwoBodyForce>());
    if (input.use_j2) prop.add_force(std::make_unique<J2Force>());
    if (input.use_drag) {
        auto atm = std::make_unique<HarrisPriester>(input.f107, 2);
        prop.add_force(std::make_unique<DragForce>(std::move(atm)));
    }
    if (input.use_srp) prop.add_force(std::make_unique<SRPForce>(ShadowType::CONICAL));
    if (input.use_sun) prop.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::SUN));
    if (input.use_moon) prop.add_force(std::make_unique<ThirdBodyForce>(ThirdBody::MOON));

    prop.set_spacecraft(input.spacecraft);

    PropagatorConfig cfg;
    cfg.integrator = IntegratorType::RKF78;
    cfg.step_size = 30.0;
    cfg.abs_tol = 1e-10;
    cfg.rel_tol = 1e-12;
    prop.set_config(cfg);

    std::vector<EphemerisPoint> ephemeris;
    double duration = 0.1;  // 2.4 hours
    State sf = prop.propagate(input.initial_state, input.epoch_jd,
                               input.epoch_jd + duration, 60.0, &ephemeris);

    printf("  Propagated %zu points over %.1f hours\n", ephemeris.size(), duration * 24);
    printf("  Final: r=[%.1f, %.1f, %.1f] km\n", sf.r[0], sf.r[1], sf.r[2]);

    CHECK(ephemeris.size() > 100, "Generated ephemeris from VCM input");
    CHECK(sf.rmag() > 6300 && sf.rmag() < 8500, "Final state physically reasonable (LEO range)");
}

int main() {
    printf("============================================================\n");
    printf("HPOP VCM Input Parser Tests\n");
    printf("============================================================\n");

    test_vcm_flatbuffer_state();
    test_vcm_keplerian();
    test_vcm_spacecraft();
    test_vcm_atmosphere();
    test_json_simple();
    test_json_vcm_schema();
    test_auto_detect();
    test_full_pipeline();

    printf("\n============================================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("============================================================\n");

    return tests_failed > 0 ? 1 : 0;
}
