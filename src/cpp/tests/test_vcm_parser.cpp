/**
 * Raw Text VCM → OCM FlatBuffer Parser Tests
 *
 * Tests parsing of 18 SPCS fixed-format ASCII VCM messages
 * into spacedatastandards.org OCM FlatBuffer binary.
 *
 * Uses the sample VCM from the spacedatastandards.org survey.
 */

#include "hpop/vcm_parser.h"
#include "hpop/vcm_writer.h"
#include "sds/ocm_generated.h"

#include <cstdio>
#include <cmath>
#include <cstring>

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
    if (!(cond)) { printf("  FAIL: %s\n", msg); tests_failed++; } \
    else { printf("  PASS: %s\n", msg); tests_passed++; } \
} while(0)

#define CHECK_TOL(val, expected, tol, msg) do { \
    double _v = (val), _e = (expected), _t = (tol); \
    if (std::abs(_v - _e) > _t) { \
        printf("  FAIL: %s (got %.10e, expected %.10e, diff %.2e)\n", msg, _v, _e, std::abs(_v - _e)); \
        tests_failed++; \
    } else { \
        printf("  PASS: %s\n", msg); \
        tests_passed++; \
    } \
} while(0)

// ── Sample VCM from spacedatastandards.org survey ──

static const char* SAMPLE_VCM = R"(<> SP VECTOR/COVARIANCE MESSAGE - V2.0
<>
<> MESSAGE TIME (UTC): 2023 007 (07 JAN) 22:45:06.000  CENTER: FAKE_CENTER
<> SATELLITE NUMBER: 00000                    INT. DES.: 0000-000A
<> COMMON NAME:
<> EPOCH TIME (UTC): 2023 007 (07 JAN) 21:39: 8.414  EPOCH REV: 37693
<> J2K POS (KM):       369.89521989    5103.23778139    4467.41312115
<> J2K VEL (KM/S):  -6.491206849227  -2.407500055425   3.279477293781
<> ECI POS (KM):       333.71334681    5104.90725904    4468.35490991
<> ECI VEL (KM/S):  -6.485984386597  -2.441004400229   3.265011217687
<> EFG POS (KM):      4957.72595592    1261.90176605    4468.35490991
<> EFG VEL (KM/S):  -4.235735566625   5.051150329994   3.265011217687
<> GEOPOTENTIAL: EGM-96 70Z,70T  DRAG: JAC70/MSIS90  LUNAR/SOLAR:  ON
<> SOLAR RAD PRESS: OFF  SOLID EARTH TIDES:  ON  IN-TRACK THRUST: OFF
<> BALLISTIC COEF (M2/KG):  0.826455E-02 BDOT (M2/KG-S): 0.000000E+00
<> SOLAR RAD PRESS COEFF (M2/KG):  0.000000E+00  EDR(W/KG):  0.39E-02
<> THRUST ACCEL (M/S2):  0.000000E+00  C.M. OFFSET (M):  0.000000E+00
<> SOLAR FLUX: F10: 153  AVERAGE F10: 137  AVERAGE AP:  10.0
<> TAI-UTC (S): 37  UT1-UTC (S): -0.01719  UT1 RATE (MS/DAY):  0.384
<> POLAR MOT X,Y (ARCSEC):  0.0505  0.2109 IAU 1980 NUTAT:   4 TERMS
<> TIME CONST LEAP SECOND TIME (UTC): 2049 365 (31 DEC) 23:59:59.999
<> INTEGRATOR MODE: ASW          COORD SYS: J2000  PARTIALS: FAST NUM
<> STEP MODE: AUTO  FIXED STEP: OFF  STEP SIZE SELECTION: MANUAL
<> INITIAL STEP SIZE (S):   20.000  ERROR CONTROL: 0.100E-13
<> VECTOR U,V,W SIGMAS (KM):           0.0084     0.0402     0.0074
<> VECTOR UD,VD,WD SIGMAS (KM/S):      0.0000     0.0000     0.0000
<> COVARIANCE MATRIX (EQUINOCTIAL ELS): ( 7x 7) WTD RMS:  0.11498E+01
<>  0.20458E-11  0.10230E-11  0.13484E-11  0.13797E-11 -0.47633E-12
<>  0.13642E-10  0.13871E-12  0.91700E-13  0.70998E-12  0.93201E-13
<>  0.50454E-12  0.26963E-12 -0.10331E-12  0.40296E-14  0.39682E-12
<> -0.38952E-12 -0.24005E-12  0.73298E-14 -0.10121E-13 -0.17845E-12
<>  0.40019E-12  0.91367E-08  0.53586E-08  0.36570E-07  0.50105E-08
<>  0.93145E-09 -0.57204E-09  0.13665E-02)";

// VCM with named satellite
static const char* ISS_VCM = R"(<> SP VECTOR/COVARIANCE MESSAGE - V2.0
<> REAL
<> MESSAGE TIME (UTC): 2024 180 (28 JUN) 15:30:00.000  CENTER: 18SPCS
<> SATELLITE NUMBER: 25544                    INT. DES.: 1998-067A
<> COMMON NAME: ISS (ZARYA)
<> EPOCH TIME (UTC): 2024 180 (28 JUN) 14:22:33.100  EPOCH REV: 46100
<> J2K POS (KM):     -2517.23200000   -663.78600000    6691.33900000
<> J2K VEL (KM/S):  -1.768000000000   7.315000000000   1.439000000000
<> ECI POS (KM):     -2500.00000000   -700.00000000    6680.00000000
<> ECI VEL (KM/S):  -1.770000000000   7.310000000000   1.440000000000
<> EFG POS (KM):      5000.00000000   1200.00000000    6680.00000000
<> EFG VEL (KM/S):  -4.200000000000   5.000000000000   1.440000000000
<> GEOPOTENTIAL: EGM-96 36Z,36T  DRAG: JB2008      LUNAR/SOLAR:  ON
<> SOLAR RAD PRESS:  ON  SOLID EARTH TIDES:  ON  IN-TRACK THRUST: OFF
<> BALLISTIC COEF (M2/KG):  3.810000E-03 BDOT (M2/KG-S): 1.200000E-06
<> SOLAR RAD PRESS COEFF (M2/KG):  5.952381E-03  EDR(W/KG):  0.00E+00
<> THRUST ACCEL (M/S2):  0.000000E+00  C.M. OFFSET (M):  0.000000E+00
<> SOLAR FLUX: F10: 180  AVERAGE F10: 165  AVERAGE AP:  15.0
<> TAI-UTC (S): 37  UT1-UTC (S):  0.02500  UT1 RATE (MS/DAY):  0.100
<> POLAR MOT X,Y (ARCSEC):  0.1000  0.3000 IAU 1980 NUTAT: 106 TERMS
<> TIME CONST LEAP SECOND TIME (UTC): 2049 365 (31 DEC) 23:59:59.999
<> INTEGRATOR MODE: ASW          COORD SYS: J2000  PARTIALS: ANALYTIC
<> STEP MODE: AUTO  FIXED STEP: OFF  STEP SIZE SELECTION: AUTO
<> INITIAL STEP SIZE (S):   30.000  ERROR CONTROL: 0.100E-13
<> VECTOR U,V,W SIGMAS (KM):           0.0150     0.0800     0.0120
<> VECTOR UD,VD,WD SIGMAS (KM/S):      0.0001     0.0002     0.0001
<> COVARIANCE MATRIX (EQUINOCTIAL ELS): ( 9x 9) WTD RMS:  0.95000E+00
<>  1.00000E-10  2.00000E-10  3.00000E-10  4.00000E-10  5.00000E-10
<>  6.00000E-10  7.00000E-10  8.00000E-10  9.00000E-10  1.10000E-09
<>  1.20000E-09  1.30000E-09  1.40000E-09  1.50000E-09  1.60000E-09
<>  1.70000E-09  1.80000E-09  1.90000E-09  2.00000E-09  2.10000E-09
<>  2.20000E-09  2.30000E-09  2.40000E-09  2.50000E-09  2.60000E-09
<>  2.70000E-09  2.80000E-09  2.90000E-09  3.00000E-09  3.10000E-09
<>  3.20000E-09  3.30000E-09  3.40000E-09  3.50000E-09  3.60000E-09
<>  3.70000E-09  3.80000E-09  3.90000E-09  4.00000E-09  4.10000E-09
<>  4.20000E-09  4.30000E-09  4.40000E-09  4.50000E-09  4.50000E-09)";

// ── Test 1: Parse sample VCM ──

void test_parse_sample() {
    printf("\n=== Test 1: Parse Sample VCM ===\n");

    auto result = hpop::parse_text_vcm(SAMPLE_VCM);
    CHECK(result.valid, "VCM parsed successfully");
    CHECK(!result.ocm_buffer.empty(), "OCM buffer generated");

    // Verify file identifier
    CHECK(result.ocm_buffer.size() >= 8, "Buffer large enough for FlatBuffer");
    const char* fid = reinterpret_cast<const char*>(result.ocm_buffer.data() + 4);
    CHECK(std::strncmp(fid, "$OCM", 4) == 0, "File identifier is $OCM");

    printf("  NORAD: %u, Designator: %s\n", result.norad_id, result.intl_designator.c_str());
    CHECK(result.norad_id == 0, "NORAD ID = 00000");
    CHECK(result.intl_designator == "0000-000A", "International designator");
}

// ── Test 2: Read back OCM fields ──

void test_ocm_fields() {
    printf("\n=== Test 2: Read Back OCM Fields ===\n");

    auto result = hpop::parse_text_vcm(SAMPLE_VCM);
    CHECK(result.valid, "Parsed OK");

    auto* ocm = flatbuffers::GetRoot<OCM>(result.ocm_buffer.data());
    CHECK(ocm != nullptr, "OCM root accessible");

    // Header
    auto* hdr = ocm->HEADER();
    CHECK(hdr != nullptr, "Header present");
    if (hdr) {
        CHECK(hdr->ORIGINATOR() && std::string(hdr->ORIGINATOR()->c_str()) == "FAKE_CENTER",
              "Originator = FAKE_CENTER");
        printf("  Creation date: %s\n", hdr->CREATION_DATE() ? hdr->CREATION_DATE()->c_str() : "null");
    }

    // Metadata
    auto* meta = ocm->METADATA();
    CHECK(meta != nullptr, "Metadata present");
    if (meta) {
        printf("  Epoch: %s\n", meta->EPOCH_TZERO() ? meta->EPOCH_TZERO()->c_str() : "null");
        CHECK(meta->INTERNATIONAL_DESIGNATOR() &&
              std::string(meta->INTERNATIONAL_DESIGNATOR()->c_str()) == "0000-000A",
              "International designator = 0000-000A");
    }
}

// ── Test 3: State vector (J2K) ──

void test_state_vector() {
    printf("\n=== Test 3: J2K State Vector in OCM ===\n");

    auto result = hpop::parse_text_vcm(SAMPLE_VCM);
    auto* ocm = flatbuffers::GetRoot<OCM>(result.ocm_buffer.data());
    CHECK(ocm != nullptr, "OCM accessible");

    auto* state_data = ocm->STATE_DATA();
    CHECK(state_data != nullptr, "STATE_DATA present");
    CHECK(state_data && state_data->size() == 6, "6 components (pos + vel)");

    if (state_data && state_data->size() >= 6) {
        CHECK_TOL(state_data->Get(0), 369.89521989, 1e-6, "X = 369.895 km");
        CHECK_TOL(state_data->Get(1), 5103.23778139, 1e-6, "Y = 5103.238 km");
        CHECK_TOL(state_data->Get(2), 4467.41312115, 1e-6, "Z = 4467.413 km");
        CHECK_TOL(state_data->Get(3), -6.491206849227, 1e-10, "Vx = -6.4912 km/s");
        CHECK_TOL(state_data->Get(4), -2.407500055425, 1e-10, "Vy = -2.4075 km/s");
        CHECK_TOL(state_data->Get(5), 3.279477293781, 1e-10, "Vz = 3.2795 km/s");
    }
}

// ── Test 4: Perturbation models ──

void test_perturbations() {
    printf("\n=== Test 4: Perturbation Models ===\n");

    auto result = hpop::parse_text_vcm(SAMPLE_VCM);
    auto* ocm = flatbuffers::GetRoot<OCM>(result.ocm_buffer.data());
    auto* pert = ocm->PERTURBATIONS();
    CHECK(pert != nullptr, "Perturbations present");

    if (pert) {
        // Gravity model
        CHECK(pert->GRAVITY_MODEL() &&
              std::string(pert->GRAVITY_MODEL()->c_str()).find("EGM-96") != std::string::npos,
              "Gravity model = EGM-96");
        CHECK(pert->GRAVITY_DEGREE() == 70, "Gravity degree = 70");
        CHECK(pert->GRAVITY_ORDER() == 70, "Gravity order = 70");

        // Atmospheric model
        auto* atm = pert->ATMOSPHERIC_MODEL();
        CHECK(atm != nullptr, "ATM present");
        if (atm) {
            printf("  ATM model: %d, year: %d\n", atm->MODEL(), atm->YEAR());
            CHECK(atm->MODEL() == AtmosphericModelFamily_JXX, "ATM = JXX (Jacchia)");
            CHECK(atm->YEAR() == 1970, "ATM year = 1970");
        }

        // N-body
        auto* nbody = pert->N_BODY_PERTURBATIONS();
        CHECK(nbody != nullptr && nbody->size() == 2, "Lunar/Solar ON → 2 bodies");

        // SRP
        CHECK(pert->SOLAR_RAD_PRESSURE() &&
              std::string(pert->SOLAR_RAD_PRESSURE()->c_str()) == "OFF",
              "SRP = OFF");

        // F10.7
        CHECK_TOL(pert->FIXED_F10P7(), 153.0, 1.0, "F10 = 153");
        CHECK_TOL(pert->FIXED_F10P7_MEAN(), 137.0, 1.0, "Avg F10 = 137");
        CHECK_TOL(pert->FIXED_GEOMAG_KP(), 10.0, 0.5, "Avg AP = 10.0");
    }
}

// ── Test 5: Physical properties ──

void test_physical_properties() {
    printf("\n=== Test 5: Physical Properties ===\n");

    auto result = hpop::parse_text_vcm(SAMPLE_VCM);
    auto* ocm = flatbuffers::GetRoot<OCM>(result.ocm_buffer.data());
    auto* phys = ocm->PHYSICAL_PROPERTIES();
    CHECK(phys != nullptr, "PhysicalProperties present");

    if (phys) {
        CHECK_TOL(phys->DRAG_COEFF_NOM(), 0.826455e-02, 1e-8, "Ballistic coeff = 0.826455E-02");
        CHECK_TOL(phys->SOLAR_RAD_COEFF(), 0.0, 1e-10, "SRP coeff = 0");
    }
}

// ── Test 6: Covariance matrix ──

void test_covariance() {
    printf("\n=== Test 6: Covariance Matrix ===\n");

    auto result = hpop::parse_text_vcm(SAMPLE_VCM);
    auto* ocm = flatbuffers::GetRoot<OCM>(result.ocm_buffer.data());

    auto* cov = ocm->COVARIANCE_DATA();
    CHECK(cov != nullptr, "Covariance data present");
    if (cov) {
        printf("  Covariance elements: %u\n", cov->size());
        // 7x7 lower triangular = 7*8/2 = 28 elements
        CHECK(cov->size() == 28, "7x7 lower triangular = 28 elements");

        if (cov->size() >= 28) {
            CHECK_TOL(cov->Get(0), 0.20458e-11, 1e-16, "Cov[0] = 0.20458E-11");
            CHECK_TOL(cov->Get(1), 0.10230e-11, 1e-16, "Cov[1] = 0.10230E-11");
            CHECK_TOL(cov->Get(27), 0.13665e-02, 1e-7, "Cov[27] = 0.13665E-02");
        }
    }
}

// ── Test 7: ISS VCM with named satellite and JB2008 ──

void test_iss_vcm() {
    printf("\n=== Test 7: ISS VCM (JB2008, SRP ON) ===\n");

    auto result = hpop::parse_text_vcm(ISS_VCM);
    CHECK(result.valid, "ISS VCM parsed");
    CHECK(result.norad_id == 25544, "NORAD ID = 25544");
    CHECK(result.object_name == "ISS (ZARYA)", "Name = ISS (ZARYA)");
    CHECK(result.intl_designator == "1998-067A", "Designator = 1998-067A");

    auto* ocm = flatbuffers::GetRoot<OCM>(result.ocm_buffer.data());
    auto* pert = ocm->PERTURBATIONS();
    CHECK(pert != nullptr, "Perturbations present");

    if (pert) {
        auto* atm = pert->ATMOSPHERIC_MODEL();
        CHECK(atm != nullptr && atm->MODEL() == AtmosphericModelFamily_JB08, "ATM = JB2008");

        CHECK(pert->SOLAR_RAD_PRESSURE() &&
              std::string(pert->SOLAR_RAD_PRESSURE()->c_str()) == "ON",
              "SRP = ON");

        CHECK_TOL(pert->FIXED_F10P7(), 180.0, 1.0, "F10 = 180");
        CHECK_TOL(pert->GRAVITY_DEGREE(), 36, 0, "Gravity degree = 36");
        CHECK_TOL(pert->GRAVITY_ORDER(), 36, 0, "Gravity order = 36");
    }

    // State vector
    auto* state = ocm->STATE_DATA();
    if (state && state->size() >= 6) {
        CHECK_TOL(state->Get(0), -2517.232, 0.001, "ISS X position");
        CHECK_TOL(state->Get(3), -1.768, 0.001, "ISS Vx velocity");
    }

    // Physical properties
    auto* phys = ocm->PHYSICAL_PROPERTIES();
    if (phys) {
        CHECK_TOL(phys->DRAG_COEFF_NOM(), 3.81e-03, 1e-6, "ISS ballistic coeff");
        CHECK_TOL(phys->SOLAR_RAD_COEFF(), 5.952381e-03, 1e-8, "ISS SRP coeff");
    }

    // Covariance (9x9 = 45 elements)
    auto* cov = ocm->COVARIANCE_DATA();
    CHECK(cov != nullptr, "ISS covariance present");
    if (cov) {
        printf("  ISS covariance elements: %u\n", cov->size());
        CHECK(cov->size() == 45, "9x9 lower triangular = 45 elements");
    }
}

// ── Test 8: User-defined parameters ──

void test_user_params() {
    printf("\n=== Test 8: User-Defined Parameters ===\n");

    auto result = hpop::parse_text_vcm(SAMPLE_VCM);
    auto* ocm = flatbuffers::GetRoot<OCM>(result.ocm_buffer.data());

    auto* params = ocm->USER_DEFINED_PARAMETERS();
    CHECK(params != nullptr, "User-defined parameters present");
    if (params) {
        printf("  Total user params: %u\n", params->size());
        CHECK(params->size() > 0, "At least one user parameter");

        // Find specific params
        bool found_eci = false, found_efg = false, found_uvw = false;
        bool found_integrator = false, found_norad = false;
        for (unsigned i = 0; i < params->size(); i++) {
            auto* p = params->Get(i);
            std::string name = p->PARAM_NAME()->str();
            if (name == "ECI_POS_KM") found_eci = true;
            if (name == "EFG_POS_KM") found_efg = true;
            if (name == "UVW_POS_SIGMA_KM") found_uvw = true;
            if (name == "INTEGRATOR_MODE") found_integrator = true;
            if (name == "NORAD_CAT_ID") found_norad = true;
        }
        CHECK(found_eci, "ECI state preserved in user params");
        CHECK(found_efg, "EFG state preserved in user params");
        CHECK(found_uvw, "UVW sigmas preserved in user params");
        CHECK(found_integrator, "Integrator mode preserved");
        CHECK(found_norad, "NORAD ID preserved");
    }
}

// ── Test 9: Round-trip (text VCM → OCM → text VCM) ──

void test_roundtrip() {
    printf("\n=== Test 9: Round-Trip (text → OCM → text) ===\n");

    // Parse original
    auto result = hpop::parse_text_vcm(SAMPLE_VCM);
    CHECK(result.valid, "Initial parse OK");

    // Convert back to text
    std::string text_back = hpop::ocm_to_text_vcm(
        result.ocm_buffer.data(), result.ocm_buffer.size());
    CHECK(!text_back.empty(), "Reverse conversion produced output");
    printf("  Output text: %zu bytes\n", text_back.size());

    // Verify key fields survive the round-trip
    CHECK(text_back.find("SP VECTOR/COVARIANCE MESSAGE") != std::string::npos,
          "Header preserved");
    CHECK(text_back.find("FAKE_CENTER") != std::string::npos,
          "Center preserved");
    CHECK(text_back.find("0000-000A") != std::string::npos,
          "Int designator preserved");
    CHECK(text_back.find("369.89521989") != std::string::npos,
          "J2K X position preserved");
    CHECK(text_back.find("EGM-96") != std::string::npos,
          "Geopotential model preserved");
    CHECK(text_back.find("70Z,70T") != std::string::npos,
          "Geopotential degree/order preserved");

    // Re-parse the round-tripped text
    auto result2 = hpop::parse_text_vcm(text_back);
    CHECK(result2.valid, "Re-parsed round-tripped text");

    // Compare state vectors
    auto* ocm1 = flatbuffers::GetRoot<OCM>(result.ocm_buffer.data());
    auto* ocm2 = flatbuffers::GetRoot<OCM>(result2.ocm_buffer.data());
    auto* s1 = ocm1->STATE_DATA();
    auto* s2 = ocm2->STATE_DATA();

    if (s1 && s2 && s1->size() >= 6 && s2->size() >= 6) {
        CHECK_TOL(s2->Get(0), s1->Get(0), 1e-4, "X survives round-trip");
        CHECK_TOL(s2->Get(3), s1->Get(3), 1e-8, "Vx survives round-trip");
    }
}

// ── Test 10: Edge case — empty/malformed VCM ──

void test_edge_cases() {
    printf("\n=== Test 10: Edge Cases ===\n");

    // Empty
    auto r1 = hpop::parse_text_vcm("");
    CHECK(!r1.valid, "Empty VCM → invalid");

    // Minimal (just header + state)
    auto r2 = hpop::parse_text_vcm(
        "<> SP VECTOR/COVARIANCE MESSAGE - V2.0\n"
        "<> EPOCH TIME (UTC): 2024 001 (01 JAN) 00:00:00.000  EPOCH REV: 1\n"
        "<> J2K POS (KM):      6778.00000000       0.00000000       0.00000000\n"
        "<> J2K VEL (KM/S):   0.000000000000   7.669000000000   0.000000000000\n"
    );
    CHECK(r2.valid, "Minimal VCM parses OK");
    if (r2.valid) {
        auto* ocm = flatbuffers::GetRoot<OCM>(r2.ocm_buffer.data());
        auto* state = ocm->STATE_DATA();
        CHECK(state && state->size() == 6, "Minimal VCM has state data");
        if (state) CHECK_TOL(state->Get(0), 6778.0, 0.01, "Minimal VCM X = 6778");
    }
}

int main() {
    printf("============================================================\n");
    printf("HPOP Raw Text VCM → OCM Parser Tests\n");
    printf("============================================================\n");

    test_parse_sample();
    test_ocm_fields();
    test_state_vector();
    test_perturbations();
    test_physical_properties();
    test_covariance();
    test_iss_vcm();
    test_user_params();
    test_roundtrip();
    test_edge_cases();

    printf("\n============================================================\n");
    printf("Results: %d passed, %d failed\n", tests_passed, tests_failed);
    printf("============================================================\n");

    return tests_failed > 0 ? 1 : 0;
}
