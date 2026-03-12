cmake_minimum_required(VERSION 3.15)
project(hpop_sdn_plugin LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ============================================================================
# Locate spacedatastandards.org dependencies (for OEM FlatBuffers)
# ============================================================================

if(NOT DEFINED SDS_ROOT)
    if(DEFINED ENV{SDS_ROOT})
        set(SDS_ROOT "$ENV{SDS_ROOT}")
    else()
        foreach(_candidate
            "${CMAKE_SOURCE_DIR}/../../../../space-data-network-project"
            "$ENV{HOME}/software/space-data-network-project"
        )
            if(EXISTS "${_candidate}/packages/flatbuffers/include")
                get_filename_component(SDS_ROOT "${_candidate}" ABSOLUTE)
                break()
            endif()
        endforeach()
    endif()
endif()

if(SDS_ROOT)
    message(STATUS "SDS_ROOT: ${SDS_ROOT}")
    set(FLATBUFFERS_INCLUDE "${SDS_ROOT}/packages/flatbuffers/include")
    set(OEM_GENERATED "${SDS_ROOT}/packages/spacedatastandards.org/lib/cpp/OEM")
else()
    message(WARNING "SDS_ROOT not found. FlatBuffers headers may be missing.")
    set(FLATBUFFERS_INCLUDE "")
    set(OEM_GENERATED "")
endif()

# ── Performance flags ──
if(EMSCRIPTEN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msimd128 -flto -ffast-math -DNDEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mrelaxed-simd")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -flto --closure 1")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -flto -ffast-math -march=native -mtune=native")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# ============================================================================
# HPOP core library
# ============================================================================

set(HPOP_SOURCES
    src/propagator.cpp
    src/gravity.cpp
    src/atmosphere.cpp
    src/srp.cpp
    src/thirdbody.cpp
)

add_library(hpop_lib STATIC ${HPOP_SOURCES})
target_include_directories(hpop_lib PUBLIC include)

# ============================================================================
# Emscripten WASM build
# ============================================================================

if(EMSCRIPTEN)
    add_executable(hpop_wasm wasm_api.cpp)
    target_link_libraries(hpop_wasm hpop_lib)
    target_include_directories(hpop_wasm PRIVATE include)

    if(FLATBUFFERS_INCLUDE)
        target_include_directories(hpop_wasm PRIVATE ${FLATBUFFERS_INCLUDE})
    endif()
    if(OEM_GENERATED)
        target_include_directories(hpop_wasm PRIVATE ${OEM_GENERATED})
    endif()

    set_target_properties(hpop_wasm PROPERTIES
        SUFFIX ".js"
        LINK_FLAGS "\
            -s WASM=1 \
            -s MODULARIZE=1 \
            -s EXPORT_ES6=1 \
            -s EXPORT_NAME='HPOP_WASM' \
            -s EXPORTED_RUNTIME_METHODS=['ccall','cwrap','UTF8ToString','stringToUTF8','lengthBytesUTF8'] \
            -s EXPORTED_FUNCTIONS=['_parse','_convert','_malloc','_free'] \
            -s ALLOW_MEMORY_GROWTH=1 \
            -s NO_EXIT_RUNTIME=1 \
            -s ENVIRONMENT='web,node' \
            -s DISABLE_EXCEPTION_CATCHING=0 \
            --bind \
            -O2 \
        "
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/../../wasm/node"
    )

# ============================================================================
# Native build (for testing)
# ============================================================================

else()
    add_executable(test_twobody tests/test_twobody.cpp)
    target_link_libraries(test_twobody hpop_lib)

    add_executable(test_geopotential tests/test_geopotential.cpp)
    target_link_libraries(test_geopotential hpop_lib)

    add_executable(test_drag tests/test_drag.cpp)
    target_link_libraries(test_drag hpop_lib)

    add_executable(test_srp tests/test_srp.cpp)
    target_link_libraries(test_srp hpop_lib)

    add_executable(test_thirdbody tests/test_thirdbody.cpp)
    target_link_libraries(test_thirdbody hpop_lib)

    # Native WASM API test (propagation + OEM FlatBuffer output)
    if(FLATBUFFERS_INCLUDE AND OEM_GENERATED)
        add_executable(test_wasm_api tests/test_wasm_api.cpp)
        target_link_libraries(test_wasm_api hpop_lib)
        target_include_directories(test_wasm_api PRIVATE include ${FLATBUFFERS_INCLUDE} ${OEM_GENERATED})
    endif()
endif()
