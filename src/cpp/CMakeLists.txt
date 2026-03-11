cmake_minimum_required(VERSION 3.15)
project(hpop_sdn_plugin LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Performance flags ──
if(EMSCRIPTEN)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -msimd128 -flto -ffast-math -DNDEBUG")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mrelaxed-simd")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O3 -flto --closure 1")
else()
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -flto -ffast-math -march=native -mtune=native")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# HPOP library
set(HPOP_SOURCES
    src/propagator.cpp
    src/gravity.cpp
    src/atmosphere.cpp
    src/srp.cpp
    src/thirdbody.cpp
)

add_library(hpop_lib STATIC ${HPOP_SOURCES})
target_include_directories(hpop_lib PUBLIC include)

# Test executables
if(NOT EMSCRIPTEN)
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
endif()
