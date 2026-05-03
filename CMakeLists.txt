cmake_minimum_required(VERSION 3.16)
project(qepp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Matplot++ REQUIRED)
find_package(Eigen3 REQUIRED)

# Optional: spglib for space-group detection
# For Release builds the static archive is preferred; if not found on the system
# spglib is automatically downloaded and compiled as a static library via FetchContent
# so the resulting binary has no runtime dependency on libsymspg.so.
# Set -DQEPP_USE_SPGLIB=OFF to disable spglib entirely.
option(QEPP_USE_SPGLIB "Use spglib for space-group detection" ON)

if(QEPP_USE_SPGLIB)
    find_library(SPGLIB_LIB_STATIC NAMES libsymspg.a)
    find_path(SPGLIB_INC spglib.h)

    if(SPGLIB_LIB_STATIC AND SPGLIB_INC)
        # System static archive found — use it directly
        set(SPGLIB_TARGET "${SPGLIB_LIB_STATIC}")
        set(SPGLIB_TARGET_INC "${SPGLIB_INC}")
        set(SPGLIB_USE_IMPORTED ON)
        message(STATUS "spglib static library: ${SPGLIB_LIB_STATIC}")
    else()
        # No system static archive — build spglib from source (always static)
        include(FetchContent)
        set(SPGLIB_SHARED_LIBS OFF CACHE BOOL "Build spglib as static library" FORCE)
        set(BUILD_SHARED_LIBS OFF)
        FetchContent_Declare(
            spglib
            GIT_REPOSITORY https://github.com/spglib/spglib.git
            GIT_TAG        v2.3.1
            GIT_SHALLOW    ON
        )
        FetchContent_MakeAvailable(spglib)
        set(SPGLIB_TARGET Spglib::symspg)
        set(SPGLIB_TARGET_INC "")  # propagated through the imported target
        set(SPGLIB_USE_IMPORTED OFF)
        message(STATUS "spglib built from source as static library (FetchContent)")
    endif()
endif()

add_executable(
	qepp
	src/main.cpp
	src/io/cli/dispatch.cpp
	src/io/cli/cif.cpp
	src/io/cli/dos.cpp
	src/io/cli/band.cpp
	src/io/cli/elastic.cpp
	src/io/cli/charge.cpp
	src/io/cli/kpath.cpp
	src/io/help.cpp
	src/algorithm/dos.cpp
	src/algorithm/struct.cpp
	src/io/dos.cpp
	src/io/charge.cpp
	src/algorithm/utils.cpp
	src/algorithm/cif.cpp
	src/io/qe_input.cpp
	src/plot/dos.cpp
	src/io/band.cpp
	src/io/elastic.cpp
	src/algorithm/band.cpp
	src/plot/band.cpp
	src/algorithm/elastic.cpp
	src/plot/charge.cpp
	src/io/mag.cpp
	src/io/stm.cpp
	src/plot/stm.cpp
	src/io/bader.cpp
	src/io/cli/mag.cpp
	src/io/cli/stm.cpp
	src/io/cli/bader.cpp
	src/io/conv.cpp
	src/plot/conv.cpp
	src/io/struct.cpp
	src/io/parse.cpp
	src/io/cli/conv.cpp
	src/io/cli/struct.cpp
	src/io/cli/parse.cpp
	src/algorithm/qha.cpp
	src/io/qha.cpp
	src/io/cli/qha.cpp
	src/algorithm/harmonic.cpp
	src/algorithm/dfpt.cpp
	src/io/dfpt.cpp
	src/io/phonon.cpp
	src/plot/phonon.cpp
	src/io/cli/phonon.cpp
)
target_include_directories(qepp PRIVATE include)
target_link_libraries(qepp PRIVATE Matplot++::matplot Eigen3::Eigen)

if(QEPP_USE_SPGLIB AND SPGLIB_TARGET)
    target_compile_definitions(qepp PRIVATE HAVE_SPGLIB)
    if(SPGLIB_USE_IMPORTED AND SPGLIB_TARGET_INC)
        target_include_directories(qepp PRIVATE ${SPGLIB_TARGET_INC})
    endif()
    target_link_libraries(qepp PRIVATE ${SPGLIB_TARGET})
elseif(NOT QEPP_USE_SPGLIB)
    message(STATUS "spglib disabled — crystal family detection will use ibrav fallback")
else()
    message(STATUS "spglib NOT found — crystal family detection will use ibrav fallback")
endif()
