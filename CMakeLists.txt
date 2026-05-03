cmake_minimum_required(VERSION 3.16)
project(qepp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Build Matplot++ from the bundled source as a static library so the binary
# has no dependency on the system's shared libmatplot / MKL / FFTW.
# Force ZLIB/PNG to use static archives.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "")
set(MATPLOTPP_BUILD_EXAMPLES OFF CACHE BOOL "")
set(MATPLOTPP_BUILD_TESTS OFF CACHE BOOL "")
# Point CMake directly at the static archives so FindZLIB/FindPNG don't pick
# up the .so symlinks that would break a -static link.
set(ZLIB_LIBRARY        /usr/lib/x86_64-linux-gnu/libz.a   CACHE FILEPATH "")
set(PNG_LIBRARY         /usr/lib/x86_64-linux-gnu/libpng.a CACHE FILEPATH "")
set(PNG_LIBRARY_RELEASE /usr/lib/x86_64-linux-gnu/libpng.a CACHE FILEPATH "")
add_subdirectory(third_party/matplotplusplus EXCLUDE_FROM_ALL)

# Strip any remaining shared-library paths (.so) from the cimg INTERFACE target
# — these come from the FFTW/MKL auto-detection inside Matplot++'s 3rd_party and
# would prevent a fully static link.
get_target_property(_cimg_iface cimg INTERFACE_LINK_LIBRARIES)
if(_cimg_iface)
    list(FILTER _cimg_iface EXCLUDE REGEX "\\.so")
    set_target_properties(cimg PROPERTIES INTERFACE_LINK_LIBRARIES "${_cimg_iface}")
endif()

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

# Fully static binary: bakes in glibc, libstdc++, libgcc_s, libm, libz, libpng.
# The resulting binary runs on any x86-64 Linux regardless of glibc version.
target_link_options(qepp PRIVATE -static -static-libstdc++ -static-libgcc)

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
