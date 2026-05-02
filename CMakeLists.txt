cmake_minimum_required(VERSION 3.16)
project(qepp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Matplot++ REQUIRED)
find_package(Eigen3 REQUIRED)

# Optional: spglib for space-group detection
find_library(SPGLIB_LIB symspg)
find_path(SPGLIB_INC spglib.h)

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
	src/io/dos.cpp
	src/io/charge.cpp
	src/algorithm/utils.cpp
	src/algorithm/cif.cpp
	src/io/qe_input.cpp
	src/plot/dos.cpp
	src/io/band.cpp
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
)
target_include_directories(qepp PRIVATE include)
target_link_libraries(qepp PRIVATE Matplot++::matplot Eigen3::Eigen)

if(SPGLIB_LIB AND SPGLIB_INC)
    target_compile_definitions(qepp PRIVATE HAVE_SPGLIB)
    target_include_directories(qepp PRIVATE ${SPGLIB_INC})
    target_link_libraries(qepp PRIVATE ${SPGLIB_LIB})
    message(STATUS "spglib found: ${SPGLIB_LIB}")
else()
    message(STATUS "spglib NOT found — crystal family detection will use ibrav fallback")
endif()
