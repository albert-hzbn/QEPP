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
	src/utils.cpp
	src/cif.cpp
	src/qe_input.cpp
	src/dos.cpp
	src/band.cpp
	src/elastic.cpp
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
