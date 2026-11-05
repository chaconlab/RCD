cmake_minimum_required(VERSION 3.16)
project(RCD VERSION 1.06 LANGUAGES CXX)

# ---------------------------------------------------------------------------
# Flavors
#   RCD_MPI    : OFF (default) = single processor, ON = MPI build (defines MPI_RELEASE)
#   RCD_INTEL  : OFF (default) = use the compiler CMake finds (gcc, ...),
#                ON = use Intel oneAPI icpx (+ MKL instead of LAPACK)
#
# Examples:
#   cmake -B build-gcc                          && cmake --build build-gcc -j
#   cmake -B build-mpi   -DRCD_MPI=ON           && cmake --build build-mpi -j
#   cmake -B build-intel -DRCD_INTEL=ON         && cmake --build build-intel -j
#   cmake -B build-impi  -DRCD_INTEL=ON -DRCD_MPI=ON
# ---------------------------------------------------------------------------
option(RCD_MPI   "Build the MPI version of rcd"        OFF)
option(RCD_INTEL "Build with the Intel compiler (icpx)" OFF)

if(RCD_INTEL AND NOT CMAKE_CXX_COMPILER)
  find_program(_ICPX icpx REQUIRED)
  set(CMAKE_CXX_COMPILER ${_ICPX})
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Output layout identical to the old makefiles
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Flavor suffix, as in the old rcd_mpi_gcc / lib*_gcc naming
if(RCD_INTEL)
  set(_toolchain intel)
else()
  set(_toolchain gcc)
endif()

set(SRC ${CMAKE_SOURCE_DIR}/src)

# --- compiler flags ---------------------------------------------------------
if(RCD_INTEL)
  set(RCD_OPT   -O3 -g -inline-level=2)
else()
  set(RCD_OPT   -O3 -g -Wall -fmessage-length=0)
  list(APPEND RCD_OPT -fpermissive -U_FORTIFY_SOURCE)
endif()

# --- helper to define a static library --------------------------------------
function(rcd_lib name dir)
  file(GLOB srcs CONFIGURE_DEPENDS ${SRC}/${dir}/*.cpp)
  add_library(${name} STATIC ${srcs})
  target_include_directories(${name} PUBLIC ${SRC}/${dir}/include)
  target_include_directories(${name} PRIVATE ${SRC})
  target_compile_options(${name} PRIVATE ${RCD_OPT})
  set_target_properties(${name} PROPERTIES POSITION_INDEPENDENT_CODE ON)
endfunction()

rcd_lib(tools  libtools)
rcd_lib(pdb    libpdb)
rcd_lib(volume libvolume)
rcd_lib(energy libenergy)
rcd_lib(ga     libga)

# Sources include headers as <libpdb/...> style via src/ and includes/
foreach(t tools pdb volume energy ga)
  target_include_directories(${t} PRIVATE ${CMAKE_SOURCE_DIR}/includes)
endforeach()


target_compile_definitions(pdb PRIVATE VOLUME_INCLUDED)
target_compile_definitions(ga  PRIVATE VOLUME_INCLUDED)

# --- linear algebra ---------------------------------------------------------
if(RCD_INTEL)
  target_compile_options(volume PRIVATE -qmkl=sequential)
  set(RCD_LAPACK -qmkl=sequential)
else()
  find_package(LAPACK REQUIRED)
  set(RCD_LAPACK ${LAPACK_LIBRARIES})
endif()

# --- executable -------------------------------------------------------------
add_executable(rcd ${SRC}/rcd/rcd.cpp ${SRC}/rcd/rt_timer.cpp)
target_include_directories(rcd PRIVATE ${SRC}/rcd/include ${SRC} ${CMAKE_SOURCE_DIR}/includes)
target_compile_options(rcd PRIVATE ${RCD_OPT})
target_compile_definitions(rcd PRIVATE VOLUME_INCLUDED)
target_link_libraries(rcd PRIVATE ga energy pdb volume tools ${RCD_LAPACK} rt)

if(RCD_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
  target_compile_definitions(rcd PRIVATE MPI_RELEASE)
  target_link_libraries(rcd PRIVATE MPI::MPI_CXX)
  set(_exe rcd_mpi_${_toolchain})
else()
  set(_exe rcd_${_toolchain})
endif()
set_target_properties(rcd PROPERTIES OUTPUT_NAME ${_exe})

# Data files the program expects next to the sources
install(TARGETS rcd RUNTIME DESTINATION bin)
install(FILES ${SRC}/rcd/CDRS.bin ${SRC}/rcd/dunbrackTCB.bin ${SRC}/rcd/dunbrackTCBIG.bin
              ${SRC}/libenergy/dunbrack.bin ${SRC}/libenergy/loco.score
        DESTINATION share/rcd)
