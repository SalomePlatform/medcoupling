# Copyright (C) 2015-2026  CEA, EDF
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this library; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email :
# webmaster.salome@opencascade.com
#

cmake_minimum_required(VERSION 2.8.11 FATAL_ERROR)

# Project name
# ============
# original
project(MEDCoupling C CXX)
# upper case
string(TOUPPER ${PROJECT_NAME} PROJECT_NAME_UC)

if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.3)
  cmake_policy(SET CMP0057 NEW)
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13)
    cmake_policy(SET CMP0078 OLD)
  endif()
  if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.14)
    cmake_policy(SET CMP0086 OLD)
  endif()
endif()

# Common CMake macros
set(CONFIGURATION_ROOT_DIR
    $ENV{CONFIGURATION_ROOT_DIR}
    CACHE PATH "Path to the Salome CMake files")
if(EXISTS ${CONFIGURATION_ROOT_DIR})
  list(APPEND CMAKE_MODULE_PATH "${CONFIGURATION_ROOT_DIR}/cmake")
  include(SalomeMacros NO_POLICY_SCOPE)
else()
  message(
    FATAL_ERROR
      "We absolutely need the Salome CMake configuration files, please define CONFIGURATION_ROOT_DIR !"
  )
endif()

# Versioning
# ===========
salome_setup_version(9.15.0)
message(
  STATUS
    "Building ${PROJECT_NAME_UC} ${${PROJECT_NAME_UC}_VERSION} from \"${${PROJECT_NAME_UC}_GIT_SHA1}\""
)

# Platform setup
# ==============
include(SalomeSetupPlatform)

# [ABN]: use the below for aggressive code quality check: if( EXISTS "aaaa")
# add_definitions(-Weverything)
# add_definitions(-Wno-inconsistent-missing-override)
# add_definitions(-Wno-c++98-compat -Wno-c++98-compat-pedantic)
# add_definitions(-Wsign-conversion) add_definitions(-Wno-switch-enum)
#
# add_definitions(-Wno-documentation-unknown-command) # \image in Doxygen
# add_definitions(-Wno-gnu-statement-expression)      # assert(...)
#
# add_definitions(-Wno-reserved-id-macro -Wno-padded -Wno-weak-vtables
# -Wdouble-promotion -Wno-comma -Wno-unused-parameter)
# add_definitions(-Wno-unreachable-code-break -Wno-old-style-cast
# -Wno-deprecated -Wno-double-promotion)
# add_definitions(-Wno-exit-time-destructors -Wno-global-constructors
# -Wno-missing-prototypes) add_definitions(-Wno-covered-switch-default
# -Wno-weak-template-vtables -Wno-undefined-func-template
# -Wno-used-but-marked-unused) add_definitions(-Wno-unreachable-code-return
# -Wno-missing-noreturn -Wno-unused-member-function -Wno-header-hygiene)
# add_definitions(-Wno-ignored-qualifiers) # because of MED-file
# add_definitions(-Wno-unused-function -Wno-unused-macros)
# add_definitions(-Wno-missing-field-initializers) add_definitions(-Wno-shadow)
# # SWIG generates a lot of those ...
#
# add_definitions(-Wno-float-equal)  # OUCH endif()

#
# User options
# ============
include(CMakeDependentOption)
option(MEDCOUPLING_MICROMED "Build MED without MED file dependency." OFF)
option(MEDCOUPLING_ENABLE_PYTHON "Build PYTHON bindings." ON)
option(MEDCOUPLING_ENABLE_PARTITIONER "Build MEDPartitioner." ON)
option(MEDCOUPLING_ENABLE_RENUMBER "Build Renumber." ON)
option(MEDCOUPLING_ENABLE_SHAPERECOGN "Build Shape Recognition module." OFF)
option(
  MEDCOUPLING_WITH_FILE_EXAMPLES
  "Install examples of files containing meshes and fields of different formats."
  ON)
option(MEDCOUPLING_USE_MPI
       "(Use MPI containers) - For MED this triggers the build of ParaMEDMEM."
       OFF)
option(MEDCOUPLING_BUILD_TESTS "Build MEDCoupling C++ tests." ON)
option(MEDCOUPLING_BUILD_PY_TESTS "Build MEDCoupling Python tests." ON)
option(MEDCOUPLING_BUILD_DOC "Build MEDCoupling doc." ON)
option(MEDCOUPLING_BUILD_STATIC "Build MEDCoupling library in static mode." OFF)
option(
  MEDCOUPLING_USE_64BIT_IDS
  "Size of IDs to refer cells and nodes. 32 bits when OFF (default), 64 bits when ON."
  ON)

if(${MEDCOUPLING_USE_MPI})
  set(USE_METIS_NOT_PARMETIS OFF)
else()
  set(USE_METIS_NOT_PARMETIS ON)
endif()
cmake_dependent_option(
  MEDCOUPLING_PARTITIONER_METIS "Enable metis graph library in MEDPartitioner."
  ${USE_METIS_NOT_PARMETIS} "MEDCOUPLING_ENABLE_PARTITIONER" OFF)
cmake_dependent_option(
  MEDCOUPLING_PARTITIONER_SCOTCH
  "Enable scotch graph library in MEDPartitioner." ON
  "MEDCOUPLING_ENABLE_PARTITIONER" OFF)
cmake_dependent_option(
  MEDCOUPLING_PARTITIONER_PARMETIS
  "Enable parmetis graph library in MEDPartitioner." ON
  "MEDCOUPLING_ENABLE_PARTITIONER;MEDCOUPLING_USE_MPI" OFF)
cmake_dependent_option(
  MEDCOUPLING_PARTITIONER_PTSCOTCH
  "Enable ptscotch graph library in MEDPartitioner." ON
  "MEDCOUPLING_ENABLE_PARTITIONER;MEDCOUPLING_USE_MPI" OFF)
set(XDR_DEFAULT_OPTION ON)
if(WIN32)
  set(XDR_DEFAULT_OPTION OFF)
endif()
cmake_dependent_option(
  MEDCOUPLING_MEDLOADER_USE_XDR "Enable use of XDR for SauvReader."
  ${XDR_DEFAULT_OPTION} "NOT MEDCOUPLING_MICROMED" OFF)
cmake_dependent_option(
  MEDCOUPLING_BUILD_FRENCH_DOC "Generate MEDCOUPLING French documentation" OFF
  "MEDCOUPLING_BUILD_DOC" OFF)

if(MEDCOUPLING_BUILD_STATIC)
  set(BUILD_SHARED_LIBS 0)
  set(CMAKE_POSITION_INDEPENDENT_CODE 1) # -fPIC option
else(MEDCOUPLING_BUILD_STATIC)
  set(BUILD_SHARED_LIBS 1)
endif(MEDCOUPLING_BUILD_STATIC)

if(${MEDCOUPLING_PARTITIONER_PARMETIS} AND ${MEDCOUPLING_PARTITIONER_METIS})
  message(
    FATAL_ERROR
      "ParMetis and Metis are mutually exclusive! Make a choice (options MEDCOUPLING_PARTITIONER_METIS and MEDCOUPLING_PARTITIONER_PARMETIS)."
  )
endif(${MEDCOUPLING_PARTITIONER_PARMETIS} AND ${MEDCOUPLING_PARTITIONER_METIS})

if(${MEDCOUPLING_PARTITIONER_PTSCOTCH} AND ${MEDCOUPLING_PARTITIONER_SCOTCH})
  message(
    FATAL_ERROR
      "PTScotch and Scotch are mutually exclusive! Make a choice (options MEDCOUPLING_PARTITIONER_SCOTCH and MEDCOUPLING_PARTITIONER_PTSCOTCH)."
  )
endif(${MEDCOUPLING_PARTITIONER_PTSCOTCH} AND ${MEDCOUPLING_PARTITIONER_SCOTCH})

#
# Set list of prerequisites
# =========================

if(NOT MEDCOUPLING_MICROMED)
  find_package(SalomeHDF5 REQUIRED)
  find_package(SalomeMEDFile REQUIRED)
  # XDR stuff
  if(MEDCOUPLING_MEDLOADER_USE_XDR)
    find_package(SalomeXDR REQUIRED)
  endif(MEDCOUPLING_MEDLOADER_USE_XDR)
  # End of XDR Stuff
  if(MEDCOUPLING_ENABLE_PARTITIONER)
    find_package(SalomeLibXml2)
    salome_log_optional_package(LibXml2 MEDCOUPLING_ENABLE_PARTITIONER)
    if(MEDCOUPLING_PARTITIONER_METIS)
      find_package(SalomeMetis)
      salome_log_optional_package(Metis MEDCOUPLING_PARTITIONER_METIS)
      if(SalomeMetis_FOUND)
        add_definitions("-DMED_ENABLE_METIS")
      endif()
    endif(MEDCOUPLING_PARTITIONER_METIS)
    if(MEDCOUPLING_PARTITIONER_SCOTCH)
      find_package(SalomeScotch)
      salome_log_optional_package(Scotch MEDCOUPLING_PARTITIONER_SCOTCH)
      if(SalomeScotch_FOUND)
        add_definitions("-DMED_ENABLE_SCOTCH")
      endif()
    endif(MEDCOUPLING_PARTITIONER_SCOTCH)
  endif(MEDCOUPLING_ENABLE_PARTITIONER)
endif(NOT MEDCOUPLING_MICROMED)

enable_testing(
  ) # let it outsite because even if MEDCOUPLING_BUILD_TESTS is OFF, python
    # tests that not need additional compilation can be run.

if(MEDCOUPLING_BUILD_TESTS)
  find_package(SalomeCppUnit)
  salome_log_optional_package(CppUnit SALOME_BUILD_TESTS)
endif(MEDCOUPLING_BUILD_TESTS)

if(MEDCOUPLING_USE_MPI)
  find_package(SalomeMPI REQUIRED)
  # Used when exporting configuration:
  set(MEDCOUPLING_CXX_COMPILER "${CMAKE_CXX_COMPILER}")
  set(MEDCOUPLING_MPI_CXX_COMPILER "${MPI_CXX_COMPILER}") # Might not always be
                                                          # set even if MPI
                                                          # there ...
  set(MEDCOUPLING_MPI_CXX_LIBRARIES "${MPI_CXX_LIBRARIES}")

  add_definitions("-DHAVE_MPI")
  salome_add_mpi_to_hdf5()
  if(MEDCOUPLING_PARTITIONER_PARMETIS)
    find_package(SalomeParMetis)
    salome_log_optional_package(ParMetis MEDCOUPLING_PARTITIONER_PARMETIS)
    if(SalomeParMetis_FOUND)
      add_definitions("-DMED_ENABLE_PARMETIS")
    endif()
  endif(MEDCOUPLING_PARTITIONER_PARMETIS)
  if(MEDCOUPLING_PARTITIONER_PTSCOTCH)
    find_package(SalomePTScotch)
    salome_log_optional_package(PTScotch MEDCOUPLING_PARTITIONER_PTSCOTCH)
    if(SalomePTScotch_FOUND)
      add_definitions("-DMED_ENABLE_PTSCOTCH")
    endif(SalomePTScotch_FOUND)
  endif(MEDCOUPLING_PARTITIONER_PTSCOTCH)
endif(MEDCOUPLING_USE_MPI)

if(MEDCOUPLING_ENABLE_RENUMBER)
  find_package(SalomeBoost)
  salome_log_optional_package(Boost MEDCOUPLING_ENABLE_RENUMBER)
endif(MEDCOUPLING_ENABLE_RENUMBER)

if(MEDCOUPLING_ENABLE_SHAPERECOGN)
  if(UNIX)
    set(CBLAS_ROOT_DIR
        $ENV{CBLAS_ROOT_DIR}
        CACHE PATH "Path to the LAPACK/CBLAS.")
    if(CBLAS_ROOT_DIR)
      list(APPEND CMAKE_PREFIX_PATH "${CBLAS_ROOT_DIR}")
      find_package(CBLAS REQUIRED)
    else(CBLAS_ROOT_DIR)
      find_package(BLAS REQUIRED)
    endif(CBLAS_ROOT_DIR)
    set(LAPACK_ROOT_DIR
        $ENV{LAPACK_ROOT_DIR}
        CACHE PATH "Path to the LAPACKE")
    if(LAPACK_ROOT_DIR)
      set(LAPACK_ROOT_DIR
          $ENV{LAPACK_ROOT_DIR}
          CACHE PATH "Path to the LAPACK.")
      list(APPEND CMAKE_PREFIX_PATH "${LAPACK_ROOT_DIR}")
    endif(LAPACK_ROOT_DIR)
    find_package(LAPACK REQUIRED)
    find_library(LAPACKE_LIB NAMES lapacke REQUIRED)
    set(LAPACK_LIBRARIES ${LAPACKE_LIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES}
                         ${CBLAS_LIBRARIES})
    find_path(
      LAPACKE_INCLUDE_DIRS
      NAMES lapacke.h
      HINTS ${LAPACK_LIBRARIES})
    if(LAPACK_FOUND)
      message(STATUS "Lapacke libraries: ${LAPACK_LIBRARIES}")
      message(STATUS "Lapacke include dirs: ${LAPACKE_INCLUDE_DIRS}")
    else()
      message(FATAL_ERROR "Error in Lapacke detection ! lapacke not found !")
    endif(LAPACK_FOUND)
    salome_log_optional_package(Lapack MEDCOUPLING_ENABLE_SHAPERECOGN)
  elseif(MSVC)
    set(OPENBLAS_ROOT_DIR
        $ENV{OPENBLAS_ROOT_DIR}
        CACHE PATH "Path to the OpenBLAS.")
    if(OPENBLAS_ROOT_DIR)
      list(APPEND CMAKE_PREFIX_PATH "${OPENBLAS_ROOT_DIR}")
      find_package(OpenBLAS REQUIRED)
    else(OPENBLAS_ROOT_DIR)
      message(FATAL "Could not find OpenBLAS")
    endif(OPENBLAS_ROOT_DIR)
    set(LAPACK_ROOT_DIR
        $ENV{LAPACK_ROOT_DIR}
        CACHE PATH "Path to the LAPACKE")
    if(LAPACK_ROOT_DIR)
      set(LAPACK_ROOT_DIR
          $ENV{LAPACK_ROOT_DIR}
          CACHE PATH "Path to the LAPACK.")
      list(APPEND CMAKE_PREFIX_PATH "${LAPACK_ROOT_DIR}")
    endif(LAPACK_ROOT_DIR)
    find_package(LAPACK REQUIRED)
    find_library(LAPACKE_LIB NAMES openblas REQUIRED)
    set(LAPACK_LIBRARIES ${LAPACKE_LIB} ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
    find_path(
      LAPACKE_INCLUDE_DIRS
      NAMES lapacke.h
      HINTS ${LAPACK_LIBRARIES} ${OpenBLAS_INCLUDE_DIRS})
    if(LAPACK_FOUND)
      message(STATUS "Lapacke libraries: ${LAPACK_LIBRARIES}")
      message(STATUS "Lapacke include dirs: ${LAPACKE_INCLUDE_DIRS}")
    else()
      message(FATAL_ERROR "Error in Lapacke detection ! lapacke not found !")
    endif(LAPACK_FOUND)
    salome_log_optional_package(Lapack MEDCOUPLING_ENABLE_SHAPERECOGN)

  else()
    message(FATAL "Not implemented")
  endif()
endif(MEDCOUPLING_ENABLE_SHAPERECOGN)

if(MEDCOUPLING_ENABLE_PYTHON)
  find_package(SalomePythonInterp)
  find_package(SalomePythonLibs)
  find_package(SalomeNumPySciPy)
  find_package(SalomeSWIG)
  salome_log_optional_package(PythonInterp MEDCOUPLING_ENABLE_PYTHON)
  salome_log_optional_package(PythonLibs MEDCOUPLING_ENABLE_PYTHON)
  salome_log_optional_package(SWIG MEDCOUPLING_ENABLE_PYTHON)
  if("${PYTHON_VERSION_MAJOR}" STREQUAL "2")
    if("${PYTHON_VERSION_MINOR}" LESS "7")
      message(
        FATAL_ERROR
          "MEDCoupling's Python requires at least Python 2.7 (you seem to have ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}). Upgrade your Python, or turn off MEDCOUPLING_ENABLE_PYTHON"
      )
    endif()
  endif()
  if("${PYTHON_VERSION_MAJOR}" STREQUAL "3")
    message(STATUS "Using Python 3")
  endif()
  if("${SWIG_VERSION}" VERSION_GREATER_EQUAL "4")
    set(MEDCOUPLING_SWIG4_COMPAT
        TRUE
        CACHE BOOL "SWIG 4 compatibility (experimental)")
  else()
    set(MEDCOUPLING_SWIG4_COMPAT FALSE)
  endif()
  if(MEDCOUPLING_SWIG4_COMPAT)
    list(APPEND CMAKE_SWIG_FLAGS "-DMEDCOUPLING_SWIG4_COMPAT")
  endif()
endif(MEDCOUPLING_ENABLE_PYTHON)

if(MEDCOUPLING_BUILD_DOC)
  find_package(SalomeDoxygen)
  find_package(SalomeGraphviz)
  find_package(SalomeSphinx)
  salome_log_optional_package(Doxygen MEDCOUPLING_BUILD_DOC)
  salome_log_optional_package(Graphviz MEDCOUPLING_BUILD_DOC)
  salome_log_optional_package(Sphinx MEDCOUPLING_BUILD_DOC)
endif(MEDCOUPLING_BUILD_DOC)

# Detection report
salome_package_report_and_check()

# Directories
#
# Directories have to be given after prerequisites (to be able to use Python
# version string for example).
# ===========
set(MEDCOUPLING_INSTALL_BINS
    bin
    CACHE PATH "Install path: MEDCoupling binaries")
set(MEDCOUPLING_INSTALL_LIBS
    lib
    CACHE PATH "Install path: MEDCoupling libs")
set(MEDCOUPLING_INSTALL_HEADERS
    include
    CACHE PATH "Install path: MEDCoupling headers")
set(MEDCOUPLING_INSTALL_SCRIPT_SCRIPTS
    ${MEDCOUPLING_INSTALL_BINS}
    CACHE PATH "Install path: MEDCoupling scripts")
set(MEDCOUPLING_INSTALL_TESTS
    bin/salome/test
    CACHE PATH "Install path: MEDCoupling tests")
# SET(MEDCOUPLING_INSTALL_SCRIPT_DATA ${MEDCOUPLING_INSTALL_BINS} CACHE PATH
# "Install path: MEDCoupling script data")
set(MEDCOUPLING_INSTALL_SCRIPT_PYTHON
    ${MEDCOUPLING_INSTALL_BINS}
    CACHE PATH "Install path: MEDCoupling Python scripts")
set(MEDCOUPLING_INSTALL_CMAKE_LOCAL
    cmake_files
    CACHE PATH "Install path: local MEDCoupling CMake files")

if(MEDCOUPLING_ENABLE_PYTHON)
  set(_pydir
      lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages)
  set(MEDCOUPLING_INSTALL_PYTHON
      ${_pydir}
      CACHE PATH "Install path: MEDCoupling Python stuff")
  set(MEDCOUPLING_INSTALL_PYTHON_SHARED
      ${MEDCOUPLING_INSTALL_PYTHON}/shared_modules
      CACHE PATH "Install path: MEDCoupling Python shared modules")
endif(MEDCOUPLING_ENABLE_PYTHON)

set(MEDCOUPLING_INSTALL_RES
    share/resources
    CACHE PATH "Install path: MEDCoupling resources")
set(MEDCOUPLING_INSTALL_DOC
    share/doc
    CACHE PATH "Install path: MEDCoupling documentation")

# Med specific:

set(MEDCOUPLING_INSTALL_RES_DATA
    "${MEDCOUPLING_INSTALL_RES}/med"
    CACHE PATH "Install path: MEDCoupling specific data")
# SET(MEDCOUPLING_INSTALL_RES_SCRIPTS "${MEDCOUPLING_INSTALL_RES}/med" CACHE
# PATH "Install path: MEDCouplng specific scripts")

mark_as_advanced(MEDCOUPLING_INSTALL_BINS MEDCOUPLING_INSTALL_LIBS
                 MEDCOUPLING_INSTALL_IDLS MEDCOUPLING_INSTALL_HEADERS)
mark_as_advanced(
  MEDCOUPLING_INSTALL_SCRIPT_SCRIPTS MEDCOUPLING_INSTALL_SCRIPT_DATA
  MEDCOUPLING_INSTALL_SCRIPT_PYTHON)
mark_as_advanced(
  MEDCOUPLING_INSTALL_APPLISKEL_SCRIPTS MEDCOUPLING_INSTALL_APPLISKEL_PYTHON
  MEDCOUPLING_INSTALL_CMAKE MEDCOUPLING_INSTALL_CMAKE_LOCAL
  MEDCOUPLING_INSTALL_RES)
mark_as_advanced(MEDCOUPLING_INSTALL_PYTHON MEDCOUPLING_INSTALL_PYTHON_SHARED
                 MEDCOUPLING_INSTALL_RES_DATA MEDCOUPLING_INSTALL_DOC)

# Header configuration
# ====================
salome_configure_file(MEDCoupling_version.h.in MEDCoupling_version.h INSTALL
                      ${MEDCOUPLING_INSTALL_HEADERS})

# Accumulate environment variables for MED module
# SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH NOCHECK
# ${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_BINS}
# ${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_PYTHON})
# SALOME_ACCUMULATE_ENVIRONMENT(LD_LIBRARY_PATH NOCHECK
# ${CMAKE_INSTALL_PREFIX}/${MEDCOUPLING_INSTALL_LIBS})

# Sources
# ========
if(WIN32)
  add_definitions("-D_USE_MATH_DEFINES")
endif(WIN32)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

if(MEDCOUPLING_ENABLE_PYTHON AND NUMPY_FOUND)
  add_definitions("-DNPY_NO_DEPRECATED_API")
endif(MEDCOUPLING_ENABLE_PYTHON AND NUMPY_FOUND)

add_subdirectory(src)

if(MEDCOUPLING_BUILD_DOC)
  add_subdirectory(doc)
endif(MEDCOUPLING_BUILD_DOC)

add_subdirectory(resources)
if(MEDCOUPLING_ENABLE_PYTHON)
  add_subdirectory(v8_work)
endif()

# Configuration export
# ====================
include(CMakePackageConfigHelpers)

# List of targets in this project we want to make visible to the rest of the
# world. They all have to be INSTALL'd with the option "EXPORT
# ${PROJECT_NAME}TargetGroup"
set(_${PROJECT_NAME}_exposed_targets interpkernel medcouplingcpp
                                     medcouplingremapper medicoco)

if(NOT MEDCOUPLING_MICROMED)
  list(APPEND _${PROJECT_NAME}_exposed_targets medloader)
  if(MEDCOUPLING_ENABLE_RENUMBER)
    list(APPEND _${PROJECT_NAME}_exposed_targets renumbercpp)
  endif()
  if(MEDCOUPLING_ENABLE_PARTITIONER)
    list(APPEND _${PROJECT_NAME}_exposed_targets medpartitionercpp)
    if(MEDCOUPLING_BUILD_TESTS)
      list(APPEND _${PROJECT_NAME}_exposed_targets MEDPARTITIONERTest)
    endif()
  endif()
  if(MEDCOUPLING_BUILD_TESTS)
    list(APPEND _${PROJECT_NAME}_exposed_targets InterpKernelTest
         InterpKernelTestUtils)
  endif()
endif()

if(MEDCOUPLING_USE_MPI)
  list(APPEND _${PROJECT_NAME}_exposed_targets paramedmem)
  if(NOT MEDCOUPLING_MICROMED)
    list(APPEND _${PROJECT_NAME}_exposed_targets paramedloader)
  endif()
  if(MEDCOUPLING_BUILD_TESTS)
    if(NOT MEDCOUPLING_MICROMED)
      list(APPEND _${PROJECT_NAME}_exposed_targets ParaMEDMEMTest)
    endif()
  endif()
endif()

if(MEDCOUPLING_ENABLE_SHAPERECOGN)
  list(APPEND _${PROJECT_NAME}_exposed_targets shaperecogn)
endif(MEDCOUPLING_ENABLE_SHAPERECOGN)

# Add all targets to the build-tree export set

export(TARGETS ${_${PROJECT_NAME}_exposed_targets}
       FILE ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake)

# Create the configuration files:

# Ensure the variables are always defined for the configure (even if empty):
set(MEDFILE_ROOT_DIR "${MEDFILE_ROOT_DIR}")
set(MPI_ROOT_DIR "${MPI_ROOT_DIR}")
set(BOOST_ROOT_DIR "${BOOST_ROOT_DIR}")
set(SWIG_ROOT_DIR "${SWIG_ROOT_DIR}")
set(PYTHON_ROOT_DIR "${PYTHON_ROOT_DIR}")
set(CPPUNIT_ROOT_DIR "${CPPUNIT_ROOT_DIR}")
set(GRAPHVIZ_ROOT_DIR "${GRAPHVIZ_ROOT_DIR}")
set(DOXYGEN_ROOT_DIR "${DOXYGEN_ROOT_DIR}")
set(SPHINX_ROOT_DIR "${SPHINX_ROOT_DIR}")

set(METIS_ROOT_DIR "${METIS_ROOT_DIR}")
set(PARMETIS_ROOT_DIR "${PARMETIS_ROOT_DIR}")
set(SCOTCH_ROOT_DIR "${SCOTCH_ROOT_DIR}")
set(XDR_ROOT_DIR "${XDR_ROOT_DIR}")

set(CONF_INCLUDE_DIRS "${CMAKE_INSTALL_PREFIX}/${INSTALL_INCLUDE_DIR}")

# Build variables that will be expanded when configuring
# Salome<MODULE>Config.cmake:
salome_configure_prepare(
  MEDFile
  MPI
  Boost
  Swig
  Python
  CppUnit
  Graphviz
  Doxygen
  Sphinx
  Metis
  ParMetis
  Scotch
  XDR)

configure_package_config_file(
  ${PROJECT_NAME}Config.cmake.in
  ${PROJECT_BINARY_DIR}/to_install/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION "${MEDCOUPLING_INSTALL_CMAKE_LOCAL}"
  PATH_VARS CONF_INCLUDE_DIRS
            MEDCOUPLING_INSTALL_CMAKE_LOCAL
            CMAKE_INSTALL_PREFIX
            MEDFILE_ROOT_DIR
            MPI_ROOT_DIR
            BOOST_ROOT_DIR
            SWIG_ROOT_DIR
            PYTHON_ROOT_DIR
            CPPUNIT_ROOT_DIR
            GRAPHVIZ_ROOT_DIR
            DOXYGEN_ROOT_DIR
            SPHINX_ROOT_DIR
            METIS_ROOT_DIR
            PARMETIS_ROOT_DIR
            SCOTCH_ROOT_DIR
            XDR_ROOT_DIR)

write_basic_package_version_file(
  ${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
  VERSION ${${PROJECT_NAME_UC}_VERSION}
  COMPATIBILITY AnyNewerVersion)

# Install the CMake configuration files:
install(FILES "${PROJECT_BINARY_DIR}/to_install/${PROJECT_NAME}Config.cmake"
              "${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake"
        DESTINATION "${MEDCOUPLING_INSTALL_CMAKE_LOCAL}")

# Install the export set for use with the install-tree
install(
  EXPORT ${PROJECT_NAME}TargetGroup
  DESTINATION "${MEDCOUPLING_INSTALL_CMAKE_LOCAL}"
  FILE ${PROJECT_NAME}Targets.cmake)
