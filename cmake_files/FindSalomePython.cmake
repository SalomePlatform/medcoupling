# Copyright (C) 2013-2015  CEA/DEN, EDF R&D, OPEN CASCADE
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
# Author: Adrien Bruneton
#

# Python libraries and interpreter detection for SALOME
#
#  !! Please read the generic detection procedure in SalomeMacros.cmake !!
#
# The interpreter is found first, and if OK, the corresponding libraries are searched.
# We ensure the version of the libraries matches the one of the interpreter.
#
# We also look for an installation of NumPy, and if found the following variables are set
#   NUMPY_INCLUDE_DIR  - NumPy header location
#   NUMPY_DEFINITIONS  - compiler flag
# and are automatically appended to PYTHON_INCLUDE_DIRS (and PYTHON_DEFINITIONS resp.)    
#

# 1. Load environment or any previously detected Python
IF(DEFINED ENV{PYTHON_ROOT_DIR})
  FILE(TO_CMAKE_PATH "$ENV{PYTHON_ROOT_DIR}" _PYTHON_ROOT_DIR_ENV)
  SET(_dflt_value "${_PYTHON_ROOT_DIR_ENV}")
ELSE()
  # will be blank if no Python was previously loaded
  SET(_dflt_value "${PYTHON_ROOT_DIR_EXP}")
ENDIF()

#   Make cache entry 
SET(PYTHON_ROOT_DIR "${_dflt_value}" CACHE PATH "Path to Python directory (interpreter and libs)")

# 2. Find package - config mode first (i.e. looking for XYZ-config.cmake)
IF(WIN32)
 IF(CMAKE_BUILD_TYPE STREQUAL Debug)
  SET(PythonInterp_FIND_VERSION _d)
  SET(PYTHON_DEFINITIONS "-DHAVE_DEBUG_PYTHON")
 ENDIF(CMAKE_BUILD_TYPE STREQUAL Debug)
ENDIF(WIN32)
IF(EXISTS "${PYTHON_ROOT_DIR}")
  # Hope to find direclty a CMake config file there
  SET(_CONF_DIR "${PYTHON_ROOT_DIR}/share/cmake") 

  # Try find_package in config mode with a hard-coded guess. This
  # has the priority.
  FIND_PACKAGE(Python CONFIG QUIET PATHS "${_CONF_DIR}")
  MARK_AS_ADVANCED(Python_DIR)
    
  IF (NOT PYTHON_FOUND)  
    LIST(APPEND CMAKE_PREFIX_PATH "${PYTHON_ROOT_DIR}")
  ELSE()
    MESSAGE(STATUS "Found Python in CONFIG mode!")
  ENDIF()
ENDIF()

# Otherwise try the standard way (module mode, with the standard CMake Find*** macro):
SALOME_FIND_PACKAGE(SalomePython PythonInterp MODULE)
SET(_found1 ${PYTHONINTERP_FOUND})

IF (PYTHONINTERP_FOUND)
  # Now ensure we find the Python libraries matching the interpreter:
  # This uses the variable PYTHON_EXECUTABLE

  GET_FILENAME_COMPONENT(_python_bin "${PYTHON_EXECUTABLE}" NAME )
  SET(PYTHONBIN "${_python_bin}" CACHE STRING "Name of Python interpreter")

  GET_FILENAME_COMPONENT(_python_dir "${PYTHON_EXECUTABLE}" PATH)
  GET_FILENAME_COMPONENT(CMAKE_INCLUDE_PATH "${_python_dir}/../include/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}" ABSOLUTE)
  GET_FILENAME_COMPONENT(CMAKE_LIBRARY_PATH "${_python_dir}/../lib" ABSOLUTE)
  # For a Windows install, this might look more like this:
  IF(WIN32)
    LIST(APPEND CMAKE_LIBRARY_PATH "${_python_dir}/libs" ABSOLUTE)
    LIST(APPEND CMAKE_INCLUDE_PATH "${_python_dir}/include" ABSOLUTE)
  ENDIF()
  # Override the EXACT and VERSION settings of the SalomePython module
  # to force the next call to SALOME_FIND_PACKAGE() to find the exact matching
  # version:
  SET(_old_EXACT ${SalomePython_FIND_VERSION_EXACT})
  SET(_old_VERSION "${SalomePython_FIND_VERSION}")
  SET(SalomePython_FIND_VERSION_EXACT TRUE)
  SET(SalomePython_FIND_VERSION "${PYTHON_VERSION_STRING}")
  # Prepare call to FIND_PACKAGE(PythonLibs) and ensure priority is given to 
  # the location found for the interpreter:
  GET_FILENAME_COMPONENT(_tmp "${_python_dir}" PATH)
#  SET(PYTHON_LIBRARY ${_tmp}/lib)
#  SET(PYTHON_INCLUDE_DIR ${_tmp}/include)
  SALOME_FIND_PACKAGE(SalomePython PythonLibs MODULE)
  # Restore variables:
  SET(SalomePython_FIND_VERSION_EXACT ${_old_EXACT})
  SET(SalomePython_FIND_VERSION "${_old_VERSION}")
ENDIF()

# Set the FOUND flag for SalomePython and Python:
SET(SALOMEPYTHON_FOUND FALSE)
IF (_found1 AND PYTHONLIBS_FOUND)
  
  # 24.03.2015 ANA: Fix problem on Windows in  Debug mode
  # If you have Python, installed by Windows MSI Installer, 
  # PYTHON_LIBRARIES variable contains redundant release libraries...
  IF(WIN32 AND CMAKE_BUILD_TYPE STREQUAL Debug)
    SET (PYTHON_LIBRARIES ${PYTHON_DEBUG_LIBRARIES})
  ENDIF()

  SET(SALOMEPYTHON_FOUND TRUE)
  SET(Python_FOUND TRUE)
ELSE()
SET(SALOMEPYTHON_FOUND FALSE)
  SET(Python_FOUND FALSE)
ENDIF()

IF (SALOMEPYTHON_FOUND)
  MESSAGE(STATUS "Python interpreter and Python libraries found:")
  MESSAGE(STATUS "Python libraries: ${PYTHON_LIBRARIES}")
  MESSAGE(STATUS "Python include dir: ${PYTHON_INCLUDE_DIR}")

  # 3. Set the root dir which was finally retained 
  # For Python this is the grand-parent of the
  # include directory:
  GET_FILENAME_COMPONENT(_tmp_ROOT_DIR "${PYTHON_INCLUDE_DIR}" PATH)
  IF(NOT WIN32)
    GET_FILENAME_COMPONENT(_tmp_ROOT_DIR "${_tmp_ROOT_DIR}" PATH)
  ENDIF()

  # 4. Warn if CMake found something not located under ENV(XYZ_ROOT_DIR)
  IF(DEFINED ENV{PYTHON_ROOT_DIR})
    SALOME_CHECK_EQUAL_PATHS(_res "${_tmp_ROOT_DIR}" "${_PYTHON_ROOT_DIR_ENV}")
    IF(NOT _res)
      MESSAGE(WARNING "Python was found, but not a the path given by the "
"environment PYTHON_ROOT_DIR! Is the variable correctly set?")
    ELSE()
      MESSAGE(STATUS "Python found directory matches what was specified in the PYTHON_ROOT_DIR, all good!")    
    ENDIF()
  ENDIF()

  # 5. Conflict detection
  # 5.1  From another prerequisite using Python
  IF(PYTHON_ROOT_DIR_EXP)
      SALOME_CHECK_EQUAL_PATHS(_res "${_tmp_ROOT_DIR}" "${PYTHON_ROOT_DIR_EXP}") 
      IF(NOT _res)
         MESSAGE(WARNING "Warning: Python: detected version conflicts with a previously found Python!"
                          "The two paths are " ${_tmp_ROOT_DIR} " vs " ${PYTHON_ROOT_DIR_EXP})
      ELSE()
          MESSAGE(STATUS "Python directory matches what was previously exposed by another prereq, all good!")
      ENDIF()        
  ENDIF()

  ##
  ## 6. Save the final detected installation
  ##
  SET(PYTHON_ROOT_DIR "${_tmp_ROOT_DIR}")
  SET(PYTHON_PYTHONPATH "${_tmp_ROOT_DIR}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")

  ## 7. Specifics
  ##
  
  # NumPy detection 
  EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy ; import sys ; sys.stdout.write(numpy.get_include())" OUTPUT_VARIABLE NUMPY_INCLUDE_DIR ERROR_QUIET )
  IF(NUMPY_INCLUDE_DIR)
    SET(NUMPY_FOUND TRUE)
  ENDIF(NUMPY_INCLUDE_DIR)
  IF(NUMPY_FOUND)
    SET(PYTHON_INCLUDE_DIRS ${NUMPY_INCLUDE_DIR} ${PYTHON_INCLUDE_DIRS})
    SET(PYTHON_DEFINITIONS "${PYTHON_DEFINITIONS} -DWITH_NUMPY")
    MESSAGE(STATUS "NumPy found : ${NUMPY_INCLUDE_DIR}")
  ELSE(NUMPY_FOUND)
    MESSAGE(STATUS "NumPy not found.")
  ENDIF(NUMPY_FOUND)
  # SciPy detection
  EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import scipy ; import sys ; sys.stdout.write(scipy.version.version)" OUTPUT_VARIABLE SCIPY_VERSION ERROR_QUIET )
  IF(SCIPY_VERSION)
    SET(SCIPY_FOUND TRUE)
  ENDIF(SCIPY_VERSION)
  IF(SCIPY_FOUND)
    MESSAGE(STATUS "Scipy found : Version ${SCIPY_VERSION}")
  ENDIF(SCIPY_FOUND)
  ## None here    
ELSE()
  MESSAGE(STATUS "Python was only partially (or not at all) found .")
ENDIF()

IF(SALOMEPYTHON_FOUND) 
  SALOME_ACCUMULATE_HEADERS(PYTHON_INCLUDE_DIR)
  SALOME_ACCUMULATE_ENVIRONMENT(PATH ${PYTHON_EXECUTABLE})
  SALOME_ACCUMULATE_ENVIRONMENT(LD_LIBRARY_PATH ${PYTHON_LIBRARIES})
  SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH ${PYTHON_PYTHONPATH})
ENDIF()
