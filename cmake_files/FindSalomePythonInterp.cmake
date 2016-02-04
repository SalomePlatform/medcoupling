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

# Python interpreter detection for SALOME
#
#  !! Please read the generic detection procedure in SalomeMacros.cmake !!
#
# We also look for an installation of NumPy, and if found the following variables are set
#   NUMPY_INCLUDE_DIR  - NumPy header location
#   NUMPY_DEFINITIONS  - compiler flag
# and are automatically appended to PYTHON_INCLUDE_DIRS (and PYTHON_DEFINITIONS resp.)    
#

# Make sure the detection of both libs and interpreter (if both needed) occur in the correct order:
IF(SALOMEPYTHONLIBS_FOUND AND NOT SALOMEPYTHONINTERP_FOUND)
   MESSAGE(FATAL_ERROR "Developer error -> Python interpreter should be detected/required before Python libs!")
ENDIF()

# Use the PYTHON_ROOT_DIR if PYTHONINTERP_ROOT_DIR is not defined:
SET(PYTHON_ROOT_DIR "$ENV{PYTHON_ROOT_DIR}" CACHE PATH "Path to the Python installation (libs+interpreter)")
IF(EXISTS "${PYTHON_ROOT_DIR}" AND (NOT PYTHONINTERP_ROOT_DIR))
  # Extract sub-directory "paraview-x.xx":
  MESSAGE(STATUS "Setting PYTHONINTERP_ROOT_DIR to: ${PYTHON_ROOT_DIR}")
  SET(PYTHONINTERP_ROOT_DIR "${PYTHON_ROOT_DIR}" CACHE PATH "Path to PythonInterp directory")
ENDIF()
SALOME_FIND_PACKAGE_AND_DETECT_CONFLICTS(PythonInterp PYTHON_EXECUTABLE 1)

IF(SALOMEPYTHONINTERP_FOUND) 
  SET(PYTHON_PYTHONPATH "${PYTHON_ROOT_DIR}/lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}")
  
  ##
  ## Specifics -- NumPy/SciPy detection
  ##
  
  # Numpy
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

  SALOME_ACCUMULATE_ENVIRONMENT(PATH ${PYTHON_EXECUTABLE})
  SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH ${PYTHON_PYTHONPATH})
ENDIF()

