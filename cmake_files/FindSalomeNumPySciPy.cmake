# Copyright (C) 2013-2016  CEA/DEN, EDF R&D, OPEN CASCADE
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
# Author: Roman NIKOLAEV
#
# Looking for an installation of NumPy and SciPy, and if found the following variables are set
#   NUMPY_INCLUDE_DIR  - NumPy header location
#   NUMPY_DEFINITIONS  - NumPy compiler flags
#   SCIPY_DEFINITIONS  - SciPy compiler flags
#   SCIPY_VERSION      - SciPy version
#

IF(SALOMEPYTHONINTERP_FOUND)   
  # Numpy
  EXECUTE_PROCESS(COMMAND ${PYTHON_EXECUTABLE} -c "import numpy ; import sys ; sys.stdout.write(numpy.get_include())" OUTPUT_VARIABLE NUMPY_INCLUDE_DIR ERROR_QUIET )
  IF(NUMPY_INCLUDE_DIR)
    SET(NUMPY_FOUND TRUE)
  ENDIF(NUMPY_INCLUDE_DIR)
  IF(NUMPY_FOUND)
    SET(NUMPY_DEFINITIONS -DWITH_NUMPY)
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
    SET(SCIPY_DEFINITIONS -DWITH_SCIPY)
  ELSE(SCIPY_FOUND)
    MESSAGE(STATUS "SciPy not found.")
  ENDIF(SCIPY_FOUND)
ENDIF()