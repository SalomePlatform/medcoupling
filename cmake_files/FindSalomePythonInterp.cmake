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
# Author: Adrien Bruneton
#

# Python interpreter detection for SALOME
#
#  !! Please read the generic detection procedure in SalomeMacros.cmake !!
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
  GET_FILENAME_COMPONENT(_python_bin "${PYTHON_EXECUTABLE}" NAME )
  SET(PYTHONBIN "${_python_bin}" CACHE STRING "Name of Python interpreter")
  SALOME_ACCUMULATE_ENVIRONMENT(PATH ${PYTHON_EXECUTABLE})
  SALOME_ACCUMULATE_ENVIRONMENT(PYTHONPATH ${PYTHON_PYTHONPATH})
ENDIF()

