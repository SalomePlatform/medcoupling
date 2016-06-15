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

# Python libraries detection for SALOME
#
#  !! Please read the generic detection procedure in SalomeMacros.cmake !!
#

# Use the PYTHON_ROOT_DIR if PYTHONLIBS_ROOT_DIR is not defined:
SET(PYTHON_ROOT_DIR "$ENV{PYTHON_ROOT_DIR}" CACHE PATH "Path to the Python installation (libs+interpreter)")
IF(EXISTS "${PYTHON_ROOT_DIR}" AND (NOT PYTHONLIBS_ROOT_DIR))
  MESSAGE(STATUS "Setting PYTHONLIBS_ROOT_DIR to: ${PYTHON_ROOT_DIR}")
  SET(PYTHONLIBS_ROOT_DIR "${PYTHON_ROOT_DIR}" CACHE PATH "Path to PythonLibs directory")
ENDIF()
IF (SALOMEPYTHONINTERP_FOUND AND NOT "${PYTHON_VERSION_STRING}" STREQUAL "") 
   # Trying to search libraries with same version as an interpreter version
   SET(PythonLibs_FIND_VERSION ${PYTHON_VERSION_STRING})
   SET(PythonLibs_FIND_VERSION_MAJOR ${PYTHON_VERSION_MAJOR})
ENDIF()
IF(WIN32)
  set(CMAKE_LIBRARY_PATH "${PYTHON_ROOT_DIR}/libs")
ENDIF()
SALOME_FIND_PACKAGE_AND_DETECT_CONFLICTS(PythonLibs PYTHON_INCLUDE_DIR 2)

IF(SALOMEPYTHONLIBS_FOUND) 
  SALOME_ACCUMULATE_HEADERS(PYTHON_INCLUDE_DIR)
  SALOME_ACCUMULATE_ENVIRONMENT(LD_LIBRARY_PATH ${PYTHON_LIBRARIES})
ENDIF()

## Specifics -- check matching version with Interpreter if already detected:
IF (SALOMEPYTHONLIBS_FOUND AND SALOMEPYTHONINTERP_FOUND)
  # Now ensure versions are matching
  IF("${PYTHONLIBS_VERSION_STRING}" STREQUAL "${PYTHON_VERSION_STRING}")
    MESSAGE(STATUS "Python libs and interpreter versions are matching: ${PYTHONLIBS_VERSION_STRING}")
  ELSE()
    MESSAGE(FATAL_ERROR "Python libs and interpreter versions are NOT matching: ${PYTHONLIBS_VERSION_STRING} vs ${PYTHON_VERSION_STRING}")
  ENDIF()
ENDIF()
