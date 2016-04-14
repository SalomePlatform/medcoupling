# Copyright (C) 2007-2016  CEA/DEN, EDF R&D, OPEN CASCADE
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

# ------

MESSAGE(STATUS "Check for metis ...")

SET(METIS_ROOT_DIR $ENV{METIS_ROOT_DIR} CACHE PATH "Path to the METIS.")
IF(METIS_ROOT_DIR)
  LIST(APPEND CMAKE_LIBRARY_PATH "${METIS_ROOT_DIR}")
  LIST(APPEND CMAKE_LIBRARY_PATH "${METIS_ROOT_DIR}/lib")
  LIST(APPEND CMAKE_INCLUDE_PATH "${METIS_ROOT_DIR}/Lib")
  LIST(APPEND CMAKE_INCLUDE_PATH "${METIS_ROOT_DIR}/include")
ENDIF(METIS_ROOT_DIR)

FIND_LIBRARY(METIS_LIBRARIES metis)
FIND_PATH(METIS_INCLUDE_DIRS metis.h)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Metis REQUIRED_VARS METIS_INCLUDE_DIRS METIS_LIBRARIES)
FILE(READ ${METIS_INCLUDE_DIRS}/metis.h metis_h_content)
STRING(REPLACE "\n" ";" list_metis_h_content ${metis_h_content})
FOREACH(ln ${list_metis_h_content})
  IF("${ln}" MATCHES "^#define METIS_VER_MAJOR")
    STRING(REPLACE "#define METIS_VER_MAJOR" "" metis_major_version "${ln}")
    STRING(STRIP "${metis_major_version}" metis_major_version)
  ENDIF("${ln}" MATCHES "^#define METIS_VER_MAJOR")
ENDFOREACH(ln ${list_metis_h_content})
IF(metis_major_version STREQUAL 5)
  SET(MEDCOUPLING_METIS_V5 1)
  MESSAGE(STATUS "Metis maj version 5 detected.")
ELSE(metis_major_version STREQUAL 5)
   MESSAGE(STATUS "Metis maj version 4 detected.")
ENDIF(metis_major_version STREQUAL 5)
