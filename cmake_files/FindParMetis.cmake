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

MESSAGE(STATUS "Check for parmetis ...")

SET(PARMETIS_ROOT_DIR $ENV{PARMETIS_ROOT_DIR} CACHE PATH "Path to the PARMETIS.")
IF(PARMETIS_ROOT_DIR)
  LIST(APPEND CMAKE_LIBRARY_PATH "${PARMETIS_ROOT_DIR}")
  LIST(APPEND CMAKE_INCLUDE_PATH "${PARMETIS_ROOT_DIR}/Lib")
ENDIF(PARMETIS_ROOT_DIR)

FIND_LIBRARY(PARMETIS_LIBRARIES parmetis)
FIND_LIBRARY(PARMETIS_SEQ_LIBRARIES metis)
SET(PARMETIS_LIBRARIES ${PARMETIS_LIBRARIES} ${PARMETIS_SEQ_LIBRARIES})
FIND_PATH(PARMETIS_INCLUDE_DIRS parmetis.h)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(ParMetis REQUIRED_VARS PARMETIS_INCLUDE_DIRS PARMETIS_LIBRARIES)
