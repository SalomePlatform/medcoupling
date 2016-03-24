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

MESSAGE(STATUS "Check for scotch ...")

SET(SCOTCH_ROOT_DIR $ENV{SCOTCH_ROOT_DIR} CACHE PATH "Path to the SCOTCH.")
IF(SCOTCH_ROOT_DIR)
  LIST(APPEND CMAKE_PREFIX_PATH "${SCOTCH_ROOT_DIR}")
ENDIF(SCOTCH_ROOT_DIR)

FIND_LIBRARY(SCOTCH_LIBRARIES scotch)
FIND_LIBRARY(SCOTCH_ERR_LIBRARIES scotcherr)
SET(SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES} ${SCOTCH_ERR_LIBRARIES})
FIND_PATH(SCOTCH_INCLUDE_DIRS scotch.h PATH_SUFFIXES "/scotch")

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Scotch REQUIRED_VARS SCOTCH_INCLUDE_DIRS SCOTCH_LIBRARIES)
