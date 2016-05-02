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

MESSAGE(STATUS "Check for XDR ...")

FIND_PATH(XDR_INCLUDE_DIRS rpc/xdr.h)
IF(XDR_INCLUDE_DIRS)
   SET(XDR_DEFINITIONS "-DHAS_XDR")
ENDIF()

IF(WIN32)
  FIND_LIBRARY(XDR_LIBRARIES xdr)                 # To get the .lib file from XDR
  FIND_PATH(XDR_INCLUDE_DIRS2 stdint.h PATH_SUFFIXES src/msvc)  # To get the stdint.h from XDR (needed by types.h)
  IF(XDR_INCLUDE_DIRS)
    IF(XDR_INCLUDE_DIRS2)
      LIST(APPEND XDR_INCLUDE_DIRS "${XDR_INCLUDE_DIRS2}")
    ELSE()
      SET(XDR_INCLUDE_DIRS "${XDR_INCLUDE_DIRS2}")  # Make the detection fail
    ENDIF()
  ENDIF()
ENDIF(WIN32)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XDR REQUIRED_VARS XDR_INCLUDE_DIRS)
