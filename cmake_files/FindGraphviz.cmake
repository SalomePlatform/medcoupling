# - Graphviz detection
#
# Output variables: GRAPHVIZ_EXECUTABLE   - where is executable 'dot' takes place.
#		    GRAPHVIZ_INCLUDE_DIRS - where to find headers.
# 		    GRAPHVIZ_LIBRARIES    - where to get libraries.
#		    GRAPHVIZ_VERSION      - Graphviz version
#		    GRAPHVIZ_DEFINITIONS  - Graphviz definitions
# 		    GRAPHVIZ_FOUND        - True if Graphviz was found.
#
###########################################################################
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D, OPEN CASCADE
#
# Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
# CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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

FIND_PROGRAM(GRAPHVIZ_EXECUTABLE dot)

FIND_PATH(GRAPHVIZ_INCLUDE_DIR NAMES graphviz/cgraph.h)
SET(GRAPHVIZ_INCLUDE_DIRS ${GRAPHVIZ_INCLUDE_DIR} ${GRAPHVIZ_INCLUDE_DIR}/graphviz)

IF(GRAPHVIZ_EXECUTABLE)
  IF(NOT CMAKE_PREFIX_PATH) # Case where dot is in the path and -DGRAPHVIZ_ROOT_DIR not defined.
    GET_FILENAME_COMPONENT(GRAPHVIZ_EXECUTABLE_DIR ${GRAPHVIZ_EXECUTABLE} PATH)
    GET_FILENAME_COMPONENT(GRAPHVIZ_EXECUTABLE_DIR ${GRAPHVIZ_EXECUTABLE_DIR} PATH)
    SET(CMAKE_PREFIX_PATH "${GRAPHVIZ_EXECUTABLE_DIR}")
  ENDIF(NOT CMAKE_PREFIX_PATH)
ENDIF(GRAPHVIZ_EXECUTABLE)

FIND_LIBRARY(GRAPHVIZ_cdt_LIBRARY      NAMES cdt PATH_SUFFIXES bin)
FIND_LIBRARY(GRAPHVIZ_graph_LIBRARY    NAMES cgraph PATH_SUFFIXES bin)
FIND_LIBRARY(GRAPHVIZ_gvc_LIBRARY      NAMES gvc PATH_SUFFIXES bin)
FIND_LIBRARY(GRAPHVIZ_pathplan_LIBRARY NAMES pathplan PATH_SUFFIXES bin)

SET(GRAPHVIZ_LIBRARIES 
  ${GRAPHVIZ_cdt_LIBRARY}
  ${GRAPHVIZ_graph_LIBRARY} 
  ${GRAPHVIZ_gvc_LIBRARY}
  ${GRAPHVIZ_pathplan_LIBRARY}
  )

IF(GRAPHVIZ_EXECUTABLE)
  EXECUTE_PROCESS(COMMAND ${GRAPHVIZ_EXECUTABLE} "-V" ERROR_VARIABLE GRAPHVIZ_VERSION ERROR_STRIP_TRAILING_WHITESPACE)
  STRING(REGEX REPLACE ".* ([0-9.]+) .*" "\\1" GRAPHVIZ_VERSION "${GRAPHVIZ_VERSION}")
ENDIF()

SET(GRAPHVIZ_DEFINITIONS)
IF("${GRAPHVIZ_VERSION}" VERSION_LESS "2.36.0")
  SET(GRAPHVIZ_DEFINITIONS -DWITH_CGRAPH)
ENDIF()

## Don't detect cgraph on Windows
#IF(NOT WIN32)
#  FIND_LIBRARY(GRAPHVIZ_cgraph_LIBRARY   NAMES cgraph PATH_SUFFIXES bin)
#  SET(GRAPHVIZ_LIBRARIES ${GRAPHVIZ_cgraph_LIBRARY})
#ENDIF()

# Handle the standard arguments of the find_package() command:
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Graphviz REQUIRED_VARS 
					   GRAPHVIZ_EXECUTABLE 
					   GRAPHVIZ_LIBRARIES 
					   GRAPHVIZ_INCLUDE_DIRS)
