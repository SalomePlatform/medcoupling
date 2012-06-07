# Copyright (C) 2007-2012  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
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

INCLUDE_DIRECTORIES(
  ${BOOST_INCLUDE_DIRS}
  ${METIS_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(renumbercpp_HEADERS_HXX
  RENUMBER_Renumbering.hxx
  RenumberingFactory.hxx
  RENUMBERDefines.hxx
  )

SET(renumbercpp_SOURCES
  RENUMBER_Renumbering.cxx
  RenumberingFactory.cxx
  )

SET(renumber_SOURCES
  renumbering.cxx
  )

SET(renumbercpp_LDFLAGS medmem)
SET(renumber_DEFINITIONS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${BOOST_DEFINITIONS}")

IF(MED_ENABLE_METIS)
  SET(renumbercpp_SOURCES ${renumbercpp_SOURCES} RENUMBER_METISRenumbering.cxx)
  SET(renumbercpp_HEADERS_HXX ${renumbercpp_HEADERS_HXX} RENUMBER_METISRenumbering.hxx)
  SET(renumbercpp_LDFLAGS ${renumbercpp_LDFLAGS} ${METIS_LIBS})
  SET(renumber_DEFINITIONS "${renumber_DEFINITIONS} ${METIS_DEFINITIONS}")
ENDIF(MED_ENABLE_METIS)

IF(BOOST_STATUS)
  SET(renumbercpp_SOURCES ${renumbercpp_SOURCES} RENUMBER_BOOSTRenumbering.cxx)
  SET(renumbercpp_HEADERS_HXX ${renumbercpp_HEADERS_HXX} RENUMBER_BOOSTRenumbering.hxx)
  SET(renumbercpp_LDFLAGS ${renumbercpp_LDFLAGS} ${BOOST_LIBS})
  SET(renumber_DEFINITIONS "${renumber_DEFINITIONS} ${BOOST_DEFINITIONS} -DENABLE_BOOST")
ENDIF(BOOST_STATUS)

ADD_EXECUTABLE(renumber ${renumber_SOURCES})
SET_TARGET_PROPERTIES(renumber PROPERTIES COMPILE_FLAGS "${renumber_DEFINITIONS}")
TARGET_LINK_LIBRARIES(renumber renumbercpp)
INSTALL(TARGETS renumber DESTINATION ${MED_salomebin_BINS})

ADD_LIBRARY(renumbercpp SHARED ${renumbercpp_SOURCES})
SET_TARGET_PROPERTIES(renumbercpp PROPERTIES COMPILE_FLAGS "${renumber_DEFINITIONS}")
TARGET_LINK_LIBRARIES(renumbercpp ${renumbercpp_LDFLAGS})
INSTALL(TARGETS renumbercpp DESTINATION ${MED_salomelib_LIBS})

INSTALL(FILES ${renumbercpp_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
