# Copyright (C) 2012-2016  CEA/DEN, EDF R&D
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

ADD_DEFINITIONS(${HDF5_DEFINITIONS} ${BOOST_DEFINITIONS} ${MEDFILE_DEFINITIONS})

IF(Boost_FOUND)
  ADD_DEFINITIONS("-DENABLE_BOOST")
ENDIF(Boost_FOUND)

# Include directories
INCLUDE_DIRECTORIES(
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

IF(MEDCOUPLING_USE_MPI)
  ADD_DEFINITIONS(${MPI_DEFINITIONS})
  INCLUDE_DIRECTORIES(${MPI_INCLUDE_DIRS})
ENDIF(MEDCOUPLING_USE_MPI)

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

SET(renumbercpp_LDFLAGS medloader)

IF(PARMETIS_FOUND)
  INCLUDE_DIRECTORIES(${PARMETIS_INCLUDE_DIRS})
ENDIF(PARMETIS_FOUND)

IF(METIS_FOUND)
  INCLUDE_DIRECTORIES(${METIS_INCLUDE_DIRS})

  SET(renumbercpp_SOURCES ${renumbercpp_SOURCES} RENUMBER_METISRenumbering.cxx)
  SET(renumbercpp_HEADERS_HXX ${renumbercpp_HEADERS_HXX} RENUMBER_METISRenumbering.hxx)
  SET(renumbercpp_LDFLAGS ${renumbercpp_LDFLAGS} ${METIS_LIBRARIES})
  SET(renumber_DEFINITIONS "${renumber_DEFINITIONS} ${METIS_DEFINITIONS}")
ENDIF(METIS_FOUND)

IF(Boost_FOUND)
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

  SET(renumbercpp_SOURCES ${renumbercpp_SOURCES} RENUMBER_BOOSTRenumbering.cxx)
  SET(renumbercpp_HEADERS_HXX ${renumbercpp_HEADERS_HXX} RENUMBER_BOOSTRenumbering.hxx)
  SET(renumbercpp_LDFLAGS ${renumbercpp_LDFLAGS} ${BOOST_LIBRARIES})
ENDIF(Boost_FOUND)

IF(MEDCOUPLING_USE_MPI)
  SET(renumbercpp_LDFLAGS ${renumbercpp_LDFLAGS} ${MPI_LIBRARIES})
ENDIF(MEDCOUPLING_USE_MPI)

ADD_EXECUTABLE(renumber ${renumber_SOURCES})
SET_TARGET_PROPERTIES(renumber PROPERTIES COMPILE_FLAGS "${renumber_DEFINITIONS}")
TARGET_LINK_LIBRARIES(renumber renumbercpp)
INSTALL(TARGETS renumber DESTINATION ${MEDCOUPLING_INSTALL_BINS})

ADD_LIBRARY(renumbercpp ${renumbercpp_SOURCES})
SET_TARGET_PROPERTIES(renumbercpp PROPERTIES COMPILE_FLAGS "${renumber_DEFINITIONS}")
TARGET_LINK_LIBRARIES(renumbercpp ${renumbercpp_LDFLAGS})
INSTALL(TARGETS renumbercpp DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

INSTALL(FILES ${renumbercpp_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
