# Copyright (C) 2012-2025  CEA, EDF
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this library; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email :
# webmaster.salome@opencascade.com
#

add_definitions(${HDF5_DEFINITIONS} ${BOOST_DEFINITIONS} ${MEDFILE_DEFINITIONS})

if(Boost_FOUND)
  add_definitions("-DENABLE_BOOST")
endif(Boost_FOUND)

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

# Include directories
include_directories(
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints)

if(MEDCOUPLING_USE_MPI)
  add_definitions(${MPI_DEFINITIONS})
  include_directories(${MPI_INCLUDE_DIRS})
endif(MEDCOUPLING_USE_MPI)

set(renumbercpp_HEADERS_HXX RENUMBER_Renumbering.hxx RenumberingFactory.hxx
                            RENUMBERDefines.hxx)

set(renumbercpp_SOURCES RENUMBER_Renumbering.cxx RenumberingFactory.cxx)

set(renumber_SOURCES renumbering.cxx)

set(renumbercpp_LDFLAGS medloader)

if(PARMETIS_FOUND)
  include_directories(${PARMETIS_INCLUDE_DIRS})
endif(PARMETIS_FOUND)

if(METIS_FOUND)
  include_directories(${METIS_INCLUDE_DIRS})

  set(renumbercpp_SOURCES ${renumbercpp_SOURCES} RENUMBER_METISRenumbering.cxx)
  set(renumbercpp_HEADERS_HXX ${renumbercpp_HEADERS_HXX}
                              RENUMBER_METISRenumbering.hxx)
  set(renumbercpp_LDFLAGS ${renumbercpp_LDFLAGS} ${METIS_LIBRARIES})
  set(renumber_DEFINITIONS "${renumber_DEFINITIONS} ${METIS_DEFINITIONS}")
endif(METIS_FOUND)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})

  set(renumbercpp_SOURCES ${renumbercpp_SOURCES} RENUMBER_BOOSTRenumbering.cxx)
  set(renumbercpp_HEADERS_HXX ${renumbercpp_HEADERS_HXX}
                              RENUMBER_BOOSTRenumbering.hxx)
  set(renumbercpp_LDFLAGS ${renumbercpp_LDFLAGS} ${BOOST_LIBRARIES})
endif(Boost_FOUND)

if(MEDCOUPLING_USE_MPI)
  set(renumbercpp_LDFLAGS ${renumbercpp_LDFLAGS} ${MPI_LIBRARIES})
endif(MEDCOUPLING_USE_MPI)

add_executable(renumber ${renumber_SOURCES})
set_target_properties(renumber PROPERTIES COMPILE_FLAGS
                                          "${renumber_DEFINITIONS}")
target_link_libraries(renumber renumbercpp)
install(TARGETS renumber DESTINATION ${MEDCOUPLING_INSTALL_BINS})

add_library(renumbercpp ${renumbercpp_SOURCES})
set_target_properties(renumbercpp PROPERTIES COMPILE_FLAGS
                                             "${renumber_DEFINITIONS}")
target_link_libraries(renumbercpp ${renumbercpp_LDFLAGS})
install(TARGETS renumbercpp DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

install(FILES ${renumbercpp_HEADERS_HXX}
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
