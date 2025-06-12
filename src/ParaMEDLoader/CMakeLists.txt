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
# Author : Anthony Geay (CEA/DEN)

add_definitions(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS} ${MPI_DEFINITIONS})

if(HDF5_ENABLE_PARALLEL OR HDF5_IS_PARALLEL)
  add_definitions("-DHDF5_IS_PARALLEL")
endif(HDF5_ENABLE_PARALLEL OR HDF5_IS_PARALLEL)

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

include_directories(
  ${MPI_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../ParaMEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D)

set(paramedloader_SOURCES ParaMEDLoader.cxx ParaMEDFileMesh.cxx)

add_library(paramedloader ${paramedloader_SOURCES})
target_link_libraries(paramedloader paramedmem medloader ${PLATFORM_LIBS})
install(
  TARGETS paramedloader
  EXPORT ${PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

file(GLOB paramedloader_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
install(FILES ${paramedloader_HEADERS_HXX}
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
