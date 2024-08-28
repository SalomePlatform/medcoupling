# Copyright (C) 2012-2024  CEA, EDF
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

# add_definitions(${MPI_DEFINITIONS})

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

# include_directories(
#   # ${MPI_INCLUDE_DIRS}
#   ${CMAKE_CURRENT_SOURCE_DIR}
#   ${CMAKE_CURRENT_SOURCE_DIR}/MPIAccess
#   ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
#   ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
#   ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
#   ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
#   ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
#   ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
#   ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints)

set(paramedmem_SOURCES
  BlockTopology.cxx
  CommInterface.cxx
  ComponentTopology.cxx
  DEC.cxx
  DisjointDEC.cxx
  ElementLocator.cxx
  ExplicitCoincidentDEC.cxx
  ExplicitMapping.cxx
  ExplicitTopology.cxx
  InterpKernelDEC.cxx
  InterpolationMatrix.cxx
  LinearTimeInterpolator.cxx
  MPIProcessorGroup.cxx
  ByStringMPIProcessorGroup.cxx
  MxN_Mapping.cxx
  OverlapDEC.cxx
  OverlapElementLocator.cxx
  OverlapInterpolationMatrix.cxx
  OverlapMapping.cxx
  ParaDataArray.cxx
  ParaFIELD.cxx
  ParaGRID.cxx
  ParaMESH.cxx
  ParaSkyLineArray.cxx
  ParaUMesh.cxx
  ProcessorGroup.cxx
  StructuredCoincidentDEC.cxx
  TimeInterpolator.cxx
  MPIAccess.cxx
  MPIAccessDEC.cxx)

add_library(paramedmem ${paramedmem_SOURCES})
target_link_libraries(paramedmem medcouplingcpp medicoco MPI::MPI_CXX)
target_include_directories(paramedmem PUBLIC "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/para_medmem>")

install(
  TARGETS paramedmem
  EXPORT ${CMAKE_PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

file(GLOB paramedmem_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
file(GLOB mpiaccess_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/MPIAccess/*.hxx")
list(APPEND paramedmem_HEADERS_HXX ${mpiaccess_HEADERS_HXX})
install(FILES ${paramedmem_HEADERS_HXX}
  DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

file(GLOB paramedmem_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
install(FILES ${paramedmem_HEADERS_TXX}
  DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# To allow usage as SWIG dependencies:
# set(paramedmem_HEADERS_HXX PARENT_SCOPE)
