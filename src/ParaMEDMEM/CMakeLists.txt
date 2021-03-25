# Copyright (C) 2012-2021  CEA/DEN, EDF R&D
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
# Author : Anthony Geay (CEA/DEN)

ADD_DEFINITIONS(${MPI_DEFINITIONS})

IF (NOT DEFINED MSVC)
  ADD_DEFINITIONS(-Wsign-compare -Wconversion)
ENDIF()

INCLUDE_DIRECTORIES(
  ${MPI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/MPIAccess
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(paramedmem_SOURCES
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
  MPIAccess/MPIAccess.cxx
  MPIAccess/MPIAccessDEC.cxx
)

ADD_LIBRARY(paramedmem ${paramedmem_SOURCES})
TARGET_LINK_LIBRARIES(paramedmem medcouplingcpp ${MPI_LIBRARIES})
INSTALL(TARGETS paramedmem EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

FILE(GLOB paramedmem_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
FILE(GLOB mpiaccess_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/MPIAccess/*.hxx")
LIST(APPEND paramedmem_HEADERS_HXX ${mpiaccess_HEADERS_HXX})
INSTALL(FILES ${paramedmem_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

FILE(GLOB paramedmem_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${paramedmem_HEADERS_TXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# To allow usage as SWIG dependencies:
SET(paramedmem_HEADERS_HXX PARENT_SCOPE)
