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
# Author : Anthony Geay (CEA/DEN)

ADD_DEFINITIONS(${MPI_DEFINITIONS})

INCLUDE_DIRECTORIES(
  ${MPI_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(paramedmem_SOURCES
  ProcessorGroup.cxx
  MPIProcessorGroup.cxx
  ParaMESH.cxx
  ComponentTopology.cxx
  MPIAccess.cxx
  InterpolationMatrix.cxx
  OverlapInterpolationMatrix.cxx
  StructuredCoincidentDEC.cxx
  ExplicitCoincidentDEC.cxx
  InterpKernelDEC.cxx
  ElementLocator.cxx
  OverlapElementLocator.cxx
  MPIAccessDEC.cxx
  TimeInterpolator.cxx
  LinearTimeInterpolator.cxx
  DEC.cxx
  DisjointDEC.cxx
  OverlapDEC.cxx
  ExplicitTopology.cxx
  MxN_Mapping.cxx
  OverlapMapping.cxx
  ParaFIELD.cxx
  ParaGRID.cxx
  BlockTopology.cxx
  ExplicitMapping.cxx
  )

ADD_LIBRARY(paramedmem ${paramedmem_SOURCES})
TARGET_LINK_LIBRARIES(paramedmem medcoupling ${MPI_LIBRARIES})
INSTALL(TARGETS paramedmem EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

FILE(GLOB paramedmem_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${paramedmem_HEADERS_HXX} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# To allow usage as SWIG dependencies:
SET(paramedmem_HEADERS_HXX PARENT_SCOPE)
