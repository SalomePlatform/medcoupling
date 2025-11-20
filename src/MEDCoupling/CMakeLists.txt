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

if(MEDCOUPLING_BUILD_TESTS)
  add_subdirectory(Test)
endif(MEDCOUPLING_BUILD_TESTS)

# KILLER WARNINGS TO HELP INT64 MANAGEMENT ADD_DEFINITIONS(-Wextra -Wsynth
# -Wno-unused-parameter -pedantic -Wall -Wshadow  -Wpointer-arith -Wcast-qual
# -Wcast-align -Wwrite-strings -Wctor-dtor-privacy -Wnon-virtual-dtor
# -Woverloaded-virtual -Wreturn-type  -Wno-endif-labels  -Wsign-compare
# -Wmissing-format-attribute -Wno-multichar -Wno-deprecated-declarations
# -Wpacked -Wredundant-decls -Wlong-long -Wdisabled-optimization
# -Wunknown-pragmas -Wundef  -Wreorder)

if(NOT DEFINED MSVC)
  add_definitions(-Wsign-compare -Wconversion)
endif()

include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}/../..
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/LinearAlgebra)

set(medcoupling_SOURCES
    MEDCouplingField.cxx
    MEDCouplingFieldFloat.cxx
    MEDCouplingFieldDouble.cxx
    MEDCouplingFieldInt32.cxx
    MEDCouplingFieldInt64.cxx
    MEDCouplingUMesh.cxx
    MEDCouplingUMesh_internal.cxx
    MEDCouplingUMesh_intersection.cxx
    MEDCoupling1GTUMesh.cxx
    MEDCouplingMemArray.cxx
    MEDCouplingMemArrayFloat.cxx
    MEDCouplingMemArrayChar.cxx
    MEDCouplingMap.cxx
    MEDCouplingTraits.cxx
    MEDCouplingTimeLabel.cxx
    MEDCouplingCMesh.cxx
    MEDCouplingIMesh.cxx
    MEDCouplingCurveLinearMesh.cxx
    MEDCouplingStructuredMesh.cxx
    MEDCouplingTimeDiscretization.cxx
    MEDCouplingFieldDiscretization.cxx
    MEDCouplingFieldDiscretizationOnNodesFE.cxx
    MEDCouplingRefCountObject.cxx
    MEDCouplingPointSet.cxx
    MEDCouplingFieldTemplate.cxx
    MEDCouplingMappedExtrudedMesh.cxx
    MEDCouplingMesh.cxx
    MEDCouplingGaussLocalization.cxx
    MEDCouplingNatureOfField.cxx
    MEDCouplingMultiFields.cxx
    MEDCouplingDefinitionTime.cxx
    MEDCouplingFieldOverTime.cxx
    MEDCouplingCartesianAMRMesh.cxx
    MEDCouplingAMRAttribute.cxx
    MEDCouplingMatrix.cxx
    MEDCouplingPartDefinition.cxx
    MEDCouplingSkyLineArray.cxx
    MEDCouplingVoronoi.cxx)

set(medcouplingremapper_SOURCES MEDCouplingRemapper.cxx)

add_library(medcouplingcpp ${medcoupling_SOURCES})
set_target_properties(medcouplingcpp PROPERTIES OUTPUT_NAME "medcoupling")
target_link_libraries(medcouplingcpp interpkernel)
install(
  TARGETS medcouplingcpp
  EXPORT ${PROJECT_NAME}TargetGroup
  ARCHIVE DESTINATION ${MEDCOUPLING_INSTALL_LIBS}
  RUNTIME DESTINATION ${MEDCOUPLING_INSTALL_BINS}
  LIBRARY DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

add_library(medcouplingremapper ${medcouplingremapper_SOURCES})
target_link_libraries(medcouplingremapper medcouplingcpp)
install(
  TARGETS medcouplingremapper
  EXPORT ${PROJECT_NAME}TargetGroup
  ARCHIVE DESTINATION ${MEDCOUPLING_INSTALL_LIBS}
  RUNTIME DESTINATION ${MEDCOUPLING_INSTALL_BINS}
  LIBRARY DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

file(GLOB medcoupling_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
file(GLOB medcoupling_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
install(FILES ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX}
              MEDCouplingNatureOfFieldEnum
        DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# To allow usage as SWIG dependencies:
set(medcoupling_HEADERS_HXX PARENT_SCOPE)
set(medcoupling_HEADERS_TXX PARENT_SCOPE)
