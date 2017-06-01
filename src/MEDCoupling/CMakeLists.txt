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

IF(MEDCOUPLING_BUILD_TESTS)
  ADD_SUBDIRECTORY(Test)
ENDIF(MEDCOUPLING_BUILD_TESTS)

### KILLER WARNINGS TO HELP INT64 MANAGEMENT
#ADD_DEFINITIONS(-Wextra -Wsynth -Wno-unused-parameter -pedantic -Wall -Wshadow  -Wpointer-arith -Wcast-qual -Wcast-align -Wwrite-strings -Wctor-dtor-privacy -Wnon-virtual-dtor -Woverloaded-virtual -Wreturn-type  -Wno-endif-labels  -Wsign-compare -Wmissing-format-attribute -Wno-multichar -Wno-deprecated-declarations -Wpacked -Wredundant-decls -Wlong-long -Wdisabled-optimization -Wunknown-pragmas -Wundef  -Wreorder)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}/../..
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(medcoupling_SOURCES
  MEDCouplingField.cxx
  MEDCouplingFieldFloat.cxx
  MEDCouplingFieldDouble.cxx
  MEDCouplingFieldInt.cxx
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
  MEDCouplingVoronoi.cxx
  ../ICoCo/ICoCoField.cxx
  ../ICoCo/ICoCoMEDField.cxx
  )

SET(medcouplingremapper_SOURCES
  MEDCouplingRemapper.cxx
  )

ADD_LIBRARY(medcoupling ${medcoupling_SOURCES})
TARGET_LINK_LIBRARIES(medcoupling interpkernel)
INSTALL(TARGETS medcoupling EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

ADD_LIBRARY(medcouplingremapper ${medcouplingremapper_SOURCES})
TARGET_LINK_LIBRARIES(medcouplingremapper medcoupling)
INSTALL(TARGETS medcouplingremapper EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

FILE(GLOB medcoupling_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
FILE(GLOB medcoupling_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
FILE(GLOB icoco_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/../ICoCo/*.hxx")
INSTALL(FILES ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX} ${icoco_HEADERS_HXX} MEDCouplingNatureOfFieldEnum DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# To allow usage as SWIG dependencies:
SET(medcoupling_HEADERS_HXX PARENT_SCOPE)
SET(medcoupling_HEADERS_TXX PARENT_SCOPE)
