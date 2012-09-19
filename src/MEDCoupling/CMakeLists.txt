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
# Author : Anthony Geay (CEA/DEN)

ADD_SUBDIRECTORY(Test)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(medcoupling_SOURCES
  MEDCouplingField.cxx
  MEDCouplingFieldDouble.cxx
  MEDCouplingUMesh.cxx
  MEDCouplingMemArray.cxx
  MEDCouplingTimeLabel.cxx
  MEDCouplingCMesh.cxx
  MEDCouplingTimeDiscretization.cxx
  MEDCouplingFieldDiscretization.cxx
  MEDCouplingRefCountObject.cxx
  MEDCouplingPointSet.cxx
  MEDCouplingUMeshDesc.cxx
  MEDCouplingFieldTemplate.cxx
  MEDCouplingExtrudedMesh.cxx
  MEDCouplingMesh.cxx
  MEDCouplingGaussLocalization.cxx
  MEDCouplingNatureOfField.cxx
  MEDCouplingMultiFields.cxx
  MEDCouplingDefinitionTime.cxx
  MEDCouplingFieldOverTime.cxx
  )

SET(medcouplingremapper_SOURCES
  MEDCouplingRemapper.cxx
  )

ADD_LIBRARY(medcoupling SHARED ${medcoupling_SOURCES})
SET_TARGET_PROPERTIES(medcoupling PROPERTIES COMPILE_FLAGS "${PLATFORM_DEFINITIONS}")
TARGET_LINK_LIBRARIES(medcoupling interpkernel ${PLATFORM_LIBS})
INSTALL(TARGETS medcoupling DESTINATION ${MED_salomelib_LIBS})

ADD_LIBRARY(medcouplingremapper SHARED ${medcouplingremapper_SOURCES})
TARGET_LINK_LIBRARIES(medcouplingremapper medcoupling ${PLATFORM_LIBS})
SET_TARGET_PROPERTIES(medcouplingremapper PROPERTIES COMPILE_FLAGS "${PLATFORM_DEFINITIONS}")
INSTALL(TARGETS medcouplingremapper DESTINATION ${MED_salomelib_LIBS})

FILE(GLOB medcoupling_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
FILE(GLOB medcoupling_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${medcoupling_HEADERS_HXX} ${medcoupling_HEADERS_TXX} DESTINATION ${MED_salomeinclude_HEADERS})
