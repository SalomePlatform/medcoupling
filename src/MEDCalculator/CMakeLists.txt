# Copyright (C) 2012-2014  CEA/DEN, EDF R&D
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

ADD_DEFINITIONS(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS} ${OMNIORB_DEFINITIONS})

IF(SALOME_MED_ENABLE_PYTHON)
  ADD_SUBDIRECTORY(Swig)
  IF(SALOME_BUILD_TESTS)
    ADD_SUBDIRECTORY(Test)
  ENDIF(SALOME_BUILD_TESTS)
ENDIF(SALOME_MED_ENABLE_PYTHON)

INCLUDE_DIRECTORIES(
  ${MEDFILE_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${OMNIORB_INCLUDE_DIR} ${OMNIORB_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}/../../idl
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCouplingCorba
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  )

SET(medcalculator_SOURCES
  MEDCalculatorBrowserMesh.cxx
  MEDCalculatorBrowserField.cxx
  MEDCalculatorBrowserStep.cxx
  MEDCalculatorBrowserLiteStruct.cxx
  MEDCalculatorDBRangeSelection.cxx
  MEDCalculatorDBSliceField.cxx
  MEDCalculatorDBField.cxx
)

ADD_LIBRARY(medcalculator SHARED ${medcalculator_SOURCES})
TARGET_LINK_LIBRARIES(medcalculator medloader medcouplingcorba ${KERNEL_LDFLAGS} 
  ${KERNEL_SALOMELocalTrace} ${KERNEL_SalomeNS} ${KERNEL_OpUtil} 
  ${OMNIORB_LIBRARIES} ${PLATFORM_LIBRARIES})
INSTALL(TARGETS medcalculator DESTINATION ${SALOME_INSTALL_LIBS})

FILE(GLOB medcalculator_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${medcalculator_HEADERS_HXX} DESTINATION ${SALOME_INSTALL_HEADERS})
