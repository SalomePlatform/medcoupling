# Copyright (C) 2012-2013  CEA/DEN, EDF R&D
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

ADD_SUBDIRECTORY(Swig)

IF(CPPUNIT_IS_OK)
  ADD_SUBDIRECTORY(Test)
ENDIF(CPPUNIT_IS_OK)

INCLUDE_DIRECTORIES(
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${OMNIORB_INCLUDE_DIRS}
  ${PTHREAD_INCLUDE_DIRS}
  ${KERNEL_ROOT_DIR}/include/salome
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
SET_TARGET_PROPERTIES(medcalculator PROPERTIES COMPILE_FLAGS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${OMNIORB_DEFINITIONS} ${PLATFORM_DEFINITIONS}")
TARGET_LINK_LIBRARIES(medcalculator medloader medcouplingcorba ${KERNEL_LDFLAGS} ${SALOMELocalTrace} ${SalomeNS} ${OpUtil} ${OMNIORB_LIBS} ${PLATFORM_LIBS})
INSTALL(TARGETS medcalculator DESTINATION ${MED_salomelib_LIBS})

FILE(GLOB medcalculator_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${medcalculator_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
