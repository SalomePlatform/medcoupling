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

INCLUDE(${SALOME_MACROS_DIR}/InstallAndCompilePythonFile.cmake)
FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDLoader/Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling_Swig
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  )

SET(medcalculatorspython_SOURCES
  SPythonInterpreter.cxx
  SPythonParser.cxx
  )

ADD_LIBRARY(medcalculatorspython SHARED ${medcalculatorspython_SOURCES})
TARGET_LINK_LIBRARIES(medcalculatorspython interpkernel ${PYTHON_LIBS})
INSTALL(TARGETS medcalculatorspython DESTINATION ${MED_salomelib_LIBS})

SET_SOURCE_FILES_PROPERTIES(MEDCalculator.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDCalculator.i PROPERTIES SWIG_DEFINITIONS "-shadow")
SWIG_ADD_MODULE(MEDCalculator python MEDCalculator.i)
SWIG_LINK_LIBRARIES(MEDCalculator ${PYTHON_LIBS} medcalculator)
SET_TARGET_PROPERTIES(_MEDCalculator PROPERTIES COMPILE_FLAGS "${PLATFORM_DEFINITIONS}")

IF(${MACHINE} STREQUAL WINDOWS)
  SET_TARGET_PROPERTIES(_MEDCalculator PROPERTIES DEBUG_OUTPUT_NAME _MEDCalculator_d)
ENDIF(${MACHINE} STREQUAL WINDOWS)

FILE(GLOB medcalculatorspython_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${medcalculatorspython_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
INSTALL(TARGETS _MEDCalculator DESTINATION ${MED_salomepythondir})
INSTALL(FILES MEDCalculator.i MEDCalculatorTypemaps.i DESTINATION ${MED_salomeinclude_HEADERS})
INSTALL_AND_COMPILE_PYTHON_FILE(${CMAKE_CURRENT_BINARY_DIR}/MEDCalculator.py ${MED_salomescript_PYTHON})
