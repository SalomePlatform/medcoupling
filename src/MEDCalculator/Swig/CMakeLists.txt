# Copyright (C) 2012-2015  CEA/DEN, EDF R&D
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

ADD_DEFINITIONS(${PYTHON_DEFINITIONS})

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
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
TARGET_LINK_LIBRARIES(medcalculatorspython interpkernel ${PYTHON_LIBRARIES})
INSTALL(TARGETS medcalculatorspython DESTINATION ${SALOME_INSTALL_LIBS})

SET_SOURCE_FILES_PROPERTIES(MEDCalculator.i PROPERTIES CPLUSPLUS ON)
SET_SOURCE_FILES_PROPERTIES(MEDCalculator.i PROPERTIES SWIG_DEFINITIONS "-shadow")
SET(SWIG_MODULE_MEDCalculator_EXTRA_FLAGS ${SWIG_EXTRA_FLAGS_FOR_NUMPYANDSCIPY})
SWIG_ADD_MODULE(MEDCalculator python MEDCalculator.i)
SWIG_LINK_LIBRARIES(MEDCalculator ${PYTHON_LIBRARIES} medcalculator)

IF(WIN32)
  SET_TARGET_PROPERTIES(_MEDCalculator PROPERTIES DEBUG_OUTPUT_NAME _MEDCalculator_d)
ENDIF(WIN32)

FILE(GLOB medcalculatorspython_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${medcalculatorspython_HEADERS_HXX} DESTINATION ${SALOME_INSTALL_HEADERS})
INSTALL(TARGETS _MEDCalculator DESTINATION ${SALOME_INSTALL_PYTHON})
INSTALL(FILES MEDCalculator.i MEDCalculatorTypemaps.i DESTINATION ${SALOME_INSTALL_HEADERS})
SALOME_INSTALL_SCRIPTS(${CMAKE_CURRENT_BINARY_DIR}/MEDCalculator.py ${SALOME_INSTALL_PYTHON})
