# Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

FIND_PACKAGE(SWIG REQUIRED)
INCLUDE(${SWIG_USE_FILE})

INCLUDE_DIRECTORIES(
  ${PYTHON_INCLUDE_PATH}
  ${MED3_INCLUDES_DIR}
  ${HDF5_INCLUDES_DIR}
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
SET_SOURCE_FILES_PROPERTIES(MEDCalculator.i PROPERTIES SWIG_FLAGS "-shadow")

SWIG_ADD_MODULE(MEDCalculator python MEDCalculator.i)
SWIG_LINK_LIBRARIES(MEDCalculator ${PYTHON_LIBRARIES} medcalculator)

FILE(GLOB medcalculatorspython_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
INSTALL(FILES ${medcalculatorspython_HEADERS_HXX} DESTINATION ${MED_salomeinclude_HEADERS})
INSTALL(TARGETS _MEDCalculator DESTINATION ${MED_salomepythondir})
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/MEDCalculator.py DESTINATION ${MED_salomescript_PYTHON})
