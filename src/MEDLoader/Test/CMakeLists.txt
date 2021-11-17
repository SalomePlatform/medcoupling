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

INCLUDE_DIRECTORIES(
  ${CPPUNIT_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNELTest # For common CppUnitTest.hxx file and TestIKUtils.hxx
  )

SET(TestMEDLoader_SOURCES
  TestMEDLoader.cxx
  MEDLoaderTest.cxx
  )

SALOME_ACCUMULATE_ENVIRONMENT(MEDCOUPLING_RESOURCE_DIR "${CMAKE_BINARY_DIR}/resources")
SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

ADD_EXECUTABLE(TestMEDLoader ${TestMEDLoader_SOURCES})
TARGET_LINK_LIBRARIES(TestMEDLoader medloader InterpKernelTestUtils ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS})

INSTALL(TARGETS TestMEDLoader DESTINATION ${MEDCOUPLING_INSTALL_BINS})

SET(BASE_TESTS TestMEDLoader)

FOREACH(test ${BASE_TESTS})
  ADD_TEST(NAME ${test}
           COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MCTestLauncher.py ${CMAKE_CURRENT_BINARY_DIR}/${test})
  SET_TESTS_PROPERTIES(${test} PROPERTIES ENVIRONMENT "${tests_env}")
ENDFOREACH()

# Application tests

SET(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDLoader)
INSTALL(TARGETS TestMEDLoader DESTINATION ${TEST_INSTALL_DIRECTORY})
INSTALL(FILES MCTestLauncher.py DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)
