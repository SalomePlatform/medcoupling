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

SET(TestSauvLoader_SOURCES
  TestSauvLoader.cxx
  SauvLoaderTest.cxx
  )

SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

ADD_EXECUTABLE(TestMEDLoader ${TestMEDLoader_SOURCES})
TARGET_LINK_LIBRARIES(TestMEDLoader medloader InterpKernelTestUtils ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS})
ADD_TEST(TestMEDLoader TestMEDLoader)
SET_TESTS_PROPERTIES(TestMEDLoader PROPERTIES ENVIRONMENT "${tests_env}")

ADD_EXECUTABLE(TestSauvLoader ${TestSauvLoader_SOURCES})

TARGET_LINK_LIBRARIES(TestSauvLoader medloader InterpKernelTestUtils ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS})
ADD_TEST(TestSauvLoader TestSauvLoader)
SET_TESTS_PROPERTIES(TestSauvLoader PROPERTIES ENVIRONMENT "${tests_env}")

INSTALL(TARGETS TestMEDLoader TestSauvLoader DESTINATION ${MEDCOUPLING_INSTALL_BINS})

# Application tests

SET(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDLoader)
INSTALL(TARGETS TestMEDLoader TestSauvLoader DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)
