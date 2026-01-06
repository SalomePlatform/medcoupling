# Copyright (C) 2012-2026  CEA, EDF
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

include_directories(
  ${CPPUNIT_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNELTest # For common
                                                      # CppUnitTest.hxx file and
                                                      # TestIKUtils.hxx
)

set(TestMEDLoader_SOURCES TestMEDLoader.cxx CrackAlgoTest.cxx MEDLoaderTest.cxx)

salome_accumulate_environment(MEDCOUPLING_RESOURCE_DIR
                              "${CMAKE_BINARY_DIR}/resources")
salome_generate_tests_environment(tests_env)

add_executable(TestMEDLoader ${TestMEDLoader_SOURCES})
target_link_libraries(TestMEDLoader medloader InterpKernelTestUtils
                      ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS})

install(TARGETS TestMEDLoader DESTINATION ${MEDCOUPLING_INSTALL_BINS})

set(BASE_TESTS TestMEDLoader)

foreach(test ${BASE_TESTS})
  add_test(
    NAME ${test}
    COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/MCTestLauncher.py
            ${CMAKE_CURRENT_BINARY_DIR}/${test})
  set_tests_properties(${test} PROPERTIES ENVIRONMENT "${tests_env}")
endforeach()

# Application tests

set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDLoader)
install(TARGETS TestMEDLoader DESTINATION ${TEST_INSTALL_DIRECTORY})
install(FILES MCTestLauncher.py DESTINATION ${TEST_INSTALL_DIRECTORY})

install(
  FILES CTestTestfileInstall.cmake
  DESTINATION ${TEST_INSTALL_DIRECTORY}
  RENAME CTestTestfile.cmake)
