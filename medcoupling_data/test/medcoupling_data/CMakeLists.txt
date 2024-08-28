# Copyright (C) 2012-2024  CEA, EDF
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
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/../interp_kernel)

set(TestMEDCoupling_SOURCES
  MEDCouplingBasicsTestInterp.cxx
  TestMEDCoupling.cxx
  MEDCouplingBasicsTest0.cxx
  MEDCouplingBasicsTest1.cxx
  MEDCouplingBasicsTest2.cxx
  MEDCouplingBasicsTest3.cxx
  MEDCouplingBasicsTest4.cxx
  MEDCouplingBasicsTest5.cxx)

set(TestMEDCouplingRemapper_SOURCES
  TestMEDCouplingRemapper.cxx
  MEDCouplingRemapperTest.cxx
  MEDCouplingBasicsTest0.cxx)

set(TestMEDCouplingExamples_SOURCES MEDCouplingExamplesTest.cxx
  MEDCouplingBasicsTest0.cxx)

# salome_generate_tests_environment(tests_env)

add_executable(TestMEDCoupling ${TestMEDCoupling_SOURCES})

target_link_libraries(TestMEDCoupling medcouplingcpp cppunit::cppunit)

add_test(TestMEDCoupling TestMEDCoupling)

# set_tests_properties(TestMEDCoupling PROPERTIES ENVIRONMENT "${tests_env}")

add_executable(TestMEDCouplingRemapper ${TestMEDCouplingRemapper_SOURCES})

target_link_libraries(TestMEDCouplingRemapper medcouplingremapper cppunit::cppunit)

add_test(TestMEDCouplingRemapper TestMEDCouplingRemapper)

# set_tests_properties(TestMEDCouplingRemapper PROPERTIES ENVIRONMENT "${tests_env}")

add_executable(TestMEDCouplingExamples ${TestMEDCouplingExamples_SOURCES})

target_link_libraries(TestMEDCouplingExamples medcouplingcpp cppunit::cppunit)

add_test(TestMEDCouplingExamples TestMEDCouplingExamples)

# set_tests_properties(TestMEDCouplingExamples PROPERTIES ENVIRONMENT "${tests_env}")

install(TARGETS TestMEDCoupling TestMEDCouplingRemapper TestMEDCouplingExamples
  DESTINATION ${MEDCOUPLING_INSTALL_BINS})

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES UnitTestsResult)

# Application tests

# set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDCoupling)
# install(
#   TARGETS TestMEDCoupling TestMEDCouplingRemapper TestMEDCouplingExamples
#   DESTINATION ${TEST_INSTALL_DIRECTORY})

# install(
#   FILES CTestTestfileInstall.cmake
#   DESTINATION ${TEST_INSTALL_DIRECTORY}
#   RENAME CTestTestfile.cmake)
