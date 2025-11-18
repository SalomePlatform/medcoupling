# Copyright (C) 2024  CEA, EDF
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

include_directories(
  ${CPPUNIT_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${LAPACKE_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNELTest # For common
                                                      # CppUnitTest.hxx file and
                                                      # TestIKUtils.hxx
)

set(TestShapeRecogn_SOURCES
    TestShapeRecogn.cxx
    MathOpsTest.cxx
    ShapeRecognTest.cxx
    PlaneTest.cxx
    CylinderTest.cxx
    ConeTest.cxx
    SphereTest.cxx
    TorusTest.cxx)

salome_accumulate_environment(MEDCOUPLING_RESOURCE_DIR
                              "${CMAKE_BINARY_DIR}/resources")
salome_generate_tests_environment(tests_env)

if(NOT MEDCOUPLING_MICROMED)
  set(TESTSHAPE_RECOGN0 TestShapeRecogn)
  add_executable(${TESTSHAPE_RECOGN0} ${TestShapeRecogn_SOURCES})
  target_link_libraries(
    ${TESTSHAPE_RECOGN0}
    shaperecogn
    InterpKernelTestUtils
    medloader
    medcouplingcpp
    ${CPPUNIT_LIBRARIES}
    ${PLATFORM_LIBS}
    ${CBLAS_LIBRARIES})

  install(TARGETS ${TESTSHAPE_RECOGN0} DESTINATION ${MEDCOUPLING_INSTALL_BINS})

  # Test section

  set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/ShapeRecogn)
  set(BASE_TESTS ${TESTSHAPE_RECOGN0})

  foreach(test ${BASE_TESTS})
    add_test(NAME ${test} COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${test})
    set_tests_properties(${test} PROPERTIES ENVIRONMENT "${tests_env}")
  endforeach()
  install(TARGETS ${TESTSHAPE_RECOGN0} DESTINATION ${TEST_INSTALL_DIRECTORY})

  set(TESTSHAPE_RECOGN1 "TestShapeRecogn.py")
  set(test ${TESTSHAPE_RECOGN1})
  set(testname "TestShapeRecognPy")
  add_test(NAME ${testname} COMMAND ${PYTHON_EXECUTABLE}
                                    ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  set_tests_properties(${testname} PROPERTIES ENVIRONMENT "${tests_env}")
  install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${TESTSHAPE_RECOGN1}
          DESTINATION ${TEST_INSTALL_DIRECTORY})

  install(
    FILES CTestTestfileInstall.cmake
    DESTINATION ${TEST_INSTALL_DIRECTORY}
    RENAME CTestTestfile.cmake)

  # End Test section
endif(NOT MEDCOUPLING_MICROMED)
