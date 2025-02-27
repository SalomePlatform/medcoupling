# Copyright (C) 2024  CEA, EDF
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

INCLUDE_DIRECTORIES(
  ${CPPUNIT_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDLoader
  ${CMAKE_CURRENT_SOURCE_DIR}/../../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNELTest # For common CppUnitTest.hxx file and TestIKUtils.hxx
  )

SET(TestShapeRecogn_SOURCES
  TestShapeRecogn.cxx
  MathOpsTest.cxx
  ShapeRecognTest.cxx
  PlaneTest.cxx
  CylinderTest.cxx
  ConeTest.cxx
  SphereTest.cxx
  TorusTest.cxx
)

SALOME_ACCUMULATE_ENVIRONMENT(MEDCOUPLING_RESOURCE_DIR "${CMAKE_BINARY_DIR}/resources")
SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)


IF(NOT MEDCOUPLING_MICROMED)
  SET(TESTSHAPE_RECOGN0 TestShapeRecogn)
  ADD_EXECUTABLE(${TESTSHAPE_RECOGN0} ${TestShapeRecogn_SOURCES})
  TARGET_LINK_LIBRARIES(${TESTSHAPE_RECOGN0} shaperecogn InterpKernelTestUtils medloader medcouplingcpp ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS} ${CBLAS_LIBRARIES})

  INSTALL(TARGETS ${TESTSHAPE_RECOGN0} DESTINATION ${MEDCOUPLING_INSTALL_BINS})

  #### Test section

  SET(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/ShapeRecogn)
  SET(BASE_TESTS ${TESTSHAPE_RECOGN0})

  FOREACH(test ${BASE_TESTS})
    ADD_TEST(NAME ${test}
            COMMAND ${CMAKE_CURRENT_BINARY_DIR}/${test})
    SET_TESTS_PROPERTIES(${test} PROPERTIES ENVIRONMENT "${tests_env}")
  ENDFOREACH()
  INSTALL(TARGETS ${TESTSHAPE_RECOGN0} DESTINATION ${TEST_INSTALL_DIRECTORY})


  SET(TESTSHAPE_RECOGN1 "TestShapeRecogn.py")
  SET(test ${TESTSHAPE_RECOGN1})
  SET(testname "TestShapeRecognPy")
  ADD_TEST(NAME ${testname} COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test})
  SET_TESTS_PROPERTIES(${testname} PROPERTIES ENVIRONMENT "${tests_env}")
  INSTALL(FILES ${CMAKE_CURRENT_SOURCE_DIR}/${TESTSHAPE_RECOGN1} DESTINATION ${TEST_INSTALL_DIRECTORY})

  INSTALL(FILES CTestTestfileInstall.cmake
          DESTINATION ${TEST_INSTALL_DIRECTORY}
          RENAME CTestTestfile.cmake)

  #### End Test section
ENDIF(NOT MEDCOUPLING_MICROMED)
