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

ADD_DEFINITIONS(${CPPUNIT_DEFINITIONS})

INCLUDE_DIRECTORIES(
  ${CPPUNIT_INCLUDE_DIRS}
  ${PTHREADS_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/..
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNEL/GaussPoints
  ${CMAKE_CURRENT_SOURCE_DIR}/../../INTERP_KERNELTest
  )

SET(TestMEDCoupling_SOURCES
  MEDCouplingBasicsTestInterp.cxx
  TestMEDCoupling.cxx
  MEDCouplingBasicsTest0.cxx
  MEDCouplingBasicsTest1.cxx
  MEDCouplingBasicsTest2.cxx
  MEDCouplingBasicsTest3.cxx
  MEDCouplingBasicsTest4.cxx
  MEDCouplingBasicsTest5.cxx
  )

SET(TestMEDCouplingRemapper_SOURCES
  TestMEDCouplingRemapper.cxx
  MEDCouplingRemapperTest.cxx
  MEDCouplingBasicsTest0.cxx
  )

SET(TestMEDCouplingExamples_SOURCES
  MEDCouplingExamplesTest.cxx
  MEDCouplingBasicsTest0.cxx
  )

SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)

ADD_EXECUTABLE(TestMEDCoupling ${TestMEDCoupling_SOURCES})
TARGET_LINK_LIBRARIES(TestMEDCoupling medcoupling ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS})
ADD_TEST(TestMEDCoupling TestMEDCoupling)
SET_TESTS_PROPERTIES(TestMEDCoupling PROPERTIES ENVIRONMENT "${tests_env}")

ADD_EXECUTABLE(TestMEDCouplingRemapper ${TestMEDCouplingRemapper_SOURCES})
TARGET_LINK_LIBRARIES(TestMEDCouplingRemapper medcouplingremapper ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS})
ADD_TEST(TestMEDCouplingRemapper TestMEDCouplingRemapper)
SET_TESTS_PROPERTIES(TestMEDCouplingRemapper PROPERTIES ENVIRONMENT "${tests_env}")

ADD_EXECUTABLE(TestMEDCouplingExamples ${TestMEDCouplingExamples_SOURCES})
TARGET_LINK_LIBRARIES(TestMEDCouplingExamples medcoupling ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS})
ADD_TEST(TestMEDCouplingExamples TestMEDCouplingExamples)
SET_TESTS_PROPERTIES(TestMEDCouplingExamples PROPERTIES ENVIRONMENT "${tests_env}")

INSTALL(TARGETS TestMEDCoupling TestMEDCouplingRemapper TestMEDCouplingExamples DESTINATION ${MEDCOUPLING_INSTALL_BINS})

SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES UnitTestsResult)

# Application tests

SET(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/MEDCoupling)
INSTALL(TARGETS TestMEDCoupling TestMEDCouplingRemapper TestMEDCouplingExamples DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)
