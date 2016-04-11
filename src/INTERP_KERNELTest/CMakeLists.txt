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

ADD_DEFINITIONS(${HDF5_DEFINITIONS} ${MEDFILE_DEFINITIONS} ${XDR_DEFINITIONS} ${CPPUNIT_DEFINITIONS})

IF(MEDCOUPLING_MICROMED)
  ADD_DEFINITIONS(-DMEDCOUPLING_MICROMED)
ENDIF()

INCLUDE_DIRECTORIES(
  ${CPPUNIT_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MEDFILE_INCLUDE_DIRS}
  ${XDR_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDCoupling
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  )

IF(NOT MEDCOUPLING_MICROMED)
  INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../MEDLoader)
ENDIF()

SET(InterpKernelTest_SOURCES
  BBTreeTest.cxx
  CppUnitTest.cxx
  ExprEvalInterpTest.cxx
  QuadraticPlanarInterpTest.cxx
  QuadraticPlanarInterpTest2.cxx
  QuadraticPlanarInterpTest3.cxx
  QuadraticPlanarInterpTest4.cxx
  QuadraticPlanarInterpTest5.cxx
  SingleElementPlanarTests.cxx
  TransformedTriangleIntersectTest.cxx
  TransformedTriangleTest.cxx
  UnitTetra3D2DIntersectionTest.cxx
  UnitTetraIntersectionBaryTest.cxx
  ThreeDSurfProjectionTest.cxx
)

SET(InterpKernelTestUtils_SOURCES
    TestInterpKernelUtils.cxx
)

SET(TestINTERP_KERNEL_SOURCES
  TestInterpKernel.cxx
  )

SET(PerfTest_SOURCES
  PerfTest.cxx
  )

IF(NOT MEDCOUPLING_MICROMED)
  SET(InterpKernelTest_SOURCES
    ${InterpKernelTest_SOURCES}
    InterpolationOptionsTest.cxx
    MEDMeshMaker.cxx
    )

  SET(PerfTest_SOURCES
    PerfTest.cxx
    )
  ADD_EXECUTABLE(PerfTest ${PerfTest_SOURCES})
  TARGET_LINK_LIBRARIES(PerfTest InterpKernelTest ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS})
  INSTALL(TARGETS PerfTest DESTINATION ${MEDCOUPLING_INSTALL_BINS})
ENDIF(NOT MEDCOUPLING_MICROMED)

ADD_LIBRARY(InterpKernelTestUtils ${InterpKernelTestUtils_SOURCES})
TARGET_LINK_LIBRARIES(InterpKernelTestUtils interpkernel)
ADD_LIBRARY(InterpKernelTest ${InterpKernelTest_SOURCES})

SET(InterpKernelTest_LIBS medcoupling interpkernel ${CPPUNIT_LIBRARIES})
IF(NOT MEDCOUPLING_MICROMED)
  SET(InterpKernelTest_LIBS medloader ${InterpKernelTest_LIBS})
ENDIF(NOT MEDCOUPLING_MICROMED)
TARGET_LINK_LIBRARIES(InterpKernelTest ${InterpKernelTest_LIBS} InterpKernelTestUtils)

ADD_EXECUTABLE(TestINTERP_KERNEL ${TestINTERP_KERNEL_SOURCES})
TARGET_LINK_LIBRARIES(TestINTERP_KERNEL InterpKernelTest ${CPPUNIT_LIBRARIES} ${PLATFORM_LIBS})
ADD_TEST(TestINTERP_KERNEL TestINTERP_KERNEL)
SALOME_GENERATE_TESTS_ENVIRONMENT(tests_env)
SET_TESTS_PROPERTIES(TestINTERP_KERNEL PROPERTIES ENVIRONMENT "${tests_env}")

INSTALL(TARGETS TestINTERP_KERNEL DESTINATION ${MEDCOUPLING_INSTALL_BINS})
INSTALL(TARGETS InterpKernelTest DESTINATION ${MEDCOUPLING_INSTALL_LIBS})
INSTALL(TARGETS InterpKernelTestUtils DESTINATION ${MEDCOUPLING_INSTALL_LIBS})
INSTALL(FILES CppUnitTest.hxx InterpKernelTestExport.hxx TestInterpKernelUtils.hxx DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# Application tests

SET(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/INTERP_KERNELTest)
INSTALL(TARGETS TestINTERP_KERNEL InterpKernelTest DESTINATION ${TEST_INSTALL_DIRECTORY})

INSTALL(FILES CTestTestfileInstall.cmake
        DESTINATION ${TEST_INSTALL_DIRECTORY}
        RENAME CTestTestfile.cmake)
