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


if(MEDCOUPLING_MICROMED)
  add_definitions(-DMEDCOUPLING_MICROMED)
else()
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../medfile_api)
endif()

set(InterpKernelTest_SOURCES
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
    ThreeDSurfProjectionTest.cxx)

set(InterpKernelTestUtils_SOURCES TestInterpKernelUtils.cxx)

set(TestINTERP_KERNEL_SOURCES TestInterpKernel.cxx)

set(PerfTest_SOURCES PerfTest.cxx)

if(NOT MEDCOUPLING_MICROMED)
  set(InterpKernelTest_SOURCES
    ${InterpKernelTest_SOURCES}
    InterpolationOptionsTest.cxx
    MEDMeshMaker.cxx)

  set(PerfTest_SOURCES PerfTest.cxx)
  add_executable(PerfTest ${PerfTest_SOURCES})
  target_link_libraries(PerfTest InterpKernelTest cppunit::cppunit)
  install(TARGETS PerfTest DESTINATION ${MEDCOUPLING_INSTALL_BINS})
endif(NOT MEDCOUPLING_MICROMED)

add_library(InterpKernelTestUtils ${InterpKernelTestUtils_SOURCES})
target_link_libraries(InterpKernelTestUtils interpkernel)
target_include_directories(InterpKernelTestUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(InterpKernelTest ${InterpKernelTest_SOURCES})

set(InterpKernelTest_LIBS medcouplingcpp interpkernel cppunit::cppunit)
if(NOT MEDCOUPLING_MICROMED)
  set(InterpKernelTest_LIBS medloader ${InterpKernelTest_LIBS})
endif(NOT MEDCOUPLING_MICROMED)
target_link_libraries(InterpKernelTest ${InterpKernelTest_LIBS}
                      InterpKernelTestUtils)

add_executable(TestINTERP_KERNEL ${TestINTERP_KERNEL_SOURCES})
target_link_libraries(TestINTERP_KERNEL InterpKernelTest cppunit::cppunit)
add_test(TestINTERP_KERNEL TestINTERP_KERNEL)

# salome_generate_tests_environment(tests_env)
# set_tests_properties(TestINTERP_KERNEL PROPERTIES ENVIRONMENT "${tests_env}")

install(TARGETS TestINTERP_KERNEL DESTINATION ${MEDCOUPLING_INSTALL_BINS})
install(TARGETS InterpKernelTest DESTINATION ${MEDCOUPLING_INSTALL_LIBS})
install(TARGETS InterpKernelTestUtils DESTINATION ${MEDCOUPLING_INSTALL_LIBS})
# install(FILES CppUnitTest.hxx InterpKernelTestExport.hxx
#               TestInterpKernelUtils.hxx
#         DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# Application tests

# set(TEST_INSTALL_DIRECTORY ${MEDCOUPLING_INSTALL_TESTS}/INTERP_KERNELTest)
# install(TARGETS TestINTERP_KERNEL InterpKernelTest
#         DESTINATION ${TEST_INSTALL_DIRECTORY})

# install(
#   FILES CTestTestfileInstall.cmake
#   DESTINATION ${TEST_INSTALL_DIRECTORY}
#   RENAME CTestTestfile.cmake)
