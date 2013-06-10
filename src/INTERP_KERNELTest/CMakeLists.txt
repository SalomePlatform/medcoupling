# Copyright (C) 2012-2013  CEA/DEN, EDF R&D
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
# Author : Anthony Geay (CEA/DEN)

INCLUDE_DIRECTORIES(
  ${CPPUNIT_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${MED3_INCLUDE_DIRS}
  ${XDR_INCLUDE_DIRS}
  ${CMAKE_CURRENT_SOURCE_DIR}/../MEDMEM
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/GaussPoints
  ${CMAKE_CURRENT_SOURCE_DIR}/../INTERP_KERNEL/Geometric2D
  )

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
  TestInterpKernelUtils.cxx
)

SET(TestINTERP_KERNEL_SOURCES
  TestInterpKernel.cxx
  )

SET(PerfTest_SOURCES
  PerfTest.cxx
  )

IF(NOT MED_ENABLE_MICROMED)
  SET(InterpKernelTest_SOURCES
    ${InterpKernelTest_SOURCES}
    InterpolationOptionsTest.cxx
    MEDMeshMaker.cxx
    PointLocatorTest.cxx
    )
  
  SET(PerfTest_SOURCES
    PerfTest.cxx
    )
  ADD_EXECUTABLE(PerfTest ${PerfTest_SOURCES})
  SET_TARGET_PROPERTIES(PerfTest PROPERTIES COMPILE_FLAGS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${XDR_DEFINITIONS} ${CPPUNIT_DEFINITIONS} ${PLATFORM_DEFINITIONS}")
  TARGET_LINK_LIBRARIES(PerfTest InterpKernelTest ${CPPUNIT_LIBS} ${PLATFORM_LIBS})
  INSTALL(TARGETS PerfTest DESTINATION ${MED_salomebin_BINS})
ENDIF(NOT MED_ENABLE_MICROMED)

ADD_LIBRARY(InterpKernelTest SHARED ${InterpKernelTest_SOURCES})
SET_TARGET_PROPERTIES(InterpKernelTest PROPERTIES COMPILE_FLAGS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${XDR_DEFINITIONS} ${CPPUNIT_DEFINITIONS} ${PLATFORM_DEFINITIONS}")
TARGET_LINK_LIBRARIES(InterpKernelTest medmem interpkernel ${CPPUNIT_LIBS})

ADD_EXECUTABLE(TestINTERP_KERNEL ${TestINTERP_KERNEL_SOURCES})
SET_TARGET_PROPERTIES(TestINTERP_KERNEL PROPERTIES COMPILE_FLAGS "${HDF5_DEFINITIONS} ${MED3_DEFINITIONS} ${XDR_DEFINITIONS} ${CPPUNIT_DEFINITIONS} ${PLATFORM_DEFINITIONS}")
TARGET_LINK_LIBRARIES(TestINTERP_KERNEL InterpKernelTest ${CPPUNIT_LIBS} ${PLATFORM_LIBS})
ADD_TEST(TestINTERP_KERNEL TestINTERP_KERNEL)

INSTALL(TARGETS TestINTERP_KERNEL DESTINATION ${MED_salomebin_BINS})
INSTALL(TARGETS InterpKernelTest DESTINATION ${MED_salomelib_LIBS})
