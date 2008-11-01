#  SALOMELocalTrace : log on local machine
#
#  Copyright (C) 2003  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
#  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
#  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

include $(top_srcdir)/adm_local/unix/make_common_starter.am

lib_LTLIBRARIES = libInterpKernelTest.la

salomeinclude_HEADERS = CppUnitTest.hxx \
	TransformedTriangleTest.hxx \
	TransformedTriangleIntersectTest.hxx \
	InterpolationTestSuite.hxx \
	SingleElementTetraTests.hxx \
	MultiElementTetraTests.hxx \
	HexaTests.hxx \
	MeshTestToolkit.hxx \
	MeshTestToolkit.txx \
	BBTreeTest.hxx \
	RemapperTest.hxx \
	SingleElementPlanarTests.hxx \
	MultiElement2DTests.hxx \
	InterpolationPlanarTestSuite.hxx \
	PointLocatorTest.hxx \
	MEDMeshMaker.hxx \
	QuadraticPlanarInterpTest.hxx

EXTRA_DIST += \
	BasicMainTest.hxx \
	Interpolation3DTest.cxx \
	Interpolation3DTest.hxx \
	MultiElement3DSurfTests.hxx \
	TestingUtils.hxx \
	perf_test.sh

CLEANFILES = \
        UnitTestsResult

dist_libInterpKernelTest_la_SOURCES= \
	CppUnitTest.cxx \
	TransformedTriangleTest.cxx \
	TransformedTriangleIntersectTest.cxx \
	BBTreeTest.cxx \
	RemapperTest.cxx \
	SingleElementPlanarTests.cxx \
	PointLocatorTest.cxx \
	MEDMeshMaker.cxx \
	InterpolationOptionsTest.hxx \
	InterpolationOptionsTest.cxx \
	QuadraticPlanarInterpTest.cxx \
	QuadraticPlanarInterpTest2.cxx \
	QuadraticPlanarInterpTest3.cxx \
	QuadraticPlanarInterpTest4.cxx

libInterpKernelTest_la_CPPFLAGS = @CPPUNIT_INCLUDES@ $(MED2_INCLUDES) $(HDF5_INCLUDES) \
	-I$(srcdir)/.. -I$(srcdir)/../../MEDWrapper/V2_1/Core -I$(srcdir)/../../MEDMEM \
	-I$(srcdir)/../Geometric2D -DOPTIMIZE -DLOG_LEVEL=0

libInterpKernelTest_la_LDFLAGS  = @CPPUNIT_LIBS@ ../../MEDWrapper/V2_1/Core/libmed_V2_1.la  \
				  ../libinterpkernel.la ../Geometric2D/libInterpGeometric2DAlg.la \
				  ../../MEDMEM/libmedmem.la  -lm

if MED_ENABLE_KERNEL
    libInterpKernelTest_la_CPPFLAGS += ${KERNEL_CXXFLAGS}
    libInterpKernelTest_la_LDFLAGS  += ${KERNEL_LDFLAGS} -lSALOMELocalTrace
endif

# Executables targets
bin_PROGRAMS = TestInterpKernel \
	PerfTest

dist_TestInterpKernel_SOURCES = TestInterpKernel.cxx
dist_PerfTest_SOURCES = PerfTest.cxx

LDADD = $(MED2_LIBS) $(libInterpKernelTest_la_LDFLAGS) \
	libInterpKernelTest.la
if MED_ENABLE_KERNEL
     LDADD += -lSALOMEBasics
endif
AM_CPPFLAGS = $(libInterpKernelTest_la_CPPFLAGS)

UNIT_TEST_PROG = TestInterpKernel
