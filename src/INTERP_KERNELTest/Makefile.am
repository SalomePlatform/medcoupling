#  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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

# Library targets
lib_LTLIBRARIES = libInterpKernelTest.la

UNUSED_FILES =                               \
	Interpolation3DTest.cxx              \
	Interpolation3DTest.hxx              \
	MultiElement3DSurfTests.hxx          \
	TestingUtils.hxx

EXTRA_DIST +=                                \
	perf_test.sh                         \
        $(UNUSED_FILES)

CLEANFILES =                                 \
        UnitTestsResult

libInterpKernelTest_la_SOURCES=              \
	InterpKernelTestExport.hxx           \
	BBTreeTest.hxx                       \
	BBTreeTest.cxx                       \
	CppUnitTest.hxx                      \
	CppUnitTest.cxx                      \
	InterpolationPlanarTestSuite.hxx     \
	ExprEvalInterpTest.hxx               \
	ExprEvalInterpTest.cxx               \
	QuadraticPlanarInterpTest.hxx        \
	QuadraticPlanarInterpTest.cxx        \
	QuadraticPlanarInterpTest2.cxx       \
	QuadraticPlanarInterpTest3.cxx       \
	QuadraticPlanarInterpTest4.cxx       \
	QuadraticPlanarInterpTest5.cxx       \
	SingleElementPlanarTests.hxx         \
	SingleElementPlanarTests.cxx         \
	TransformedTriangleIntersectTest.hxx \
	TransformedTriangleIntersectTest.cxx \
	TransformedTriangleTest.hxx          \
	TransformedTriangleTest.cxx          \
	UnitTetraIntersectionBaryTest.hxx    \
	UnitTetraIntersectionBaryTest.cxx    \
	TestInterpKernelUtils.hxx            \
	TestInterpKernelUtils.cxx

libInterpKernelTest_la_CPPFLAGS =                \
	@CPPUNIT_INCLUDES@                       \
	-I$(srcdir)/../INTERP_KERNEL             \
	-I$(srcdir)/../INTERP_KERNEL/Geometric2D \
	-I$(srcdir)/../INTERP_KERNEL/Bases       \
	-I$(srcdir)/../INTERP_KERNEL/ExprEval    \
	-DOPTIMIZE -DLOG_LEVEL=0

libInterpKernelTest_la_LDFLAGS  =                \
	@CPPUNIT_LIBS@                           \
	../INTERP_KERNEL/libinterpkernel.la      \
	../INTERP_KERNEL/Geometric2D/libInterpGeometric2DAlg.la

if !MED_ENABLE_MICROMED
   libInterpKernelTest_la_SOURCES +=         \
	InterpolationOptionsTest.hxx         \
	InterpolationOptionsTest.cxx         \
	InterpolationTestSuite.hxx           \
	MEDMeshMaker.hxx                     \
	MEDMeshMaker.cxx                     \
	MeshTestToolkit.hxx                  \
	MeshTestToolkit.txx                  \
	MultiElement2DTests.hxx              \
	MultiElementTetraTests.hxx           \
	PointLocatorTest.hxx                 \
	PointLocatorTest.cxx                 \
	SingleElementTetraTests.hxx

   libInterpKernelTest_la_CPPFLAGS +=            \
	-DDISABLE_MICROMED                       \
	$(MED2_INCLUDES) $(HDF5_INCLUDES)        \
	-I$(srcdir)/../MEDWrapper/V2_1/Core      \
	-I$(srcdir)/../MEDMEM

   libInterpKernelTest_la_LDFLAGS +=             \
	../MEDWrapper/V2_1/Core/libmed_V2_1.la   \
	../MEDMEM/libmedmem.la -lm 
endif

if MED_ENABLE_KERNEL
    libInterpKernelTest_la_CPPFLAGS += ${KERNEL_CXXFLAGS}
    libInterpKernelTest_la_LDFLAGS  += ${KERNEL_LDFLAGS} -lSALOMELocalTrace
endif

# Executables targets
bin_PROGRAMS = TestINTERP_KERNEL

if !MED_ENABLE_MICROMED
   bin_PROGRAMS += PerfTest
endif

TestINTERP_KERNEL_SOURCES =      \
	BasicMainTest.hxx        \
	TestInterpKernel.cxx

TestINTERP_KERNEL_CPPFLAGS = $(libInterpKernelTest_la_CPPFLAGS)
TestINTERP_KERNEL_LDADD = $(libInterpKernelTest_la_LDFLAGS) libInterpKernelTest.la

if !MED_ENABLE_MICROMED
   TestINTERP_KERNEL_SOURCES += HexaTests.hxx
endif

PerfTest_SOURCES =               \
	PerfTest.cxx             \
	MeshTestToolkit.hxx      \
	MeshTestToolkit.txx

PerfTest_CPPFLAGS = $(libInterpKernelTest_la_CPPFLAGS)
PerfTest_LDADD = $(MED2_LIBS) $(libInterpKernelTest_la_LDFLAGS) libInterpKernelTest.la

if MED_ENABLE_KERNEL
     TestINTERP_KERNEL_LDADD += -lSALOMEBasics
     PerfTest_LDADD   += -lSALOMEBasics
endif

UNIT_TEST_PROG = TestINTERP_KERNEL

check : tests

clean-local:
	rm -rf *.fig
