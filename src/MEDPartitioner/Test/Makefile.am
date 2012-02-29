# Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#  SALOMELocalTrace : log on local machine
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

lib_LTLIBRARIES = libMEDPARTITIONERTest.la

dist_libMEDPARTITIONERTest_la_SOURCES= \
	MEDPARTITIONERTest.cxx
#	MEDPARTITIONERTest_Utils.cxx

#if MED_ENABLE_PARMETIS
#  dist_libMEDPARTITIONERTest_la_SOURCES+= \
#	ParaMEDPARTITIONERTest.hxx
#endif

salomeinclude_HEADERS = \
	MEDPARTITIONERTest.hxx
#	MEDPARTITIONERTest_Utils.hxx

libMEDPARTITIONERTest_la_CPPFLAGS= \
	@CPPUNIT_INCLUDES@ \
	$(MED2_INCLUDES) \
	$(HDF5_INCLUDES) \
	@CXXTMPDPTHFLAGS@ \
	$(BOOST_CPPFLAGS) \
	$(LIBXML_CXXFLAGS) \
	-I$(srcdir)/.. \
	-I$(srcdir)/../../MEDWrapper/V2_1/Core \
	-I$(srcdir)/../../MEDMEM \
	-I$(srcdir)/../../MEDMEMCppTest \
	-I$(srcdir)/../../MEDCoupling \
	-I$(srcdir)/../../MEDLoader \
	-I$(srcdir)/../../INTERP_KERNEL/Bases \
	-I$(srcdir)/../../INTERP_KERNEL/GaussPoints \
	-I$(srcdir)/../../INTERP_KERNEL

libMEDPARTITIONERTest_la_LDFLAGS= $(MED2_LIBS) @CPPUNIT_LIBS@

if MED_ENABLE_KERNEL
    libMEDPARTITIONERTest_la_CPPFLAGS+= ${KERNEL_CXXFLAGS}
    libMEDPARTITIONERTest_la_LDFLAGS+= ${KERNEL_LDFLAGS} -lSALOMELocalTrace 
endif

if MED_ENABLE_PARMETIS
  libMEDPARTITIONERTest_la_CPPFLAGS+= $(PARMETIS_CPPFLAGS)
  libMEDPARTITIONERTest_la_LDFLAGS+= $(PARMETIS_LIBS)
endif
if MED_ENABLE_METIS
  libMEDPARTITIONERTest_la_CPPFLAGS+= $(METIS_CPPFLAGS)
  libMEDPARTITIONERTest_la_LDFLAGS+= $(METIS_LIBS)
endif
if MED_ENABLE_SCOTCH
  libMEDPARTITIONERTest_la_CPPFLAGS+= $(SCOTCH_CPPFLAGS)
  libMEDPARTITIONERTest_la_LDFLAGS+= $(SCOTCH_LIBS)
endif

# Executables targets
bin_PROGRAMS= TestMEDPARTITIONER

TestMEDPARTITIONER_SOURCES= TestMEDPARTITIONER.cxx

TestMEDPARTITIONER_CPPFLAGS= $(libMEDPARTITIONERTest_la_CPPFLAGS) 
TestMEDPARTITIONER_LDADD= \
	$(MED2_LIBS) \
	$(libMEDPARTITIONERTest_la_LDFLAGS) -lm \
	libMEDPARTITIONERTest.la \
	../../MEDWrapper/V2_1/Core/libmed_V2_1.la \
	../../MEDMEM/libmedmem.la \
	../../INTERP_KERNEL/libinterpkernel.la \
	../libmedpartitioner.la

if MED_ENABLE_KERNEL
     TestMEDPARTITIONER_LDADD+= -lSALOMEBasics
endif

dist_salomescript_DATA = 
#testMEDPARTITIONER.py

UNIT_TEST_PROG = TestMEDPARTITIONER

# parallel exe

if MED_ENABLE_PARMETIS
#  bin_PROGRAMS+= TestParaMEDPARTITIONER
#  TestParaMEDPARTITIONER_SOURCES= TestParaMEDPARTITIONER.cxx
#  TestParaMEDPARTITIONER_CPPFLAGS = $(TestMEDPARTITIONER_CPPFLAGS)
#  TestParaMEDPARTITIONER_LDADD = $(TestMEDPARTITIONER_LDADD) $(HDF5_LIBS)
#  dist_salomescript_DATA += testParaMEDPARTITIONER.py
#  TESTS = TestParaMEDPARTITIONER
endif

check: tests

CLEANFILES =            \
        UnitTestsResult
