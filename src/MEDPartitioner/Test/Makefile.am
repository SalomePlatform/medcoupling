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

dist_libMEDPARTITIONERTest_la_SOURCES = \
	MEDPARTITIONERTest.cxx

if MPI_IS_OK
  dist_libMEDPARTITIONERTest_la_SOURCES += MEDPARTITIONERTestPara.cxx
endif

salomeinclude_HEADERS = \
	MEDPARTITIONERTest.hxx

libMEDPARTITIONERTest_la_CPPFLAGS = \
	@CPPUNIT_INCLUDES@ \
	$(MED3_INCLUDES) \
	$(HDF5_INCLUDES) \
	$(LIBXML_CXXFLAGS) \
	$(MPI_INCLUDES) \
	-I$(srcdir)/.. \
	-I$(srcdir)/../../MEDCoupling \
	-I$(srcdir)/../../MEDLoader \
	-I$(srcdir)/../../INTERP_KERNEL/Bases \
	-I$(srcdir)/../../INTERP_KERNEL/GaussPoints \
	-I$(srcdir)/../../INTERP_KERNEL

libMEDPARTITIONERTest_la_LDFLAGS = \
	$(MED3_LIBS_C_ONLY) \
	$(HDF5_LIBS) \
	$(MPI_LIBS) \
	@CPPUNIT_LIBS@ \
	../../MEDLoader/libmedloader.la \
	../libmedpartitioner.la 

# Executables targets

dist_salomescript_DATA = 
UNIT_TEST_PROG =

bin_PROGRAMS = TestMEDPARTITIONER
TestMEDPARTITIONER_SOURCES = TestMEDPARTITIONER.cxx
TestMEDPARTITIONER_CPPFLAGS = $(libMEDPARTITIONERTest_la_CPPFLAGS) \
	-I$(srcdir)/../../INTERP_KERNELTest

if MPI_IS_OK
  TestMEDPARTITIONER_CPPFLAGS +=
	$(MPI_INCLUDES) \
	$(PARMETIS_CPPFLAGS)
else
if MED_ENABLE_METIS
  TestMEDPARTITIONER_CPPFLAGS += $(METIS_CPPFLAGS)
endif
if MED_ENABLE_SCOTCH
  TestMEDPARTITIONER_CPPFLAGS += $(SCOTCH_CPPFLAGS)
endif
endif

TestMEDPARTITIONER_LDADD = \
	$(libMEDPARTITIONERTest_la_LDFLAGS) -lm \
	libMEDPARTITIONERTest.la \
	../../INTERP_KERNEL/libinterpkernel.la \
	../libmedpartitioner.la

UNIT_TEST_PROG += TestMEDPARTITIONER

check: tests

CLEANFILES = UnitTestsResult

clean-local:
        rm -rf *.med *.xml
