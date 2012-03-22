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

include $(top_srcdir)/adm_local/unix/make_common_starter.am

EXTRA_DIST += CMakeLists.txt

bin_PROGRAMS= TestMEDLoader TestSauvLoader

TestMEDLoader_CPPFLAGS=@CPPUNIT_INCLUDES@ @PTHREAD_CFLAGS@ -I$(srcdir)/.. -I$(srcdir)/../../INTERP_KERNEL/Bases -I$(srcdir)/../../INTERP_KERNELTest -I$(srcdir)/../../INTERP_KERNEL -I$(srcdir)/../../INTERP_KERNEL/Geometric2D -I$(srcdir)/../../MEDCoupling
TestSauvLoader_CPPFLAGS=$(TestMEDLoader_CPPFLAGS) $(MED3_INCLUDES) $(HDF5_INCLUDES)

TestMEDLoader_LDFLAGS = @CPPUNIT_LIBS@ ../libmedloader.la ../../MEDCoupling/libmedcoupling.la ../../INTERP_KERNEL/libinterpkernel.la
TestSauvLoader_LDFLAGS=$(TestMEDLoader_LDFLAGS)

dist_TestMEDLoader_SOURCES = TestMEDLoader.cxx MEDLoaderTest.cxx MEDLoaderTest.hxx
dist_TestSauvLoader_SOURCES = TestSauvLoader.cxx SauvLoaderTest.cxx SauvLoaderTest.hxx

UNIT_TEST_PROG = TestMEDLoader TestSauvLoader

check : tests

CLEANFILES = \
        UnitTestsResult

clean-local:
	rm -rf *.med
