#  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

bin_PROGRAMS= TestMEDCoupling TestMEDCouplingRemapper

TestMEDCoupling_CPPFLAGS=@CPPUNIT_INCLUDES@ @PTHREAD_CFLAGS@ -I$(srcdir)/.. -I$(srcdir)/../../INTERP_KERNEL/Bases -I$(srcdir)/../../INTERP_KERNELTest -I$(srcdir)/../../INTERP_KERNEL -I$(srcdir)/../../INTERP_KERNEL/Geometric2D

TestMEDCoupling_LDFLAGS = @CPPUNIT_LIBS@ ../libmedcoupling.la ../../INTERP_KERNEL/libinterpkernel.la

dist_TestMEDCoupling_SOURCES = TestMEDCoupling.cxx MEDCouplingBasicsTest.hxx MEDCouplingBasicsTest0.cxx MEDCouplingBasicsTest1.cxx \
	MEDCouplingBasicsTestInterp.cxx MEDCouplingBasicsTestData1.hxx

TestMEDCouplingRemapper_CPPFLAGS=@CPPUNIT_INCLUDES@ @PTHREAD_CFLAGS@ -I$(srcdir)/.. -I$(srcdir)/../../INTERP_KERNEL/Bases -I$(srcdir)/../../INTERP_KERNELTest -I$(srcdir)/../../INTERP_KERNEL -I$(srcdir)/../../INTERP_KERNEL/Geometric2D

TestMEDCouplingRemapper_LDFLAGS = @CPPUNIT_LIBS@ ../libmedcouplingremapper.la ../../INTERP_KERNEL/libinterpkernel.la ../libmedcoupling.la

dist_TestMEDCouplingRemapper_SOURCES = TestMEDCouplingRemapper.cxx MEDCouplingRemapperTest.hxx MEDCouplingRemapperTest.cxx MEDCouplingBasicsTest0.cxx
