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

#  MED MEDCoupling_Swig : binding of C++ implementation and Python
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

lib_LTLIBRARIES = _MEDCoupling.la _MEDCouplingRemapper.la

salomeinclude_HEADERS =  \
	MEDCoupling.i MEDCouplingRemapper.i MEDCouplingTypemaps.i

BUILT_SOURCES = MEDCouplingWRAP.cxx MEDCouplingRemapperWRAP.cxx

SWIG_FLAGS = @SWIG_FLAGS@ -I$(srcdir) -I$(srcdir)/../MEDCoupling -I$(srcdir)/../INTERP_KERNEL/Bases -I$(srcdir)/../INTERP_KERNEL

dist__MEDCoupling_la_SOURCES   = MEDCoupling.i MEDCouplingTypemaps.i
nodist__MEDCoupling_la_SOURCES = MEDCouplingWRAP.cxx
MEDCoupling.py: MEDCouplingWRAP.cxx

MEDCouplingWRAP.cxx: MEDCoupling.i MEDCouplingTypemaps.i
	$(SWIG) $(SWIG_FLAGS) -o $@  $<

_MEDCoupling_la_CPPFLAGS = $(PYTHON_INCLUDES) \
	@CXXTMPDPTHFLAGS@ -I$(srcdir)/../MEDCoupling -I$(srcdir)/../INTERP_KERNEL -I$(srcdir)/../INTERP_KERNEL/Bases

_MEDCoupling_la_LDFLAGS  = -module $(PYTHON_LIBS) \
	../MEDCoupling/libmedcoupling.la ../INTERP_KERNEL/libinterpkernel.la

dist__MEDCouplingRemapper_la_SOURCES = MEDCouplingRemapper.i
nodist__MEDCouplingRemapper_la_SOURCES = MEDCouplingRemapperWRAP.cxx
MEDCouplingRemapper.py: MEDCouplingRemapperWRAP.cxx

MEDCouplingRemapperWRAP.cxx: MEDCouplingRemapper.i MEDCouplingTypemaps.i
	$(SWIG) $(SWIG_FLAGS) -o $@  $<

_MEDCouplingRemapper_la_CPPFLAGS = $(PYTHON_INCLUDES) \
	@CXXTMPDPTHFLAGS@ -I$(srcdir)/../MEDCoupling -I$(srcdir)/../INTERP_KERNEL -I$(srcdir)/../INTERP_KERNEL/Bases

_MEDCouplingRemapper_la_LDFLAGS  = -module $(PYTHON_LIBS) \
	../MEDCoupling/libmedcouplingremapper.la ../INTERP_KERNEL/libinterpkernel.la

CLEANFILES = MEDCouplingWRAP.cxx MEDCoupling.py MEDCouplingRemapperWRAP.cxx MEDCouplingRemapper.py

dist_salomescript_DATA= MEDCouplingBasicsTest.py MEDCouplingRemapperTest.py MEDCoupling.py MEDCouplingRemapper.py MEDCouplingDataForTest.py libMEDCoupling_Swig.py

UNIT_TEST_PROG = MEDCouplingBasicsTest.py MEDCouplingRemapperTest.py

#check : tests
