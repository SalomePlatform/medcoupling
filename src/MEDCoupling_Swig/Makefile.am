#  MED MEDMEM_SWIG : binding of C++ implementation and Python
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
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#

include $(top_srcdir)/adm_local/unix/make_common_starter.am

lib_LTLIBRARIES = _libMEDCoupling_Swig.la

salomeinclude_HEADERS =  \
	libMEDCoupling_Swig.i MEDCouplingTypemaps.i

SWIG_DEF = libMEDCoupling_Swig.i MEDCouplingTypemaps.i

SWIG_FLAGS = @SWIG_FLAGS@ -I$(srcdir) -I$(srcdir)/../MEDCoupling -I$(srcdir)/../INTERP_KERNEL/Bases

dist__libMEDCoupling_Swig_la_SOURCES   = libMEDCoupling_Swig.i MEDCouplingTypemaps.i
nodist__libMEDCoupling_Swig_la_SOURCES = libMEDCoupling_Swig_wrap.cxx
libMEDCoupling_Swig.py: libMEDCoupling_Swig_wrap.cxx

libMEDCoupling_Swig_wrap.cxx: $(SWIG_DEF)
	$(SWIG) $(SWIG_FLAGS) -o $@  $<

_libMEDCoupling_Swig_la_CPPFLAGS = $(CORBA_CXXFLAGS) $(CORBA_INCLUDES) $(PYTHON_INCLUDES) \
	@CXXTMPDPTHFLAGS@ -I$(srcdir)/../MEDCoupling -I$(srcdir)/../INTERP_KERNEL -I$(srcdir)/../INTERP_KERNEL/Bases

_libMEDCoupling_Swig_la_LDFLAGS  = -module $(PYTHON_LIBS) \
	../MEDCoupling/libmedcoupling.la ../INTERP_KERNEL/libinterpkernel.la

CLEANFILES = libMEDCoupling_Swig_wrap.cxx libMEDCoupling_Swig.py

dist_salomescript_DATA= MEDCouplingBasicsTest.py libMEDCoupling_Swig.py
