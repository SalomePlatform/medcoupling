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

lib_LTLIBRARIES = _libMEDLoader_Swig.la

salomeinclude_HEADERS =  \
	libMEDLoader_Swig.i

SWIG_DEF = libMEDLoader_Swig.i

SWIG_FLAGS = @SWIG_FLAGS@ -I$(srcdir) -I$(srcdir)/.. -I$(srcdir)/../../MEDCoupling -I$(srcdir)/../../MEDCoupling_Swig \
	-I$(srcdir)/../../INTERP_KERNEL -I$(srcdir)/../../INTERP_KERNEL/Bases

dist__libMEDLoader_Swig_la_SOURCES   = $(SWIG_DEF)
nodist__libMEDLoader_Swig_la_SOURCES = libMEDLoader_Swig_wrap.cxx
libMEDLoader_Swig.py: libMEDLoader_Swig_wrap.cxx

libMEDLoader_Swig_wrap.cxx: $(SWIG_DEF)
	$(SWIG) $(SWIG_FLAGS) -o $@  $<

_libMEDLoader_Swig_la_CPPFLAGS = $(PYTHON_INCLUDES) \
	$(MED2_INCLUDES) $(HDF5_INCLUDES) @CXXTMPDPTHFLAGS@ \
	-I$(srcdir)/../../INTERP_KERNEL \
	-I$(srcdir)/.. -I$(srcdir)/../../MEDCoupling_Swig -I$(srcdir)/../../INTERP_KERNEL/Bases \
	-I$(srcdir)/../../MEDCoupling

_libMEDLoader_Swig_la_LDFLAGS  = -module $(MED2_LIBS) $(HDF5_LIBS) $(PYTHON_LIBS) \
	../../MEDCoupling/libmedcoupling.la ../../INTERP_KERNEL/libinterpkernel.la \
	../libmedloader.la

CLEANFILES = libMEDLoader_Swig_wrap.cxx libMEDLoader_Swig.py

nodist_salomescript_DATA= libMEDLoader_Swig.py

