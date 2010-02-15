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

lib_LTLIBRARIES = _libParaMEDMEM_Swig.la

salomeinclude_HEADERS =  \
	libParaMEDMEM_Swig.i

SWIG_DEF = libParaMEDMEM_Swig.i libParaMEDMEM_Swig.typemap

SWIG_FLAGS = @SWIG_FLAGS@ -I$(srcdir) $(MPI_INCLUDES)  -I$(srcdir)/../ParaMEDMEM -I$(srcdir)/../MEDCoupling -I$(srcdir)/../MEDCoupling_Swig \
	-I$(srcdir)/../MEDLoader/Swig -I$(srcdir)/../INTERP_KERNEL -I$(srcdir)/../INTERP_KERNEL/Bases -I$(srcdir)/../ParaMEDLoader        \
	-I$(srcdir)/../MEDLoader

dist__libParaMEDMEM_Swig_la_SOURCES   = $(SWIG_DEF)
nodist__libParaMEDMEM_Swig_la_SOURCES = libParaMEDMEM_Swig_wrap.cxx
libParaMEDMEM_Swig.py: libParaMEDMEM_Swig_wrap.cxx

libParaMEDMEM_Swig_wrap.cxx: $(SWIG_DEF)
	$(SWIG) $(SWIG_FLAGS) -o $@  $<

_libParaMEDMEM_Swig_la_CPPFLAGS = $(PYTHON_INCLUDES) \
	$(MED2_INCLUDES) $(HDF5_INCLUDES) @CXXTMPDPTHFLAGS@ \
	-I$(srcdir)/../INTERP_KERNEL \
	$(MPI_INCLUDES) -I$(srcdir)/../ParaMEDMEM -I$(srcdir)/../MEDCoupling_Swig -I$(srcdir)/../INTERP_KERNEL/Bases \
	-I$(srcdir)/../MEDCoupling -I$(srcdir)/../ParaMEDLoader -I$(srcdir)/../MEDLoader

_libParaMEDMEM_Swig_la_LDFLAGS  = -module $(MED2_LIBS) $(HDF5_LIBS) $(PYTHON_LIBS) $(MPI_LIBS) \
	../MEDCoupling/libmedcoupling.la ../INTERP_KERNEL/libinterpkernel.la \
	../ParaMEDMEM/libparamedmem.la ../ParaMEDLoader/libparamedloader.la

if MED_ENABLE_KERNEL
  _libParaMEDMEM_Swig_la_CPPFLAGS += ${KERNEL_CXXFLAGS}
  _libParaMEDMEM_Swig_la_LDFLAGS  += ${KERNEL_LDFLAGS} -lSALOMELocalTrace
endif

CLEANFILES = libParaMEDMEM_Swig_wrap.cxx libParaMEDMEM_Swig.py

dist_salomescript_DATA= test_InterpKernelDEC.py \
			test_NonCoincidentDEC.py \
			test_StructuredCoincodentDEC.py \
			ParaMEDMEM.py libParaMEDMEM_Swig.py
