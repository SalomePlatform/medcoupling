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

#  MED MEDMEM_SWIG : binding of C++ implementation and Python
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

salomepython_PYTHON = ParaMEDMEM.py
salomepyexec_LTLIBRARIES = _ParaMEDMEM.la

salomeinclude_HEADERS =  \
	ParaMEDMEM.i

SWIG_DEF = ParaMEDMEM.i ParaMEDMEM.typemap

SWIG_FLAGS = @SWIG_FLAGS@ -I$(srcdir) $(MPI_INCLUDES)  -I$(srcdir)/../ParaMEDMEM -I$(srcdir)/../MEDCoupling -I$(srcdir)/../MEDCoupling_Swig \
	-I$(srcdir)/../MEDLoader/Swig -I$(srcdir)/../INTERP_KERNEL -I$(srcdir)/../INTERP_KERNEL/Bases -I$(srcdir)/../ParaMEDLoader        \
	-I$(srcdir)/../MEDLoader

dist__ParaMEDMEM_la_SOURCES   = $(SWIG_DEF)
nodist__ParaMEDMEM_la_SOURCES = ParaMEDMEM_wrap.cxx
ParaMEDMEM.py : ParaMEDMEM_wrap.cxx

ParaMEDMEM_wrap.cxx: $(SWIG_DEF)
	$(SWIG) $(SWIG_FLAGS) -o $@  $<

_ParaMEDMEM_la_CPPFLAGS = $(PYTHON_INCLUDES) \
	$(MED2_INCLUDES) $(HDF5_INCLUDES) @CXXTMPDPTHFLAGS@ \
	-I$(srcdir)/../INTERP_KERNEL \
	$(MPI_INCLUDES) -I$(srcdir)/../ParaMEDMEM -I$(srcdir)/../MEDCoupling_Swig -I$(srcdir)/../INTERP_KERNEL/Bases \
	-I$(srcdir)/../MEDCoupling -I$(srcdir)/../ParaMEDLoader -I$(srcdir)/../MEDLoader -I$(srcdir)/../MEDLoader/Swig

_ParaMEDMEM_la_LDFLAGS  = -module $(MED2_LIBS) $(HDF5_LIBS) $(PYTHON_LIBS) $(MPI_LIBS) \
	../MEDCoupling/libmedcoupling.la ../INTERP_KERNEL/libinterpkernel.la \
	../ParaMEDMEM/libparamedmem.la ../ParaMEDLoader/libparamedloader.la

if MED_ENABLE_KERNEL
  _ParaMEDMEM_la_CPPFLAGS += ${KERNEL_CXXFLAGS}
  _ParaMEDMEM_la_LDFLAGS  += ${KERNEL_LDFLAGS} -lSALOMELocalTrace
endif

CLEANFILES = ParaMEDMEM_wrap.cxx ParaMEDMEM.py

dist_salomescript_DATA= test_InterpKernelDEC.py \
			test_NonCoincidentDEC.py \
			test_StructuredCoincidentDEC.py \
			ParaMEDMEM.py
