# Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
# Author : Anthony Geay (CEA/DEN)

#  MED MEDMEM_SWIG : binding of C++ implementation and Python
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

EXTRA_DIST += SauvLoaderTest.py MEDLoaderExamplesTest.py

lib_LTLIBRARIES = _MEDLoader.la

salomeinclude_HEADERS =  \
	MEDLoader.i MEDLoaderTypemaps.i MEDLoaderCommon.i

SWIG_DEF = MEDLoader.i MEDLoaderTypemaps.i MEDLoaderCommon.i

SWIG_FLAGS = @SWIG_FLAGS@ -I$(srcdir) -I$(srcdir)/.. -I$(srcdir)/../../MEDCoupling -I$(srcdir)/../../MEDCoupling_Swig \
	-I$(srcdir)/../../INTERP_KERNEL -I$(srcdir)/../../INTERP_KERNEL/Bases

dist__MEDLoader_la_SOURCES   = $(SWIG_DEF)
nodist__MEDLoader_la_SOURCES = MEDLoader_wrap.cxx
MEDLoader.py: MEDLoader_wrap.cxx

MEDLoader_wrap.cxx: $(SWIG_DEF)
	$(SWIG) $(SWIG_FLAGS) -o $@  $<

_MEDLoader_la_CPPFLAGS = $(PYTHON_INCLUDES) \
	$(MED3_INCLUDES) $(HDF5_INCLUDES) @CXXTMPDPTHFLAGS@ \
	-I$(srcdir)/../../INTERP_KERNEL \
	-I$(srcdir)/.. -I$(srcdir)/../../MEDCoupling_Swig -I$(srcdir)/../../INTERP_KERNEL/Bases \
	-I$(srcdir)/../../MEDCoupling

_MEDLoader_la_LDFLAGS  = -module $(MED3_LIBS_C_ONLY) $(HDF5_LIBS) $(PYTHON_LIBS) \
	../../MEDCoupling/libmedcoupling.la ../../INTERP_KERNEL/libinterpkernel.la \
	../libmedloader.la

CLEANFILES = MEDLoader_wrap.cxx MEDLoader.py

dist_salomescript_SCRIPTS= med2sauv sauv2med case2med med2case

dist_salomescript_DATA= MEDLoaderDataForTest.py MEDLoaderTest.py MEDLoaderTest2.py MEDLoaderTest3.py MEDLoaderCouplingTrainingSession.py CaseIO.py CaseReader.py CaseWriter.py VTKReader.py MEDLoaderSplitter.py libMEDLoader_Swig.py medutilities.py
nodist_salomescript_DATA = MEDLoader.py

UNIT_TEST_PROG = MEDLoaderTest.py MEDLoaderTest2.py

install-data-hook:
	@for f in $(dist_salomescript_SCRIPTS) ; do \
	  chmod -f a+x $(DESTDIR)$(salomescriptdir)/$$f ; \
	done
