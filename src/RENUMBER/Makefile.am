# Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#  MED MEDMEM : MED files in memory

include $(top_srcdir)/adm_local/unix/make_common_starter.am

# this directory must be recompiled before Test folder

if CPPUNIT_IS_OK
# SUBDIRS=. Test
endif

lib_LTLIBRARIES= librenumber.la

salomeinclude_HEADERS= \
RENUMBER_Renumbering.hxx \
RenumberingFactory.hxx \
RENUMBERDefines.hxx

if MED_ENABLE_METIS
    salomeinclude_HEADERS+= RENUMBER_METISRenumbering.hxx
endif
if BOOST_IS_OK
    salomeinclude_HEADERS+= RENUMBER_BOOSTRenumbering.hxx
endif

dist_librenumber_la_SOURCES= \
RENUMBER_Renumbering.cxx \
RenumberingFactory.cxx

if MED_ENABLE_METIS
    dist_librenumber_la_SOURCES+= RENUMBER_METISRenumbering.cxx
endif
if BOOST_IS_OK
    dist_librenumber_la_SOURCES+= RENUMBER_BOOSTRenumbering.cxx 
endif

librenumber_la_CPPFLAGS= $(MED3_INCLUDES) $(HDF5_INCLUDES) @CXXTMPDPTHFLAGS@ \
	$(BOOST_CPPFLAGS) \
	-I$(srcdir)/../MEDMEM -I$(srcdir)/../INTERP_KERNEL/Bases \
	-I$(srcdir)/../INTERP_KERNEL/GaussPoints \
	-I$(srcdir)/../INTERP_KERNEL

librenumber_la_LDFLAGS= 
#libmedsplitter_la_LDFLAGS= $(MED3_LIBS_C_ONLY) $(HDF5_LIBS) $(STDLIB) \
#	../MEDMEM/libmedmem.la

if MED_ENABLE_METIS
  librenumber_la_CPPFLAGS+= $(METIS_CPPFLAGS)
  librenumber_la_LDFLAGS+= $(METIS_LIBS)
endif 
if BOOST_IS_OK
  librenumber_la_CPPFLAGS+= $(BOOST_CPPFLAGS) -DENABLE_BOOST
  librenumber_la_LDFLAGS+= $(BOOST_LIBS)
endif
if MED_ENABLE_KERNEL
  librenumber_la_CPPFLAGS+= ${KERNEL_CXXFLAGS}
  librenumber_la_LDFLAGS+= ${KERNEL_LDFLAGS} -lSALOMELocalTrace 
endif

librenumber_la_LDFLAGS+= $(MED3_LIBS_C_ONLY) $(HDF5_LIBS) $(STDLIB) \
	../MEDMEM/libmedmem.la ../INTERP_KERNEL/libinterpkernel.la

# Executables targets
bin_PROGRAMS= renumber

dist_renumber_SOURCES= renumbering.cxx

renumber_CPPFLAGS= $(librenumber_la_CPPFLAGS)
renumber_LDADD= $(librenumber_la_LDFLAGS) -lm $(BOOST_LIBS) librenumber.la
if MED_ENABLE_KERNEL
  renumber_LDADD+= -lSALOMEBasics
endif

OBSOLETE_FILES = 
#	MEDSPLITTER_SequentialTopology.cxx \
#	test_HighLevelAPI.cxx

EXTRA_DIST += $(OBSOLETE_FILES) testRenumbering.py

