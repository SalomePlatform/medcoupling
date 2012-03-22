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

# this directory must be recompiled before Test folder

if CPPUNIT_IS_OK
 SUBDIRS=. Test
endif

lib_LTLIBRARIES= libmedpartitioner.la

salomeinclude_HEADERS= \
MEDPARTITIONER_MeshCollection.hxx \
MEDPARTITIONER_MeshCollectionDriver.hxx \
MEDPARTITIONER_MeshCollectionMedXmlDriver.hxx \
MEDPARTITIONER_MeshCollectionMedAsciiDriver.hxx \
MEDPARTITIONER_ParallelTopology.hxx \
MEDPARTITIONER_JointFinder.hxx \
MEDPARTITIONER_Graph.hxx \
MEDPARTITIONER_UserGraph.hxx\
MEDPARTITIONER_Utils.hxx \
MEDPARTITIONER.hxx \
MEDPARTITIONER_ParaDomainSelector.hxx \
MEDPARTITIONER_ConnectZone.hxx \
MEDPARTITIONER_SkyLineArray.hxx \
MEDPARTITIONER_Topology.hxx


if MED_ENABLE_PARMETIS
    salomeinclude_HEADERS += MEDPARTITIONER_MetisGraph.hxx
endif
if MED_ENABLE_METIS
    salomeinclude_HEADERS += MEDPARTITIONER_MetisGraph.hxx
endif
if MED_ENABLE_SCOTCH
    salomeinclude_HEADERS += MEDPARTITIONER_ScotchGraph.hxx
endif

dist_libmedpartitioner_la_SOURCES = \
MEDPARTITIONER_MeshCollection.cxx \
MEDPARTITIONER_MeshCollectionDriver.cxx \
MEDPARTITIONER_MeshCollectionMedXmlDriver.cxx \
MEDPARTITIONER_MeshCollectionMedAsciiDriver.cxx \
MEDPARTITIONER_ParallelTopology.cxx \
MEDPARTITIONER_Graph.cxx\
MEDPARTITIONER_UserGraph.cxx\
MEDPARTITIONER_Utils.cxx \
MEDPARTITIONER_ParaDomainSelector.cxx \
MEDPARTITIONER_ConnectZone.cxx \
MEDPARTITIONER_SkyLineArray.cxx

libmedpartitioner_la_CPPFLAGS = \
	$(MPI_INCLUDES) \
	$(MED3_INCLUDES) \
	$(HDF5_INCLUDES) \
	$(LIBXML_INCLUDES) \
	-I$(srcdir)/../INTERP_KERNEL/Bases \
	-I$(srcdir)/../MEDCoupling \
	-I$(srcdir)/../MEDLoader \
	-I$(srcdir)/../INTERP_KERNEL

libmedpartitioner_la_LDFLAGS = 

if MPI_IS_OK
  dist_libmedpartitioner_la_SOURCES += MEDPARTITIONER_UtilsPara.cxx \
                                       MEDPARTITIONER_JointFinder.cxx

if MED_ENABLE_PARMETIS
    dist_libmedpartitioner_la_SOURCES += MEDPARTITIONER_MetisGraph.cxx
    libmedpartitioner_la_CPPFLAGS += $(PARMETIS_CPPFLAGS)
    libmedpartitioner_la_LDFLAGS += $(PARMETIS_LIBS)
endif

else !MPI_IS_OK

if MED_ENABLE_METIS
    dist_libmedpartitioner_la_SOURCES += MEDPARTITIONER_MetisGraph.cxx
    libmedpartitioner_la_CPPFLAGS += $(METIS_CPPFLAGS)
    libmedpartitioner_la_LDFLAGS += $(METIS_LIBS)
endif

endif

if MED_ENABLE_SCOTCH
    dist_libmedpartitioner_la_SOURCES += MEDPARTITIONER_ScotchGraph.cxx
    libmedpartitioner_la_CPPFLAGS += $(SCOTCH_CPPFLAGS)
    libmedpartitioner_la_LDFLAGS += $(SCOTCH_LIBS)
endif

libmedpartitioner_la_LDFLAGS += \
	$(MED3_LIBS_C_ONLY) \
	$(HDF5_LIBS) \
	$(STDLIB) \
	$(LIBXML_LIBS) \
	$(MPI_LIBS) \
	../INTERP_KERNEL/libinterpkernel.la \
	../MEDCoupling/libmedcoupling.la \
	../MEDLoader/libmedloader.la

# Executables targets
bin_PROGRAMS = medpartitioner
dist_medpartitioner_SOURCES = medpartitioner.cxx
medpartitioner_CPPFLAGS = \
	-I$(srcdir)/../INTERP_KERNEL/Bases \
	-I$(srcdir)/../MEDCoupling \
	-I$(srcdir)/../MEDLoader \
	-I$(srcdir)/../INTERP_KERNEL

if MED_ENABLE_METIS
  medpartitioner_CPPFLAGS += $(METIS_CPPFLAGS)
endif
if MED_ENABLE_SCOTCH
  medpartitioner_CPPFLAGS += $(SCOTCH_CPPFLAGS)
endif

medpartitioner_LDADD = libmedpartitioner.la

if MPI_IS_OK
 
  bin_PROGRAMS += medpartitioner_para
  dist_medpartitioner_para_SOURCES = medpartitioner_para.cxx
  medpartitioner_para_CPPFLAGS = \
	-I$(srcdir)/../INTERP_KERNEL/Bases \
	-I$(srcdir)/../MEDCoupling \
	-I$(srcdir)/../MEDLoader \
	-I$(srcdir)/../INTERP_KERNEL \
	$(MPI_INCLUDES) \
	$(PARMETIS_CPPFLAGS)

  medpartitioner_para_LDADD = libmedpartitioner.la

endif

OBSOLETE_FILES =

EXTRA_DIST += $(OBSOLETE_FILES)
