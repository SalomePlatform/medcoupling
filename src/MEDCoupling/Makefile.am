#  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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

if CPPUNIT_IS_OK
  SUBDIRS = . Test
endif

# first '.' says that this folder must be compiled before Test (MEDMEM/Test uses MEDMEM)

lib_LTLIBRARIES = libmedcoupling.la


salomeinclude_HEADERS = MEDCoupling.hxx \
MEDCouplingFieldDouble.hxx  MEDCouplingMesh.hxx MEDCouplingUMesh.hxx TimeLabel.hxx \
MEDCouplingField.hxx MEDCouplingNormalizedUnstructuredMesh.hxx MemArray.hxx        \
MEDCouplingNormalizedUnstructuredMesh.txx  MemArray.txx RefCountObject.hxx         \
MEDCouplingRMesh.hxx MEDCouplingTimeDiscretization.hxx                             \
MEDCouplingFieldDiscretization.hxx MEDCouplingPointSet.hxx                         \
MEDCouplingUMeshDesc.hxx

# Libraries targets

dist_libmedcoupling_la_SOURCES = \
	MEDCouplingField.cxx  MEDCouplingFieldDouble.cxx       \
	MEDCouplingUMesh.cxx  MemArray.cxx  TimeLabel.cxx      \
	MEDCouplingRMesh.cxx MEDCouplingTimeDiscretization.cxx \
	MEDCouplingFieldDiscretization.cxx                     \
	MEDCouplingPointSet.cxx MEDCouplingUMeshDesc.cxx

libmedcoupling_la_LDFLAGS= 

libmedcoupling_la_CPPFLAGS=  @CXXTMPDPTHFLAGS@ \
	-I$(srcdir) -I$(srcdir)/../INTERP_KERNEL/Bases -I$(srcdir)/../INTERP_KERNEL

# the geom2D library is included in the interpkernel one
libmedcoupling_la_LIBADD= ../INTERP_KERNEL/libinterpkernel.la

AM_CPPFLAGS= $(libinterpkernel_la_CPPFLAGS)
LDADD= $(libinterpkernel_la_LDFLAGS)
if MED_ENABLE_KERNEL
     LDADD+=-lSALOMEBasics
endif

EXTRA_DIST += \
	MEDCouplingFieldDouble.hxx                \
	MEDCouplingTimeDiscretization.hxx         \
	MEDCouplingFieldDiscretization.hxx        \
        MEDCouplingMesh.hxx                       \
        MEDCouplingUMesh.hxx                      \
        MEDCouplingRMesh.hxx                      \
        TimeLabel.hxx                             \
	MEDCouplingField.hxx                      \
	MEDCouplingNormalizedUnstructuredMesh.hxx \
	MemArray.hxx                              \
	MEDCouplingNormalizedUnstructuredMesh.txx \
	MemArray.txx                              \
	RefCountObject.hxx                        \
	MEDCouplingPointSet.hxx                   \
	MEDCouplingUMeshDesc.hxx
