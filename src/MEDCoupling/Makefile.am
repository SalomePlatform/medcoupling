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

if CPPUNIT_IS_OK
  SUBDIRS = . Test
endif

# first '.' says that this folder must be compiled before Test (MEDMEM/Test uses MEDMEM)

lib_LTLIBRARIES = libmedcoupling.la libmedcouplingremapper.la


salomeinclude_HEADERS = MEDCoupling.hxx \
MEDCouplingFieldDouble.hxx  MEDCouplingMesh.hxx MEDCouplingUMesh.hxx MEDCouplingTimeLabel.hxx     \
MEDCouplingField.hxx MEDCouplingNormalizedUnstructuredMesh.hxx MEDCouplingMemArray.hxx            \
MEDCouplingNormalizedUnstructuredMesh.txx  MEDCouplingMemArray.txx MEDCouplingRefCountObject.hxx  \
MEDCouplingCMesh.hxx MEDCouplingTimeDiscretization.hxx                                            \
MEDCouplingFieldDiscretization.hxx MEDCouplingPointSet.hxx MEDCouplingPointSet.txx                \
MEDCouplingUMeshDesc.hxx MEDCouplingNatureOfField.hxx                                             \
MEDCouplingNormalizedCartesianMesh.hxx MEDCouplingNormalizedCartesianMesh.txx                     \
MEDCouplingRemapper.hxx MEDCouplingExtrudedMesh.hxx MEDCouplingGaussLocalization.hxx              \
MEDCouplingAutoRefCountObjectPtr.hxx

# Libraries targets

dist_libmedcoupling_la_SOURCES = \
	MEDCouplingField.cxx  MEDCouplingFieldDouble.cxx                                    \
	MEDCouplingUMesh.cxx  MEDCouplingMemArray.cxx  MEDCouplingTimeLabel.cxx             \
	MEDCouplingCMesh.cxx MEDCouplingTimeDiscretization.cxx                              \
	MEDCouplingFieldDiscretization.cxx MEDCouplingRefCountObject.cxx                    \
	MEDCouplingPointSet.cxx MEDCouplingUMeshDesc.cxx                                    \
	MEDCouplingExtrudedMesh.cxx MEDCouplingMesh.cxx MEDCouplingGaussLocalization.cxx    \
	MEDCouplingNatureOfField.cxx

libmedcoupling_la_LDFLAGS= 

libmedcoupling_la_CPPFLAGS=  @CXXTMPDPTHFLAGS@ \
	-I$(srcdir) -I$(srcdir)/../INTERP_KERNEL/Bases \
	-I$(srcdir)/../INTERP_KERNEL -I$(srcdir)/../INTERP_KERNEL/Geometric2D \
	-I$(srcdir)/../INTERP_KERNEL/ExprEval

# the geom2D library is included in the interpkernel one
libmedcoupling_la_LIBADD= ../INTERP_KERNEL/libinterpkernel.la

dist_libmedcouplingremapper_la_SOURCES = MEDCouplingRemapper.cxx

libmedcouplingremapper_la_CPPFLAGS =  @CXXTMPDPTHFLAGS@ \
	-I$(srcdir) -I$(srcdir)/../INTERP_KERNEL/Bases -I$(srcdir)/../INTERP_KERNEL \
	-I$(srcdir)/../INTERP_KERNEL/Geometric2D

libmedcouplingremapper_la_LDFLAGS = 

libmedcouplingremapper_la_LIBADD = libmedcoupling.la

AM_CPPFLAGS += $(libinterpkernel_la_CPPFLAGS)
LDADD= $(libinterpkernel_la_LDFLAGS)

EXTRA_DIST += \
	MEDCouplingFieldDouble.hxx                \
	MEDCouplingTimeDiscretization.hxx         \
	MEDCouplingFieldDiscretization.hxx        \
        MEDCouplingMesh.hxx                       \
        MEDCouplingUMesh.hxx                      \
        MEDCouplingCMesh.hxx                      \
        MEDCouplingTimeLabel.hxx                  \
	MEDCouplingField.hxx                      \
	MEDCouplingNormalizedUnstructuredMesh.hxx \
	MEDCouplingMemArray.hxx                   \
	MEDCouplingNormalizedUnstructuredMesh.txx \
	MEDCouplingMemArray.txx                   \
	MEDCouplingRefCountObject.hxx             \
	MEDCouplingPointSet.hxx                   \
	MEDCouplingPointSet.txx                   \
	MEDCouplingUMeshDesc.hxx                  \
	MEDCouplingNatureOfField.hxx              \
	MEDCouplingRemapper.hxx                   \
	MEDCouplingExtrudedMesh.hxx               \
	MEDCouplingGaussLocalization.hxx 
