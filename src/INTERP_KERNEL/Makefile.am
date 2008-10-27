#  MED MEDMEM : MED files in memory
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
#
#
#  File   : Makefile.am
#  Author : Vincent BERGEAUD (CEA/DEN/DANS/DM2S/SFME/LGLS)
#  Module : MED

include $(top_srcdir)/adm_local/unix/make_common_starter.am

# first '.' says that this folder must be compiled before Test (MEDMEM/Test uses MEDMEM)
if CPPUNIT_IS_OK
 SUBDIRS = Geometric2D . Test 
else
 SUBDIRS = Geometric2D .
endif

DIST_SUBDIRS = Geometric2D Test

lib_LTLIBRARIES = libinterpkernel.la


salomeinclude_HEADERS = \
BoundingBox.hxx             IntersectorHexa.hxx                RegionNode.hxx\
ConvexIntersector.hxx       IntersectorTetra.hxx               Remapper.hxx\
Geometric2DIntersector.hxx  Log.hxx                            TargetIntersector.hxx\
INTERPKERNEL_defines.hxx    MEDNormalizedUnstructuredMesh.hxx  TetraAffineTransform.hxx\
InterpKernelMatrix.hxx      MeshElement.hxx                    TransformedTriangle.hxx\
Interpolation2D.hxx         MeshRegion.hxx                     TransformedTriangle_inline.hxx\
Interpolation3D.hxx         MeshUtils.hxx                      TranslationRotationMatrix.hxx\
Interpolation3DSurf.hxx     NormalizedUnstructuredMesh.hxx     TriangulationIntersector.hxx\
Interpolation.hxx           PlanarIntersector.hxx              VectorUtils.hxx\
InterpolationPlanar.hxx     PointLocator.hxx                   VTKNormalizedUnstructuredMesh.hxx\
InterpolationUtils.hxx      PolygonAlgorithms.hxx\
BBTree.txx                  InterpolationPlanar.txx            PlanarIntersector.txx\
ConvexIntersector.txx       IntersectorHexa.txx                PointLocatorAlgos.txx\
Geometric2DIntersector.txx  IntersectorTetra.txx               PolygonAlgorithms.txx\
Interpolation2D.txx         MEDNormalizedUnstructuredMesh.txx  TriangulationIntersector.txx\
Interpolation3DSurf.txx     MeshElement.txx                    VTKNormalizedUnstructuredMesh.txx\
Interpolation3D.txx         MeshRegion.txx

# Libraries targets

dist_libinterpkernel_la_SOURCES = \
	TransformedTriangle.cxx\
	TransformedTriangle_intersect.cxx\
	TransformedTriangle_math.cxx\
	BoundingBox.cxx \
	TetraAffineTransform.cxx\
	CellModel.cxx \
	Remapper.cxx \
	PointLocator.cxx

libinterpkernel_la_LDFLAGS= $(MED2_LIBS) $(HDF5_LIBS) ../MEDWrapper/V2_1/Core/libmed_V2_1.la $(STDLIB) ../MEDMEM/libmedmem.la -lutil -lm

libinterpkernel_la_CPPFLAGS=  @CXXTMPDPTHFLAGS@ \
	-I$(srcdir) -I$(srcdir)/../MEDMEM $(BOOST_CPPFLAGS) -DOPTIMIZE -DLOG_LEVEL=0 \
	-I$(srcdir)/../MEDWrapper/Base \
	-I$(srcdir)/../MEDWrapper/V2_1/Core \
	-I$(srcdir)/Geometric2D \
	$(MED2_INCLUDES)

# the geom2D library is included in the interpkernel one
libinterpkernel_la_LIBADD= ./Geometric2D/libInterpGeometric2DAlg.la

if MED_ENABLE_KERNEL
    libinterpkernel_la_CPPFLAGS+= ${KERNEL_CXXFLAGS}
    libinterpkernel_la_LDFLAGS+= ${KERNEL_LDFLAGS} -lSALOMELocalTrace 
endif

AM_CPPFLAGS= $(libinterpkernel_la_CPPFLAGS)
LDADD= $(libinterpkernel_la_LDFLAGS)
if MED_ENABLE_KERNEL
     LDADD+=-lSALOMEBasics
endif

EXTRA_DIST += \
	BBTree.txx \
	BoundingBox.hxx \
	ConvexIntersector.hxx \
	ConvexIntersector.txx \
	Geometric2DIntersector.hxx \
	Geometric2DIntersector.txx \
	Geometric2DIntersector.hxx \
	Geometric2DIntersector.txx \
	InterpKernelMatrix.hxx \
	Interpolation2D.txx \
	Interpolation3D.txx \
	Interpolation3DSurf.txx \
	InterpolationPlanar.hxx \
	InterpolationPlanar.txx \
	InterpolationUtils.hxx \
	IntersectorHexa.hxx \
	IntersectorHexa.txx \
	IntersectorTetra.hxx \
	IntersectorTetra.txx \
	Log.hxx \
	MEDNormalizedUnstructuredMesh.hxx \
	MEDNormalizedUnstructuredMesh.txx \
	MeshElement.hxx \
	MeshElement.txx \
	MeshRegion.hxx \
	MeshRegion.txx \
	MeshUtils.hxx \
	NormalizedUnstructuredMesh.hxx \
	PlanarIntersector.hxx \
	PlanarIntersector.txx \
	PointLocatorAlgos.txx \
	PolygonAlgorithms.hxx \
	PolygonAlgorithms.txx \
	RegionNode.hxx \
	TargetIntersector.hxx \
	TestInterpKernel.cpp \
	TetraAffineTransform.hxx \
	TransformedTriangle.hxx \
	TransformedTriangle_inline.hxx \
	TranslationRotationMatrix.hxx \
	TriangulationIntersector.hxx \
	TriangulationIntersector.txx \
	VTKNormalizedUnstructuredMesh.hxx \
	VTKNormalizedUnstructuredMesh.txx \
	VectorUtils.hxx
