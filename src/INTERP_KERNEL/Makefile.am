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
#  MED files in memory
#  File   : Makefile.am
#  Author : Vincent BERGEAUD (CEA/DEN/DANS/DM2S/SFME/LGLS)
#  Module : MED
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

SUBDIRS = Bases Geometric2D .

DIST_SUBDIRS = Bases Geometric2D

lib_LTLIBRARIES = libinterpkernel.la


salomeinclude_HEADERS = \
INTERPKERNELDefines.hxx \
BoundingBox.hxx                 PolyhedronIntersector.hxx          RegionNode.hxx\
ConvexIntersector.hxx           SplitterTetra.hxx                  Intersector3DP0P1.txx \
Geometric2DIntersector.hxx      Log.hxx                            TargetIntersector.hxx\
CellModel.hxx                   TetraAffineTransform.hxx           Intersector3DP1P0.hxx \
InterpKernelMatrix.hxx          MeshElement.hxx                    TransformedTriangle.hxx\
Interpolation2D.hxx             MeshRegion.hxx                     TransformedTriangleInline.hxx\
Interpolation3D.hxx             MeshUtils.hxx                      TranslationRotationMatrix.hxx\
Interpolation3DSurf.hxx         MeshRegion.txx                     TriangulationIntersector.hxx\
Interpolation.hxx               PlanarIntersector.hxx              VectorUtils.hxx\
InterpolationPlanar.hxx         PolyhedronIntersectorP1P0.txx      VTKNormalizedUnstructuredMesh.hxx\
InterpolationUtils.hxx          PolygonAlgorithms.hxx              InterpolationOptions.hxx\
BBTree.txx                      InterpolationPlanar.txx            PlanarIntersector.txx\
ConvexIntersector.txx           PolyhedronIntersector.txx          Intersector3DP1P0.txx\
Geometric2DIntersector.txx      SplitterTetra.txx                  PolygonAlgorithms.txx\
Interpolation2D.txx             TriangulationIntersector.txx       PolyhedronIntersectorP1P0.hxx\
Interpolation3DSurf.txx         MeshElement.txx                    VTKNormalizedUnstructuredMesh.txx\
Interpolation3D.txx             PlanarIntersectorP0P0.hxx          PlanarIntersectorP0P0.txx\
PolyhedronIntersectorP0P1.hxx   PolyhedronIntersectorP0P1.txx      Intersector3D.hxx\
Intersector3D.txx               Intersector3DP0P1.hxx 


# Libraries targets

dist_libinterpkernel_la_SOURCES = \
	TransformedTriangle.cxx\
	TransformedTriangleIntersect.cxx\
	TransformedTriangleMath.cxx\
	BoundingBox.cxx \
	TetraAffineTransform.cxx\
	CellModel.cxx\
	UnitTetraIntersectionBary.cxx

libinterpkernel_la_CPPFLAGS=-I$(srcdir)/Geometric2D -I$(srcdir)/Bases

libinterpkernel_la_LDFLAGS=

# the geom2D library is included in the interpkernel one
libinterpkernel_la_LIBADD= ./Geometric2D/libInterpGeometric2DAlg.la Bases/libinterpkernelbases.la

AM_CPPFLAGS= $(libinterpkernel_la_CPPFLAGS)
LDADD= $(libinterpkernel_la_LDFLAGS)

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
	Intersector3D.hxx \
	Intersector3D.txx \
	Intersector3DP0P1.hxx \
	Intersector3DP0P1.txx \
	Intersector3DP1P0.hxx \
	Intersector3DP1P0.txx \
	PolyhedronIntersector.hxx \
	PolyhedronIntersector.txx \
	PolyhedronIntersectorP0P1.hxx \
	PolyhedronIntersectorP0P1.txx \
	SplitterTetra.hxx \
	SplitterTetra.txx \
	Log.hxx \
	MEDNormalizedUnstructuredMesh.hxx \
	MEDNormalizedUnstructuredMesh.txx \
	MeshElement.hxx \
	MeshElement.txx \
	MeshRegion.hxx \
	MeshRegion.txx \
	MeshUtils.hxx \
	PlanarIntersector.hxx \
	PlanarIntersector.txx \
	PlanarIntersectorP0P0.hxx \
	PlanarIntersectorP0P0.txx \
	PointLocatorAlgos.txx \
	PolygonAlgorithms.hxx \
	PolygonAlgorithms.txx \
	RegionNode.hxx \
	TargetIntersector.hxx \
	TestInterpKernel.cpp \
	TetraAffineTransform.hxx \
	TransformedTriangle.hxx \
	TransformedTriangleInline.hxx \
	TranslationRotationMatrix.hxx \
	TriangulationIntersector.hxx \
	TriangulationIntersector.txx \
	VTKNormalizedUnstructuredMesh.hxx \
	VTKNormalizedUnstructuredMesh.txx \
	VectorUtils.hxx
