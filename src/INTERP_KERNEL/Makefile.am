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
#

#  MED files in memory
#  File   : Makefile.am
#  Author : Vincent BERGEAUD (CEA/DEN/DANS/DM2S/SFME/LGLS)
#  Module : MED
#
include $(top_srcdir)/adm_local/unix/make_common_starter.am

EXTRA_DIST += CMakeLists.txt

lib_LTLIBRARIES = libinterpkernel.la

salomeinclude_HEADERS = 		\
BBTree.txx				\
BoundingBox.hxx				\
CellModel.hxx				\
ConvexIntersector.hxx			\
ConvexIntersector.txx			\
Geometric2DIntersector.hxx		\
Geometric2DIntersector.txx		\
PointLocator2DIntersector.hxx           \
PointLocator2DIntersector.txx           \
INTERPKERNELDefines.hxx			\
InterpKernelMatrix.hxx			\
InterpKernelMatrixTools.hxx		\
InterpKernelMeshQuality.hxx             \
InterpKernelCellSimplify.hxx            \
Interpolation.hxx			\
Interpolation.txx                       \
Interpolation2D.hxx			\
Interpolation2D.txx			\
Interpolation3D.hxx			\
Interpolation3D.txx			\
Interpolation2D1D.hxx   	        \
Interpolation2D1D.txx           	\
Interpolation3D2D.hxx			\
Interpolation3D2D.txx			\
Interpolation3DSurf.hxx			\
InterpolationOptions.hxx		\
InterpolationPlanar.hxx			\
InterpolationPlanar.txx			\
InterpolationCC.hxx                     \
InterpolationCC.txx                     \
InterpolationCU.hxx                     \
InterpolationCU.txx                     \
InterpolationUtils.hxx			\
IntersectorCU.hxx			\
IntersectorCU.txx			\
IntersectorCU1D.hxx			\
IntersectorCU1D.txx			\
IntersectorCU2D.hxx			\
IntersectorCU2D.txx			\
IntersectorCU3D.hxx			\
IntersectorCU3D.txx			\
Intersector3D.hxx			\
Intersector3D.txx			\
Intersector3DP0P0.hxx	          	\
Intersector3DP0P0.txx	          	\
Intersector3DP0P1.hxx			\
Intersector3DP0P1.txx			\
Intersector3DP1P0.hxx			\
Intersector3DP1P0.txx			\
Intersector3DP1P0Bary.hxx		\
Intersector3DP1P0Bary.txx		\
Intersector3DP1P1.hxx			\
Intersector3DP1P1.txx			\
Log.hxx					\
MeshElement.hxx				\
MeshElement.txx				\
MeshRegion.hxx				\
MeshRegion.txx				\
MeshUtils.hxx				\
PointLocatorAlgos.txx                   \
PlanarIntersector.hxx			\
PlanarIntersector.txx			\
Planar2D1DIntersectorP0P0.hxx           \
Planar2D1DIntersectorP0P0.txx           \
PlanarIntersectorP0P0.hxx		\
PlanarIntersectorP0P0.txx		\
PlanarIntersectorP0P1.hxx		\
PlanarIntersectorP0P1.txx		\
PlanarIntersectorP1P0.hxx		\
PlanarIntersectorP1P0.txx		\
PlanarIntersectorP1P0Bary.hxx		\
PlanarIntersectorP1P0Bary.txx		\
PlanarIntersectorP0P1PL.hxx		\
PlanarIntersectorP0P1PL.txx		\
PlanarIntersectorP1P0PL.hxx		\
PlanarIntersectorP1P0PL.txx		\
PlanarIntersectorP1P1PL.hxx             \
PlanarIntersectorP1P1PL.txx             \
PlanarIntersectorP1P1.hxx		\
PlanarIntersectorP1P1.txx		\
PolygonAlgorithms.hxx			\
PolygonAlgorithms.txx			\
Polyhedron3D2DIntersectorP0P0.hxx       \
Polyhedron3D2DIntersectorP0P0.txx       \
PolyhedronIntersectorP0P0.hxx		\
PolyhedronIntersectorP0P0.txx		\
PointLocator3DIntersectorP0P0.hxx       \
PointLocator3DIntersectorP0P0.txx       \
PolyhedronIntersectorP0P1.hxx		\
PolyhedronIntersectorP0P1.txx		\
PointLocator3DIntersectorP0P1.hxx       \
PointLocator3DIntersectorP0P1.txx       \
PolyhedronIntersectorP1P0.hxx		\
PolyhedronIntersectorP1P0.txx		\
PointLocator3DIntersectorP1P0.hxx       \
PointLocator3DIntersectorP1P0.txx       \
PolyhedronIntersectorP1P0Bary.hxx	\
PolyhedronIntersectorP1P0Bary.txx	\
PolyhedronIntersectorP1P1.hxx		\
PolyhedronIntersectorP1P1.txx		\
PointLocator3DIntersectorP1P1.hxx       \
PointLocator3DIntersectorP1P1.txx       \
RegionNode.hxx				\
SplitterTetra.hxx			\
SplitterTetra.txx			\
TargetIntersector.hxx			\
TetraAffineTransform.hxx		\
TransformedTriangle.hxx			\
TransformedTriangleInline.hxx		\
TranslationRotationMatrix.hxx		\
TriangulationIntersector.hxx		\
TriangulationIntersector.txx		\
IntegralUniformIntersector.hxx          \
IntegralUniformIntersector.txx          \
UnitTetraIntersectionBary.hxx		\
VTKNormalizedUnstructuredMesh.hxx	\
VTKNormalizedUnstructuredMesh.txx	\
VectorUtils.hxx                         \
VolSurfFormulae.hxx                     \
GenMathFormulae.hxx                     \
VolSurfUser.hxx                         \
VolSurfUser.txx                         \
CurveIntersector.hxx                    \
CurveIntersector.txx                    \
CurveIntersectorP0P0.hxx                \
CurveIntersectorP0P0.txx                \
CurveIntersectorP0P1.hxx                \
CurveIntersectorP0P1.txx                \
CurveIntersectorP1P0.hxx                \
CurveIntersectorP1P0.txx                \
CurveIntersectorP1P1.hxx                \
CurveIntersectorP1P1.txx                \
Interpolation1D.hxx	                \
Interpolation1D.txx	                \
Interpolation2DCurve.hxx                \
InterpolationCurve.hxx	                \
InterpolationCurve.txx			\
DirectedBoundingBox.hxx			\
Bases/InterpKernelException.hxx         \
Bases/NormalizedUnstructuredMesh.hxx    \
Bases/InterpKernelStlExt.hxx            \
Bases/InterpKernelHashMap.hxx           \
Bases/InterpKernelHashTable.hxx         \
Bases/InterpKernelHashFun.hxx           \
Bases/InterpKernelAutoPtr.hxx           \
Geometric2D/InterpKernelGeo2DAbstractEdge.hxx      \
Geometric2D/InterpKernelGeo2DBounds.hxx            \
Geometric2D/InterpKernelGeo2DPrecision.hxx         \
Geometric2D/InterpKernelGeo2DComposedEdge.hxx      \
Geometric2D/InterpKernelGeo2DEdgeArcCircle.hxx     \
Geometric2D/InterpKernelGeo2DEdge.hxx              \
Geometric2D/InterpKernelGeo2DEdgeInfLin.hxx        \
Geometric2D/InterpKernelGeo2DEdgeLin.hxx           \
Geometric2D/InterpKernelGeo2DElementaryEdge.hxx    \
Geometric2D/InterpKernelGeo2DNode.hxx              \
Geometric2D/InterpKernelGeo2DQuadraticPolygon.hxx  \
Geometric2D/InterpKernelGeo2DEdge.txx              \
ExprEval/InterpKernelExprParser.hxx              \
ExprEval/InterpKernelFunction.hxx                \
ExprEval/InterpKernelUnit.hxx                    \
ExprEval/InterpKernelValue.hxx                   \
ExprEval/InterpKernelAsmX86.hxx                  \
GaussPoints/InterpKernelGaussCoords.hxx


EXTRA_DIST +=			\
InterpKernelUtilities.hxx	\
Intersector3DP0P0.hxx		\
Intersector3DP0P0.txx		\
PlanarIntersectorP0P1.hxx	\
PlanarIntersectorP0P1.txx	\
PlanarIntersectorP1P0.hxx	\
PlanarIntersectorP1P0.txx	\
VolSurfFormulae.hxx             \
VolSurfUser.hxx                 \
VolSurfUser.txx                 \
ExprEval/InterpKernelExprParser.hxx              \
ExprEval/InterpKernelFunction.hxx                \
ExprEval/InterpKernelUnit.hxx                    \
ExprEval/InterpKernelValue.hxx

# Libraries targets

dist_libinterpkernel_la_SOURCES = \
	TransformedTriangle.cxx \
	TransformedTriangleIntersect.cxx \
	TransformedTriangleMath.cxx \
	BoundingBox.cxx \
	TranslationRotationMatrix.cxx \
	TetraAffineTransform.cxx \
	CellModel.cxx \
	UnitTetraIntersectionBary.cxx \
	InterpolationOptions.cxx \
	DirectedBoundingBox.cxx \
	Interpolation2DCurve.cxx \
	Interpolation3DSurf.cxx \
	Interpolation3D.cxx \
	Interpolation3D2D.cxx \
	MeshElement.cxx \
	InterpKernelMeshQuality.cxx \
	InterpKernelCellSimplify.cxx \
	InterpKernelMatrixTools.cxx \
	Bases/InterpKernelException.cxx \
	Geometric2D/InterpKernelGeo2DAbstractEdge.cxx         \
	Geometric2D/InterpKernelGeo2DBounds.cxx               \
	Geometric2D/InterpKernelGeo2DPrecision.cxx            \
	Geometric2D/InterpKernelGeo2DComposedEdge.cxx         \
	Geometric2D/InterpKernelGeo2DEdgeArcCircle.cxx        \
	Geometric2D/InterpKernelGeo2DEdge.cxx                 \
	Geometric2D/InterpKernelGeo2DEdgeInfLin.cxx           \
	Geometric2D/InterpKernelGeo2DEdgeLin.cxx              \
	Geometric2D/InterpKernelGeo2DElementaryEdge.cxx       \
	Geometric2D/InterpKernelGeo2DNode.cxx                 \
	Geometric2D/InterpKernelGeo2DQuadraticPolygon.cxx     \
	ExprEval/InterpKernelExprParser.cxx        \
	ExprEval/InterpKernelFunction.cxx          \
	ExprEval/InterpKernelUnit.cxx              \
	ExprEval/InterpKernelValue.cxx             \
	ExprEval/InterpKernelAsmX86.cxx            \
	GaussPoints/InterpKernelGaussCoords.cxx

libinterpkernel_la_CPPFLAGS= -I$(srcdir)/Bases -I$(srcdir)/Geometric2D -I$(srcdir)/ExprEval -I$(srcdir)/GaussPoints

AM_CPPFLAGS += $(libinterpkernel_la_CPPFLAGS)
