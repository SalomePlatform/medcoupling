# Copyright (C) 2012-2016  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
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

SET(interpkernel_SOURCES
  TransformedTriangle.cxx
  TransformedTriangleIntersect.cxx
  TransformedTriangleMath.cxx
  BoundingBox.cxx
  TranslationRotationMatrix.cxx
  TetraAffineTransform.cxx
  CellModel.cxx
  DiameterCalculator.cxx
  OrientationInverter.cxx
  UnitTetraIntersectionBary.cxx
  InterpolationOptions.cxx
  BoxSplittingOptions.cxx
  DirectedBoundingBox.cxx
  Interpolation2DCurve.cxx
  Interpolation3DSurf.cxx
  Interpolation3D.cxx
  Interpolation2D3D.cxx
  Interpolation3D1D.cxx
  MeshElement.cxx
  InterpKernelMeshQuality.cxx
  InterpKernelCellSimplify.cxx
  InterpKernelMatrixTools.cxx
  VolSurfUser.cxx
  SplitterTetra.cxx
  Bases/InterpKernelException.cxx
  Geometric2D/InterpKernelGeo2DAbstractEdge.cxx
  Geometric2D/InterpKernelGeo2DBounds.cxx
  Geometric2D/InterpKernelGeo2DPrecision.cxx
  Geometric2D/InterpKernelGeo2DComposedEdge.cxx
  Geometric2D/InterpKernelGeo2DEdgeArcCircle.cxx
  Geometric2D/InterpKernelGeo2DEdge.cxx
  Geometric2D/InterpKernelGeo2DEdgeInfLin.cxx
  Geometric2D/InterpKernelGeo2DEdgeLin.cxx
  Geometric2D/InterpKernelGeo2DElementaryEdge.cxx
  Geometric2D/InterpKernelGeo2DNode.cxx
  Geometric2D/InterpKernelGeo2DQuadraticPolygon.cxx
  ExprEval/InterpKernelExprParser.cxx
  ExprEval/InterpKernelFunction.cxx
  ExprEval/InterpKernelUnit.cxx
  ExprEval/InterpKernelValue.cxx
  ExprEval/InterpKernelAsmX86.cxx
  GaussPoints/InterpKernelGaussCoords.cxx
  )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/Bases
  ${CMAKE_CURRENT_SOURCE_DIR}/Geometric2D
  ${CMAKE_CURRENT_SOURCE_DIR}/ExprEval
  ${CMAKE_CURRENT_SOURCE_DIR}/GaussPoints
  )

SET(PLATFORM_MMAP)
IF(NOT WIN32)
  SET(PLATFORM_MMAP "-D_POSIX_MAPPED_FILES")
ENDIF(NOT WIN32)

ADD_LIBRARY(interpkernel ${interpkernel_SOURCES})
SET_TARGET_PROPERTIES(interpkernel PROPERTIES COMPILE_FLAGS "${PLATFORM_MMAP}")
TARGET_LINK_LIBRARIES(interpkernel ${PLATFORM_LIBS})
INSTALL(TARGETS interpkernel EXPORT ${PROJECT_NAME}TargetGroup DESTINATION ${MEDCOUPLING_INSTALL_LIBS})

FILE(GLOB_RECURSE interpkernel_HEADERS_HXX "${CMAKE_CURRENT_SOURCE_DIR}/*.hxx")
FILE(GLOB_RECURSE interpkernel_HEADERS_TXX "${CMAKE_CURRENT_SOURCE_DIR}/*.txx")
INSTALL(FILES ${interpkernel_HEADERS_HXX} ${interpkernel_HEADERS_TXX} Bases/NormalizedGeometricTypes DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

# Will be used for SWIG dependencies:
SET (interpkernel_HEADERS_HXX PARENT_SCOPE)
SET (interpkernel_HEADERS_TXX PARENT_SCOPE)
