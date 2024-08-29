# Copyright (C) 2012-2024  CEA, EDF
#
# This library is free software; you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by the Free
# Software Foundation; either version 2.1 of the License, or (at your option)
# any later version.
#
# This library is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more
# details.
#
# You should have received a copy of the GNU Lesser General Public License along
# with this library; if not, write to the Free Software Foundation, Inc., 59
# Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email :
# webmaster.salome@opencascade.com
#
# Author : Anthony Geay (CEA/DEN)

set(INTERPKERNEL_SOURCES
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
    Interpolation1D0D.cxx
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
    LinearAlgebra/InterpKernelDenseMatrix.cxx
    LinearAlgebra/InterpKernelLUDecomp.cxx
    LinearAlgebra/InterpKernelQRDecomp.cxx)

add_library(interpkernel ${INTERPKERNEL_SOURCES})
add_library(MEDCoupling::interpkernel ALIAS interpkernel)

if(MEDCOUPLING_USE_64BIT_IDS)
  target_compile_definitions(interpkernel PUBLIC "MEDCOULPING_USE_64BIT_IDS")
endif()

target_include_directories(
  interpkernel
  PUBLIC
  "\$<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>" # MEDCouplingVersion.h
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/interp_kernel>"
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/interp_kernel/Bases>"
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/interp_kernel/Geometric2D>"
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/interp_kernel/ExprEval>"
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/interp_kernel/GaussPoints>"
    "\$<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include/interp_kernel/LinearAlgebra>"
)

if(NOT DEFINED MSVC)
  target_compile_options(interpkernel PRIVATE -Wsign-compare -Wconversion)
endif()

if(NOT WIN32)
  set_target_properties(interpkernel PROPERTIES COMPILE_FLAGS "-D_POSIX_MAPPED_FILES")
endif(NOT WIN32)

install(
  TARGETS interpkernel
  EXPORT ${CMAKE_PROJECT_NAME}TargetGroup
  DESTINATION ${MEDCOUPLING_INSTALL_LIBS}
  INCLUDES DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})

file(GLOB_RECURSE interpkernel_HEADERS "${PROJECT_SOURCE_DIR}/include/interp_kernel/*")
install(FILES ${interpkernel_HEADERS} DESTINATION ${MEDCOUPLING_INSTALL_HEADERS})
