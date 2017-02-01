// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// Author : Anthony Geay (CEA/DEN)

#ifndef __VOLSURFUSER_HXX__
#define __VOLSURFUSER_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelException.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{
  template<class ConnType, NumberingPolicy numPolConn, int SPACEDIM>
  double computeVolSurfOfCell(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords);

  template<class ConnType, NumberingPolicy numPolConn>
  double computeVolSurfOfCell2(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, int spaceDim);

  template<class ConnType, NumberingPolicy numPolConn, int SPACEDIM>
  void computeBarycenter(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, double *res);

  template<class ConnType, NumberingPolicy numPolConn>
  void computeBarycenter2(NormalizedCellType type, const ConnType *connec, int lgth, const double *coords, int spaceDim, double *res);

  double INTERPKERNEL_EXPORT OrthoDistanceFromPtToPlaneInSpaceDim3(const double *p, const double *p1, const double *p2, const double *p3);

  double INTERPKERNEL_EXPORT SquareDistanceFromPtToSegInSpaceDim2(const double *pt, const double *pt0Seg2, const double *pt1Seg2, std::size_t &nbOfHint);

  double INTERPKERNEL_EXPORT DistanceFromPtToTriInSpaceDim3(const double *pt, const double *pt0Tri3, const double *pt1Tri3, const double *pt2Tri3);

  double INTERPKERNEL_EXPORT DistanceFromPtToPolygonInSpaceDim3(const double *pt, const int *connOfPolygonBg, const int *connOfPolygonEnd, const double *coords);

  bool ComputeRotTranslationMatrixToPut3PointsOnOXY(const double *pt0Tri3, const double *pt1Tri3, const double *pt2Tri3, double *matrix);
}

#endif
