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

#ifndef __CONVEXINTERSECTOR_HXX__
#define __CONVEXINTERSECTOR_HXX__

#include "PlanarIntersectorP0P0.hxx"
#include "PlanarIntersectorP0P1.hxx"
#include "PlanarIntersectorP1P0.hxx"
#include "PlanarIntersectorP1P1.hxx"
#include "PlanarIntersectorP1P0Bary.hxx"
#include "PlanarIntersectorP0P1Bary.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType >
  class ConvexIntersector : public InterpType<MyMeshType,MyMatrix,ConvexIntersector<MyMeshType,MyMatrix,InterpType> >
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    ConvexIntersector(const MyMeshType& meshT, const MyMeshType& meshS, 
                      double dimCaracteristic, double precision, double md3DSurf, double minDot3DSurf, double medianPlane,
                      bool doRotate, int orientation, int printLevel);
    double intersectGeometry(ConnType icellT, ConnType icellS, ConnType nbNodesT, ConnType nbNodesS);
    double intersectGeometryWithQuadrangle(const double *quadrangle, const std::vector<double>& sourceCoords, bool isSourceQuad);
    double intersectGeometryGeneral(const std::vector<double>& targetCoords, const std::vector<double>& sourceCoords);
    double intersectGeoBary(const std::vector<double>& targetCell, bool targetCellQuadratic, const double *sourceCell, std::vector<double>& res);
  private :
    double _epsilon;
  };
}

#endif
