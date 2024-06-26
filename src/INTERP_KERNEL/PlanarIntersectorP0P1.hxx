// Copyright (C) 2007-2024  CEA, EDF
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

#ifndef __PLANARINTERSECTORP0P1_HXX__
#define __PLANARINTERSECTORP0P1_HXX__

#include "PlanarIntersector.hxx"
#include "InterpKernelUtilities.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, class ConcreteP0P1Intersector>
  class PlanarIntersectorP0P1 : public PlanarIntersector<MyMeshType,MyMatrix>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  protected:
    PlanarIntersectorP0P1(const MyMeshType& meshT, const MyMeshType& meshS, double dimCaracteristic, double precision, double md3DSurf, double minDot3DSurf, double medianPlane, bool doRotate, int orientation, int printLevel);
  public:
    void intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res);
    ConnType getNumberOfRowsOfResMatrix() const;
    ConnType getNumberOfColsOfResMatrix() const;
    /*!
     * Contrary to intersectCells method here icellS and icellT are \b not in \b C mode but in mode of MyMeshType.
     */
    double intersectGeometryWithQuadrangle(const double *quadrangle, const std::vector<double>& sourceCoords, bool isSourceQuad) { return asLeaf().intersectGeometryWithQuadrangle(quadrangle,sourceCoords,isSourceQuad); }
  protected:
    ConcreteP0P1Intersector& asLeaf() { return static_cast<ConcreteP0P1Intersector&>(*this); }
  };
}

#endif
