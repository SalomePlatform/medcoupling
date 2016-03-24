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

#ifndef __PLANARINTERSECTORP0P0_HXX__
#define __PLANARINTERSECTORP0P0_HXX__

#include "PlanarIntersector.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, class ConcreteP0P0Intersector>
  class PlanarIntersectorP0P0 : public PlanarIntersector<MyMeshType,MyMatrix>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  protected:
    PlanarIntersectorP0P0(const MyMeshType& meshT, const MyMeshType& meshS, double dimCaracteristic, double precision, double md3DSurf, double minDot3DSurf, double medianPlane, bool doRotate, int orientation, int printLevel);
  public:
    int getNumberOfRowsOfResMatrix() const;
    int getNumberOfColsOfResMatrix() const;
    void intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res);
    /*!
     * Contrary to intersectCells method here icellS and icellT are \b not in \b C mode but in mode of MyMeshType.
     */
    double intersectGeometry(ConnType icellT, ConnType icellS, ConnType nbNodesT, ConnType nbNodesS) { return asLeaf().intersectGeometry(icellT,icellS,nbNodesT,nbNodesS); }
  protected:
    ConcreteP0P0Intersector& asLeaf() { return static_cast<ConcreteP0P0Intersector&>(*this); }
  };
}

#endif
