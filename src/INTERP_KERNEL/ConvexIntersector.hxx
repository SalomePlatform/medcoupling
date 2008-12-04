//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef __CONVEXINTERSECTOR_HXX__
#define __CONVEXINTERSECTOR_HXX__

#include "PlanarIntersector.hxx"
#include "InterpolationUtils.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType>
  class INTERPKERNEL_EXPORT ConvexIntersector : public PlanarIntersector<MyMeshType>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    ConvexIntersector(const MyMeshType& mesh_A, const MyMeshType& mesh_B, 
                      double dimCaracteristic, double precision, double medianPlane,
                      bool doRotate, int printLevel);
    double intersectCells(ConnType icell_A, ConnType icell_B, int nb_NodesA, int nb_NodesB);
  private :
    double _epsilon;
    const ConnType *_connectA;
    const ConnType *_connectB;
    const ConnType *_connIndexA;
    const ConnType *_connIndexB;
    const double *_coordsA;
    const double *_coordsB;
  };
}

#endif
