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
#ifndef __PLANARINTERSECTORP1P0_TXX__
#define __PLANARINTERSECTORP1P0_TXX__

#include "PlanarIntersectorP1P0.hxx"
#include "InterpolationUtils.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, class ConcreteP1P0Intersector>
  PlanarIntersectorP1P0<MyMeshType,MyMatrix,ConcreteP1P0Intersector>::PlanarIntersectorP1P0(const MyMeshType& meshT, const MyMeshType& meshS,
                                                                                            double dimCaracteristic, double precision, double md3DSurf, double minDot3DSurf, double medianPlane,
                                                                                            bool doRotate, int orientation, int printLevel):
    PlanarIntersector<MyMeshType,MyMatrix>(meshT,meshS,dimCaracteristic,precision,md3DSurf,minDot3DSurf,medianPlane,doRotate,orientation,printLevel)
  {
  }

  template<class MyMeshType, class MyMatrix, class ConcreteP1P0Intersector>
  int PlanarIntersectorP1P0<MyMeshType,MyMatrix,ConcreteP1P0Intersector>::getNumberOfRowsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getNumberOfElements();
  }

  template<class MyMeshType, class MyMatrix, class ConcreteP1P0Intersector>
  int PlanarIntersectorP1P0<MyMeshType,MyMatrix,ConcreteP1P0Intersector>::getNumberOfColsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getNumberOfNodes();
  }

  /*!
   * This methods split on the fly, into triangles in order to compute dual mesh of target cell (with icellT id in target mesh in C mode).
   */
  template<class MyMeshType, class MyMatrix, class ConcreteP1P0Intersector>
  void PlanarIntersectorP1P0<MyMeshType,MyMatrix,ConcreteP1P0Intersector>::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    double triangle[9];
    double quadrangle[12];
    std::vector<double> targetCellCoords;
    int orientation=1;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(icellT),targetCellCoords);
    NormalizedCellType tT=PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getTypeOfElement(OTT<ConnType,numPol>::indFC(icellT));
    bool isTargetQuad=CellModel::GetCellModel(tT).isQuadratic();
    typename MyMatrix::value_type& resRow=res[icellT];
    for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++)
      {
        int iS=*iter;
        int nbNodesS=PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[iS+1]-PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[iS];
        const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectS+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[iS]);
        for(int nodeIdS=0;nodeIdS<nbNodesS;nodeIdS++)
          {
            ConnType curNodeSInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[nodeIdS]);
            std::copy(PlanarIntersector<MyMeshType,MyMatrix>::_coordsS+curNodeSInCmode*SPACEDIM,
                      PlanarIntersector<MyMeshType,MyMatrix>::_coordsS+curNodeSInCmode*SPACEDIM+SPACEDIM,triangle);
            for(int subTriS=1;subTriS<=nbNodesS-2;subTriS++)
              {
                std::copy(PlanarIntersector<MyMeshType,MyMatrix>::_coordsS+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[(nodeIdS+subTriS)%nbNodesS])*SPACEDIM,
                          PlanarIntersector<MyMeshType,MyMatrix>::_coordsS+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[(nodeIdS+subTriS)%nbNodesS])*SPACEDIM+SPACEDIM,
                          triangle+SPACEDIM);
                std::copy(PlanarIntersector<MyMeshType,MyMatrix>::_coordsS+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[(nodeIdS+subTriS+1)%nbNodesS])*SPACEDIM,
                          PlanarIntersector<MyMeshType,MyMatrix>::_coordsS+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[(nodeIdS+subTriS+1)%nbNodesS])*SPACEDIM+SPACEDIM,
                          triangle+2*SPACEDIM);
                fillDualCellOfTri<SPACEDIM>(triangle,quadrangle);
                std::vector<double> targetCellCoordsTmp(targetCellCoords);
                if(SPACEDIM==3)
                  orientation=PlanarIntersector<MyMeshType,MyMatrix>::projectionThis(&targetCellCoordsTmp[0],quadrangle,targetCellCoords.size()/SPACEDIM,4);
                double surf=orientation*intersectGeometryWithQuadrangle(quadrangle,targetCellCoordsTmp,isTargetQuad);
                surf=PlanarIntersector<MyMeshType,MyMatrix>::getValueRegardingOption(surf);
                if(surf!=0.)
                  {
                    typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(curNodeSInCmode));
                    if(iterRes==resRow.end())
                      resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(curNodeSInCmode),surf));
                    else
                      {
                        double val=(*iterRes).second+surf;
                        resRow.erase(OTT<ConnType,numPol>::indFC(curNodeSInCmode));
                        resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(curNodeSInCmode),val));
                      }
                  }
              }
          }
      }
  }
}

#endif
