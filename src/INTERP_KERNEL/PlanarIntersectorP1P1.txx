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
#ifndef __PLANARINTERSECTORP1P1_TXX__
#define __PLANARINTERSECTORP1P1_TXX__

#include "PlanarIntersectorP1P1.hxx"
#include "InterpolationUtils.hxx"
#include "CellModel.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, class ConcreteP1P1Intersector>
  PlanarIntersectorP1P1<MyMeshType,MyMatrix,ConcreteP1P1Intersector>::PlanarIntersectorP1P1(const MyMeshType& meshT, const MyMeshType& meshS,
                                                                                            double dimCaracteristic, double precision, double md3DSurf, double minDot3DSurf, double medianPlane,
                                                                                            bool doRotate, int orientation, int printLevel):
    PlanarIntersector<MyMeshType,MyMatrix>(meshT,meshS,dimCaracteristic,precision,md3DSurf,minDot3DSurf,medianPlane,doRotate,orientation,printLevel)
  {
  }

  template<class MyMeshType, class MyMatrix, class ConcreteP1P1Intersector>
  int PlanarIntersectorP1P1<MyMeshType,MyMatrix,ConcreteP1P1Intersector>::getNumberOfRowsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getNumberOfNodes();
  }

  template<class MyMeshType, class MyMatrix, class ConcreteP1P1Intersector>
  int PlanarIntersectorP1P1<MyMeshType,MyMatrix,ConcreteP1P1Intersector>::getNumberOfColsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getNumberOfNodes();
  }

  /*!
   * This methods split on the fly, into triangles in order to compute dual mesh of target cell (with icellT id in target mesh in C mode).
   */
  template<class MyMeshType, class MyMatrix, class ConcreteP1P1Intersector>
  void PlanarIntersectorP1P1<MyMeshType,MyMatrix,ConcreteP1P1Intersector>::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    int nbNodesT=PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[icellT+1]-PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[icellT];
    int orientation=1;
    const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectT+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[icellT]);
    std::vector<double> polygT;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(icellT),polygT);
    for(int nodeIdT=0;nodeIdT<nbNodesT;nodeIdT++)
      {
        ConnType curNodeTInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[nodeIdT]);
        PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinatesPermute(OTT<ConnType,numPol>::indFC(icellT),nodeIdT,polygT);
        std::vector<double> polygDualT(SPACEDIM*2*(nbNodesT-1));
        fillDualCellOfPolyg<SPACEDIM>(&polygT[0],polygT.size()/SPACEDIM,&polygDualT[0]);
        typename MyMatrix::value_type& resRow=res[curNodeTInCmode];
        for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++)
          {
            int iS=*iter;
            int nbNodesS=PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[iS+1]-PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[iS];
            const ConnType *startOfCellNodeConnS=PlanarIntersector<MyMeshType,MyMatrix>::_connectS+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[iS]);
            for(int nodeIdS=0;nodeIdS<nbNodesS;nodeIdS++)
              {
                ConnType curNodeSInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConnS[nodeIdS]);
                std::vector<double> polygS;
                PlanarIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinatesPermute(OTT<ConnType,numPol>::indFC(iS),nodeIdS,polygS);
                std::vector<double> polygDualS(SPACEDIM*2*(nbNodesS-1));
                fillDualCellOfPolyg<SPACEDIM>(&polygS[0],polygS.size()/SPACEDIM,&polygDualS[0]);
                std::vector<double> polygDualTTmp(polygDualT);
                if(SPACEDIM==3)
                  orientation=PlanarIntersector<MyMeshType,MyMatrix>::projectionThis(&polygDualS[0],&polygDualTTmp[0],polygDualS.size()/SPACEDIM,polygDualT.size()/SPACEDIM);
                double surf=orientation*intersectGeometryGeneral(polygDualTTmp,polygDualS);
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
