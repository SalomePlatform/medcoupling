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
#ifndef __PLANARINTERSECTORP1P0PL_TXX__
#define __PLANARINTERSECTORP1P0PL_TXX__

#include "PlanarIntersectorP1P0PL.hxx"
#include "PlanarIntersector.txx"
#include "CellModel.hxx"

#include "PointLocatorAlgos.txx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "MeshUtils.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  PlanarIntersectorP1P0PL<MyMeshType,MyMatrix>::PlanarIntersectorP1P0PL(const MyMeshType& meshT, const MyMeshType& meshS,
                                                                        double dimCaracteristic, double md3DSurf, double minDot3DSurf,
                                                                        double medianPlane, double precision, int orientation):
    PlanarIntersector<MyMeshType,MyMatrix>(meshT,meshS,dimCaracteristic,precision,md3DSurf,minDot3DSurf,medianPlane,true,orientation,0)
  {
  }

  template<class MyMeshType, class MyMatrix>
  void PlanarIntersectorP1P0PL<MyMeshType,MyMatrix>::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    std::vector<double> CoordsT;
    typename MyMatrix::value_type& resRow=res[icellT];
    PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(icellT),CoordsT);
    double baryT[SPACEDIM];
    double baryTTmp[SPACEDIM];
    calculateBarycenterDyn2<SPACEDIM>(&CoordsT[0],CoordsT.size()/SPACEDIM,baryT);
    for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++)
      {
        NormalizedCellType tS=PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getTypeOfElement(OTT<ConnType,numPol>::indFC(*iter));
        if(tS!=NORM_TRI3)
          throw INTERP_KERNEL::Exception("Invalid source cell detected for meshdim==2. Only TRI3 supported !");
        std::vector<double> CoordsS;
        PlanarIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinates(OTT<ConnType,numPol>::indFC(*iter),CoordsS);
        if(SPACEDIM==2)
          {
            std::copy(baryT,baryT+SPACEDIM,baryTTmp);
          }
        else
          {
            double littleTargetCell[9];
            std::copy(baryT,baryT+SPACEDIM,littleTargetCell);
            std::copy(CoordsT.begin(),CoordsT.begin()+3,littleTargetCell+3);
            std::copy(CoordsT.begin()+3,CoordsT.begin()+6,littleTargetCell+6);
            PlanarIntersector<MyMeshType,MyMatrix>::projectionThis(&CoordsS[0],littleTargetCell,3,3);
            std::copy(littleTargetCell,littleTargetCell+3,baryTTmp);
          }
        if(PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg2D(baryTTmp,&CoordsS[0],3,PlanarIntersector<MyMeshType,MyMatrix>::_precision))
          {
            double resLoc[3];
            barycentric_coords<SPACEDIM>(&CoordsS[0],baryTTmp,resLoc);
            const ConnType *startOfCellNodeConnS=PlanarIntersector<MyMeshType,MyMatrix>::_connectS+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[*iter]);
            for(int nodeIdS=0;nodeIdS<3;nodeIdS++)
              {
                if(fabs(resLoc[nodeIdS])>PlanarIntersector<MyMeshType,MyMatrix>::_precision)
                  {
                    ConnType curNodeSInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConnS[nodeIdS]);
                    typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(curNodeSInCmode));
                    if(iterRes==resRow.end())
                      resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(curNodeSInCmode),resLoc[nodeIdS]));
                    else
                      {
                        double val=(*iterRes).second+resLoc[nodeIdS];
                        resRow.erase(OTT<ConnType,numPol>::indFC(curNodeSInCmode));
                        resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(curNodeSInCmode),val));
                      }
                  }
              }
          }
      }
  }

  template<class MyMeshType, class MyMatrix>
  int PlanarIntersectorP1P0PL<MyMeshType,MyMatrix>::getNumberOfRowsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getNumberOfElements();
  }
  
  template<class MyMeshType, class MyMatrix>
  int PlanarIntersectorP1P0PL<MyMeshType,MyMatrix>::getNumberOfColsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getNumberOfNodes();
  }
}

#endif
