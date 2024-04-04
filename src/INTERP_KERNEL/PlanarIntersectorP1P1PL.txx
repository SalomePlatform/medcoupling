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
#ifndef __PLANARINTERSECTORP1P1PL_TXX__
#define __PLANARINTERSECTORP1P1PL_TXX__

#include "PlanarIntersector.txx"
#include "InterpolationUtils.hxx"
#include "NormalizedGeometricTypes"
#include "PlanarIntersectorP1P1PL.hxx"
#include "InterpKernelUtilities.hxx"

#include "PointLocatorAlgos.txx"
#include <vector>
#include <cmath>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  PlanarIntersectorP1P1PL<MyMeshType,MyMatrix>::PlanarIntersectorP1P1PL(const MyMeshType& meshT, const MyMeshType& meshS,
                                                                        double dimCaracteristic, double md3DSurf, double minDot3DSurf,
                                                                        double medianPlane, double precision, int orientation):
    PlanarIntersector<MyMeshType,MyMatrix>(meshT,meshS,dimCaracteristic,precision,md3DSurf,minDot3DSurf,medianPlane,true,orientation,0)
  {
  }

  template<class MyMeshType, class MyMatrix>
  void PlanarIntersectorP1P1PL<MyMeshType,MyMatrix>::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    std::vector<double> CoordsT;
    this->getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(icellT),CoordsT);
    ConnType nbOfNodesT=ToConnType(CoordsT.size())/SPACEDIM;
    for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++)
      {
        NormalizedCellType tS=this->_meshS.getTypeOfElement(OTT<ConnType,numPol>::indFC(*iter));
        if(tS!=NORM_TRI3)
          throw INTERP_KERNEL::Exception("Invalid source cell detected for meshdim==2. Only TRI3 supported !");
        std::vector<double> CoordsS;
        this->getRealSourceCoordinates(OTT<ConnType,numPol>::indFC(*iter),CoordsS);
        std::vector<double> CoordsTTmp(CoordsT);
        if(SPACEDIM==3)
          this->projectionThis(CoordsS.data(),CoordsTTmp.data(),ToConnType(CoordsS.size())/SPACEDIM,nbOfNodesT);
        const ConnType *startOfCellNodeConnT=this->_connectT+OTT<ConnType,numPol>::conn2C(this->_connIndexT[icellT]);
        for(int nodeIdT=0;nodeIdT<nbOfNodesT;nodeIdT++)
          {
            typename MyMatrix::value_type& resRow=res[OTT<ConnType,numPol>::ind2C(startOfCellNodeConnT[nodeIdT])];
            if( PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg2DSimple(CoordsTTmp.data()+nodeIdT*SPACEDIM,CoordsS.data(),3,this->_precision) )
              {
                double resLoc[3];
                barycentric_coords<SPACEDIM>(&CoordsS[0],&CoordsTTmp[nodeIdT*SPACEDIM],resLoc);
                const ConnType *startOfCellNodeConnS=this->_connectS+OTT<ConnType,numPol>::conn2C(this->_connIndexS[*iter]);
                for(int nodeIdS=0;nodeIdS<3;nodeIdS++)
                  {
                    if(fabs(resLoc[nodeIdS])>this->_precision)
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
  }

  template<class MyMeshType, class MyMatrix>
  typename MyMeshType::MyConnType PlanarIntersectorP1P1PL<MyMeshType,MyMatrix>::getNumberOfRowsOfResMatrix() const
  {
    return this->_meshT.getNumberOfNodes();
  }
  
  template<class MyMeshType, class MyMatrix>
  typename MyMeshType::MyConnType PlanarIntersectorP1P1PL<MyMeshType,MyMatrix>::getNumberOfColsOfResMatrix() const
  {
    return this->_meshS.getNumberOfNodes();
  }
}

#endif
