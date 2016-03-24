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
#ifndef __PLANARINTERSECTORP0P1PL_TXX__
#define __PLANARINTERSECTORP0P1PL_TXX__

#include "PlanarIntersectorP0P1PL.hxx"
#include "PlanarIntersector.txx"
#include "CellModel.hxx"

#include "PointLocatorAlgos.txx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  PlanarIntersectorP0P1PL<MyMeshType,MyMatrix>::PlanarIntersectorP0P1PL(const MyMeshType& meshT, const MyMeshType& meshS,
                                                                        double dimCaracteristic, double md3DSurf, double minDot3DSurf,
                                                                        double medianPlane, double precision, int orientation):
    PlanarIntersector<MyMeshType,MyMatrix>(meshT,meshS,dimCaracteristic,precision,md3DSurf,minDot3DSurf,medianPlane,true,orientation,0)
  {
  }

  template<class MyMeshType, class MyMatrix>
  void PlanarIntersectorP0P1PL<MyMeshType,MyMatrix>::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    std::vector< std::vector<double> > coordsOfSources(icellsS.size());
    int ii=0;
    for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++,ii++)
      PlanarIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinates(OTT<ConnType,numPol>::indFC(*iter),coordsOfSources[ii]);
    const ConnType *startOfCellNodeConnT=PlanarIntersector<MyMeshType,MyMatrix>::_connectT+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[icellT]);
    std::vector<double> coordsTarget;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(icellT),coordsTarget);
    int nbNodesT=coordsTarget.size()/SPACEDIM;
    ii=0;
    for(typename std::vector<ConnType>::const_iterator iter2=icellsS.begin();iter2!=icellsS.end();iter2++,ii++)
      {
        std::vector<double> tmpSource(coordsOfSources[ii]);
        std::vector<double> tmpTarget(coordsTarget);
        if(SPACEDIM==3)
          PlanarIntersector<MyMeshType,MyMatrix>::projectionThis(&tmpSource[0],&tmpTarget[0],tmpSource.size()/SPACEDIM,nbNodesT);
        for(int nodeIdT=0;nodeIdT<nbNodesT;nodeIdT++)
          {
            if(PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg2D(&tmpTarget[0]+nodeIdT*SPACEDIM,&tmpSource[0],tmpSource.size()/SPACEDIM,PlanarIntersector<MyMeshType,MyMatrix>::_precision))
              {
                ConnType curNodeTInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConnT[nodeIdT]);
                typename MyMatrix::value_type& resRow=res[curNodeTInCmode];
                typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(*iter2));
                if(iterRes==resRow.end())
                  resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(*iter2),1.));
              }
          }
      }
  }

  template<class MyMeshType, class MyMatrix>
  int PlanarIntersectorP0P1PL<MyMeshType,MyMatrix>::getNumberOfRowsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getNumberOfNodes();
  }
  
  template<class MyMeshType, class MyMatrix>
  int PlanarIntersectorP0P1PL<MyMeshType,MyMatrix>::getNumberOfColsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getNumberOfElements();
  }
}

#endif
