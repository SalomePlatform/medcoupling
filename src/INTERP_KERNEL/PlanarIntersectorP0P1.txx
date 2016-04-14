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
#ifndef __PLANARINTERSECTORP0P1_TXX__
#define __PLANARINTERSECTORP0P1_TXX__

#include "PlanarIntersectorP0P1.hxx"
#include "InterpolationUtils.hxx"
#include "CellModel.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, class ConcreteP0P1Intersector>
  PlanarIntersectorP0P1<MyMeshType,MyMatrix,ConcreteP0P1Intersector>::PlanarIntersectorP0P1(const MyMeshType& meshT, const MyMeshType& meshS,
                                                                                            double dimCaracteristic, double precision, double md3DSurf, double minDot3DSurf, double medianPlane,
                                                                                            bool doRotate, int orientation, int printLevel):
    PlanarIntersector<MyMeshType,MyMatrix>(meshT,meshS,dimCaracteristic,precision,md3DSurf,minDot3DSurf,medianPlane,doRotate,orientation,printLevel)
  {
  }

  template<class MyMeshType, class MyMatrix, class ConcreteP0P1Intersector>
  int PlanarIntersectorP0P1<MyMeshType,MyMatrix,ConcreteP0P1Intersector>::getNumberOfRowsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getNumberOfNodes();
  }

  template<class MyMeshType, class MyMatrix, class ConcreteP0P1Intersector>
  int PlanarIntersectorP0P1<MyMeshType,MyMatrix,ConcreteP0P1Intersector>::getNumberOfColsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getNumberOfElements();
  }

  /*!
   * This methods split on the fly, into triangles in order to compute dual mesh of target cell (with icellT id in target mesh in C mode).
   */
  template<class MyMeshType, class MyMatrix, class ConcreteP0P1Intersector>
  void PlanarIntersectorP0P1<MyMeshType,MyMatrix,ConcreteP0P1Intersector>::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    int nbNodesT=PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[icellT+1]-PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[icellT];
    double triangle[9];
    double quadrangle[12];
    std::vector<double> sourceCellCoords;
    int orientation=1;
    const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectT+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[icellT]);
    for(int nodeIdT=0;nodeIdT<nbNodesT;nodeIdT++)
      {
        ConnType curNodeTInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[nodeIdT]);
        std::copy(PlanarIntersector<MyMeshType,MyMatrix>::_coordsT+curNodeTInCmode*SPACEDIM,
                  PlanarIntersector<MyMeshType,MyMatrix>::_coordsT+curNodeTInCmode*SPACEDIM+SPACEDIM,triangle);
        typename MyMatrix::value_type& resRow=res[curNodeTInCmode];
        for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++)
          {
            int iS=*iter;
            PlanarIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinates(OTT<ConnType,numPol>::indFC(iS),sourceCellCoords);
            for(int subTriT=1;subTriT<=nbNodesT-2;subTriT++)
              {
                std::copy(PlanarIntersector<MyMeshType,MyMatrix>::_coordsT+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[(nodeIdT+subTriT)%nbNodesT])*SPACEDIM,
                          PlanarIntersector<MyMeshType,MyMatrix>::_coordsT+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[(nodeIdT+subTriT)%nbNodesT])*SPACEDIM+SPACEDIM,
                          triangle+SPACEDIM);
                std::copy(PlanarIntersector<MyMeshType,MyMatrix>::_coordsT+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[(nodeIdT+subTriT+1)%nbNodesT])*SPACEDIM,
                          PlanarIntersector<MyMeshType,MyMatrix>::_coordsT+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[(nodeIdT+subTriT+1)%nbNodesT])*SPACEDIM+SPACEDIM,
                          triangle+2*SPACEDIM);
                fillDualCellOfTri<SPACEDIM>(triangle,quadrangle);
                std::vector<double> sourceCellCoordsTmp(sourceCellCoords);
                if(SPACEDIM==3)
                  orientation=PlanarIntersector<MyMeshType,MyMatrix>::projectionThis(&sourceCellCoordsTmp[0],quadrangle,sourceCellCoords.size()/SPACEDIM,4);
                NormalizedCellType tS=PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getTypeOfElement(OTT<ConnType,numPol>::indFC(iS));
                double surf=orientation*intersectGeometryWithQuadrangle(quadrangle,sourceCellCoordsTmp,CellModel::GetCellModel(tS).isQuadratic());
                surf=PlanarIntersector<MyMeshType,MyMatrix>::getValueRegardingOption(surf);
                if(surf!=0.)
                  {
                    typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(iS));
                    if(iterRes==resRow.end())
                      resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(iS),surf));
                    else
                      {
                        double val=(*iterRes).second+surf;
                        resRow.erase(OTT<ConnType,numPol>::indFC(iS));
                        resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(iS),val));
                      }
                  }
              }
          }
      }
  }
}

#endif
