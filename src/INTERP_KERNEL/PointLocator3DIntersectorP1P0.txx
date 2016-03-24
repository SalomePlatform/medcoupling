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
#ifndef __POINTLOCATOR3DINTERSECTORP1P0_TXX__
#define __POINTLOCATOR3DINTERSECTORP1P0_TXX__

#include "PointLocator3DIntersectorP1P0.hxx"
#include "Intersector3DP1P0.txx"
#include "MeshUtils.hxx"

namespace INTERP_KERNEL
{
  /**
   * @param targetMesh  mesh containing the target elements
   * @param srcMesh     mesh containing the source elements
   * @param policy      splitting policy to be used
   *
   * WARNING : in _split attribute, sourceMesh and targetMesh are switched in order to fit intersectCells feature.
   */
  template<class MyMeshType, class MyMatrix>
  PointLocator3DIntersectorP1P0<MyMeshType,MyMatrix>::PointLocator3DIntersectorP1P0(const MyMeshType& targetMesh, const MyMeshType& srcMesh, double precision):Intersector3DP1P0<MyMeshType,MyMatrix>(targetMesh,srcMesh),_precision(precision)
  {
  }

  template<class MyMeshType, class MyMatrix>
  PointLocator3DIntersectorP1P0<MyMeshType,MyMatrix>::~PointLocator3DIntersectorP1P0()
  {
  }

  /**
   * @param targetCell in C mode.
   * @param srcCells in C mode.
   *
   * WARNING : for all methods on _split object source and target are switched !
   */
  template<class MyMeshType, class MyMatrix>
  void PointLocator3DIntersectorP1P0<MyMeshType,MyMatrix>::intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res)
  {
    std::vector<double> CoordsT;
    typename MyMatrix::value_type& resRow=res[targetCell];
    const double *coordsS=Intersector3DP1P0<MyMeshType,MyMatrix>::_src_mesh.getCoordinatesPtr();
    Intersector3DP1P0<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(targetCell),CoordsT);
    double baryT[SPACEDIM];
    calculateBarycenterDyn2<SPACEDIM>(&CoordsT[0],CoordsT.size()/SPACEDIM,baryT);
    for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
      {
        NormalizedCellType tS=Intersector3DP1P0<MyMeshType,MyMatrix>::_src_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(*iterCellS));
        if(tS!=NORM_TETRA4)
          throw INTERP_KERNEL::Exception("Invalid source cell detected for meshdim==3. Only TETRA4 supported !");
        const CellModel& cmTypeS=CellModel::GetCellModel(tS);
        std::vector<ConnType> connOfCurCellS;
        Intersector3DP1P0<MyMeshType,MyMatrix>::getConnOfSourceCell(OTT<ConnType,numPol>::indFC(*iterCellS),connOfCurCellS);
        if( PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg3D(baryT,&connOfCurCellS[0],connOfCurCellS.size(),coordsS,cmTypeS,_precision) )
          {
            double resLoc[4];
            std::vector<double> srcCell;
            Intersector3DP1P0<MyMeshType,MyMatrix>::getRealSourceCoordinates(OTT<ConnType,numPol>::indFC(*iterCellS),srcCell);
            std::vector<const double*> eap(4);
            eap[0]=&srcCell[0]; eap[1]=&srcCell[3]; eap[2]=&srcCell[6]; eap[3]=&srcCell[9];
            barycentric_coords(eap,baryT,resLoc);
            const ConnType *startOfCellNodeConn=Intersector3DP1P0<MyMeshType,MyMatrix>::getStartConnOfSourceCell(*iterCellS);
            for(int nodeIdS=0;nodeIdS<4;nodeIdS++)
              {
                if(fabs(resLoc[nodeIdS])>_precision)
                  {
                    ConnType curNodeSInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[nodeIdS]);
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

#endif
