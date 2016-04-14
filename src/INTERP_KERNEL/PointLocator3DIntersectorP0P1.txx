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
#ifndef __POINTLOCATOR3DINTERSECTORP0P1_TXX__
#define __POINTLOCATOR3DINTERSECTORP0P1_TXX__

#include "PointLocator3DIntersectorP0P1.hxx"
#include "Intersector3DP0P1.txx"
#include "MeshUtils.hxx"

namespace INTERP_KERNEL
{

  /** 
   * @param targetMesh  mesh containing the target elements
   * @param srcMesh     mesh containing the source elements
   * @param policy      splitting policy to be used
   */
  template<class MyMeshType, class MyMatrix>
  PointLocator3DIntersectorP0P1<MyMeshType,MyMatrix>::PointLocator3DIntersectorP0P1(const MyMeshType& targetMesh, const MyMeshType& srcMesh, double precision):
    Intersector3DP0P1<MyMeshType,MyMatrix>(targetMesh,srcMesh),_precision(precision)
  {
  }

  template<class MyMeshType, class MyMatrix>
  PointLocator3DIntersectorP0P1<MyMeshType,MyMatrix>::~PointLocator3DIntersectorP0P1()
  {
  }

  /**
   * 
   * @param targetCell in C mode.
   * @param srcCells in C mode.
   *
   */
  template<class MyMeshType, class MyMatrix>
  void PointLocator3DIntersectorP0P1<MyMeshType,MyMatrix>::intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res)
  {
    std::vector<double> coordsTarget;
    Intersector3DP0P1<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(targetCell),coordsTarget);
    int nbNodesT=coordsTarget.size()/SPACEDIM;
    const double *coordsS=Intersector3DP0P1<MyMeshType,MyMatrix>::_src_mesh.getCoordinatesPtr();
    const ConnType *startOfCellNodeConnT=Intersector3DP0P1<MyMeshType,MyMatrix>::getStartConnOfTargetCell(targetCell);
    for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
      {
        NormalizedCellType tS=Intersector3DP0P1<MyMeshType,MyMatrix>::_src_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(*iterCellS));
        const CellModel& cmTypeS=CellModel::GetCellModel(tS);
        std::vector<ConnType> connOfCurCellS;
        Intersector3DP0P1<MyMeshType,MyMatrix>::getConnOfSourceCell(OTT<ConnType,numPol>::indFC(*iterCellS),connOfCurCellS);
        for(int nodeIdT=0;nodeIdT<nbNodesT;nodeIdT++)
          {
            if(PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg3D(&coordsTarget[nodeIdT*SPACEDIM],&connOfCurCellS[0],connOfCurCellS.size(),coordsS,cmTypeS,_precision))
              {
                ConnType curNodeTInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConnT[nodeIdT]);
                typename MyMatrix::value_type& resRow=res[curNodeTInCmode];
                typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(*iterCellS));
                if(iterRes==resRow.end())
                  resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(*iterCellS),1.));
              }
          }
      }
  }
}

#endif
