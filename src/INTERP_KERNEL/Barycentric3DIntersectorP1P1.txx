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

#ifndef __BARYCENTRIC3DINTERSECTORP1P1_TXX__
#define __BARYCENTRIC3DINTERSECTORP1P1_TXX__

#include "Barycentric3DIntersectorP1P1.hxx"
#include "Intersector3DP1P1.txx"
#include "MeshUtils.hxx"

namespace INTERP_KERNEL
{

  /**
   * Constructor creating object from target cell global number 
   * 
   * @param targetMesh  mesh containing the target elements
   * @param srcMesh     mesh containing the source elements
   * @param policy      splitting policy to be used
   */
  template<class MyMeshType, class MyMatrix>
  Barycentric3DIntersectorP1P1<MyMeshType,MyMatrix>::Barycentric3DIntersectorP1P1(const MyMeshType& targetMesh, const MyMeshType& srcMesh, double precision):
    Intersector3DP1P1<MyMeshType,MyMatrix>(targetMesh,srcMesh),_precision(precision)
  {
  }

  template<class MyMeshType, class MyMatrix>
  Barycentric3DIntersectorP1P1<MyMeshType,MyMatrix>::~Barycentric3DIntersectorP1P1()
  {
  }

  /**
   * @param targetCell in C mode.
   * @param srcCells in C mode.
   */
  template<class MyMeshType, class MyMatrix>
  void Barycentric3DIntersectorP1P1<MyMeshType,MyMatrix>::intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res)
  {
    std::vector<double> CoordsT;
    const ConnType *startOfCellNodeConnT=Intersector3DP1P1<MyMeshType,MyMatrix>::getStartConnOfTargetCell(targetCell);
    Intersector3DP1P1<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(targetCell),CoordsT);
    int nbOfNodesT=CoordsT.size()/SPACEDIM;
    const double *coordsS=Intersector3DP1P1<MyMeshType,MyMatrix>::_src_mesh.getCoordinatesPtr();
    for(int nodeIdT=0;nodeIdT<nbOfNodesT;nodeIdT++)
      {
        typename MyMatrix::value_type& resRow=res[OTT<ConnType,numPol>::ind2C(startOfCellNodeConnT[nodeIdT])];
        if(!resRow.empty())
          continue;
        for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
          {
            NormalizedCellType tS=Intersector3DP1P1<MyMeshType,MyMatrix>::_src_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(*iterCellS));
            if(tS!=NORM_TETRA4)
              throw INTERP_KERNEL::Exception("Invalid source cell detected for meshdim==3. Only TETRA4 supported !");
            const CellModel& cmTypeS=CellModel::GetCellModel(tS);
            //
            std::vector<ConnType> connOfCurCellS;
            Intersector3DP1P1<MyMeshType,MyMatrix>::getConnOfSourceCell(OTT<ConnType,numPol>::indFC(*iterCellS),connOfCurCellS);
            if( PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg3D(&CoordsT[nodeIdT*SPACEDIM],&connOfCurCellS[0],connOfCurCellS.size(),coordsS,cmTypeS,_precision) )
              {
                double resLoc[4];
                std::vector<double> localCoordsS;
                Intersector3DP1P1<MyMeshType,MyMatrix>::getRealSourceCoordinates(OTT<ConnType,numPol>::indFC(*iterCellS),localCoordsS);
                std::vector<const double*> eap(4);
                eap[0]=&localCoordsS[0]; eap[1]=&localCoordsS[3]; eap[2]=&localCoordsS[6]; eap[3]=&localCoordsS[9];
                barycentric_coords(eap,&CoordsT[nodeIdT*SPACEDIM],resLoc);
                const ConnType *startOfCellNodeConnS=Intersector3DP1P1<MyMeshType,MyMatrix>::getStartConnOfSourceCell(*iterCellS);
                for(int nodeIdS=0;nodeIdS<4;nodeIdS++)
                  {
                    if(fabs(resLoc[nodeIdS])>_precision)
                      {
                        ConnType curNodeSInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConnS[nodeIdS]);
                        resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(curNodeSInCmode),resLoc[nodeIdS]));
                      }
                  }
              }
          }
      }
  }
}

#endif
