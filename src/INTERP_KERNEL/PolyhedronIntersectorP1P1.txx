//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#ifndef __PolyhedronIntersectorP1P1_TXX__
#define __PolyhedronIntersectorP1P1_TXX__

#include "PolyhedronIntersectorP1P1.hxx"
#include "Intersector3DP1P1.txx"
#include "MeshUtils.hxx"

#include "SplitterTetra.txx"

namespace INTERP_KERNEL
{

  /**
   * Constructor creating object from target cell global number 
   * The constructor first calculates the necessary nodes, 
   * (depending on the splitting policy) and then splits the hexahedron into 
   * tetrahedra, placing these in the internal vector _tetra.
   * 
   * @param targetMesh  mesh containing the target elements
   * @param srcMesh     mesh containing the source elements
   * @param policy      splitting policy to be used
   *
   * WARNING : in _split attribute, sourceMesh and targetMesh are switched in order to fit intersectCells feature.
   */
  template<class MyMeshType, class MyMatrix>
  PolyhedronIntersectorP1P1<MyMeshType,MyMatrix>::PolyhedronIntersectorP1P1(const MyMeshType& targetMesh, const MyMeshType& srcMesh, SplittingPolicy policy):Intersector3DP1P1<MyMeshType,MyMatrix>(targetMesh,srcMesh)
  {
  }

  /**
   * Destructor.
   * Liberates the SplitterTetra objects and potential sub-node points that have been allocated.
   *
   */
  template<class MyMeshType, class MyMatrix>
  PolyhedronIntersectorP1P1<MyMeshType,MyMatrix>::~PolyhedronIntersectorP1P1()
  {
  }

  /**
   * Calculates the volume of intersection of an element in the source mesh and the target element
   * represented by the object.
   * The calculation is performed by calling the corresponding method for
   * each SplitterTetra object created by the splitting.
   * 
   * @param targetCell in C mode.
   * @param srcCells in C mode.
   *
   * WARNING : for all methods on _split object source and target are switched !
   */
  template<class MyMeshType, class MyMatrix>
  void PolyhedronIntersectorP1P1<MyMeshType,MyMatrix>::intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res)
  {
    // split the targetCell into dual cells
    std::pair< int, double[12] > subTetraNodes[24]; // a node of sub tetra and its coordinates
    const double* nodes[4]; int conn[4];
    for(int node = 0; node < 4 ; ++node)
      nodes[node]=getCoordsOfNode2(node, OTT<ConnType,numPol>::indFC(targetCell),
                                   Intersector3D<MyMeshType,MyMatrix>::_target_mesh,conn[node]);
    SplitterTetra<MyMeshType> tgtTetra(Intersector3D<MyMeshType,MyMatrix>::_src_mesh, nodes, conn);
    for (int i=0; i<24; i++)
      tgtTetra.splitMySelfForDual(subTetraNodes[i].second,i,subTetraNodes[i].first);

    // intersect each source tetrahedron with each of target dual cells
    SplitterTetra<MyMeshType>* subTetrasS[24];
    for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
    {
      // split a source cell into dual cells
      for(int node = 0; node < 4 ; ++node)
        nodes[node]=getCoordsOfNode2(node, OTT<ConnType,numPol>::indFC(*iterCellS),
                                     Intersector3D<MyMeshType,MyMatrix>::_src_mesh,conn[node]);

      SplitterTetra<MyMeshType> srcTetra(Intersector3D<MyMeshType,MyMatrix>::_target_mesh, nodes, conn);
      srcTetra.splitIntoDualCells(subTetrasS);

      // intersect each target subTetra with each source one
      for(int i=0;i<24;i++)
      {
        SplitterTetra<MyMeshType> *tmp=subTetrasS[i];
        ConnType sourceNode=OTT<ConnType,numPol>::indFC(tmp->getId(0));
        for(int j=0;j<24;j++)
        {
          const double* tetraNodes12 = subTetraNodes[j].second;
          const double* tetraNodesT[4]={ tetraNodes12, tetraNodes12+3, tetraNodes12+6, tetraNodes12+9 };
          double volume = tmp->intersectTetra( tetraNodesT );
          if(volume!=0.)
          {
            ConnType tgtNode=OTT<ConnType,numPol>::indFC( subTetraNodes[j].first );
            typename MyMatrix::value_type& resRow = res[tgtNode];
            typename MyMatrix::value_type::const_iterator iterRes=resRow.find( sourceNode );
            if(iterRes!=resRow.end())
            {
              volume += (*iterRes).second;
              resRow.erase(sourceNode);
            }
            resRow.insert(std::make_pair(sourceNode,volume));
          }
        }
        delete tmp;
      }
    }
  }
}

#endif
