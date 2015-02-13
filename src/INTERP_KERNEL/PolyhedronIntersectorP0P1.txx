// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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
#ifndef __POLYHEDRONINTERSECTORP0P1_TXX__
#define __POLYHEDRONINTERSECTORP0P1_TXX__

#include "PolyhedronIntersectorP0P1.hxx"
#include "Intersector3DP0P1.txx"
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
   */
  template<class MyMeshType, class MyMatrix>
  PolyhedronIntersectorP0P1<MyMeshType,MyMatrix>::PolyhedronIntersectorP0P1(const MyMeshType& targetMesh, const MyMeshType& srcMesh, SplittingPolicy policy):Intersector3DP0P1<MyMeshType,MyMatrix>(targetMesh,srcMesh),_split(targetMesh,srcMesh,policy)
  {
  }

  /**
   * Destructor.
   * Liberates the SplitterTetra objects and potential sub-node points that have been allocated.
   *
   */
  template<class MyMeshType, class MyMatrix>
  PolyhedronIntersectorP0P1<MyMeshType,MyMatrix>::~PolyhedronIntersectorP0P1()
  {
    releaseArrays();
  }
    
  template<class MyMeshType, class MyMatrix>
  void PolyhedronIntersectorP0P1<MyMeshType,MyMatrix>::releaseArrays()
  {
    for(typename std::vector< SplitterTetra<MyMeshType>* >::iterator iter = _tetra.begin(); iter != _tetra.end(); ++iter)
      delete *iter;
    _split.releaseArrays();
    _tetra.clear();
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
   */
  template<class MyMeshType, class MyMatrix>
  void PolyhedronIntersectorP0P1<MyMeshType,MyMatrix>::intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res)
  {
    SplitterTetra<MyMeshType>* subTetras[24];
    releaseArrays();
    _split.splitTargetCell2(targetCell,_tetra);
    for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
      {
        for(typename std::vector<SplitterTetra<MyMeshType>*>::iterator iter = _tetra.begin(); iter != _tetra.end(); ++iter)
          {
            (*iter)->splitIntoDualCells(subTetras);
            for(int i=0;i<24;i++)
              {
                SplitterTetra<MyMeshType> *tmp=subTetras[i];
                double volume = tmp->intersectSourceCell(*iterCellS);
                if(volume!=0.)
                  {
                    int targetNodeId(tmp->getId(0));
                    if(targetNodeId<0)
                      {
                        std::ostringstream oss; oss << "PolyhedronIntersectorP0P1::intersectCells : On target cell #" <<  targetCell << " the splitting into tetra4 leads to the creation of an additional point that interacts with source cell Id #" << *iterCellS << " !";
                        throw INTERP_KERNEL::Exception(oss.str().c_str());
                      }
                    typename MyMatrix::value_type& resRow=res[targetNodeId];
                    typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(*iterCellS));
                    if(iterRes==resRow.end())
                      resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(*iterCellS),volume));
                    else
                      {
                        double val=(*iterRes).second+volume;
                        resRow.erase(OTT<ConnType,numPol>::indFC(*iterCellS));
                        resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(*iterCellS),val));
                      }
                  }
                delete tmp;
              }
          }
      }
  }
}

#endif
