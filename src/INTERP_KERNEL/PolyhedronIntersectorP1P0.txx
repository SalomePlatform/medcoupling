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
#ifndef __POLYHEDRONINTERSECTORP1P0_TXX__
#define __POLYHEDRONINTERSECTORP1P0_TXX__

#include "PolyhedronIntersectorP1P0.hxx"
#include "Intersector3DP1P0.txx"
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
  PolyhedronIntersectorP1P0<MyMeshType,MyMatrix>::PolyhedronIntersectorP1P0(const MyMeshType& targetMesh, const MyMeshType& srcMesh, SplittingPolicy policy):Intersector3DP1P0<MyMeshType,MyMatrix>(targetMesh,srcMesh),_split(srcMesh,targetMesh,policy)
  {
  }

  /**
   * Destructor.
   * Liberates the SplitterTetra objects and potential sub-node points that have been allocated.
   *
   */
  template<class MyMeshType, class MyMatrix>
  PolyhedronIntersectorP1P0<MyMeshType,MyMatrix>::~PolyhedronIntersectorP1P0()
  {
    releaseArrays();
  }
    
  template<class MyMeshType, class MyMatrix>
  void PolyhedronIntersectorP1P0<MyMeshType,MyMatrix>::releaseArrays()
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
   * WARNING : for all methods on _split object source and target are switched !
   */
  template<class MyMeshType, class MyMatrix>
  void PolyhedronIntersectorP1P0<MyMeshType,MyMatrix>::intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res)
  {
    SplitterTetra<MyMeshType>* subTetras[24];
    typename MyMatrix::value_type& resRow=res[targetCell];
    for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
      {
        releaseArrays();
        int nbOfNodesS=Intersector3D<MyMeshType,MyMatrix>::_src_mesh.getNumberOfNodesOfElement(OTT<ConnType,numPol>::indFC(*iterCellS));
        _split.splitTargetCell(*iterCellS,nbOfNodesS,_tetra);
        for(typename std::vector<SplitterTetra<MyMeshType>*>::iterator iter = _tetra.begin(); iter != _tetra.end(); ++iter)
          {
            (*iter)->splitIntoDualCells(subTetras);
            for(int i=0;i<24;i++)
              {
                SplitterTetra<MyMeshType> *tmp=subTetras[i];
                double volume = tmp->intersectSourceCell(targetCell);
                ConnType sourceNode=tmp->getId(0);
                if(volume!=0.)
                  {
                    typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(sourceNode));
                    if(iterRes==resRow.end())
                      resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(sourceNode),volume));
                    else
                      {
                        double val=(*iterRes).second+volume;
                        resRow.erase(OTT<ConnType,numPol>::indFC(sourceNode));
                        resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(sourceNode),val));
                      }
                  }
                delete tmp;
              }
          }
      }
  }
}

#endif
