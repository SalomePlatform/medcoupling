// Copyright (C) 2007-2023  CEA, EDF
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

  template<class RowType, class ConnType>
  void AddContributionInRow(RowType& row, ConnType colId, double value)
  {
    if(value != 0.)
    {
      typename RowType::const_iterator iterRes=row.find(colId);
      if(iterRes==row.end())
        row.insert(std::make_pair(colId,value));
      else
      {
        double val=(*iterRes).second+value;
        row.erase(colId);
        row.insert(std::make_pair(colId,val));
      }
    }
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
    typename MyMatrix::value_type& resRow=res[targetCell];
    INTERP_KERNEL::SplittingPolicy sp( _split.getSplittingPolicy() );
    if( sp == GENERAL_48 )
      THROW_IK_EXCEPTION("GENERAL_28 spliting is not supported for P1P0 interpolation");
    SplitterTetra<MyMeshType>* subTetras[24];
    for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
    {
      releaseArrays();
      ConnType nbOfNodesS=this->_src_mesh.getNumberOfNodesOfElement(OTT<ConnType,numPol>::indFC(*iterCellS));
      _split.splitTargetCell(*iterCellS,nbOfNodesS,_tetra);
      INTERP_KERNEL::NormalizedCellType srcType = this->_src_mesh.getTypeOfElement( OTT<ConnType,numPol>::indFC(*iterCellS) );
      if( srcType == NORM_TETRA4 || (srcType == NORM_HEXA8 && sp != GENERAL_24 ))
      {
        for(typename std::vector<SplitterTetra<MyMeshType>*>::const_iterator iter = _tetra.cbegin(); iter != _tetra.cend(); ++iter)
          {
            (*iter)->splitIntoDualCells(subTetras);
            double vol2 = 0.;
            for(int i=0;i<24;i++)
              {
                SplitterTetra<MyMeshType> *tmp=subTetras[i];
                double volume = tmp->intersectSourceCell(targetCell);
                vol2 += volume;
                ConnType sourceNode=tmp->getId(0);
                AddContributionInRow(resRow,OTT<ConnType,numPol>::indFC(sourceNode),volume);
                delete tmp;
              }
          }
      }
      else
      {// for HEXA and GENERAL_24 no need to use subsplitting into dual mesh
        for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
          {
            releaseArrays();
            ConnType nbOfNodesS=Intersector3D<MyMeshType,MyMatrix>::_src_mesh.getNumberOfNodesOfElement(OTT<ConnType,numPol>::indFC(*iterCellS));
            _split.splitTargetCell2(*iterCellS,_tetra);
            for(typename std::vector<SplitterTetra<MyMeshType>*>::const_iterator iter = _tetra.cbegin(); iter != _tetra.cend(); ++iter)
              {
                double volume = std::abs( (*iter)->intersectSourceCell(targetCell) );
                // node #0 is for internal node node #1 is for the node at the middle of the face
                ConnType sourceNode0( (*iter)->getId(0) ), sourceNode1( (*iter)->getId(1) );
                AddContributionInRow(resRow,OTT<ConnType,numPol>::indFC(sourceNode0),volume/2.);
                AddContributionInRow(resRow,OTT<ConnType,numPol>::indFC(sourceNode1),volume/2.);
              }
          }
      }
    }
  }
}

#endif
