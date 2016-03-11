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
#ifndef __PolyhedronIntersectorP1P0Bary_TXX__
#define __PolyhedronIntersectorP1P0Bary_TXX__

#include "PolyhedronIntersectorP1P0Bary.hxx"
#include "Intersector3DP1P0Bary.txx"
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
  PolyhedronIntersectorP1P0Bary<MyMeshType,MyMatrix>::PolyhedronIntersectorP1P0Bary(const MyMeshType& targetMesh,
                                                                                    const MyMeshType& srcMesh,
                                                                                    SplittingPolicy policy)
    :Intersector3DP1P0Bary<MyMeshType,MyMatrix>(targetMesh,srcMesh),_split(targetMesh,srcMesh,policy)
  {
    // SPEC:
    // "Limitation. For the P1P0 barycentric improvement only triangle source cells in 2D and
    // tetrahedrons in 3D will be supported by interpolators. If a non
    // triangle/tetrahedron source cell is detected an INTERP_KERNEL::Exception should be thrown."

    // Check types of source elements here rather than in intersectCells() since a wrong type can be
    // found late after a long time of calculation.

    const unsigned long numSrcElems = srcMesh.getNumberOfElements();
    for(unsigned long i = 0 ; i < numSrcElems ; ++i)
      if ( srcMesh.getTypeOfElement( OTT<ConnType,numPol>::indFC(i) ) != NORM_TETRA4 )
        throw INTERP_KERNEL::Exception("P1P0 barycentric algorithm works only with tetrahedral source meshes");
  }

  /**
   * Destructor.
   * Liberates the SplitterTetra objects and potential sub-node points that have been allocated.
   *
   */
  template<class MyMeshType, class MyMatrix>
  PolyhedronIntersectorP1P0Bary<MyMeshType,MyMatrix>::~PolyhedronIntersectorP1P0Bary()
  {
    releaseArrays();
  }

  template<class MyMeshType, class MyMatrix>
  void PolyhedronIntersectorP1P0Bary<MyMeshType,MyMatrix>::releaseArrays()
  {
    for(typename std::vector< SplitterTetra<MyMeshType>* >::iterator iter = _tetra.begin(); iter != _tetra.end(); ++iter)
      delete *iter;
    _split.releaseArrays();
    _tetra.clear();
  }

  //================================================================================
  /*!
   * \brief This method computes a value per each node of source cell for each target cell.
   *  \param srcCell - a source tetrahedron
   *  \param tgtCells - target elements
   *  \param res - matrix to fill in 
   */
  //================================================================================

  template<class MyMeshType, class MyMatrix>
  void PolyhedronIntersectorP1P0Bary<MyMeshType,MyMatrix>::intersectCells(ConnType                     tgtCell,
                                                                          const std::vector<ConnType>& srcCells,
                                                                          MyMatrix&                    res)
  {
    typename MyMatrix::value_type& resRow=res[tgtCell];

    int nbOfNodesT=Intersector3D<MyMeshType,MyMatrix>::_target_mesh.getNumberOfNodesOfElement(OTT<ConnType,numPol>::indFC(tgtCell));
    releaseArrays();
    _split.splitTargetCell(tgtCell,nbOfNodesT,_tetra);

    for(typename std::vector<ConnType>::const_iterator iterCellS=srcCells.begin();iterCellS!=srcCells.end();iterCellS++)
    {
      // intersect a source tetrahedron with each target tetrahedron: get intersection volume and barycenter
      double baryCentre[SPACEDIM], total_baryCentre[3] = { 0., 0., 0.};
      double interVolume = 0;
      for(typename std::vector<SplitterTetra<MyMeshType>*>::iterator iterTetraT = _tetra.begin(); iterTetraT != _tetra.end(); ++iterTetraT)
      {
        SplitterTetra<MyMeshType> *tmp=*iterTetraT;
        tmp->clearVolumesCache();
        double volume = tmp->intersectSourceCell(*iterCellS, baryCentre);
        if ( volume > 0 )
        {
          interVolume += volume;
          for ( int i = 0; i < SPACEDIM; ++i )
            total_baryCentre[i] += baryCentre[i]*volume;
        }
      }
      if(interVolume!=0)
      {
        for ( int i = 0; i < SPACEDIM; ++i )
          total_baryCentre[i] /= interVolume;

        // coordinates of the source tetrahedron
        std::vector<const double*> srcCellCoords(4);
        for ( int n = 0; n < 4; ++n )
          srcCellCoords[ n ] = getCoordsOfNode( n, *iterCellS, Intersector3D<MyMeshType,MyMatrix>::_src_mesh );

        // compute barycentric coordinates
        double baryCoords[4];
        barycentric_coords( srcCellCoords, total_baryCentre, baryCoords);

        // store coeffs of each node of the source tetrahedron
        const ConnType *srcCellNodes=Intersector3D<MyMeshType,MyMatrix>::_src_mesh.getConnectivityPtr()+OTT<ConnType,numPol>::conn2C(Intersector3D<MyMeshType,MyMatrix>::_src_mesh.getConnectivityIndexPtr()[*iterCellS]);
        for ( int n = 0; n < 4; ++n )
        {
          double val = baryCoords[n] * interVolume;
          ConnType curNodeS = srcCellNodes[n];
          typename MyMatrix::value_type::const_iterator iterRes=resRow.find(curNodeS);
          if(iterRes!=resRow.end())
          {
            val += iterRes->second;
            resRow.erase( curNodeS );
          }
          resRow.insert(std::make_pair(curNodeS,val));
        }
      }
    }
  }
}

#endif
