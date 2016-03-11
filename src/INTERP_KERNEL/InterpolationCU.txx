// Copyright (C) 2009-2016  OPEN CASCADE
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
// File      : InterpolationCU.txx
// Created   : Mon Dec 14 17:30:25 2009
// Author    : Edward AGAPOV (eap)

#ifndef __InterpolationCU_TXX__
#define __InterpolationCU_TXX__

#include "InterpolationCU.hxx"

#include "Interpolation.txx"
#include "IntersectorCU1D.txx"
#include "IntersectorCU2D.txx"
#include "IntersectorCU3D.txx"

#include <map>

// // convert index "From Mesh Index"
#define _FMIU(i) OTT<typename MyUMeshType::MyConnType,MyUMeshType::My_numPol>::ind2C((i))
#define _FMIC(i) OTT<typename MyCMeshType::MyConnType,MyCMeshType::My_numPol>::ind2C((i))
// convert index "To Mesh Index"
#define _TMIU(i) OTT<typename MyUMeshType::MyConnType,MyUMeshType::My_numPol>::indFC((i))
#define _TMIC(i) OTT<typename MyCMeshType::MyConnType,MyCMeshType::My_numPol>::indFC((i))
// convert coord "From Mesh Coord"
#define _FMCOO(i) OTT<typename MyUMeshType::MyConnType,MyUMeshType::My_numPol>::coo2C((i))
// convert connectivity "From Mesh Connectivity"
#define _FMCON(i) OTT<typename MyUMeshType::MyConnType,MyUMeshType::My_numPol>::conn2C((i))

namespace INTERP_KERNEL
{
  /**
   * \defgroup InterpolationCU InterpolationCU
   * \class InterpolationCU
   * \brief Class used to calculate the volumes of intersection between the elements of a cartesian and an unstructured  meshes.
   * 
   */
  //================================================================================
  /**
   * Default constructor
   * 
   */
  //================================================================================

  InterpolationCU::InterpolationCU()
  {
  }

  InterpolationCU::InterpolationCU(const InterpolationOptions & io)
    :Interpolation<InterpolationCU>(io)
  {
  }

  //================================================================================
  /**
   * Calculates the matrix of volumes of intersection between the elements of srcMesh and the elements of targetMesh.
   * The calculation is done in two steps. First a filtering process reduces the number of pairs of elements for which the
   * calculation must be carried out by eliminating pairs that do not intersect based on their bounding boxes. Then, the 
   * volume of intersection is calculated for the remaining pairs, and entered into the
   * intersection matrix. 
   * 
   * The matrix is partially sparse : it is a vector of maps of integer - double pairs. 
   * It can also be an INTERP_KERNEL::Matrix object.
   * The length of the vector is equal to the number of target elements - for each target element there is a map, regardless
   * of whether the element intersects any source elements or not. But in the maps there are only entries for those source elements
   * which have a non-zero intersection volume with the target element. The vector has indices running from 
   * 0 to (nb target elements - 1), meaning that the map for target element i is stored at index i - 1. In the maps, however,
   * the indexing is more natural : the intersection volume of the target element i with source element j is found at matrix[i-1][j].
   * 

   * @param srcMesh     cartesian source mesh
   * @param targetMesh  unstructured target mesh
   * @param result      matrix in which the result is stored 
   * @param method      interpolation method
   */
  //================================================================================

  template<class MyCMeshType, class MyUMeshType, class MatrixType>
  int InterpolationCU::interpolateMeshes(const MyCMeshType& src_mesh,
                                         const MyUMeshType& tgt_mesh,
                                         MatrixType&        result,
                                         const char *       method)
  {
    typedef typename MyCMeshType::MyConnType CConnType;

    if ( std::string("P0P0") != method )
      throw Exception("Only P0P0 method is implemented so far");
    if ( MyCMeshType::MY_SPACEDIM != MyUMeshType::MY_SPACEDIM ||
         MyCMeshType::MY_SPACEDIM != MyUMeshType::MY_MESHDIM )
      throw Exception("InterpolationCU::interpolateMeshes(): dimension of meshes must be same");

    const double eps = getPrecision();
    const int dim = MyCMeshType::MY_SPACEDIM;

    TargetIntersector<MyCMeshType, MatrixType>* intersector = 0;
    switch( dim )
      {
      case 1: intersector = new IntersectorCU1D<MyCMeshType, MyUMeshType, MatrixType>( src_mesh, tgt_mesh ); break;
      case 2: intersector = new IntersectorCU2D<MyCMeshType, MyUMeshType, MatrixType>( src_mesh, tgt_mesh ); break;
      case 3: intersector = new IntersectorCU3D<MyCMeshType, MyUMeshType, MatrixType>( src_mesh, tgt_mesh, getSplittingPolicy() ); break;
      }
    // create empty maps for all target elements
    result.resize( intersector->getNumberOfRowsOfResMatrix() );
    const int ret = intersector->getNumberOfColsOfResMatrix();

    const double* src_coords[ dim ];
    int        src_nb_coords[ dim ];
    std::map< double, int> src_coord_to_index[ dim ];
    for ( int j = 0; j < dim; ++j )
      {
        src_coords   [j] = src_mesh.getCoordsAlongAxis( _TMIC( j ));
        src_nb_coords[j] = src_mesh.nbCellsAlongAxis  ( _TMIC( j )) + 1;
        for (int i = 0; i < src_nb_coords[j]; ++i )
          src_coord_to_index[j].insert( std::make_pair( src_coords[j][i], i ));
      }

    const unsigned long tgtu_nb_cells = tgt_mesh.getNumberOfElements();

    IntersectorCU<MyCMeshType, MyUMeshType, MatrixType> bbHelper(src_mesh, tgt_mesh);
    double bb[2*dim];

    // loop on unstructured tgt cells

    for(unsigned int iT=0; iT<tgtu_nb_cells; iT++)
      {
        result[ iT ].clear();

        // get bounding box of target cell
        bbHelper.getUElemBB( bb, _TMIU(iT));

        bool doItersect = true;
        for ( int j = 0; j < dim && doItersect; ++j )
          doItersect =
            bb[j*2]   < src_coords[j][ src_nb_coords[j]-1 ] - eps &&
            bb[j*2+1] > src_coords[j][0] + eps;
        if ( !doItersect )
          continue; // no intersection

        // find structured src cells intersecting iT cell
        std::vector< std::vector< CConnType > > structIndices(1);
        std::map< double, int>::iterator coo_ind;
        for ( int j = 0; j < dim; ++j )
          {
            coo_ind = src_coord_to_index[j].lower_bound( bb[2*j+1] - eps );
            if ( coo_ind == src_coord_to_index[j].end() )
              --coo_ind;
            int max_i = coo_ind->second;

            coo_ind = src_coord_to_index[j].upper_bound( bb[2*j  ] + eps );
            if ( coo_ind != src_coord_to_index[j].begin() )
              --coo_ind;
            int min_i = coo_ind->second;

            std::vector< std::vector< CConnType > > newStructIndices;
            for ( unsigned int iInd = 0; iInd < structIndices.size(); ++iInd )
              {
                for ( int i = min_i; i < max_i; ++i )
                  {
                    std::vector< CConnType > index = structIndices[iInd];
                    index.push_back( i );
                    newStructIndices.push_back( index );
                  }
              }
            structIndices.swap( newStructIndices );
          }

        // perform intersection

        for ( unsigned int iInd = 0; iInd < structIndices.size(); ++iInd )
          intersector->intersectCells( iT, structIndices[iInd], result );
      }
    delete intersector;
    return ret;
  }

  //================================================================================
  /**
   * Calculates the matrix of volumes of intersection between the elements of srcMesh and the elements of targetMesh.
   * The calculation is done in two steps. First a filtering process reduces the number of pairs of elements for which the
   * calculation must be carried out by eliminating pairs that do not intersect based on their bounding boxes. Then, the 
   * volume of intersection is calculated for the remaining pairs, and entered into the
   * intersection matrix. 
   * 
   * The matrix is partially sparse : it is a vector of maps of integer - double pairs. 
   * It can also be an INTERP_KERNEL::Matrix object.
   * The length of the vector is equal to the number of target elements - for each target element there is a map, regardless
   * of whether the element intersects any source elements or not. But in the maps there are only entries for those source elements
   * which have a non-zero intersection volume with the target element. The vector has indices running from 
   * 0 to (nb target elements - 1), meaning that the map for target element i is stored at index i - 1. In the maps, however,
   * the indexing is more natural : the intersection volume of the target element i with source element j is found at matrix[i-1][j].
   * 

   * @param srcMesh     2-dimesional unstructured target mesh
   * @param targetMesh  2-dimensional cartesian source mesh
   * @param result      matrix in which the result is stored 
   * @param method      interpolation method
   */
  //================================================================================

  template<class MyUMeshType, class MyCMeshType, class MatrixType>
  int InterpolationCU::interpolateMeshesRev(const MyUMeshType& meshS, const MyCMeshType& meshT, MatrixType& result, const char *method)
  {
    MatrixType revResult;
    int sizeT = interpolateMeshes( meshT, meshS, revResult, method );
    int sizeS = revResult.size();
    result.resize( sizeT );

    for ( int iS = 0; iS < sizeS; ++iS )
      {
        typename MatrixType::value_type & row = revResult[iS];
        typename MatrixType::value_type::iterator iT_surf = row.begin();
        for ( ; iT_surf != row.end(); ++iT_surf )
          result[ iT_surf->first ][ iS ] = iT_surf->second;
      }
    return sizeS;
  }

}

#endif
