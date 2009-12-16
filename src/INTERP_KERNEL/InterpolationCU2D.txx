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
// File      : InterpolationCU2D.txx
// Created   : Mon Dec 14 17:30:25 2009
// Author    : Edward AGAPOV (eap)

#ifndef __InterpolationCU2D_TXX__
#define __InterpolationCU2D_TXX__

#include "InterpolationCU2D.hxx"

#include "Interpolation.txx"
#include "PlanarIntersector.txx"
#include "TriangulationIntersector.txx"

#include <map>

// convert index "From Mesh Index"
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
   * \defgroup InterpolationCU2D InterpolationCU2D
   * \class InterpolationCU2D
   * \brief Class used to calculate the volumes of intersection between the elements of a cartesian and an unstructured 2D meshes.
   * 
   */
  //================================================================================
  /**
   * Default constructor
   * 
   */
  //================================================================================

  InterpolationCU2D::InterpolationCU2D()
  {
  }

  InterpolationCU2D::InterpolationCU2D(const InterpolationOptions & io)
    :Interpolation<InterpolationCU2D>(io)
  {
  }

  //================================================================================
  /*!
   * \brief A stub to anable inheriting from Interpolation
   */
  //================================================================================

  template<class MyMeshType, class MatrixType>
  int interpolateMeshes(const MyMeshType& meshS, const MyMeshType& meshT, MatrixType& result, const char *method)
  {
    throw INTERP_KERNEL::Exception("InterpolationCU2D class is intended to interpolate meshes of different nature: cartesian and unstructured");
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

   * @param srcMesh     2-dimensional cartesian source mesh
   * @param targetMesh  2-dimesional unstructured target mesh
   * @param result      matrix in which the result is stored 
   * @param method      interpolation method
   */
  //================================================================================

  template<class MyCMeshType, class MyUMeshType, class MatrixType>
  int InterpolationCU2D::interpolateMeshes(const MyCMeshType& src_mesh,
                                           const MyUMeshType& tgt_mesh,
                                           MatrixType&        result,
                                           const char *       method)
  {
    typedef typename MyUMeshType::MyConnType ConnType;

    if ( std::string("P0P0") != method )
      throw Exception("Only P0P0 method is implemented so far");

    // create empty maps for all target elements
    result.resize( tgt_mesh.getNumberOfElements() );

    const int ret = src_mesh.getNumberOfElements();

    const double eps = getPrecision();
    const int dim = 2;

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
    const ConnType *        tgtu_conn = tgt_mesh.getConnectivityPtr();
    const ConnType *    tgtu_conn_ind = tgt_mesh.getConnectivityIndexPtr();
    const double *         tgtu_coord = tgt_mesh.getCoordinatesPtr();

    double bb[2*dim];
    TriangulationIntersector<MyUMeshType,MatrixType,PlanarIntersectorP0P0>
      intersector (tgt_mesh,tgt_mesh, 0,0,0,0,0,0 );

    // loop on unstructured tgt cells

    for(int iT=0; iT<tgtu_nb_cells; iT++)
      {
        result[ iT ].clear();

        // get bounding box of target cell
        ConnType tgtu_nb_nodes = tgtu_conn_ind[iT+1] - tgtu_conn_ind[iT];
        intersector.getElemBB( bb, tgt_mesh, _TMIU(iT), tgtu_nb_nodes);

        if ( bb[0] > src_coords[0][ src_nb_coords[0]-1 ] - eps ||
             bb[2] > src_coords[1][ src_nb_coords[1]-1 ] - eps ||
             bb[1] < src_coords[0][0] + eps ||
             bb[3] < src_coords[1][0] + eps )
          continue; // no intersection

        // find structured src cells intersecting iT cell
        int min_i[dim], max_i[dim];
        std::map< double, int>::iterator coo_ind;
        for ( int j = 0; j < dim; ++j )
          {
            coo_ind = src_coord_to_index[j].lower_bound( bb[2*j+1] - eps );
            if ( coo_ind == src_coord_to_index[j].end() )
              --coo_ind;
            max_i[j] = coo_ind->second;

            coo_ind = src_coord_to_index[j].upper_bound( bb[2*j  ] + eps );
            if ( coo_ind != src_coord_to_index[j].begin() )
              --coo_ind;
            min_i[j] = coo_ind->second;
          }

        // get tgt cell coords
        std::vector<double> tgtu_coords(2*tgtu_nb_nodes);
        for (ConnType i=0; i<tgtu_nb_nodes; i++)
          for(int j=0; j<2; j++)
            tgtu_coords[2*i+j]=tgtu_coord[2*_FMCOO( tgtu_conn[_FMCON (tgtu_conn_ind[_FMIU(iT)]+i)])+j];
        
        NormalizedCellType tT = tgt_mesh.getTypeOfElement( _TMIU(iT));
        bool is_tgt_quad = CellModel::getCellModel(tT).isQuadratic();

        // perform intersection

        for ( int iX = min_i[0]; iX < max_i[0]; ++iX )
          for ( int iY = min_i[1]; iY < max_i[1]; ++iY )
            {
              double src_quad[8] = { src_coords[0][iX],   src_coords[1][iY],
                                     src_coords[0][iX+1], src_coords[1][iY], 
                                     src_coords[0][iX+1], src_coords[1][iY+1], 
                                     src_coords[0][iX],   src_coords[1][iY+1] };

              double surf = intersector.intersectGeometryWithQuadrangle( src_quad,
                                                                         tgtu_coords,
                                                                         is_tgt_quad);
              //surf = intersector.getValueRegardingOption(surf);
              if ( surf > eps )
                result[ iT ][ _TMIC( iX + iY * ( src_nb_coords[1]-1 )) ] = surf;
            }
      }

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
  int InterpolationCU2D::interpolateMeshesRev(const MyUMeshType& meshS, const MyCMeshType& meshT, MatrixType& result, const char *method)
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
