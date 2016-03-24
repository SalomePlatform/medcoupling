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
// File      : InterpolationCC.txx
// Created   : Fri Aug 14 11:39:27 2009
// Author    : Edward AGAPOV (eap)
//

#include "InterpolationCC.hxx"
#include "InterpolationUtils.hxx"

// convert index "From Mesh Index"
#define _FMI(i) OTT<typename MyMeshType::MyConnType,MyMeshType::My_numPol>::ind2C((i))
// convert index "To Mesh Index"
#define _TMI(i) OTT<typename MyMeshType::MyConnType,MyMeshType::My_numPol>::indFC((i))

namespace INTERP_KERNEL
{
  //================================================================================
  /*!
   * \brief Constructor does nothing
   */
  //================================================================================
  InterpolationCC::InterpolationCC()
  {
  }

  InterpolationCC::InterpolationCC(const InterpolationOptions& io):Interpolation<InterpolationCC>(io)
  {
  }

  //================================================================================
  /*!
   * \brief An 1D intersection result
   */
  //================================================================================

  struct Interference
  {
    int _src_index; // source cell index along an axis
    int _tgt_index; // target cell index along an axis
    double _length; // interference length
    Interference(int is = -1, int it = -1, double l = 0):_src_index(is),_tgt_index(it),_length(l){}
  };

  //================================================================================
  /*!
   * \brief Fills the matrix by precomputed cell interferences along axes
   *  \param inter_of_axis - cell/cell interferences along each axis
   *  \param result - matrix to fill in
   *  \param src_nb_cells[] - nb of cells along each of axes in the source mesh
   *  \param tgt_nb_cells[] - nb of cells along each of axes in the target mesh
   *  \param src_i_cell - source cell number accumulated by previous axes
   *  \param tgt_i_cell - target cell number accumulated by previous axes
   *  \param src_prev_area - factor by which this axis icreases cell number
   *  \param tgt_prev_area - factor by which this axis icreases cell number
   *  \param axis - the axis to treat
   *  \param prev_value - intersection size computed by previous axes
   */
  //================================================================================

  template <class MyMeshType, class MatrixType, int dim>
  void fillMatrix(const std::list< Interference >  inter_of_axis[dim],
                  MatrixType&                      result,
                  const int                        src_nb_cells[dim],
                  const int                        tgt_nb_cells[dim],
                  const int                        src_i_cell = 0,
                  const int                        tgt_i_cell = 0,
                  const int                        src_prev_area = 1,
                  const int                        tgt_prev_area = 1,
                  const int                        axis = 0,
                  const double                     prev_value = 1.0)
  {
    typedef std::list < Interference >::const_iterator TIntIterator;

    if ( axis + 1 == dim )
    {
      for ( TIntIterator i = inter_of_axis[axis].begin(); i != inter_of_axis[axis].end(); ++i )
      {
        double value = i->_length * prev_value;
        int src_i    = i->_src_index * src_prev_area + src_i_cell;
        int tgt_i    = i->_tgt_index * tgt_prev_area + tgt_i_cell;

        result[ tgt_i ].insert( std::make_pair( _TMI( src_i ), value ));
      }
    }
    else
    {
      int src_prev_area_next = src_prev_area * src_nb_cells[ axis ];
      int tgt_prev_area_next = tgt_prev_area * tgt_nb_cells[ axis ];

      for ( TIntIterator i = inter_of_axis[axis].begin(); i != inter_of_axis[axis].end(); ++i )
      {
        double value = i->_length * prev_value;
        int src_i    = i->_src_index * src_prev_area + src_i_cell;
        int tgt_i    = i->_tgt_index * tgt_prev_area + tgt_i_cell;

        // call for the next axis
        fillMatrix<MyMeshType, MatrixType, dim>(inter_of_axis, result,
                                                src_nb_cells, tgt_nb_cells, src_i, tgt_i,
                                                src_prev_area_next, tgt_prev_area_next,
                                                axis+1, value );
      }
    }
  }

  //================================================================================
  /*!
   * \brief Calculates the matrix of volumes of intersection between the elements of
   *        src_mesh and the elements of targetMesh
   *  \param src_mesh - source mesh
   *  \param tgt_mesh - target mesh
   *  \param result - matrix in which the result is stored 
   *  \param method - interpolation method, not used as only "P0P0" is implemented so far
   * 
   * The matrix is partially sparse : it is a vector of maps of integer - double pairs. 
   * It can also be an INTERP_KERNEL::Matrix object.
   * The length of the vector is equal to the number of target elements - for each target
   * element there is a map, regardless of whether the element intersects any source
   * elements or not. But in the maps there are only entries for those source elements
   * which have a non-zero intersection volume with the target element. The vector has
   * indices running from 0 to (nb target elements - 1), meaning that the map for target
   * element i is stored at index i - 1. In the maps, however, the indexing is more natural:
   * the intersection volume of the target element i with source element j is found at matrix[i-1][j]
   */
  //================================================================================

  template<class MyMeshType, class MatrixType>
  int InterpolationCC::interpolateMeshes(const MyMeshType& src_mesh,
                                         const MyMeshType& tgt_mesh,
                                         MatrixType&       result,
                                         const char *      method)
  {
    if ( std::string("P0P0") != method )
      throw Exception("Only P0P0 method is implemented so far");

    // create empty maps for all target elements
    result.resize( tgt_mesh.getNumberOfElements() );

    const int ret = src_mesh.getNumberOfElements();

    const double eps = getPrecision();
    const int dim = MyMeshType::MY_MESHDIM;
    //const NumberingPolicy numPol = MyMeshType::My_numPol;

    const double* src_coords[ dim ];
    const double* tgt_coords[ dim ];
    int src_nb_cells[ dim ];
    int tgt_nb_cells[ dim ];
    for ( int j = 0; j < dim; ++j )
    {
      src_coords[ j ] = src_mesh.getCoordsAlongAxis( _TMI( j ));
      tgt_coords[ j ] = tgt_mesh.getCoordsAlongAxis( _TMI( j ));
      src_nb_cells[ j ] = src_mesh.nbCellsAlongAxis( _TMI( j ));
      tgt_nb_cells[ j ] = tgt_mesh.nbCellsAlongAxis( _TMI( j ));
    }
    
    // ============================================
    // Calculate cell interferences along the axes
    // ============================================

    std::list < Interference > interferences[ dim ];

    for ( int j = 0; j < dim; ++j ) // loop on axes of castesian space
    {
      std::list < Interference >& axis_interferences = interferences[j];

      int it = 0, is = 0;
      double x1t, x2t, x1s, x2s; // left and right ordinates of target and source cells

      // look for the first interference
      // --------------------------------
      bool intersection = false;
      while ( !intersection && it < tgt_nb_cells[j] && is < src_nb_cells[j] )
      {
        x1s = src_coords[ j ][ is ];
        x2t = tgt_coords[ j ][ it+1 ];
        if ( x2t < x1s+eps )
        {
          it++; // source lays on the right of target
          continue;
        }
        x1t = tgt_coords[ j ][ it ];
        x2s = src_coords[ j ][ is+1 ];
        if ( x2s < x1t+eps )
        {
          is++; // source lays on the left of target
          continue;
        }
        intersection = true;
      }
      if ( !intersection ) return ret; // no intersections

      // get all interferences
      // ----------------------
      while ( intersection )
      {
        x1s = src_coords[ j ][ is ];
        x1t = tgt_coords[ j ][ it ];
        x2t = tgt_coords[ j ][ it+1 ];
        x2s = src_coords[ j ][ is+1 ];

        double x1 = std::max( x1s ,x1t );
        double x2 = std::min( x2s ,x2t );
        axis_interferences.push_back( Interference( is, it, x2 - x1 ));

        // to the next target and/or source cell
        double diff2 = x2s - x2t;
        if ( diff2 > -eps )
          intersection = ( ++it < tgt_nb_cells[j] );
        if ( diff2 < eps )
          intersection = ( ++is < src_nb_cells[j] && intersection);
      }
    }

    // ================
    // Fill the matrix
    // ================

    switch ( dim )
    {
    case 3:
      fillMatrix<MyMeshType,MatrixType,3>( interferences, result, src_nb_cells,tgt_nb_cells );
      break;

    case 2:
      fillMatrix<MyMeshType,MatrixType,2>( interferences, result, src_nb_cells,tgt_nb_cells );
      break;

    case 1:
      fillMatrix<MyMeshType,MatrixType,1>( interferences, result, src_nb_cells,tgt_nb_cells );
      break;
    }

    return ret;
  }
}
