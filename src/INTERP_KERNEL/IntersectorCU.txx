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
// File      : IntersectorCU.txx
// Created   : Thu Dec 17 12:30:17 2009
// Author    : Edward AGAPOV (eap)
//
#ifndef __IntersectorCU_TXX__
#define __IntersectorCU_TXX__

#include "IntersectorCU.hxx"

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


#define _CU_TEMPLATE \
template<class MyCMeshType, class MyUMeshType, class MyMatrix, class ConcreteIntersector>
#define _INTERSECTOR_CU_ \
IntersectorCU<MyCMeshType, MyUMeshType, MyMatrix, ConcreteIntersector>

namespace INTERP_KERNEL
{
  //================================================================================
  /*!
   * \brief Constructor
   */
  //================================================================================

  _CU_TEMPLATE
  _INTERSECTOR_CU_::IntersectorCU(const MyCMeshType& meshS, const MyUMeshType& meshT):
    _meshU(meshT), _meshC(meshS)
  {
    _connectU  =meshT.getConnectivityPtr();
    _connIndexU=meshT.getConnectivityIndexPtr();
    _coordsU   =meshT.getCoordinatesPtr();

    for ( int j = 0; j < SPACEDIM; ++j )
      {
        _coordsC [ j ] = _meshC.getCoordsAlongAxis( _TMIC( j ));
        _nbCellsC[ j ] = _meshC.nbCellsAlongAxis  ( _TMIC( j ));
      }
  }

  //================================================================================
  /*!
   * \brief Destructor
   */
  //================================================================================

  _CU_TEMPLATE
  _INTERSECTOR_CU_::~IntersectorCU()
  {
  }

  //================================================================================
  /*!
   * \brief Return bounding box of an unstructured element
   */
  //================================================================================

  _CU_TEMPLATE
  void _INTERSECTOR_CU_::getUElemBB(double* bb, UConnType icell)
  {
    //initializing bounding box limits
    for(int idim=0; idim<SPACEDIM; idim++)
      {
        bb[2*idim  ] =  std::numeric_limits<double>::max();
        bb[2*idim+1] = -std::numeric_limits<double>::max();
      }

    UConnType nb_nodes = _connIndexU[icell+1] - _connIndexU[icell];
    for (UConnType i=0; i<nb_nodes; i++)
      {
        //const double* coord_node=_coordsU+SPACEDIM*(OTT<ConnType,numPol>::coo2C(conn[OTT<ConnType,numPol>::conn2C(conn_index[OTT<ConnType,numPol>::ind2C(iP)]+i)]));
        const double* coord_node=_coordsU+SPACEDIM*(_FMCOO( _connectU[_FMCON (_connIndexU[_FMIU(icell)]+i)]));
        for(int idim=0; idim<SPACEDIM; idim++)
          {            
            double x = *(coord_node+idim);
            bb[2*idim  ] = (x<bb[2*idim  ])?x:bb[2*idim  ];
            bb[2*idim+1] = (x>bb[2*idim+1])?x:bb[2*idim+1];
          }
      }
  }

  //================================================================================
  /*!
   * \brief Return coordinates of nodes of an unstructured element
   */
  //================================================================================

  _CU_TEMPLATE
  void _INTERSECTOR_CU_::getUCoordinates(UConnType icell, std::vector<double>& coords)
  {
    UConnType nb_nodes = _connIndexU[icell+1] - _connIndexU[icell];
    coords.resize( SPACEDIM * nb_nodes );
    for (UConnType i=0; i<nb_nodes; i++)
      for(int j=0; j<SPACEDIM; j++)
        coords[SPACEDIM*i+j]=_coordsU[SPACEDIM*_FMCOO( _connectU[_FMCON (_connIndexU[_FMIU(icell)]+i)])+j];
  }

  _CU_TEMPLATE
  int _INTERSECTOR_CU_::getNumberOfRowsOfResMatrix() const
  {
    return _meshU.getNumberOfElements();
  }
  _CU_TEMPLATE
  int _INTERSECTOR_CU_::getNumberOfColsOfResMatrix() const
  {
    return _meshC.getNumberOfElements();
  }

  //================================================================================
  /*!
   * \brief Intersects unstructured target cell with structured source cell given by
   *        its indices
   *  \param icellU - index of an unstructured element
   *  \param icellC - i,j,k of a structured element
   *  \param res - matrix to fill in
   */
  //================================================================================

  _CU_TEMPLATE
  void _INTERSECTOR_CU_::intersectCells(CConnType                     icellU,
                                        const std::vector<CConnType>& icellC,
                                        MyMatrix&                     res)
  {
    double v = intersectGeometry(icellU, icellC);

    CConnType iC = icellC[0], area = _nbCellsC[0];
    for ( int j = 1; j < SPACEDIM; ++j )
      {
        iC += icellC[j] * area;
        area *= _nbCellsC[j];
      }
    res[ icellU ][ iC ] = v;
  }
}

#endif
