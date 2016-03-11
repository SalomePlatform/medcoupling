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
// File      : IntersectorCU1D.txx
// Created   : Thu Dec 17 14:17:49 2009
// Author    : Edward AGAPOV (eap)

#ifndef __IntersectorCU1D_TXX__
#define __IntersectorCU1D_TXX__

#include "IntersectorCU1D.hxx"
#include "IntersectorCU.txx"

#define  IntersectorCU1D_TEMPLATE template<class MyCMeshType, class MyUMeshType, class MyMatrix>
#define  INTERSECTOR_CU1D IntersectorCU1D<MyCMeshType,MyUMeshType,MyMatrix >
#define _INTER_CU         IntersectorCU  <MyCMeshType,MyUMeshType,MyMatrix,IntersectorCU1D<MyCMeshType,MyUMeshType,MyMatrix> >

namespace INTERP_KERNEL
{
  //================================================================================
  /*!
   * \brief intersector of the unstructured mesh and the cartesian mesh in 1D
   */
  //================================================================================

  IntersectorCU1D_TEMPLATE
  INTERSECTOR_CU1D::IntersectorCU1D(const MyCMeshType& meshS,
                                    const MyUMeshType& meshT):
    _INTER_CU( meshS, meshT )
  {
    if ( MyCMeshType::MY_SPACEDIM != 1 || MyCMeshType::MY_MESHDIM != 1 ||
         MyUMeshType::MY_SPACEDIM != 1 || MyUMeshType::MY_MESHDIM != 1 )
      throw Exception("IntersectorCU1D(): Invalid mesh dimension, it must be 1");
  }

  //================================================================================
  /*!
   * \brief destructor
   */
  //================================================================================

  IntersectorCU1D_TEMPLATE
  INTERSECTOR_CU1D::~IntersectorCU1D()
  {
  }

  //================================================================================
  /*!
   * \brief Calculate length of intersection of an unstructured cell and a cartesian one.
   * The cartesian cell is given by its [i,j,k] indices
   */
  //================================================================================

  IntersectorCU1D_TEMPLATE
  double INTERSECTOR_CU1D::intersectGeometry(UConnType                     icellT,
                                             const std::vector<CConnType>& icellS)
  {
    std::vector<double> coordsU;
    _INTER_CU::getUCoordinates(icellT, coordsU);

    const double* coordsC = & _INTER_CU::_coordsC[0][ _FMIC(icellS[0]) ];

    double res = std::min( coordsU[1], coordsC[1] ) - std::max( coordsU[0], coordsC[0] );
    return res;
  }
}
#endif
