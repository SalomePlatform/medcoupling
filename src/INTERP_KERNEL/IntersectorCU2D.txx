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
// File      : IntersectorCU2D.txx
// Created   : Thu Dec 17 14:17:49 2009
// Author    : Edward AGAPOV (eap)

#ifndef __IntersectorCU2D_TXX__
#define __IntersectorCU2D_TXX__

#include "IntersectorCU2D.hxx"
#include "IntersectorCU.txx"

#define IntersectorCU2D_TEMPLATE template<class MyCMeshType, class MyUMeshType, class MyMatrix>
#define INTERSECTOR_CU2D IntersectorCU2D<MyCMeshType, MyUMeshType, MyMatrix >
#define INTER_CU IntersectorCU<MyCMeshType,MyUMeshType,MyMatrix,IntersectorCU2D<MyCMeshType,MyUMeshType,MyMatrix> >


namespace INTERP_KERNEL
{
  IntersectorCU2D_TEMPLATE
  INTERSECTOR_CU2D::IntersectorCU2D(const MyCMeshType& meshS,
                                    const MyUMeshType& meshT):
    IntersectorCU<MyCMeshType, MyUMeshType, MyMatrix, IntersectorCU2D<MyCMeshType,MyUMeshType,MyMatrix> >( meshS, meshT ),
    _intersector(meshT, meshT, 0,0,0,0,0,0,0 )
  {
    if ( MyCMeshType::MY_SPACEDIM != 2 || MyCMeshType::MY_MESHDIM != 2 ||
         MyUMeshType::MY_SPACEDIM != 2 || MyUMeshType::MY_MESHDIM != 2 )
      throw Exception("IntersectorCU2D(): Invalid mesh dimension, it must be 2");
  }


  //================================================================================
  /*!
   * \brief Calculate area of intersection of an unstructured cell and a cartesian one.
   * The cartesian cell is given by its [i,j] indices
   */
  //================================================================================

  IntersectorCU2D_TEMPLATE
  double INTERSECTOR_CU2D::intersectGeometry(UConnType                     icellT,
                                             const std::vector<CConnType>& icellS)
  {
    std::vector<double> uCoords;
    this->getUCoordinates( icellT, uCoords );

    NormalizedCellType tT = INTER_CU::_meshU.getTypeOfElement( _TMIU(icellT));
    bool is_tgt_quad = CellModel::GetCellModel(tT).isQuadratic();

    double quad[8] = { INTER_CU::_coordsC[0][icellS[0]],   INTER_CU::_coordsC[1][icellS[1]],
                       INTER_CU::_coordsC[0][icellS[0]+1], INTER_CU::_coordsC[1][icellS[1]], 
                       INTER_CU::_coordsC[0][icellS[0]+1], INTER_CU::_coordsC[1][icellS[1]+1], 
                       INTER_CU::_coordsC[0][icellS[0]],   INTER_CU::_coordsC[1][icellS[1]+1] };

    double surf = _intersector.intersectGeometryWithQuadrangle( quad,
                                                                uCoords,
                                                                is_tgt_quad);
    return surf;
  }
}
#endif
