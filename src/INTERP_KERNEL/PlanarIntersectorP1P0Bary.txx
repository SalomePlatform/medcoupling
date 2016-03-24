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
#ifndef __PlanarIntersectorP1P0Bary_TXX__
#define __PlanarIntersectorP1P0Bary_TXX__

#include "PlanarIntersectorP1P0Bary.hxx"
#include "InterpolationUtils.hxx"

#define PLAN_INTERSECTOR PlanarIntersectorP1P0Bary<MyMeshType,MyMatrix,ConcreteP1P0Intersector>
#define PLAN_INTER_TEMPLATE template<class MyMeshType, class MyMatrix, class ConcreteP1P0Intersector>

namespace INTERP_KERNEL
{
  PLAN_INTER_TEMPLATE
  PLAN_INTERSECTOR::PlanarIntersectorP1P0Bary(const MyMeshType& meshT, const MyMeshType& meshS,
                                              double dimCaracteristic, double precision,
                                              double md3DSurf, double minDot3DSurf, double medianPlane,
                                              bool doRotate, int orientation, int printLevel):
    PlanarIntersector<MyMeshType,MyMatrix>(meshT,meshS,dimCaracteristic,precision,md3DSurf,minDot3DSurf,
                                           medianPlane,doRotate,orientation,printLevel)
  {
    // SPEC:
    // "Limitation. For the P1P0 barycentric improvement only triangle source cells in 2D and
    // tetrahedrons in 3D will be supported by interpolators. If a non
    // triangle/tetrahedron source cell is detected an INTERP_KERNEL::Exception should be thrown."

    // Check types of source elements here rather than in intersectCells() since a wrong type can be
    // found late after a long time of calculation.

    const unsigned long numSrcElems = meshS.getNumberOfElements();
    for(unsigned long i = 0 ; i < numSrcElems ; ++i)
      if ( meshS.getTypeOfElement( OTT<ConnType,numPol>::indFC( i )) != NORM_TRI3 )
        throw INTERP_KERNEL::Exception("P1P0 barycentric algorithm works only with triangular source meshes");
  }

  PLAN_INTER_TEMPLATE
  int PLAN_INTERSECTOR::getNumberOfRowsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getNumberOfElements();
  }

  PLAN_INTER_TEMPLATE
  int PLAN_INTERSECTOR::getNumberOfColsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getNumberOfNodes();
  }

  /*!
   * This method computes a value per each node of each source triangle for target.
   */
  PLAN_INTER_TEMPLATE
  void PLAN_INTERSECTOR::intersectCells(ConnType                     icellT,
                                        const std::vector<ConnType>& icellsS,
                                        MyMatrix&                    res)
  {
    int orientation=1;
    std::vector<double> srcTriaCoords, tgtCellCoords, tgtCellCoordsTmp, nodeCeffs;

    // target cell data
    PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(icellT),tgtCellCoords);
    std::vector<double> * tgtCoords = & tgtCellCoords;
    int tgtNbNodes = tgtCellCoords.size()/SPACEDIM;
    NormalizedCellType tT=PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getTypeOfElement(OTT<ConnType,numPol>::indFC(icellT));
    bool isTargetQuad=CellModel::GetCellModel(tT).isQuadratic();

    typename MyMatrix::value_type& resRow=res[icellT];

    // treat each source triangle
    for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++)
    {
      int iS=*iter;
      PlanarIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinates(OTT<ConnType,numPol>::indFC(iS),srcTriaCoords);
      const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectS+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[iS]);
      if(SPACEDIM==3)
      {
        tgtCellCoordsTmp = tgtCellCoords;
        tgtCoords = & tgtCellCoordsTmp;
        orientation=PlanarIntersector<MyMeshType,MyMatrix>::projectionThis(&tgtCellCoordsTmp[0], &srcTriaCoords[0],
                                                                           tgtNbNodes, 3);
      }
      //double surf=orientation*intersectGeometryWithQuadrangle(quadrangle,targetCellCoordsTmp,isTargetQuad);
      double surf=orientation*intersectGeoBary( *tgtCoords, isTargetQuad, &srcTriaCoords[0], nodeCeffs );
      surf=PlanarIntersector<MyMeshType,MyMatrix>::getValueRegardingOption(surf);
      if(surf!=0.)
      {
        for(int nodeIdS=0;nodeIdS<3;nodeIdS++)
        {
          ConnType curNodeS=startOfCellNodeConn[nodeIdS];
          typename MyMatrix::value_type::const_iterator iterRes=resRow.find(curNodeS);
          if(iterRes!=resRow.end())
          {
            nodeCeffs[nodeIdS] += iterRes->second;
            resRow.erase( curNodeS );
          }
          resRow.insert(std::make_pair(curNodeS,nodeCeffs[nodeIdS]));
        }
      }
    }
  }
}
#endif
