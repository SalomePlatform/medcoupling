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

#ifndef __PlanarIntersectorP0P1Bary_TXX__
#define __PlanarIntersectorP0P1Bary_TXX__

#include "PlanarIntersectorP0P1Bary.hxx"
#include "InterpolationUtils.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, class ConcreteP0P1Intersector>
  PlanarIntersectorP0P1Bary<MyMeshType,MyMatrix,ConcreteP0P1Intersector>::PlanarIntersectorP0P1Bary(const MyMeshType& meshT, const MyMeshType& meshS,
                                                                                                    double dimCaracteristic, double precision,
                                                                                                    double md3DSurf, double minDot3DSurf, double medianPlane,
                                                                                                    bool doRotate, int orientation, int printLevel):
  PlanarIntersector<MyMeshType,MyMatrix>(meshT,meshS,dimCaracteristic,precision,md3DSurf,minDot3DSurf,
                                         medianPlane,doRotate,orientation,printLevel)
  {
    // SPEC:
    // "Limitation. For the P0P1 barycentric improvement only triangle target cells in 2D and
    // tetrahedrons in 3D will be supported by interpolators. If a non
    // triangle/tetrahedron source cell is detected an INTERP_KERNEL::Exception should be thrown."

    // Check types of source elements here rather than in intersectCells() since a wrong type can be
    // found late after a long time of calculation.

    const unsigned long numTrgElems = meshT.getNumberOfElements();
    for(unsigned long i = 0 ; i < numTrgElems ; ++i)
      if ( meshT.getTypeOfElement( OTT<ConnType,numPol>::indFC( i )) != NORM_TRI3 )
        throw INTERP_KERNEL::Exception("P0P1 barycentric algorithm works only with triangular target meshes");
  }

  template<class MyMeshType, class MyMatrix, class ConcreteP0P1Intersector>
  int PlanarIntersectorP0P1Bary<MyMeshType,MyMatrix,ConcreteP0P1Intersector>::getNumberOfRowsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getNumberOfNodes();
  }

  template<class MyMeshType, class MyMatrix, class ConcreteP0P1Intersector>
  int PlanarIntersectorP0P1Bary<MyMeshType,MyMatrix,ConcreteP0P1Intersector>::getNumberOfColsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getNumberOfElements();
  }

  /*!
   * This method computes a value per each node of each source triangle for target.
   */
  template<class MyMeshType, class MyMatrix, class ConcreteP0P1Intersector>
  void PlanarIntersectorP0P1Bary<MyMeshType,MyMatrix,ConcreteP0P1Intersector>::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    int orientation=1;
    std::vector<double> trgTriaCoords,trgTriaCoordsTmp;
    // target cell data
    PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(icellT),trgTriaCoords);
    std::vector<double> *tgtCoords(&trgTriaCoords);
    const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectT+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[icellT]);
    // treat each source cells
    for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++)
    {
      std::vector<double> srcCellCoords,srcCellCoordsTmp,nodeCeffs;
      int iS=*iter;
      NormalizedCellType tS=PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getTypeOfElement(OTT<ConnType,numPol>::indFC(iS));
      bool isSourceQuad=CellModel::GetCellModel(tS).isQuadratic();
      PlanarIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinates(OTT<ConnType,numPol>::indFC(iS),srcCellCoords);
      std::vector<double> *srcCoords(&srcCellCoords);
      int srcNbNodes = srcCellCoords.size()/SPACEDIM;
      if(SPACEDIM==3)
        {
          srcCellCoordsTmp=srcCellCoords;
          trgTriaCoordsTmp=trgTriaCoords;
          srcCoords=&srcCellCoordsTmp;
          tgtCoords=&trgTriaCoordsTmp;
          orientation=PlanarIntersector<MyMeshType,MyMatrix>::projectionThis(&trgTriaCoordsTmp[0],&srcCellCoordsTmp[0],
                                                                             3,srcNbNodes);
        }
      //double surf=orientation*intersectGeometryWithQuadrangle(quadrangle,targetCellCoordsTmp,isTargetQuad);
      double surf=orientation*intersectGeoBary(*srcCoords,isSourceQuad,&((*tgtCoords)[0]),nodeCeffs);
      surf=PlanarIntersector<MyMeshType,MyMatrix>::getValueRegardingOption(surf);
      if(surf!=0.)
      {
        for(int nodeIdT=0;nodeIdT<3;nodeIdT++)
        {
          ConnType curNodeT=startOfCellNodeConn[nodeIdT];
          typename MyMatrix::value_type& resRow=res[curNodeT];
          typename MyMatrix::value_type::const_iterator iterRes=resRow.find(*iter);
          if(iterRes!=resRow.end())
          {
            nodeCeffs[*iter] += iterRes->second;
            resRow.erase(*iter);
          }
          resRow.insert(std::make_pair(*iter,nodeCeffs[nodeIdT]));
        }
      }
    }
  }
}
#endif
