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
// Author : Adrien Bruneton (CEA/DEN)

#ifndef __MappedBarycentric2DIntersectorP1P1_TXX__
#define __MappedBarycentric2DIntersectorP1P1_TXX__

#include "MappedBarycentric2DIntersectorP1P1.hxx"
#include "PlanarIntersector.txx"
#include "CellModel.hxx"

#include "PointLocatorAlgos.txx"
#include "MeshUtils.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  MappedBarycentric2DIntersectorP1P1<MyMeshType,MyMatrix>::MappedBarycentric2DIntersectorP1P1(const MyMeshType& meshT, const MyMeshType& meshS,
                                                                        double dimCaracteristic, double md3DSurf, double minDot3DSurf,
                                                                        double medianPlane, double precision, int orientation):
    PlanarIntersector<MyMeshType,MyMatrix>(meshT,meshS,dimCaracteristic,precision,md3DSurf,minDot3DSurf,medianPlane,true,orientation,0)
  {
  }

  template<class MyMeshType, class MyMatrix>
  void MappedBarycentric2DIntersectorP1P1<MyMeshType,MyMatrix>::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    std::vector<double> CoordsT;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(OTT<ConnType,numPol>::indFC(icellT),CoordsT);
    int nbOfNodesT=CoordsT.size()/SPACEDIM;
    for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++)
      {
        NormalizedCellType tS=PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getTypeOfElement(OTT<ConnType,numPol>::indFC(*iter));
        if(tS!=NORM_QUAD4)
          throw INTERP_KERNEL::Exception("Invalid source cell detected for meshdim==2. Only QUAD4 supported !");
        std::vector<double> CoordsS;
        PlanarIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinates(OTT<ConnType,numPol>::indFC(*iter),CoordsS);
        std::vector<double> CoordsTTmp(CoordsT);
        if(SPACEDIM==3)
          PlanarIntersector<MyMeshType,MyMatrix>::projectionThis(&CoordsS[0],&CoordsTTmp[0],CoordsS.size()/SPACEDIM,nbOfNodesT);
        const ConnType *startOfCellNodeConnT=PlanarIntersector<MyMeshType,MyMatrix>::_connectT+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[icellT]);
        for(int nodeIdT=0;nodeIdT<nbOfNodesT;nodeIdT++)
          {
            typename MyMatrix::value_type& resRow=res[OTT<ConnType,numPol>::ind2C(startOfCellNodeConnT[nodeIdT])];
            if( PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg2D(&CoordsTTmp[nodeIdT*SPACEDIM],&CoordsS[0],4,PlanarIntersector<MyMeshType,MyMatrix>::_precision) )
              {
                double mco[2];  // mapped coordinates in the quad4
                std::vector<const double*> coo(4);
                coo[0]=&CoordsS[0]; coo[1]=&CoordsS[SPACEDIM]; coo[2]=&CoordsS[2*SPACEDIM]; coo[3]=&CoordsS[3*SPACEDIM];
                quad_mapped_coords(coo,&CoordsTTmp[nodeIdT*SPACEDIM],mco);

                // Now use the form function of the QUAD4 to map the field values
                double resLoc[4];
                // See QUAD4 standard connectivity and cuboid_mapped_coords() convention:
                resLoc[0] = (1.-mco[0]) * (1.-mco[1]);
                resLoc[1] = (1.-mco[0]) *   mco[1]   ;
                resLoc[2] =  mco[0]     *   mco[1]   ;
                resLoc[3] =  mco[0]     * (1.-mco[1]);

                const ConnType *startOfCellNodeConnS=PlanarIntersector<MyMeshType,MyMatrix>::_connectS+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[*iter]);
                for(int nodeIdS=0;nodeIdS<4;nodeIdS++)
                  {
                    if(fabs(resLoc[nodeIdS])>PlanarIntersector<MyMeshType,MyMatrix>::_precision)
                      {
                        ConnType curNodeSInCmode=OTT<ConnType,numPol>::coo2C(startOfCellNodeConnS[nodeIdS]);
                        typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(curNodeSInCmode));
                        if(iterRes==resRow.end())
                          resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(curNodeSInCmode),resLoc[nodeIdS]));
                        else
                          {
                            double val=(*iterRes).second+resLoc[nodeIdS];
                            resRow.erase(OTT<ConnType,numPol>::indFC(curNodeSInCmode));
                            resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(curNodeSInCmode),val));
                          }
                      }
                  }
              }
          }
      }
  }

  template<class MyMeshType, class MyMatrix>
  int MappedBarycentric2DIntersectorP1P1<MyMeshType,MyMatrix>::getNumberOfRowsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getNumberOfNodes();
  }

  template<class MyMeshType, class MyMatrix>
  int MappedBarycentric2DIntersectorP1P1<MyMeshType,MyMatrix>::getNumberOfColsOfResMatrix() const
  {
    return PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getNumberOfNodes();
  }
}

#endif
