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
// Author : Anthony Geay (EDF R&D)

#ifndef __CURVEINTERSECTORP1P1PL_TXX__
#define __CURVEINTERSECTORP1P1PL_TXX__

#include "CurveIntersectorP1P1PL.hxx"
#include "CurveIntersector.txx"

#include <cassert>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  CurveIntersectorP1P1PL<MyMeshType,MyMatrix>::CurveIntersectorP1P1PL(const MyMeshType& meshT, const MyMeshType& meshS, double precision, double tolerance, double medianLine, int printLevel):CurveIntersector<MyMeshType,MyMatrix>(meshT, meshS, precision, tolerance, medianLine, printLevel)
  {
  }

  template<class MyMeshType, class MyMatrix>
  int CurveIntersectorP1P1PL<MyMeshType,MyMatrix>::getNumberOfRowsOfResMatrix() const
  {
    return CurveIntersector<MyMeshType,MyMatrix>::_meshT.getNumberOfNodes();
  }

  template<class MyMeshType, class MyMatrix>
  int CurveIntersectorP1P1PL<MyMeshType,MyMatrix>::getNumberOfColsOfResMatrix() const
  {
    return CurveIntersector<MyMeshType,MyMatrix>::_meshS.getNumberOfNodes();
  }

  template<class MyMeshType, class MyMatrix>
  void CurveIntersectorP1P1PL<MyMeshType,MyMatrix>::AppendValueInMatrix2(typename MyMatrix::value_type& resRow, ConnType nodeIdS0, double val0)
  {
    typename MyMatrix::value_type::const_iterator iterRes(resRow.find(OTT<ConnType,numPol>::indFC(nodeIdS0)));    
    if(iterRes==resRow.end())
      {
        resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(nodeIdS0),val0));
      }
    else
      {
        double val((*iterRes).second+val0);
        resRow.erase(OTT<ConnType,numPol>::indFC(nodeIdS0));
        resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(nodeIdS0),val));
      }
  }

  template<class MyMeshType, class MyMatrix>
  void CurveIntersectorP1P1PL<MyMeshType,MyMatrix>::AppendValueInMatrix(MyMatrix& res, ConnType nodeIdT, ConnType nodeIdS0, double val0, ConnType nodeIdS1, double val1)
  {
    typename MyMatrix::value_type& resRow(res[nodeIdT]);
    AppendValueInMatrix2(resRow,nodeIdS0,val0);
    AppendValueInMatrix2(resRow,nodeIdS1,val1);
  }

  template<class MyMeshType, class MyMatrix>
  void CurveIntersectorP1P1PL<MyMeshType,MyMatrix>::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    std::vector<double> coordsT;
    if(CurveIntersector<MyMeshType,MyMatrix>::getRealTargetCoordinates(icellT,coordsT))
      throw INTERP_KERNEL::Exception("Invalid target cell detected for meshdim==1. Only SEG2 supported !");
    assert(coordsT.size()/SPACEDIM==2);
    for(typename std::vector<ConnType>::const_iterator iter=icellsS.begin();iter!=icellsS.end();iter++)
      {
        std::vector<double> coordsS;
        if(CurveIntersector<MyMeshType,MyMatrix>::getRealSourceCoordinates(*iter,coordsS))
          throw INTERP_KERNEL::Exception("Invalid source cell detected for meshdim==1. Only SEG2 supported !");
        assert(coordsS.size()/SPACEDIM==2);
        double xs0,xs1,xt0,xt1;
        double lgth(CurveIntersector<MyMeshType,MyMatrix>::intersectSegmentsInternal(&coordsT[0],&coordsS[0],xs0,xs1,xt0,xt1));
        ConnType nodeIdS0(CurveIntersector<MyMeshType,MyMatrix>::getNodeIdOfSourceCellAt(*iter,0));
        ConnType nodeIdS1(CurveIntersector<MyMeshType,MyMatrix>::getNodeIdOfSourceCellAt(*iter,1));
        if(lgth>0.)
          {
            double a,b;
            // for first
            ConnType nodeIdT0(CurveIntersector<MyMeshType,MyMatrix>::getNodeIdOfTargetCellAt(icellT,0));
            if(CurveIntersector<MyMeshType,MyMatrix>::ComputeBaryCoordsOf(xs0,xs1,xt0,a,b))
              {
                a*=lgth; b*=lgth;
                AppendValueInMatrix(res,nodeIdT0,nodeIdS0,a,nodeIdS1,b);
              }
            //
            ConnType nodeIdT1(CurveIntersector<MyMeshType,MyMatrix>::getNodeIdOfTargetCellAt(icellT,1));
            typename MyMatrix::value_type& resRow1=res[nodeIdT1];
            if(CurveIntersector<MyMeshType,MyMatrix>::ComputeBaryCoordsOf(xs0,xs1,xt1,a,b))
              {
                a*=lgth; b*=lgth;
                AppendValueInMatrix(res,nodeIdT1,nodeIdS0,a,nodeIdS1,b);
              }
          }
      }
  }
}

#endif
