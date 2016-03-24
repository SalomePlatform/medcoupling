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
#ifndef __CurveIntersectorP1P1_TXX__
#define __CurveIntersectorP1P1_TXX__

#include "CurveIntersectorP1P1.hxx"
#include "CurveIntersector.txx"

#define BASE_INTERSECTOR  CurveIntersector<MyMeshType,MyMatrix>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  CurveIntersectorP1P1<MyMeshType,MyMatrix>
    ::CurveIntersectorP1P1(const MyMeshType& meshT, const MyMeshType& meshS,
                           double precision, double tolerance,
                           double medianLine, int printLevel):
    BASE_INTERSECTOR (meshT, meshS, precision, tolerance, medianLine, printLevel)
  {
  }

  template<class MyMeshType, class MyMatrix>
  int CurveIntersectorP1P1<MyMeshType,MyMatrix>
  ::getNumberOfRowsOfResMatrix() const
  {
    return BASE_INTERSECTOR::_meshT.getNumberOfNodes();
  }

  template<class MyMeshType, class MyMatrix>
  int CurveIntersectorP1P1<MyMeshType,MyMatrix>
  ::getNumberOfColsOfResMatrix() const
  {
    return BASE_INTERSECTOR::_meshS.getNumberOfNodes();
  }

  template<class MyMeshType, class MyMatrix>
  void CurveIntersectorP1P1<MyMeshType,MyMatrix>
  ::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    std::vector<typename BASE_INTERSECTOR::TDualSegment> segmentsT, segmentsS;
    BASE_INTERSECTOR::getDualSegments( icellT, BASE_INTERSECTOR::_meshT, segmentsT);
    for ( int t = 0; t < (int)segmentsT.size(); ++t )
      {
        typename MyMatrix::value_type& resRow = res[ OTT<ConnType,numPol>::ind2C( segmentsT[t]._nodeId )];
        for(typename std::vector<ConnType>::const_iterator
              iter=icellsS.begin(); iter!=icellsS.end(); iter++)
          {
            int iS = *iter;
            BASE_INTERSECTOR::getDualSegments( OTT<ConnType,numPol>::ind2C(iS),
                                               BASE_INTERSECTOR::_meshS, segmentsS);
            for ( int s = 0; s < (int)segmentsS.size(); ++s )
              {
                double surf = BASE_INTERSECTOR::intersectSegments(&segmentsT[t]._coords[0],
                                                                  &segmentsS[s]._coords[0]);
                if(surf!=0.)
                  {
                    int nS = segmentsS[s]._nodeId;
                    typename MyMatrix::value_type::const_iterator iterRes=resRow.find(nS);
                    if(iterRes==resRow.end())
                      resRow.insert(std::make_pair(nS,surf));
                    else
                      {
                        surf+=(*iterRes).second;
                        resRow.erase(nS);
                        resRow.insert(std::make_pair(nS,surf));
                      }
                  }
              }
          }
      }
  }
}
#undef BASE_INTERSECTOR

#endif
