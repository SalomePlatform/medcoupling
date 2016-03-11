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
#ifndef __CurveIntersectorP1P0_TXX__
#define __CurveIntersectorP1P0_TXX__

#include "CurveIntersectorP1P0.hxx"
#include "CurveIntersector.txx"

#define BASE_INTERSECTOR CurveIntersector<MyMeshType,MyMatrix>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  CurveIntersectorP1P0<MyMeshType,MyMatrix>
    ::CurveIntersectorP1P0(const MyMeshType& meshT, const MyMeshType& meshS,
                           double precision, double tolerance,
                           double medianLine, int printLevel):
    BASE_INTERSECTOR (meshT, meshS, precision, tolerance, medianLine, printLevel)
  {
  }

  template<class MyMeshType, class MyMatrix>
  int CurveIntersectorP1P0<MyMeshType,MyMatrix>
    ::getNumberOfRowsOfResMatrix() const
  {
    return BASE_INTERSECTOR::_meshT.getNumberOfElements();
  }

  template<class MyMeshType, class MyMatrix>
  int CurveIntersectorP1P0<MyMeshType,MyMatrix>
    ::getNumberOfColsOfResMatrix() const
  {
    return BASE_INTERSECTOR::_meshS.getNumberOfNodes();
  }

  //================================================================================
  /*!
   * \brief Project from source nodes to target segments
   */
  //================================================================================

  template<class MyMeshType, class MyMatrix>
  void CurveIntersectorP1P0<MyMeshType,MyMatrix>
  ::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    typename MyMatrix::value_type& resRow = res[ icellT ];
    std::vector<typename BASE_INTERSECTOR::TDualSegment> segmentsS;
    std::vector<double> coordsT;
    int t, nbSegT = 1 + BASE_INTERSECTOR::getRealTargetCoordinates(icellT,coordsT);
    for ( t = 0; t < nbSegT; ++t )
      for(typename std::vector<ConnType>::const_iterator
            iter=icellsS.begin(); iter!=icellsS.end(); iter++)
        {
          int iS = *iter;
          BASE_INTERSECTOR::getDualSegments( OTT<ConnType,numPol>::ind2C(iS),
                                             BASE_INTERSECTOR::_meshS, segmentsS);
          for ( int s = 0; s < (int)segmentsS.size(); ++s )
            {
              double surf = BASE_INTERSECTOR::intersectSegments(&segmentsS[s]._coords[0],
                                                                &coordsT[0] + t*SPACEDIM);
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
#undef BASE_INTERSECTOR

#endif
