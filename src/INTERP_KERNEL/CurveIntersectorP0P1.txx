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
#ifndef __CurveIntersectorP0P1_TXX__
#define __CurveIntersectorP0P1_TXX__

#include "CurveIntersectorP0P1.hxx"
#include "CurveIntersector.txx"

#define BASE_INTERSECTOR CurveIntersector<MyMeshType,MyMatrix>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  CurveIntersectorP0P1<MyMeshType,MyMatrix>
    ::CurveIntersectorP0P1(const MyMeshType& meshT, const MyMeshType& meshS,
                           double precision, double tolerance,
                           double medianLine, int printLevel):
      BASE_INTERSECTOR(meshT, meshS, precision, tolerance, medianLine, printLevel)
  {
  }

  template<class MyMeshType, class MyMatrix>
  int CurveIntersectorP0P1<MyMeshType,MyMatrix>
    ::getNumberOfRowsOfResMatrix() const
  {
    return BASE_INTERSECTOR::_meshT.getNumberOfNodes();
  }

  template<class MyMeshType, class MyMatrix>
  int CurveIntersectorP0P1<MyMeshType,MyMatrix>
    ::getNumberOfColsOfResMatrix() const
  {
    return BASE_INTERSECTOR::_meshS.getNumberOfElements();
  }

  //================================================================================
  /*!
   * \brief Project from segments to nodes
   */
  //================================================================================

  template<class MyMeshType, class MyMatrix>
  void CurveIntersectorP0P1<MyMeshType,MyMatrix>
  ::intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res)
  {
    std::vector<typename BASE_INTERSECTOR::TDualSegment> segmentsT;
    BASE_INTERSECTOR::getDualSegments( icellT, BASE_INTERSECTOR::_meshT, segmentsT);
    for ( int t = 0; t < (int)segmentsT.size(); ++t )
      {
        typename MyMatrix::value_type& resRow = res[ OTT<ConnType,numPol>::ind2C( segmentsT[t]._nodeId )];
        for(typename std::vector<ConnType>::const_iterator
              iter=icellsS.begin(); iter!=icellsS.end(); iter++)
          {
            int iS = *iter;
            std::vector<double> coordsS;
            int s, nbSegS = 1 + BASE_INTERSECTOR::getRealSourceCoordinates(iS,coordsS);
            for ( s = 0; s < nbSegS; ++s )
              {
                double surf = BASE_INTERSECTOR::intersectSegments(&segmentsT[t]._coords[0],
                                                                  &coordsS[0] + s*SPACEDIM);
                if(surf!=0.)
                  {
                    typename MyMatrix::value_type::const_iterator iterRes=resRow.find(OTT<ConnType,numPol>::indFC(iS));
                    if(iterRes==resRow.end())
                      resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(iS),surf));
                    else
                      {
                        surf+=(*iterRes).second;
                        resRow.erase(OTT<ConnType,numPol>::indFC(iS));
                        resRow.insert(std::make_pair(OTT<ConnType,numPol>::indFC(iS),surf));
                      }
                  }
              }
          }
      }
  }
}
#undef BASE_INTERSECTOR

#endif
