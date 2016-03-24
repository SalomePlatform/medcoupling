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

#ifndef __CURVEINTERSECTORP1P1PL_HXX__
#define __CURVEINTERSECTORP1P1PL_HXX__

#include "CurveIntersector.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  class CurveIntersectorP1P1PL : public CurveIntersector<MyMeshType,MyMatrix>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;

    CurveIntersectorP1P1PL(const MyMeshType& meshT, const MyMeshType& meshS,
                           double  precision, double tolerance,
                           double medianLine, int printLevel);
  public:
    void intersectCells(ConnType icellT, const std::vector<ConnType>& icellsS, MyMatrix& res);
    int getNumberOfRowsOfResMatrix() const;
    int getNumberOfColsOfResMatrix() const;
  private:
    static void AppendValueInMatrix(MyMatrix& res, ConnType nodeIdT, ConnType nodeIdS0, double val0, ConnType nodeIdS1, double val1);
    static void AppendValueInMatrix2(typename MyMatrix::value_type& resRow, ConnType nodeIdS0, double val0);
  };
}

#endif
