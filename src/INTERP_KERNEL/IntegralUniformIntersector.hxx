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

#ifndef __INTEGRALUNIFORMINTERSECTOR_HXX__
#define __INTEGRALUNIFORMINTERSECTOR_HXX__

#include "TargetIntersector.hxx"

#include <cmath>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  class IntegralUniformIntersector : public TargetIntersector<MyMeshType,MyMatrix>
  {
  public:
    typedef typename MyMeshType::MyConnType ConnType;
  public:
    IntegralUniformIntersector(const MyMeshType& mesh, bool isAbs);
    double performNormalization(double val) const { if(_is_abs) return fabs(val); else return val; }
    void setFromTo(bool val) { _from_to=val; }
    void putValueIn(ConnType i, double val, MyMatrix& res) const;
  protected:
    const MyMeshType& _mesh;
    //! if false means fromIntegralUniform if true means toIntegralUniform
    bool _from_to;
    bool _is_abs;
  };

  template<class MyMeshType, class MyMatrix>
  class IntegralUniformIntersectorP0 : public IntegralUniformIntersector<MyMeshType,MyMatrix>
  {
  public:
    typedef typename MyMeshType::MyConnType ConnType;
  public:
    IntegralUniformIntersectorP0(const MyMeshType& mesh, bool isAbs);
    int getNumberOfRowsOfResMatrix() const;
    int getNumberOfColsOfResMatrix() const;
    void intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res);
  };

  template<class MyMeshType, class MyMatrix>
  class IntegralUniformIntersectorP1 : public IntegralUniformIntersector<MyMeshType,MyMatrix>
  {
  public:
    typedef typename MyMeshType::MyConnType ConnType;
  public:
    IntegralUniformIntersectorP1(const MyMeshType& mesh, bool isAbs);
    int getNumberOfRowsOfResMatrix() const;
    int getNumberOfColsOfResMatrix() const;
    void intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res);
  };
}

#endif
