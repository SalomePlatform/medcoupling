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

#ifndef __TARGETINTERSECTOR__HXX__
#define __TARGETINTERSECTOR__HXX__

#include "INTERPKERNELDefines.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  /**
   * \brief Abstract base class of Intersector classes. 
   * These classes represent a target element and calculate its intersection
   * with the source elements.
   */
  template<class MyMeshType, class MyMatrix>
  class TargetIntersector
  {
  public:
    typedef typename MyMeshType::MyConnType ConnType;
  public:
    /*!
     * Tool for cell intersection, result is always positive.
     * @param icellT id of cell in target mesh in \b C \b mode.
     * @param icellsS ids of cells in source mesh in \b C \b mode.
     * @param res is an IN/OUT parameter that represents the icellTth row in final matrix, fed with at most icellsS elements. 
     */
    virtual void intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res) = 0;

    virtual int getNumberOfRowsOfResMatrix() const = 0;
    virtual int getNumberOfColsOfResMatrix() const = 0;
    virtual ~TargetIntersector() { }
  };
}

#endif
