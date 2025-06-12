// Copyright (C) 2022-2025  CEA, EDF
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

#pragma once

#include "BBTree.txx"
#include <memory>

/*!
 * Wrapper over BBTree to deal with ownership of bbox double array.
 */
template <int dim, class ConnType>
class BBTreeStandAlone
{
   private:
    std::unique_ptr<double[]> _bbox;
    BBTree<dim, ConnType> _effective;

   public:
    const static int dimension = dim;
    BBTreeStandAlone(std::unique_ptr<double[]> &&bbs, ConnType nbelems, double epsilon = BBTREE_DFT_EPSILON)
        : _bbox(std::move(bbs)), _effective(_bbox.get(), nullptr, 0, nbelems, epsilon)
    {
    }
    void getIntersectingElems(const double *bb, std::vector<ConnType> &elems) const
    {
        _effective.getIntersectingElems(bb, elems);
    }
    void getElementsAroundPoint(const double *xx, std::vector<ConnType> &elems) const
    {
        _effective.getElementsAroundPoint(xx, elems);
    }
};
