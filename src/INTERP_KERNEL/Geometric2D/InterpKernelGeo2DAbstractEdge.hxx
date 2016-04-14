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

#ifndef __INTERPKERNELGEO2DABSTRACTEDGE_HXX__
#define __INTERPKERNELGEO2DABSTRACTEDGE_HXX__

#include "INTERPKERNELDefines.hxx"

#include <set>
#include <list>
#include <fstream>

namespace INTERP_KERNEL
{
  class Edge;
  class Node;
  class Bounds;

  class ComposedEdge;
  class ElementaryEdge;

  /*!
   * Asumption is done with this iterator that we iterate on a container containing more than one edge.
   */
  class IteratorOnComposedEdge
  {
    friend class ComposedEdge;
    friend class ElementaryEdge;
    friend class QuadraticPolygon;
  public:
    INTERPKERNEL_EXPORT IteratorOnComposedEdge();
    INTERPKERNEL_EXPORT IteratorOnComposedEdge(ComposedEdge *compEdges);
    INTERPKERNEL_EXPORT bool isValid() const { return _list_handle!=0; } 
    INTERPKERNEL_EXPORT void operator=(const IteratorOnComposedEdge& other);
    INTERPKERNEL_EXPORT void first() { _deep_it=_list_handle->begin(); }
    INTERPKERNEL_EXPORT void next() { _deep_it++; }
    INTERPKERNEL_EXPORT void last();
    INTERPKERNEL_EXPORT void nextLoop();
    INTERPKERNEL_EXPORT void previousLoop();
    INTERPKERNEL_EXPORT bool finished() const { return _deep_it==_list_handle->end(); }
    INTERPKERNEL_EXPORT bool goToNextInOn(bool direction, int& i, int nbMax);
    INTERPKERNEL_EXPORT ElementaryEdge *current() { return *_deep_it; }
    INTERPKERNEL_EXPORT void assignMySelfToAllElems(ComposedEdge *elems);
    INTERPKERNEL_EXPORT void insertElemEdges(ComposedEdge *elems, bool changeMySelf);
  private:
    std::list<ElementaryEdge *>::iterator _deep_it;
    std::list<ElementaryEdge *>* _list_handle;
  };
}

#endif
