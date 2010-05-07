//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef __INTERPKERNELGEO2DABSTRACTEDGE_HXX__
#define __INTERPKERNELGEO2DABSTRACTEDGE_HXX__

#include "INTERPKERNELGEOMETRIC2DDefines.hxx"

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
  class INTERPKERNELGEOMETRIC2D_EXPORT IteratorOnComposedEdge
  {
    friend class ComposedEdge;
    friend class ElementaryEdge;
    friend class QuadraticPolygon;
  public:
    IteratorOnComposedEdge();
    IteratorOnComposedEdge(ComposedEdge *compEdges);
    bool isValid() const { return _list_handle!=0; } 
    void operator=(const IteratorOnComposedEdge& other);
    void first() { _deep_it=_list_handle->begin(); }
    void next() { _deep_it++; }
    void last();
    void nextLoop();
    void previousLoop();
    bool finished() const { return _deep_it==_list_handle->end(); }
    bool goToNextInOn(bool direction, int& i, int nbMax);
    ElementaryEdge *current() { return *_deep_it; }
    void assignMySelfToAllElems(ComposedEdge *elems);
    void insertElemEdges(ComposedEdge *elems, bool changeMySelf);
  private:
    std::list<ElementaryEdge *>::iterator _deep_it;
    std::list<ElementaryEdge *>* _list_handle;
  };
}

#endif
