#ifndef __ABSTRACTEDGE_HXX__
#define __ABSTRACTEDGE_HXX__

#include "Geometric2D_defines.hxx"

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
  class GEOMETRIC2D_EXPORT IteratorOnComposedEdge
  {
    friend class ComposedEdge;
    friend class ElementaryEdge;
    friend class QuadraticPolygon;
  public:
    IteratorOnComposedEdge();
    IteratorOnComposedEdge(ComposedEdge *compEdges);
    bool isValid() const { return _listHandle!=0; } 
    void operator=(const IteratorOnComposedEdge& other);
    void first() { _deepIt=_listHandle->begin(); }
    void next() { _deepIt++; }
    void last();
    void nextLoop();
    void previousLoop();
    bool finished() const { return _deepIt==_listHandle->end(); }
    bool goToNextInOn(bool direction, int& i, int nbMax);
    ElementaryEdge *current() { return *_deepIt; }
    void assignMySelfToAllElems(ComposedEdge *elems);
    void insertElemEdges(ComposedEdge *elems, bool changeMySelf);
  private:
    std::list<ElementaryEdge *>::iterator _deepIt;
    std::list<ElementaryEdge *>* _listHandle;
  };
}

#endif
