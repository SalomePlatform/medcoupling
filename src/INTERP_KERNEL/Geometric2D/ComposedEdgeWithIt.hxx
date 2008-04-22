#ifndef __COMPOSEDNODEWITHIT_HXX__
#define __COMPOSEDNODEWITHIT_HXX__

#include "ComposedEdge.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  class ComposedEdgeWithIt : public ComposedEdge
  {
  public:
    ComposedEdgeWithIt():_iterator(this) { }
    const IteratorOnComposedEdge& getIterator() const { return _iterator; }
    void setIterator(const IteratorOnComposedEdge& it) { _iterator=it; }
  protected:
    IteratorOnComposedEdge _iterator;
  };
}

#endif
