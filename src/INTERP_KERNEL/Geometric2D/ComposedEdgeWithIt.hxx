#ifndef __COMPOSEDNODEWITHIT_HXX__
#define __COMPOSEDNODEWITHIT_HXX__

#include "Geometric2D_defines.hxx"

#include "ComposedEdge.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  class GEOMETRIC2D_EXPORT ComposedEdgeWithIt : public ComposedEdge
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
