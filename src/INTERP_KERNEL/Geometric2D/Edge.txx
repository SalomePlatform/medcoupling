#ifndef __EDGE_TXX__
#define __EDGE_TXX__

#include "EdgeArcCircle.hxx"

template<INTERP_KERNEL::TypeOfMod4QuadEdge type>
INTERP_KERNEL::Edge *INTERP_KERNEL::Edge::buildEdgeFrom(Node *start, Node *middle, Node *end)
{
  return new INTERP_KERNEL::EdgeArcCircle(start,middle,end);
}

#endif
