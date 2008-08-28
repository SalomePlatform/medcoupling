#ifndef __EDGEINFLIN_HXX__
#define __EDGEINFLIN_HXX__

#include "Geometric2D_defines.hxx"

#include "EdgeLin.hxx"

namespace INTERP_KERNEL
{
  class GEOMETRIC2D_EXPORT EdgeInfLin : public EdgeLin
  {
  public:
    EdgeInfLin(Node *start, Node *end):EdgeLin(start,end,true) { }
    EdgeInfLin(Node *pointPassingThrough, double slope);
    bool isIn(double characterVal) const { return true; }
    void dynCastFunction(const EdgeLin * &seg,
                         const EdgeArcCircle * &arcSeg) const { seg=this; }
    void dumpInXfigFile(std::ostream& stream) const { }
  private:
    ~EdgeInfLin() { }
  };
}

#endif
