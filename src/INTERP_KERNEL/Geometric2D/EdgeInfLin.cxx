#include "EdgeInfLin.hxx"

using namespace INTERP_KERNEL;

EdgeInfLin::EdgeInfLin(Node *pointPassingThrough, double slope)
{
  _start=pointPassingThrough;
  _start->incrRef();
  _end=new Node((*_start)[0]+cos(slope),(*_start)[1]+sin(slope));
}
