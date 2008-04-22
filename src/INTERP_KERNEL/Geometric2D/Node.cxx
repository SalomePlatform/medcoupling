#include "Node.hxx"
#include "EdgeArcCircle.hxx"

using namespace std;
using namespace INTERP_KERNEL;

Node::Node(double x, double y):_isToDel(true),_cnt(1),_loc(UNKNOWN)
{
  const unsigned SPACEDIM=2;
  _coords=new double[SPACEDIM];
  _coords[0]=x; _coords[1]=y;
}

Node::Node(const double *coords):_isToDel(false),_cnt(1),_loc(UNKNOWN),_coords((double *)coords)
{
}

Node::Node(std::istream& stream):_isToDel(true),_cnt(1),_loc(UNKNOWN)
{
  const unsigned SPACEDIM=2;
  _coords=new double[SPACEDIM];
  for(unsigned i=0;i<SPACEDIM;i++)
    {
      int tmp;
      stream >> tmp;
      _coords[i]=((double) tmp)/1e4;
    }
}

Node::~Node()
{
  if(_isToDel)
    delete [] _coords;
}

bool Node::decrRef()
{
  bool ret=(--_cnt==0);
  if(ret)
    delete this;
  return ret;
}

bool Node::isEqual(const Node& other) const
{
  const unsigned SPACEDIM=2;
  bool ret=true;
  for(unsigned i=0;i<SPACEDIM;i++)
    ret&=areDoubleEquals((*this)[i],other[i]);
  return ret;
}

double Node::getSlope(const Node& other) const
{
  double x=other[0]-(*this)[0];
  double y=other[1]-(*this)[1];
  double norm=sqrt(x*x+y*y);
  double ret=EdgeArcCircle::safeAcos(fabs(x)/norm);
  if( (x>=0. && y>=0.) || (x<0. && y<0.) )
    return ret;
  else
    return M_PI-ret;
}

/*!
 * Convenient method. Equivalent to isEqual method. In case of true is returned, '&other' is
 * added in 'track' container.
 */
bool Node::isEqualAndKeepTrack(const Node& other, std::vector<Node *>& track) const
{
  bool ret=isEqual(other);
  if(ret)
    track.push_back((Node *)&other);
  return ret;
}

void Node::dumpInXfigFile(std::ostream& stream, int resolution, const Bounds& box) const
{
  stream << box.fitXForXFig(_coords[0],resolution) << " " << box.fitYForXFig(_coords[1],resolution) << " ";
}

double Node::distanceWithSq(const Node& other) const
{
  return (_coords[0]-other._coords[0])*(_coords[0]-other._coords[0])+(_coords[1]-other._coords[1])*(_coords[1]-other._coords[1]);
}
