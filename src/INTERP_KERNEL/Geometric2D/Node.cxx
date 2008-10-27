#include "Node.hxx"
#include "EdgeArcCircle.hxx"

using namespace std;
using namespace INTERP_KERNEL;

Node::Node(double x, double y):_cnt(1),_loc(UNKNOWN)
{
  _coords[0]=x; _coords[1]=y;
}

Node::Node(const double *coords):_cnt(1),_loc(UNKNOWN)
{
  _coords[0]=coords[0];
  _coords[1]=coords[1];
}

Node::Node(std::istream& stream):_cnt(1),_loc(UNKNOWN)
{
  int tmp;
  stream >> tmp;
  _coords[0]=((double) tmp)/1e4;
  stream >> tmp;
  _coords[1]=((double) tmp)/1e4;
}

Node::~Node()
{
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
  const int SPACEDIM=2;
  bool ret=true;
  for(int i=0;i<SPACEDIM;i++)
    ret&=areDoubleEquals((*this)[i],other[i]);
  return ret;
}

double Node::getSlope(const Node& other) const
{
  return computeSlope(*this, other);
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

/*!
 * WARNING different from 'computeAngle' method ! The returned value are not in the same interval !
 * Here in -Pi/2; Pi/2. Typically this method returns the same value by exchanging pt1 and pt2.
 * Use in process of detection of a point in or not in polygon.
 */
double Node::computeSlope(const double *pt1, const double *pt2)
{
  double x=pt2[0]-pt1[0];
  double y=pt2[1]-pt1[1];
  double norm=sqrt(x*x+y*y);
  double ret=EdgeArcCircle::safeAcos(fabs(x)/norm);
  if( (x>=0. && y>=0.) || (x<0. && y<0.) )
    return ret;
  else
    return M_PI-ret;
}

/*!
 * WARNING different from 'computeSlope' method. Here angle in -Pi;Pi is returned.
 * This method is anti-symetric.
 */
double Node::computeAngle(const double *pt1, const double *pt2)
{
  double x=pt2[0]-pt1[0];
  double y=pt2[1]-pt1[1];
  double norm=sqrt(x*x+y*y);
  double ret=EdgeArcCircle::safeAcos(x/norm);
  if(y>=0)
    return ret;
  else
    return -ret;
}

/*!
 * apply a Similarity transformation on this.
 * @param xBary is the opposite of the X translation to do.
 * @param yBary is the opposite of the Y translation to do.
 * @param dimChar is the reduction factor.
 */
void Node::applySimilarity(double xBary, double yBary, double dimChar)
{
  _coords[0]=(_coords[0]-xBary)/dimChar;
  _coords[1]=(_coords[1]-yBary)/dimChar;
}
