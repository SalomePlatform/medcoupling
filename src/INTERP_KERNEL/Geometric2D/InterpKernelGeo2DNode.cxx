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

#include "InterpKernelGeo2DNode.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"

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
  const unsigned SPACEDIM=2;
  bool ret=true;
  for(unsigned i=0;i<SPACEDIM;i++)
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
    track.push_back(const_cast<Node *>(&other));
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
 * Here in [0; Pi). Typically this method returns the same value by exchanging pt1 and pt2.
 * Use in process of detection of a point in or not in polygon.
 */
double Node::computeSlope(const double *pt1, const double *pt2)
{
  double x=pt2[0]-pt1[0];
  double y=pt2[1]-pt1[1];
  double norm=sqrt(x*x+y*y);
  double ret=EdgeArcCircle::SafeAcos(fabs(x)/norm);
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
  return EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect(x/norm,y/norm);
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

/*!
 * apply the reverse Similarity transformation on this.
 * This method is the opposite of Node::applySimilarity method to retrieve the initial state.
 * @param xBary is the opposite of the X translation to do.
 * @param yBary is the opposite of the Y translation to do.
 * @param dimChar is the reduction factor.
 */
void Node::unApplySimilarity(double xBary, double yBary, double dimChar)
{
  _coords[0]=_coords[0]*dimChar+xBary;
  _coords[1]=_coords[1]*dimChar+yBary;
}

/*!
 * Called by QuadraticPolygon::splitAbs method.
 */
void Node::fillGlobalInfoAbs(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                             std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo, int *nodeId) const
{
  std::map<INTERP_KERNEL::Node *,int>::const_iterator it=mapThis.find(const_cast<Node *>(this));
  if(it!=mapThis.end())
    {
      *nodeId=(*it).second;
      return;
    }
  it=mapOther.find(const_cast<Node *>(this));
  if(it!=mapOther.end())
    {
      *nodeId=(*it).second+offset1;
      return;
    }
  it=mapAddCoo.find(const_cast<Node *>(this));
  if(it!=mapAddCoo.end())
    {
      *nodeId=(*it).second;
      return;
    }
  int id=(int)addCoo.size()/2;
  addCoo.push_back(fact*_coords[0]+baryX);
  addCoo.push_back(fact*_coords[1]+baryY);
  *nodeId=offset2+id;
  mapAddCoo[const_cast<Node *>(this)]=offset2+id;
}

/*!
 * Called by QuadraticPolygon::splitAbs method.
 */
void Node::fillGlobalInfoAbs2(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                              std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo, std::vector<int>& pointsOther) const
{
  int tmp;
  std::size_t sz1=addCoo.size();
  fillGlobalInfoAbs(mapThis,mapOther,offset1,offset2,fact,baryX,baryY,addCoo,mapAddCoo,&tmp);
  if(sz1!=addCoo.size())
    {
      pointsOther.push_back(tmp);
      return ;
    }
  std::vector<int>::const_iterator it=std::find(pointsOther.begin(),pointsOther.end(),tmp);
  if(it!=pointsOther.end())
    return ;
  if(tmp<offset1)
    pointsOther.push_back(tmp);
}
