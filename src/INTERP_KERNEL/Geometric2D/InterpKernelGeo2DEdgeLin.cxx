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

#include "InterpKernelGeo2DEdgeLin.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "InterpKernelException.hxx"
#include "NormalizedUnstructuredMesh.hxx"

using namespace INTERP_KERNEL;

namespace INTERP_KERNEL
{
  extern const unsigned MAX_SIZE_OF_LINE_XFIG_FILE=1024;
}

SegSegIntersector::SegSegIntersector(const EdgeLin& e1, const EdgeLin& e2):SameTypeEdgeIntersector(e1,e2)
{
  _matrix[0]=(*(e2.getStartNode()))[0]-(*(e2.getEndNode()))[0];
  _matrix[1]=(*(e1.getEndNode()))[0]-(*(e1.getStartNode()))[0];
  _matrix[2]=(*(e2.getStartNode()))[1]-(*(e2.getEndNode()))[1];
  _matrix[3]=(*(e1.getEndNode()))[1]-(*(e1.getStartNode()))[1];
  _col[0]=_matrix[3]*(*(e1.getStartNode()))[0]-_matrix[1]*(*(e1.getStartNode()))[1];
  _col[1]=-_matrix[2]*(*(e2.getStartNode()))[0]+_matrix[0]*(*(e2.getStartNode()))[1];
  //Little trick to avoid problems if 'e1' and 'e2' are colinears and along Ox or Oy axes.
  if(fabs(_matrix[3])>fabs(_matrix[1]))
    _ind=0;
  else
    _ind=1;
}

/*!
 * Must be called when 'this' and 'other' have been detected to be at least colinear. Typically they are overlapped.
 * Must be called after call of areOverlappedOrOnlyColinears.
 */
bool SegSegIntersector::haveTheySameDirection() const
{
  return (_matrix[3]*_matrix[1]+_matrix[2]*_matrix[0])>0.;
  //return (_matrix[_ind?1:0]>0. && _matrix[_ind?3:2]>0.) || (_matrix[_ind?1:0]<0. && _matrix[_ind?3:2]<0.);
}

/*!
 * Precondition start and end must be so that there predecessor was in the same direction than 'e1'
 */
void SegSegIntersector::getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const
{
  getCurveAbscisse(start,whereStart,commonNode);
  getCurveAbscisse(end,whereEnd,commonNode);
}

void SegSegIntersector::getCurveAbscisse(Node *node, TypeOfLocInEdge& where, MergePoints& commonNode) const
{
  bool obvious;
  obviousCaseForCurvAbscisse(node,where,commonNode,obvious);
  if(obvious)
    return ;
  double ret=((*node)[!_ind]-(*_e1.getStartNode())[!_ind])/((*_e1.getEndNode())[!_ind]-(*_e1.getStartNode())[!_ind]);
  if(ret>0. && ret <1.)
    where=INSIDE;
  else if(ret<0.)
    where=OUT_BEFORE;
  else
    where=OUT_AFTER;
}

/*!
 * areColinears method should be called before with a returned colinearity equal to false to avoid bad news.
 */
std::list< IntersectElement > SegSegIntersector::getIntersectionsCharacteristicVal() const
{
  std::list< IntersectElement > ret;
  double x=_matrix[0]*_col[0]+_matrix[1]*_col[1];
  double y=_matrix[2]*_col[0]+_matrix[3]*_col[1];
  //Only one intersect point possible
  Node *node=new Node(x,y);
  node->declareOn();
  bool i_1S=_e1.getStartNode()->isEqual(*node);
  bool i_1E=_e1.getEndNode()->isEqual(*node);
  bool i_2S=_e2.getStartNode()->isEqual(*node);
  bool i_2E=_e2.getEndNode()->isEqual(*node);
  ret.push_back(IntersectElement(_e1.getCharactValue(*node),
      _e2.getCharactValue(*node),
      i_1S,i_1E,i_2S,i_2E,node,_e1,_e2,keepOrder()));
  return ret;
}

/*!
 * retrieves if segs are colinears.
 * WARNING !!! Contrary to areOverlappedOrOnlyColinears method, this method use an
 * another precision to detect colinearity !
 */
bool SegSegIntersector::areColinears() const
{
  double determinant=_matrix[0]*_matrix[3]-_matrix[1]*_matrix[2];
  return fabs(determinant)<QuadraticPlanarArcDetectionPrecision::getArcDetectionPrecision();
}

/*!
 * Should be called \b once ! non const method.
 * \param whereToFind specifies the box where final seek should be done. Essentially it is used for caracteristic reason.
 * \param colinearity returns if regarding QuadraticPlanarPrecision::getPrecision() ; e1 and e2 are colinears
 *                    If true 'this' is modified ! So this method be called once above all if true is returned for this parameter.
 * \param areOverlapped if colinearity if true, this parameter looks if e1 and e2 are overlapped.
 */
void SegSegIntersector::areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& colinearity, bool& areOverlapped)
{
  double determinant=_matrix[0]*_matrix[3]-_matrix[1]*_matrix[2];
  if(fabs(determinant)>2.*QuadraticPlanarPrecision::getPrecision())//2*_precision due to max of offset on _start and _end
    {
      colinearity=false; areOverlapped=false;
      _matrix[0]/=determinant; _matrix[1]/=determinant; _matrix[2]/=determinant; _matrix[3]/=determinant;
    }
  else
    {
      colinearity=true;
      //retrieving initial matrix
      double tmp=_matrix[0]; _matrix[0]=_matrix[3]; _matrix[3]=tmp;
      _matrix[1]=-_matrix[1]; _matrix[2]=-_matrix[2];
      //
      double deno=sqrt(_matrix[0]*_matrix[0]+_matrix[1]*_matrix[1]);
      double x=(*(_e1.getStartNode()))[0]-(*(_e2.getStartNode()))[0];
      double y=(*(_e1.getStartNode()))[1]-(*(_e2.getStartNode()))[1];
      areOverlapped=fabs((_matrix[1]*y+_matrix[0]*x)/deno)<QuadraticPlanarPrecision::getPrecision();
    }
}

EdgeLin::EdgeLin(std::istream& lineInXfig)
{
  char currentLine[MAX_SIZE_OF_LINE_XFIG_FILE];
  lineInXfig.getline(currentLine,MAX_SIZE_OF_LINE_XFIG_FILE);
  _start=new Node(lineInXfig);
  _end=new Node(lineInXfig);
  updateBounds();
}

EdgeLin::EdgeLin(Node *start, Node *end, bool direction):Edge(start,end,direction)
{
  updateBounds();
}

EdgeLin::EdgeLin(double sX, double sY, double eX, double eY):Edge(sX,sY,eX,eY)
{
  updateBounds();
}

EdgeLin::~EdgeLin()
{
}

/*!
 * Characteristic for edges is relative position btw 0.;1.
 */
bool EdgeLin::isIn(double characterVal) const
{
  return characterVal>0. && characterVal<1.;
}

Node *EdgeLin::buildRepresentantOfMySelf() const
{
  return new Node(((*(_start))[0]+(*(_end))[0])/2.,((*(_start))[1]+(*(_end))[1])/2.);
}

double EdgeLin::getCharactValue(const Node& node) const
{
  return getCharactValueEng(node);
}

double EdgeLin::getCharactValueBtw0And1(const Node& node) const
{
  return getCharactValueEng(node);
}

double EdgeLin::getDistanceToPoint(const double *pt) const
{
  double loc=getCharactValueEng(pt);
  if(loc>0. && loc<1.)
    {
      double tmp[2];
      tmp[0]=(*_start)[0]*(1-loc)+loc*(*_end)[0];
      tmp[1]=(*_start)[1]*(1-loc)+loc*(*_end)[1];
      return Node::distanceBtw2Pt(pt,tmp);
    }
  else
    {
      double dist1=Node::distanceBtw2Pt(*_start,pt);
      double dist2=Node::distanceBtw2Pt(*_end,pt);
      return std::min(dist1,dist2);
    }
}

bool EdgeLin::isNodeLyingOn(const double *coordOfNode) const
{
  double dBase=sqrt(_start->distanceWithSq(*_end));
  double d1=Node::distanceBtw2Pt(*_start,coordOfNode);
  d1+=Node::distanceBtw2Pt(*_end,coordOfNode);
  return Node::areDoubleEquals(dBase,d1);
}

void EdgeLin::dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const
{
  stream << "2 1 0 1 ";
  fillXfigStreamForLoc(stream);
  stream << " 7 50 -1 -1 0.000 0 0 -1 1 0 2" << std::endl << "1 1 1.00 60.00 120.00" << std::endl;
  direction?_start->dumpInXfigFile(stream,resolution,box):_end->dumpInXfigFile(stream,resolution,box);
  direction?_end->dumpInXfigFile(stream,resolution,box):_start->dumpInXfigFile(stream,resolution,box);
  stream << std::endl;
}

void EdgeLin::update(Node *m)
{
  updateBounds();
}

double EdgeLin::getNormSq() const
{
  return _start->distanceWithSq(*_end);
}

/*!
 * This methods computes :
 * \f[
 * \int_{Current Edge} -ydx
 * \f]
 */
double EdgeLin::getAreaOfZone() const
{
  return ((*_start)[0]-(*_end)[0])*((*_start)[1]+(*_end)[1])/2.;
}

void EdgeLin::getBarycenter(double *bary) const
{
  bary[0]=((*_start)[0]+(*_end)[0])/2.;
  bary[1]=((*_start)[1]+(*_end)[1])/2.;
}

/*!
 * \f[
 * bary[0]=\int_{Current Edge} -yxdx
 * \f]
 * \f[
 * bary[1]=\int_{Current Edge} -\frac{y^{2}}{2}dx
 * \f]
 * To compute these 2 expressions in this class we have :
 * \f[
 * y=y_{1}+\frac{y_{2}-y_{1}}{x_{2}-x_{1}}(x-x_{1})
 * \f]
 */
void EdgeLin::getBarycenterOfZone(double *bary) const
{
  double x1=(*_start)[0];
  double y1=(*_start)[1];
  double x2=(*_end)[0];
  double y2=(*_end)[1];
  bary[0]=(x1-x2)*(y1*(2.*x1+x2)+y2*(2.*x2+x1))/6.;
  //bary[0]+=(y1-y2)*(x2*x2/3.-(x1*x2+x1*x1)/6.)+y1*(x1*x1-x2*x2)/2.;
  //bary[0]+=(y1-y2)*((x2*x2+x1*x2+x1*x1)/3.-(x2+x1)*x1/2.)+y1*(x1*x1-x2*x2)/2.;
  bary[1]=(x1-x2)*(y1*(y1+y2)+y2*y2)/6.;
}

/*!
 * Here \a this is not used (contrary to EdgeArcCircle class).
 */
void EdgeLin::getMiddleOfPoints(const double *p1, const double *p2, double *mid) const
{
  mid[0]=(p1[0]+p2[0])/2.;
  mid[1]=(p1[1]+p2[1])/2.;
}

double EdgeLin::getCurveLength() const
{
  double x=(*_start)[0]-(*_end)[0];
  double y=(*_start)[1]-(*_end)[1];
  return sqrt(x*x+y*y);
}

Edge *EdgeLin::buildEdgeLyingOnMe(Node *start, Node *end, bool direction) const
{
  return new EdgeLin(start,end,direction);
}

/*!
 * No precision should be introduced here. Just think as if precision was perfect.
 */
void EdgeLin::updateBounds()
{
  _bounds.setValues(std::min((*_start)[0],(*_end)[0]),std::max((*_start)[0],(*_end)[0]),std::min((*_start)[1],(*_end)[1]),std::max((*_start)[1],(*_end)[1]));
}

double EdgeLin::getCharactValueEng(const double *node) const
{
  double car1_1x=node[0]-(*(_start))[0]; double car1_2x=(*(_end))[0]-(*(_start))[0];
  double car1_1y=node[1]-(*(_start))[1]; double car1_2y=(*(_end))[1]-(*(_start))[1];
  return (car1_1x*car1_2x+car1_1y*car1_2y)/(car1_2x*car1_2x+car1_2y*car1_2y);
}

void EdgeLin::fillGlobalInfoAbs(bool direction, const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                                std::vector<int>& edgesThis, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int> mapAddCoo) const
{
  int tmp[2];
  _start->fillGlobalInfoAbs(mapThis,mapOther,offset1,offset2,fact,baryX,baryY,addCoo,mapAddCoo,tmp);
  _end->fillGlobalInfoAbs(mapThis,mapOther,offset1,offset2,fact,baryX,baryY,addCoo,mapAddCoo,tmp+1);
  if(direction)
    {
      edgesThis.push_back(tmp[0]);
      edgesThis.push_back(tmp[1]);
    }
  else
    {
      edgesThis.push_back(tmp[1]);
      edgesThis.push_back(tmp[0]);
    }
}

void EdgeLin::fillGlobalInfoAbs2(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                                 std::vector<int>& edgesOther, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo) const
{
  _start->fillGlobalInfoAbs2(mapThis,mapOther,offset1,offset2,fact,baryX,baryY,addCoo,mapAddCoo,edgesOther);
  _end->fillGlobalInfoAbs2(mapThis,mapOther,offset1,offset2,fact,baryX,baryY,addCoo,mapAddCoo,edgesOther);
}
