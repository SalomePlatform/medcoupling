#include "EdgeArcCircle.hxx"
#include "EdgeLin.hxx"
#include "InterpolationUtils.hxx"
#include "Node.hxx"

#include <sstream>

using namespace std;
using namespace INTERP_KERNEL;

ArcCArcCIntersector::ArcCArcCIntersector(const EdgeArcCircle& e1, const EdgeArcCircle& e2):SameTypeIntersector(e1,e2),_dist(0.)
{
}

bool ArcCArcCIntersector::haveTheySameDirection() const
{
  return (getE1().getAngle()>0. &&  getE2().getAngle()>0.) || (getE1().getAngle()<0. &&  getE2().getAngle()<0.);
}

/*!
 * Precondition 'start' and 'end' are on the same curve than this.
 */
void ArcCArcCIntersector::getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const
{
  bool obvious1,obvious2;
  obviousCaseForCurvAbscisse(start,whereStart,commonNode,obvious1);
  obviousCaseForCurvAbscisse(end,whereEnd,commonNode,obvious2);
  if(obvious1 && obvious2)
    return ;
  double angleInRadStart=getAngle(start);
  double angleInRadEnd=getAngle(end);
  if(obvious1 || obvious2)
    {
      if(obvious1)
        {
          if(isIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadEnd))
            whereEnd=INSIDE;
          else
            whereEnd=OUT_AFTER;
          return ;
        }
      else
        {
          if(isIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadStart))
            whereStart=INSIDE;
          else
            whereStart=OUT_BEFORE;
          return ;
        }
    }
  if(isIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadStart))
    {
      whereStart=INSIDE;
      if(isIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadEnd))
        whereEnd=INSIDE;
      else
        whereEnd=OUT_AFTER;
    }
  else
    {//we are out in start.
      if(isIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadEnd))
        {
          whereStart=OUT_BEFORE;
          whereEnd=INSIDE;
        }
      else
        {
          if(isIn2Pi(getE2().getAngle0(),getE2().getAngle(),getE1().getAngle0()))
            {//_e2 contains stictly _e1
              whereStart=OUT_BEFORE;
              whereEnd=OUT_AFTER;
            }
          else
            {//_e2 is outside from _e1
              whereStart=OUT_BEFORE;
              whereEnd=OUT_BEFORE;
            }
        }
    }
}

/*!
 * Return angle between ]-Pi;Pi[
 */
double ArcCArcCIntersector::getAngle(Node *node) const
{
  double ret=EdgeArcCircle::safeAcos(((*node)[0]-getE1().getCenter()[0])/getE1().getRadius());
  if(((*node)[1]<getE1().getCenter()[1]))
    ret=-ret;
  return ret;
}

void ArcCArcCIntersector::areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& obviousNoIntersection, bool& areOverlapped)
{
  _dist=Node::distanceBtw2Pt(getE1().getCenter(),getE2().getCenter());
  double radius1=getE1().getRadius(); double radius2=getE2().getRadius();
  if(_dist>radius1+radius2+QUADRATIC_PLANAR::_precision || _dist+fmin(radius1,radius2)+QUADRATIC_PLANAR::_precision<fmax(radius1,radius2))
    {
      obviousNoIntersection=true;
      areOverlapped=false;
      return ;
    }
  if(Node::areDoubleEquals(_dist,0.) && Node::areDoubleEquals(radius1,radius2))
    {
      obviousNoIntersection=false;
      areOverlapped=true;
    }
  else
    {
      obviousNoIntersection=false;
      areOverlapped=false;
    }
}

std::list< IntersectElement > ArcCArcCIntersector::getIntersectionsCharacteristicVal() const
{
  std::list< IntersectElement > ret;
  const double *center1=getE1().getCenter();
  const double *center2=getE2().getCenter();
  double radius1=getE1().getRadius(); double radius2=getE2().getRadius();
  double d1_1=(_dist*_dist-radius2*radius2+radius1*radius1)/(2.*_dist);
  double u[2];//u is normalized vector from center1 to center2.
  u[0]=(center2[0]-center1[0])/_dist; u[1]=(center2[1]-center1[1])/_dist;
  double angle0_1=EdgeArcCircle::safeAcos(u[0]);
  if(u[1]<0.)
    angle0_1=-angle0_1;
  double d1_1y=EdgeArcCircle::safeSqrt(radius1*radius1-d1_1*d1_1);
  double angleE1=normalizeAngle(getE1().getAngle0()+getE1().getAngle());
  double angleE2=normalizeAngle(getE2().getAngle0()+getE2().getAngle());
  if(!Node::areDoubleEquals(d1_1y,0))
    {
      //2 intersections
      double v1[2],v2[2];
      v1[0]=u[0]*d1_1-u[1]*d1_1y; v1[1]=u[1]*d1_1+u[0]*d1_1y;
      v2[0]=u[0]*d1_1+u[1]*d1_1y; v2[1]=u[1]*d1_1-u[0]*d1_1y;
      Node *node1=new Node(center1[0]+v1[0],center1[1]+v1[1]); node1->declareOn();
      Node *node2=new Node(center1[0]+v2[0],center1[1]+v2[1]); node2->declareOn();
      double angle1_1=EdgeArcCircle::safeAcos(v1[0]/radius1);
      if(v1[1]<0.)
	angle1_1=-angle1_1;
      double angle2_1=EdgeArcCircle::safeAcos(v2[0]/radius1);
      if(v2[1]<0.)
	angle2_1=-angle2_1;
      double v3[2],v4[2];
      v3[0]=center1[0]-center2[0]+v1[0]; v3[1]=center1[1]-center2[1]+v1[1];
      v4[0]=center1[0]-center2[0]+v2[0]; v4[1]=center1[1]-center2[1]+v2[1];
      double angle1_2=EdgeArcCircle::safeAcos(v3[0]/radius2);
      if(v3[1]<0.)
	angle1_2=-angle1_2;
      double angle2_2=EdgeArcCircle::safeAcos(v4[0]/radius2);
      if(v4[1]<0.)
	angle2_2=-angle2_2;
      //
      bool e1_1S=Node::areDoubleEqualsWP(angle1_1,getE1().getAngle0(),radius1);
      bool e1_1E=Node::areDoubleEqualsWP(angle1_1,angleE1,radius1);
      bool e1_2S=Node::areDoubleEqualsWP(angle1_2,getE2().getAngle0(),radius1);
      bool e1_2E=Node::areDoubleEqualsWP(angle1_2,angleE2,radius1);
      //
      bool e2_1S=Node::areDoubleEqualsWP(angle2_1,getE1().getAngle0(),radius2);
      bool e2_1E=Node::areDoubleEqualsWP(angle2_1,angleE1,radius2);
      bool e2_2S=Node::areDoubleEqualsWP(angle2_2,getE2().getAngle0(),radius2);
      bool e2_2E=Node::areDoubleEqualsWP(angle2_2,angleE2,radius2);
      ret.push_back(IntersectElement(angle1_1,angle1_2,e1_1S,e1_1E,e1_2S,e1_2E,node1,_e1,_e2,keepOrder()));
      ret.push_back(IntersectElement(angle2_1,angle2_2,e2_1S,e2_1E,e2_2S,e2_2E,node2,_e1,_e2,keepOrder()));
    }
  else
    {
      //tangent intersection
      double v1[2],v2[2];
      v1[0]=d1_1*u[0]; v1[1]=d1_1*u[1];
      v2[0]=center1[0]-center2[0]+v1[0]; v2[1]=center1[1]-center2[1]+v1[1];
      double angle0_1=EdgeArcCircle::safeAcos(v1[0]/radius1);
      if(v1[1]<0.)
	angle0_1=-angle0_1;
      double angle0_2=EdgeArcCircle::safeAcos(v2[0]/radius2);
      if(v2[1]<0.)
	angle0_2=-angle0_2;
      bool e0_1S=Node::areDoubleEqualsWP(angle0_1,getE1().getAngle0(),radius1);
      bool e0_1E=Node::areDoubleEqualsWP(angle0_1,angleE1,radius1);
      bool e0_2S=Node::areDoubleEqualsWP(angle0_2,getE2().getAngle0(),radius2);
      bool e0_2E=Node::areDoubleEqualsWP(angle0_2,angleE2,radius2);
      Node *node=new Node(center1[0]+d1_1*u[0],center1[1]+d1_1*u[1]); node->declareOnTangent();
      ret.push_back(IntersectElement(angle0_1,angle0_2,e0_1S,e0_1E,e0_2S,e0_2E,node,_e1,_e2,keepOrder()));
    }
  return ret;
}
  /*double angle0_2;
  double signDeltaAngle2;
  double d1_2;
  if(u[1]<0.)
    angle0_1=-angle0_1;
  if(d1_1>=0.)
    {
      if(_dist>radius1)
        {
          angle0_2=angle0_1+M_PI;
          signDeltaAngle2=-1.;
        }
      else
        {
          angle0_2=angle0_1;
          signDeltaAngle2=1.;
        }
    }
  else
    {
      angle0_1+=M_PI;
      angle0_2=angle0_1;
      signDeltaAngle2=1.;
    }
  angle0_1=normalizeAngle(angle0_1);
  angle0_2=normalizeAngle(angle0_2);
  double angleE1=normalizeAngle(getE1().getAngle0()+getE1().getAngle());
  double angleE2=normalizeAngle(getE2().getAngle0()+getE2().getAngle());
  if(!(Node::areDoubleEquals(d1_1,radius1) || Node::areDoubleEquals(d1_1,-radius1)) )
    {
      //2 intersections   
      double deltaAngle1=EdgeArcCircle::safeAcos(fabs(d1_1)/radius1); //owns to 0;Pi/2 by construction
      double deltaAngle2=EdgeArcCircle::safeAcos(fabs(d1_2)/radius2); //owns to 0;Pi/2 by construction
      double angle1_1=normalizeAngle(angle0_1+deltaAngle1);// Intersection 1 seen for _e1
      double angle2_1=normalizeAngle(angle0_1-deltaAngle1);// Intersection 2 seen for _e1
      double angle1_2=normalizeAngle(angle0_2+signDeltaAngle2*deltaAngle2);// Intersection 1 seen for _e2
      double angle2_2=normalizeAngle(angle0_2-signDeltaAngle2*deltaAngle2);// Intersection 2 seen for _e2
      //
      bool e1_1S=Node::areDoubleEqualsWP(angle1_1,getE1().getAngle0(),radius1);
      bool e1_1E=Node::areDoubleEqualsWP(angle1_1,angleE1,radius1);
      bool e1_2S=Node::areDoubleEqualsWP(angle1_2,getE2().getAngle0(),radius1);
      bool e1_2E=Node::areDoubleEqualsWP(angle1_2,angleE2,radius1);
      //
      bool e2_1S=Node::areDoubleEqualsWP(angle2_1,getE1().getAngle0(),radius2);
      bool e2_1E=Node::areDoubleEqualsWP(angle2_1,angleE1,radius2);
      bool e2_2S=Node::areDoubleEqualsWP(angle2_2,getE2().getAngle0(),radius2);
      bool e2_2E=Node::areDoubleEqualsWP(angle2_2,angleE2,radius2);
      Node *node1=new Node(center1[0]+radius1*cos(angle1_1),center1[0]+radius1*sin(angle1_1)); node1->declareOn();
      Node *node2=new Node(center1[0]+radius1*cos(angle2_1),center1[0]+radius1*sin(angle2_1)); node2->declareOn();
      ret.push_back(IntersectElement(angle1_1,angle1_2,e1_1S,e1_1E,e1_2S,e1_2E,node1,_e1,_e2,keepOrder()));
      ret.push_back(IntersectElement(angle2_1,angle2_2,e2_1S,e2_1E,e2_2S,e2_2E,node2,_e1,_e2,keepOrder()));
    }
  else
    //tangent intersection
    {
      bool e0_1S=Node::areDoubleEqualsWP(angle0_1,getE1().getAngle0(),radius1);
      bool e0_1E=Node::areDoubleEqualsWP(angle0_1,angleE1,radius1);
      bool e0_2S=Node::areDoubleEqualsWP(angle0_2,getE2().getAngle0(),radius2);
      bool e0_2E=Node::areDoubleEqualsWP(angle0_2,angleE2,radius2);
      Node *node=new Node(center1[0]+radius1*cos(angle0_1),center1[0]+radius1*sin(angle0_1)); node->declareOnTangent();
      ret.push_back(IntersectElement(angle0_1,angle0_2,e0_1S,e0_1E,e0_2S,e0_2E,node,_e1,_e2,keepOrder()));
    }
    return ret;*/

/*!
 * Idem isAngleNotIn except that here 'start' in ]-Pi;Pi[ and delta in ]-2*Pi;2Pi[. 
 * @param angleIn in ]-Pi;Pi[.
 */
bool ArcCArcCIntersector::isIn2Pi(double start, double delta, double angleIn)
{
  double myDelta=angleIn-start;
  if(delta>0.)
    {
      myDelta=myDelta>=0.?myDelta:myDelta+2.*M_PI;
      return myDelta>0. && myDelta<delta;
    }
  else
    {
      myDelta=myDelta<=0.?myDelta:myDelta-2.*M_PI;
      return myDelta<0. && myDelta>delta;
    }
}

/*!
 * Given the arc 'a' defined by 'start' angle and a 'delta' [-Pi;Pi] states for the angle 'angleIn' [-Pi;Pi] if it owns or not 'a'.
 */
bool ArcCArcCIntersector::isAngleNotIn(double start, double delta, double angleIn)
{
  double tmp=start;
  if(tmp<0.)
    tmp+=2*M_PI;
  double tmp2=angleIn;
  if(tmp2<0.)
    tmp2+=2*M_PI;
  if(tmp+delta>=2.*M_PI)
    return (tmp2<tmp) && (tmp2>tmp+delta-2*M_PI);
  else if(tmp+delta>=0.)
    return (tmp2<fmin(tmp,tmp+delta) || tmp2>fmax(tmp,tmp+delta));
  else
    return (tmp2>tmp) && (tmp2<(tmp+delta+2.*M_PI));
}

ArcCSegIntersector::ArcCSegIntersector(const EdgeArcCircle& e1, const EdgeLin& e2, bool reverse):CrossTypeIntersector(e1,e2,reverse)
{
}

void ArcCSegIntersector::areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& obviousNoIntersection, bool& areOverlapped)
{
  areOverlapped=false;//No overlapping by contruction
  const double *center=getE1().getCenter();
  _dx=(*(_e2.getEndNode()))[0]-(*(_e2.getStartNode()))[0];
  _dy=(*(_e2.getEndNode()))[1]-(*(_e2.getStartNode()))[1];
  _drSq=_dx*_dx+_dy*_dy;
  _cross=
    ((*(_e2.getStartNode()))[0]-center[0])*((*(_e2.getEndNode()))[1]-center[1])-
    ((*(_e2.getStartNode()))[1]-center[1])*((*(_e2.getEndNode()))[0]-center[0]);
  _determinant=getE1().getRadius()*getE1().getRadius()*_drSq-_cross*_cross;
  if(_determinant>-QUADRATIC_PLANAR::_precision*10.)//QUADRATIC_PLANAR::_precision*QUADRATIC_PLANAR::_precision*_drSq*_drSq/(2.*_dx*_dx))
    obviousNoIntersection=false;
  else
    obviousNoIntersection=true;   
}

void ArcCSegIntersector::getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const
{
  throw Exception("Internal error. Should never been called : no overlapping possible between arc of circle and a segment.");
}

std::list< IntersectElement > ArcCSegIntersector::getIntersectionsCharacteristicVal() const
{
  std::list< IntersectElement > ret;
  const double *center=getE1().getCenter();
  if(!(fabs(_determinant)<(QUADRATIC_PLANAR::_precision*10.)))//QUADRATIC_PLANAR::_precision*QUADRATIC_PLANAR::_precision*_drSq*_drSq/(2.*_dx*_dx))
    {
      double determinant=EdgeArcCircle::safeSqrt(_determinant);
      double x1=(_cross*_dy+Node::sign(_dy)*_dx*determinant)/_drSq+center[0];
      double y1=(-_cross*_dx+fabs(_dy)*determinant)/_drSq+center[1];
      Node *intersect1=new Node(x1,y1); intersect1->declareOn();
      bool i1_1S=_e1.getStartNode()->isEqual(*intersect1);
      bool i1_1E=_e1.getEndNode()->isEqual(*intersect1);
      bool i1_2S=_e2.getStartNode()->isEqual(*intersect1);
      bool i1_2E=_e2.getEndNode()->isEqual(*intersect1);
      ret.push_back(IntersectElement(getE1().getCharactValue(*intersect1),getE2().getCharactValue(*intersect1),i1_1S,i1_1E,i1_2S,i1_2E,intersect1,_e1,_e2,keepOrder()));
      //
      double x2=(_cross*_dy-Node::sign(_dy)*_dx*determinant)/_drSq+center[0];
      double y2=(-_cross*_dx-fabs(_dy)*determinant)/_drSq+center[1];
      Node *intersect2=new Node(x2,y2); intersect2->declareOn();
      bool i2_1S=_e1.getStartNode()->isEqual(*intersect2);
      bool i2_1E=_e1.getEndNode()->isEqual(*intersect2);
      bool i2_2S=_e2.getStartNode()->isEqual(*intersect2);
      bool i2_2E=_e2.getEndNode()->isEqual(*intersect2);
      ret.push_back(IntersectElement(getE1().getCharactValue(*intersect2),getE2().getCharactValue(*intersect2),i2_1S,i2_1E,i2_2S,i2_2E,intersect2,_e1,_e2,keepOrder()));
    }
  else//tangent intersection
    {
      double x=(_cross*_dy)/_drSq+center[0];
      double y=(-_cross*_dx)/_drSq+center[1];
      Node *intersect=new Node(x,y); intersect->declareOnTangent();
      bool i_1S=_e1.getStartNode()->isEqual(*intersect);
      bool i_1E=_e1.getEndNode()->isEqual(*intersect);
      bool i_2S=_e2.getStartNode()->isEqual(*intersect);
      bool i_2E=_e2.getEndNode()->isEqual(*intersect);
      ret.push_back(IntersectElement(_e1.getCharactValue(*intersect),_e2.getCharactValue(*intersect),i_1S,i_1E,i_2S,i_2E,intersect,_e1,_e2,keepOrder()));
    }
  return ret;
}

EdgeArcCircle::EdgeArcCircle(std::istream& lineInXfig)
{
  const unsigned NB_OF_SKIP_FIELDS=15;
  std::string tmpS;
  for(unsigned i=0;i<NB_OF_SKIP_FIELDS;i++)
    lineInXfig >> tmpS;
  _start=new Node(lineInXfig);
  Node *middle=new Node(lineInXfig);
  _end=new Node(lineInXfig);
  getArcOfCirclePassingThru(*_start,*middle,*_end,_center,_radius,_angle,_angle0);
  middle->decrRef();
  updateBounds();
}

EdgeArcCircle::EdgeArcCircle(Node *start, Node *middle, Node *end, bool direction):Edge(start,end, direction)
{
  getArcOfCirclePassingThru(*_start,*middle,*_end,_center,_radius,_angle,_angle0);
  updateBounds();
}

EdgeArcCircle::EdgeArcCircle(double sX, double sY, double mX, double mY, double eX, double eY):Edge(sX,sY,eX,eY)
{
  double middle[2]; middle[0]=mX; middle[1]=mY;
  getArcOfCirclePassingThru(*_start,middle,*_end,_center,_radius,_angle,_angle0);
  updateBounds();
}

/*!
 * @param angle0 in ]-Pi;Pi[
 * @param deltaAngle in ]-2.*Pi;2.*Pi[
 */
EdgeArcCircle::EdgeArcCircle(Node *start, Node *end, const double *center, double radius, double angle0, double deltaAngle, bool direction):Edge(start,end,direction),_angle(deltaAngle),
                                                                                                                                      _angle0(angle0),_radius(radius)
{
  _center[0]=center[0];
  _center[1]=center[1];
  updateBounds();
}

void EdgeArcCircle::changeMiddle(Node *newMiddle)
{
  getArcOfCirclePassingThru(*_start,*newMiddle,*_end,_center,_radius,_angle,_angle0);
  updateBounds();
}

Edge *EdgeArcCircle::buildEdgeLyingOnMe(Node *start, Node *end, bool direction) const
{
  double sx=((*start)[0]-_center[0])/_radius;
  double sy=((*start)[1]-_center[1])/_radius;
  double ex=((*end)[0]-_center[0])/_radius;
  double ey=((*end)[1]-_center[1])/_radius;
  double angle0=safeAcos(direction?sx:ex);
  angle0=(direction?sy:ey)>0.?angle0:-angle0;
  double deltaAngle=safeAcos(sx*ex+sy*ey);
  deltaAngle=sx*ey-sy*ex>0.?deltaAngle:-deltaAngle;
  if(deltaAngle>0. && _angle<0.)
    deltaAngle-=2.*M_PI;
  else if(deltaAngle<0. && _angle>0.)
    deltaAngle+=2.*M_PI;
  deltaAngle=direction?deltaAngle:-deltaAngle;
  return new EdgeArcCircle(start,end,_center,_radius,angle0,deltaAngle,direction);
}

void EdgeArcCircle::getArcOfCirclePassingThru(const double *start, const double *middle, const double *end, 
                                              double *center, double& radius, double& angleInRad, double& angleInRad0)
{
  double delta=(middle[0]-start[0])*(end[1]-middle[1])-(end[0]-middle[0])*(middle[1]-start[1]);
  double b1=(middle[1]*middle[1]+middle[0]*middle[0]-start[0]*start[0]-start[1]*start[1])/2;
  double b2=(end[1]*end[1]+end[0]*end[0]-middle[0]*middle[0]-middle[1]*middle[1])/2;
  center[0]=((end[1]-middle[1])*b1+(start[1]-middle[1])*b2)/delta;
  center[1]=((middle[0]-end[0])*b1+(middle[0]-start[0])*b2)/delta;
  radius=safeSqrt((start[0]-center[0])*(start[0]-center[0])+(start[1]-center[1])*(start[1]-center[1]));
  angleInRad0=safeAcos((start[0]-center[0])/radius);
  if((start[1]-center[1])<0.)
    angleInRad0=-angleInRad0;
  double angleInRadM=safeAcos((middle[0]-center[0])/radius);
  if((middle[1]-center[1])<0.)
    angleInRadM=-angleInRadM;
  angleInRad=safeAcos(((start[0]-center[0])*(end[0]-center[0])+(start[1]-center[1])*(end[1]-center[1]))/(radius*radius));
  bool signOfAngle=((start[0]-center[0])*(end[1]-center[1])-(start[1]-center[1])*(end[0]-center[0]))>0.;
  angleInRad=signOfAngle?angleInRad:-angleInRad;
  if(ArcCArcCIntersector::isAngleNotIn(angleInRad0,angleInRad,angleInRadM))
    angleInRad=angleInRad<0?2*M_PI+angleInRad:angleInRad-2*M_PI;
}

void EdgeArcCircle::dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const
{
  stream << "5 1 0 1 ";
  fillXfigStreamForLoc(stream);
  stream << " 7 50 -1 -1 0.000 0 ";
  if( (direction && _angle>=0) || (!direction && _angle<0))
    stream << '0';//'0'
  else
    stream << '1';//'1'
  stream << " 0 0 ";
  stream << box.fitXForXFigD(_center[0],resolution) << " " << box.fitYForXFigD(_center[1],resolution) << " ";
  direction?_start->dumpInXfigFile(stream,resolution,box):_end->dumpInXfigFile(stream,resolution,box);
  Node *middle=buildRepresentantOfMySelf();
  middle->dumpInXfigFile(stream,resolution,box);
  middle->decrRef();
  direction?_end->dumpInXfigFile(stream,resolution,box):_start->dumpInXfigFile(stream,resolution,box);
  stream << endl;
}

void EdgeArcCircle::update(Node *m)
{
  getArcOfCirclePassingThru(*_start,*m,*_end,_center,_radius,_angle,_angle0);
  updateBounds();
}

double EdgeArcCircle::getAreaOfZone() const
{
  return -_radius*_radius*(sin(_angle)-_angle)/2.+((*_start)[0]-(*_end)[0])*((*_start)[1]+(*_end)[1])/2.;
}

double EdgeArcCircle::getCurveLength() const
{
  return fabs(_angle*_radius);
}

void EdgeArcCircle::getBarycenter(double *bary) const
{
  bary[0]=((*_start)[0]+(*_end)[0])/2.;
  bary[1]=((*_start)[1]+(*_end)[1])/2.;
}

/*!
 * Characteristic value used is angle in ]_Pi;Pi[ from axe 0x.
 */
bool EdgeArcCircle::isIn(double characterVal) const
{
  return ArcCArcCIntersector::isIn2Pi(_angle0,_angle,characterVal);
}

Node *EdgeArcCircle::buildRepresentantOfMySelf() const
{
  return new Node(_center[0]+_radius*cos(_angle0+_angle/2.),_center[1]+_radius*sin(_angle0+_angle/2.));
}

/*!
 * Characteristic value used is angle in ]_Pi;Pi[ from axe 0x.
 * 'val1' and 'val2' have been detected previously as owning to this.
 */
bool EdgeArcCircle::isLower(double val1, double val2) const
{
  double myDelta1=val1-_angle0;
  double myDelta2=val2-_angle0;
  if(_angle>0.)
    {
      myDelta1=myDelta1>-(_radius*QUADRATIC_PLANAR::_precision)?myDelta1:myDelta1+2.*M_PI;//in some cases val1 or val2 are so close to angle0 that myDelta is close to 0. but negative.
      myDelta2=myDelta2>-(_radius*QUADRATIC_PLANAR::_precision)?myDelta2:myDelta2+2.*M_PI;
      return myDelta1<myDelta2;
    }
  else
    {
      myDelta1=myDelta1<(_radius*QUADRATIC_PLANAR::_precision)?myDelta1:myDelta1-2.*M_PI;
      myDelta2=myDelta2<(_radius*QUADRATIC_PLANAR::_precision)?myDelta2:myDelta2-2.*M_PI;
      return myDelta2<myDelta1;
    }
}

/*!
 * For Arc circle the caract value is angle with Ox between -Pi and Pi.
 */
double EdgeArcCircle::getCharactValue(const Node& node) const
{
  double dx=(node[0]-_center[0])/_radius;
  double dy=(node[1]-_center[1])/_radius;
  double angle0=safeAcos(dx);
  angle0=dy>=0.?angle0:-angle0;
  return angle0;
}

void EdgeArcCircle::updateBounds()
{
  _bounds.setValues(fmin((*_start)[0],(*_end)[0]),fmax((*_start)[0],(*_end)[0]),fmin((*_start)[1],(*_end)[1]),fmax((*_start)[1],(*_end)[1]));
  if(ArcCArcCIntersector::isIn2Pi(_angle0,_angle,M_PI/2))
    _bounds[3]=_center[1]+_radius;
  if(ArcCArcCIntersector::isIn2Pi(_angle0,_angle,-M_PI/2))
    _bounds[2]=_center[1]-_radius;
  if(ArcCArcCIntersector::isIn2Pi(_angle0,_angle,0.))
    _bounds[1]=_center[0]+_radius;
  if(ArcCArcCIntersector::isIn2Pi(_angle0,_angle,M_PI))
  _bounds[0]=_center[0]-_radius;
}
