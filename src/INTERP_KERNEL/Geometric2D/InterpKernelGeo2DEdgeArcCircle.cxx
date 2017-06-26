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

#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"
#include "InterpKernelException.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "NormalizedUnstructuredMesh.hxx"

#include <sstream>
#include <algorithm>

using namespace INTERP_KERNEL;

ArcCArcCIntersector::ArcCArcCIntersector(const EdgeArcCircle& e1, const EdgeArcCircle& e2):SameTypeEdgeIntersector(e1,e2),_dist(0.)
{
}

bool ArcCArcCIntersector::haveTheySameDirection() const
{
  return (getE1().getAngle()>0. &&  getE2().getAngle()>0.) || (getE1().getAngle()<0. &&  getE2().getAngle()<0.);
}

bool ArcCArcCIntersector::areColinears() const
{
  double radiusL,radiusB;
  double centerL[2],centerB[2];
  double tmp,cst;
  return internalAreColinears(getE1(),getE2(),tmp,cst,radiusL,centerL,radiusB,centerB);
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
          if(EdgeArcCircle::IsIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadEnd))
            whereEnd=INSIDE;
          else
            whereEnd=OUT_AFTER;
          return ;
        }
      else
        {
          if(EdgeArcCircle::IsIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadStart))
            whereStart=INSIDE;
          else
            whereStart=OUT_BEFORE;
          return ;
        }
    }
  if(EdgeArcCircle::IsIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadStart))
    {
      whereStart=INSIDE;
      if(EdgeArcCircle::IsIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadEnd))
        whereEnd=INSIDE;
      else
        whereEnd=OUT_AFTER;
    }
  else
    {//we are out in start.
      if(EdgeArcCircle::IsIn2Pi(getE1().getAngle0(),getE1().getAngle(),angleInRadEnd))
        {
          whereStart=OUT_BEFORE;
          whereEnd=INSIDE;
        }
      else
        {
          if(EdgeArcCircle::IsIn2Pi(getE2().getAngle0(),getE2().getAngle(),getE1().getAngle0()))
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
  return EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect(((*node)[0]-getE1().getCenter()[0])/getE1().getRadius(),((*node)[1]-getE1().getCenter()[1])/getE1().getRadius());
}

bool ArcCArcCIntersector::internalAreColinears(const EdgeArcCircle& a1, const EdgeArcCircle& a2, double& distBetweenCenters, double& cst,
                                               double& radiusL, double centerL[2], double& radiusB, double centerB[2])
{
  double lgth1=fabs(a1.getAngle()*a1.getRadius());
  double lgth2=fabs(a2.getAngle()*a2.getRadius());
  if(lgth1<lgth2)
    {//a1 is the little one ('L') and a2 the big one ('B')
      a1.getCenter(centerL); radiusL=a1.getRadius();
      a2.getCenter(centerB); radiusB=a2.getRadius();
    }
  else
    {
      a2.getCenter(centerL); radiusL=a2.getRadius();
      a1.getCenter(centerB); radiusB=a1.getRadius();
    }
  // dividing from the begining by radiusB^2 to keep precision
  distBetweenCenters=Node::distanceBtw2PtSq(centerL,centerB);
  cst=distBetweenCenters/(radiusB*radiusB);
  cst+=radiusL*radiusL/(radiusB*radiusB);
  return Node::areDoubleEqualsWP(cst,1.,2.);
}

bool ArcCArcCIntersector::areArcsOverlapped(const EdgeArcCircle& a1, const EdgeArcCircle& a2)
{
  double radiusL,radiusB;
  double centerL[2],centerB[2];
  double tmp(0.),cst(0.);
  if(!internalAreColinears(a1,a2,tmp,cst,radiusL,centerL,radiusB,centerB))
    return false;
  //
  double angle0L,angleL;
  Bounds *merge=a1.getBounds().nearlyAmIIntersectingWith(a2.getBounds());
  merge->getInterceptedArc(centerL,radiusL,angle0L,angleL);
  delete merge;
  //
  tmp=sqrt(tmp);
  if(Node::areDoubleEqualsWP(tmp,0.,1/(10*std::max(radiusL,radiusB))))
    return Node::areDoubleEquals(radiusL,radiusB);
  double phi=EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect((centerL[0]-centerB[0])/tmp,(centerL[1]-centerB[1])/tmp);
  double cst2=2*radiusL*tmp/(radiusB*radiusB);
  double cmpContainer[4];
  int sizeOfCmpContainer=2;
  cmpContainer[0]=cst+cst2*cos(phi-angle0L);
  cmpContainer[1]=cst+cst2*cos(phi-angle0L+angleL);
  double a=EdgeArcCircle::NormalizeAngle(phi-angle0L);
  if(EdgeArcCircle::IsIn2Pi(angle0L,angleL,a))
    cmpContainer[sizeOfCmpContainer++]=cst+cst2;
  a=EdgeArcCircle::NormalizeAngle(phi-angle0L+M_PI);
  if(EdgeArcCircle::IsIn2Pi(angle0L,angleL,a))
    cmpContainer[sizeOfCmpContainer++]=cst-cst2;
  a=*std::max_element(cmpContainer,cmpContainer+sizeOfCmpContainer);
  return Node::areDoubleEqualsWP(a,1.,2.);
}

void ArcCArcCIntersector::areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& obviousNoIntersection, bool& areOverlapped)
{
  _dist=Node::distanceBtw2Pt(getE1().getCenter(),getE2().getCenter());
  double radius1=getE1().getRadius(); double radius2=getE2().getRadius();
  if(_dist>radius1+radius2+QuadraticPlanarPrecision::getPrecision() || _dist+std::min(radius1,radius2)+QuadraticPlanarPrecision::getPrecision()<std::max(radius1,radius2))
    {
      obviousNoIntersection=true;
      areOverlapped=false;
      return ;
    }
  if(areArcsOverlapped(getE1(),getE2()))//(Node::areDoubleEquals(_dist,0.) && Node::areDoubleEquals(radius1,radius2))
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
  double d1_1y=EdgeArcCircle::SafeSqrt(radius1*radius1-d1_1*d1_1);
  double angleE1=EdgeArcCircle::NormalizeAngle(getE1().getAngle0()+getE1().getAngle());
  double angleE2=EdgeArcCircle::NormalizeAngle(getE2().getAngle0()+getE2().getAngle());
  if(!Node::areDoubleEquals(d1_1y,0))
    {
      //2 intersections
      double v1[2],v2[2];
      v1[0]=u[0]*d1_1-u[1]*d1_1y; v1[1]=u[1]*d1_1+u[0]*d1_1y;
      v2[0]=u[0]*d1_1+u[1]*d1_1y; v2[1]=u[1]*d1_1-u[0]*d1_1y;
      Node *node1=new Node(center1[0]+v1[0],center1[1]+v1[1]); node1->declareOn();
      Node *node2=new Node(center1[0]+v2[0],center1[1]+v2[1]); node2->declareOn();
      double angle1_1=EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect(v1[0]/radius1,v1[1]/radius1);
      double angle2_1=EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect(v2[0]/radius1,v2[1]/radius1);
      double v3[2],v4[2];
      v3[0]=center1[0]-center2[0]+v1[0]; v3[1]=center1[1]-center2[1]+v1[1];
      v4[0]=center1[0]-center2[0]+v2[0]; v4[1]=center1[1]-center2[1]+v2[1];
      double angle1_2=EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect(v3[0]/radius2,v3[1]/radius2);
      double angle2_2=EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect(v4[0]/radius2,v4[1]/radius2);
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
      double angle0_1=EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect(v1[0]/radius1,v1[1]/radius1);
      double angle0_2=EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect(v2[0]/radius2,v2[1]/radius2);
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
  angle0_1=NormalizeAngle(angle0_1);
  angle0_2=NormalizeAngle(angle0_2);
  double angleE1=NormalizeAngle(getE1().getAngle0()+getE1().getAngle());
  double angleE2=NormalizeAngle(getE2().getAngle0()+getE2().getAngle());
  if(!(Node::areDoubleEquals(d1_1,radius1) || Node::areDoubleEquals(d1_1,-radius1)) )
  {
  //2 intersections   
  double deltaAngle1=EdgeArcCircle::SafeAcos(fabs(d1_1)/radius1); //owns to 0;Pi/2 by construction
  double deltaAngle2=EdgeArcCircle::SafeAcos(fabs(d1_2)/radius2); //owns to 0;Pi/2 by construction
  double angle1_1=NormalizeAngle(angle0_1+deltaAngle1);// Intersection 1 seen for _e1
  double angle2_1=NormalizeAngle(angle0_1-deltaAngle1);// Intersection 2 seen for _e1
  double angle1_2=NormalizeAngle(angle0_2+signDeltaAngle2*deltaAngle2);// Intersection 1 seen for _e2
  double angle2_2=NormalizeAngle(angle0_2-signDeltaAngle2*deltaAngle2);// Intersection 2 seen for _e2
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

ArcCSegIntersector::ArcCSegIntersector(const EdgeArcCircle& e1, const EdgeLin& e2, bool reverse):CrossTypeEdgeIntersector(e1,e2,reverse)
{
}

void ArcCSegIntersector::areOverlappedOrOnlyColinears(const Bounds *whereToFind, bool& obviousNoIntersection, bool& areOverlapped)
{
  areOverlapped=false;//No overlapping by construction
  const double *center=getE1().getCenter();
  _dx=(*(_e2.getEndNode()))[0]-(*(_e2.getStartNode()))[0];
  _dy=(*(_e2.getEndNode()))[1]-(*(_e2.getStartNode()))[1];
  _drSq=_dx*_dx+_dy*_dy;
  _cross=
      ((*(_e2.getStartNode()))[0]-center[0])*((*(_e2.getEndNode()))[1]-center[1])-
      ((*(_e2.getStartNode()))[1]-center[1])*((*(_e2.getEndNode()))[0]-center[0]);
  _determinant=getE1().getRadius()*getE1().getRadius()/_drSq-_cross*_cross/(_drSq*_drSq);
  if(_determinant>-2*QuadraticPlanarPrecision::getPrecision())//QuadraticPlanarPrecision::getPrecision()*QuadraticPlanarPrecision::getPrecision()*_drSq*_drSq/(2.*_dx*_dx))
    obviousNoIntersection=false;
  else
    obviousNoIntersection=true;   
}

/*!
 * By construction, no chance that an arc of circle and line to be colinear.
 */
bool ArcCSegIntersector::areColinears() const
{
  return false;
}

void ArcCSegIntersector::getPlacements(Node *start, Node *end, TypeOfLocInEdge& whereStart, TypeOfLocInEdge& whereEnd, MergePoints& commonNode) const
{
  throw Exception("Internal error. Should never been called : no overlapping possible between arc of circle and a segment.");
}

std::list< IntersectElement > ArcCSegIntersector::getIntersectionsCharacteristicVal() const
{
  std::list< IntersectElement > ret;
  const double *center=getE1().getCenter();
  if(!(fabs(_determinant)<(2.*QuadraticPlanarPrecision::getPrecision())))//QuadraticPlanarPrecision::getPrecision()*QuadraticPlanarPrecision::getPrecision()*_drSq*_drSq/(2.*_dx*_dx))
    {
      double determinant=EdgeArcCircle::SafeSqrt(_determinant);
      double x1=(_cross*_dy/_drSq+Node::sign(_dy)*_dx*determinant)+center[0];
      double y1=(-_cross*_dx/_drSq+fabs(_dy)*determinant)+center[1];
      Node *intersect1=new Node(x1,y1); intersect1->declareOn();
      bool i1_1S=_e1.getStartNode()->isEqual(*intersect1);
      bool i1_1E=_e1.getEndNode()->isEqual(*intersect1);
      bool i1_2S=_e2.getStartNode()->isEqual(*intersect1);
      bool i1_2E=_e2.getEndNode()->isEqual(*intersect1);
      ret.push_back(IntersectElement(getE1().getCharactValue(*intersect1),getE2().getCharactValue(*intersect1),i1_1S,i1_1E,i1_2S,i1_2E,intersect1,_e1,_e2,keepOrder()));
      //
      double x2=(_cross*_dy/_drSq-Node::sign(_dy)*_dx*determinant)+center[0];
      double y2=(-_cross*_dx/_drSq-fabs(_dy)*determinant)+center[1];
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
      Node *intersect3=new Node(x,y); intersect3->declareOnTangent();
      bool i_1S=_e1.getStartNode()->isEqual(*intersect3);
      bool i_1E=_e1.getEndNode()->isEqual(*intersect3);
      bool i_2S=_e2.getStartNode()->isEqual(*intersect3);
      bool i_2E=_e2.getEndNode()->isEqual(*intersect3);
      ret.push_back(IntersectElement(_e1.getCharactValue(*intersect3),_e2.getCharactValue(*intersect3),i_1S,i_1E,i_2S,i_2E,intersect3,_e1,_e2,keepOrder()));
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
  GetArcOfCirclePassingThru(*_start,*middle,*_end,_center,_radius,_angle,_angle0);
  middle->decrRef();
  updateBounds();
}

EdgeArcCircle::EdgeArcCircle(Node *start, Node *middle, Node *end, bool direction):Edge(start,end, direction)
{
  GetArcOfCirclePassingThru(*_start,*middle,*_end,_center,_radius,_angle,_angle0);
  updateBounds();
}

EdgeArcCircle::EdgeArcCircle(double sX, double sY, double mX, double mY, double eX, double eY):Edge(sX,sY,eX,eY)
{
  double middle[2]; middle[0]=mX; middle[1]=mY;
  GetArcOfCirclePassingThru(*_start,middle,*_end,_center,_radius,_angle,_angle0);
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
  GetArcOfCirclePassingThru(*_start,*newMiddle,*_end,_center,_radius,_angle,_angle0);
  updateBounds();
}

Edge *EdgeArcCircle::buildEdgeLyingOnMe(Node *start, Node *end, bool direction) const
{
  double sx=((*start)[0]-_center[0])/_radius;
  double sy=((*start)[1]-_center[1])/_radius;
  double ex=((*end)[0]-_center[0])/_radius;
  double ey=((*end)[1]-_center[1])/_radius;
  double angle0=GetAbsoluteAngleOfNormalizedVect(direction?sx:ex,direction?sy:ey);
  double deltaAngle=GetAbsoluteAngleOfNormalizedVect(sx*ex+sy*ey,sx*ey-sy*ex);
  if(deltaAngle>0. && _angle<0.)
    deltaAngle-=2.*M_PI;
  else if(deltaAngle<0. && _angle>0.)
    deltaAngle+=2.*M_PI;
  deltaAngle=direction?deltaAngle:-deltaAngle;
  return new EdgeArcCircle(start,end,_center,_radius,angle0,deltaAngle,direction);
}

void EdgeArcCircle::applySimilarity(double xBary, double yBary, double dimChar)
{
  Edge::applySimilarity(xBary,yBary,dimChar);
  _radius/=dimChar;
  _center[0]=(_center[0]-xBary)/dimChar;
  _center[1]=(_center[1]-yBary)/dimChar;
}

void EdgeArcCircle::unApplySimilarity(double xBary, double yBary, double dimChar)
{
  Edge::unApplySimilarity(xBary,yBary,dimChar);
  _radius*=dimChar;
  _center[0]=_center[0]*dimChar+xBary;
  _center[1]=_center[1]*dimChar+yBary;
}

/*!
 * 'eps' is expected to be > 0.
 * 'conn' is of size 3. conn[0] is start id, conn[1] is end id and conn[2] is middle id.
 * 'offset' is typically the number of nodes already existing in global 2D curve mesh. Additionnal coords 'addCoo' ids will be put after the already existing.
 */
void EdgeArcCircle::tesselate(const int *conn, int offset, double eps, std::vector<int>& newConn, std::vector<double>& addCoo) const
{
  newConn.push_back(INTERP_KERNEL::NORM_POLYL);
  int nbOfSubDiv=(int)(fabs(_angle)/eps);
  if(nbOfSubDiv<=2)
    {
      newConn.push_back(conn[0]); newConn.push_back(conn[2]); newConn.push_back(conn[1]);
      return ;
    }
  double signOfAngle=_angle>0.?1.:-1.;
  int offset2=offset+((int)addCoo.size())/2;
  newConn.push_back(conn[0]);
  for(int i=1;i<nbOfSubDiv;i++,offset2++)
    {
      double angle=_angle0+i*eps*signOfAngle;
      newConn.push_back(offset2);
      addCoo.push_back(_center[0]+_radius*cos(angle)); addCoo.push_back(_center[1]+_radius*sin(angle));
    }
  newConn.push_back(conn[1]);
}

EdgeArcCircle *EdgeArcCircle::BuildFromNodes(Node *start, Node *middle, Node *end)
{
  EdgeLin *e1,*e2;
  e1=new EdgeLin(start,middle);
  e2=new EdgeLin(middle,end);
  SegSegIntersector inters(*e1,*e2);
  bool colinearity=inters.areColinears();
  delete e1; delete e2;
  if(colinearity)
    {
      start->decrRef(); middle->decrRef(); end->decrRef();
      return 0;
    }
  else
    {
      EdgeArcCircle *ret=new EdgeArcCircle(start,middle,end);
      start->decrRef(); middle->decrRef(); end->decrRef();
      return ret;
    }
}

/*!
 * Given an \b NON normalized vector 'vect', returns its norm 'normVect' and its
 * angle in ]-Pi,Pi] relative to Ox axe.
 */
double EdgeArcCircle::GetAbsoluteAngle(const double *vect, double& normVect)
{
  normVect=Node::norm(vect);
  return GetAbsoluteAngleOfNormalizedVect(vect[0]/normVect,vect[1]/normVect);
}

/*!
 * Given a \b normalized vector defined by (ux,uy) returns its angle in ]-Pi;Pi].
 * So before using this method ux*ux+uy*uy should as much as possible close to 1.
 * This methods is quite time consuming in order to keep as much as possible precision.
 * It is NOT ALWAYS possible to do that only in one call of acos. Sometimes call to asin is necessary
 * due to imperfection of acos near 0. and Pi (cos x ~ 1-x*x/2.)
 */
double EdgeArcCircle::GetAbsoluteAngleOfNormalizedVect(double ux, double uy)
{
  return atan2(uy, ux);
}

void EdgeArcCircle::GetArcOfCirclePassingThru(const double *start, const double *middle, const double *end, 
                                              double *center, double& radius, double& angleInRad, double& angleInRad0)
{
  double delta=(middle[0]-start[0])*(end[1]-middle[1])-(end[0]-middle[0])*(middle[1]-start[1]);
  double b1=(middle[1]*middle[1]+middle[0]*middle[0]-start[0]*start[0]-start[1]*start[1])/2;
  double b2=(end[1]*end[1]+end[0]*end[0]-middle[0]*middle[0]-middle[1]*middle[1])/2;
  center[0]=((end[1]-middle[1])*b1+(start[1]-middle[1])*b2)/delta;
  center[1]=((middle[0]-end[0])*b1+(middle[0]-start[0])*b2)/delta;
  radius=SafeSqrt((start[0]-center[0])*(start[0]-center[0])+(start[1]-center[1])*(start[1]-center[1]));
  angleInRad0=GetAbsoluteAngleOfNormalizedVect((start[0]-center[0])/radius,(start[1]-center[1])/radius);
  double angleInRadM=GetAbsoluteAngleOfNormalizedVect((middle[0]-center[0])/radius,(middle[1]-center[1])/radius);
  angleInRad=GetAbsoluteAngleOfNormalizedVect(((start[0]-center[0])*(end[0]-center[0])+(start[1]-center[1])*(end[1]-center[1]))/(radius*radius),
      ((start[0]-center[0])*(end[1]-center[1])-(start[1]-center[1])*(end[0]-center[0]))/(radius*radius));
  if(IsAngleNotIn(angleInRad0,angleInRad,angleInRadM))
    angleInRad=angleInRad<0?2*M_PI+angleInRad:angleInRad-2*M_PI;
}

void EdgeArcCircle::dumpInXfigFile(std::ostream& stream, bool direction, int resolution, const Bounds& box) const
{
  stream << "5 1 0 1 ";
  fillXfigStreamForLoc(stream);
  stream << " 7 50 -1 -1 0.000 0 ";
  if( (direction && (-_angle)>=0) || (!direction && (-_angle)<0))
    stream << '0';//'0'
  else
    stream << '1';//'1'
  stream << " 1 0 ";
  stream << box.fitXForXFigD(_center[0],resolution) << " " << box.fitYForXFigD(_center[1],resolution) << " ";
  direction?_start->dumpInXfigFile(stream,resolution,box):_end->dumpInXfigFile(stream,resolution,box);
  Node *middle=buildRepresentantOfMySelf();
  middle->dumpInXfigFile(stream,resolution,box);
  middle->decrRef();
  direction?_end->dumpInXfigFile(stream,resolution,box):_start->dumpInXfigFile(stream,resolution,box);
  stream << std::endl << "1 1 2.00 120.00 180.00" << std::endl;
}

void EdgeArcCircle::update(Node *m)
{
  GetArcOfCirclePassingThru(*_start,*m,*_end,_center,_radius,_angle,_angle0);
  updateBounds();
}

/*!
 * This methods computes :
 * \f[
 * \int_{Current Edge} -ydx
 * \f]
 */
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
  bary[0]=_center[0]+_radius*cos(_angle0+_angle/2.);
  bary[1]=_center[1]+_radius*sin(_angle0+_angle/2.);
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
 * x=x_{0}+Radius \cdot cos(\theta)
 * \f]
 * \f[
 * y=y_{0}+Radius \cdot sin(\theta)
 * \f]
 * \f[
 * dx=-Radius \cdot sin(\theta) \cdot d\theta
 * \f]
 */
void EdgeArcCircle::getBarycenterOfZone(double *bary) const
{
  double x0=_center[0];
  double y0=_center[1];
  double angle1=_angle0+_angle;
  double tmp1=sin(angle1);
  double tmp0=sin(_angle0);
  double tmp2=_radius*_radius*_radius;
  double tmp3=cos(angle1);
  double tmp4=cos(_angle0);
  bary[0]=_radius*x0*y0*(tmp4-tmp3)+_radius*_radius*(y0*(cos(2*_angle0)-cos(2*angle1))/4.+
      x0*(_angle/2.+(sin(2.*_angle0)-sin(2.*angle1))/4.))
      +tmp2*(tmp1*tmp1*tmp1-tmp0*tmp0*tmp0)/3.;
  bary[1]=y0*y0*_radius*(tmp4-tmp3)/2.+_radius*_radius*y0*(_angle/2.+(sin(2.*_angle0)-sin(2.*angle1))/4.)
        +tmp2*(tmp4-tmp3+(tmp3*tmp3*tmp3-tmp4*tmp4*tmp4)/3.)/2.;
}

/**
 * Compute the "middle" of two points on the arc of circle.
 * The order (p1,p2) or (p2,p1) doesn't matter. p1 and p2 have to be localized on the edge defined by this.
 * \param[out] mid the point located half-way between p1 and p2 on the arc defined by this.
 * \sa getMiddleOfPointsOriented() a generalisation working also when p1 and p2 are not on the arc.
 */
void EdgeArcCircle::getMiddleOfPoints(const double *p1, const double *p2, double *mid) const
{
  double dx1((p1[0]-_center[0])/_radius),dy1((p1[1]-_center[1])/_radius),dx2((p2[0]-_center[0])/_radius),dy2((p2[1]-_center[1])/_radius);
  double angle1(GetAbsoluteAngleOfNormalizedVect(dx1,dy1)),angle2(GetAbsoluteAngleOfNormalizedVect(dx2,dy2));
  //
  double myDelta1(angle1-_angle0),myDelta2(angle2-_angle0);
  if(_angle>0.)
    { myDelta1=myDelta1>-QuadraticPlanarPrecision::getPrecision()?myDelta1:myDelta1+2.*M_PI; myDelta2=myDelta2>-QuadraticPlanarPrecision::getPrecision()?myDelta2:myDelta2+2.*M_PI; }
  else
    { myDelta1=myDelta1<QuadraticPlanarPrecision::getPrecision()?myDelta1:myDelta1-2.*M_PI; myDelta2=myDelta2<QuadraticPlanarPrecision::getPrecision()?myDelta2:myDelta2-2.*M_PI; }
  ////
  mid[0]=_center[0]+_radius*cos(_angle0+(myDelta1+myDelta2)/2.);
  mid[1]=_center[1]+_radius*sin(_angle0+(myDelta1+myDelta2)/2.);
}

/**
 * Compute the "middle" of two points on the arc of circle.
 * Walk on the circle from p1 to p2 using the rotation direction indicated by this->_angle (i.e. by the orientation of the arc).
 * This function is sensitive to the ordering of p1 and p2.
 * \param[out] mid the point located half-way between p1 and p2
 * \sa getMiddleOfPoints() to be used when the order of p1 and p2 is not relevant.
 */
void EdgeArcCircle::getMiddleOfPointsOriented(const double *p1, const double *p2, double *mid) const
{
  double dx1((p1[0]-_center[0])/_radius),dy1((p1[1]-_center[1])/_radius),dx2((p2[0]-_center[0])/_radius),dy2((p2[1]-_center[1])/_radius);
  double angle1(GetAbsoluteAngleOfNormalizedVect(dx1,dy1)),angle2(GetAbsoluteAngleOfNormalizedVect(dx2,dy2));

  if (angle1 <= 0.0)
    angle1 += 2.*M_PI;
  if (angle2 <= 0.0)
    angle2 += 2.*M_PI;

  double avg;
  if((_angle>0. && angle1 <= angle2) || (_angle<=0. && angle1 >= angle2))
    avg = (angle1+angle2)/2.;
  else
    avg = (angle1+angle2)/2. - M_PI;

  mid[0]=_center[0]+_radius*cos(avg);
  mid[1]=_center[1]+_radius*sin(avg);
}


/*!
 * Characteristic value used is angle in ]_Pi;Pi[ from axe 0x.
 */
bool EdgeArcCircle::isIn(double characterVal) const
{
  return IsIn2Pi(_angle0,_angle,characterVal);
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
      myDelta1=myDelta1>-(_radius*QuadraticPlanarPrecision::getPrecision())?myDelta1:myDelta1+2.*M_PI;//in some cases val1 or val2 are so close to angle0 that myDelta is close to 0. but negative.
      myDelta2=myDelta2>-(_radius*QuadraticPlanarPrecision::getPrecision())?myDelta2:myDelta2+2.*M_PI;
      return myDelta1<myDelta2;
    }
  else
    {
      myDelta1=myDelta1<(_radius*QuadraticPlanarPrecision::getPrecision())?myDelta1:myDelta1-2.*M_PI;
      myDelta2=myDelta2<(_radius*QuadraticPlanarPrecision::getPrecision())?myDelta2:myDelta2-2.*M_PI;
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
  return GetAbsoluteAngleOfNormalizedVect(dx,dy);
}

double EdgeArcCircle::getCharactValueBtw0And1(const Node& node) const
{
  double dx=(node[0]-_center[0])/_radius;
  double dy=(node[1]-_center[1])/_radius;
  double angle=GetAbsoluteAngleOfNormalizedVect(dx,dy);
  //
  double myDelta=angle-_angle0;
  if(_angle>0.)
    myDelta=myDelta>=0.?myDelta:myDelta+2.*M_PI;
  else
    myDelta=myDelta<=0.?myDelta:myDelta-2.*M_PI;
  return myDelta/_angle;
}

double EdgeArcCircle::getDistanceToPoint(const double *pt) const
{
  double angle=Node::computeAngle(_center,pt);
  if(IsIn2Pi(_angle0,_angle,angle))
    return fabs(Node::distanceBtw2Pt(_center,pt)-_radius);
  else
    {
      double dist1=Node::distanceBtw2Pt(*_start,pt);
      double dist2=Node::distanceBtw2Pt(*_end,pt);
      return std::min(dist1,dist2);
    }
}

bool EdgeArcCircle::isNodeLyingOn(const double *coordOfNode) const
{
  double dist=Node::distanceBtw2Pt(_center,coordOfNode);
  if(Node::areDoubleEquals(dist,_radius))
    {
      double angle=Node::computeAngle(_center,coordOfNode);
      return IsIn2Pi(_angle0,_angle,angle);
    }
  else
    return false;
}

/*!
 * Idem IsAngleNotIn except that here 'start' in ]-Pi;Pi[ and delta in ]-2*Pi;2Pi[. 
 * @param angleIn in ]-Pi;Pi[.
 */
bool EdgeArcCircle::IsIn2Pi(double start, double delta, double angleIn)
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
bool EdgeArcCircle::IsAngleNotIn(double start, double delta, double angleIn)
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
    return (tmp2<std::min(tmp,tmp+delta) || tmp2>std::max(tmp,tmp+delta));
  else
    return (tmp2>tmp) && (tmp2<(tmp+delta+2.*M_PI));
}

void EdgeArcCircle::updateBounds()
{
  _bounds.setValues(std::min((*_start)[0],(*_end)[0]),std::max((*_start)[0],(*_end)[0]),std::min((*_start)[1],(*_end)[1]),std::max((*_start)[1],(*_end)[1]));
  if(IsIn2Pi(_angle0,_angle,M_PI/2))
    _bounds[3]=_center[1]+_radius;
  if(IsIn2Pi(_angle0,_angle,-M_PI/2))
    _bounds[2]=_center[1]-_radius;
  if(IsIn2Pi(_angle0,_angle,0.))
    _bounds[1]=_center[0]+_radius;
  if(IsIn2Pi(_angle0,_angle,M_PI))
    _bounds[0]=_center[0]-_radius;
}

void EdgeArcCircle::fillGlobalInfoAbs(bool direction, const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
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

void EdgeArcCircle::fillGlobalInfoAbs2(const std::map<INTERP_KERNEL::Node *,int>& mapThis, const std::map<INTERP_KERNEL::Node *,int>& mapOther, int offset1, int offset2, double fact, double baryX, double baryY,
                                       std::vector<int>& edgesOther, std::vector<double>& addCoo, std::map<INTERP_KERNEL::Node *,int>& mapAddCoo) const
{
  _start->fillGlobalInfoAbs2(mapThis,mapOther,offset1,offset2,fact,baryX,baryY,addCoo,mapAddCoo,edgesOther);
  _end->fillGlobalInfoAbs2(mapThis,mapOther,offset1,offset2,fact,baryX,baryY,addCoo,mapAddCoo,edgesOther);
}
