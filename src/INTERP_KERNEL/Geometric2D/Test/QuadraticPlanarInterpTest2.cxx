#include "QuadraticPlanarInterpTest.hxx"
#include "QuadraticPolygon.hxx"
#include "EdgeArcCircle.hxx"
#include "EdgeLin.hxx"

#include <cmath>
#include <sstream>
#include <iostream>

using namespace std;
using namespace INTERP_KERNEL;

static const double ADMISSIBLE_ERROR = 1.e-14;

void QuadraticPlanarInterpTest::IntersectArcCircleBase()
{
  double center[2]={0.5,0.5};
  double radius=0.3;
  EdgeArcCircle *e1=buildArcOfCircle(center,radius,M_PI/4.,M_PI/3.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[0],center[0]+radius*cos(M_PI/3),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[1],center[0]+radius*cos(M_PI/4),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[2],center[1]+radius*sin(M_PI/4),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[3],center[1]+radius*sin(M_PI/3),ADMISSIBLE_ERROR);
  e1->decrRef();
  //
  e1=buildArcOfCircle(center,radius,M_PI/3.,M_PI/2.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[0],center[0]+radius*cos(M_PI/2),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[1],center[0]+radius*cos(M_PI/3),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[2],center[1]+radius*sin(M_PI/3),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[3],center[1]+radius*sin(M_PI/2),ADMISSIBLE_ERROR);
  e1->decrRef();
  //
  e1=buildArcOfCircle(center,radius,M_PI/3.,3.*M_PI/4.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[0],center[0]+radius*cos(3*M_PI/4),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[1],center[0]+radius*cos(M_PI/3),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[2],center[1]+radius*sin(3*M_PI/4),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[3],center[1]+radius*sin(M_PI/2),ADMISSIBLE_ERROR);//<<
  e1->decrRef();
  //
  e1=buildArcOfCircle(center,radius,3*M_PI/4,7*M_PI/8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[0],center[0]+radius*cos(7*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[1],center[0]+radius*cos(3*M_PI/4),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[2],center[1]+radius*sin(7*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[3],center[1]+radius*sin(3*M_PI/4),ADMISSIBLE_ERROR);
  e1->decrRef();
  //
  e1=buildArcOfCircle(center,radius,7.*M_PI/8.,9.*M_PI/8.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[0],center[0]+radius*cos(M_PI),ADMISSIBLE_ERROR);//<<
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[1],center[0]+radius*cos(7*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[2],center[1]+radius*sin(9*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[3],center[1]+radius*sin(7*M_PI/8),ADMISSIBLE_ERROR);
  e1->decrRef();
  //
  e1=buildArcOfCircle(center,radius,9.*M_PI/8.,11.*M_PI/8.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[0],center[0]+radius*cos(9*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[1],center[0]+radius*cos(11*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[2],center[1]+radius*sin(11*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[3],center[1]+radius*sin(9*M_PI/8),ADMISSIBLE_ERROR);
  e1->decrRef();
  //
  e1=buildArcOfCircle(center,radius,11.*M_PI/8.,7.*M_PI/4.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[0],center[0]+radius*cos(11*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[1],center[0]+radius*cos(7*M_PI/4),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[2],center[1]+radius*sin(3*M_PI/2),ADMISSIBLE_ERROR);//<<
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[3],center[1]+radius*sin(7*M_PI/4),ADMISSIBLE_ERROR);
  e1->decrRef();
  //
   e1=buildArcOfCircle(center,radius,7.*M_PI/4.,15.*M_PI/8.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[0],center[0]+radius*cos(7*M_PI/4),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[1],center[0]+radius*cos(15*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[2],center[1]+radius*sin(7*M_PI/4),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[3],center[1]+radius*sin(15*M_PI/8),ADMISSIBLE_ERROR);
  e1->decrRef();
  //
  e1=buildArcOfCircle(center,radius,-M_PI/8.,M_PI/4.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[0],center[0]+radius*cos(M_PI/4),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[1],center[0]+radius*cos(0.),ADMISSIBLE_ERROR);      //<<
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[2],center[1]+radius*sin(15*M_PI/8),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getBounds()[3],center[1]+radius*sin(M_PI/4),ADMISSIBLE_ERROR);
  e1->decrRef();
  //
  // ArcCArcCIntersector
  //
  TypeOfLocInEdge where1,where2;
  vector<Node *> v4;
  MergePoints v3;
  EdgeArcCircle *e2;
  ArcCArcCIntersector *intersector=0;
  for(unsigned k=0;k<8;k++)
    {
      e1=buildArcOfCircle(center,radius,M_PI/4.+k*M_PI/4.,M_PI/3.+k*M_PI/4.);
      e2=buildArcOfCircle(center,radius,M_PI/4.+k*M_PI/4.,M_PI/3.+k*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      intersector->getPlacements(e2->getStartNode(),e2->getEndNode(),where1,where2,v3);
      CPPUNIT_ASSERT(where1==START && where2==END);
      delete intersector; v3.clear(); e2->decrRef();
      //
      e2=buildArcOfCircle(center,radius,7*M_PI/24.+k*M_PI/4.,M_PI/3.+k*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      intersector->getPlacements(e2->getStartNode(),e2->getEndNode(),where1,where2,v3);
      CPPUNIT_ASSERT(where1==INSIDE && where2==END);
      delete intersector; v3.clear(); e2->decrRef();
      //
      e2=buildArcOfCircle(center,radius,M_PI/4.+k*M_PI/4.,7*M_PI/24.+k*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      intersector->getPlacements(e2->getStartNode(),e2->getEndNode(),where1,where2,v3);
      CPPUNIT_ASSERT(where1==START && where2==INSIDE);
      delete intersector; v3.clear(); e2->decrRef();
      //
      e2=buildArcOfCircle(center,radius,13.*M_PI/48.+k*M_PI/4.,15*M_PI/48.+k*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      intersector->getPlacements(e2->getStartNode(),e2->getEndNode(),where1,where2,v3);
      CPPUNIT_ASSERT(where1==INSIDE && where2==INSIDE);
      delete intersector; v3.clear(); e2->decrRef();
      //
      e2=buildArcOfCircle(center,radius,-M_PI/4.+k*M_PI/4.,M_PI/6.+k*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      intersector->getPlacements(e2->getStartNode(),e2->getEndNode(),where1,where2,v3);
      CPPUNIT_ASSERT(where1==OUT_BEFORE && where2==OUT_BEFORE);
      delete intersector; v3.clear(); e2->decrRef();
      //
      e2=buildArcOfCircle(center,radius,0+k*M_PI/4.,5*M_PI/6.+k*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      intersector->getPlacements(e2->getStartNode(),e2->getEndNode(),where1,where2,v3);
      CPPUNIT_ASSERT(where1==OUT_BEFORE && where2==OUT_AFTER);
      delete intersector; v3.clear(); e2->decrRef();
      e1->decrRef();
    }
  // Ok now let's see intersection only. 2 intersections R1 > R2 ; dist(circle1,circle2)>R1; Opposite order.
  for(unsigned k=0;k<8;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,3.,(k-1)*M_PI/4.,(k+1)*M_PI/4.);
      e2=buildArcOfCircle(center2,1.,M_PI+(k-1)*M_PI/4.,M_PI+(k+1)*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(!order);
      CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT(!v4[0]->isEqual(*v4[1]));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(btw2NodesAndACenter(*v4[0],*v4[1],e1->getCenter()),0.35587863972199624,1e-10);
      for(vector<Node *>::iterator iter=v4.begin();iter!=v4.end();iter++)
        (*iter)->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // Ok now let's see intersection only. 2 intersections R1 > R2 ; dist(circle1,circle2)>R1; Same order.
  for(unsigned k=0;k<7;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,3.,(k-1)*M_PI/4.,(k+1)*M_PI/4.);
      e2=buildArcOfCircle(center2,1.,M_PI+(k+1)*M_PI/4.,M_PI+(k-1)*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order);
      CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT(!v4[0]->isEqual(*v4[1]));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(btw2NodesAndACenter(*v4[0],*v4[1],e1->getCenter()),0.35587863972199624,1e-10);
      for(vector<Node *>::iterator iter=v4.begin();iter!=v4.end();iter++)
        (*iter)->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // 2 intersections R1>R2 ; dist(circle1,circle2)<R1; Same order.
  for(unsigned k=0;k<8;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=2.8*cos(k*M_PI/4.); center2[1]=2.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,3.,(k-1)*M_PI/4.,(k+1)*M_PI/4.);
      e2=buildArcOfCircle(center2,1.,(k)*M_PI/4.-M_PI/2.,(k)*M_PI/4.+M_PI/2.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order);
      CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT(!v4[0]->isEqual(*v4[1]));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(btw2NodesAndACenter(*v4[0],*v4[1],e1->getCenter()),0.6793851523346941,1e-10);
      for(vector<Node *>::iterator iter=v4.begin();iter!=v4.end();iter++)
        (*iter)->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // 2 intersections R1>R2 ; dist(circle1,circle2)<R1; Opp order.
  for(unsigned k=0;k<8;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=2.8*cos(k*M_PI/4.); center2[1]=2.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,3.,(k-1)*M_PI/4.,(k+1)*M_PI/4.);
      e2=buildArcOfCircle(center2,1.,(k)*M_PI/4.+M_PI/2.,(k)*M_PI/4.-M_PI/2.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(!order);
      CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT(!v4[0]->isEqual(*v4[1]));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(btw2NodesAndACenter(*v4[0],*v4[1],e1->getCenter()),0.6793851523346941,1e-10);
      for(vector<Node *>::iterator iter=v4.begin();iter!=v4.end();iter++)
	(*iter)->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // Ok now let's see intersection only. 2 intersections R1 < R2 ; dist(circle1,circle2)>R2; Opposite order.
  for(unsigned k=0;k<1;k++)
    {
      double center2[2]; center[0]=0.; center[1]=0.;
      center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,1.,(k-1)*M_PI/4.,(k+1)*M_PI/4.);
      e2=buildArcOfCircle(center2,3.,M_PI+(k-1)*M_PI/4.,M_PI+(k+1)*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(!order);
      CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT(!v4[0]->isEqual(*v4[1]));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.1195732971845034,btw2NodesAndACenter(*v4[0],*v4[1],e1->getCenter()),1e-10);
      for(vector<Node *>::iterator iter=v4.begin();iter!=v4.end();iter++)
        (*iter)->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // Ok now let's see intersection only. 2 intersections R1 < R2 ; dist(circle1,circle2)>R2; same order.
  for(unsigned k=0;k<8;k++)
    {
      double center2[2]; center[0]=0.; center[1]=0.;
      center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,1.,(k+1)*M_PI/4.,(k-1)*M_PI/4.);
      e2=buildArcOfCircle(center2,3.,M_PI+(k-1)*M_PI/4.,M_PI+(k+1)*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order);
      CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT(!v4[0]->isEqual(*v4[1]));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.1195732971845034,btw2NodesAndACenter(*v4[0],*v4[1],e1->getCenter()),1e-10);
      for(vector<Node *>::iterator iter=v4.begin();iter!=v4.end();iter++)
        (*iter)->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // Ok now let's see intersection only. 2 intersections R1 < R2 ; dist(circle1,circle2)<R2; same order.
  for(unsigned k=0;k<8;k++)
    {
      double center2[2]; center[0]=0.; center[1]=0.;
      center2[0]=-2.8*cos(k*M_PI/4.); center2[1]=-2.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,1.,(k)*M_PI/4.+M_PI/2.,(k)*M_PI/4.-M_PI/2.);
      e2=buildArcOfCircle(center2,3.,(k+1)*M_PI/4.,(k-1)*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order);
      CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT(!v4[0]->isEqual(*v4[1]));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0844420190512074,btw2NodesAndACenter(*v4[0],*v4[1],e1->getCenter()),1e-10);
      for(vector<Node *>::iterator iter=v4.begin();iter!=v4.end();iter++)
        (*iter)->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // Ok now let's see intersection only. 2 intersections R1 < R2 ; dist(circle1,circle2)<R2; opp. order.
  for(unsigned k=0;k<8;k++)
    {
      double center2[2]; center[0]=0.; center[1]=0.;
      center2[0]=-2.8*cos(k*M_PI/4.); center2[1]=-2.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,1.,(k)*M_PI/4.+M_PI/2.,(k)*M_PI/4.-M_PI/2.);
      e2=buildArcOfCircle(center2,3.,(k-1)*M_PI/4.,(k+1)*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(!order);
      CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[1]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT(!v4[0]->isEqual(*v4[1]));
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.0844420190512074,btw2NodesAndACenter(*v4[0],*v4[1],e1->getCenter()),1e-10);
      for(vector<Node *>::iterator iter=v4.begin();iter!=v4.end();iter++)
        (*iter)->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // Tangent intersection
  for(unsigned k=0;k<8;k++)
    {
      double center2[2]; center[0]=0.; center[1]=0.;
      center2[0]=4.*cos(k*M_PI/4.); center2[1]=4.*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,1.,(k+1)*M_PI/4.,(k-1)*M_PI/4.);
      e2=buildArcOfCircle(center2,3.,M_PI+(k-1)*M_PI/4.,M_PI+(k+1)*M_PI/4.);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order); // order has no sence here because v4.size() expected to 1 but for valgrind serenity test.
      CPPUNIT_ASSERT_EQUAL(1,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getRadius(),Node::distanceBtw2Pt(e1->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(e2->getRadius(),Node::distanceBtw2Pt(e2->getCenter(),(*(v4[0]))),ADMISSIBLE_ERROR);
      for(vector<Node *>::iterator iter=v4.begin();iter!=v4.end();iter++)
        (*iter)->decrRef();
      v4.clear(); v4.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // Extremities # 1
  for(unsigned k=0;k<8;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,3.,k*M_PI/4.-0.17793931986099812,k*M_PI/4.+0.17793931986099812);
      e2=buildArcOfCircle(center2,1.,M_PI+k*M_PI/4.-0.55978664859225125,M_PI+k*M_PI/4.+0.55978664859225125);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(!intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT_EQUAL(0,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(2,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT(e1->getStartNode()==e2->getEndNode()); CPPUNIT_ASSERT(e2->getStartNode()==e1->getEndNode());
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  for(unsigned k=0;k<8;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,3.,k*M_PI/4.-0.17793931986099812,k*M_PI/4.+0.17793931986099812);
      e2=buildArcOfCircle(center2,1.,M_PI+k*M_PI/4.+0.55978664859225125,M_PI+k*M_PI/4.-0.55978664859225125);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(!intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT_EQUAL(0,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(2,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(e2->getEndNode()==e1->getEndNode());
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // Extremities # 2
  for(unsigned k=0;k<8;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,3.,k*M_PI/4.-0.17793931986099812,k*M_PI/4.+0.17793931986099812);
      e2=buildArcOfCircle(center2,1.,M_PI+k*M_PI/4.+0.55978664859225125,M_PI+k*M_PI/4.-0.7);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3));
      CPPUNIT_ASSERT(order); CPPUNIT_ASSERT_EQUAL(1,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v4[0]);
      v4[0]->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
      }
  // Extremities # 3
  for(unsigned k=0;k<8;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      e1=buildArcOfCircle(center,3.,k*M_PI/4.-0.17793931986099812,k*M_PI/4.+0.17793931986099812);
      e2=buildArcOfCircle(center2,1.,M_PI+k*M_PI/4.+0.7,M_PI+k*M_PI/4.-0.7);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order); CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT(e1->getStartNode()==v4[0]); CPPUNIT_ASSERT(e1->getEndNode()==v4[1]);
      v4[0]->decrRef(); v4[1]->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  // Extremities # 4
  for(unsigned k=0;k<8;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      Node *nodeS=new Node(center[0]+3.*cos(k*M_PI/4.-0.17793931986099812),center[1]+3.*sin(k*M_PI/4.-0.17793931986099812));
      Node *nodeE=new Node(center[0]+3.*cos(k*M_PI/4.),center[1]+3.*sin(k*M_PI/4.));
      double angle=k*M_PI/4.-0.17793931986099812;
      angle=angle>M_PI?angle-2.*M_PI:angle;
      e1=new EdgeArcCircle(nodeS,nodeE,//Problem of precision 1e-14 to easily reached.
			   center,3.,angle,0.17793931986099812);
      nodeS->decrRef(); nodeE->decrRef();
      e2=buildArcOfCircle(center2,1.,M_PI+k*M_PI/4.+0.7,M_PI+k*M_PI/4.-0.7);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order); CPPUNIT_ASSERT_EQUAL(1,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT(e1->getStartNode()==v4[0]);
      v4[0]->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
  //Extremities # 5
  for(unsigned k=0;k<8;k++)
    {
      center[0]=0.; center[1]=0.;
      double center2[2]; center2[0]=3.8*cos(k*M_PI/4.); center2[1]=3.8*sin(k*M_PI/4.);
      Node *nodeS=new Node(center[0]+3.*cos(k*M_PI/4.-0.17793931986099812),center[1]+3.*sin(k*M_PI/4.-0.17793931986099812));
      Node *nodeE=new Node(center[0]+3.*cos(k*M_PI/4.)+0.5,center[1]+3.*sin(k*M_PI/4.));
      double angle=k*M_PI/4.-0.17793931986099812;
      angle=angle>M_PI?angle-2.*M_PI:angle;
      e1=new EdgeArcCircle(nodeS,nodeE,//Problem of precision 1e-14 to easily reached.
			   center,3.,angle,0.67793931986099812);
      nodeS->decrRef(); nodeE->decrRef();
      e2=buildArcOfCircle(center2,1.,M_PI+k*M_PI/4.+0.7,M_PI+k*M_PI/4.-0.7);
      intersector=new ArcCArcCIntersector(*e1,*e2);
      bool order;
      bool obvious,areOverlapped;
      intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
      CPPUNIT_ASSERT(!obvious && !areOverlapped);
      CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order); CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
      CPPUNIT_ASSERT(e1->getStartNode()==v4[0]);
      v4[0]->decrRef(); v4[1]->decrRef();
      v4.clear(); v3.clear();
      delete intersector; e2->decrRef(); e1->decrRef();
    }
}

void QuadraticPlanarInterpTest::IntersectArcCircleFull()
{
  double center1[2]; center1[0]=0.;   center1[1]=0.;   double radius1=3.;
  double center2[2]; center2[0]=0.75; center2[1]=-2.6; double radius2=1.;
  EdgeArcCircle *e1=buildArcOfCircle(center1,radius1,-M_PI/3.,4.*M_PI/3.);
  EdgeArcCircle *e2=buildArcOfCircle(center2,radius2,0.,M_PI/2.);
  MergePoints commonNode;
  QuadraticPolygon pol1; QuadraticPolygon pol2;
  QuadraticPolygon pol3; QuadraticPolygon pol4;
  pol3.pushBack(e1); pol4.pushBack(e2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol3.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5707963267949,pol4.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.6648305849,pol3.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.8146018366,pol4.getAreaFast(),1e-6);
  CPPUNIT_ASSERT(e1->intersectWith(e2,commonNode,pol1,pol2));
  CPPUNIT_ASSERT_EQUAL(2,pol1.size());
  CPPUNIT_ASSERT_EQUAL(2,pol2.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(19.6648305849,pol1.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.8146018366,pol2.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol1.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5707963267949,pol2.getPerimeterFast(),1e-6);
  //
  e1=buildArcOfCircle(center1,radius1,-2*M_PI/3.,-7.*M_PI/3.);
  e2=buildArcOfCircle(center2,radius2,0.,M_PI/2.);
  commonNode.clear();
  QuadraticPolygon pol5; QuadraticPolygon pol6;
  QuadraticPolygon pol7; QuadraticPolygon pol8;
  pol7.pushBack(e1); pol8.pushBack(e2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol7.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5707963267949,pol8.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-19.6648305849,pol7.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.8146018366,pol8.getAreaFast(),1e-6);
  CPPUNIT_ASSERT(e1->intersectWith(e2,commonNode,pol5,pol6));
  CPPUNIT_ASSERT_EQUAL(2,pol5.size());
  CPPUNIT_ASSERT_EQUAL(2,pol6.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-19.6648305849,pol5.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.8146018366,pol6.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol5.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5707963267949,pol6.getPerimeterFast(),1e-6);
  //
  center2[0]=3.5; center2[1]=0.;
  e1=buildArcOfCircle(center1,radius1,-2*M_PI/3.,-7.*M_PI/3.);
  e2=buildArcOfCircle(center2,radius2,M_PI/2.,3*M_PI/2.);
  commonNode.clear();
  QuadraticPolygon pol9; QuadraticPolygon pol10;
  QuadraticPolygon pol11; QuadraticPolygon pol12;
  pol11.pushBack(e1); pol12.pushBack(e2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol11.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.1415926535897931,pol12.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-19.6648305849,pol11.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5707963267949,pol12.getAreaFast(),1e-6);
  CPPUNIT_ASSERT(e1->intersectWith(e2,commonNode,pol9,pol10));
  CPPUNIT_ASSERT_EQUAL(3,pol9.size());
  CPPUNIT_ASSERT_EQUAL(3,pol10.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol9.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.1415926535897931,pol10.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-19.6648305849,pol9.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5707963267949,pol10.getAreaFast(),1e-6);
  //
  center2[0]=0.; center2[1]=0.; radius2=radius1;
  e1=buildArcOfCircle(center1,radius1,-2*M_PI/3.,-7.*M_PI/3.);
  e2=buildArcOfCircle(center2,radius2,M_PI/3.,2*M_PI/3.);
  commonNode.clear();
  QuadraticPolygon pol13; QuadraticPolygon pol14;
  QuadraticPolygon pol15; QuadraticPolygon pol16;
  pol15.pushBack(e1); pol16.pushBack(e2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol15.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.1415926535897931,pol16.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-19.6648305849,pol15.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.6095032974147,pol16.getAreaFast(),1e-6);
  CPPUNIT_ASSERT(e1->intersectWith(e2,commonNode,pol13,pol14));
  CPPUNIT_ASSERT_EQUAL(3,pol13.size());
  CPPUNIT_ASSERT_EQUAL(1,pol14.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol13.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-19.6648305849,pol13.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.1415926535897931,pol14.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.6095032974147,pol14.getAreaFast(),1e-6);
  //
  e1=buildArcOfCircle(center1,radius1,-2*M_PI/3.,-7.*M_PI/3.);
  e2=buildArcOfCircle(center2,radius2,2*M_PI/3.,M_PI/3.);
  commonNode.clear();
  QuadraticPolygon pol17; QuadraticPolygon pol18;
  QuadraticPolygon pol19; QuadraticPolygon pol20;
  pol19.pushBack(e1); pol20.pushBack(e2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol19.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.1415926535897931,pol20.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-19.6648305849,pol19.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-8.6095032974147,pol20.getAreaFast(),1e-6);
  CPPUNIT_ASSERT(e1->intersectWith(e2,commonNode,pol17,pol18));
  CPPUNIT_ASSERT_EQUAL(3,pol17.size());
  CPPUNIT_ASSERT_EQUAL(1,pol18.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15.707963267948966,pol17.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-19.6648305849,pol17.getAreaFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.1415926535897931,pol18.getPerimeterFast(),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-8.6095032974147,pol18.getAreaFast(),1e-6);
  //no intersection #1
  center2[0]=4.277; center2[1]=-4.277;
  e1=buildArcOfCircle(center1,radius1,-2*M_PI/3.,-7.*M_PI/3.);
  e2=buildArcOfCircle(center2,radius2,M_PI/4.,5*M_PI/4.);
  QuadraticPolygon polTemp1; QuadraticPolygon polTemp2;
  CPPUNIT_ASSERT(!e1->intersectWith(e2,commonNode,polTemp1,polTemp2));
  e1->decrRef(); e2->decrRef();
  //no intersection #2
  center2[0]=1.; center2[1]=-1.; radius2=0.2;
  e1=buildArcOfCircle(center1,radius1,-2*M_PI/3.,-7.*M_PI/3.);
  e2=buildArcOfCircle(center2,radius2,M_PI/4.,5*M_PI/4.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,commonNode,polTemp1,polTemp2));
  e1->decrRef(); e2->decrRef();
}

void QuadraticPlanarInterpTest::IntersectArcCircleSegumentBase()
{
  double center[2]={2.,2.};
  EdgeArcCircle *e1=buildArcOfCircle(center,2.3,M_PI/4.,5.*M_PI/4.);
  EdgeLin *e2=new EdgeLin(-1.3,1.,3.,5.3);
  Intersector *intersector=new ArcCSegIntersector(*e1,*e2);
  bool order;
  bool obvious,areOverlapped;
  intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped);
  CPPUNIT_ASSERT(!obvious && !areOverlapped);
  vector<Node *> v4;
  MergePoints v3;
  CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(!order); CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,(*v4[0])[0],1e-10); CPPUNIT_ASSERT_DOUBLES_EQUAL(4.3,(*v4[0])[1],1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.3,(*v4[1])[0],1e-10); CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,(*v4[1])[1],1e-10);
  v4[0]->decrRef(); v4[1]->decrRef(); e2->decrRef(); v3.clear(); v4.clear(); delete intersector;
  //
  e2=new EdgeLin(3.,5.3,-1.3,1.);
  intersector=new ArcCSegIntersector(*e1,*e2);
  intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped); CPPUNIT_ASSERT(!obvious && !areOverlapped);
  CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order); CPPUNIT_ASSERT_EQUAL(2,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,(*v4[0])[0],1e-10); CPPUNIT_ASSERT_DOUBLES_EQUAL(4.3,(*v4[0])[1],1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.3,(*v4[1])[0],1e-10); CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,(*v4[1])[1],1e-10);
  v4[0]->decrRef(); v4[1]->decrRef(); e2->decrRef(); v3.clear(); v4.clear(); delete intersector;
  // tangent intersection
  e2=new EdgeLin(-1.,4.3,3.,4.3);
  intersector=new ArcCSegIntersector(*e1,*e2);
  intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped); CPPUNIT_ASSERT(!obvious && !areOverlapped);
  CPPUNIT_ASSERT(intersector->intersect(0,v4,order,v3)); CPPUNIT_ASSERT(order); CPPUNIT_ASSERT_EQUAL(1,(int)v4.size()); CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.,(*v4[0])[0],1e-10); CPPUNIT_ASSERT_DOUBLES_EQUAL(4.3,(*v4[0])[1],1e-10);
  v4[0]->decrRef(); e2->decrRef(); v3.clear(); delete intersector;
  // no intersection
  e2=new EdgeLin(-2.,-2.,-1.,-3.);
  intersector=new ArcCSegIntersector(*e1,*e2);
  intersector->areOverlappedOrOnlyColinears(0,obvious,areOverlapped); CPPUNIT_ASSERT(obvious && !areOverlapped);
  e2->decrRef(); v3.clear(); delete intersector;
  //
  e1->decrRef();
}

EdgeArcCircle *QuadraticPlanarInterpTest::buildArcOfCircle(const double *center, double radius, double alphaStart, double alphaEnd)
{
  double alphaM=(alphaStart+alphaEnd)/2;
  return new EdgeArcCircle(center[0]+cos(alphaStart)*radius,center[1]+sin(alphaStart)*radius,
			   center[0]+cos(alphaM)*radius,center[1]+sin(alphaM)*radius,
			   center[0]+cos(alphaEnd)*radius,center[1]+sin(alphaEnd)*radius);
}

double QuadraticPlanarInterpTest::btw2NodesAndACenter(const Node& n1, const Node& n2, const double *center)
{
  const double *n1Pt=n1;
  const double *n2Pt=n2;
  double tmp1[2],tmp2[2];
  tmp1[0]=n1Pt[0]-center[0]; tmp1[1]=n1Pt[1]-center[1];
  tmp2[0]=n2Pt[0]-center[0]; tmp2[1]=n2Pt[1]-center[1];
  double distTmp1=sqrt(tmp1[0]*tmp1[0]+tmp1[1]*tmp1[1]);
  double distTmp2=sqrt(tmp2[0]*tmp2[0]+tmp2[1]*tmp2[1]);
  double ret=acos((tmp1[0]*tmp2[0]+tmp1[1]*tmp2[1])/(distTmp1*distTmp2));
  if(tmp1[0]*tmp2[1]-tmp1[1]*tmp2[0]<0)
    ret=-ret;
  return ret;
}
