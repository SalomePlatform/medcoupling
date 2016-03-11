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

#include "QuadraticPlanarInterpTest.hxx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DElementaryEdge.hxx"
#include "InterpKernelGeo2DComposedEdge.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"
#include "TestInterpKernelUtils.hxx"

#include <sstream>
#include <iostream>

using namespace INTERP_KERNEL;

namespace INTERP_TEST
{

static const double ADMISSIBLE_ERROR = 1.e-14;

void QuadraticPlanarInterpTest::setUp()
{
}

void QuadraticPlanarInterpTest::tearDown()
{
}

void QuadraticPlanarInterpTest::cleanUp()
{
}

void QuadraticPlanarInterpTest::ReadWriteInXfigElementary()
{
  //Testing bounds calculation. For Seg2
  std::istringstream stream("2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2\n3200 3400 4500 4700");
  EdgeLin *e1=new EdgeLin(stream);
  Bounds bound=e1->getBounds();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.32,bound[0],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.45,bound[1],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.34,bound[2],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.47,bound[3],ADMISSIBLE_ERROR);
  e1->decrRef();
  std::istringstream stream2("2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2\n4500 4700 3200 3400");
  e1=new EdgeLin(stream2);
  bound=e1->getBounds();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.32,bound[0],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.45,bound[1],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.34,bound[2],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.47,bound[3],ADMISSIBLE_ERROR);
  e1->decrRef();
  //Testing bounds calculation For Arc of circle.
  
}

void QuadraticPlanarInterpTest::ReadWriteInXfigGlobal()
{
  QuadraticPolygon pol1(INTERP_TEST::getResourceFile("Pol1.fig").c_str());
  pol1.dumpInXfigFile("Pol1_gen.fig");
  QuadraticPolygon pol2(INTERP_TEST::getResourceFile("Pol2.fig").c_str());
  pol2.dumpInXfigFile("Pol2_gen.fig");
  QuadraticPolygon pol3(INTERP_TEST::getResourceFile("Pol3.fig").c_str());
  pol3.dumpInXfigFile("Pol3_gen.fig");
  QuadraticPolygon pol4(INTERP_TEST::getResourceFile("Pol4.fig").c_str());
  CPPUNIT_ASSERT_EQUAL(1,pol4.size());
  ElementaryEdge *edge1=dynamic_cast<ElementaryEdge *>(pol4[0]);
  CPPUNIT_ASSERT(edge1);
  Edge *edge2=edge1->getPtr();
  EdgeArcCircle *edge=dynamic_cast<EdgeArcCircle *>(edge2);
  CPPUNIT_ASSERT(edge);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.24375,edge->getRadius(),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.7857653289925404,edge->getAngle(),ADMISSIBLE_ERROR);
  double center[2];
  edge->getCenter(center);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.48,center[0],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.48375,center[1],ADMISSIBLE_ERROR);
  const double *start=*edge->getStartNode();
  Node *n1=new Node(start[0]+2*(center[0]-start[0]),start[1]+2*(center[1]-start[1]));
  edge->changeMiddle(n1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.24375,edge->getRadius(),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.7857653289925404,edge->getAngle(),ADMISSIBLE_ERROR);
  n1->decrRef();
  n1=new Node(center[0],center[1]+0.24375);
  edge->changeMiddle(n1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.24375,edge->getRadius(),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.49741997818704586,edge->getAngle(),ADMISSIBLE_ERROR);//5.7857653289925404 + 2*PI
  n1->decrRef();
  //A half circle.
  EdgeArcCircle *e=new EdgeArcCircle(0.84,0.54,0.78,0.6,0.84,0.66);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.06,e->getRadius(),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-3.1415925921507317,e->getAngle(),1e-5);
  e->decrRef();
  e=new EdgeArcCircle(0.84,0.54,0.9,0.6,0.84,0.66);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.06,e->getRadius(),ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.1415925921507317,e->getAngle(),1e-5);
  e->decrRef();
}

void QuadraticPlanarInterpTest::BasicGeometricTools()
{
  Node *n1=new Node(1.,1.);
  Node *n2=new Node(4.,2.);
  EdgeLin *e1=new EdgeLin(n1,n2);
  double tmp[2];
  e1->getNormalVector(tmp);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.94868329805051377,tmp[1],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.31622776601683794,tmp[0],ADMISSIBLE_ERROR);
  e1->decrRef();
  n1->decrRef(); n2->decrRef();
  n1=new Node(1.,1.);
  n2=new Node(0.,4.);
  e1=new EdgeLin(n1,n2);
  double tmp2[2];
  e1->getNormalVector(tmp2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,Node::dot(tmp,tmp2),1e-10);
  tmp[0]=0.5; tmp[1]=2.5;
  CPPUNIT_ASSERT(e1->isNodeLyingOn(tmp));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,e1->getDistanceToPoint(tmp),1e-12);
  tmp[1]=2.55; CPPUNIT_ASSERT(!e1->isNodeLyingOn(tmp));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.0158113883008418,e1->getDistanceToPoint(tmp),1e-12);
  tmp[0]=0.; tmp[1]=5.;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,e1->getDistanceToPoint(tmp),1e-12);
  EdgeArcCircle *e=new EdgeArcCircle(4.,3.,0.,5.,-5.,0.);
  tmp[0]=-4.; tmp[1]=3.;
  CPPUNIT_ASSERT(e->isNodeLyingOn(tmp));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,e->getDistanceToPoint(tmp),1e-12);
  tmp[1]=3.1; CPPUNIT_ASSERT(!e->isNodeLyingOn(tmp));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.0632371551998077e-2,e->getDistanceToPoint(tmp),1e-12);
  tmp[0]=-4.; tmp[1]=-3.;
  CPPUNIT_ASSERT(!e->isNodeLyingOn(tmp));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.1622776601683795,e->getDistanceToPoint(tmp),1e-12);
  e->decrRef();
  e1->decrRef();
  n1->decrRef(); n2->decrRef();
}

void QuadraticPlanarInterpTest::IntersectionBasics()
{
  //Testing intersection of Bounds.
  std::istringstream stream1("2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2\n3200 3400 4500 4800");
  EdgeLin *e1=new EdgeLin(stream1);
  std::istringstream stream2("2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2\n3200 3400 4500 4800");
  EdgeLin *e2=new EdgeLin(stream2);
  Bounds *bound=e1->getBounds().amIIntersectingWith(e2->getBounds()); CPPUNIT_ASSERT(bound);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.32,(*bound)[0],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.45,(*bound)[1],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.34,(*bound)[2],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.48,(*bound)[3],ADMISSIBLE_ERROR);
  delete bound;
  e2->decrRef(); e1->decrRef();
  //
  std::istringstream stream3("2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2\n3000 7200 6000 3700");
  EdgeLin *e3=new EdgeLin(stream3);
  std::istringstream stream4("2 1 0 1 0 7 50 -1 -1 0.000 0 0 -1 0 0 2\n4800 6600 7200 4200");
  EdgeLin *e4=new EdgeLin(stream4);
  bound=e3->getBounds().amIIntersectingWith(e4->getBounds()); CPPUNIT_ASSERT(bound);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.48,(*bound)[0],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.6,(*bound)[1],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.42,(*bound)[2],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.66,(*bound)[3],ADMISSIBLE_ERROR);
  delete bound;
  e3->decrRef(); e4->decrRef();
}

void QuadraticPlanarInterpTest::EdgeLinUnitary()
{
  EdgeLin *e1=new EdgeLin(0.5,0.5,3.7,4.1);
  Node *n=new Node(2.1,2.3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getCharactValue(*n),0.5,1e-8);
  n->decrRef();
  n=new Node(3.7,4.1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getCharactValue(*n),1.,1e-8);
  n->decrRef();
  n=new Node(0.5,0.5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getCharactValue(*n),0.,1e-8);
  n->decrRef();
  n=new Node(-1.1,-1.3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getCharactValue(*n),-0.5,1e-8);
  n->decrRef();
  n=new Node(5.3,5.9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getCharactValue(*n),1.5,1e-8);
  n->decrRef(); e1->decrRef();
}

/*!
 * Here two things are tested. 
 * 1 ) One the overlapping calculation capability of edge/edge intersector.
 * 2 ) Then the capability to handle the case where 2 segs (whatever their type) are overlapped.
 * All the configuration of full or part overlapping have been tested.
 */
void QuadraticPlanarInterpTest::IntersectionEdgeOverlapUnitarySegSeg()
{
  ComposedEdge& v1=*(new ComposedEdge);
  ComposedEdge& v2=*(new ComposedEdge);
  MergePoints v3;
  //Testing merge of geometric equals seg2.
  Edge *e1=new EdgeLin(0.5,0.5,1.,1.); Edge *e2=new EdgeLin(0.5,0.5,1.,1.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(2,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size()); CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && v2[0]->getDirection());
  v1.clear(); v2.clear(); v3.clear();
  //  - testing by adding some noise
  e1->decrRef(); e1=new EdgeLin(0.5+5.e-15,0.5-5.e-15,1.,1.+7.e-15);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(2,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size()); CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && v2[0]->getDirection());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Testing merge of geometric equals seg2 but now with opposite direction
  e1=new EdgeLin(0.5,0.5,0.7,0.7); e2=new EdgeLin(0.7+6.e-15,0.7-2.e-15,0.5+3.e-15,0.5-4.e-15);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(2,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size()); CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && !v2[0]->getDirection());//compared 8 lines above !v2[0]->getDirection()
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 0
  //Test 1 - OUT_AFTER - OUT_AFTER | same dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(1.5,0.,2.,0.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(0,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(0,(int)v2.size());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 2 - INSIDE - OUT_AFTER | same dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(0.5,0.,1.5,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresicEqualDirSensitive(v2[0]));
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode()); CPPUNIT_ASSERT(e2->getStartNode()==v2[0]->getStartNode()); CPPUNIT_ASSERT(e2->getEndNode()==v2[1]->getEndNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 2 - INSIDE - OUT_AFTER | same dir. - 90°
  e1=new EdgeLin(0.,0.,0.,1.); e2=new EdgeLin(0.,0.5,0.,1.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresicEqualDirSensitive(v2[0]));
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode()); CPPUNIT_ASSERT(e2->getStartNode()==v2[0]->getStartNode()); CPPUNIT_ASSERT(e2->getEndNode()==v2[1]->getEndNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 2 - INSIDE - OUT_AFTER | same dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(0.5,0.5,1.5,1.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresicEqualDirSensitive(v2[0]));
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode()); CPPUNIT_ASSERT(e2->getStartNode()==v2[0]->getStartNode()); CPPUNIT_ASSERT(e2->getEndNode()==v2[1]->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 2 - INSIDE - OUT_AFTER | opp. dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(1.5,1.5,0.5,0.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(!v1[1]->intresicEqualDirSensitive(v2[1]) && v1[1]->intresicEqual(v2[1]));
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode()); CPPUNIT_ASSERT(e2->getStartNode()==v2[0]->getStartNode()); CPPUNIT_ASSERT(e2->getEndNode()==v2[1]->getEndNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 3 - INSIDE - INSIDE | same dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(0.25,0.,0.75,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(3,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresincEqCoarse(e2) && v1[1]->getDirection());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v1[1]->getEndNode()==v1[2]->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()== e1->getStartNode()); CPPUNIT_ASSERT(v1[2]->getEndNode()== e1->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==e2->getStartNode()); CPPUNIT_ASSERT(v1[1]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 3 - INSIDE - INSIDE | same dir. - 90°
  e1=new EdgeLin(0.,0.,0.,1.); e2=new EdgeLin(0.,0.25,0.,0.75);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(3,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresincEqCoarse(e2) && v1[1]->getDirection());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v1[1]->getEndNode()==v1[2]->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()== e1->getStartNode()); CPPUNIT_ASSERT(v1[2]->getEndNode()== e1->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==e2->getStartNode()); CPPUNIT_ASSERT(v1[1]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 3 - INSIDE - INSIDE | same dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(0.25,0.25,0.75,0.75);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(3,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresincEqCoarse(e2) && v1[1]->getDirection());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v1[1]->getEndNode()==v1[2]->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()== e1->getStartNode()); CPPUNIT_ASSERT(v1[2]->getEndNode()== e1->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==e2->getStartNode()); CPPUNIT_ASSERT(v1[1]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 3 - INSIDE - INSIDE | opp dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(0.75,0.75,0.25,0.25);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(3,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresincEqCoarse(e2) && !v1[1]->getDirection());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v1[1]->getEndNode()==v1[2]->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()== e1->getStartNode()); CPPUNIT_ASSERT(v1[2]->getEndNode()== e1->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==e2->getEndNode()); CPPUNIT_ASSERT(v1[1]->getEndNode()==e2->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 3bis - INSIDE - INSIDE - Bis | opp dir.
  double center[2]={0.,0.};
  double radius=1.;
  e1=buildArcOfCircle(center,radius,-M_PI,0); e2=buildArcOfCircle(center,radius,-2*M_PI/3.+2*M_PI,-M_PI/3.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,e1->getCurveLength(),1e-12); CPPUNIT_ASSERT_DOUBLES_EQUAL(5.*M_PI/3.,e2->getCurveLength(),1e-12);// To check that in the previous line +2.M_PI has done its job.
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(3,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(v2[0]->getPtr())); CPPUNIT_ASSERT(v1[0]->getDirection()); CPPUNIT_ASSERT(!v2[0]->getDirection());
  CPPUNIT_ASSERT(v1[2]->intresincEqCoarse(v2[2]->getPtr())); CPPUNIT_ASSERT(v1[2]->getDirection()); CPPUNIT_ASSERT(!v2[2]->getDirection());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,v2[1]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI/3.,v1[1]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT(v2[1]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v2[1]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==e2->getStartNode());
  CPPUNIT_ASSERT(v1[1]->getStartNode()==e2->getStartNode());
  CPPUNIT_ASSERT(v1[1]->getEndNode()==e2->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 3bis - INSIDE - INSIDE - Bis | same dir.
  e1=buildArcOfCircle(center,radius,-M_PI,0); e2=buildArcOfCircle(center,radius,-M_PI/3.,-2*M_PI/3.+2*M_PI);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,e1->getCurveLength(),1e-12); CPPUNIT_ASSERT_DOUBLES_EQUAL(5.*M_PI/3.,e2->getCurveLength(),1e-12);// To check that in the previous line +2.M_PI has done its job.
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(3,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(v2[2]->getPtr())); CPPUNIT_ASSERT(v1[0]->getDirection()); CPPUNIT_ASSERT(v2[2]->getDirection());
  CPPUNIT_ASSERT(v1[2]->intresincEqCoarse(v2[0]->getPtr())); CPPUNIT_ASSERT(v1[2]->getDirection()); CPPUNIT_ASSERT(v2[0]->getDirection());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,v2[1]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI/3.,v1[1]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT(v2[1]->getStartNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v2[1]->getEndNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v1[1]->getStartNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v1[1]->getEndNode()==e2->getStartNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 3bis - INSIDE - INSIDE - Bis | opp dir. | e1<->e2 to test symetry
  e1=buildArcOfCircle(center,radius,-M_PI,0); e2=buildArcOfCircle(center,radius,-2*M_PI/3.+2*M_PI,-M_PI/3.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,e1->getCurveLength(),1e-12); CPPUNIT_ASSERT_DOUBLES_EQUAL(5.*M_PI/3.,e2->getCurveLength(),1e-12);// To check that in the previous line +2.M_PI has done its job.
  CPPUNIT_ASSERT(e2->intersectWith(e1,v3,v2,v1));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(3,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(v2[0]->getPtr())); CPPUNIT_ASSERT(!v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->getDirection());
  CPPUNIT_ASSERT(v1[2]->intresincEqCoarse(v2[2]->getPtr())); CPPUNIT_ASSERT(!v1[2]->getDirection()); CPPUNIT_ASSERT(v2[2]->getDirection());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,v2[1]->getCurveLength(),1.e-5); // << not maximal precision because node switching 
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI/3.,v1[1]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT(v2[1]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v2[1]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==e2->getStartNode());
  CPPUNIT_ASSERT(v1[1]->getStartNode()==e2->getStartNode());
  CPPUNIT_ASSERT(v1[1]->getEndNode()==e2->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 3bis - INSIDE - INSIDE - Bis | same dir. | e1<->e2 to test symetry
  e1=buildArcOfCircle(center,radius,-M_PI,0); e2=buildArcOfCircle(center,radius,-M_PI/3.,-2*M_PI/3.+2*M_PI);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,e1->getCurveLength(),1e-12); CPPUNIT_ASSERT_DOUBLES_EQUAL(5.*M_PI/3.,e2->getCurveLength(),1e-12);// To check that in the previous line +2.M_PI has done its job.
  CPPUNIT_ASSERT(e2->intersectWith(e1,v3,v2,v1));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(3,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(v2[2]->getPtr())); CPPUNIT_ASSERT(v1[0]->getDirection()); CPPUNIT_ASSERT(v2[2]->getDirection());
  CPPUNIT_ASSERT(v1[2]->intresincEqCoarse(v2[0]->getPtr())); CPPUNIT_ASSERT(v1[2]->getDirection()); CPPUNIT_ASSERT(v2[0]->getDirection());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,v2[1]->getCurveLength(),1.e-5); // << not maximal precision because node switching 
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI/3.,v1[1]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT(v2[1]->getStartNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v2[1]->getEndNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v1[1]->getStartNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v1[1]->getEndNode()==e2->getStartNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 4 - OUT_BEFORE - OUT_BEFORE | same dir. - 0 °
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(-1.,0.,-0.5,0.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(0,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(0,(int)v2.size());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 5 - OUT_BEFORE - INSIDE | same dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(-0.5,0.,0.5,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresicEqualDirSensitive(v2[1]));
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 5 - OUT_BEFORE - INSIDE | same dir. - 90°
  e1=new EdgeLin(0.,0.,0.,1.); e2=new EdgeLin(0,-0.5,0.,0.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresicEqualDirSensitive(v2[1]));
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 5 - OUT_BEFORE - INSIDE | same dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(-0.5,-0.5,0.5,0.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresicEqualDirSensitive(v2[1]));
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 5 - OUT_BEFORE - INSIDE | opp dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(0.5,0.5,-0.5,-0.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(!v1[0]->intresicEqualDirSensitive(v2[0]) && v1[0]->intresicEqual(v2[0]) );
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 6 - OUT_BEFORE - OUT_AFTER | same dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(-0.5,0.,1.5,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection());
  CPPUNIT_ASSERT(v2[1]->intresincEqCoarse(e1) && v2[1]->getDirection());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode() && v2[1]->getEndNode()==v2[2]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 6 - OUT_BEFORE - OUT_AFTER | same dir. - 90°
  e1=new EdgeLin(0.,0.,0.,1.); e2=new EdgeLin(0.,-0.5,0.,1.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection());
  CPPUNIT_ASSERT(v2[1]->intresincEqCoarse(e1) && v2[1]->getDirection());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode() && v2[1]->getEndNode()==v2[2]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 6 - OUT_BEFORE - OUT_AFTER | same dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(-0.5,-0.5,1.5,1.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection());
  CPPUNIT_ASSERT(v2[1]->intresincEqCoarse(e1) && v2[1]->getDirection());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode() && v2[1]->getEndNode()==v2[2]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 6 - OUT_BEFORE - OUT_AFTER | opp dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(1.5,1.5,-0.5,-0.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection());
  CPPUNIT_ASSERT(v2[1]->intresincEqCoarse(e1) && !v2[1]->getDirection());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode() && v2[1]->getEndNode()==v2[2]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 7 - END - OUT_AFTER | same dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(1.,0.,1.5,0.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(0,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(0,(int)v2.size());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 7 - END - OUT_AFTER | opp dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(1.5,0.,1.,0.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(0,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(0,(int)v2.size());
  CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 8 - START - END | same dir. - 0°
  e1=new EdgeLin(0.,0.,0.7,0.); e2=new EdgeLin(0.,0.,0.7,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(2,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && v2[0]->getDirection());
  CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 8 - START - END | same dir. - 90°
  e1=new EdgeLin(0.,0.,0.,0.7); e2=new EdgeLin(0.,0.,0.,0.7);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(2,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && v2[0]->getDirection());
  CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 8 - START - END | same dir. - 45°
  e1=new EdgeLin(0.,0.,0.7,0.7); e2=new EdgeLin(0.,0.,0.7,0.7);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(2,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && v2[0]->getDirection());
  CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 8 - START - END | opp. dir. - 45°
  e1=new EdgeLin(0.,0.,0.7,0.7); e2=new EdgeLin(0.7,0.7,0.,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(2,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && !v2[0]->getDirection());
  CPPUNIT_ASSERT(e1->getStartNode()==e2->getEndNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 9 - OUT_BEFORE - START | same dir.
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(-0.5,0.,0.,0.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(0,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(0,(int)v2.size());
  CPPUNIT_ASSERT(e2->getEndNode()==e1->getStartNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 10 - START - OUT_AFTER | same dir. - 0°
  e1=new EdgeLin(0.,0.,0.7,0.); e2=new EdgeLin(0.,0.,1.,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && v2[0]->getDirection()); 
  CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(v2[1]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 10 - START - OUT_AFTER | same dir. - 90°
  e1=new EdgeLin(0.,0.,0.,0.7); e2=new EdgeLin(0.,0.,0.,1.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && v2[0]->getDirection()); 
  CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(v2[1]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 10 - START - OUT_AFTER | same dir. - 45°
  e1=new EdgeLin(0.,0.,0.7,0.7); e2=new EdgeLin(0.,0.,1.,1.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && v2[0]->getDirection()); 
  CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(v2[1]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 10 - START - OUT_AFTER | opp dir. - 45°
  e1=new EdgeLin(0.,0.,0.7,0.7); e2=new EdgeLin(1.,1.,0.,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[1]->intresincEqCoarse(e1) && !v2[1]->getDirection()); 
  CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==e2->getEndNode()); CPPUNIT_ASSERT(v2[1]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 11 - INSIDE - END | same dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(0.7,0.,1.,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresincEqCoarse(e2) && v1[1]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection()); 
  CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 11 - INSIDE - END | same dir. - 90°
  e1=new EdgeLin(0.,0.,0.,1.); e2=new EdgeLin(0.,0.7,0.,1.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresincEqCoarse(e2) && v1[1]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection()); 
  CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 11 - INSIDE - END | same dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(0.7,0.7,1.,1.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresincEqCoarse(e2) && v1[1]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection()); 
  CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 11 - INSIDE - END | opp dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(1.,1.,0.7,0.7);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(e1->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getStartNode());
  CPPUNIT_ASSERT(v1[1]->intresincEqCoarse(e2) && !v1[1]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection()); 
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 12 - OUT_BEFORE - END | same dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(-0.5,0.,1.,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[1]->intresincEqCoarse(e1) && v2[1]->getDirection());
  CPPUNIT_ASSERT(e2->getStartNode()==v2[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode()); CPPUNIT_ASSERT(e2->getEndNode()==v2[1]->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 12 - OUT_BEFORE - END | same dir. - 90°
  e1=new EdgeLin(0.,0.,0.,1.); e2=new EdgeLin(0.,-0.5,0.,1.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[1]->intresincEqCoarse(e1) && v2[1]->getDirection());
  CPPUNIT_ASSERT(e2->getStartNode()==v2[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode()); CPPUNIT_ASSERT(e2->getEndNode()==v2[1]->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 12 - OUT_BEFORE - END | same dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(-0.5,-0.5,1.,1.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[1]->intresincEqCoarse(e1) && v2[1]->getDirection());
  CPPUNIT_ASSERT(e2->getStartNode()==v2[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getEndNode()); CPPUNIT_ASSERT(e2->getEndNode()==v2[1]->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 12 - OUT_BEFORE - END | opp dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(1.,1.,-0.5,-0.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e1) && !v2[0]->getDirection());
  CPPUNIT_ASSERT(e2->getStartNode()==v2[0]->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==e2->getStartNode()); CPPUNIT_ASSERT(e2->getEndNode()==v2[1]->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 13 - START - INSIDE | same dir. - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(0.,0.,0.5,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e2) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection());
  CPPUNIT_ASSERT(e2->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 13 - START - INSIDE | same dir. - 90°
  e1=new EdgeLin(0.,0.,0.,1.); e2=new EdgeLin(0.,0.,0.,0.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e2) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection());
  CPPUNIT_ASSERT(e2->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 13 - START - INSIDE | same dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(0.,0.,0.5,0.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e2) && v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection());
  CPPUNIT_ASSERT(e2->getStartNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==e2->getStartNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 13 - START - INSIDE | opp dir. - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(0.5,0.5,0.,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e2) && !v1[0]->getDirection()); CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection());
  CPPUNIT_ASSERT(e2->getEndNode()==v1[0]->getStartNode()); CPPUNIT_ASSERT(e1->getStartNode()==e2->getEndNode()); CPPUNIT_ASSERT(e1->getEndNode()==v1[1]->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 14 - INSIDE - START | same dir.
  e1=buildArcOfCircle(center,radius,-M_PI,2.*M_PI); e2=buildArcOfCircle(center,radius,M_PI/3.,-M_PI);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,e1->getCurveLength(),1e-12); CPPUNIT_ASSERT_DOUBLES_EQUAL(4.*M_PI/3.,e2->getCurveLength(),1e-12);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresicEqual(v2[0]));
  CPPUNIT_ASSERT(v2[1]->getEndNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[1]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v2[1]->getStartNode()==e1->getEndNode());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2*M_PI/3.,v1[0]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI/3.,v2[0]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,v2[1]->getCurveLength(),1.e-12);
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 14 - INSIDE - START | opp dir.
  e1=buildArcOfCircle(center,radius,-M_PI,2.*M_PI); e2=buildArcOfCircle(center,radius,-M_PI,M_PI/3.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,e1->getCurveLength(),1e-12); CPPUNIT_ASSERT_DOUBLES_EQUAL(4.*M_PI/3.,e2->getCurveLength(),1e-12);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[1]->intresincEqCoarse(v2[1]->getPtr()) && !v2[1]->getDirection() && v1[1]->getDirection());
  CPPUNIT_ASSERT(v2[0]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[1]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v2[1]->getStartNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2*M_PI/3.,v1[0]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI/3.,v2[1]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,v2[0]->getCurveLength(),1.e-12);
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 15 - END - INSIDE | same dir.
  e1=buildArcOfCircle(center,radius,-M_PI,2.*M_PI); e2=buildArcOfCircle(center,radius,0.,-4.*M_PI/3);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,e1->getCurveLength(),1e-12); CPPUNIT_ASSERT_DOUBLES_EQUAL(4.*M_PI/3.,e2->getCurveLength(),1e-12);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresicEqual(v2[1]));
  CPPUNIT_ASSERT(v2[0]->getEndNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[1]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getStartNode()==e1->getEndNode());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI/3.,v1[0]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,v2[0]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.*M_PI/3.,v1[1]->getCurveLength(),1.e-12);
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 15 - END - INSIDE | opp dir.
  e1=buildArcOfCircle(center,radius,-M_PI,2.*M_PI); e2=buildArcOfCircle(center,radius,-4.*M_PI/3,0.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,e1->getCurveLength(),1e-12); CPPUNIT_ASSERT_DOUBLES_EQUAL(4.*M_PI/3.,e2->getCurveLength(),1e-12);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(v2[0]->getPtr()) && !v2[0]->getDirection() && v1[0]->getDirection());
  CPPUNIT_ASSERT(v2[0]->getEndNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode());
  CPPUNIT_ASSERT(v1[1]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v2[1]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI/3.,v1[0]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(M_PI,v2[1]->getCurveLength(),1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.*M_PI/3.,v1[1]->getCurveLength(),1.e-12);
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //
  ComposedEdge::Delete(&v1);
  ComposedEdge::Delete(&v2);
}

/*!
 * Here there is test of cases where between 2 edges intersects only in points not on edge.
 */
void QuadraticPlanarInterpTest::IntersectionPointOnlyUnitarySegSeg()
{
  // 0° - classical
  EdgeLin *e1=new EdgeLin(0.,0.,1.,0.);
  EdgeLin *e2=new EdgeLin(0.3,0.3,0.5,-0.3);
  ComposedEdge& v1=*(new ComposedEdge);
  ComposedEdge& v2=*(new ComposedEdge); MergePoints v3;
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4,(*v1[0]->getEndNode())[0],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,(*v1[0]->getEndNode())[1],ADMISSIBLE_ERROR);
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  // 90° - classical
  e1=new EdgeLin(0.,0.,0.,1.);
  e2=new EdgeLin(-0.3,0.3,0.3,0.5);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT(v1[0]->getEndNode()==v1[1]->getStartNode()); CPPUNIT_ASSERT(v2[0]->getEndNode()==v2[1]->getStartNode());
  CPPUNIT_ASSERT(e1->getStartNode()==v1.front()->getStartNode() && e1->getEndNode()==v1.back()->getEndNode());
  CPPUNIT_ASSERT(e2->getStartNode()==v2.front()->getStartNode() && e2->getEndNode()==v2.back()->getEndNode());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,(*v1[0]->getEndNode())[0],ADMISSIBLE_ERROR);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.4,(*v1[0]->getEndNode())[1],ADMISSIBLE_ERROR);
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 1 - 0°
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(0.,0.,0.,1.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v3.isStart1(0)); CPPUNIT_ASSERT(v3.isStart2(0));
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 1 - 90°
  e1=new EdgeLin(0.,0.,0.,1.); e2=new EdgeLin(0.,0.,1.,0.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v3.isStart1(0)); CPPUNIT_ASSERT(v3.isStart2(0));
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 1 - 45°
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(0.,0.,1.,-1.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v3.isStart1(0)); CPPUNIT_ASSERT(v3.isStart2(0));
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 2
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(1.,1.,1.,0.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v3.isEnd1(0)); CPPUNIT_ASSERT(v3.isEnd2(0));
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 3
  e1=new EdgeLin(0.,0.,1.,0.); e2=new EdgeLin(1.,0.,1.,1.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v3.isEnd1(0)); CPPUNIT_ASSERT(v3.isStart2(0));
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Test 4
  e1=new EdgeLin(0.,0.,1.,1.); e2=new EdgeLin(1.,-1.,0.,0.);
  CPPUNIT_ASSERT(!e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v3.isStart1(0)); CPPUNIT_ASSERT(v3.isEnd2(0));
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Intersection extremity of one edge and inside of other edge. 2 End.
  e1=new EdgeLin(0.,0.,1.,0.);
  e2=new EdgeLin(0.5,1.,0.5,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode() && v1[0]->getEndNode()==e2->getEndNode() && v1[1]->getStartNode()==e2->getEndNode() && v1[1]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getDirection() && v1[1]->getDirection());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Intersection extremity of one edge and inside of other edge. 2 Start.
  e1=new EdgeLin(0.,0.,1.,0.);
  e2=new EdgeLin(0.5,0.,0.5,1.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(2,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)v2.size());
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v2[0]->intresincEqCoarse(e2) && v2[0]->getDirection());
  CPPUNIT_ASSERT(v1[0]->getStartNode()==e1->getStartNode() && v1[0]->getEndNode()==e2->getStartNode() && v1[1]->getStartNode()==e2->getStartNode() && v1[1]->getEndNode()==e1->getEndNode());
  CPPUNIT_ASSERT(v1[0]->getDirection() && v1[1]->getDirection());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Intersection extremity of one edge and inside of other edge. 1 Start.
  e1=new EdgeLin(0.5,0.,0.5,1.);
  e2=new EdgeLin(0.,0.,1.,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection());
  CPPUNIT_ASSERT(v2[0]->getStartNode()==e2->getStartNode() && v2[0]->getEndNode()==e1->getStartNode() && v2[1]->getStartNode()==e1->getStartNode() && v2[1]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getDirection() && v2[1]->getDirection());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  //Intersection extremity of one edge and inside of other edge. 1 End.
  e1=new EdgeLin(0.5,1.,0.5,0.);
  e2=new EdgeLin(0.,0.,1.,0.);
  CPPUNIT_ASSERT(e1->intersectWith(e2,v3,v1,v2));
  CPPUNIT_ASSERT_EQUAL(1,(int)v1.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)v2.size());
  CPPUNIT_ASSERT_EQUAL(0,(int)v3.getNumberOfAssociations());
  CPPUNIT_ASSERT(v1[0]->intresincEqCoarse(e1) && v1[0]->getDirection());
  CPPUNIT_ASSERT(v2[0]->getStartNode()==e2->getStartNode() && v2[0]->getEndNode()==e1->getEndNode() && v2[1]->getStartNode()==e1->getEndNode() && v2[1]->getEndNode()==e2->getEndNode());
  CPPUNIT_ASSERT(v2[0]->getDirection() && v2[1]->getDirection());
  e2->decrRef(); e1->decrRef();
  v1.clear(); v2.clear(); v3.clear();
  ComposedEdge::Delete(&v2);
  ComposedEdge::Delete(&v1);
}

}
