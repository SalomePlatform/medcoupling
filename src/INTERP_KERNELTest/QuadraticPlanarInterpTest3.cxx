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
#include "InterpKernelGeo2DElementaryEdge.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"

#include <cmath>
#include <sstream>
#include <iostream>

using namespace INTERP_KERNEL;

namespace INTERP_TEST
{

void QuadraticPlanarInterpTest::checkInOutDetection()
{
  Node *n1=new Node(0.,0.);
  Node *n2=new Node(1.,0.);
  Node *n3=new Node(0.5,1.);
  EdgeLin *e1=new EdgeLin(n1,n2);
  EdgeLin *e2=new EdgeLin(n2,n3);
  EdgeLin *e3=new EdgeLin(n3,n1);
  ComposedEdge *tri=new ComposedEdge;
  tri->pushBack(e1); tri->pushBack(e2); tri->pushBack(e3);
  //
  Node *where=new Node(0.4,0.1);
  CPPUNIT_ASSERT(tri->isInOrOut(where)); where->decrRef();
  where=new Node(-0.1,1.);
  CPPUNIT_ASSERT(!tri->isInOrOut(where)); where->decrRef();
  where=new Node(0.6,-0.1);
  CPPUNIT_ASSERT(!tri->isInOrOut(where)); where->decrRef();
  //Clean-up
  n1->decrRef(); n2->decrRef(); n3->decrRef();
  ComposedEdge::Delete(tri);
}

/*!
 * Check Iterators mechanism.
 */
void QuadraticPlanarInterpTest::checkAssemblingBases1()
{
  Node *n1=new Node(0.,0.);
  Node *n2=new Node(0.1,0.); EdgeLin *e1_2=new EdgeLin(n1,n2);
  Node *n3=new Node(0.2,0.); EdgeLin *e2_3=new EdgeLin(n2,n3);
  Node *n4=new Node(0.3,0.); EdgeLin *e3_4=new EdgeLin(n3,n4);
  Node *n5=new Node(0.4,0.); EdgeLin *e4_5=new EdgeLin(n4,n5);
  Node *n6=new Node(0.5,0.); EdgeLin *e5_6=new EdgeLin(n5,n6);
  Node *n7=new Node(0.6,0.); EdgeLin *e6_7=new EdgeLin(n6,n7);
  Node *n8=new Node(0.7,0.); EdgeLin *e7_8=new EdgeLin(n7,n8);
  Node *n9=new Node(0.8,0.); EdgeLin *e8_9=new EdgeLin(n8,n9);
  Node *n10=new Node(0.9,0.); EdgeLin *e9_10=new EdgeLin(n9,n10);
  Node *n11=new Node(1.,0.); EdgeLin *e10_11=new EdgeLin(n10,n11);
  Node *n12=new Node(0.5,1.); EdgeLin *e11_12=new EdgeLin(n11,n12);
  EdgeLin *e12_1=new EdgeLin(n12,n1);
  //Only one level
  e1_2->incrRef(); e2_3->incrRef(); e3_4->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_7->incrRef(); 
  e7_8->incrRef(); e8_9->incrRef(); e9_10->incrRef(); e10_11->incrRef(); e11_12->incrRef(); e12_1->incrRef();
  ComposedEdge *c=new ComposedEdge;
  c->pushBack(e1_2); c->pushBack(e2_3); c->pushBack(e3_4); c->pushBack(e4_5); c->pushBack(e5_6); c->pushBack(e6_7);
  c->pushBack(e7_8); c->pushBack(e8_9); c->pushBack(e9_10); c->pushBack(e10_11); c->pushBack(e11_12); c->pushBack(e12_1);
  CPPUNIT_ASSERT_EQUAL(12,c->recursiveSize());
  IteratorOnComposedEdge it(c);
  CPPUNIT_ASSERT(it.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it.finished());
  it.next(); CPPUNIT_ASSERT(it.current()->getPtr()==e2_3); CPPUNIT_ASSERT(!it.finished());
  it.next(); it.next(); CPPUNIT_ASSERT(it.current()->getPtr()==e4_5); CPPUNIT_ASSERT(!it.finished());
  it.previousLoop(); CPPUNIT_ASSERT(it.current()->getPtr()==e3_4); CPPUNIT_ASSERT(!it.finished());
  it.previousLoop(); CPPUNIT_ASSERT(it.current()->getPtr()==e2_3); CPPUNIT_ASSERT(!it.finished());
  it.previousLoop(); CPPUNIT_ASSERT(it.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it.finished());
  it.previousLoop(); CPPUNIT_ASSERT(it.current()->getPtr()==e12_1); CPPUNIT_ASSERT(!it.finished());
  it.next(); CPPUNIT_ASSERT(it.finished());
  it.first(); CPPUNIT_ASSERT(it.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it.finished());
  it.previousLoop(); CPPUNIT_ASSERT(it.current()->getPtr()==e12_1); CPPUNIT_ASSERT(!it.finished());
  it.nextLoop(); CPPUNIT_ASSERT(it.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it.finished());
  it.last(); CPPUNIT_ASSERT(it.current()->getPtr()==e12_1); CPPUNIT_ASSERT(!it.finished());
  //Multi-Level
  ComposedEdge::Delete(c);
  //(e1_2, (e2_3,(e3_4, e4_5, e5_6, e6_7, (e7_8, e8_9 ), ( e9_10 , e10_11 ), e11_12 ),e12_1 ) )
  e1_2->incrRef(); e2_3->incrRef(); e3_4->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_7->incrRef(); 
  e7_8->incrRef(); e8_9->incrRef(); e9_10->incrRef(); e10_11->incrRef(); e11_12->incrRef(); e12_1->incrRef();
  ComposedEdge *c2_2_4=new ComposedEdge; c2_2_4->pushBack(e7_8); c2_2_4->pushBack(e8_9);
  ComposedEdge *c2_2_5=new ComposedEdge; c2_2_5->pushBack(e9_10); c2_2_5->pushBack(e10_11);
  ComposedEdge *c2_2=new ComposedEdge; c2_2->pushBack(e3_4); c2_2->pushBack(e4_5); c2_2->pushBack(e5_6);  c2_2->pushBack(e6_7); c2_2->pushBack(c2_2_4); c2_2->pushBack(c2_2_5); c2_2->pushBack(e11_12);
  ComposedEdge *c2=new ComposedEdge; c2->pushBack(e2_3); c2->pushBack(c2_2); c2->pushBack(e12_1);
  c=new ComposedEdge; c->pushBack(e1_2); c->pushBack(c2); CPPUNIT_ASSERT_EQUAL(12,c->recursiveSize());
  IteratorOnComposedEdge it2(c);
  CPPUNIT_ASSERT(it2.current()->getPtr()==e1_2);
  it2.next(); CPPUNIT_ASSERT(it2.current()->getPtr()==e2_3); CPPUNIT_ASSERT(!it2.finished());
  it2.next(); CPPUNIT_ASSERT(it2.current()->getPtr()==e3_4); CPPUNIT_ASSERT(!it2.finished());
  it2.next(); CPPUNIT_ASSERT(it2.current()->getPtr()==e4_5); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e3_4); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e2_3); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e12_1); CPPUNIT_ASSERT(!it2.finished());
  it2.next(); CPPUNIT_ASSERT(it2.finished());
  it2.first(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e12_1); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it2.finished());
  it2.last(); CPPUNIT_ASSERT(it2.current()->getPtr()==e12_1); CPPUNIT_ASSERT(!it2.finished());
  it2.first(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e2_3); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e3_4); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e4_5); CPPUNIT_ASSERT(!it2.finished());
  //  substitutions.
  /*it2.first(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it2.finished());
  ElementaryEdge *&tmp=it2.current(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_2); CPPUNIT_ASSERT(!it2.finished());
  ComposedEdge *c1=new ComposedEdge;  Node *n1_bis=new Node(0.,0.05); EdgeLin *e1_1bis=new EdgeLin(n1,n1_bis); EdgeLin *e1bis_2=new EdgeLin(n1_bis,n2); e1_1bis->incrRef(); e1bis_2->incrRef();
  c1->pushBack(e1_1bis); c1->pushBack(e1bis_2); delete tmp; tmp=(ElementaryEdge *)c1; CPPUNIT_ASSERT_EQUAL(13,c->recursiveSize());
  CPPUNIT_ASSERT(it2.current()->getPtr()==e1_1bis); CPPUNIT_ASSERT(!it2.finished());// here testing capability of Iterator.'current' method to deal with change of hierarchy.
  it2.next(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1bis_2); CPPUNIT_ASSERT(!it2.finished());
  it2.next(); CPPUNIT_ASSERT(it2.current()->getPtr()==e2_3); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1bis_2); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_1bis); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e12_1); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e11_12); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e10_11); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e9_10); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e8_9); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e7_8); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e6_7); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e5_6); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e4_5); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e3_4); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e2_3); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1bis_2); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_1bis); CPPUNIT_ASSERT(!it2.finished());
  it2.previousLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e12_1); CPPUNIT_ASSERT(!it2.finished());
  //go forward
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_1bis); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1bis_2); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e2_3); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e3_4); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e4_5); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e5_6); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e6_7); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e7_8); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e8_9); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e9_10); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e10_11); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e11_12); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e12_1); CPPUNIT_ASSERT(!it2.finished());
  it2.nextLoop(); CPPUNIT_ASSERT(it2.current()->getPtr()==e1_1bis); CPPUNIT_ASSERT(!it2.finished());*/
  ComposedEdge::SoftDelete(c2_2_4);
  ComposedEdge::SoftDelete(c2_2_5);
  ComposedEdge::SoftDelete(c2_2);
  ComposedEdge::SoftDelete(c2);
  ComposedEdge::Delete(c);
  //clean-up
  //e1_1bis->decrRef(); e1bis_2->decrRef();
  e1_2->decrRef(); e2_3->decrRef(); e3_4->decrRef(); e4_5->decrRef(); e5_6->decrRef(); e6_7->decrRef(); 
  e7_8->decrRef(); e8_9->decrRef(); e9_10->decrRef(); e10_11->decrRef(); e11_12->decrRef(); e12_1->decrRef(); 
  //n1_bis->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();
  n7->decrRef(); n8->decrRef(); n9->decrRef(); n10->decrRef(); n11->decrRef(); n12->decrRef();
}

/*!
 * Check splitting of 2 polygons. After this operation, all ElementaryEdge are either in/out/on.
 */
void QuadraticPlanarInterpTest::checkAssemblingBases2()
{
  //The "most" basic test1
  Node *n1=new Node(0.,0.);                Node *n4=new Node(0.,-0.3);   
  Node *n2=new Node(1.,0.);                Node *n5=new Node(1.,-0.3);
  Node *n3=new Node(0.5,1.);               Node *n6=new Node(0.5,0.7);
  EdgeLin *e1_2=new EdgeLin(n1,n2);        EdgeLin *e4_5=new EdgeLin(n4,n5);
  EdgeLin *e2_3=new EdgeLin(n2,n3);        EdgeLin *e5_6=new EdgeLin(n5,n6);
  EdgeLin *e3_1=new EdgeLin(n3,n1);        EdgeLin *e6_4=new EdgeLin(n6,n4);
  //
  e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_4->incrRef(); 
  QuadraticPolygon pol1; pol1.pushBack(e1_2); pol1.pushBack(e2_3); pol1.pushBack(e3_1);
  QuadraticPolygon pol2; pol2.pushBack(e4_5); pol2.pushBack(e5_6); pol2.pushBack(e6_4);
  QuadraticPolygon cpyPol1(pol1); int nbOfSplits=0;
  cpyPol1.SplitPolygonsEachOther(pol1,pol2,nbOfSplits);
  CPPUNIT_ASSERT_EQUAL(5,pol1.recursiveSize());
  CPPUNIT_ASSERT_EQUAL(5,pol2.recursiveSize());CPPUNIT_ASSERT_EQUAL(15,nbOfSplits);
  checkBasicsOfPolygons(pol1,pol2,true);
  CPPUNIT_ASSERT(pol2[1]->getEndNode()==pol1[1]->getEndNode());
  CPPUNIT_ASSERT(pol2[1]->getEndNode()->getLoc()==ON_1);
  CPPUNIT_ASSERT(pol2[3]->getEndNode()==pol1[0]->getEndNode());
  CPPUNIT_ASSERT(pol2[3]->getEndNode()->getLoc()==ON_1);
  cpyPol1.performLocatingOperation(pol2);
  ElementaryEdge *tmp=dynamic_cast<ElementaryEdge *>(pol2[0]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e4_5);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol2[1]); CPPUNIT_ASSERT(tmp);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol2[2]); CPPUNIT_ASSERT(tmp);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_IN_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol2[3]); CPPUNIT_ASSERT(tmp);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_IN_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol2[4]); CPPUNIT_ASSERT(tmp);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  //clean-up for test1
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e4_5->decrRef(); e5_6->decrRef(); e6_4->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();

  //Deeper test some extremities of pol2 are on edges of pol1.

  n1=new Node(0.,0.);                n4=new Node(1.5,-0.5);   
  n2=new Node(1.,0.);                n5=new Node(0.5,0.);
  n3=new Node(0.5,1.);               n6=new Node(0.75,0.5); Node *n7=new Node(2.,0.5);
  e1_2=new EdgeLin(n1,n2); e2_3=new EdgeLin(n2,n3); e3_1=new EdgeLin(n3,n1);
  EdgeLin *e5_4=new EdgeLin(n5,n4); EdgeLin *e4_7=new EdgeLin(n4,n7); EdgeLin *e7_6=new EdgeLin(n7,n6); EdgeLin *e6_5=new EdgeLin(n6,n5);
  //
  e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e5_4->incrRef(); e4_7->incrRef(); e7_6->incrRef(); e6_5->incrRef();
  QuadraticPolygon pol3; pol3.pushBack(e1_2); pol3.pushBack(e2_3); pol3.pushBack(e3_1);
  QuadraticPolygon pol4; pol4.pushBack(e5_4); pol4.pushBack(e4_7); pol4.pushBack(e7_6); pol4.pushBack(e6_5);
  QuadraticPolygon cpyPol3(pol3); nbOfSplits=0;
  cpyPol3.SplitPolygonsEachOther(pol3,pol4,nbOfSplits);
  CPPUNIT_ASSERT_EQUAL(5,pol3.recursiveSize());
  CPPUNIT_ASSERT_EQUAL(4,pol4.recursiveSize());CPPUNIT_ASSERT_EQUAL(16,nbOfSplits);
  checkBasicsOfPolygons(pol3,pol4,true);
  CPPUNIT_ASSERT(pol4[0]->getStartNode()==pol3[0]->getEndNode()); CPPUNIT_ASSERT(pol4[0]->getStartNode()==n5);
  CPPUNIT_ASSERT(n5->getLoc()==ON_LIM_1);
  CPPUNIT_ASSERT(pol4[2]->getEndNode()==pol3[2]->getEndNode()); CPPUNIT_ASSERT(pol4[2]->getEndNode()==n6);
  CPPUNIT_ASSERT(n6->getLoc()==ON_LIM_1);
  cpyPol3.performLocatingOperation(pol4);
  tmp=dynamic_cast<ElementaryEdge *>(pol4[1]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e4_7);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol4[3]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e6_5);
  tmp=dynamic_cast<ElementaryEdge *>(pol4[0]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e5_4);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol4[2]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e7_6);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol4[3]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e6_5);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_IN_1);
  //clean-up for test2
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e5_4->decrRef(); e4_7->decrRef(); e7_6->decrRef(); e6_5->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef(); n7->decrRef();

  //Test with one edge of pol2 is included in pol1.

  n1=new Node(0.,0.);                n4=new Node(-0.5,0.);   
  n2=new Node(1.,0.);                n5=new Node(0.,-1.);
  n3=new Node(0.5,1.);               n6=new Node(0.5,0.);
  e1_2=new EdgeLin(n1,n2); e2_3=new EdgeLin(n2,n3); e3_1=new EdgeLin(n3,n1);
  e4_5=new EdgeLin(n4,n5); e5_6=new EdgeLin(n5,n6); e6_4=new EdgeLin(n6,n4);
  e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_4->incrRef();
  QuadraticPolygon pol5; pol5.pushBack(e1_2); pol5.pushBack(e2_3); pol5.pushBack(e3_1);
  QuadraticPolygon pol6; pol6.pushBack(e4_5); pol6.pushBack(e5_6); pol6.pushBack(e6_4);
  QuadraticPolygon cpyPol5(pol5); nbOfSplits=0;
  cpyPol5.SplitPolygonsEachOther(pol5,pol6,nbOfSplits);
  CPPUNIT_ASSERT_EQUAL(4,pol5.recursiveSize());
  CPPUNIT_ASSERT_EQUAL(4,pol6.recursiveSize()); CPPUNIT_ASSERT_EQUAL(13,nbOfSplits);
  checkBasicsOfPolygons(pol5,pol6,false);
  CPPUNIT_ASSERT(pol6[2]->getStartNode()==pol5[0]->getEndNode()); CPPUNIT_ASSERT(pol6[2]->getStartNode()==n6);
  CPPUNIT_ASSERT(n6->getLoc()==ON_LIM_1);
  CPPUNIT_ASSERT(pol6[2]->getEndNode()==pol5[0]->getStartNode()); CPPUNIT_ASSERT(pol5[0]->getStartNode()==n1);
  CPPUNIT_ASSERT(n1->getLoc()==ON_LIM_1);
  cpyPol5.performLocatingOperation(pol6);
  tmp=dynamic_cast<ElementaryEdge *>(pol6[0]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e4_5);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol6[1]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e5_6);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol6[2]); CPPUNIT_ASSERT(tmp);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_ON_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol6[3]); CPPUNIT_ASSERT(tmp);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_OUT_1);
  //clean-up test3
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e4_5->decrRef(); e5_6->decrRef(); e6_4->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();

  //Test of full overlapped polygons.

  n1=new Node(0.,0.);                n4=new Node(0.,0.);   
  n2=new Node(1.,0.);                n5=new Node(1.,0.);
  n3=new Node(0.5,1.);               n6=new Node(0.5,1.);
  e1_2=new EdgeLin(n1,n2); e2_3=new EdgeLin(n2,n3); e3_1=new EdgeLin(n3,n1);
  e4_5=new EdgeLin(n4,n5); e5_6=new EdgeLin(n5,n6); e6_4=new EdgeLin(n6,n4);
  e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_4->incrRef();
  QuadraticPolygon pol7; pol7.pushBack(e1_2); pol7.pushBack(e2_3); pol7.pushBack(e3_1);
  QuadraticPolygon pol8; pol8.pushBack(e4_5); pol8.pushBack(e5_6); pol8.pushBack(e6_4);
  QuadraticPolygon cpyPol7(pol7); nbOfSplits=0;
  cpyPol7.SplitPolygonsEachOther(pol7,pol8,nbOfSplits);
  tmp=dynamic_cast<ElementaryEdge *>(pol8[0]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e1_2);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_ON_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol8[1]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e2_3);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_ON_1);
  tmp=dynamic_cast<ElementaryEdge *>(pol8[2]); CPPUNIT_ASSERT(tmp); CPPUNIT_ASSERT(tmp->getPtr()==e3_1);
  CPPUNIT_ASSERT(tmp->getLoc()==FULL_ON_1);
  //clean-up test4
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e4_5->decrRef(); e5_6->decrRef(); e6_4->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();
}

void QuadraticPlanarInterpTest::checkBasicsOfPolygons(QuadraticPolygon& pol1, QuadraticPolygon& pol2, bool checkDirection)
{
  IteratorOnComposedEdge it1(&pol1),it2(&pol2); it1.previousLoop(); it2.previousLoop();
  Node *nIter1=it1.current()->getEndNode(); Node *nIter2=it2.current()->getEndNode();
  for(it2.first();!it2.finished();it2.next())
    {
      CPPUNIT_ASSERT(nIter2==it2.current()->getStartNode());
      if(checkDirection)
        CPPUNIT_ASSERT(it2.current()->getDirection());
      nIter2=it2.current()->getEndNode();
    }
  for(it1.first();!it1.finished();it1.next())
    {
      CPPUNIT_ASSERT(nIter1==it1.current()->getStartNode());
      if(checkDirection)
        CPPUNIT_ASSERT(it1.current()->getDirection());
      nIter1=it1.current()->getEndNode();
    }
}

}
