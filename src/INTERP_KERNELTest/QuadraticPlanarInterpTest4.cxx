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
#include <iterator>

using namespace INTERP_KERNEL;

namespace INTERP_TEST
{

void QuadraticPlanarInterpTest::checkPolygonsIntersection1()
{
  //The "most" basic test1
  Node *n1=new Node(0.,0.);                Node *n4=new Node(0.,-0.3);   
  Node *n2=new Node(1.,0.);                Node *n5=new Node(1.,-0.3);
  Node *n3=new Node(0.5,1.);               Node *n6=new Node(0.5,0.7);
  EdgeLin *e1_2=new EdgeLin(n1,n2);        EdgeLin *e4_5=new EdgeLin(n4,n5);
  EdgeLin *e2_3=new EdgeLin(n2,n3);        EdgeLin *e5_6=new EdgeLin(n5,n6);
  EdgeLin *e3_1=new EdgeLin(n3,n1);        EdgeLin *e6_4=new EdgeLin(n6,n4);
  //
  std::vector<QuadraticPolygon *> result;
  for(int k=0;k<2;k++)
    for(int i=0;i<3;i++)
      {
        for(int j=0;j<1;j++)
          {
            e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_4->incrRef(); 
            QuadraticPolygon pol1; pol1.circularPermute(); pol1.pushBack(e1_2); pol1.pushBack(e2_3); pol1.pushBack(e3_1);
            for(int i1=0;i1<i;i1++) pol1.circularPermute(); if(k==1) pol1.reverse();
            QuadraticPolygon pol2; pol2.pushBack(e4_5); pol2.pushBack(e5_6); pol2.pushBack(e6_4);
            for(int j1=0;j1<j;j1++) pol2.circularPermute();
            result=pol1.intersectMySelfWith(pol2);
            CPPUNIT_ASSERT_EQUAL(1,(int)result.size()); checkBasicsOfPolygons(*result[0],*result[0],false);
            CPPUNIT_ASSERT_EQUAL(3,result[0]->recursiveSize());
            double tmp1=0.,tmp2=0.,tmp3=0.;
            pol1.intersectForPerimeter(pol2,tmp1,tmp2,tmp3);
            std::vector<double> v1,v2;
            std::vector<int> v3;
            pol1.intersectForPerimeterAdvanced(pol2,v1,v2);//no common edge
            pol1.intersectForPoint(pol2,v3);
            CPPUNIT_ASSERT_EQUAL(3,(int)v1.size());
            CPPUNIT_ASSERT_EQUAL(3,(int)v2.size());
            CPPUNIT_ASSERT_EQUAL(3,(int)v3.size());
            if(k==0)
              {
                CPPUNIT_ASSERT_EQUAL(2,v3[(3-i)%3]);
                CPPUNIT_ASSERT_EQUAL(0,v3[(4-i)%3]);
                CPPUNIT_ASSERT_EQUAL(0,v3[(5-i)%3]);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7,v1[(3-i)%3],1.e-14);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,v1[(4-i)%3],1.e-14);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,v1[(5-i)%3],1.e-14);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,v2[0],1.e-14);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.78262379212492639,v2[1],1.e-14);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(0.78262379212492639,v2[2],1.e-14);
              }
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7,tmp1,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5652475842498528,tmp2,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,tmp3,1.e-14);//no common edge
            delete result[0];
          }
      }
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
  for(int k=0;k<2;k++)
    for(int i=0;i<3;i++)
      {
        for(int j=0;j<4;j++)
          {
            e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e5_4->incrRef(); e4_7->incrRef(); e7_6->incrRef(); e6_5->incrRef();
            QuadraticPolygon pol3; pol3.pushBack(e1_2); pol3.pushBack(e2_3); pol3.pushBack(e3_1);
            for(int i1=0;i1<i;i1++) pol3.circularPermute(); if(k==1) pol3.reverse();
            QuadraticPolygon pol4; pol4.pushBack(e5_4); pol4.pushBack(e4_7); pol4.pushBack(e7_6); pol4.pushBack(e6_5);
            for(int j1=0;j1<j;j1++) pol4.circularPermute();
            result=pol3.intersectMySelfWith(pol4);
            CPPUNIT_ASSERT_EQUAL(1,(int)result.size()); checkBasicsOfPolygons(*result[0],*result[0],false);
            CPPUNIT_ASSERT_EQUAL(3,result[0]->recursiveSize());
            delete result[0];          
          }
      }
  //clean-up for test2
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e5_4->decrRef(); e4_7->decrRef(); e7_6->decrRef(); e6_5->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef(); n7->decrRef();

  //Test with one edge of pol2 is included in pol1.

  n1=new Node(0.,0.);                n4=new Node(-0.5,0.);   
  n2=new Node(1.,0.);                n5=new Node(0.,-1.);
  n3=new Node(0.5,1.);               n6=new Node(0.5,0.);
  e1_2=new EdgeLin(n1,n2); e2_3=new EdgeLin(n2,n3); e3_1=new EdgeLin(n3,n1);
  e4_5=new EdgeLin(n4,n5); e5_6=new EdgeLin(n5,n6); e6_4=new EdgeLin(n6,n4);
  for(int k=0;k<2;k++)
    for(int i=0;i<3;i++)
      {
        for(int j=0;j<3;j++)
          {
            e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_4->incrRef();
            QuadraticPolygon pol5; pol5.pushBack(e1_2); pol5.pushBack(e2_3); pol5.pushBack(e3_1);
            for(int i1=0;i1<i;i1++) pol5.circularPermute(); if(k==1) pol5.reverse();
            QuadraticPolygon pol6; pol6.pushBack(e4_5); pol6.pushBack(e5_6); pol6.pushBack(e6_4);
            for(int j1=0;j1<j;j1++) pol6.circularPermute();
            result=pol5.intersectMySelfWith(pol6);
            CPPUNIT_ASSERT_EQUAL(0,(int)result.size());
          }
      }
  //clean-up test3
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e4_5->decrRef(); e5_6->decrRef(); e6_4->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();

  //Test of full overlapped polygons.

  n1=new Node(0.,0.);                n4=new Node(0.,0.);   
  n2=new Node(1.,0.);                n5=new Node(1.,0.);
  n3=new Node(0.5,1.);               n6=new Node(0.5,1.);
  e1_2=new EdgeLin(n1,n2); e2_3=new EdgeLin(n2,n3); e3_1=new EdgeLin(n3,n1);
  e4_5=new EdgeLin(n4,n5); e5_6=new EdgeLin(n5,n6); e6_4=new EdgeLin(n6,n4);
  for(int k=0;k<2;k++)
    for(int i=0;i<3;i++)
      {
        for(int j=0;j<3;j++)
          {
            e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_4->incrRef();
            QuadraticPolygon pol7; pol7.pushBack(e1_2); pol7.pushBack(e2_3); pol7.pushBack(e3_1);
            for(int i1=0;i1<i;i1++) pol7.circularPermute(); if(k==1) pol7.reverse();
            QuadraticPolygon pol8; pol8.pushBack(e4_5); pol8.pushBack(e5_6); pol8.pushBack(e6_4);
            for(int j1=0;j1<j;j1++) pol8.circularPermute();
            result=pol7.intersectMySelfWith(pol8);
            CPPUNIT_ASSERT_EQUAL(1,(int)result.size()); checkBasicsOfPolygons(*result[0],*result[0],false);
            CPPUNIT_ASSERT_EQUAL(3,result[0]->recursiveSize());
            delete result[0];
            double tmp1=0.,tmp2=0.,tmp3=0.;
            pol7.intersectForPerimeter(pol8,tmp1,tmp2,tmp3);
            std::vector<double> v1,v2;
            pol7.intersectForPerimeterAdvanced(pol8,v1,v2);//only common edges.
            CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2360679774997898,v1[0]+v1[1]+v1[2],1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2360679774997898,v2[0]+v2[1]+v2[2],1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,tmp1,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,tmp2,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2360679774997898,tmp3,1.e-14);
          }
      }
  //clean-up test4
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e4_5->decrRef(); e5_6->decrRef(); e6_4->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();

  //Test of closing process
  
  n1=new Node(0.,0.);                n4=new Node(0.539,-0.266);   
  n2=new Node(1.,0.);                n5=new Node(1.039,0.6);
  n3=new Node(0.5,1.);               n6=new Node(-0.077,0.667);
  e1_2=new EdgeLin(n1,n2); e2_3=new EdgeLin(n2,n3); e3_1=new EdgeLin(n3,n1);
  e4_5=new EdgeLin(n4,n5); e5_6=new EdgeLin(n5,n6); e6_4=new EdgeLin(n6,n4);
  for(int k=0;k<2;k++)
    for(int i=0;i<3;i++)
      {
        for(int j=0;j<3;j++)
          {
            e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_4->incrRef();
            QuadraticPolygon pol9; pol9.pushBack(e1_2); pol9.pushBack(e2_3); pol9.pushBack(e3_1);
            for(int i1=0;i1<i;i1++) pol9.circularPermute(); if(k==1) pol9.reverse();
            QuadraticPolygon pol10; pol10.pushBack(e5_6); pol10.pushBack(e6_4); pol10.pushBack(e4_5);
            for(int j1=0;j1<j;j1++) pol10.circularPermute();
            result=pol9.intersectMySelfWith(pol10);
            CPPUNIT_ASSERT_EQUAL(1,(int)result.size()); checkBasicsOfPolygons(*result[0],*result[0],false);
            CPPUNIT_ASSERT_EQUAL(6,result[0]->recursiveSize());
            delete result[0];
          }
      }
  //clean-up test5
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e4_5->decrRef(); e5_6->decrRef(); e6_4->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();

  // Full in case

  n1=new Node(0.,0.);                n4=new Node(0.3,0.1);   
  n2=new Node(1.,0.);                n5=new Node(0.7,0.1);
  n3=new Node(0.5,1.);               n6=new Node(0.5,0.7);
  e1_2=new EdgeLin(n1,n2); e2_3=new EdgeLin(n2,n3); e3_1=new EdgeLin(n3,n1);
  e4_5=new EdgeLin(n4,n5); e5_6=new EdgeLin(n5,n6); e6_4=new EdgeLin(n6,n4);
  for(int k=0;k<2;k++)
    for(int i=0;i<3;i++)
      {
        for(int j=0;j<3;j++)
          {
            e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_4->incrRef();
            QuadraticPolygon pol11; pol11.pushBack(e1_2); pol11.pushBack(e2_3); pol11.pushBack(e3_1);
            for(int i1=0;i1<i;i1++) pol11.circularPermute(); if(k==1) pol11.reverse();
            QuadraticPolygon pol12; pol12.pushBack(e5_6); pol12.pushBack(e6_4); pol12.pushBack(e4_5);
            for(int j1=0;j1<j;j1++) pol12.circularPermute();
            result=pol11.intersectMySelfWith(pol12);
            CPPUNIT_ASSERT_EQUAL(1,(int)result.size()); checkBasicsOfPolygons(*result[0],*result[0],false);
            CPPUNIT_ASSERT_EQUAL(3,result[0]->recursiveSize());
            delete result[0];
          }
      }
  //clean-up test6
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e4_5->decrRef(); e5_6->decrRef(); e6_4->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();

  // Full out case

  n1=new Node(0.,0.);                n4=new Node(-2,0.);   
  n2=new Node(1.,0.);                n5=new Node(-1.,0.);
  n3=new Node(0.5,1.);               n6=new Node(-1.5,1.);
  e1_2=new EdgeLin(n1,n2); e2_3=new EdgeLin(n2,n3); e3_1=new EdgeLin(n3,n1);
  e4_5=new EdgeLin(n4,n5); e5_6=new EdgeLin(n5,n6); e6_4=new EdgeLin(n6,n4);
  for(int k=0;k<2;k++)
    for(int i=0;i<3;i++)
      {
        for(int j=0;j<3;j++)
          {
            e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef(); e4_5->incrRef(); e5_6->incrRef(); e6_4->incrRef();
            QuadraticPolygon pol13; pol13.pushBack(e1_2); pol13.pushBack(e2_3); pol13.pushBack(e3_1);
            for(int i1=0;i1<i;i1++) pol13.circularPermute(); if(k==1) pol13.reverse();
            QuadraticPolygon pol14; pol14.pushBack(e5_6); pol14.pushBack(e6_4); pol14.pushBack(e4_5);
            for(int j1=0;j1<j;j1++) pol14.circularPermute();
            result=pol13.intersectMySelfWith(pol14);
            CPPUNIT_ASSERT_EQUAL(0,(int)result.size());
          }
      }
  //clean-up test7
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef(); e4_5->decrRef(); e5_6->decrRef(); e6_4->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();

  //Multi polygons
  
  n1=new Node(0.,0.);
  n2=new Node(1.,0.);
  n3=new Node(1.,1.);
  n4=new Node(0.,1.);
  //
  n5=new Node(0.2,0.7);
  n6=new Node(0.4,0.7);
  n7=new Node(0.4,1.3);
  Node *n8=new Node(0.6,1.3);
  Node *n9=new Node(0.6,0.7);
  Node *n10=new Node(0.9,0.7);
  Node *n11=new Node(0.9,2.);
  Node *n12=new Node(0.2,2.);
  //
  e1_2=new EdgeLin(n1,n2); e2_3=new EdgeLin(n2,n3); Edge *e3_4=new EdgeLin(n3,n4); Edge *e4_1=new EdgeLin(n4,n1);
  e5_6=new EdgeLin(n5,n6); Edge *e6_7=new EdgeLin(n6,n7); Edge *e7_8=new EdgeLin(n7,n8); Edge *e8_9=new EdgeLin(n8,n9); Edge *e9_10=new EdgeLin(n9,n10); Edge *e10_11=new EdgeLin(n10,n11);
  Edge *e11_12=new EdgeLin(n11,n12); Edge *e12_1=new EdgeLin(n12,n5);
  //
  for(int k=0;k<2;k++)
    for(int i=0;i<4;i++)
      {
        for(int j=0;j<8;j++)
          {
            e1_2->incrRef(); e2_3->incrRef(); e3_4->incrRef(); e4_1->incrRef(); e5_6->incrRef(); e6_7->incrRef(); e7_8->incrRef(); e8_9->incrRef(); e9_10->incrRef(); e10_11->incrRef(); e11_12->incrRef(); e12_1->incrRef();
            QuadraticPolygon pol15; pol15.pushBack(e1_2); pol15.pushBack(e2_3); pol15.pushBack(e3_4); pol15.pushBack(e4_1);
            for(int i1=0;i1<i;i1++) pol15.circularPermute(); if(k==1) pol15.reverse();
            QuadraticPolygon pol16; pol16.pushBack(e5_6); pol16.pushBack(e6_7); pol16.pushBack(e7_8); pol16.pushBack(e8_9); pol16.pushBack(e9_10); pol16.pushBack(e10_11); pol16.pushBack(e11_12); pol16.pushBack(e12_1);
            for(int j1=0;j1<j;j1++) pol16.circularPermute();
            result=pol15.intersectMySelfWith(pol16);
            CPPUNIT_ASSERT_EQUAL(2,(int)result.size());
            checkBasicsOfPolygons(*result[0],*result[1],false);
            CPPUNIT_ASSERT_EQUAL(4,result[0]->recursiveSize()); CPPUNIT_ASSERT_EQUAL(4,result[1]->recursiveSize());
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.15,result[0]->getArea()+result[1]->getArea(),1e-10);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.03,fabs(result[0]->getArea()-result[1]->getArea()),1e-10);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.15,pol15.intersectWith(pol16),1e-10);
            delete result[0]; delete result[1];
          }
      }
  //clean-up test8
  e1_2->decrRef(); e2_3->decrRef(); e3_4->decrRef(); e4_1->decrRef(); e5_6->decrRef(); e6_7->decrRef(); e7_8->decrRef(); e8_9->decrRef(); e9_10->decrRef(); e10_11->decrRef(); e11_12->decrRef(); e12_1->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef(); n7->decrRef(); n8->decrRef(); n9->decrRef(); n10->decrRef(); n11->decrRef(); n12->decrRef();
}

/*!
 * Testing case where a polygon pol1 is included in an onother polygon pol2.
 */
void QuadraticPlanarInterpTest::checkPolygonsIntersection2()
{
  Node *n1=new Node(0.,0.);          Node *n4=new Node(0.2,0.2);
  Node *n2=new Node(1.,0.);          Node *n5=new Node(0.8,0.2);
  Node *n3=new Node(0.5,1.);         Node *n6=new Node(0.5,0.8);
  Edge *e1_2=new EdgeLin(n1,n2);     Edge *e4_5=new EdgeLin(n4,n5);
  Edge *e2_3=new EdgeLin(n2,n3);     Edge *e5_6=new EdgeLin(n5,n6);
  Edge *e3_1=new EdgeLin(n3,n1);     Edge *e6_4=new EdgeLin(n6,n4);
  //
  QuadraticPolygon pol1; pol1.pushBack(e1_2); pol1.pushBack(e2_3); pol1.pushBack(e3_1);
  QuadraticPolygon pol2; pol2.pushBack(e4_5); pol2.pushBack(e5_6); pol2.pushBack(e6_4);
  std::vector<QuadraticPolygon *> result=pol1.intersectMySelfWith(pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)result.size());
  CPPUNIT_ASSERT_EQUAL(3,result[0]->recursiveSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.18,result[0]->getArea(),1e-10);
  delete result[0];
  result.clear();
  pol1.initLocations();
  pol2.initLocations();
  result=pol2.intersectMySelfWith(pol1);
  CPPUNIT_ASSERT_EQUAL(1,(int)result.size());
  CPPUNIT_ASSERT_EQUAL(3,result[0]->recursiveSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.18,result[0]->getArea(),1e-10);
  delete result[0];
  //clean-up
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();
}

void QuadraticPlanarInterpTest::checkAreasCalculations()
{
  Node *n1=new Node(0.,0.);
  Node *n2=new Node(1.,0.);
  Node *n3=new Node(0.5,1.);
  Edge *e1_2=new EdgeLin(n1,n2);
  Edge *e2_3=new EdgeLin(n2,n3);
  Edge *e3_1=new EdgeLin(n3,n1);
  //
  e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef();
  QuadraticPolygon pol1; pol1.pushBack(e1_2); pol1.pushBack(e2_3); pol1.pushBack(e3_1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,pol1.getArea(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2360679774997898,pol1.getPerimeter(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.61803398874989479,pol1.getHydraulicDiameter(),1e-10);
  pol1.reverse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5,pol1.getArea(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2360679774997898,pol1.getPerimeter(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.61803398874989479,pol1.getHydraulicDiameter(),1e-10);
  //clean-up
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef();

  //case 2

  n1=new Node(0.,0.);
  n2=new Node(1.,0.);
  Node *n3m=new Node(1.5,0.5);
  n3=new Node(1.,1.);
  Node *n4=new Node(0.,1.);
  e1_2=new EdgeLin(n1,n2);
  e2_3=new EdgeArcCircle(n2,n3m,n3);
  Edge *e3_4=new EdgeLin(n3,n4);
  Edge *e4_1=new EdgeLin(n4,n1);
  //
  for(int k=0;k<8;k++)
    {
      n2->setNewCoords(cos(k*M_PI/4),sin(k*M_PI/4));
      n3->setNewCoords(sqrt(2.)*cos((k+1)*M_PI/4),sqrt(2.)*sin((k+1)*M_PI/4));
      n3m->setNewCoords(1.5811388300841898*cos(0.3217505543966423+k*M_PI/4),1.5811388300841898*sin(0.3217505543966423+k*M_PI/4));
      n4->setNewCoords(cos(k*M_PI/4+M_PI/2),sin(k*M_PI/4+M_PI/2));
      e1_2->update(n3m); e2_3->update(n3m); e3_4->update(n3m); e4_1->update(n3m);
      e1_2->incrRef(); e2_3->incrRef(); e3_4->incrRef(); e4_1->incrRef();
      QuadraticPolygon pol2; pol2.pushBack(e1_2); pol2.pushBack(e2_3); pol2.pushBack(e3_4); pol2.pushBack(e4_1);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.3926990816987241,pol2.getArea(),1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(4.5707963267948966,pol2.getPerimeter(),1e-6);
      pol2.reverse();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.3926990816987241,pol2.getArea(),1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(4.5707963267948966,pol2.getPerimeter(),1e-6);
    }
  //clean-up case2
  e1_2->decrRef(); e2_3->decrRef(); e3_4->decrRef(); e4_1->decrRef(); 
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n3m->decrRef(); n4->decrRef();
  
  //case 3

  const double radius1=0.7;
  const double radius2=0.9;
  n1=new Node(1.+radius1*cos(-2.*M_PI/3.),1.+radius1*sin(-2.*M_PI/3.));
  n2=new Node(1.+radius1*cos(-M_PI/3.),1.+radius1*sin(-M_PI/3.));
  Node *n2m=new Node(1.+radius1*cos(M_PI/2.),1.+radius1*sin(M_PI/2.));
  n3=new Node(1.+radius2*cos(-M_PI/3.),1.+radius2*sin(-M_PI/3.));
  n3m=new Node(1.+radius2*cos(M_PI/2.),1.+radius2*sin(M_PI/2.));
  n4=new Node(1.+radius2*cos(-2.*M_PI/3.),1.+radius2*sin(-2.*M_PI/3.));
  e1_2=new EdgeArcCircle(n1,n2m,n2);
  e2_3=new EdgeLin(n2,n3);
  e3_4=new EdgeArcCircle(n3,n3m,n4);
  e4_1=new EdgeLin(n4,n1);
  //
  e1_2->incrRef(); e2_3->incrRef(); e3_4->incrRef(); e4_1->incrRef();
  QuadraticPolygon pol3; pol3.pushBack(e1_2); pol3.pushBack(e2_3); pol3.pushBack(e3_4); pol3.pushBack(e4_1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.83775804095727857,pol3.getArea(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.7775804095727832,pol3.getPerimeter(),1e-10);
  pol3.reverse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.83775804095727857,pol3.getArea(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.7775804095727832,pol3.getPerimeter(),1e-10);
  //clean-up case3
  e1_2->decrRef(); e2_3->decrRef(); e3_4->decrRef(); e4_1->decrRef(); 
  n1->decrRef(); n2->decrRef(); n2m->decrRef(); n3->decrRef(); n3m->decrRef(); n4->decrRef();
}

void QuadraticPlanarInterpTest::checkBarycenterCalculations()
{
  Node *n1=new Node(3.,7.);
  Node *n2=new Node(5.,7.);
  Node *n3=new Node(4.,8.);
  Edge *e1_2=new EdgeLin(n1,n2);
  Edge *e2_3=new EdgeLin(n2,n3);
  Edge *e3_1=new EdgeLin(n3,n1);
  //
  double bary[2];
  e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef();
  QuadraticPolygon pol1; pol1.pushBack(e1_2); pol1.pushBack(e2_3); pol1.pushBack(e3_1);
  bary[0]=0.; bary[1]=0.;
  e1_2->getBarycenterOfZone(bary);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-56.,bary[0],1.e-10);
  bary[0]=0.; bary[1]=0.;
  e2_3->getBarycenterOfZone(bary);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(33.66666666666667,bary[0],1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(28.16666666666667,bary[1],1.e-10);
  bary[0]=0.; bary[1]=0.;
  e3_1->getBarycenterOfZone(bary);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(26.333333333333336,bary[0],1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(28.1666666666667,bary[1],1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,pol1.getArea(),1e-10);
  pol1.getBarycenter(bary);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.,bary[0],1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.333333333333333,bary[1],1.e-10);
  //
  e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef();
  QuadraticPolygon pol4; pol4.pushBack(e3_1,false); pol4.pushBack(e2_3,false); pol4.pushBack(e1_2,false);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.,pol4.getArea(),1e-10);
  pol4.getBarycenter(bary);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.,bary[0],1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.333333333333333,bary[1],1.e-10);
  //clean-up
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef();
  //Inverting polygon
  n1=new Node(3.,7.);
  n2=new Node(5.,7.);
  n3=new Node(4.,8.);
  e1_2=new EdgeLin(n1,n3);
  e2_3=new EdgeLin(n3,n2);
  e3_1=new EdgeLin(n2,n1);
  e1_2->incrRef(); e2_3->incrRef(); e3_1->incrRef();
  QuadraticPolygon pol3; pol3.pushBack(e1_2); pol3.pushBack(e2_3); pol3.pushBack(e3_1);
  bary[0]=0.; bary[1]=0.;
  pol3.getBarycenter(bary);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.,pol3.getArea(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.,bary[0],1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.333333333333333,bary[1],1.e-10);
  //clean-up
  e1_2->decrRef(); e2_3->decrRef(); e3_1->decrRef();
  n1->decrRef(); n2->decrRef(); n3->decrRef();
  //
  double center[2]={3.,7.};
  e1_2=buildArcOfCircle(center,4.,M_PI/3.,4.*M_PI/3.);
  bary[0]=0.; bary[1]=0.;
  e1_2->getBarycenterOfZone(bary);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(131.685410765053,bary[0],1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(303.262521934362,bary[1],1.e-10);
  n1=new Node(0.99999999999999822,3.5358983848622465);
  n2=new Node(5.,10.4641016151377544);
  Edge *e2_1=new EdgeLin(n1,n2);
  //
  e1_2->incrRef(); e2_1->incrRef();
  QuadraticPolygon pol2; pol2.pushBack(e1_2); pol2.pushBack(e2_1);
  pol2.getBarycenter(bary);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(25.132741228718345,pol2.getArea(),1e-10);
  //4*radius/(3.*pi)
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5297896122085546,bary[0],1.e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.8488263631567756,bary[1],1.e-10);
  //clean-up
  e1_2->decrRef(); e2_1->decrRef();
  n1->decrRef(); n2->decrRef();
}

/*!
 * Testing user interface high level function.
 */
void QuadraticPlanarInterpTest::checkHighLevelFunctionTest1()
{
  QuadraticPlanarPrecision::setPrecision(1e-12);
  QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-9);
  double coords[]={
    8.8334591186000004, 5.0999999999999996,
    7.1014083111000001, 6.0999999999999996,
    7.8334591186000004, 6.8320508074999999,
    7.9674337149000003, 5.5999999999999996,
    7.4192455562999999, 6.5142135623000001,
    8.3334591186000004, 5.9660254036999998
  };
  std::vector<Node *> nodes;
  nodes.push_back(new Node(coords));
  nodes.push_back(new Node(coords+2));
  nodes.push_back(new Node(coords+4));
  nodes.push_back(new Node(coords+6));
  nodes.push_back(new Node(coords+8));
  nodes.push_back(new Node(coords+10));
  QuadraticPolygon *pol=QuadraticPolygon::BuildArcCirclePolygon(nodes);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.04719755,pol->getArea(),1e-5);
  CPPUNIT_ASSERT_EQUAL(3,pol->size());
  ElementaryEdge *e0=dynamic_cast<ElementaryEdge *>((*pol)[0]);
  ElementaryEdge *e1=dynamic_cast<ElementaryEdge *>((*pol)[1]);
  ElementaryEdge *e2=dynamic_cast<ElementaryEdge *>((*pol)[0]);
  CPPUNIT_ASSERT(e0); CPPUNIT_ASSERT(e1); CPPUNIT_ASSERT(e2);
  CPPUNIT_ASSERT(dynamic_cast<EdgeLin *>(e0->getPtr()));//<- testing detection of colinearity
  CPPUNIT_ASSERT(dynamic_cast<EdgeArcCircle *>(e1->getPtr()));
  CPPUNIT_ASSERT(dynamic_cast<EdgeLin *>(e2->getPtr()));//<- testing detection of colinearity
  nodes.clear();
  delete pol;
  nodes.push_back(new Node(coords));
  nodes.push_back(new Node(coords+4));
  nodes.push_back(new Node(coords+2));
  nodes.push_back(new Node(coords+10));
  nodes.push_back(new Node(coords+8));
  nodes.push_back(new Node(coords+6));
  pol=QuadraticPolygon::BuildArcCirclePolygon(nodes);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.04719755,pol->getArea(),1e-5);
  CPPUNIT_ASSERT_EQUAL(3,pol->size());
  e0=dynamic_cast<ElementaryEdge *>((*pol)[0]);
  e1=dynamic_cast<ElementaryEdge *>((*pol)[1]);
  e2=dynamic_cast<ElementaryEdge *>((*pol)[0]);
  CPPUNIT_ASSERT(e0); CPPUNIT_ASSERT(e1); CPPUNIT_ASSERT(e2);
  CPPUNIT_ASSERT(dynamic_cast<EdgeLin *>(e0->getPtr()));//<- testing detection of colinearity
  CPPUNIT_ASSERT(dynamic_cast<EdgeArcCircle *>(e1->getPtr()));
  CPPUNIT_ASSERT(dynamic_cast<EdgeLin *>(e2->getPtr()));//<- testing detection of colinearity
  delete pol;
  const double coords2[]={
    0.,0.,
    1.5,0.,
    1.5,1.,
    0.,1.
  };
  nodes.clear();
  nodes.push_back(new Node(coords2));
  nodes.push_back(new Node(coords2+2));
  nodes.push_back(new Node(coords2+4));
  nodes.push_back(new Node(coords2+6));
  pol=QuadraticPolygon::BuildLinearPolygon(nodes);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5,pol->getArea(),1e-12);
  double tmp[2],tmp2;
  pol->getBarycenter(tmp,tmp2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.75,tmp[0],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,tmp[1],1e-12);
  delete pol;
  const double coords3[]={
    1.0999999999000001, -1.9052558882999999,
    1.9052558881999999, -1.0999999999000001,
    1.7320508075000001, -0.99999999989999999,
    0.99999999989999999, -1.7320508075000001,
    1.5556349186, -1.5556349185,
    1.8186533478, -1.0499999999,
    1.4142135623000001, -1.4142135623000001,
    1.0499999999, -1.8186533479
  };
  nodes.clear();
  nodes.push_back(new Node(coords3));
  nodes.push_back(new Node(coords3+2));
  nodes.push_back(new Node(coords3+4));
  nodes.push_back(new Node(coords3+6));
  nodes.push_back(new Node(coords3+8));
  nodes.push_back(new Node(coords3+10));
  nodes.push_back(new Node(coords3+12));
  nodes.push_back(new Node(coords3+14));
  pol=QuadraticPolygon::BuildArcCirclePolygon(nodes);
  pol->getBarycenter(tmp,tmp2);
  delete pol;
  QuadraticPlanarPrecision::setPrecision(1e-14);
}

void QuadraticPlanarInterpTest::check1DInterpLin()
{
  QuadraticPlanarPrecision::setPrecision(1e-7);
  QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-9);
  const int NB_OF_CELL_AXIAL_1=30;
  static const double Z_VALS_1[NB_OF_CELL_AXIAL_1+1]=
    { -0.1550 , -0.1356, -0.1162, -0.0969, -0.0775 ,-0.0581, -0.0387, -0.0194,  0.0000 , 0.0500, 
      0.1000 , 0.1500 , 0.2000 , 0.2500,  0.3000,  0.3500,  0.4000,  0.4500,  0.5000,  0.5500, 
      0.6000,  0.6500,  0.7000,  0.7194,  0.7388,  0.7581,  0.7775,  0.7969,  0.8163,  0.8356, 
      0.8550};
  std::vector<double> zLev1(Z_VALS_1,Z_VALS_1+NB_OF_CELL_AXIAL_1+1);

  const int NB_OF_CELL_AXIAL_2=46;
  static const double Z_VALS_2[NB_OF_CELL_AXIAL_2+1]=
    { -0.3050 ,-0.2863,-0.2675,-0.2488,-0.2300,-0.2113,-0.1925,-0.1738,-0.1550,-0.1356  
      , -0.1162,-0.0969,-0.0775,-0.0581,-0.0387,-0.0194,0.0000, 0.0500, 0.1 ,0.15 
      ,  0.20,  0.25, 0.30, 0.350 ,0.40 ,0.450 ,0.500 , 0.550, 0.600 ,0.650 ,0.700
      , 0.7194 ,0.7388 ,0.7581 ,0.7775 ,0.7969 ,0.8163 ,0.8356, 0.8550
      ,  0.8738 ,0.8925 ,0.9113 ,0.9300 ,0.9488 ,0.9675 ,0.9863, 1.0050};
  std::vector<double> zLev2(Z_VALS_2,Z_VALS_2+NB_OF_CELL_AXIAL_2+1);
  std::map<int,std::map<int,double> > m;
  Edge::Interpolate1DLin(zLev1,zLev2,m);
  CPPUNIT_ASSERT_EQUAL(30,(int)m.size());
  double ret=0;
  for(int i=0;i<30;i++)
    {
      CPPUNIT_ASSERT_EQUAL(1,(int)m[i].size());
      CPPUNIT_ASSERT(m[i][8+i] > 0.15);
      ret+=m[i][8+i];
    }
  CPPUNIT_ASSERT_DOUBLES_EQUAL(ret,30.,1e-12);
  //
  m.clear();
  const int NB_OF_CELL_AXIAL_3=13;
  static const double Z_VALS_3[NB_OF_CELL_AXIAL_3+1]={
    0.,0.01,0.05,0.10,0.15,0.20,0.25,0.30,
    0.35,0.40,0.45,0.50,0.55,0.60 };
  std::vector<double> zLev3(Z_VALS_3,Z_VALS_3+NB_OF_CELL_AXIAL_3+1);
  Edge::Interpolate1DLin(zLev3,zLev1,m);
  CPPUNIT_ASSERT_EQUAL(13,(int)m.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,m[0][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,m[1][8],1e-12);
  for(int i=0;i<11;i++)
    {
      CPPUNIT_ASSERT_EQUAL(1,(int)m[i+2].size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,m[i+2][i+9],1e-12);
    }
  QuadraticPlanarPrecision::setPrecision(1e-14);
}

/*!
 * This test looks if intersectors are in coherency.
 */
void QuadraticPlanarInterpTest::checkEpsilonCoherency1()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-12);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-5);

  const double pol1[]={
    -2.1083388455000001, 1.2172499999999999,
    -1.7320508075000001, 1,
    -1.9201948265, 1.108625
  };

  const double pol2[]={
    -2.2379999998, 0,
    -1.9381648534, 1.1189999998,
    -2.1617419990000002, 0.57923702298000002,
    -1.9381648534, 1.1189999998,
    -1.9909924031999999, 1.1494999999,
    -1.9645786283, 1.1342499998
  };
  //
  Node *n1=new Node(pol1[0],pol1[1]);
  Node *n2=new Node(pol1[2],pol1[3]);
  Node *n3;
  //
  Edge *e1=new EdgeLin(n1,n2); n1->decrRef(); n2->decrRef();
  n1=new Node(pol2[0],pol2[1]);
  n2=new Node(pol2[4],pol2[5]);
  n3=new Node(pol2[2],pol2[3]);
  Edge *e2=new EdgeArcCircle(n1,n2,n3); n1->decrRef(); n2->decrRef(); n3->decrRef();
  e2->decrRef();
  e1->decrRef();
}

/*!
 * Tests to avoid regressions : Basic one.
 */
void QuadraticPlanarInterpTest::checkNonRegression1()
{
  const double coords1[]=
    {
      16.1732057215, -25.110999999800001,
      16.02555485246479, -25.340997988918762
    };
  Node *nS1=new Node(coords1);
  Node *nE1=new Node(coords1+2);
  const double radius1=2.902;
  const double angleS1=-0.49999999950907054; const double angleL1=-0.0942156629996692;
  const double center1[2]={13.66, -23.66};
  EdgeArcCircle *e1=new EdgeArcCircle(nS1,nE1,center1,radius1,angleS1,angleL1);
  //
  const double coords2[]=
    {
      16.041579804000001, -25.350249998999999,
      16.367740958999999, -24.132999999999999
    };
  Node *nS2=new Node(coords2);
  Node *nE2=new Node(coords2+2);
  const double radius2=2.4345;
  const double angleS2=-0.523598776190207; const double angleL2=0.5235987755846041;
  const double center2[]={ 13.933240960547204, -24.132999998525658 };
  EdgeArcCircle *e2=new EdgeArcCircle(nS2,nE2,center2,radius2,angleS2,angleL2);
  MergePoints merge;
  QuadraticPolygon c1,c2;
  e1->intersectWith(e2,merge,c1,c2);
  CPPUNIT_ASSERT_EQUAL(2,c1.size()); CPPUNIT_ASSERT_EQUAL(2,c2.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(e1->getCurveLength(),c1.getPerimeter(),1e-5);
  //clean-up
  nS1->decrRef(); nE1->decrRef(); nS2->decrRef(); nE2->decrRef(); e1->decrRef(); e2->decrRef();
}

void QuadraticPlanarInterpTest::checkNonRegression2()
{
  QuadraticPlanarPrecision::setPrecision(1e-12);
  QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-9);
  double coords1[]=
    {
      15.141499999899999, -26.226033271399999,
      16.226033271199999, -25.141499999800001,
      16.1732057215, -25.110999999800001,
      15.110999999899999, -26.1732057217,
      15.755157392699999, -25.755157392499999,
      16.199619496299999, -25.126249999799999,
      15.7120238788, -25.712023879099998,
      15.126249999899999, -26.199619496499999
    };
  double coords2[]=
    {
      15.933240959000001, -24.132999999999999,
      15.665291765999999, -25.132999998999999,
      16.041579804000001, -25.350249998999999,
      16.367740958999999, -24.132999999999999,
      15.865092611, -24.650638091000001,
      15.853435785, -25.241624998999999,
      16.284787383000001, -24.763094964,
      16.150490958999999, -24.132999999999999
    };
  std::vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::BuildArcCirclePolygon(nodes1);
  std::vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
  std::vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00173945,v[0]->getArea(),1e-7);
  delete v[0];
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00173945,pol1->intersectWith(*pol2),1e-7);
  delete pol1;
  delete pol2;
}

/*!
 * Tests to avoid regressions : Basic one.
 */
void QuadraticPlanarInterpTest::checkNonRegression3()
{
  const double coords1[]=
    {
      10.962340811000001, -22.417749999000002,
      12.217990959, -21.162099852000001
    };
  Node *nS1=new Node(coords1);
  Node *nE1=new Node(coords1+2);
  const double radius1=3.4304999897666599;
  const double angleS1=2.6179938783536514; const double angleL1=-0.52359877711901204;
  const double center1[2]={13.933240950441375, -24.132999992807399};
  EdgeArcCircle *e1=new EdgeArcCircle(nS1,nE1,center1,radius1,angleS1,angleL1);
  //
  const double coords2[]=
    {
      11.1467942784, -22.2090000002,
      11.0939667286, -22.178500000099998
    };
  Node *nS2=new Node(coords2);
  Node *nE2=new Node(coords2+2);
  EdgeLin *e2=new EdgeLin(nS2,nE2);
  MergePoints merge;
  QuadraticPolygon c1,c2;
  CPPUNIT_ASSERT(e1->intersectWith(e2,merge,c1,c2));
  CPPUNIT_ASSERT_EQUAL(2,c1.size());
  CPPUNIT_ASSERT_EQUAL(2,c2.size());
  ElementaryEdge *tmp1=dynamic_cast<ElementaryEdge *>(c1.front()); CPPUNIT_ASSERT(tmp1);
  EdgeArcCircle *tmp2=dynamic_cast<EdgeArcCircle *>(tmp1->getPtr()); CPPUNIT_ASSERT(tmp2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.6179938783536514,tmp2->getAngle0(),1e-14);
  //clean-up
  nS1->decrRef(); nE1->decrRef(); nS2->decrRef(); nE2->decrRef(); e1->decrRef(); e2->decrRef();
}

void QuadraticPlanarInterpTest::checkNonRegression4()
{
  QuadraticPlanarPrecision::setPrecision(1e-12);
  QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-9);
  double coords1[]=
    {
      10.962340811000001, -22.417749999000002,
      12.217990959, -21.162099852000001,
      12.051990958999999, -20.874579418,
      10.674820377, -22.251749999000001,
      11.507511146000001, -21.707270185999999,
      12.134990959, -21.018339635,
      11.272751694, -21.472510735,
      10.818580594, -22.334749999
    };

  double coords2[]=
    {
      10.758000000199999, -23.66,
      11.1467942784, -22.2090000002,
      11.0939667286, -22.178500000099998,
      10.696999999999999, -23.66,
      10.856883252299999, -22.908907131159999,
      11.1203805035, -22.1937500001,
      10.797961776699999, -22.893119169449999,
      10.727500000099999, -23.66
    };
  std::vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::BuildArcCirclePolygon(nodes1);
  std::vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
  std::vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00164773941455998,v[0]->getArea(),1e-7);
  delete v[0];
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00164773941455998,pol1->intersectWith(*pol2),1e-7);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression5()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-12);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-5);
  double coords1[]=
    {
      -1.7320508075000001, 1,
      -1, 1.7320508075000001 ,
      -1.2172499999999999, 2.1083388455000001,
      -2.1083388455000001, 1.2172499999999999,
      -1.4142135623000001, 1.4142135623000001,
      -1.108625, 1.9201948265,
      -1.7214514588000001, 1.7214514588000001,
      -1.9201948265, 1.108625};

  double coords2[]=
    {
      -2.2379999998, 0,
      -1.9381648534, 1.1189999998,
      -1.9909924031999999, 1.1494999999,
      -2.2989999998999999, 0,
      -2.1617419990000002, 0.57923702298000002,
      -1.9645786283, 1.1342499998,
      -2.2206634745999998, 0.59502498461999997,
      -2.2684999997999999, 0};
  //Edge1_of_pol2 inter Edge4_of_pol1 = {-1.9381648533711939, 1.1189999998498941}
  //Edge4_of_pol1 _angle = -0.523598775922546, _angle0 = -3.1415926535897931, _radius = 2.2379999983074721, _center = {-1.4925279436059493e-09, 1.3300635705141101e-10}}
  std::vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::BuildArcCirclePolygon(nodes1);
  std::vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
  std::vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(0,(int)v.size());
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00164773941455998,v[0]->getArea(),1e-7);
  //delete v[0];
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00164773941455998,pol1->intersectWith(*pol2),1e-7);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression6()
{
  QuadraticPlanarPrecision::setPrecision(1e-12);
  QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-5);
  double coords1[]=
    {
      10.962340811000001, -22.417749999000002,
      12.217990959, -21.162099852000001,
      12.051990958999999, -20.874579418,
      10.674820377, -22.251749999000001,
      11.507511146000001, -21.707270185999999,
      12.134990959, -21.018339635,
      11.272751694, -21.472510735,
      10.818580594, -22.334749999
    };
  double coords2[]=
    { 10.426, -23.66,
      10.859273844199999, -22.043000000100001,
      10.806446294799999, -22.012500000199999,
      10.3650000002, -23.66,
      10.536195877799999, -22.822979208099998,
      10.832860069499999, -22.027750000200001,
      10.477274402499999, -22.80719124657,
      10.3955000001, -23.66};
  std::vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::BuildArcCirclePolygon(nodes1);
  std::vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
  std::vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v[0]->getArea(),0.0150659,1e-7);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression7()
{
  QuadraticPlanarPrecision::setPrecision(1e-5);
  QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-5);
  double coords1[]=
    {
      -2., 0,
      -1.7320508075000001, 1,
      -2.1083388455000001, 1.2172499999999999,
      -2.4344999999999999, 0,
      -1.9318516525603098, 0.51763809027157182,
      -1.9201948265, 1.108625,
      -2.3515464241024469, 0.63009496529570408,
      -2.2172499999999999, 0
    };
  double coords2[]=
    { -2.3369999999000002, 0,
      -2.0239013684999998, 1.1684999999000001,
      -2.1927763221999998, 1.2659999998,
      -2.5319999998, 0,
      -2.2573686559260442, 0.60486010843437632,
      -2.1083388453499996, 1.2172499998499999,
      -2.445724191994314, 0.65532982205982326,
      -2.4344999998499999, 0 };
  std::vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::BuildArcCirclePolygon(nodes1);
  std::vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
  std::vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.121795,v[0]->getArea(),1.e-6);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression8()
{
  QuadraticPlanarPrecision::setPrecision(1e-3);
  QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-5);
  double coords1[]=
    {
      -13.933240959000001, -28.559499999,
      -16.146490959000001, -27.966461449000001,
      -16.383240958999998, -28.376524478,
      -13.933240959000001, -29.032999999000001,
      -15.078903461873765, -28.408670669106311,
      -16.264865958999998, -28.1714929635,
      -15.201454280317435, -28.866036547696734,
      -13.933240959000001, -28.796249999 };
  double coords2[]=
    { -16.382999999950002, -28.376524478457149,
      -13.933000000014729, -29.03299999982551,
      -13.93300000006697, -28.793999999915993,
      -16.263500000000001, -28.169544407039268,
      -15.201213320921273, -28.866036548734634,
      -13.933000000040851, -28.913499999870751,
      -15.139355569325469, -28.635180276305853,
      -16.323249999975001, -28.273034442748209 };
  std::vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::BuildArcCirclePolygon(nodes1);
  std::vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
  std::vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.598232,v[0]->getArea(),1.e-6);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression9()
{
  QuadraticPlanarPrecision::setPrecision(1e-7);
  QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-8);
  double coords1[]=
    {
      -0.04476229252902969, -0.085118027765365603,
      -0.046952683430894329, -0.085704941238358354,
      -0.046952683430894329, -0.088063823748058725,
      -0.043582851274179504, -0.087160879944491371,
      -0.045818853668170414, -0.085555669718918592,
      -0.046952683430894329, -0.086884382493208526,
      -0.045208329947517549, -0.087834175256748526,
      -0.044172571901604597, -0.086139453854928494 };

  double coords2[]=
    { -0.05065868681155701, -0.087744551996665671,
      -0.046951871439587615, -0.088737790182236015,
      -0.046951871439683469, -0.088063823751059062,
      -0.050321703596054014, -0.087160879946116557,
      -0.0488706602695924, -0.08848517684025306,
      -0.046951871439635542, -0.088400806966647538,
      -0.048696224921445964, -0.087834175258503858,
      -0.050490195203805516, -0.087452715971391121};

  std::vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::BuildArcCirclePolygon(nodes1);
  std::vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
  std::vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(0,(int)v.size());
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression10()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords1[]=
    { -0.002269581957210453, -0.09851030343724453,
      -0.004268022334182935, -0.1059685844580936,
      -0.002777851483521377, -0.1023709937816271};
  double coords2[]=
    { -0.004114727297178323, -0.1049870239624718,
      -0.003544545103522544, -0.1053162188055505};
  Node *n1_1=new Node(coords1);
  Node *n2_1=new Node(coords1+2);
  Node *n3_1=new Node(coords1+4);
  Node *n1_2=new Node(coords2);
  Node *n2_2=new Node(coords2+2);
  EdgeArcCircle *e1=new EdgeArcCircle(n1_1,n3_1,n2_1);
  EdgeLin *e2=new EdgeLin(n1_2,n2_2);
  MergePoints merge;
  ComposedEdge *c1=new ComposedEdge;
  ComposedEdge *c2=new ComposedEdge;
  CPPUNIT_ASSERT(e1->intersectWith(e2,merge,*c1,*c2));
  CPPUNIT_ASSERT_EQUAL(2,c1->size());
  CPPUNIT_ASSERT_EQUAL(2,c2->size());
  ComposedEdge::Delete(c1); ComposedEdge::Delete(c2);
  n1_1->decrRef(); n2_1->decrRef(); n3_1->decrRef();
  n1_2->decrRef(); n2_2->decrRef();
  e1->decrRef(); e2->decrRef();
}

void QuadraticPlanarInterpTest::checkNonRegression11()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords1[]=
    { -0.002269581957210453, -0.09851030343724453,
      -0.004268022334182935, -0.1059685844580936,
      -0.002886178753789801, -0.1067663922211958,
      -0.0006739664310059821, -0.09851030343724453,
      -0.002777851483521377, -0.1023709937816271,
      -0.003577100543986368, -0.1063674883396447,
      -0.001236605237717319, -0.1027839694676665,
      -0.001471774194108217, -0.09851030343724453};
  double coords2[]=
    { -0.003544545103522544, -0.1053162188055505,
      -0.001941023322604723, -0.09851030343724451,
      -0.002598140593501099, -0.09851030343724451,
      -0.004114727297178323, -0.1049870239624718,
      -0.002347317802266182, -0.1020064358043286,
      -0.002269581958052911, -0.09851030343724451,
      -0.002982346712452072, -0.1018362598405457,
      -0.003829636200350435, -0.1051516213840111};
  
  std::vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::BuildArcCirclePolygon(nodes1);
  std::vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
  std::vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.28973e-06,v[0]->getArea(),1.e-11);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression12()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-6);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords1[]=
    { -0.5032251558760915, -0.8716087994449138,
      -0.4695268343089433, -0.8806382374805872,
      -0.4695268343089433, -0.8570494123835835,
      -0.4914307433275896, -0.8511802776536561,
      -0.4869703691141082, -0.8783417525751493,
      -0.4695268343089433, -0.8688438249320853,
      -0.480865131947653, -0.8555566971861125,
      -0.4973279496018406, -0.8613945385492849};

  double coords2[]=
    { -0.5065868681155701, -0.8774455199666568,
      -0.4695187143958762, -0.8873779018223601,
      -0.4695187143968347, -0.8806382375105907,
      -0.5032170359605401, -0.8716087994611657,
      -0.488706602695924, -0.8848517684025307,
      -0.4695187143963554, -0.8840080696664754,
      -0.4869622492144596, -0.8783417525850385,
      -0.5049019520380551, -0.8745271597139112};

  std::vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::BuildArcCirclePolygon(nodes1);
  std::vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
  std::vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,v[0]->getArea(),1.e-6);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression13()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-6);

  double coords_1[194]={ 
    0, 0, 0.304375, -7.454791178893722e-17, 0.2152256265236553, -0.2152256265236555, -5.591093384170291e-17, -0.304375, 
    -0.2152256265236555, -0.2152256265236554, -0.304375, 3.727395589446861e-17, -0.2152256265236554, 0.2152256265236554, 1.86369779472343e-17, 0.304375, 
    0.2152256265236554, 0.2152256265236554, 0.60875, -1.490958235778744e-16, 0.5624116654162459, -0.2329585394522483, 0.4304512530473107, -0.4304512530473109, 
    0.2329585394522485, -0.5624116654162458, -1.118218676834058e-16, -0.60875, -0.2329585394522482, -0.5624116654162459, -0.4304512530473109, -0.4304512530473108, 
    -0.5624116654162459, -0.2329585394522483, -0.60875, 7.454791178893722e-17, -0.5624116654162458, 0.2329585394522485, -0.4304512530473108, 0.4304512530473109, 
    -0.2329585394522484, 0.5624116654162458, 3.727395589446861e-17, 0.60875, 0.2329585394522485, 0.5624116654162458, 0.4304512530473109, 0.4304512530473108, 
    0.5624116654162458, 0.2329585394522484, 0.913125, -2.236437353668116e-16, 0.645676879570966, -0.6456768795709663, -1.677328015251087e-16, -0.913125, 
    -0.6456768795709663, -0.6456768795709661, -0.913125, 1.118218676834058e-16, -0.6456768795709661, 0.6456768795709662, 5.591093384170291e-17, 0.913125, 
    0.6456768795709662, 0.6456768795709661, 1.2175, -2.981916471557489e-16, 1.124823330832492, -0.4659170789044966, 0.8609025060946214, -0.8609025060946218, 
    0.4659170789044971, -1.124823330832492, -2.236437353668116e-16, -1.2175, -0.4659170789044965, -1.124823330832492, -0.8609025060946218, -0.8609025060946216, 
    -1.124823330832492, -0.4659170789044967, -1.2175, 1.490958235778744e-16, -1.124823330832492, 0.465917078904497, -0.8609025060946216, 0.8609025060946217, 
    -0.4659170789044967, 1.124823330832492, 7.454791178893722e-17, 1.2175, 0.4659170789044969, 1.124823330832492, 0.8609025060946217, 0.8609025060946216, 
    1.124823330832492, 0.4659170789044968, 1.521875, -3.727395589446861e-16, 1.076128132618277, -1.076128132618277, -2.795546692085146e-16, -1.521875, 
    -1.076128132618277, -1.076128132618277, -1.521875, 1.86369779472343e-16, -1.076128132618277, 1.076128132618277, 9.318488973617152e-17, 1.521875, 
    1.076128132618277, 1.076128132618277, 1.82625, -4.472874707336233e-16, 1.687234996248738, -0.6988756183567448, 1.291353759141932, -1.291353759141933, 
    0.6988756183567456, -1.687234996248737, -3.354656030502175e-16, -1.82625, -0.6988756183567447, -1.687234996248738, -1.291353759141933, -1.291353759141932, 
    -1.687234996248738, -0.6988756183567449, -1.82625, 2.236437353668116e-16, -1.687234996248737, 0.6988756183567454, -1.291353759141932, 1.291353759141932, 
    -0.6988756183567451, 1.687234996248737, 1.118218676834058e-16, 1.82625, 0.6988756183567453, 1.687234996248737, 1.291353759141932, 1.291353759141932, 
    1.687234996248737, 0.6988756183567452, 2.130625, -5.218353825225606e-16, 1.506579385665588, -1.506579385665588, -3.913765368919204e-16, -2.130625, 
    -1.506579385665588, -1.506579385665588, -2.130625, 2.609176912612803e-16, -1.506579385665588, 1.506579385665588, 1.304588456306401e-16, 2.130625, 
    1.506579385665588, 1.506579385665588, 2.435, -5.963832943114977e-16, 2.249646661664984, -0.9318341578089931, 1.721805012189243, -1.721805012189244, 
    0.9318341578089941, -2.249646661664983, -4.472874707336233e-16, -2.435, -0.9318341578089929, -2.249646661664984, -1.721805012189244, -1.721805012189243, 
    -2.249646661664984, -0.9318341578089934, -2.435, 2.981916471557489e-16, -2.249646661664983, 0.9318341578089939, -1.721805012189243, 1.721805012189243, 
    -0.9318341578089935, 2.249646661664983, 1.490958235778744e-16, 2.435, 0.9318341578089938, 2.249646661664983, 1.721805012189243, 1.721805012189243, 
    2.249646661664983, 0.9318341578089936 };

  int tab6_1[48]={ 
    0, 9, 11, 1, 10, 2, 0, 11, 13, 2, 12, 3, 0, 13, 15, 3, 14, 4, 0, 15, 
    17, 4, 16, 5, 0, 17, 19, 5, 18, 6, 0, 19, 21, 6, 20, 7, 0, 21, 23, 7, 
    22, 8, 0, 23, 9, 8, 24, 1 };

  int tab8_1[192]={ 
    9, 33, 35, 11, 25, 34, 26, 10, 11, 35, 37, 13, 26, 36, 27, 12, 13, 37, 39, 15, 
    27, 38, 28, 14, 15, 39, 41, 17, 28, 40, 29, 16, 17, 41, 43, 19, 29, 42, 30, 18, 
    19, 43, 45, 21, 30, 44, 31, 20, 21, 45, 47, 23, 31, 46, 32, 22, 23, 47, 33, 9, 
    32, 48, 25, 24, 33, 57, 59, 35, 49, 58, 50, 34, 35, 59, 61, 37, 50, 60, 51, 36, 
    37, 61, 63, 39, 51, 62, 52, 38, 39, 63, 65, 41, 52, 64, 53, 40, 41, 65, 67, 43, 
    53, 66, 54, 42, 43, 67, 69, 45, 54, 68, 55, 44, 45, 69, 71, 47, 55, 70, 56, 46, 
    47, 71, 57, 33, 56, 72, 49, 48, 57, 81, 83, 59, 73, 82, 74, 58, 59, 83, 85, 61, 
    74, 84, 75, 60, 61, 85, 87, 63, 75, 86, 76, 62, 63, 87, 89, 65, 76, 88, 77, 64, 
    65, 89, 91, 67, 77, 90, 78, 66, 67, 91, 93, 69, 78, 92, 79, 68, 69, 93, 95, 71, 
    79, 94, 80, 70, 71, 95, 81, 57, 80, 96, 73, 72 };

  double coords_2[20]={ 
    0.5159941860137611, 0, 0, -0.5159941860137611, -0.5159941860137611, 0, 0, 0.5159941860137611, 
    0.6684941860137611, 0, 0, -0.6684941860137611, -0.6684941860137611, 0, 0, 0.6684941860137611, 
    0.5922441860137611, 0, -0.5922441860137611, 0 };
  
  int tab8_2[16]={ 
    0, 4, 6, 2, 8, 5, 9, 1, 2, 6, 4, 0, 9, 7, 8, 3 };
  
  double perimeterFromPol1,perimeterFromPol2,perimeterFromPol1AndPol2;

  const int *work1=tab6_1;
  for(int i=0;i<8;i++,work1+=6)
    {
      QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords_1,work1,6);
      const int *work2=tab8_2;
      for(int j=0;j<2;j++,work2+=8)
        {
          QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords_2,work2,8);
          //std::vector<int> tmp;
          //pol1->intersectForPoint(*pol2,tmp);
          pol1->intersectForPerimeter(*pol2,perimeterFromPol1,perimeterFromPol2,perimeterFromPol1AndPol2);
          //pol1->intersectMySelfWith(*pol2);
          delete pol2;
        }
      delete pol1;
    }
  work1=tab8_1;
  for(int i=0;i<24;i++,work1+=8)
    {
      QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords_1,work1,8);
      const int *work2=tab8_2;
      for(int j=0;j<2;j++,work2+=8)
        {
          
          QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords_2,work2,8);
          //std::vector<int> tmp;
          //pol1->intersectForPoint(*pol2,tmp);
          pol1->intersectForPerimeter(*pol2,perimeterFromPol1,perimeterFromPol2,perimeterFromPol1AndPol2);
          delete pol2;
        }
      delete pol1;
    }
}

/*!
  Some overlapping cases for intersectForPoint.
*/
void QuadraticPlanarInterpTest::checkNonRegression14()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-6);

  double coords[72]={
    1.,0.,1.3,0.,-1.3,0.,-1.,0.,1.15,0.,0.,1.3,-1.15,0.,0.,1.,
    -0.91923881554251186,-0.91923881554251186,-0.91923881554251186,0.91923881554251186,-1.0606601717798214,1.0606601717798214,-1.0606601717798214,-1.0606601717798214,-1.5,0.,
    -0.98994949366116658,-0.98994949366116658,-0.98994949366116658,0.98994949366116658,
    0.91923881554251186,0.91923881554251186,1.0606601717798214,1.0606601717798214,0.98994949366116658,0.98994949366116658, 0., 1.5,
    -0.83562389259250125,0.99585777605467141, -0.65, 1.1258330249197703, -1.2216004070216808, 0.44462618632336953, -1.1258330249197703, 0.65,
    -0.74564936725635955, 1.0648976575756897, -1.6770646146510724, 1.4072242996141826, -1.1782001231476449, 0.54940374026290939, -1.5873847317707279, 0.74020965686300877,
    -1.1782001231476449, 0.54940374026290939, -1.0648976575756894, 0.74564936725635977, -1.2950531075192693, -0.11330246557195534, -1.2950531075192693, 0.11330246557195565,
    -1.1258330249197703, 0.65, -2.1146554070041046, 0.56662020857685746, -1.6918048488667423, 0.45331774300490169,
    0.,-1.3,0.,-1.5
  };
  int tab[48]={
    0,1,2,3,4,5,6,7,
    8,9,10,11,2,14,12,13,
    9,15,16,10,5,17,18,14,
    9,15,16,10,34,17,35,14,
    19,20,21,22,23,24,25,26,
    27,28,29,30,31,32,2,33
  };
  QuadraticPolygon *pol1,*pol2;
  std::vector<int> goalOfTest;
  //
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab,8);
  // Level 1
  pol2=buildQuadraticPolygonCoarseInfo(coords,tab+8,8);
  pol1->intersectForPoint(*pol2,goalOfTest);
  const int res1[4]={0,1,0,0};
  CPPUNIT_ASSERT_EQUAL(4,(int)goalOfTest.size());
  CPPUNIT_ASSERT(equal(goalOfTest.begin(),goalOfTest.end(),res1));
  delete pol2;
  // Level 2
  pol2=buildQuadraticPolygonCoarseInfo(coords,tab+16,8);
  pol1->intersectForPoint(*pol2,goalOfTest);
  const int res2[4]={0,2,0,0};
  CPPUNIT_ASSERT_EQUAL(4,(int)goalOfTest.size());
  CPPUNIT_ASSERT(equal(goalOfTest.begin(),goalOfTest.end(),res2));
  delete pol2;
  //Level 2 bis
  pol2=buildQuadraticPolygonCoarseInfo(coords,tab+24,8);
  pol1->intersectForPoint(*pol2,goalOfTest);
  const int res2Bis[4]={0,2,0,0};
  CPPUNIT_ASSERT_EQUAL(4,(int)goalOfTest.size());
  CPPUNIT_ASSERT(equal(goalOfTest.begin(),goalOfTest.end(),res2Bis));
  delete pol2;
  // Level 3
  pol2=buildQuadraticPolygonCoarseInfo(coords,tab+40,8);
  pol1->intersectForPoint(*pol2,goalOfTest);
  const int res3[4]={0,3,0,0};
  CPPUNIT_ASSERT_EQUAL(4,(int)goalOfTest.size());
  CPPUNIT_ASSERT(equal(goalOfTest.begin(),goalOfTest.end(),res3));
  delete pol2;
  // Level 4
  pol2=buildQuadraticPolygonCoarseInfo(coords,tab+32,8);
  pol1->intersectForPoint(*pol2,goalOfTest);
  const int res4[4]={0,4,0,0};
  CPPUNIT_ASSERT_EQUAL(4,(int)goalOfTest.size());
  CPPUNIT_ASSERT(equal(goalOfTest.begin(),goalOfTest.end(),res4));
  delete pol2;
  //
  delete pol1;
}

/*!
 * This test is one of the most complicated intersection configuration.
 */
void QuadraticPlanarInterpTest::checkNonRegression15()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-6);

  double coords[72]={
    1.,0.,1.3,0.,-1.3,0.,-1.,0.,1.15,0.,0.,1.3,-1.15,0.,0.,1.,
    -0.91923881554251186,-0.91923881554251186,-0.91923881554251186,0.91923881554251186,-1.0606601717798214,1.0606601717798214,-1.0606601717798214,-1.0606601717798214,-1.5,0.,
    -0.98994949366116658,-0.98994949366116658,-0.98994949366116658,0.98994949366116658,
    0.91923881554251186,0.91923881554251186,1.0606601717798214,1.0606601717798214,0.98994949366116658,0.98994949366116658, 0., 1.5,
    -0.83562389259250125,0.99585777605467141, -0.65, 1.1258330249197703, -1.2216004070216808, 0.44462618632336953, -1.1258330249197703, 0.65,
    -0.74564936725635955, 1.0648976575756897, -1.6770646146510724, 1.4072242996141826, -1.1782001231476449, 0.54940374026290939, -1.5873847317707279, 0.74020965686300877,
    -1.1782001231476449, 0.54940374026290939, -1.0648976575756894, 0.74564936725635977, -1.2950531075192693, -0.11330246557195534, -1.2950531075192693, 0.11330246557195565,
    -1.1258330249197703, 0.65, -2.1146554070041046, 0.56662020857685746, -1.6918048488667423, 0.45331774300490169,
    0.,-1.3,0.,-1.5
  };

  int tab[24]={
    0,1,2,3,4,5,6,7,
    9,15,16,10,7,17,5,14,
    9,10,16,15,14,5,17,7
  };

  const double RefLgth=3.88995883524451;
  const double RefArea=0.383185168001075;
  //
  QuadraticPolygon *pol1,*pol2;
  //pol1 and pol2 in same orientation
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords,tab+8,8);
  std::vector<QuadraticPolygon *> res=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_EQUAL(4,res[0]->recursiveSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(RefLgth,res[0]->getPerimeter(),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(RefArea,res[0]->getArea(),1e-12);
  delete res[0];
  //pol1 and pol2 in same orientation but inversing intersection call pol1<->pol2
  res=pol2->intersectMySelfWith(*pol1);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_EQUAL(4,res[0]->recursiveSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(RefLgth,res[0]->getPerimeter(),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(RefArea,res[0]->getArea(),1e-12);
  delete res[0];
  delete pol2;
  //pol1 and pol2 in opposite orientation
  pol2=buildQuadraticPolygonCoarseInfo(coords,tab+16,8);
  res=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_EQUAL(4,res[0]->recursiveSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(RefLgth,res[0]->getPerimeter(),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-RefArea,res[0]->getArea(),1e-12);
  delete res[0];
  //pol1 and pol2 in opposite orientation but inversing intersection call pol1<->pol2
  res=pol2->intersectMySelfWith(*pol1);
  CPPUNIT_ASSERT_EQUAL(1,(int)res.size());
  CPPUNIT_ASSERT_EQUAL(4,res[0]->recursiveSize());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(RefLgth,res[0]->getPerimeter(),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(RefArea,res[0]->getArea(),1e-12);
  delete res[0];
  delete pol2;
  //
  delete pol1;
}

class DoubleEqual
{
public:
  DoubleEqual(double eps):_eps(eps) { }
  bool operator()(double x, double y) { return fabs(x-y)<_eps; }
private:
  double _eps;
};

/*!
 * This test is to see the reuse of a polygon in intersect* methods. initLocation needed ...
 */
void QuadraticPlanarInterpTest::checkNonRegression16()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords1[194]={ 
    0, 0, 0.304375, 0, 0.2152256265236554, 0.2152256265236554, 1.86369779472343e-17, 0.304375, 
    -0.2152256265236554, 0.2152256265236554, -0.304375, 3.727395589446861e-17, -0.2152256265236555, -0.2152256265236554, -5.591093384170291e-17, -0.304375, 
    0.2152256265236553, -0.2152256265236555, 0.60875, 0, 0.5624116654162458, 0.2329585394522484, 0.4304512530473109, 0.4304512530473108, 
    0.2329585394522485, 0.5624116654162458, 3.727395589446861e-17, 0.60875, -0.2329585394522484, 0.5624116654162458, -0.4304512530473108, 0.4304512530473109, 
    -0.5624116654162458, 0.2329585394522485, -0.60875, 7.454791178893722e-17, -0.5624116654162459, -0.2329585394522483, -0.4304512530473109, -0.4304512530473108, 
    -0.2329585394522482, -0.5624116654162459, -1.118218676834058e-16, -0.60875, 0.2329585394522485, -0.5624116654162458, 0.4304512530473107, -0.4304512530473109, 
    0.5624116654162459, -0.2329585394522483, 0.913125, 0, 0.6456768795709662, 0.6456768795709661, 5.591093384170291e-17, 0.913125, 
    -0.6456768795709661, 0.6456768795709662, -0.913125, 1.118218676834058e-16, -0.6456768795709663, -0.6456768795709661, -1.677328015251087e-16, -0.913125, 
    0.645676879570966, -0.6456768795709663, 1.2175, 0, 1.124823330832492, 0.4659170789044968, 0.8609025060946217, 0.8609025060946216, 
    0.4659170789044969, 1.124823330832492, 7.454791178893722e-17, 1.2175, -0.4659170789044967, 1.124823330832492, -0.8609025060946216, 0.8609025060946217, 
    -1.124823330832492, 0.465917078904497, -1.2175, 1.490958235778744e-16, -1.124823330832492, -0.4659170789044967, -0.8609025060946218, -0.8609025060946216, 
    -0.4659170789044965, -1.124823330832492, -2.236437353668116e-16, -1.2175, 0.4659170789044971, -1.124823330832492, 0.8609025060946214, -0.8609025060946218, 
    1.124823330832492, -0.4659170789044966, 1.521875, 0, 1.076128132618277, 1.076128132618277, 9.318488973617152e-17, 1.521875, 
    -1.076128132618277, 1.076128132618277, -1.521875, 1.86369779472343e-16, -1.076128132618277, -1.076128132618277, -2.795546692085146e-16, -1.521875, 
    1.076128132618277, -1.076128132618277, 1.82625, 0, 1.687234996248737, 0.6988756183567452, 1.291353759141932, 1.291353759141932, 
    0.6988756183567453, 1.687234996248737, 1.118218676834058e-16, 1.82625, -0.6988756183567451, 1.687234996248737, -1.291353759141932, 1.291353759141932, 
    -1.687234996248737, 0.6988756183567454, -1.82625, 2.236437353668116e-16, -1.687234996248738, -0.6988756183567449, -1.291353759141933, -1.291353759141932, 
    -0.6988756183567447, -1.687234996248738, -3.354656030502175e-16, -1.82625, 0.6988756183567456, -1.687234996248737, 1.291353759141932, -1.291353759141933, 
    1.687234996248738, -0.6988756183567448, 2.130625, 0, 1.506579385665588, 1.506579385665588, 1.304588456306401e-16, 2.130625, 
    -1.506579385665588, 1.506579385665588, -2.130625, 2.609176912612803e-16, -1.506579385665588, -1.506579385665588, -3.913765368919204e-16, -2.130625, 
    1.506579385665588, -1.506579385665588, 2.435, 0, 2.249646661664983, 0.9318341578089936, 1.721805012189243, 1.721805012189243, 
    0.9318341578089938, 2.249646661664983, 1.490958235778744e-16, 2.435, -0.9318341578089935, 2.249646661664983, -1.721805012189243, 1.721805012189243, 
    -2.249646661664983, 0.9318341578089939, -2.435, 2.981916471557489e-16, -2.249646661664984, -0.9318341578089934, -1.721805012189244, -1.721805012189243, 
    -0.9318341578089929, -2.249646661664984, -4.472874707336233e-16, -2.435, 0.9318341578089941, -2.249646661664983, 1.721805012189243, -1.721805012189244, 
    2.249646661664984, -0.9318341578089931, };

  int tab1_8[192]={ 
    11, 35, 33, 9, 26, 34, 25, 10, 13, 37, 35, 11, 27, 36, 26, 12, 15, 39, 37, 13, 
    28, 38, 27, 14, 17, 41, 39, 15, 29, 40, 28, 16, 19, 43, 41, 17, 30, 42, 29, 18, 
    21, 45, 43, 19, 31, 44, 30, 20, 23, 47, 45, 21, 32, 46, 31, 22, 9, 33, 47, 23, 
    25, 48, 32, 24, 35, 59, 57, 33, 50, 58, 49, 34, 37, 61, 59, 35, 51, 60, 50, 36, 
    39, 63, 61, 37, 52, 62, 51, 38, 41, 65, 63, 39, 53, 64, 52, 40, 43, 67, 65, 41, 
    54, 66, 53, 42, 45, 69, 67, 43, 55, 68, 54, 44, 47, 71, 69, 45, 56, 70, 55, 46, 
    33, 57, 71, 47, 49, 72, 56, 48, 59, 83, 81, 57, 74, 82, 73, 58, 61, 85, 83, 59, 
    75, 84, 74, 60, 63, 87, 85, 61, 76, 86, 75, 62, 65, 89, 87, 63, 77, 88, 76, 64, 
    67, 91, 89, 65, 78, 90, 77, 66, 69, 93, 91, 67, 79, 92, 78, 68, 71, 95, 93, 69, 
    80, 94, 79, 70, 57, 81, 95, 71, 73, 96, 80, 72, };

  double coords2[20]={ 
    2.435, 0, 0, -2.435, -2.435, 0, 0, 2.435, 
    2.6925, 0, 0, -2.6925, -2.6925, 0, 0, 2.6925, 
    2.56375, 0, -2.56375, 0, };

  int tab2_8[16]={ 0, 4, 6, 2, 8, 5, 9, 1, 2, 6, 4, 0, 9, 7, 8, 3 };

  QuadraticPolygon *pol1,*pol2;
  //pol1 and pol2 in same orientation
  std::vector<double> test1,test2;
  for(int ii=0;ii<24;ii++)
    {
      pol1=buildQuadraticPolygonCoarseInfo(coords1,tab1_8+8*ii,8);
      for(int jj=0;jj<2;jj++)
        {
          pol2=buildQuadraticPolygonCoarseInfo(coords2,tab2_8+jj*8,8);
          //
          std::vector<double> v1,v2;
          pol1->initLocations();
          pol1->intersectForPerimeterAdvanced(*pol2,v1,v2);
          if(ii==16 && jj==1)
            test1=v1;
          if(ii==20 && jj==1)
            test2=v1;
          delete pol2;
        }
      delete pol1;
    }
  const double test1_res[4]={0.,1.9124445278727873,0.,0.};
  CPPUNIT_ASSERT(std::equal(test1.begin(),test1.end(),test1_res,DoubleEqual(1e-10)));
  const double test2_res[4]={0.,0.,0.,0.};
  CPPUNIT_ASSERT(std::equal(test2.begin(),test2.end(),test2_res,DoubleEqual(1e-10)));
}

/*!
 * This test checks overlapped intersections END-INSIDE and INSIDE-START with same and opposite orientation.
 */
void QuadraticPlanarInterpTest::checkNonRegression17()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-7);
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(1e-7);
  double coords[16]={
    -1., 0., 1., 0. , 1.5, 0., -1.5, 0., 
    0. , 1., 1.25, 0., 0., 1.5, -1.25, 0.};
  
  double coords2[16]={
    0.70710678118654757, 0.70710678118654757, -1., 0., -1.25, 0.,  0.88388347648318444, 0.88388347648318444,
    0., -1., -1.125, 0., 0., -1.25, 0.79549512883486606, 0.79549512883486606 };

  double coords3[16]={
    0.70710678118654757, 0.70710678118654757, 0.88388347648318444, 0.88388347648318444, -1.25, 0., -1., 0.,
    0.79549512883486606, 0.79549512883486606, 0., -1.25, -1.125, 0., 0., -1. };

  int tab8[8]={
    0, 1, 2, 3, 4, 5, 6, 7 };
  QuadraticPolygon *pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  QuadraticPolygon *pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.22089323345553233,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords2,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.22089323345553233,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords3,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.22089323345553233,pol1->intersectWith(*pol2),1.e-13);
  delete pol1;
  delete pol2;
  pol1=buildQuadraticPolygonCoarseInfo(coords,tab8,8);
  pol2=buildQuadraticPolygonCoarseInfo(coords3,tab8,8);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.22089323345553233,pol2->intersectWith(*pol1),1.e-13);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNormalize()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-14);
  Node *n1=new Node(0.,0.);                Node *n4=new Node(0.,-3.);
  Node *n2=new Node(10.,0.);               Node *n5=new Node(10.,-3.);
  Node *n3=new Node(5.,10.);               Node *n6=new Node(5.,7.);
  EdgeLin *e1_2=new EdgeLin(n1,n2);        EdgeLin *e4_5=new EdgeLin(n4,n5);
  EdgeLin *e2_3=new EdgeLin(n2,n3);        EdgeLin *e5_6=new EdgeLin(n5,n6);
  EdgeLin *e3_1=new EdgeLin(n3,n1);        EdgeLin *e6_4=new EdgeLin(n6,n4);
  //
  QuadraticPolygon pol1; pol1.pushBack(e1_2); pol1.pushBack(e2_3); pol1.pushBack(e3_1);
  QuadraticPolygon pol2; pol2.pushBack(e4_5); pol2.pushBack(e5_6); pol2.pushBack(e6_4);
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();
  double area1Start=pol1.getArea();
  double xb,yb;
  double fact=pol1.normalize(&pol2,xb,yb);
  double area1End=pol1.getArea();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(area1Start,area1End*fact*fact,1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(13.,fact,1.e-14);
  double area=pol1.intersectWith(pol2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(24.5,area*fact*fact,1e-14);
  //
  n1=new Node(0.,0.);  n4=new Node(0.,-3.);
  n2=new Node(10.,0.); n5=new Node(10.,-3.);
  n3=new Node(5.,10.); n6=new Node(5.,7.);
  e1_2=new EdgeLin(n1,n2);        e4_5=new EdgeLin(n4,n5);
  e2_3=new EdgeLin(n2,n3);        e5_6=new EdgeLin(n5,n6);
  e3_1=new EdgeLin(n3,n1);        e6_4=new EdgeLin(n6,n4);
  QuadraticPolygon pol3; pol3.pushBack(e1_2); pol3.pushBack(e2_3); pol3.pushBack(e3_1);
  QuadraticPolygon pol4; pol4.pushBack(e4_5); pol4.pushBack(e5_6); pol4.pushBack(e6_4);
  n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(24.5,pol3.intersectWithAbs(pol4),1.e-14);
  // Ok testing EdgeArcCircle update.
  double center[2]={5.,5.};
  double radius=300.;
  EdgeArcCircle *e1=buildArcOfCircle(center,radius,M_PI/4.,M_PI/3.);
  const Bounds& b=e1->getBounds();
  double x,y,fact2;
  fact2=b.getCaracteristicDim();
  b.getBarycenter(x,y);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(78.539816339744817,e1->getCurveLength(),1e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(15106.061037591669,e1->getAreaOfZone(),1e-10);
  e1->getStartNode()->applySimilarity(x,y,fact2);
  e1->getEndNode()->applySimilarity(x,y,fact2);
  e1->applySimilarity(x,y,fact2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(62.132034355964237,fact2,1e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.2640792652913602,e1->getCurveLength(),1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.034741420428165526,e1->getAreaOfZone(),1e-13);
  e1->decrRef();
}

void QuadraticPlanarInterpTest::checkMakePartitionAbs1()
{
  INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(1e-14);
  Node *n0=new Node(0.,0.);                Node *n4=new Node(0.5,0.25);
  Node *n1=new Node(0.,0.5);               Node *n5=new Node(0.3,1.2);
  Node *n2=new Node(1.,0.5);               Node *n6=new Node(1.1,1.3);
  Node *n3=new Node(1.,0.);                Node *n7=new Node(-0.1,0.9);
  EdgeLin *e0_1=new EdgeLin(n0,n1);
  EdgeLin *e1_2=new EdgeLin(n1,n2);        EdgeLin *e4_5=new EdgeLin(n4,n5);
  EdgeLin *e2_3=new EdgeLin(n2,n3);        EdgeLin *e5_6=new EdgeLin(n5,n6);
  EdgeLin *e3_0=new EdgeLin(n3,n0);        EdgeLin *e6_4=new EdgeLin(n6,n4);
  EdgeLin *e4_7=new EdgeLin(n4,n7);        EdgeLin *e7_5=new EdgeLin(n7,n5);
  QuadraticPolygon pol1; pol1.pushBack(e0_1); pol1.pushBack(e1_2); pol1.pushBack(e2_3); pol1.pushBack(e3_0);
  QuadraticPolygon pol2; pol2.pushBack(e4_5); pol2.pushBack(e5_6); pol2.pushBack(e6_4);
  pol2.pushBack(e7_5); e4_5->incrRef(); pol2.pushBack(new ElementaryEdge(e4_5,false)); pol2.pushBack(e4_7);
  n0->decrRef(); n1->decrRef(); n2->decrRef(); n3->decrRef(); n4->decrRef(); n5->decrRef(); n6->decrRef(); n7->decrRef();
  pol1.dumpInXfigFileWithOther(pol2,"tony.fig");
}

}
