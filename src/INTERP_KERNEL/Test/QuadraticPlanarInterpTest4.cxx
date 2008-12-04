//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#include "QuadraticPlanarInterpTest.hxx"
#include "QuadraticPolygon.hxx"
#include "ElementaryEdge.hxx"
#include "EdgeArcCircle.hxx"
#include "EdgeLin.hxx"

#include <cmath>
#include <sstream>
#include <iostream>

using namespace std;
using namespace INTERP_KERNEL;

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
  vector<QuadraticPolygon *> result;
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
            double tmp1=0.,tmp2=0.,tmp3=0.,tmp4=0.;
            pol1.intersectForPerimeter(pol2,tmp1,tmp2,tmp3,tmp4);
            vector<double> v1,v2;
            vector<int> v3;
            double area2;
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,pol1.intersectForPerimeterAdvanced(pol2,v1,v2,area2),1.e-14);//no common edge
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
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.245,area2,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.7,tmp1,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(1.5652475842498528,tmp2,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,tmp3,1.e-14);//no common edge
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.245,tmp4,1.e-14);
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
            double tmp1=0.,tmp2=0.,tmp3=0.,tmp4=0.;
            pol7.intersectForPerimeter(pol8,tmp1,tmp2,tmp3,tmp4);
            vector<double> v1,v2;
            double area2;
            CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2360679774997898,pol7.intersectForPerimeterAdvanced(pol8,v1,v2,area2),1.e-14);//only common edges.
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,v1[0]+v1[1]+v1[2],1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,v2[0]+v2[1]+v2[2],1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,area2,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,tmp1,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,tmp2,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2360679774997898,tmp3,1.e-14);
            CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,tmp4,1.e-14);
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

/*!
 * Testing user interface high level function.
 */
void QuadraticPlanarInterpTest::checkHighLevelFunctionTest1()
{
  QUADRATIC_PLANAR::setPrecision(1e-12);
  QUADRATIC_PLANAR::setArcDetectionPrecision(1e-9);
  double coords[]={
    8.8334591186000004, 5.0999999999999996,
    7.1014083111000001, 6.0999999999999996,
    7.8334591186000004, 6.8320508074999999,
    7.9674337149000003, 5.5999999999999996,
    7.4192455562999999, 6.5142135623000001,
    8.3334591186000004, 5.9660254036999998
  };
  vector<Node *> nodes;
  nodes.push_back(new Node(coords));
  nodes.push_back(new Node(coords+2));
  nodes.push_back(new Node(coords+4));
  nodes.push_back(new Node(coords+6));
  nodes.push_back(new Node(coords+8));
  nodes.push_back(new Node(coords+10));
  QuadraticPolygon *pol=QuadraticPolygon::buildArcCirclePolygon(nodes);
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
  pol=QuadraticPolygon::buildArcCirclePolygon(nodes);
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
  pol=QuadraticPolygon::buildLinearPolygon(nodes);
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
  pol=QuadraticPolygon::buildArcCirclePolygon(nodes);
  pol->getBarycenter(tmp,tmp2);
  delete pol;
  QUADRATIC_PLANAR::setPrecision(1e-14);
}

void QuadraticPlanarInterpTest::check1DInterpLin()
{
  QUADRATIC_PLANAR::setPrecision(1e-7);
  QUADRATIC_PLANAR::setArcDetectionPrecision(1e-9);
  const int NB_OF_CELL_AXIAL_1=30;
  static const double Z_VALS_1[NB_OF_CELL_AXIAL_1+1]=
    { -0.1550 , -0.1356, -0.1162, -0.0969, -0.0775 ,-0.0581, -0.0387, -0.0194,  0.0000 , 0.0500, 
      0.1000 , 0.1500 , 0.2000 , 0.2500,  0.3000,  0.3500,  0.4000,  0.4500,  0.5000,  0.5500, 
      0.6000,  0.6500,  0.7000,  0.7194,  0.7388,  0.7581,  0.7775,  0.7969,  0.8163,  0.8356, 
      0.8550};
  vector<double> zLev1(Z_VALS_1,Z_VALS_1+NB_OF_CELL_AXIAL_1+1);

  const int NB_OF_CELL_AXIAL_2=46;
  static const double Z_VALS_2[NB_OF_CELL_AXIAL_2+1]=
    { -0.3050 ,-0.2863,-0.2675,-0.2488,-0.2300,-0.2113,-0.1925,-0.1738,-0.1550,-0.1356  
      , -0.1162,-0.0969,-0.0775,-0.0581,-0.0387,-0.0194,0.0000, 0.0500, 0.1 ,0.15 
      ,  0.20,  0.25, 0.30, 0.350 ,0.40 ,0.450 ,0.500 , 0.550, 0.600 ,0.650 ,0.700
      , 0.7194 ,0.7388 ,0.7581 ,0.7775 ,0.7969 ,0.8163 ,0.8356, 0.8550
      ,  0.8738 ,0.8925 ,0.9113 ,0.9300 ,0.9488 ,0.9675 ,0.9863, 1.0050};
  vector<double> zLev2(Z_VALS_2,Z_VALS_2+NB_OF_CELL_AXIAL_2+1);
  map<int,map<int,double> > m;
  Edge::interpolate1DLin(zLev1,zLev2,m);
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
  vector<double> zLev3(Z_VALS_3,Z_VALS_3+NB_OF_CELL_AXIAL_3+1);
  Edge::interpolate1DLin(zLev3,zLev1,m);
  CPPUNIT_ASSERT_EQUAL(13,(int)m.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,m[0][8],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,m[1][8],1e-12);
  for(int i=0;i<11;i++)
    {
      CPPUNIT_ASSERT_EQUAL(1,(int)m[i+2].size());
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,m[i+2][i+9],1e-12);
    }
  QUADRATIC_PLANAR::setPrecision(1e-14);
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
  QUADRATIC_PLANAR::setPrecision(1e-12);
  QUADRATIC_PLANAR::setArcDetectionPrecision(1e-9);
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
  vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::buildArcCirclePolygon(nodes1);
  vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
  vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
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
  QUADRATIC_PLANAR::setPrecision(1e-12);
  QUADRATIC_PLANAR::setArcDetectionPrecision(1e-9);
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
  vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::buildArcCirclePolygon(nodes1);
  vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
  vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00164773941455998,v[0]->getArea(),1e-7);
  delete v[0];
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00164773941455998,pol1->intersectWith(*pol2),1e-7);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression5()
{
  INTERP_KERNEL::QUADRATIC_PLANAR::setPrecision(1e-12);
  INTERP_KERNEL::QUADRATIC_PLANAR::setArcDetectionPrecision(1e-5);
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
  vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::buildArcCirclePolygon(nodes1);
  vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
  pol1->dumpInXfigFileWithOther(*pol2,"this.fig");
  vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(0,(int)v.size());
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00164773941455998,v[0]->getArea(),1e-7);
  //delete v[0];
  //CPPUNIT_ASSERT_DOUBLES_EQUAL(0.00164773941455998,pol1->intersectWith(*pol2),1e-7);
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression6()
{
  QUADRATIC_PLANAR::setPrecision(1e-12);
  QUADRATIC_PLANAR::setArcDetectionPrecision(1e-5);
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
  vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::buildArcCirclePolygon(nodes1);
  vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
  vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(v[0]->getArea(),0.0150659,1e-7);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression7()
{
  QUADRATIC_PLANAR::setPrecision(1e-5);
  QUADRATIC_PLANAR::setArcDetectionPrecision(1e-5);
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
  vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::buildArcCirclePolygon(nodes1);
  vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
  vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.121795,v[0]->getArea(),1.e-6);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression8()
{
  QUADRATIC_PLANAR::setPrecision(1e-3);
  QUADRATIC_PLANAR::setArcDetectionPrecision(1e-5);
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
  vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::buildArcCirclePolygon(nodes1);
  vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
  vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.598232,v[0]->getArea(),1.e-6);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression9()
{
  QUADRATIC_PLANAR::setPrecision(1e-7);
  QUADRATIC_PLANAR::setArcDetectionPrecision(1e-8);
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

  vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::buildArcCirclePolygon(nodes1);
  vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
  vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(0,(int)v.size());
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression10()
{
  INTERP_KERNEL::QUADRATIC_PLANAR::setPrecision(1e-7);
  INTERP_KERNEL::QUADRATIC_PLANAR::setArcDetectionPrecision(1e-7);
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
  INTERP_KERNEL::QUADRATIC_PLANAR::setPrecision(1e-7);
  INTERP_KERNEL::QUADRATIC_PLANAR::setArcDetectionPrecision(1e-7);
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
  
  vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::buildArcCirclePolygon(nodes1);
  vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
  vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.28973e-06,v[0]->getArea(),1.e-11);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNonRegression12()
{
  INTERP_KERNEL::QUADRATIC_PLANAR::setPrecision(1e-6);
  INTERP_KERNEL::QUADRATIC_PLANAR::setArcDetectionPrecision(1e-7);
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

  vector<Node *> nodes1;
  nodes1.push_back(new Node(coords1));
  nodes1.push_back(new Node(coords1+2));
  nodes1.push_back(new Node(coords1+4));
  nodes1.push_back(new Node(coords1+6));
  nodes1.push_back(new Node(coords1+8));
  nodes1.push_back(new Node(coords1+10));
  nodes1.push_back(new Node(coords1+12));
  nodes1.push_back(new Node(coords1+14));
  QuadraticPolygon *pol1=QuadraticPolygon::buildArcCirclePolygon(nodes1);
  vector<Node *> nodes2;
  nodes2.push_back(new Node(coords2));
  nodes2.push_back(new Node(coords2+2));
  nodes2.push_back(new Node(coords2+4));
  nodes2.push_back(new Node(coords2+6));
  nodes2.push_back(new Node(coords2+8));
  nodes2.push_back(new Node(coords2+10));
  nodes2.push_back(new Node(coords2+12));
  nodes2.push_back(new Node(coords2+14));
  QuadraticPolygon *pol2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
  vector<QuadraticPolygon *> v=pol1->intersectMySelfWith(*pol2);
  CPPUNIT_ASSERT_EQUAL(1,(int)v.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.28973e-06,v[0]->getArea(),1.e-11);
  delete v[0];
  delete pol1;
  delete pol2;
}

void QuadraticPlanarInterpTest::checkNormalize()
{
  INTERP_KERNEL::QUADRATIC_PLANAR::setPrecision(1e-14);
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
  double fact=pol1.normalize(&pol2);
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
