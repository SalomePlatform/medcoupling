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
    for(int i=0;i<1;i++)
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
             CPPUNIT_ASSERT_DOUBLES_EQUAL(0.15,result[0]->getAreaFast()+result[1]->getAreaFast(),1e-10);
             CPPUNIT_ASSERT_DOUBLES_EQUAL(0.03,fabs(result[0]->getAreaFast()-result[1]->getAreaFast()),1e-10);
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
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,pol1.getAreaFast(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2360679774997898,pol1.getPerimeterFast(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.61803398874989479,pol1.getHydroulicDiameter(),1e-10);
  pol1.reverse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.5,pol1.getAreaFast(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.2360679774997898,pol1.getPerimeterFast(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.61803398874989479,pol1.getHydroulicDiameter(),1e-10);
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
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1.3926990816987241,pol2.getAreaFast(),1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(4.5707963267948966,pol2.getPerimeterFast(),1e-6);
      pol2.reverse();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(-1.3926990816987241,pol2.getAreaFast(),1e-6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(4.5707963267948966,pol2.getPerimeterFast(),1e-6);
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
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.83775804095727857,pol3.getAreaFast(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.7775804095727832,pol3.getPerimeterFast(),1e-10);
  pol3.reverse();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.83775804095727857,pol3.getAreaFast(),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.7775804095727832,pol3.getPerimeterFast(),1e-10);
  //clean-up case3
  e1_2->decrRef(); e2_3->decrRef(); e3_4->decrRef(); e4_1->decrRef(); 
  n1->decrRef(); n2->decrRef(); n2m->decrRef(); n3->decrRef(); n3m->decrRef(); n4->decrRef();
}
