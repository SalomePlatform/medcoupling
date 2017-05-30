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

#include "MEDCouplingBasicsTest3.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingGaussLocalization.hxx"

#include <cmath>
#include <functional>
#include <iterator>

using namespace MEDCoupling;

void MEDCouplingBasicsTest3::testGetMeasureFieldCMesh1()
{
  MEDCouplingCMesh *m=MEDCouplingCMesh::New();
  DataArrayDouble *da=DataArrayDouble::New();
  const double discX[4]={2.3,3.4,5.8,10.2};
  const double discY[3]={12.3,23.4,45.8};
  const double discZ[5]={-0.7,1.2,1.25,2.13,2.67};
  da->alloc(4,1);
  std::copy(discX,discX+4,da->getPointer());
  m->setCoordsAt(0,da);
  da->decrRef();
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(4,m->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(1,m->getSpaceDimension());
  MEDCouplingFieldDouble *f=m->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(3,(int)f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)f->getNumberOfComponents());
  const double expected1[3]={1.1,2.4,4.4};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f->getIJ(i,0),1e-12);
  f->decrRef();
  DataArrayDouble *coords=m->getCoordinatesAndOwner();
  CPPUNIT_ASSERT_EQUAL(4,(int)coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)coords->getNumberOfComponents());
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(discX[i],coords->getIJ(i,0),1e-12);
  coords->decrRef();
  coords=m->computeCellCenterOfMass();
  CPPUNIT_ASSERT_EQUAL(3,(int)coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)coords->getNumberOfComponents());
  const double expected1_3[3]={2.85,4.6,8.};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1_3[i],coords->getIJ(i,0),1e-12);
  coords->decrRef();
  //
  da=DataArrayDouble::New();
  da->alloc(3,1);
  std::copy(discY,discY+3,da->getPointer());
  m->setCoordsAt(1,da);
  da->decrRef();
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(12,m->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(6,(int)m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(2,m->getSpaceDimension());
  f=m->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(6,(int)f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)f->getNumberOfComponents());
  const double expected2[6]={12.21,26.64,48.84,24.64,53.76,98.56};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f->getIJ(i,0),1e-12);
  f->decrRef();
  coords=m->getCoordinatesAndOwner();
  CPPUNIT_ASSERT_EQUAL(12,(int)coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)coords->getNumberOfComponents());
  const double expected2_2[24]={2.3,12.3,3.4,12.3,5.8,12.3,10.2,12.3, 2.3,23.4,3.4,23.4,5.8,23.4,10.2,23.4, 2.3,45.8,3.4,45.8,5.8,45.8,10.2,45.8};
  for(int i=0;i<24;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2_2[i],coords->getIJ(0,i),1e-12);
  coords->decrRef();
  coords=m->computeCellCenterOfMass();
  CPPUNIT_ASSERT_EQUAL(6,(int)coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)coords->getNumberOfComponents());
  const double expected2_3[12]={2.85,17.85,4.6,17.85,8.,17.85, 2.85,34.6,4.6,34.6,8.,34.6};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2_3[i],coords->getIJ(0,i),1e-12);
  coords->decrRef();
  //
  da=DataArrayDouble::New();
  da->alloc(5,1);
  std::copy(discZ,discZ+5,da->getPointer());
  m->setCoordsAt(2,da);
  da->decrRef();
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(60,m->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(24,(int)m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(3,m->getSpaceDimension());
  f=m->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(24,(int)f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)f->getNumberOfComponents());
  const double expected3[24]={23.199, 50.616, 92.796, 46.816, 102.144, 187.264, 0.6105, 1.332, 2.442, 1.232, 2.688, 4.928, 10.7448, 23.4432, 42.9792, 21.6832, 47.3088, 86.7328, 6.5934, 14.3856, 26.3736, 13.3056, 29.0304, 53.2224};
  for(int i=0;i<24;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],f->getIJ(i,0),1e-12);
  f->decrRef();
  coords=m->getCoordinatesAndOwner();
  CPPUNIT_ASSERT_EQUAL(60,(int)coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)coords->getNumberOfComponents());
  const double expected3_2[180]={
    2.3,12.3,-0.7, 3.4,12.3,-0.7, 5.8,12.3,-0.7, 10.2,12.3,-0.7, 2.3,23.4,-0.7, 3.4,23.4,-0.7, 5.8,23.4,-0.7, 10.2,23.4,-0.7, 2.3,45.8,-0.7, 3.4,45.8,-0.7, 5.8,45.8,-0.7, 10.2,45.8,-0.7,
    2.3,12.3,1.2, 3.4,12.3,1.2, 5.8,12.3,1.2, 10.2,12.3,1.2, 2.3,23.4,1.2, 3.4,23.4,1.2, 5.8,23.4,1.2, 10.2,23.4,1.2, 2.3,45.8,1.2, 3.4,45.8,1.2, 5.8,45.8,1.2, 10.2,45.8,1.2,
    2.3,12.3,1.25, 3.4,12.3,1.25, 5.8,12.3,1.25, 10.2,12.3,1.25, 2.3,23.4,1.25, 3.4,23.4,1.25, 5.8,23.4,1.25, 10.2,23.4,1.25, 2.3,45.8,1.25, 3.4,45.8,1.25, 5.8,45.8,1.25, 10.2,45.8,1.25,
    2.3,12.3,2.13, 3.4,12.3,2.13, 5.8,12.3,2.13, 10.2,12.3,2.13, 2.3,23.4,2.13, 3.4,23.4,2.13, 5.8,23.4,2.13, 10.2,23.4,2.13, 2.3,45.8,2.13, 3.4,45.8,2.13, 5.8,45.8,2.13, 10.2,45.8,2.13,
    2.3,12.3,2.67, 3.4,12.3,2.67, 5.8,12.3,2.67, 10.2,12.3,2.67, 2.3,23.4,2.67, 3.4,23.4,2.67, 5.8,23.4,2.67, 10.2,23.4,2.67, 2.3,45.8,2.67, 3.4,45.8,2.67, 5.8,45.8,2.67, 10.2,45.8,2.67
  };
  for(int i=0;i<180;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3_2[i],coords->getIJ(0,i),1e-12);
  coords->decrRef();
  coords=m->computeCellCenterOfMass();
  CPPUNIT_ASSERT_EQUAL(24,(int)coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)coords->getNumberOfComponents());
  const double expected3_3[72]={
    2.85,17.85,0.25,4.6,17.85,0.25,8.,17.85,0.25, 2.85,34.6,0.25,4.6,34.6,0.25,8.,34.6,0.25,
    2.85,17.85,1.225,4.6,17.85,1.225,8.,17.85,1.225, 2.85,34.6,1.225,4.6,34.6,1.225,8.,34.6,1.225,
    2.85,17.85,1.69,4.6,17.85,1.69,8.,17.85,1.69, 2.85,34.6,1.69,4.6,34.6,1.69,8.,34.6,1.69,
    2.85,17.85,2.4,4.6,17.85,2.4,8.,17.85,2.4, 2.85,34.6,2.4,4.6,34.6,2.4,8.,34.6,2.4
  };
  for(int i=0;i<72;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3_3[i],coords->getIJ(0,i),1e-12);
  coords->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest3::testFieldDoubleZipCoords1()
{
  MEDCouplingUMesh *m=build2DTargetMeshMergeNode_1();
  MEDCouplingFieldDouble *f=m->fillFromAnalytic(ON_NODES,2,"x*2.");
  f->getArray()->setInfoOnComponent(0,"titi");
  f->getArray()->setInfoOnComponent(1,"tutu");
  f->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(18,(int)f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)f->getNumberOfComponents());
  const double expected1[36]={-0.6, -0.6, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 0.4, 0.4};
  for(int i=0;i<36;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(f->zipCoords());
  f->checkConsistencyLight();
  const double expected2[30]={-0.6, -0.6, 1.4, 1.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 0.4, 0.4};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(!f->zipCoords());
  f->checkConsistencyLight();
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(std::string(f->getArray()->getInfoOnComponent(0))=="titi");
  CPPUNIT_ASSERT(std::string(f->getArray()->getInfoOnComponent(1))=="tutu");
  f->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest3::testFieldDoubleZipConnectivity1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  const int cells1[3]={2,3,4};
  MEDCouplingPointSet *m3_1=m2->buildPartOfMySelf(cells1,cells1+3,true);
  MEDCouplingUMesh *m3=dynamic_cast<MEDCouplingUMesh *>(m3_1);
  CPPUNIT_ASSERT(m3);
  m2->decrRef();
  MEDCouplingUMesh *m4=build2DSourceMesh_1();
  MEDCouplingUMesh *m5=MEDCouplingUMesh::MergeUMeshes(m1,m3);
  m1->decrRef();
  m3->decrRef();
  MEDCouplingUMesh *m6=MEDCouplingUMesh::MergeUMeshes(m5,m4);
  m4->decrRef();
  m5->decrRef();
  //
  CPPUNIT_ASSERT_EQUAL(10,(int)m6->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(22,m6->getNumberOfNodes());
  bool areNodesMerged;
  int newNbOfNodes;
  DataArrayInt *arr=m6->mergeNodes(1e-13,areNodesMerged,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(9,m6->getNumberOfNodes());
  arr->decrRef();
  MEDCouplingFieldDouble *f=m6->fillFromAnalytic(ON_CELLS,2,"x");
  MEDCouplingFieldDouble *f2=m6->fillFromAnalytic(ON_NODES,2,"x");
  CPPUNIT_ASSERT_EQUAL(10,(int)f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)f->getNumberOfComponents());
  const double expected1[20]={-0.05, -0.05, 0.3666666666666667, 0.3666666666666667, 0.53333333333333321, 0.53333333333333321,
                              -0.05, -0.05, 0.45, 0.45, 0.53333333333333321, 0.53333333333333321, -0.05, -0.05, 0.45, 0.45,
                              0.36666666666666659, 0.36666666666666659, 0.033333333333333326, 0.033333333333333326};
  for(int i=0;i<20;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f->getIJ(0,i),1e-12);
  f->getArray()->setInfoOnComponent(0,"titi");
  f->getArray()->setInfoOnComponent(1,"tutu");
  f->checkConsistencyLight();
  CPPUNIT_ASSERT(f->zipConnectivity(0));
  const double expected2[14]={-0.05, -0.05, 0.3666666666666667, 0.3666666666666667, 0.53333333333333321, 0.53333333333333321,
                              -0.05, -0.05, 0.45, 0.45, 0.36666666666666659, 0.36666666666666659, 0.033333333333333326, 0.033333333333333326};
  CPPUNIT_ASSERT_EQUAL(7,(int)f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)f->getNumberOfComponents());
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(std::string(f->getArray()->getInfoOnComponent(0))=="titi");
  CPPUNIT_ASSERT(std::string(f->getArray()->getInfoOnComponent(1))=="tutu");
  CPPUNIT_ASSERT(!f->zipConnectivity(0));
  f->decrRef();
  //
  const double expected3[18]={-0.3, -0.3, 0.2, 0.2, 0.7, 0.7, -0.3, -0.3, 0.2, 0.2, 0.7, 0.7, 
                              -0.3, -0.3, 0.2, 0.2, 0.7, 0.7};
  CPPUNIT_ASSERT_EQUAL(9,(int)f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)f2->getNumberOfComponents());
  for(int i=0;i<18;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],f2->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(f2->zipConnectivity(0));
  CPPUNIT_ASSERT_EQUAL(9,(int)f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)f2->getNumberOfComponents());
  for(int i=0;i<18;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],f2->getIJ(0,i),1e-12);
  f2->decrRef();
  //
  m6->decrRef();
}

void MEDCouplingBasicsTest3::testDaDoubleRenumber1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(7,2);
  a->setInfoOnComponent(0,"toto");
  a->setInfoOnComponent(1,"tata");
  const double arr1[14]={1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1};
  std::copy(arr1,arr1+14,a->getPointer());
  //
  const int arr2[7]={3,1,0,6,5,4,2};
  DataArrayDouble *b=a->renumber(arr2);
  CPPUNIT_ASSERT_EQUAL(7,(int)b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)b->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(1))=="tata");
  const double expected1[14]={3.1, 13.1, 2.1, 12.1, 7.1, 17.1, 1.1, 11.1, 6.1, 16.1, 5.1, 15.1, 4.1, 14.1};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],b->getIJ(0,i),1e-14);
  b->decrRef();
  a->decrRef();
  //
  DataArrayInt *c=DataArrayInt::New();
  c->alloc(7,2);
  c->setInfoOnComponent(0,"toto");
  c->setInfoOnComponent(1,"tata");
  const int arr3[14]={1,11,2,12,3,13,4,14,5,15,6,16,7,17};
  std::copy(arr3,arr3+14,c->getPointer());
  DataArrayInt *d=c->renumber(arr2);
  CPPUNIT_ASSERT_EQUAL(7,(int)d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(1))=="tata");
  const int expected2[14]={3, 13, 2, 12, 7, 17, 1, 11, 6, 16, 5, 15, 4, 14};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d->getIJ(0,i));
  c->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest3::testDaDoubleRenumberAndReduce1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(7,2);
  a->setInfoOnComponent(0,"toto");
  a->setInfoOnComponent(1,"tata");
  const double arr1[14]={1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1};
  std::copy(arr1,arr1+14,a->getPointer());
  //
  const int arr2[7]={2,-1,1,-1,0,4,3};
  DataArrayDouble *b=a->renumberAndReduce(arr2,5);
  CPPUNIT_ASSERT_EQUAL(5,(int)b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)b->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(1))=="tata");
  const double expected1[10]={5.1,15.1,3.1,13.1,1.1,11.1,7.1,17.1,6.1,16.1};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],b->getIJ(0,i),1e-14);
  b->decrRef();
  a->decrRef();
  //
  DataArrayInt *c=DataArrayInt::New();
  c->alloc(7,2);
  c->setInfoOnComponent(0,"toto");
  c->setInfoOnComponent(1,"tata");
  const int arr3[14]={1,11,2,12,3,13,4,14,5,15,6,16,7,17};
  std::copy(arr3,arr3+14,c->getPointer());
  DataArrayInt *d=c->renumberAndReduce(arr2,5);
  CPPUNIT_ASSERT_EQUAL(5,(int)d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(1))=="tata");
  const int expected2[10]={5,15,3,13,1,11,7,17,6,16};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d->getIJ(0,i));
  c->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest3::testDaDoubleRenumberInPlace1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(7,2);
  const double arr1[14]={1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1};
  std::copy(arr1,arr1+14,a->getPointer());
  //
  const int arr2[7]={3,1,0,6,5,4,2};
  a->renumberInPlace(arr2);
  CPPUNIT_ASSERT_EQUAL(7,(int)a->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)a->getNumberOfComponents());
  const double expected1[14]={3.1, 13.1, 2.1, 12.1, 7.1, 17.1, 1.1, 11.1, 6.1, 16.1, 5.1, 15.1, 4.1, 14.1};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],a->getIJ(0,i),1e-14);
  a->decrRef();
  //
  DataArrayInt *c=DataArrayInt::New();
  c->alloc(7,2);
  const int arr3[14]={1,11,2,12,3,13,4,14,5,15,6,16,7,17};
  std::copy(arr3,arr3+14,c->getPointer());
  c->renumberInPlace(arr2);
  CPPUNIT_ASSERT_EQUAL(7,(int)c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)c->getNumberOfComponents());
  const int expected2[14]={3, 13, 2, 12, 7, 17, 1, 11, 6, 16, 5, 15, 4, 14};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],c->getIJ(0,i));
  c->decrRef();
}

void MEDCouplingBasicsTest3::testDaDoubleRenumberR1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(7,2);
  a->setInfoOnComponent(0,"toto");
  a->setInfoOnComponent(1,"tata");
  const double arr1[14]={1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1};
  std::copy(arr1,arr1+14,a->getPointer());
  //
  const int arr2[7]={3,1,0,6,5,4,2};
  DataArrayDouble *b=a->renumberR(arr2);
  CPPUNIT_ASSERT_EQUAL(7,(int)b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)b->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(1))=="tata");
  const double expected1[14]={4.1, 14.1, 2.1, 12.1, 1.1, 11.1, 7.1, 17.1, 6.1, 16.1, 5.1, 15.1, 3.1, 13.1};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],b->getIJ(0,i),1e-14);
  b->decrRef();
  a->decrRef();
  //
  DataArrayInt *c=DataArrayInt::New();
  c->alloc(7,2);
  c->setInfoOnComponent(0,"toto");
  c->setInfoOnComponent(1,"tata");
  const int arr3[14]={1,11,2,12,3,13,4,14,5,15,6,16,7,17};
  std::copy(arr3,arr3+14,c->getPointer());
  DataArrayInt *d=c->renumberR(arr2);
  CPPUNIT_ASSERT_EQUAL(7,(int)d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(1))=="tata");
  const int expected2[14]={4, 14, 2, 12, 1, 11, 7, 17, 6, 16, 5, 15, 3, 13};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d->getIJ(0,i));
  c->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest3::testDaDoubleRenumberInPlaceR1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(7,2);
  const double arr1[14]={1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1};
  std::copy(arr1,arr1+14,a->getPointer());
  //
  const int arr2[7]={3,1,0,6,5,4,2};
  a->renumberInPlaceR(arr2);
  CPPUNIT_ASSERT_EQUAL(7,(int)a->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)a->getNumberOfComponents());
  const double expected1[14]={4.1, 14.1, 2.1, 12.1, 1.1, 11.1, 7.1, 17.1, 6.1, 16.1, 5.1, 15.1, 3.1, 13.1};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],a->getIJ(0,i),1e-14);
  a->decrRef();
  //
  DataArrayInt *c=DataArrayInt::New();
  c->alloc(7,2);
  const int arr3[14]={1,11,2,12,3,13,4,14,5,15,6,16,7,17};
  std::copy(arr3,arr3+14,c->getPointer());
  c->renumberInPlaceR(arr2);
  CPPUNIT_ASSERT_EQUAL(7,(int)c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)c->getNumberOfComponents());
  const int expected2[14]={4, 14, 2, 12, 1, 11, 7, 17, 6, 16, 5, 15, 3, 13};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],c->getIJ(0,i));
  c->decrRef();
}

void MEDCouplingBasicsTest3::testDaDoubleSelectByTupleId1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(7,2);
  a->setInfoOnComponent(0,"toto");
  a->setInfoOnComponent(1,"tata");
  const double arr1[14]={1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1};
  std::copy(arr1,arr1+14,a->getPointer());
  //
  const int arr2[7]={4,2,0,6,5};
  DataArrayDouble *b=a->selectByTupleId(arr2,arr2+5);
  CPPUNIT_ASSERT_EQUAL(5,(int)b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)b->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(1))=="tata");
  const double expected1[10]={5.1,15.1,3.1,13.1,1.1,11.1,7.1,17.1,6.1,16.1};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],b->getIJ(0,i),1e-14);
  b->decrRef();
  a->decrRef();
  //
  DataArrayInt *c=DataArrayInt::New();
  c->alloc(7,2);
  c->setInfoOnComponent(0,"toto");
  c->setInfoOnComponent(1,"tata");
  const int arr3[14]={1,11,2,12,3,13,4,14,5,15,6,16,7,17};
  std::copy(arr3,arr3+14,c->getPointer());
  DataArrayInt *d=c->selectByTupleId(arr2,arr2+5);
  CPPUNIT_ASSERT_EQUAL(5,(int)d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(1))=="tata");
  const int expected2[10]={5,15,3,13,1,11,7,17,6,16};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d->getIJ(0,i));
  c->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest3::testDaDoubleGetMinMaxValues1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(9,1);
  const double arr1[9]={2.34,4.56,-6.77,4.55,4.56,2.24,2.34,1.02,4.56};
  std::copy(arr1,arr1+9,a->getPointer());
  int where;
  double m=a->getMaxValue(where);
  CPPUNIT_ASSERT_EQUAL(1,where);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.56,m,1e-12);
  DataArrayInt *ws;
  m=a->getMaxValue2(ws);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.56,m,1e-12);
  CPPUNIT_ASSERT_EQUAL(3,(int)ws->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)ws->getNumberOfComponents());
  const int expected1[3]={1,4,8};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],ws->getIJ(i,0));
  ws->decrRef();
  a->decrRef();
  a=DataArrayDouble::New();
  const double arr2[9]={-2.34,-4.56,6.77,-4.55,-4.56,-2.24,-2.34,-1.02,-4.56};
  a->alloc(9,1);
  std::copy(arr2,arr2+9,a->getPointer());
  where=-2;
  m=a->getMinValue(where);
  CPPUNIT_ASSERT_EQUAL(1,where);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-4.56,m,1e-12);
  m=a->getMinValue2(ws);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-4.56,m,1e-12);
  CPPUNIT_ASSERT_EQUAL(3,(int)ws->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)ws->getNumberOfComponents());
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],ws->getIJ(i,0));
  ws->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest3::testFieldDoubleGetMinMaxValues2()
{
  MEDCouplingUMesh *m1=0;
  MEDCouplingUMesh *m2=build3DExtrudedUMesh_1(m1);
  m1->decrRef();
  CPPUNIT_ASSERT_EQUAL(18,(int)m2->getNumberOfCells());
  const double arr1[18]={8.71,4.53,-12.41,8.71,-8.71,8.7099,4.55,8.71,5.55,6.77,-1e-200,4.55,8.7099,0.,1.23,0.,2.22,8.71};
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(18,1);
  std::copy(arr1,arr1+18,a->getPointer());
  f->setArray(a);
  a->decrRef();
  f->setMesh(m2);
  //
  f->checkConsistencyLight();
  double m=f->getMaxValue();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.71,m,1e-12);
  DataArrayInt *ws;
  m=f->getMaxValue2(ws);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.71,m,1e-12);
  CPPUNIT_ASSERT_EQUAL(4,(int)ws->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)ws->getNumberOfComponents());
  const int expected1[4]={0,3,7,17};
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],ws->getIJ(i,0));
  ws->decrRef();
  //
  const double arr2[18]={-8.71,-4.53,12.41,-8.71,8.71,-8.7099,-4.55,-8.71,-5.55,-6.77,1e-200,-4.55,-8.7099,0.,-1.23,0.,-2.22,-8.71};
  std::copy(arr2,arr2+18,a->getPointer());
  f->checkConsistencyLight();
  m=f->getMinValue();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-8.71,m,1e-12);
  m=f->getMinValue2(ws);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-8.71,m,1e-12);
  CPPUNIT_ASSERT_EQUAL(4,(int)ws->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)ws->getNumberOfComponents());
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],ws->getIJ(i,0));
  ws->decrRef();
  //
  f->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest3::testBuildUnstructuredCMesh1()
{
  MEDCouplingCMesh *m=MEDCouplingCMesh::New();
  DataArrayDouble *da=DataArrayDouble::New();
  const double discX[4]={2.3,3.4,5.8,10.2};
  const double discY[3]={12.3,23.4,45.8};
  const double discZ[5]={-0.7,1.2,1.25,2.13,2.67};
  da->alloc(4,1);
  std::copy(discX,discX+4,da->getPointer());
  m->setCoordsAt(0,da);
  da->decrRef();
  m->checkConsistencyLight();
  double pos=2.4;
  CPPUNIT_ASSERT_EQUAL(0,m->getCellContainingPoint(&pos,1e-12));
  pos=3.7;
  CPPUNIT_ASSERT_EQUAL(1,m->getCellContainingPoint(&pos,1e-12));
  pos=5.9;
  CPPUNIT_ASSERT_EQUAL(2,m->getCellContainingPoint(&pos,1e-12));
  pos=10.3;
  CPPUNIT_ASSERT_EQUAL(-1,m->getCellContainingPoint(&pos,1e-12));
  pos=1.3;
  CPPUNIT_ASSERT_EQUAL(-1,m->getCellContainingPoint(&pos,1e-12));
  //
  MEDCouplingUMesh *m2=m->buildUnstructured();
  m2->checkConsistencyLight();
  MEDCouplingFieldDouble *f1=m->getMeasureField(false);
  MEDCouplingFieldDouble *f2=m2->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL((int)f1->getNumberOfTuples(),3);
  CPPUNIT_ASSERT_EQUAL((int)f2->getNumberOfTuples(),3);
  CPPUNIT_ASSERT_EQUAL(1,m2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(1,m2->getSpaceDimension());
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f1->getIJ(i,0),f2->getIJ(i,0),1e-10);
  da=DataArrayDouble::New();
  da->alloc(3,1);
  std::copy(discY,discY+3,da->getPointer());
  m->setCoordsAt(1,da);
  da->decrRef();
  m2->decrRef();
  f1->decrRef();
  f2->decrRef();
  //
  m2=m->buildUnstructured();
  m2->checkConsistencyLight();
  f1=m->getMeasureField(false);
  f2=m2->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL((int)f1->getNumberOfTuples(),6);
  CPPUNIT_ASSERT_EQUAL((int)f2->getNumberOfTuples(),6);
  CPPUNIT_ASSERT_EQUAL(2,m2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(2,m2->getSpaceDimension());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f1->getIJ(i,0),f2->getIJ(i,0),1e-10);
  f1->decrRef();
  f2->decrRef();
  m2->decrRef();
  //
  da=DataArrayDouble::New();
  da->alloc(5,1);
  std::copy(discZ,discZ+5,da->getPointer());
  m->setCoordsAt(2,da);
  da->decrRef();
  m2=m->buildUnstructured();
  m2->checkConsistencyLight();
  f1=m->getMeasureField(false);
  f2=m2->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL((int)f1->getNumberOfTuples(),24);
  CPPUNIT_ASSERT_EQUAL((int)f2->getNumberOfTuples(),24);
  CPPUNIT_ASSERT_EQUAL(3,m2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,m2->getSpaceDimension());
  for(int i=0;i<24;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f1->getIJ(i,0),f2->getIJ(i,0),1e-10);
  f1->decrRef();
  f2->decrRef();
  //
  double pos1[3]={5.,30.,2.};
  CPPUNIT_ASSERT_EQUAL(16,m->getCellContainingPoint(pos1,1e-12));
  //
  const double pt[3]={2.4,12.7,-3.4};
  m->scale(pt,3.7);
  MEDCouplingUMesh *m3=m->buildUnstructured();
  m2->scale(pt,3.7);
  CPPUNIT_ASSERT(m3->isEqual(m2,1e-12));
  m2->decrRef();
  m3->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest3::testDataArrayIntInvertO2NNO21()
{
  const int arr1[6]={2,0,4,1,5,3};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(6,1);
  std::copy(arr1,arr1+6,da->getPointer());
  DataArrayInt *da2=da->invertArrayO2N2N2O(6);
  CPPUNIT_ASSERT_EQUAL(6,(int)da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  const int expected1[6]={1,3,0,5,2,4};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(i,0));
  DataArrayInt *da3=da2->invertArrayN2O2O2N(6);
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(arr1[i],da3->getIJ(i,0));
  da3->decrRef();
  da2->decrRef();
  da->decrRef();
  //
  const int arr2[10]={3,-1,5,4,-1,0,-1,1,2,-1};
  da=DataArrayInt::New();
  da->alloc(10,1);
  std::copy(arr2,arr2+10,da->getPointer());
  da2=da->invertArrayO2N2N2O(6);
  CPPUNIT_ASSERT_EQUAL(6,(int)da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  const int expected2[10]={5,7,8,0,3,2};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],da2->getIJ(i,0));
  da3=da2->invertArrayN2O2O2N(10);
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(arr2[i],da3->getIJ(i,0));
  da3->decrRef();
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest3::testKeepSetSelectedComponent1()
{
  const double arr1[20]={1.,2.,3.,4., 11.,12.,13.,14., 21.,22.,23.,24., 31.,32.,33.,34., 41.,42.,43.,44.};
  DataArrayDouble *a1=DataArrayDouble::New();
  a1->alloc(5,4);
  std::copy(arr1,arr1+20,a1->getPointer());
  a1->setInfoOnComponent(0,"aaaa");
  a1->setInfoOnComponent(1,"bbbb");
  a1->setInfoOnComponent(2,"cccc");
  a1->setInfoOnComponent(3,"dddd");
  const int arr2[6]={1,2,1,2,0,0};
  std::vector<int> arr2V(arr2,arr2+6);
  DataArrayDouble *a2=static_cast<DataArrayDouble *>(a1->keepSelectedComponents(arr2V));
  CPPUNIT_ASSERT_EQUAL(6,(int)a2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,(int)a2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(0))=="bbbb");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(1))=="cccc");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(2))=="bbbb");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(3))=="cccc");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(4))=="aaaa");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(5))=="aaaa");
  const double expected1[30]={2.,3.,2.,3.,1.,1., 12.,13.,12.,13.,11.,11., 22.,23.,22.,23.,21.,21., 32.,33.,32.,33.,31.,31., 42.,43.,42.,43.,41.,41.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],a2->getIJ(0,i),1e-14);
  MCAuto<DataArrayInt> a3(a1->convertToIntArr());
  DataArrayInt *a4=static_cast<DataArrayInt *>(a3->keepSelectedComponents(arr2V));
  CPPUNIT_ASSERT_EQUAL(6,(int)a4->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,(int)a4->getNumberOfTuples());
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(0))=="bbbb");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(1))=="cccc");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(2))=="bbbb");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(3))=="cccc");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(4))=="aaaa");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(5))=="aaaa");
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_EQUAL(int(expected1[i]),a4->getIJ(0,i));
  // setSelectedComponents
  const int arr3[2]={3,2};
  std::vector<int> arr3V(arr3,arr3+2);
  DataArrayDouble *a5=static_cast<DataArrayDouble *>(a1->keepSelectedComponents(arr3V));
  a5->setInfoOnComponent(0,"eeee");
  a5->setInfoOnComponent(1,"ffff");
  const int arr4[2]={1,2};
  std::vector<int> arr4V(arr4,arr4+2);
  a2->setSelectedComponents(a5,arr4V);
  CPPUNIT_ASSERT_EQUAL(6,(int)a2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,(int)a2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(0))=="bbbb");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(1))=="eeee");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(2))=="ffff");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(3))=="cccc");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(4))=="aaaa");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(5))=="aaaa");
  const double expected2[30]={2.,4.,3.,3.,1.,1., 12.,14.,13.,13.,11.,11., 22.,24.,23.,23.,21.,21., 32.,34.,33.,33.,31.,31., 42.,44.,43.,43.,41.,41.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],a2->getIJ(0,i),1e-14);
  MCAuto<DataArrayInt> a6=a5->convertToIntArr();
  a6->setInfoOnComponent(0,"eeee");
  a6->setInfoOnComponent(1,"ffff");
  a4->setSelectedComponents(a6,arr4V);
  CPPUNIT_ASSERT_EQUAL(6,(int)a4->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,(int)a4->getNumberOfTuples());
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(0))=="bbbb");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(1))=="eeee");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(2))=="ffff");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(3))=="cccc");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(4))=="aaaa");
  CPPUNIT_ASSERT(std::string(a4->getInfoOnComponent(5))=="aaaa");
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_EQUAL(int(expected2[i]),a4->getIJ(0,i));
  // test of throw
  const int arr5[3]={2,3,6};
  const int arr6[3]={2,7,5};
  const int arr7[4]={2,1,4,6};
  std::vector<int> arr5V(arr5,arr5+3);
  std::vector<int> arr6V(arr6,arr6+3);
  std::vector<int> arr7V(arr7,arr7+4);
  CPPUNIT_ASSERT_THROW(a2->keepSelectedComponents(arr5V),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(a2->keepSelectedComponents(arr6V),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(a2->setSelectedComponents(a1,arr7V),INTERP_KERNEL::Exception);
  arr7V.resize(3);
  CPPUNIT_ASSERT_THROW(a2->setSelectedComponents(a1,arr7V),INTERP_KERNEL::Exception);
  //
  a5->decrRef();
  a4->decrRef();
  a2->decrRef();
  a1->decrRef();
}

void MEDCouplingBasicsTest3::testKeepSetSelectedComponent2()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  const double arr1[20]={1.,2.,3.,4., 11.,12.,13.,14., 21.,22.,23.,24., 31.,32.,33.,34., 41.,42.,43.,44.};
  DataArrayDouble *a1=DataArrayDouble::New();
  a1->alloc(5,4);
  std::copy(arr1,arr1+20,a1->getPointer());
  a1->setInfoOnComponent(0,"aaaa");
  a1->setInfoOnComponent(1,"bbbb");
  a1->setInfoOnComponent(2,"cccc");
  a1->setInfoOnComponent(3,"dddd");
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setTime(2.3,4,5);
  f1->setMesh(m1);
  f1->setName("f1");
  f1->setArray(a1);
  f1->checkConsistencyLight();
  //
  const int arr2[6]={1,2,1,2,0,0};
  std::vector<int> arr2V(arr2,arr2+6);
  MEDCouplingFieldDouble *f2=f1->keepSelectedComponents(arr2V);
  CPPUNIT_ASSERT(f2->getMesh()==f1->getMesh());
  CPPUNIT_ASSERT(f2->getTimeDiscretization()==ONE_TIME);
  int dt,it;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.3,f2->getTime(dt,it),1e-13);
  CPPUNIT_ASSERT_EQUAL(4,dt);
  CPPUNIT_ASSERT_EQUAL(5,it);
  f2->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(6,(int)f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,(int)f2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(0))=="bbbb");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(1))=="cccc");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(2))=="bbbb");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(3))=="cccc");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(4))=="aaaa");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(5))=="aaaa");
  const double expected1[30]={2.,3.,2.,3.,1.,1., 12.,13.,12.,13.,11.,11., 22.,23.,22.,23.,21.,21., 32.,33.,32.,33.,31.,31., 42.,43.,42.,43.,41.,41.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f2->getIJ(0,i),1e-14);
  //setSelectedComponents
  const int arr3[2]={3,2};
  std::vector<int> arr3V(arr3,arr3+2);
  MEDCouplingFieldDouble *f5=f1->keepSelectedComponents(arr3V);
  f5->setTime(6.7,8,9);
  f5->getArray()->setInfoOnComponent(0,"eeee");
  f5->getArray()->setInfoOnComponent(1,"ffff");
  f5->checkConsistencyLight();
  const int arr4[2]={1,2};
  std::vector<int> arr4V(arr4,arr4+2);
  f2->setSelectedComponents(f5,arr4V);
  CPPUNIT_ASSERT_EQUAL(6,(int)f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,(int)f2->getNumberOfTuples());
  f2->checkConsistencyLight();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.3,f2->getTime(dt,it),1e-13);
  CPPUNIT_ASSERT_EQUAL(4,dt);
  CPPUNIT_ASSERT_EQUAL(5,it);
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(0))=="bbbb");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(1))=="eeee");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(2))=="ffff");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(3))=="cccc");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(4))=="aaaa");
  CPPUNIT_ASSERT(std::string(f2->getArray()->getInfoOnComponent(5))=="aaaa");
  const double expected2[30]={2.,4.,3.,3.,1.,1., 12.,14.,13.,13.,11.,11., 22.,24.,23.,23.,21.,21., 32.,34.,33.,33.,31.,31., 42.,44.,43.,43.,41.,41.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f2->getIJ(0,i),1e-14);
  f5->decrRef();
  f1->decrRef();
  f2->decrRef();
  a1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest3::testElementaryDAThrowAndSpecialCases()
{
  DataArrayInt *da=DataArrayInt::New();
  CPPUNIT_ASSERT_THROW(da->checkAllocated(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->fillWithValue(1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->iota(1),INTERP_KERNEL::Exception);
  da->alloc(7,1);
  da->fillWithValue(11); //11,11,11,11...
  da->iota(10); //10,11,12,13...
  
  DataArrayInt *db=DataArrayInt::New();
  db->alloc(7,2);
  
  DataArrayDouble *dbl2=DataArrayDouble::New();
  dbl2->alloc(7,2);
  CPPUNIT_ASSERT_THROW(dbl2->isUniform(10.,1e-15),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl2->sort(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl2->iota(10.),INTERP_KERNEL::Exception);
  
  DataArrayDouble *dbl=DataArrayDouble::New();
  //DataArrayDouble not allocated yet
  CPPUNIT_ASSERT_THROW(dbl->iota(10.),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->isUniform(10.,1e-15),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->sort(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->reverse(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->fromNoInterlace(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->toNoInterlace(),INTERP_KERNEL::Exception);
  
  dbl->alloc(7,1);
  dbl->iota(10.);
  CPPUNIT_ASSERT(!dbl->isUniform(10.,1e-15));
  dbl->sort();
  CPPUNIT_ASSERT(dbl->isMonotonic(true, .99));
  CPPUNIT_ASSERT(dbl->isMonotonic(true, -.99));
  CPPUNIT_ASSERT(!dbl->isMonotonic(true, 1.1));
  CPPUNIT_ASSERT(!dbl->isMonotonic(true, -1.1));
  dbl->reverse();
  CPPUNIT_ASSERT(dbl->isMonotonic(false, .99));
  CPPUNIT_ASSERT(!dbl->isMonotonic(false, 1.1));
  CPPUNIT_ASSERT(!dbl->isMonotonic(false, -1.1));
  
  DataArrayInt *dc=DataArrayInt::New();
  dc->alloc(14,1);
  
  DataArrayDouble *dd=DataArrayDouble::New();
  CPPUNIT_ASSERT_THROW(dd->checkAllocated(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dd->fillWithValue(1.),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dd->iota(1.),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(!((dd->repr().find("No data"))==std::string::npos));
  
  dd->alloc(0,1); //Allocated but nbOfElements==0!
  CPPUNIT_ASSERT(!((dd->repr().find("Number of tuples : 0"))==std::string::npos));
  CPPUNIT_ASSERT(!((dd->repr().find("Empty Data"))==std::string::npos));
  dd->fillWithValue(11); //?!...
  dd->iota(10); //?!...
  CPPUNIT_ASSERT(dd->isMonotonic(true, 1.));
  CPPUNIT_ASSERT(dd->isMonotonic(false, 1.));
 
  CPPUNIT_ASSERT_THROW(db->copyStringInfoFrom(*da),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(db->copyStringInfoFrom(*da),INTERP_KERNEL::Exception);
  std::vector<int> cIds(2,2);
  CPPUNIT_ASSERT_THROW(da->copyPartOfStringInfoFrom(*db,cIds),INTERP_KERNEL::Exception);
  cIds[0]=1;
  cIds[0]=-1;
  CPPUNIT_ASSERT_THROW(da->copyPartOfStringInfoFrom(*db,cIds),INTERP_KERNEL::Exception);
  
  std::vector<std::string> info(2,"infoOfOneComponent");
  CPPUNIT_ASSERT_THROW(da->setInfoOnComponents(info),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->setInfoOnComponent(1,info[0].c_str()),INTERP_KERNEL::Exception);
  db->setInfoOnComponents(info);
  
  CPPUNIT_ASSERT_THROW(da->getInfoOnComponent(-1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->getInfoOnComponent(2),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT(db->getInfoOnComponent(1)==db->getInfoOnComponent(0));
  CPPUNIT_ASSERT_THROW(db->getVarOnComponent(-1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(db->getVarOnComponent(2),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(db->getUnitOnComponent(-1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(db->getUnitOnComponent(2),INTERP_KERNEL::Exception);
  
  CPPUNIT_ASSERT(da->GetVarNameFromInfo(std::string("varname unit "))==std::string("varname unit "));
  CPPUNIT_ASSERT(da->GetVarNameFromInfo(std::string("varname]unit["))==std::string("varname]unit["));
  CPPUNIT_ASSERT(da->GetVarNameFromInfo(std::string("[unit]"))==std::string());
  CPPUNIT_ASSERT(da->GetVarNameFromInfo(std::string("varname [unit]"))==std::string("varname"));
  
  CPPUNIT_ASSERT(da->GetUnitFromInfo(std::string("varname unit "))==std::string());
  CPPUNIT_ASSERT(da->GetUnitFromInfo(std::string("varname]unit["))==std::string());
  CPPUNIT_ASSERT(da->GetUnitFromInfo(std::string("[unit]"))==std::string("unit"));
  CPPUNIT_ASSERT(da->GetUnitFromInfo(std::string("varname [unit]"))==std::string("unit"));
  
  CPPUNIT_ASSERT_THROW(da->checkNbOfTuplesAndComp(*db,"theMessageInThrow"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->checkNbOfTuplesAndComp(*dc,"theMessageInThrow"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(db->checkNbOfTuplesAndComp(*dc,"theMessageInThrow"),INTERP_KERNEL::Exception);

  CPPUNIT_ASSERT_THROW(da->checkNbOfTuplesAndComp(7,2,"theMessageInThrow"),INTERP_KERNEL::Exception);
  da->checkNbOfTuplesAndComp(7,1,"theMessageInThrow");
  
  CPPUNIT_ASSERT_THROW(db->checkNbOfElems(7*2+1,"theMessageInThrow"),INTERP_KERNEL::Exception);
  db->checkNbOfElems(7*2,"theMessageInThrow");
  
  CPPUNIT_ASSERT_THROW(db->GetNumberOfItemGivenBES(10,9,1,"theMessageInThrow"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(db->GetNumberOfItemGivenBES(0,1,-1,"theMessageInThrow"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_EQUAL(10,db->GetNumberOfItemGivenBES(0,10,1,"theMessageInThrow"));
  CPPUNIT_ASSERT_EQUAL(5,db->GetNumberOfItemGivenBES(0,10,2,"theMessageInThrow"));
  CPPUNIT_ASSERT_EQUAL(6,db->GetNumberOfItemGivenBES(0,11,2,"theMessageInThrow"));
  
  //std::cout<<"\n!!!!!!!!!\n"<<dd->repr()<<"\n!!!!!!!!!\n";
  CPPUNIT_ASSERT(!((da->repr().find("Number of components : 1"))==std::string::npos));
  CPPUNIT_ASSERT(!((dd->repr().find("Number of components : 1"))==std::string::npos));
  CPPUNIT_ASSERT(!((dbl->repr().find("Number of components : 1"))==std::string::npos));
  
  CPPUNIT_ASSERT(!((da->reprZip().find("Number of components : 1"))==std::string::npos));
  CPPUNIT_ASSERT(!((dd->reprZip().find("Number of components : 1"))==std::string::npos));
  CPPUNIT_ASSERT(!((dbl->reprZip().find("Number of components : 1"))==std::string::npos));
  
  std::ostringstream ret;
  dbl->writeVTK(ret,2,"file.tmp",0);
  CPPUNIT_ASSERT(!((ret.str().find("<DataArray"))==std::string::npos));
  CPPUNIT_ASSERT(!((ret.str().find("Float32"))==std::string::npos));
  CPPUNIT_ASSERT(!((ret.str().find("16 15 14 13 12 11 10"))==std::string::npos));
  
  CPPUNIT_ASSERT_THROW(dbl->selectByTupleIdSafeSlice(0,1,-1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->subArray(-1,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->subArray(8,1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->subArray(0,8),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->meldWith(dd),INTERP_KERNEL::Exception);
  
  CPPUNIT_ASSERT_THROW(dbl->setPartOfValuesAdv(dbl2,da),INTERP_KERNEL::Exception); //dbl dbl2 not have the same number of components
  CPPUNIT_ASSERT_THROW(dbl->setPartOfValuesAdv(dd,da),INTERP_KERNEL::Exception);  //da tuple selector DataArrayInt instance not have exactly 2 components
  
  DataArrayDouble *dbl3=DataArrayDouble::New();
  dbl3->alloc(6,2);
  dbl3->fillWithValue(11.);
  int tupleId;
  //bad number of components
  CPPUNIT_ASSERT_THROW(dbl3->getMaxValue(tupleId),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dd->getMaxValue(tupleId),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->getMinValue(tupleId),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dd->getMinValue(tupleId),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->getAverageValue(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dd->getAverageValue(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dd->accumulate(100),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl->fromPolarToCart(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->fromCylToCart(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->fromSpherToCart(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->doublyContractedProduct(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->determinant(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->eigenValues(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->eigenVectors(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->inverse(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->trace(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->deviator(),INTERP_KERNEL::Exception);
 
  dbl3->setIJ(5,1,12.);
  CPPUNIT_ASSERT(dbl3->getMaxValueInArray()==12.);
  CPPUNIT_ASSERT(dbl3->getMinValueInArray()==11.);
 
  db->fillWithValue(100); //bad Ids
  CPPUNIT_ASSERT_THROW(dbl3->setPartOfValuesAdv(dbl2,db),INTERP_KERNEL::Exception);
  db->fillWithValue(-1); //bad Ids
  CPPUNIT_ASSERT_THROW(dbl3->setPartOfValuesAdv(dbl2,db),INTERP_KERNEL::Exception);
  db->fillWithValue(6); //bad Ids for dbl3
  CPPUNIT_ASSERT_THROW(dbl3->setPartOfValuesAdv(dbl2,db),INTERP_KERNEL::Exception);
  
  DataArrayDouble::SetArrayIn(dbl,dbl3); //dbl->dbl3 memLeaks?
  dbl3->checkNoNullValues();
  dbl3->setIJ(6,0,0.);
  CPPUNIT_ASSERT_THROW(dbl3->checkNoNullValues(),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(dbl3->applyInv(1.),INTERP_KERNEL::Exception);  //div by zero
  CPPUNIT_ASSERT_THROW(dbl2->findIdsInRange(1.,2.),INTERP_KERNEL::Exception);
  std::vector<const DataArrayDouble *> a(0); //input list must be NON EMPTY
  CPPUNIT_ASSERT_THROW(DataArrayDouble::Aggregate(a),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(DataArrayDouble::Meld(a),INTERP_KERNEL::Exception);
  
  a.push_back(dbl2);
  a.push_back(dbl); //Nb of components mismatch
  CPPUNIT_ASSERT_THROW(DataArrayDouble::Aggregate(a),INTERP_KERNEL::Exception);
  
  CPPUNIT_ASSERT_THROW(DataArrayDouble::Dot(dbl2,dbl),INTERP_KERNEL::Exception);
  
  CPPUNIT_ASSERT_THROW(DataArrayDouble::CrossProduct(dbl2,dbl),INTERP_KERNEL::Exception); //Nb of components mismatch
  CPPUNIT_ASSERT_THROW(DataArrayDouble::CrossProduct(dbl2,dbl2),INTERP_KERNEL::Exception); //Nb of components must be equal to 3 
  DataArrayDouble *dbl4=DataArrayDouble::New();
  dbl4->alloc(6,3);
  DataArrayDouble *dbl5=DataArrayDouble::New();
  dbl5->alloc(7,3);
  CPPUNIT_ASSERT_THROW(DataArrayDouble::CrossProduct(dbl4,dbl5),INTERP_KERNEL::Exception); //Nb of tuples mismatch
  
  a[0]=dbl4; //Nb of tuple mismatch
  a[1]=dbl5; //Nb of tuple mismatch
  CPPUNIT_ASSERT_THROW(DataArrayDouble::Meld(a),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(DataArrayDouble::Dot(dbl4,dbl5),INTERP_KERNEL::Exception);
  
  da->decrRef();
  db->decrRef();
  dbl->decrRef();
  dbl2->decrRef();
  dbl3->decrRef();
  dbl4->decrRef();
  dbl5->decrRef();
  dc->decrRef();
  dd->decrRef();
}

void MEDCouplingBasicsTest3::testDAIGetIdsEqual1()
{
  const int tab1[7]={5,-2,-4,-2,3,2,-2};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(7,1);
  std::copy(tab1,tab1+7,da->getPointer());
  DataArrayInt *da2=da->findIdsEqual(-2);
  CPPUNIT_ASSERT_EQUAL(3,(int)da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  const int expected1[3]={1,3,6};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,da2->getConstPointer()));
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest3::testDAIGetIdsEqualList1()
{
  const int tab1[7]={5,-2,-4,-2,3,2,-2};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(7,1);
  std::copy(tab1,tab1+7,da->getPointer());
  const int tab2[3]={3,-2,0};
  std::vector<int> tab2V(tab2,tab2+3);
  DataArrayInt *da2=da->findIdsEqualList(&tab2V[0],&tab2V[0]+tab2V.size());
  CPPUNIT_ASSERT_EQUAL(4,(int)da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  const int expected1[4]={1,3,4,6};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+4,da2->getConstPointer()));
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest3::testDAFromNoInterlace1()
{
  const int tab1[15]={1,11,21,31,41,2,12,22,32,42,3,13,23,33,43};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(5,3);
  std::copy(tab1,tab1+15,da->getPointer());
  DataArrayInt *da2=da->fromNoInterlace();
  const int expected1[15]={1,2,3,11,12,13,21,22,23,31,32,33,41,42,43};
  CPPUNIT_ASSERT_EQUAL(5,(int)da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)da2->getNumberOfComponents());// it's not a bug. Avoid to have 1 million components !
  CPPUNIT_ASSERT(std::equal(expected1,expected1+15,da2->getConstPointer()));
  MCAuto<DataArrayDouble> da3=da->convertToDblArr();
  DataArrayDouble *da4=da3->fromNoInterlace();
  CPPUNIT_ASSERT_EQUAL(5,(int)da4->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)da4->getNumberOfComponents());// it's not a bug. Avoid to have 1 million components !
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)expected1[i],da4->getIJ(0,i),1e-14);
  da4->decrRef();
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest3::testDAToNoInterlace1()
{
  const int tab1[15]={1,2,3,11,12,13,21,22,23,31,32,33,41,42,43};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(5,3);
  std::copy(tab1,tab1+15,da->getPointer());
  DataArrayInt *da2=da->toNoInterlace();
  const int expected1[15]={1,11,21,31,41,2,12,22,32,42,3,13,23,33,43};
  CPPUNIT_ASSERT_EQUAL(5,(int)da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)da2->getNumberOfComponents());// it's not a bug. Avoid to have 1 million components !
  CPPUNIT_ASSERT(std::equal(expected1,expected1+15,da2->getConstPointer()));
  MCAuto<DataArrayDouble> da3=da->convertToDblArr();
  DataArrayDouble *da4=da3->toNoInterlace();
  CPPUNIT_ASSERT_EQUAL(5,(int)da4->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)da4->getNumberOfComponents());// it's not a bug. Avoid to have 1 million components !
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)expected1[i],da4->getIJ(0,i),1e-14);
  da4->decrRef();
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest3::testDAIsUniform1()
{
  const int tab1[5]={1,1,1,1,1};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(5,1);
  std::copy(tab1,tab1+5,da->getPointer());
  CPPUNIT_ASSERT(da->isUniform(1));
  da->setIJ(2,0,2);
  CPPUNIT_ASSERT(!da->isUniform(1));
  da->setIJ(2,0,1);
  CPPUNIT_ASSERT(da->isUniform(1));
  MCAuto<DataArrayDouble> da2=da->convertToDblArr();
  CPPUNIT_ASSERT(da2->isUniform(1.,1e-12));
  da2->setIJ(1,0,1.+1.e-13);
  CPPUNIT_ASSERT(da2->isUniform(1.,1e-12));
  da2->setIJ(1,0,1.+1.e-11);
  CPPUNIT_ASSERT(!da2->isUniform(1.,1e-12));
  da->decrRef();
}

void MEDCouplingBasicsTest3::testDADFromPolarToCart1()
{
  const double tab1[4]={2.,0.2,2.5,0.7};
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(2,2);
  std::copy(tab1,tab1+4,da->getPointer());
  DataArrayDouble *da2=da->fromPolarToCart();
  const double expected1[4]={1.9601331556824833,0.39733866159012243, 1.9121054682112213,1.6105442180942275};
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da2->getIJ(0,i),1e-13);
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest3::testDADFromCylToCart1()
{
  const double tab1[6]={2.,0.2,4.,2.5,0.7,9.};
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(2,3);
  std::copy(tab1,tab1+6,da->getPointer());
  DataArrayDouble *da2=da->fromCylToCart();
  const double expected1[6]={1.9601331556824833,0.39733866159012243,4., 1.9121054682112213,1.6105442180942275,9.};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da2->getIJ(0,i),1e-13);
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest3::testDADFromSpherToCart1()
{
  const double tab1[6]={2.,0.2,0.3,2.5,0.7,0.8};
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(2,3);
  std::copy(tab1,tab1+6,da->getPointer());
  DataArrayDouble *da2=da->fromSpherToCart();
  const double expected1[6]={0.37959212195737485,0.11742160338765303,1.9601331556824833, 1.1220769624465328,1.1553337045129035,1.9121054682112213};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da2->getIJ(0,i),1e-13);
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest3::testUnPolyze1()
{
  const int elts[8]={0,1,2,3,4,5,6,7};
  std::vector<int> eltsV(elts,elts+8);
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  mesh->unPolyze();
  MEDCouplingUMesh *mesh2=build3DTargetMesh_1();
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT(mesh->isEqual(mesh2,1e-12));
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  CPPUNIT_ASSERT(!mesh->isEqual(mesh2,1e-12));
  mesh->getNodalConnectivity()->setIJ(0,6,10);
  mesh->getNodalConnectivity()->setIJ(0,7,9);
  mesh->getNodalConnectivity()->setIJ(0,8,12);
  mesh->getNodalConnectivity()->setIJ(0,9,13);
  mesh->unPolyze();
  CPPUNIT_ASSERT(mesh->isEqual(mesh2,1e-12));
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  mesh->getNodalConnectivity()->setIJ(0,6,12);
  mesh->getNodalConnectivity()->setIJ(0,7,13);
  mesh->getNodalConnectivity()->setIJ(0,8,10);
  mesh->getNodalConnectivity()->setIJ(0,9,9);
  mesh->unPolyze();
  CPPUNIT_ASSERT(mesh->isEqual(mesh2,1e-12));
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  mesh->getNodalConnectivity()->setIJ(0,6,12);
  mesh->getNodalConnectivity()->setIJ(0,7,10);
  mesh->getNodalConnectivity()->setIJ(0,8,13);
  mesh->getNodalConnectivity()->setIJ(0,9,9);
  mesh->unPolyze();
  CPPUNIT_ASSERT(!mesh->isEqual(mesh2,1e-12));
  mesh->decrRef();
  mesh2->decrRef();
  // Test for 2D mesh
  mesh=build2DTargetMesh_1();
  mesh2=build2DTargetMesh_1();
  eltsV.resize(5);
  mesh->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  CPPUNIT_ASSERT(!mesh->isEqual(mesh2,1e-12));
  mesh->unPolyze();
  CPPUNIT_ASSERT(mesh->isEqual(mesh2,1e-12));
  mesh->decrRef();
  mesh2->decrRef();
}

void MEDCouplingBasicsTest3::testConvertDegeneratedCells1()
{
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  int conn[32]={0,1,3,3,9,10,12,12, 0,1,3,4,9,9,9,9, 1,1,1,1,10,12,9,10, 10,11,12,9,1,1,1,1};
  mesh->allocateCells(4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+8);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+16);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+24);
  mesh->finishInsertingCells();
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(4,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,mesh->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,mesh->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,mesh->getTypeOfCell(3));
  MEDCouplingFieldDouble *f1=mesh->getMeasureField(true);
  mesh->convertDegeneratedCells();
  mesh->checkConsistencyLight();
  MEDCouplingFieldDouble *f2=mesh->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(4,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_PENTA6,mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_PYRA5,mesh->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TETRA4,mesh->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_PYRA5,mesh->getTypeOfCell(3));
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f1->getArray()->getIJ(0,i),f2->getArray()->getIJ(0,i),1e-5);
  f1->decrRef();
  f2->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest3::testGetNodeIdsNearPoints1()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  DataArrayDouble *coords=mesh->getCoords();
  DataArrayDouble *tmp=DataArrayDouble::New();
  tmp->alloc(3,2);
  const double vals[6]={0.2,0.2,0.1,0.2,0.2,0.2};
  std::copy(vals,vals+6,tmp->getPointer());
  DataArrayDouble *tmp2=DataArrayDouble::Aggregate(coords,tmp);
  tmp->decrRef();
  mesh->setCoords(tmp2);
  tmp2->decrRef();
  const double pts[6]={0.2,0.2,0.1,0.3,-0.3,0.7};
  DataArrayInt *c=mesh->getNodeIdsNearPoint(pts,1e-7);
  CPPUNIT_ASSERT_EQUAL(3,(int)c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(4,c->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(9,c->getIJ(1,0));
  CPPUNIT_ASSERT_EQUAL(11,c->getIJ(2,0));
  c->decrRef();
  DataArrayInt *cI=0;
  mesh->getNodeIdsNearPoints(pts,3,1e-7,c,cI);
  CPPUNIT_ASSERT_EQUAL(4,(int)cI->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(4,(int)c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(4,c->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(9,c->getIJ(1,0));
  CPPUNIT_ASSERT_EQUAL(11,c->getIJ(2,0));
  CPPUNIT_ASSERT_EQUAL(6,c->getIJ(3,0));
  CPPUNIT_ASSERT_EQUAL(0,cI->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(3,cI->getIJ(1,0));
  CPPUNIT_ASSERT_EQUAL(3,cI->getIJ(2,0));
  CPPUNIT_ASSERT_EQUAL(4,cI->getIJ(3,0));
  c->decrRef();
  cI->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest3::testFieldCopyTinyAttrFrom1()
{
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setName("f1");
  f1->setTimeTolerance(1.e-5);
  f1->setDescription("f1Desc");
  f1->setTime(1.23,4,5);
  MEDCouplingFieldDouble *f2=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f2->setName("f2");
  f2->setDescription("f2Desc");
  f2->setTime(6.78,9,10);
  f2->setTimeTolerance(4.556e-12);
  //
  int dt,it;
  f1->copyTinyAttrFrom(f2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.556e-12,f1->getTimeTolerance(),1e-24);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.78,f1->getTime(dt,it),1e-12);
  CPPUNIT_ASSERT_EQUAL(9,dt);
  CPPUNIT_ASSERT_EQUAL(10,it);
  CPPUNIT_ASSERT(std::string(f1->getName())=="f1");//name unchanged
  CPPUNIT_ASSERT(std::string(f1->getDescription())=="f1Desc");//description unchanged
  f1->decrRef();
  f2->decrRef();
  //
  f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setName("f1");
  f1->setTimeTolerance(1.e-5);
  f1->setDescription("f1Desc");
  f2=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f2->setName("f2");
  f2->setDescription("f2Desc");
  f2->setTimeTolerance(4.556e-12);
  //
  f1->copyTinyAttrFrom(f2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.556e-12,f1->getTimeTolerance(),1e-24);
  CPPUNIT_ASSERT(std::string(f1->getName())=="f1");//name unchanged
  CPPUNIT_ASSERT(std::string(f1->getDescription())=="f1Desc");//description unchanged
  f1->decrRef();
  f2->decrRef();
  //
  f1=MEDCouplingFieldDouble::New(ON_CELLS,CONST_ON_TIME_INTERVAL);
  f1->setName("f1");
  f1->setTimeTolerance(1.e-5);
  f1->setDescription("f1Desc");
  f1->setTime(1.23,4,5);
  f1->setEndTime(5.43,2,1);
  f2=MEDCouplingFieldDouble::New(ON_CELLS,CONST_ON_TIME_INTERVAL);
  f2->setName("f2");
  f2->setDescription("f2Desc");
  f2->setTimeTolerance(4.556e-12);
  f2->setTime(6.78,9,10);
  f2->setEndTime(10.98,7,6);
  //
  f1->copyTinyAttrFrom(f2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.556e-12,f1->getTimeTolerance(),1e-24);
  CPPUNIT_ASSERT(std::string(f1->getName())=="f1");//name unchanged
  CPPUNIT_ASSERT(std::string(f1->getDescription())=="f1Desc");//description unchanged
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.78,f1->getTime(dt,it),1e-12);
  CPPUNIT_ASSERT_EQUAL(9,dt);
  CPPUNIT_ASSERT_EQUAL(10,it);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(10.98,f1->getEndTime(dt,it),1e-12);
  CPPUNIT_ASSERT_EQUAL(7,dt);
  CPPUNIT_ASSERT_EQUAL(6,it);
  f1->decrRef();
  f2->decrRef();
  //
  f1=MEDCouplingFieldDouble::New(ON_CELLS,LINEAR_TIME);
  f1->setName("f1");
  f1->setTimeTolerance(1.e-5);
  f1->setDescription("f1Desc");
  f1->setTime(1.23,4,5);
  f1->setEndTime(5.43,2,1);
  f2=MEDCouplingFieldDouble::New(ON_CELLS,LINEAR_TIME);
  f2->setName("f2");
  f2->setDescription("f2Desc");
  f2->setTimeTolerance(4.556e-12);
  f2->setTime(6.78,9,10);
  f2->setEndTime(10.98,7,6);
  //
  f1->copyTinyAttrFrom(f2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4.556e-12,f1->getTimeTolerance(),1e-24);
  CPPUNIT_ASSERT(std::string(f1->getName())=="f1");//name unchanged
  CPPUNIT_ASSERT(std::string(f1->getDescription())=="f1Desc");//description unchanged
  CPPUNIT_ASSERT_DOUBLES_EQUAL(6.78,f1->getTime(dt,it),1e-12);
  CPPUNIT_ASSERT_EQUAL(9,dt);
  CPPUNIT_ASSERT_EQUAL(10,it);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(10.98,f1->getEndTime(dt,it),1e-12);
  CPPUNIT_ASSERT_EQUAL(7,dt);
  CPPUNIT_ASSERT_EQUAL(6,it);
  f1->decrRef();
  f2->decrRef();
}

/*!
 * 1D -> 2D extrusion with rotation
 */
void MEDCouplingBasicsTest3::testExtrudedMesh5()
{
  const double coo1[4]={0.,1.,2.,3.5};
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(4,1);
  std::copy(coo1,coo1+4,a->getPointer());
  MEDCouplingCMesh *b=MEDCouplingCMesh::New();
  b->setCoordsAt(0,a);
  MEDCouplingUMesh *c=b->buildUnstructured();
  CPPUNIT_ASSERT_EQUAL(1,c->getSpaceDimension());
  c->changeSpaceDimension(2);
  //
  DataArrayDouble *d=DataArrayDouble::New();
  d->alloc(13,1);
  d->iota();
  MEDCouplingCMesh *ee=MEDCouplingCMesh::New();
  ee->setCoordsAt(0,d);
  MEDCouplingUMesh *f=ee->buildUnstructured();
  DataArrayDouble *g=f->getCoords()->applyFunc(2,"3.5*IVec+x/6*3.14159265359*JVec");
  CPPUNIT_ASSERT_THROW(f->getCoords()->applyFunc(2,"3.5*IVec+x/6*3.14159265359*KVec"),INTERP_KERNEL::Exception); // KVec refers to component #2 and there is only 2 components !
  DataArrayDouble *h=g->fromPolarToCart();
  f->setCoords(h);
  MEDCouplingUMesh *i=c->buildExtrudedMesh(f,1);
  CPPUNIT_ASSERT_EQUAL(52,i->getNumberOfNodes());
  bool tmp2;
  int tmp3;
  DataArrayInt *tmp=i->mergeNodes(1e-9,tmp2,tmp3);
  CPPUNIT_ASSERT(tmp2);
  CPPUNIT_ASSERT_EQUAL(37,tmp3);
  tmp->decrRef();
  i->convertDegeneratedCells();
  i->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(36,(int)i->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(37,i->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(12,(int)i->getNumberOfCellsWithType(INTERP_KERNEL::NORM_TRI3));
  CPPUNIT_ASSERT_EQUAL(24,(int)i->getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD4));
  const double expected1[3]={0.25,0.75,2.0625};
  MEDCouplingFieldDouble *j=i->getMeasureField(true);
  for(int ii=0;ii<12;ii++)
    for(int k=0;k<3;k++)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[k],j->getIJ(0,ii*3+k),1e-10);
  const double expected2[72]={0.62200846792814113, 0.16666666666681595, 1.4513530918323276, 0.38888888888923495, 2.6293994326053212, 0.7045454545460802, 0.45534180126145435, 0.45534180126150181, 1.0624642029433926, 1.0624642029435025, 1.9248539780597826, 1.9248539780599816, 0.16666666666661334, 0.62200846792815856, 0.38888888888876294, 1.4513530918323678, 0.70454545454522521, 2.629399432605394, -0.16666666666674007, 0.62200846792812436, -0.38888888888906142, 1.4513530918322881, -0.70454545454576778, 2.6293994326052488, -0.45534180126154766, 0.45534180126140844, -1.0624642029436118, 1.0624642029432834, -1.9248539780601803, 1.9248539780595841, -0.62200846792817499, 0.1666666666665495, -1.451353091832408, 0.388888888888613, -2.6293994326054668, 0.70454545454495332, -0.62200846792810593, -0.16666666666680507, -1.451353091832247, -0.38888888888921297, -2.6293994326051746, -0.70454545454604123, -0.45534180126135926, -0.45534180126159562, -1.0624642029431723, -1.0624642029437235, -1.9248539780593836, -1.9248539780603811, -0.1666666666664828, -0.62200846792819242, -0.38888888888846079, -1.4513530918324489, -0.70454545454467987, -2.6293994326055397, 0.16666666666687083, -0.62200846792808862, 0.38888888888936374, -1.4513530918322073, 0.70454545454631357, -2.6293994326051022, 0.45534180126164348, -0.45534180126131207, 1.0624642029438327, -1.0624642029430627, 1.9248539780605791, -1.9248539780591853, 0.62200846792821063, -0.16666666666641802, 1.4513530918324888, -0.38888888888831086, 2.6293994326056125, -0.70454545454440853};
  DataArrayDouble *m=i->computeCellCenterOfMass();
  for(int ii=0;ii<72;ii++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[ii],m->getIJ(0,ii),1e-10);
  //
  m->decrRef();
  j->decrRef();
  i->decrRef();
  h->decrRef();
  g->decrRef();
  f->decrRef();
  ee->decrRef();
  d->decrRef();
  c->decrRef();
  b->decrRef();
  a->decrRef();
}

/*!
 * 1D -> 2D extrusion without rotation
 */
void MEDCouplingBasicsTest3::testExtrudedMesh6()
{
  const double coo1[4]={0.,1.,2.,3.5};
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(4,1);
  std::copy(coo1,coo1+4,a->getPointer());
  MEDCouplingCMesh *b=MEDCouplingCMesh::New();
  b->setCoordsAt(0,a);
  MEDCouplingUMesh *c=b->buildUnstructured();
  CPPUNIT_ASSERT_EQUAL(1,c->getSpaceDimension());
  c->changeSpaceDimension(2);
  //
  DataArrayDouble *d=DataArrayDouble::New();
  d->alloc(5,1);
  d->iota();
  MEDCouplingCMesh *e=MEDCouplingCMesh::New();
  e->setCoordsAt(0,d);
  MEDCouplingUMesh *f=e->buildUnstructured();
  DataArrayDouble *d2=f->getCoords()->applyFunc("x*x/2");
  f->setCoords(d2);
  f->changeSpaceDimension(2);
  //
  const double center[2]={0.,0.};
  f->rotate(center,0,M_PI/3);
  MEDCouplingUMesh *g=c->buildExtrudedMesh(f,0);
  g->checkConsistencyLight();
  const double expected1[]={ 0.4330127018922193, 0.4330127018922193, 0.649519052838329, 1.2990381056766578, 1.299038105676658, 1.948557158514987, 2.1650635094610955, 2.1650635094610964, 3.2475952641916446, 3.031088913245533, 3.0310889132455352, 4.546633369868303 };
  MEDCouplingFieldDouble *f1=g->getMeasureField(true);
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getIJ(0,i),1e-12);
  
  const double expected2[]={0.625, 0.21650635094610962, 1.625, 0.21650635094610959, 2.8750000000000004, 0.21650635094610965, 1.1250000000000002, 1.0825317547305482, 2.125, 1.0825317547305482, 3.3750000000000004, 1.0825317547305484, 2.125, 2.8145825622994254, 3.125, 2.8145825622994254, 4.375, 2.8145825622994254, 3.6250000000000009, 5.4126587736527414, 4.625, 5.4126587736527414, 5.875, 5.4126587736527414};
  DataArrayDouble *f2=g->computeCellCenterOfMass();
  for(int i=0;i<24;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f2->getIJ(0,i),1e-12);
  //
  f1->decrRef();
  f2->decrRef();
  g->decrRef();
  f->decrRef();
  e->decrRef();
  d->decrRef();
  d2->decrRef();
  c->decrRef();
  b->decrRef();
  a->decrRef();
}

/*!
 * 2D -> 3D extrusion with rotation
 */
void MEDCouplingBasicsTest3::testExtrudedMesh7()
{
  const double coo1[4]={0.,1.,2.,3.5};
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(4,1);
  std::copy(coo1,coo1+4,a->getPointer());
  MEDCouplingCMesh *b=MEDCouplingCMesh::New();
  b->setCoordsAt(0,a);
  MEDCouplingUMesh *c=b->buildUnstructured();
  CPPUNIT_ASSERT_EQUAL(1,c->getSpaceDimension());
  c->changeSpaceDimension(2);
  //
  DataArrayDouble *d=DataArrayDouble::New();
  d->alloc(13,1);
  d->iota();
  MEDCouplingCMesh *e=MEDCouplingCMesh::New();
  e->setCoordsAt(0,d);
  MEDCouplingUMesh *f=e->buildUnstructured();
  DataArrayDouble *g=f->getCoords()->applyFunc(2,"3.5*IVec+x/6*3.14159265359*JVec");
  DataArrayDouble *h=g->fromPolarToCart();
  f->setCoords(h);
  MEDCouplingUMesh *i=c->buildExtrudedMesh(f,1);
  CPPUNIT_ASSERT_EQUAL(52,i->getNumberOfNodes());
  bool tmp2;
  int tmp3;
  DataArrayInt *tmp=i->mergeNodes(1e-9,tmp2,tmp3);
  CPPUNIT_ASSERT(tmp2);
  CPPUNIT_ASSERT_EQUAL(37,tmp3);
  tmp->decrRef();
  i->convertDegeneratedCells();
  const double vec1[3]={10.,0.,0.};
  i->translate(vec1);
  DataArrayDouble *g2=h->applyFunc(3,"13.5/3.5*x*IVec+0*JVec+13.5/3.5*y*KVec");
  f->setCoords(g2);
  i->changeSpaceDimension(3);
  MEDCouplingUMesh *i3=i->buildExtrudedMesh(f,1);
  MEDCouplingFieldDouble *f2=i3->getMeasureField(true);
  tmp=i->mergeNodes(1e-9,tmp2,tmp3);
  CPPUNIT_ASSERT(tmp2);
  CPPUNIT_ASSERT_EQUAL(444,tmp3);
  tmp->decrRef();
  const double expected1[36]={1.327751058489274, 4.2942574094314701, 13.024068164857139, 1.3069177251569044, 4.1484240761012954, 12.297505664866796, 1.270833333332571, 3.8958333333309674, 11.039062499993179, 1.2291666666659207, 3.6041666666644425, 9.585937499993932, 1.1930822748415895, 3.3515759238941376, 8.3274943351204556, 1.1722489415082769, 3.2057425905609289, 7.6009318351210622, 1.1722489415082862, 3.2057425905609884, 7.6009318351213713, 1.1930822748416161, 3.3515759238943001, 8.3274943351212727, 1.2291666666659564, 3.6041666666646734, 9.5859374999950777, 1.2708333333326081, 3.8958333333311868, 11.039062499994293, 1.3069177251569224, 4.1484240761014384, 12.297505664867627, 1.3277510584902354, 4.2942574094346071, 13.024068164866796};
  int kk=0;
  for(int ii=0;ii<12;ii++)
    for(int jj=0;jj<36;jj++,kk++)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[jj],f2->getIJ(0,kk),1e-9);
  //
  f2->decrRef();
  i3->decrRef();
  g2->decrRef();
  i->decrRef();
  h->decrRef();
  g->decrRef();
  f->decrRef();
  e->decrRef();
  d->decrRef();
  c->decrRef();
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest3::testSimplexize1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  std::vector<int> v(1);
  v[0]=3;
  m->convertToPolyTypes(&v[0],&v[0]+v.size());
  DataArrayInt *da=m->simplexize(0);
  CPPUNIT_ASSERT_EQUAL(7,(int)da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfComponents());
  const int expected2[7]={0,0,1,2,3,4,4};
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],da->getIJ(i,0));
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(7,(int)m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(3));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_POLYGON,m->getTypeOfCell(4));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(5));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(6));
  const double expected1[7]={0.125,0.125,0.125,0.125,0.25,0.125,0.125};
  MEDCouplingFieldDouble *f=m->getMeasureField(false);
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i]*sqrt(2.),f->getIJ(i,0),1e-10);
  std::set<INTERP_KERNEL::NormalizedCellType> types=m->getAllGeoTypes();
  CPPUNIT_ASSERT_EQUAL(2,(int)types.size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*(types.begin()));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_POLYGON,*(++(types.begin())));
  f->decrRef();
  da->decrRef();
  m->decrRef();
  //
  m=build3DSurfTargetMesh_1();
  v[0]=3;
  m->convertToPolyTypes(&v[0],&v[0]+v.size());
  da=m->simplexize(1);
  CPPUNIT_ASSERT_EQUAL(7,(int)da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfComponents());
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],da->getIJ(i,0));
  m->checkConsistencyLight();
  types=m->getAllGeoTypes();
  CPPUNIT_ASSERT_EQUAL(2,(int)types.size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*(types.begin()));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_POLYGON,*(++(types.begin())));
  CPPUNIT_ASSERT_EQUAL(7,(int)m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(3));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_POLYGON,m->getTypeOfCell(4));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(5));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,m->getTypeOfCell(6));
  f=m->getMeasureField(false);
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i]*sqrt(2.),f->getIJ(i,0),1e-10);
  f->decrRef();
  da->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest3::testSimplexize2()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  std::vector<int> v(1);
  v[0]=3;
  m->convertToPolyTypes(&v[0],&v[0]+v.size());
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setMesh(m);
  DataArrayDouble *arr=DataArrayDouble::New();
  const double arr1[10]={10.,110.,20.,120.,30.,130.,40.,140.,50.,150.};
  arr->alloc(5,2);
  std::copy(arr1,arr1+10,arr->getPointer());
  f1->setArray(arr);
  arr->decrRef();
  //
  f1->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->simplexize(0));
  f1->checkConsistencyLight();
  const double expected1[14]={10.,110.,10.,110.,20.,120.,30.,130.,40.,140.,50.,150.,50.,150.};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getIJ(0,i),1e-10);
  CPPUNIT_ASSERT(!f1->simplexize(0));
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getIJ(0,i),1e-10);
  //
  f1->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest3::testDAMeld1()
{
  DataArrayDouble *da1=DataArrayDouble::New();
  da1->alloc(7,2);
  DataArrayDouble *da2=DataArrayDouble::New();
  da2->alloc(7,1);
  //
  da1->fillWithValue(7.);
  da2->iota(0.);
  MCAuto<DataArrayDouble> da3=da2->applyFunc(3,"10*x*IVec+100*x*JVec+1000*x*KVec");
  //
  da1->setInfoOnComponent(0,"c0da1");
  da1->setInfoOnComponent(1,"c1da1");
  da3->setInfoOnComponent(0,"c0da3");
  da3->setInfoOnComponent(1,"c1da3");
  da3->setInfoOnComponent(2,"c2da3");
  //
  DataArrayDouble *da1C=da1->deepCopy();
  da1->meldWith(da3);
  CPPUNIT_ASSERT_EQUAL(5,(int)da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,(int)da1->getNumberOfTuples());
  CPPUNIT_ASSERT(da1->getInfoOnComponent(0)=="c0da1");
  CPPUNIT_ASSERT(da1->getInfoOnComponent(1)=="c1da1");
  CPPUNIT_ASSERT(da1->getInfoOnComponent(2)=="c0da3");
  CPPUNIT_ASSERT(da1->getInfoOnComponent(3)=="c1da3");
  CPPUNIT_ASSERT(da1->getInfoOnComponent(4)=="c2da3");
  //
  const double expected1[35]={7.,7.,0.,0.,0., 7.,7.,10.,100.,1000., 7.,7.,20.,200.,2000., 7.,7.,30.,300.,3000., 7.,7.,40.,400.,4000.,7.,7.,50.,500.,5000.,7.,7.,60.,600.,6000.};
  for(int i=0;i<35;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da1->getIJ(0,i),1e-10);
  //
  MCAuto<DataArrayInt> dai1=da1C->convertToIntArr();
  MCAuto<DataArrayInt> dai3=da3->convertToIntArr();
  dai1->meldWith(dai3);
  CPPUNIT_ASSERT_EQUAL(5,(int)dai1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,(int)dai1->getNumberOfTuples());
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(0)=="c0da1");
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(1)=="c1da1");
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(2)=="c0da3");
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(3)=="c1da3");
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(4)=="c2da3");
  for(int i=0;i<35;i++)
    CPPUNIT_ASSERT_EQUAL((int)expected1[i],dai1->getIJ(0,i));
  // test of static method DataArrayDouble::meld
  DataArrayDouble *da4=DataArrayDouble::Meld(da1C,da3);
  CPPUNIT_ASSERT_EQUAL(5,(int)da4->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,(int)da4->getNumberOfTuples());
  CPPUNIT_ASSERT(da4->getInfoOnComponent(0)=="c0da1");
  CPPUNIT_ASSERT(da4->getInfoOnComponent(1)=="c1da1");
  CPPUNIT_ASSERT(da4->getInfoOnComponent(2)=="c0da3");
  CPPUNIT_ASSERT(da4->getInfoOnComponent(3)=="c1da3");
  CPPUNIT_ASSERT(da4->getInfoOnComponent(4)=="c2da3");
  for(int i=0;i<35;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da4->getIJ(0,i),1e-10);
  // test of static method DataArrayInt::meld
  dai1=da1C->convertToIntArr();
  DataArrayInt *dai4=DataArrayInt::Meld(dai1,dai3);
  CPPUNIT_ASSERT_EQUAL(5,(int)dai4->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,(int)dai4->getNumberOfTuples());
  CPPUNIT_ASSERT(dai4->getInfoOnComponent(0)=="c0da1");
  CPPUNIT_ASSERT(dai4->getInfoOnComponent(1)=="c1da1");
  CPPUNIT_ASSERT(dai4->getInfoOnComponent(2)=="c0da3");
  CPPUNIT_ASSERT(dai4->getInfoOnComponent(3)=="c1da3");
  CPPUNIT_ASSERT(dai4->getInfoOnComponent(4)=="c2da3");
  for(int i=0;i<35;i++)
    CPPUNIT_ASSERT_EQUAL((int)expected1[i],dai4->getIJ(0,i));
  //
  dai4->decrRef();
  da4->decrRef();
  da1C->decrRef();
  da1->decrRef();
  da2->decrRef();
}

void MEDCouplingBasicsTest3::testFieldMeld1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setMesh(m);
  DataArrayDouble *da1=DataArrayDouble::New();
  const double arr1[5]={12.,23.,34.,45.,56.};
  da1->alloc(5,1);
  std::copy(arr1,arr1+5,da1->getPointer());
  da1->setInfoOnComponent(0,"aaa");
  f1->setArray(da1);
  f1->setTime(3.4,2,1);
  f1->checkConsistencyLight();
  //
  MEDCouplingFieldDouble *f2=f1->deepCopy();
  f2->setMesh(f1->getMesh());
  f2->checkConsistencyLight();
  f2->changeNbOfComponents(2,5.);
  (*f2)=5.;
  f2->getArray()->setInfoOnComponent(0,"bbb");
  f2->getArray()->setInfoOnComponent(1,"ccc");
  f2->checkConsistencyLight();
  //
  MEDCouplingFieldDouble *f3=MEDCouplingFieldDouble::MeldFields(f2,f1);
  f3->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(5,(int)f3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)f3->getNumberOfComponents());
  CPPUNIT_ASSERT(f3->getArray()->getInfoOnComponent(0)=="bbb");
  CPPUNIT_ASSERT(f3->getArray()->getInfoOnComponent(1)=="ccc");
  CPPUNIT_ASSERT(f3->getArray()->getInfoOnComponent(2)=="aaa");
  const double expected1[15]={5.,5.,12.,5.,5.,23.,5.,5.,34.,5.,5.,45.,5.,5.,56.};
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f3->getIJ(0,i),1e-12);
  int dt,it;
  double time=f3->getTime(dt,it);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.4,time,1e-14);
  CPPUNIT_ASSERT_EQUAL(2,dt);
  CPPUNIT_ASSERT_EQUAL(1,it);
  //
  MEDCouplingFieldDouble *f4=f2->buildNewTimeReprFromThis(NO_TIME,false);
  MEDCouplingFieldDouble *f5=f1->buildNewTimeReprFromThis(NO_TIME,false);
  MEDCouplingFieldDouble *f6=MEDCouplingFieldDouble::MeldFields(f4,f5);
  f6->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(5,(int)f6->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,(int)f6->getNumberOfComponents());
  CPPUNIT_ASSERT(f6->getArray()->getInfoOnComponent(0)=="bbb");
  CPPUNIT_ASSERT(f6->getArray()->getInfoOnComponent(1)=="ccc");
  CPPUNIT_ASSERT(f6->getArray()->getInfoOnComponent(2)=="aaa");
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f6->getIJ(0,i),1e-12);
  //
  f6->decrRef();
  f4->decrRef();
  f5->decrRef();
  f3->decrRef();
  da1->decrRef();
  f2->decrRef();
  f1->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest3::testMergeNodes2()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  const double vec[2]={0.002,0.};
  m2->translate(vec);
  //
  std::vector<const MEDCouplingUMesh *> tmp(2);
  tmp[0]=m1;
  tmp[1]=m2;
  MEDCouplingUMesh *m3=MEDCouplingUMesh::MergeUMeshes(tmp);
  bool b;
  int newNbOfNodes;
  DataArrayInt *da=m3->mergeNodesCenter(0.01,b,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(9,m3->getNumberOfNodes());
  const double expected1[18]={-0.299,-0.3, 0.201,-0.3, 0.701,-0.3, -0.299,0.2, 0.201,0.2, 0.701,0.2, -0.299,0.7, 0.201,0.7, 0.701,0.7};
  for(int i=0;i<18;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],m3->getCoords()->getIJ(0,i),1e-13);
  //
  da->decrRef();
  m3->decrRef();
  m1->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest3::testMergeField2()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setMesh(m);
  DataArrayDouble *arr=DataArrayDouble::New();
  arr->alloc(5,2);
  arr->fillWithValue(2.);
  f1->setArray(arr);
  arr->decrRef();
  MEDCouplingFieldDouble *f2=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f2->setMesh(m);
  arr=DataArrayDouble::New();
  arr->alloc(5,2);
  arr->fillWithValue(5.);
  f2->setArray(arr);
  arr->decrRef();
  MEDCouplingFieldDouble *f3=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f3->setMesh(m);
  arr=DataArrayDouble::New();
  arr->alloc(5,2);
  arr->fillWithValue(7.);
  f3->setArray(arr);
  arr->decrRef();
  //
  std::vector<const MEDCouplingFieldDouble *> tmp(3);
  tmp[0]=f1; tmp[1]=f2; tmp[2]=f3;
  MEDCouplingFieldDouble *f4=MEDCouplingFieldDouble::MergeFields(tmp);
  CPPUNIT_ASSERT_EQUAL(15,(int)f4->getMesh()->getNumberOfCells());
  const double expected1[30]={2.,2.,2.,2.,2.,2.,2.,2.,2.,2., 5.,5.,5.,5.,5.,5.,5.,5.,5.,5., 7.,7.,7.,7.,7.,7.,7.,7.,7.,7.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f4->getIJ(0,i),1.e-13);
  //
  f4->decrRef();
  f1->decrRef();
  f2->decrRef();
  f3->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest3::testDAIBuildComplement1()
{
  DataArrayInt *a=DataArrayInt::New();
  const int tab[4]={3,1,7,8};
  a->alloc(4,1);
  std::copy(tab,tab+4,a->getPointer());
  DataArrayInt *b=a->buildComplement(12);
  CPPUNIT_ASSERT_EQUAL(8,(int)b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)b->getNumberOfComponents());
  const int expected1[8]={0,2,4,5,6,9,10,11};
  for(int i=0;i<8;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],b->getIJ(0,i));
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest3::testDAIBuildUnion1()
{
  DataArrayInt *a=DataArrayInt::New();
  const int tab1[4]={3,1,7,8};
  a->alloc(4,1);
  std::copy(tab1,tab1+4,a->getPointer());
  DataArrayInt *c=DataArrayInt::New();
  const int tab2[5]={5,3,0,18,8};
  c->alloc(5,1);
  std::copy(tab2,tab2+5,c->getPointer());
  DataArrayInt *b=a->buildUnion(c);
  CPPUNIT_ASSERT_EQUAL(7,(int)b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)b->getNumberOfComponents());
  const int expected1[7]={0,1,3,5,7,8,18};
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],b->getIJ(0,i));
  c->decrRef();
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest3::testDAIBuildIntersection1()
{
  DataArrayInt *a=DataArrayInt::New();
  const int tab1[4]={3,1,7,8};
  a->alloc(4,1);
  std::copy(tab1,tab1+4,a->getPointer());
  DataArrayInt *c=DataArrayInt::New();
  const int tab2[5]={5,3,0,18,8};
  c->alloc(5,1);
  std::copy(tab2,tab2+5,c->getPointer());
  DataArrayInt *b=a->buildIntersection(c);
  CPPUNIT_ASSERT_EQUAL(2,(int)b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)b->getNumberOfComponents());
  const int expected1[2]={3,8};
  for(int i=0;i<2;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],b->getIJ(0,i));
  c->decrRef();
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest3::testDAIDeltaShiftIndex1()
{
  DataArrayInt *a=DataArrayInt::New();
  const int tab[7]={1,3,6,7,7,9,15};
  a->alloc(7,1);
  std::copy(tab,tab+7,a->getPointer());
  DataArrayInt *b=a->deltaShiftIndex();
  CPPUNIT_ASSERT_EQUAL(6,(int)b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)b->getNumberOfComponents());
  const int expected1[6]={2,3,1,0,2,6};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],b->getIJ(0,i));
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest3::testDaDoubleSelectByTupleIdSafe1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(7,2);
  a->setInfoOnComponent(0,"toto");
  a->setInfoOnComponent(1,"tata");
  const double arr1[14]={1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1};
  std::copy(arr1,arr1+14,a->getPointer());
  //
  const int arr2[7]={4,2,0,6,5};
  DataArrayDouble *b=a->selectByTupleIdSafe(arr2,arr2+5);
  CPPUNIT_ASSERT_EQUAL(5,(int)b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)b->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(b->getInfoOnComponent(1))=="tata");
  const double expected1[10]={5.1,15.1,3.1,13.1,1.1,11.1,7.1,17.1,6.1,16.1};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],b->getIJ(0,i),1e-14);
  const int arr4[5]={4,-1,0,6,5};
  CPPUNIT_ASSERT_THROW(a->selectByTupleIdSafe(arr4,arr4+5),INTERP_KERNEL::Exception);
  const int arr5[5]={4,2,0,6,7};
  CPPUNIT_ASSERT_THROW(a->selectByTupleIdSafe(arr5,arr5+5),INTERP_KERNEL::Exception);
  b->decrRef();
  a->decrRef();
  //
  DataArrayInt *c=DataArrayInt::New();
  c->alloc(7,2);
  c->setInfoOnComponent(0,"toto");
  c->setInfoOnComponent(1,"tata");
  const int arr3[14]={1,11,2,12,3,13,4,14,5,15,6,16,7,17};
  std::copy(arr3,arr3+14,c->getPointer());
  DataArrayInt *d=c->selectByTupleIdSafe(arr2,arr2+5);
  CPPUNIT_ASSERT_EQUAL(5,(int)d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(1))=="tata");
  const int expected2[10]={5,15,3,13,1,11,7,17,6,16};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d->getIJ(0,i));
  CPPUNIT_ASSERT_THROW(c->selectByTupleIdSafe(arr4,arr4+5),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(c->selectByTupleIdSafe(arr5,arr5+5),INTERP_KERNEL::Exception);
  c->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest3::testAreCellsIncludedIn1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  const int pt[2]={1,3};
  MEDCouplingUMesh *m2=(MEDCouplingUMesh *)m->buildPartOfMySelf(pt,pt+2,true);
  DataArrayInt *tmp;
  CPPUNIT_ASSERT(m->areCellsIncludedIn(m2,0,tmp));
  CPPUNIT_ASSERT_EQUAL(2,(int)tmp->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)tmp->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(pt[0],tmp->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(pt[1],tmp->getIJ(0,1));
  tmp->decrRef();
  CPPUNIT_ASSERT(!m2->areCellsIncludedIn(m,0,tmp));
  tmp->decrRef();
  m2->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest3::testDAIBuildSubstraction1()
{
  DataArrayInt *a=DataArrayInt::New();
  const int aa[]={2,3,6,8,9};
  a->alloc(5,1);
  std::copy(aa,aa+5,a->getPointer());
  DataArrayInt *b=DataArrayInt::New();
  const int bb[]={1,3,5,9,11};
  b->alloc(5,1);
  std::copy(bb,bb+5,b->getPointer());
  //
  DataArrayInt *c=a->buildSubstraction(b);
  CPPUNIT_ASSERT_EQUAL(3,(int)c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)c->getNumberOfComponents());
  const int expected1[3]={2,6,8};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,c->getConstPointer()));
  //
  c->decrRef();
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest3::testBuildOrthogonalField2()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  DataArrayInt *d1=DataArrayInt::New();
  DataArrayInt *d2=DataArrayInt::New();
  DataArrayInt *d3=DataArrayInt::New();
  DataArrayInt *d4=DataArrayInt::New();
  MEDCouplingUMesh *m1=m->buildDescendingConnectivity(d1,d2,d3,d4);
  //
  MEDCouplingFieldDouble *f1=m1->buildOrthogonalField();
  DataArrayDouble *da1=f1->getArray();
  CPPUNIT_ASSERT_EQUAL(2,(int)da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(13,(int)da1->getNumberOfTuples());
  //
  const double expected1[26]={-1.,0.,0.,1.,1.,0.,0.,-1.,0.707106781186548,0.707106781186548,0.,-1.,0.,1.,1.,0.,0.,1.,1.,0.,-1.,0.,0.,1.,1.,0.};
  for(int i=0;i<26;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da1->getIJ(0,i),1e-14);
  //
  f1->decrRef();
  m1->decrRef();
  d1->decrRef();
  d2->decrRef();
  d3->decrRef();
  d4->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest3::testUMInsertNextCell1()
{
  double targetCoords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  int targetConn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->allocateCells(5);
  CPPUNIT_ASSERT_THROW(targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn),INTERP_KERNEL::Exception);
  targetMesh->setMeshDimension(2);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  CPPUNIT_ASSERT_THROW(targetMesh->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,targetConn),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(targetMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,targetConn),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT1,1,targetConn),INTERP_KERNEL::Exception);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+7);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+10);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+14);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,2);
  std::copy(targetCoords,targetCoords+18,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  targetMesh->checkConsistencyLight();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest3::testFieldOperatorDivDiffComp1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  DataArrayInt *d1=DataArrayInt::New();
  DataArrayInt *d2=DataArrayInt::New();
  DataArrayInt *d3=DataArrayInt::New();
  DataArrayInt *d4=DataArrayInt::New();
  MEDCouplingUMesh *m1=m->buildDescendingConnectivity(d1,d2,d3,d4);
  //
  MEDCouplingFieldDouble *f1=m1->buildOrthogonalField();
  const double arr1[13]={2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.};
  DataArrayDouble *arr=DataArrayDouble::New();
  arr->alloc(13,1);
  std::copy(arr1,arr1+13,arr->getPointer());
  MEDCouplingFieldDouble *f2=MEDCouplingFieldDouble::New(ON_CELLS);
  f2->setArray(arr);
  f2->setMesh(m1);
  f2->checkConsistencyLight();
  //
  MEDCouplingFieldDouble *f3=(*f1)/(*f2);
  CPPUNIT_ASSERT_THROW((*f2)/(*f1),INTERP_KERNEL::Exception);
  f3->checkConsistencyLight();
  (*f1)/=(*f2);
  CPPUNIT_ASSERT(f1->isEqual(f3,1e-10,1e-10));
  CPPUNIT_ASSERT_THROW((*f2)/=(*f1),INTERP_KERNEL::Exception);
  const double expected1[26]={-0.5, 0.0, 0.0, 0.33333333333333331, 0.25, 0.0, 0.0, -0.20000000000000001, 0.117851130197758, 0.117851130197758, 0.0, -0.14285714285714285, 0.0, 0.125, 0.1111111111111111, 0.0, 0.0, 0.10000000000000001, 0.090909090909090912, 0.0, -0.083333333333333329, 0.0, 0.0, 0.076923076923076927, 0.071428571428571425, 0.0};
  for(int i=0;i<26;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f3->getIJ(0,i),1e-10);
  //
  f3->decrRef();
  f2->decrRef();
  arr->decrRef();
  f1->decrRef();
  m1->decrRef();
  d1->decrRef();
  d2->decrRef();
  d3->decrRef();
  d4->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest3::testDARearrange1()
{
  DataArrayInt *da1=DataArrayInt::New();
  da1->alloc(12,1);
  da1->iota(0);
  const int *ptr=da1->getConstPointer();
  //
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(1,(int)da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(12,(int)da1->getNumberOfTuples());
  da1->rearrange(4);
  CPPUNIT_ASSERT(ptr==da1->getConstPointer());
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(4,(int)da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,(int)da1->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(i,da1->getIJ(0,i));
  //
  da1->rearrange(6);
  CPPUNIT_ASSERT(ptr==da1->getConstPointer());
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(6,(int)da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(2,(int)da1->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(i,da1->getIJ(0,i));
  //
  CPPUNIT_ASSERT_THROW(da1->rearrange(7),INTERP_KERNEL::Exception);
  //
  da1->rearrange(12);
  CPPUNIT_ASSERT(ptr==da1->getConstPointer());
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(12,(int)da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(1,(int)da1->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(i,da1->getIJ(0,i));
  //
  da1->rearrange(3);
  CPPUNIT_ASSERT(ptr==da1->getConstPointer());
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(3,(int)da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(4,(int)da1->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(i,da1->getIJ(0,i));
  //double
  MCAuto<DataArrayDouble> da2=da1->convertToDblArr();
  da1->decrRef();
  const double *ptr2=da2->getConstPointer();
  //
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(3,(int)da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(4,(int)da2->getNumberOfTuples());
  da2->rearrange(4);
  CPPUNIT_ASSERT(ptr2==da2->getConstPointer());
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(4,(int)da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,(int)da2->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)i,da2->getIJ(0,i),1e-14);
  //
  da2->rearrange(6);
  CPPUNIT_ASSERT(ptr2==da2->getConstPointer());
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(6,(int)da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(2,(int)da2->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)i,da2->getIJ(0,i),1e-14);
  //
  CPPUNIT_ASSERT_THROW(da2->rearrange(7),INTERP_KERNEL::Exception);
  //
  da2->rearrange(1);
  CPPUNIT_ASSERT(ptr2==da2->getConstPointer());
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(12,(int)da2->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)i,da2->getIJ(0,i),1e-14);
  //
  da2->rearrange(3);
  CPPUNIT_ASSERT(ptr2==da2->getConstPointer());
  CPPUNIT_ASSERT_EQUAL((std::size_t)12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(3,(int)da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(4,(int)da2->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)i,da2->getIJ(0,i),1e-14);
}

void MEDCouplingBasicsTest3::testGetDifferentValues1()
{
  DataArrayInt *da1=DataArrayInt::New();
  const int arr[12]={1,2,3,2,2,3,5,1,5,5,2,2};
  da1->alloc(4,3);
  std::copy(arr,arr+12,da1->getPointer());
  DataArrayInt *s=da1->getDifferentValues();
  const int expected1[4]={1,2,3,5};
  CPPUNIT_ASSERT_EQUAL(4,(int)s->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+4,s->begin()));
  da1->decrRef();
  s->decrRef();
}

void MEDCouplingBasicsTest3::testDAIBuildPermutationArr1()
{
  DataArrayInt *a=DataArrayInt::New();
  const int vala[5]={4,5,6,7,8};
  a->alloc(5,1);
  std::copy(vala,vala+5,a->getPointer());
  DataArrayInt *b=DataArrayInt::New();
  const int valb[5]={5,4,8,6,7};
  b->alloc(5,1);
  std::copy(valb,valb+5,b->getPointer());
  DataArrayInt *c=a->buildPermutationArr(*b);
  const int expect1[5]={1,0,4,2,3};
  CPPUNIT_ASSERT_EQUAL(5,(int)c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)c->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expect1,expect1+5,c->getConstPointer()));
  CPPUNIT_ASSERT(a->isEqualWithoutConsideringStrAndOrder(*b));
  b->setIJ(0,0,9);
  CPPUNIT_ASSERT(!a->isEqualWithoutConsideringStrAndOrder(*b));
  CPPUNIT_ASSERT_THROW(a->buildPermutationArr(*b),INTERP_KERNEL::Exception);
  a->setIJ(3,0,4);
  b->setIJ(0,0,5);
  b->setIJ(4,0,4);//;a==[4,5,6,4,8] and b==[5,4,8,6,4]
  CPPUNIT_ASSERT(a->isEqualWithoutConsideringStrAndOrder(*b));
  c->decrRef();
  c=a->buildPermutationArr(*b);
  const int expect2[5]={1,3,4,2,3};
  CPPUNIT_ASSERT(std::equal(expect2,expect2+5,c->getConstPointer()));
  MCAuto<DataArrayDouble> d=b->convertToDblArr();
  b->sort();
  const int expect3[5]={4,4,5,6,8};
  CPPUNIT_ASSERT(std::equal(expect3,expect3+5,b->getConstPointer()));
  d->sort();
  CPPUNIT_ASSERT_EQUAL(5,(int)d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(double(expect3[i]),d->getIJ(i,0),1e-14);
  //
  b->decrRef();
  c->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest3::testAreCellsIncludedIn2()
{
  const char myName[]="Vitoo";
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  MEDCouplingUMesh *m2=(MEDCouplingUMesh *)m->buildPartOfMySelf(0,0,true);
  CPPUNIT_ASSERT_EQUAL(0,(int)m2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(3,m2->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,m2->getMeshDimension());
  m2->setName(myName);
  DataArrayInt *tmp;
  CPPUNIT_ASSERT(m->areCellsIncludedIn(m2,0,tmp));
  CPPUNIT_ASSERT(std::string(myName)==tmp->getName());
  CPPUNIT_ASSERT_EQUAL(0,(int)tmp->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)tmp->getNumberOfComponents());
  m->decrRef();
  m2->decrRef();
  tmp->decrRef();
}

void MEDCouplingBasicsTest3::testUMeshGetPartBarycenterAndOwner1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  const int part[3]={1,0,4};
  DataArrayDouble *b=m1->getPartBarycenterAndOwner(part,part+3);
  CPPUNIT_ASSERT_EQUAL(2,(int)b->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,(int)b->getNumberOfTuples());
  const double expected1[6]={0.36666666666666665,-0.13333333333333333,-0.05,-0.05,0.45,0.45};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],b->getIJ(0,i),1e-14);
  b->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest3::testUMeshGetPartMeasureField1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  const int part[3]={1,0,4};
  DataArrayDouble *b=m1->getPartMeasureField(true,part,part+3);
  CPPUNIT_ASSERT_EQUAL(1,(int)b->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,(int)b->getNumberOfTuples());
  const double expected1[3]={0.125,0.25,0.25};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],b->getIJ(0,i),1e-14);
  b->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest3::testUMeshBuildPartOrthogonalField1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  m1->changeSpaceDimension(3);
  const int part[3]={1,0,4};
  MEDCouplingFieldDouble *b=m1->buildPartOrthogonalField(part,part+3);
  CPPUNIT_ASSERT_EQUAL(3,(int)b->getArray()->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,(int)b->getArray()->getNumberOfTuples());
  const double expected1[9]={0.,0.,-1.,0.,0.,-1.,0.,0.,-1.};
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],b->getArray()->getIJ(0,i),1e-14);
  b->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest3::testUMeshGetTypesOfPart1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  const int part1[]={0,3,4};
  std::set<INTERP_KERNEL::NormalizedCellType> s;
  s=m1->getTypesOfPart(part1,part1+3);
  CPPUNIT_ASSERT(s.size()==1);
  CPPUNIT_ASSERT(*s.begin()==INTERP_KERNEL::NORM_QUAD4);
  const int part2[]={2,2,2,1};
  s=m1->getTypesOfPart(part2,part2+4);
  CPPUNIT_ASSERT(s.size()==1);
  CPPUNIT_ASSERT(*s.begin()==INTERP_KERNEL::NORM_TRI3);
  const int part3[]={3,2,1};
  s=m1->getTypesOfPart(part3,part3+3);
  CPPUNIT_ASSERT(s.size()==2);
  CPPUNIT_ASSERT(*s.begin()==INTERP_KERNEL::NORM_TRI3);
  CPPUNIT_ASSERT(*(++s.begin())==INTERP_KERNEL::NORM_QUAD4);
  m1->decrRef();
}

void MEDCouplingBasicsTest3::testUMeshKeepCellIdsByType1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  const int part1[3]={0,3,4};
  DataArrayInt *a=m1->keepCellIdsByType(INTERP_KERNEL::NORM_TRI3,part1,part1+3);
  CPPUNIT_ASSERT_EQUAL(1,(int)a->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(0,(int)a->getNumberOfTuples());
  a->decrRef();
  //
  const int part2[5]={3,2,0,2,4};
  a=m1->keepCellIdsByType(INTERP_KERNEL::NORM_TRI3,part2,part2+5);
  CPPUNIT_ASSERT_EQUAL(1,(int)a->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(2,(int)a->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,a->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(2,a->getIJ(1,0));
  a->decrRef();
  //
  a=m1->keepCellIdsByType(INTERP_KERNEL::NORM_QUAD4,part2,part2+5);
  CPPUNIT_ASSERT_EQUAL(1,(int)a->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,(int)a->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,a->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(0,a->getIJ(1,0));
  CPPUNIT_ASSERT_EQUAL(4,a->getIJ(2,0));
  //
  a->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest3::testDAIAggregateMulti1()
{
  DataArrayInt *a=DataArrayInt::New();
  a->setName("aa");
  a->alloc(4,1);
  a->iota(0);
  a->rearrange(2);
  DataArrayInt *b=DataArrayInt::New();
  b->setName("bb");
  b->alloc(6,1);
  b->iota(0);
  b->rearrange(2);
  //
  std::vector<const DataArrayInt *> v(2);
  v[0]=a; v[1]=b;
  DataArrayInt *c=DataArrayInt::Aggregate(v);
  CPPUNIT_ASSERT_EQUAL(5,(int)c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,(int)c->getNumberOfComponents());
  CPPUNIT_ASSERT(c->getName()=="aa");
  const int expect1[10]={0,1,2,3,0,1,2,3,4,5};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(expect1[i],c->getIJ(0,i));
  //
  c->decrRef();
  a->decrRef();
  b->decrRef();
}

void MEDCouplingBasicsTest3::testMergeUMeshes2()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  MEDCouplingUMesh *m2=build3DSurfTargetMesh_1();
  MEDCouplingUMesh *m3=build3DSurfTargetMesh_1();
  //
  const int vec1[3]={0,2,3};
  MEDCouplingUMesh *m2_2=(MEDCouplingUMesh *)m2->buildPartOfMySelf(vec1,vec1+3,false);
  const int vec2[2]={1,1};
  MEDCouplingUMesh *m3_2=(MEDCouplingUMesh *)m3->buildPartOfMySelf(vec2,vec2+2,false);
  //
  std::vector<const MEDCouplingUMesh *> ms(3);
  std::vector<const MEDCouplingMesh *> ms2(3);
  ms[0]=m1; ms[1]=m2_2; ms[2]=m3_2;
  ms2[0]=m1; ms2[1]=m2_2; ms2[2]=m3_2;
  //
  MEDCouplingUMesh *m4=MEDCouplingUMesh::MergeUMeshes(ms);
  m4->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(10,(int)m4->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(20,m4->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(45,m4->getNodalConnectivityArrayLen());
  //
  MEDCouplingMesh *m4bis=MEDCouplingMesh::MergeMeshes(ms2);
  CPPUNIT_ASSERT(m4->isEqual(m4bis,1e-12));
  m4bis->decrRef();
  //
  const int vec3[5]={0,1,2,3,4};
  MEDCouplingUMesh *m4_1=(MEDCouplingUMesh *)m4->buildPartOfMySelf(vec3,vec3+5,false);
  m4_1->setName(m1->getName().c_str());
  CPPUNIT_ASSERT(m4_1->isEqual(m1,1e-12));
  m4_1->decrRef();
  //
  const int vec4[3]={5,6,7};
  MEDCouplingUMesh *m4_2=(MEDCouplingUMesh *)m4->buildPartOfMySelf(vec4,vec4+3,false);
  DataArrayInt *cellCor=0;
  DataArrayInt *nodeCor=0;
  m4_2->checkGeoEquivalWith(m2_2,10,1e-12,cellCor,nodeCor);
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  m4_2->decrRef();
  //
  const int vec5[2]={8,9};
  MEDCouplingUMesh *m4_3=(MEDCouplingUMesh *)m4->buildPartOfMySelf(vec5,vec5+2,false);
  CPPUNIT_ASSERT_EQUAL(2,(int)m4_3->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(3,m4_3->getNumberOfNodes());
  m3_2->zipCoords();
  m4_3->setName(m3_2->getName().c_str());
  CPPUNIT_ASSERT(m4_3->isEqual(m3_2,1e-12));
  m4_3->decrRef();
  //
  m4->decrRef();
  m1->decrRef();
  m2->decrRef();
  m2_2->decrRef();
  m3->decrRef();
  m3_2->decrRef();
}

void MEDCouplingBasicsTest3::testBuild0DMeshFromCoords1()
{
  const double sourceCoords[12]={-0.3,-0.3,0., 0.7,-0.3,0., -0.3,0.7,0., 0.7,0.7,0.};
  DataArrayDouble *coo=DataArrayDouble::New();
  coo->alloc(4,3);
  coo->setName("My0D");
  std::copy(sourceCoords,sourceCoords+12,coo->getPointer());
  MEDCouplingUMesh *m=MEDCouplingUMesh::Build0DMeshFromCoords(coo);
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(4,m->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(4,(int)m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(3,m->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(0,m->getMeshDimension());
  std::set<INTERP_KERNEL::NormalizedCellType> types=m->getAllGeoTypes();
  CPPUNIT_ASSERT_EQUAL(1,(int)types.size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_POINT1,*types.begin());
  for(int i=0;i<4;i++)
    {
      std::vector<int> conn;
      m->getNodeIdsOfCell(i,conn);
      CPPUNIT_ASSERT_EQUAL(1,(int)conn.size());
      CPPUNIT_ASSERT_EQUAL(i,conn[0]);
      CPPUNIT_ASSERT(INTERP_KERNEL::NORM_POINT1==m->getTypeOfCell(i));
    }
  CPPUNIT_ASSERT(std::string(m->getName())=="My0D");
  m->decrRef();
  coo->decrRef();
}

