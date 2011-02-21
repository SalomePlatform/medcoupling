//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#include "MEDCouplingBasicsTest.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingFieldOverTime.hxx"

#include <cmath>
#include <functional>
#include <iterator>

using namespace ParaMEDMEM;

void MEDCouplingBasicsTest::testDescriptionInMeshTimeUnit1()
{
  static const char text1[]="totoTTEDD";
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  m->setDescription(text1);
  CPPUNIT_ASSERT(std::string(m->getDescription())==text1);
  MEDCouplingUMesh *m2=(MEDCouplingUMesh *)m->deepCpy();
  CPPUNIT_ASSERT(m->isEqual(m2,1e-12));
  CPPUNIT_ASSERT(std::string(m2->getDescription())==text1);
  m2->setDescription("ggg");
  CPPUNIT_ASSERT(!m->isEqual(m2,1e-12));
  m2->decrRef();
  //
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f->setTimeUnit(text1);
  CPPUNIT_ASSERT(std::string(f->getTimeUnit())==text1);
  MEDCouplingFieldDouble *f2=f->deepCpy();
  CPPUNIT_ASSERT(std::string(f2->getTimeUnit())==text1);
  f2->decrRef();
  //
  f->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testMultiFields1()
{
  MEDCouplingMultiFields *mfs=buildMultiFields_1();
  std::vector<MEDCouplingMesh *> ms=mfs->getMeshes();
  std::vector<int> refs;
  std::vector<MEDCouplingMesh *> dms=mfs->getDifferentMeshes(refs);
  std::vector<DataArrayDouble *> das=mfs->getArrays();
  std::vector< std::vector<int> > refs2;
  std::vector<DataArrayDouble *> das2=mfs->getDifferentArrays(refs2);
  //
  CPPUNIT_ASSERT_EQUAL(5,(int)ms.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)dms.size());
  CPPUNIT_ASSERT_EQUAL(6,(int)das.size());
  CPPUNIT_ASSERT_EQUAL(5,(int)das2.size());
  //
  MEDCouplingMultiFields *mfs2=mfs->deepCpy();
  CPPUNIT_ASSERT(mfs->isEqual(mfs2,1e-12,1e-12));
  mfs2->decrRef();
  //
  mfs->decrRef();
}

void MEDCouplingBasicsTest::testFieldOverTime1()
{
  std::vector<MEDCouplingFieldDouble *> fs=buildMultiFields_2();
  CPPUNIT_ASSERT_THROW(MEDCouplingFieldOverTime::New(fs),INTERP_KERNEL::Exception);
  MEDCouplingFieldDouble *f4bis=fs[4]->buildNewTimeReprFromThis(ONE_TIME,false);
  fs[4]->decrRef();
  fs[4]=f4bis;
  CPPUNIT_ASSERT_THROW(MEDCouplingFieldOverTime::New(fs),INTERP_KERNEL::Exception);
  f4bis->setTime(2.7,20,21);
  MEDCouplingFieldOverTime *fot=MEDCouplingFieldOverTime::New(fs);
  MEDCouplingDefinitionTime dt=fot->getDefinitionTimeZone();
  std::vector<double> hs=dt.getHotSpotsTime();
  CPPUNIT_ASSERT_EQUAL(6,(int)hs.size());
  const double expected1[]={0.2,0.7,1.2,1.35,1.7,2.7};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],hs[i],1e-12);
  int meshId,arrId,arrIdInField,fieldId;
  dt.getIdsOnTimeRight(0.2,meshId,arrId,arrIdInField,fieldId);
  CPPUNIT_ASSERT_EQUAL(0,meshId);
  CPPUNIT_ASSERT_EQUAL(0,arrId);
  CPPUNIT_ASSERT_EQUAL(0,arrIdInField);
  CPPUNIT_ASSERT_EQUAL(0,fieldId);
  //
  dt.getIdsOnTimeRight(0.7,meshId,arrId,arrIdInField,fieldId);
  CPPUNIT_ASSERT_EQUAL(0,meshId);
  CPPUNIT_ASSERT_EQUAL(1,arrId);
  CPPUNIT_ASSERT_EQUAL(0,arrIdInField);
  CPPUNIT_ASSERT_EQUAL(1,fieldId);
  //
  dt.getIdsOnTimeLeft(1.2,meshId,arrId,arrIdInField,fieldId);//**** WARNING left here
  CPPUNIT_ASSERT_EQUAL(0,meshId);
  CPPUNIT_ASSERT_EQUAL(2,arrId);
  CPPUNIT_ASSERT_EQUAL(1,arrIdInField);
  CPPUNIT_ASSERT_EQUAL(1,fieldId);
  //
  dt.getIdsOnTimeRight(1.2,meshId,arrId,arrIdInField,fieldId);//**** WARNING right again here
  CPPUNIT_ASSERT_EQUAL(1,meshId);
  CPPUNIT_ASSERT_EQUAL(3,arrId);
  CPPUNIT_ASSERT_EQUAL(0,arrIdInField);
  CPPUNIT_ASSERT_EQUAL(2,fieldId);
  //
  dt.getIdsOnTimeRight(1.35,meshId,arrId,arrIdInField,fieldId);
  CPPUNIT_ASSERT_EQUAL(1,meshId);
  CPPUNIT_ASSERT_EQUAL(3,arrId);
  CPPUNIT_ASSERT_EQUAL(0,arrIdInField);
  CPPUNIT_ASSERT_EQUAL(2,fieldId);
  //
  dt.getIdsOnTimeRight(1.7,meshId,arrId,arrIdInField,fieldId);
  CPPUNIT_ASSERT_EQUAL(0,meshId);
  CPPUNIT_ASSERT_EQUAL(3,arrId);
  CPPUNIT_ASSERT_EQUAL(0,arrIdInField);
  CPPUNIT_ASSERT_EQUAL(3,fieldId);
  //
  dt.getIdsOnTimeRight(2.7,meshId,arrId,arrIdInField,fieldId);
  CPPUNIT_ASSERT_EQUAL(1,meshId);
  CPPUNIT_ASSERT_EQUAL(4,arrId);
  CPPUNIT_ASSERT_EQUAL(0,arrIdInField);
  CPPUNIT_ASSERT_EQUAL(4,fieldId);
  //
  MEDCouplingDefinitionTime dt2;
  CPPUNIT_ASSERT(!dt2.isEqual(dt));
  dt2.assign(dt);
  dt2.assign(dt);//to check memory management
  CPPUNIT_ASSERT(dt2.isEqual(dt));
  //
  MEDCouplingDefinitionTime dt3;
  std::vector<int> tmp1;
  std::vector<double> tmp2;
  CPPUNIT_ASSERT(!dt2.isEqual(dt3));
  dt2.getTinySerializationInformation(tmp1,tmp2);
  dt3.unserialize(tmp1,tmp2);
  CPPUNIT_ASSERT(dt2.isEqual(dt3));
  //
  for(std::vector<MEDCouplingFieldDouble *>::iterator it=fs.begin();it!=fs.end();it++)
    (*it)->decrRef();
  fot->decrRef();
}

void MEDCouplingBasicsTest::testDAICheckAndPreparePermutation1()
{
  const int vals1[]={9,10,0,6,4,11,3,7};
  const int expect1[]={5,6,0,3,2,7,1,4};
  const int vals2[]={9,10,0,6,10,11,3,7};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(8,1);
  std::copy(vals1,vals1+8,da->getPointer());
  DataArrayInt *da2=da->checkAndPreparePermutation();
  CPPUNIT_ASSERT_EQUAL(8,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  for(int i=0;i<8;i++)
    CPPUNIT_ASSERT_EQUAL(expect1[i],da2->getIJ(i,0));
  da2->decrRef();
  da->decrRef();
  //
  da=DataArrayInt::New();
  da->alloc(8,1);
  da->iota(0);
  da2=da->checkAndPreparePermutation();
  CPPUNIT_ASSERT_EQUAL(8,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  CPPUNIT_ASSERT(da2->isIdentity());
  da2->decrRef();
  da->decrRef();
  //
  da=DataArrayInt::New();
  da->alloc(8,1);
  std::copy(vals2,vals2+8,da->getPointer());
  CPPUNIT_ASSERT_THROW(da->checkAndPreparePermutation(),INTERP_KERNEL::Exception);
  da->decrRef();
}

void MEDCouplingBasicsTest::testDAIChangeSurjectiveFormat1()
{
  const int vals1[8]={0,3,2,3,2,2,1,2};
  const int expected1[5]={0,1,2,6,8};
  const int expected2[8]={0,  6,  2,4,5,7,  1,3};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(8,1);
  std::copy(vals1,vals1+8,da->getPointer());
  //
  DataArrayInt *da2,*da2I;
  da->changeSurjectiveFormat(4,da2,da2I);
  CPPUNIT_ASSERT_EQUAL(5,da2I->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(8,da2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+5,da2I->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+8,da2->getConstPointer()));
  da2->decrRef();
  da2I->decrRef();
  //
  CPPUNIT_ASSERT_THROW(da->changeSurjectiveFormat(3,da2,da2I),INTERP_KERNEL::Exception);
  //
  da->decrRef();
}

void MEDCouplingBasicsTest::testUMeshGetCellIdsLyingOnNodes1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  const int nodeIds1[5]={1,2,3,4,6};
  const int nodeIds2[2]={6,7};
  DataArrayInt *da=m->getCellIdsLyingOnNodes(nodeIds1,nodeIds1+5,true);
  CPPUNIT_ASSERT_EQUAL(1,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(1,da->getIJ(0,0));
  da->decrRef();
  da=m->getCellIdsLyingOnNodes(nodeIds2,nodeIds2+2,false);
  CPPUNIT_ASSERT_EQUAL(2,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,da->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(4,da->getIJ(1,0));
  da->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest::testUMeshFindCellsIdsOnBoundary1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  DataArrayInt *da5=m->findCellsIdsOnBoundary();
  CPPUNIT_ASSERT_EQUAL(5,da5->getNumberOfTuples());
  CPPUNIT_ASSERT(da5->isIdentity());
  //
  da5->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testMeshSetTime1()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  MEDCouplingUMesh *m2=build3DSurfTargetMesh_1();
  //
  CPPUNIT_ASSERT(m1->isEqual(m2,1e-12));
  m1->setTime(3.14,6,7);
  int tmp1,tmp2;
  double tmp3=m1->getTime(tmp1,tmp2);
  CPPUNIT_ASSERT_EQUAL(6,tmp1);
  CPPUNIT_ASSERT_EQUAL(7,tmp2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.14,tmp3,1e-12);
  CPPUNIT_ASSERT(!m1->isEqual(m2,1e-12));
  m2->setTime(3.14,6,7);
  CPPUNIT_ASSERT(m1->isEqual(m2,1e-12));
  m2->setTime(3.14,6,8);
  CPPUNIT_ASSERT(!m1->isEqual(m2,1e-12));
  m2->setTime(3.14,7,7);
  CPPUNIT_ASSERT(!m1->isEqual(m2,1e-12));
  m2->setTime(3.15,6,7);
  CPPUNIT_ASSERT(!m1->isEqual(m2,1e-12));
  //
  m1->setTime(10.34,55,12);
  MEDCouplingUMesh *m3=(MEDCouplingUMesh *)m1->deepCpy();
  CPPUNIT_ASSERT(m1->isEqual(m3,1e-12));
  tmp3=m3->getTime(tmp1,tmp2);
  CPPUNIT_ASSERT_EQUAL(55,tmp1);
  CPPUNIT_ASSERT_EQUAL(12,tmp2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(10.34,tmp3,1e-12);
  //
  m3->decrRef();
  m1->decrRef();
  m2->decrRef();
  // testing CMesh
  const double coo1[4]={0.,1.,2.,3.5};
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(4,1);
  std::copy(coo1,coo1+4,a->getPointer());
  MEDCouplingCMesh *b=MEDCouplingCMesh::New();
  b->setCoordsAt(0,a);
  a->decrRef();
  //
  b->setTime(5.67,8,100);
  tmp3=b->getTime(tmp1,tmp2);
  CPPUNIT_ASSERT_EQUAL(8,tmp1);
  CPPUNIT_ASSERT_EQUAL(100,tmp2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.67,tmp3,1e-12);
  MEDCouplingCMesh *c=(MEDCouplingCMesh *)b->deepCpy();
  CPPUNIT_ASSERT(c->isEqual(b,1e-12));
  tmp3=c->getTime(tmp1,tmp2);
  CPPUNIT_ASSERT_EQUAL(8,tmp1);
  CPPUNIT_ASSERT_EQUAL(100,tmp2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.67,tmp3,1e-12);
  c->decrRef();
  b->decrRef();
}

void MEDCouplingBasicsTest::testApplyFuncTwo1()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setMesh(m1);
  //
  const double vals[15]={1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.};
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(5,3);
  std::copy(vals,vals+15,da->getPointer());
  f1->setArray(da);
  //
  CPPUNIT_ASSERT_THROW(da->applyFunc2(1,"y+z"),INTERP_KERNEL::Exception);
  da->setInfoOnComponent(0,"x [m]");
  da->setInfoOnComponent(1,"y [mm]");
  da->setInfoOnComponent(2,"z [km]");
  DataArrayDouble *da2=da->applyFunc2(1,"y+z");
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,da2->getNumberOfTuples());
  const double expected1[5]={32.,34.,36.,38.,40.};
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da2->getIJ(0,i),1e-12);
  da2->decrRef();
  da2=da->applyFunc(1,"y+z");
  const double expected2[5]={12.,14.,16.,18.,20.};
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],da2->getIJ(0,i),1e-12);
  da2->decrRef();
  //
  CPPUNIT_ASSERT_EQUAL(3,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  f1->applyFunc2(1,"y+z");
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getArray()->getIJ(0,i),1e-12);
  //
  da->decrRef();
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest::testApplyFuncThree1()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setMesh(m1);
  //
  const double vals[15]={1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.};
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(5,3);
  std::copy(vals,vals+15,da->getPointer());
  f1->setArray(da);
  //
  std::vector<std::string> vs(3);
  vs[0]="x"; vs[1]="Y"; vs[2]="z";
  CPPUNIT_ASSERT_THROW(da->applyFunc3(1,vs,"y+z"),INTERP_KERNEL::Exception);
  vs[1]="y";
  DataArrayDouble *da2=da->applyFunc3(1,vs,"y+z");
  const double expected1[5]={32.,34.,36.,38.,40.};
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da2->getIJ(0,i),1e-12);
  da2->decrRef();
  f1->setArray(da);
  CPPUNIT_ASSERT_EQUAL(3,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  f1->applyFunc3(1,vs,"y+z");
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getArray()->getIJ(0,i),1e-12);
  //
  da->decrRef();
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest::testFillFromAnalyticTwo1()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  CPPUNIT_ASSERT_THROW(m1->fillFromAnalytic2(ON_NODES,1,"y+z"),INTERP_KERNEL::Exception);
  m1->getCoords()->setInfoOnComponent(0,"x [m]");
  m1->getCoords()->setInfoOnComponent(1,"y");
  m1->getCoords()->setInfoOnComponent(2,"z");
  MEDCouplingFieldDouble *f1=m1->fillFromAnalytic2(ON_NODES,1,"y+z");
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  const double expected1[9]={0.2, 0.7, 1.2, 0.7, 1.2, 1.7, 1.2, 1.7, 2.2};
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getArray()->getIJ(0,i),1e-12);
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest::testFillFromAnalyticThree1()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  std::vector<std::string> vs(3);
  vs[0]="x"; vs[1]="Y"; vs[2]="z";
  CPPUNIT_ASSERT_THROW(m1->fillFromAnalytic3(ON_NODES,1,vs,"y+z"),INTERP_KERNEL::Exception);
  vs[1]="y";
  MEDCouplingFieldDouble *f1=m1->fillFromAnalytic3(ON_NODES,1,vs,"y+z");
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  const double expected1[9]={0.2, 0.7, 1.2, 0.7, 1.2, 1.7, 1.2, 1.7, 2.2};
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getArray()->getIJ(0,i),1e-12);
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest::testDAUnitVar1()
{
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(1,3);
  da->setInfoOnComponent(0,"XPS [m]");
  std::string st1,st2;
  st1=da->getVarOnComponent(0);
  CPPUNIT_ASSERT(st1=="XPS");
  st2=da->getUnitOnComponent(0);
  CPPUNIT_ASSERT(st2=="m");
  //
  da->setInfoOnComponent(0,"XPS         [m]");
  st1=da->getVarOnComponent(0);
  CPPUNIT_ASSERT(st1=="XPS");
  st2=da->getUnitOnComponent(0);
  CPPUNIT_ASSERT(st2=="m");
  //
  da->setInfoOnComponent(0,"XPP         [m]");
  st1=da->getVarOnComponent(0);
  CPPUNIT_ASSERT(st1=="XPP");
  st2=da->getUnitOnComponent(0);
  CPPUNIT_ASSERT(st2=="m");
  //
  da->setInfoOnComponent(0,"XPP kdep  kefer   [ m  ]");
  st1=da->getVarOnComponent(0);
  CPPUNIT_ASSERT(st1=="XPP kdep  kefer");
  st2=da->getUnitOnComponent(0);
  CPPUNIT_ASSERT(st2==" m  ");
  //
  da->setInfoOnComponent(0,"     XPP k[  dep  k]efer   [ m^ 2/s^3*kJ  ]");
  st1=da->getVarOnComponent(0);
  CPPUNIT_ASSERT(st1=="     XPP k[  dep  k]efer");
  st2=da->getUnitOnComponent(0);
  CPPUNIT_ASSERT(st2==" m^ 2/s^3*kJ  ");
  //
  da->setInfoOnComponent(0,"     XPP kefer   ");
  st1=da->getVarOnComponent(0);
  CPPUNIT_ASSERT(st1=="     XPP kefer   ");
  st2=da->getUnitOnComponent(0);
  CPPUNIT_ASSERT(st2=="");
  //
  da->setInfoOnComponent(0,"temperature( bof)");
  st1=da->getVarOnComponent(0);
  CPPUNIT_ASSERT(st1=="temperature( bof)");
  st2=da->getUnitOnComponent(0);
  CPPUNIT_ASSERT(st2=="");
  //
  da->setInfoOnComponent(0,"kkk [m]");
  da->setInfoOnComponent(1,"ppp   [m^2/kJ]");
  da->setInfoOnComponent(2,"abcde   [MW/s]");
  //
  std::vector<std::string> vs;
  vs=da->getVarsOnComponent();
  CPPUNIT_ASSERT_EQUAL(3,(int)vs.size());
  CPPUNIT_ASSERT(vs[0]=="kkk");
  CPPUNIT_ASSERT(vs[1]=="ppp");
  CPPUNIT_ASSERT(vs[2]=="abcde");
  vs=da->getUnitsOnComponent();
  CPPUNIT_ASSERT_EQUAL(3,(int)vs.size());
  CPPUNIT_ASSERT(vs[0]=="m");
  CPPUNIT_ASSERT(vs[1]=="m^2/kJ");
  CPPUNIT_ASSERT(vs[2]=="MW/s");
  //
  da->decrRef();
}

