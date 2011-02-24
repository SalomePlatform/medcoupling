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
  m1->setTimeUnit("ms");
  CPPUNIT_ASSERT(std::string(m1->getTimeUnit())=="ms");
  m1->setTimeUnit("us");
  CPPUNIT_ASSERT(std::string(m1->getTimeUnit())=="us");
  CPPUNIT_ASSERT(!m1->isEqual(m2,1e-12));
  m2->setTimeUnit("us");
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

void MEDCouplingBasicsTest::testGaussCoordinates1()
{
  //Testing 1D cell types
  MEDCouplingUMesh *m1=build1DMultiTypes_1();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_GAUSS_PT,ONE_TIME);
  f->setMesh(m1);
  std::vector<double> wg1(1); wg1[0]=0.3;
  std::vector<double> gsCoo1(1); gsCoo1[0]=0.2;
  std::vector<double> refCoo1(2); refCoo1[0]=-1.0; refCoo1[1]=1.0;
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_SEG2,refCoo1,gsCoo1,wg1);
  std::vector<double> wg2(wg1);
  std::vector<double> gsCoo2(1); gsCoo2[0]=0.2;
  std::vector<double> refCoo2(3); refCoo2[0]=-1.0; refCoo2[1]=1.0; refCoo2[2]=0.0;
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_SEG3,refCoo2,gsCoo2,wg2);
  //
  DataArrayDouble *resToTest=f->getLocalizationOfDiscr();
  CPPUNIT_ASSERT_EQUAL(3,resToTest->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(2,resToTest->getNumberOfTuples());
  const double expected1[6]={0.6,0.6,0.6, 0.6,0.6,0.6};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],resToTest->getIJ(0,i),1e-14);
  resToTest->decrRef();
  //
  m1->decrRef();
  f->decrRef();
  //Testing 2D cell types
  MEDCouplingUMesh *m2=build2DMultiTypes_1();
  f=MEDCouplingFieldDouble::New(ON_GAUSS_PT,ONE_TIME);
  f->setMesh(m2);
  std::vector<double> wg3(2); wg3[0]=0.3; wg3[1]=0.3;
  const double tria3CooGauss[4]={ 0.1, 0.8, 0.2, 0.7 };
  std::vector<double> gsCoo3(tria3CooGauss,tria3CooGauss+4);
  const double tria3CooRef[6]={ 0.0, 0.0, 1.0 , 0.0, 0.0, 1.0 };
  std::vector<double> refCoo3(tria3CooRef,tria3CooRef+6);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_TRI3,refCoo3,gsCoo3,wg3);
  std::vector<double> wg4(3); wg4[0]=0.3; wg4[1]=0.3; wg4[2]=0.3;
  const double tria6CooGauss[6]={ 0.3, 0.2, 0.2, 0.1, 0.2, 0.4 };
  std::vector<double> gsCoo4(tria6CooGauss,tria6CooGauss+6);
  const double tria6CooRef[12]={0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5};
  std::vector<double> refCoo4(tria6CooRef,tria6CooRef+12);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_TRI6,refCoo4,gsCoo4,wg4);
  std::vector<double> wg5(4); wg5[0]=0.3; wg5[1]=0.3; wg5[2]=0.3; wg5[3]=0.3;
  const double quad4CooGauss[8]={ 0.3, 0.2, 0.2, 0.1, 0.2, 0.4, 0.15, 0.27 };
  std::vector<double> gsCoo5(quad4CooGauss,quad4CooGauss+8);
  const double quad4CooRef[8]={-1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0};
  std::vector<double> refCoo5(quad4CooRef,quad4CooRef+8);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_QUAD4,refCoo5,gsCoo5,wg5);
  std::vector<double> wg6(4); wg6[0]=0.3; wg6[1]=0.3; wg6[2]=0.3; wg6[3]=0.3;
  const double quad8CooGauss[8]={ 0.34, 0.16, 0.21, 0.3, 0.23, 0.4, 0.14, 0.37 };
  std::vector<double> gsCoo6(quad8CooGauss,quad8CooGauss+8);
  const double quad8CooRef[16]={ -1.0, -1.0, 1.0, -1.0, 1.0,  1.0, -1.0,  1.0, 0.0, -1.0, 1.0,  0.0, 0.0,  1.0, -1.0,  0.0};
  std::vector<double> refCoo6(quad8CooRef,quad8CooRef+16);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_QUAD8,refCoo6,gsCoo6,wg6);
  //
  resToTest=f->getLocalizationOfDiscr();
  CPPUNIT_ASSERT_EQUAL(3,resToTest->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(13,resToTest->getNumberOfTuples());//2+3+4+4 gauss points for resp TRI3,TRI6,QUAD4,QUAD8
  const double expected2[39]={5.1,1.55,0.0, 4.7,1.65,0.0, //TRI3
                              2.32,1.52,0.0, 1.6,1.32,0.0, 3.52,1.26,0.0,//TRI6
                              2.6,1.6,0.0, 2.4,1.8,0.0, 2.4,1.2,0.0, 2.3,1.46,0.0,//QUAD4
                              2.32,2.68,0.0, 2.6,2.42,0.0, 2.8,2.46,0.0, 2.74,2.28,0.0 };//QUAD8
  for(int i=0;i<39;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],resToTest->getIJ(0,i),1e-14);
  resToTest->decrRef();
  //
  m2->decrRef();
  f->decrRef();
  //Testing 3D cell types
  MEDCouplingUMesh *m3=build3DMultiTypes_1();
  f=MEDCouplingFieldDouble::New(ON_GAUSS_PT,ONE_TIME);
  f->setMesh(m3);
  //
  std::vector<double> wg7(1); wg7[0]=0.3;
  const double tetra4CooGauss[3]={0.34, 0.16, 0.21};
  std::vector<double> gsCoo7(tetra4CooGauss,tetra4CooGauss+3);
  const double tetra4CooRef[12]={0.0,1.0,0.0, 0.0,0.0,1.0, 0.0,0.0,0.0, 1.0,0.0,0.0};
  std::vector<double> refCoo7(tetra4CooRef,tetra4CooRef+12);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_TETRA4,refCoo7,gsCoo7,wg7);
  std::vector<double> wg8(1); wg8[0]=0.3;
  const double tetra10CooGauss[3]={0.2, 0.3, 0.1};
  std::vector<double> gsCoo8(tetra10CooGauss,tetra10CooGauss+3);
  const double tetra10CooRef[30]={0.0,1.0,0.0, 0.0,0.0,0.0, 0.0,0.0,1.0, 1.0,0.0,0.0, 0.0,0.5,0.0, 0.0,0.0,0.5, 0.0,0.5,0.5, 0.5,0.5,0.0, 0.5,0.0,0.0, 0.5,0.0,0.5};
  std::vector<double> refCoo8(tetra10CooRef,tetra10CooRef+30);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_TETRA10,refCoo8,gsCoo8,wg8);
  std::vector<double> wg9(1); wg9[0]=0.3;
  const double pyra5CooGauss[3]={0.2, 0.3, 0.1};
  std::vector<double> gsCoo9(pyra5CooGauss,pyra5CooGauss+3);
  const double pyra5CooRef[15]={1.0,0.0,0.0, 0.0,1.0,0.0, -1.0,0.0,0.0, 0.0,-1.0,0.0, 0.0,0.0,1.0};
  std::vector<double> refCoo9(pyra5CooRef,pyra5CooRef+15);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_PYRA5,refCoo9,gsCoo9,wg9);
  std::vector<double> wg10(1); wg10[0]=0.3;
  const double pyra13CooGauss[3]={0.1, 0.2, 0.7};
  std::vector<double> gsCoo10(pyra13CooGauss,pyra13CooGauss+3);
  const double pyra13CooRef[39]={1.0,0.0,0.0, 0.0,1.0,0.0,-1.0,0.0,0.0,0.0,-1.0,0.0,0.0,0.0,1.0,0.5,0.5,0.0,-0.5,0.5,0.0,-0.5,-0.5,0.0,0.5,-0.5,0.0,0.5,0.0,0.5,0.0,0.5,0.5,-0.5,0.0,0.5,0.0,-0.5,0.5};
  std::vector<double> refCoo10(pyra13CooRef,pyra13CooRef+39);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_PYRA13,refCoo10,gsCoo10,wg10);
  std::vector<double> wg11(1); wg11[0]=0.3;
  const double penta6CooGauss[3]={0.2, 0.3, 0.1};
  std::vector<double> gsCoo11(penta6CooGauss,penta6CooGauss+3);
  const double penta6CooRef[18]={-1.0,1.0,0.0,-1.0,-0.0,1.0,-1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0};
  std::vector<double> refCoo11(penta6CooRef,penta6CooRef+18);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_PENTA6,refCoo11,gsCoo11,wg11);
  std::vector<double> wg12(1); wg12[0]=0.3;
  const double penta15CooGauss[3]={0.2, 0.3,0.15};
  std::vector<double> gsCoo12(penta15CooGauss,penta15CooGauss+3);
  const double penta15CooRef[45]={-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,0.0,0.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,0.0,0.0,-1.0,0.5,0.5,-1.0,0.0,0.5,-1.0,0.5,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.5,0.5,1.0,0.0, 0.5,1.0,0.5,0.0};
  std::vector<double> refCoo12(penta15CooRef,penta15CooRef+45);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_PENTA15,refCoo12,gsCoo12,wg12);
  std::vector<double> wg13(1); wg13[0]=0.3;
  const double hexa8CooGauss[3]={0.2,0.3,0.15};
  std::vector<double> gsCoo13(hexa8CooGauss,hexa8CooGauss+3);
  const double hexa8CooRef[24]={-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0};
  std::vector<double> refCoo13(hexa8CooRef,hexa8CooRef+24);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_HEXA8,refCoo13,gsCoo13,wg13);
  std::vector<double> wg14(1); wg14[0]=0.3;
  const double hexa20CooGauss[3]={0.11,0.3,0.55};
  std::vector<double> gsCoo14(hexa20CooGauss,hexa20CooGauss+3);
  const double hexa20CooRef[60]={-1.0,-1.0,-1.0,1.0,-1.0,-1.0,1.0,1.0,-1.0,-1.0,1.0,-1.0,-1.0,-1.0,1.0,1.0,-1.0,1.0,1.0,1.0,1.0,-1.0,1.0,1.0,0.0,-1.0,-1.0,1.0,0.0,-1.0,0.0,1.0,-1.0,-1.0,0.0,-1.0,-1.0,-1.0,0.0,1.0,-1.0,0.0,1.0,1.0,0.0,-1.0,1.0,0.0,0.0,-1.0,1.0,1.0,0.0,1.0,0.0,1.0,1.0,-1.0,0.0,1.0};
  std::vector<double> refCoo14(hexa20CooRef,hexa20CooRef+60);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_HEXA20,refCoo14,gsCoo14,wg14);
  //
  resToTest=f->getLocalizationOfDiscr();
  CPPUNIT_ASSERT_EQUAL(3,resToTest->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(8,resToTest->getNumberOfTuples());//2+3+4+4 gauss points for resp TRI3,TRI6,QUAD4,QUAD8
  const double expected3[24]={1.312,3.15,1.02, 0.56,3.3,0.6, 2.18,1.1,0.2, 1.18,1.54,0.98, 1.56,0.3,3.6, 1.613,0.801,4.374, 2.6,2.4,2.3, 2.31232,2.3933985,1.553255};
  for(int i=0;i<24;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],resToTest->getIJ(0,i),1e-14);
  resToTest->decrRef();
  //
  m3->decrRef();
  f->decrRef();
}

/*!
 * Not activated test ! To be implemented ! WARNING Hexa8 should be replaced by Quad8...
 */
void MEDCouplingBasicsTest::testP2Localization1()
{
  MEDCouplingUMesh *m=buildHexa8Mesh_1();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  DataArrayDouble *da=DataArrayDouble::New();
  const double vals1[27]={1.0,3.0,4.0,1.0,3.0,4.0,3.0,2.0,5.0,1.0,3.0,4.0,1.0,3.0,4.0,3.0,2.0,5.0,1.0,3.0,4.0,1.0,3.0,4.0,3.0,2.0,5.0};
  da->alloc(27,1);
  std::copy(vals1,vals1+27,da->getPointer());
  f->setMesh(m);
  f->setArray(da);
  da->decrRef();
  //
  const double point1[3]={0.25,0.75,0.25};
  //const double points1[6]={0.25,0.75,0.25,1.0,1.0,1.0};
  double res1[3];
  f->getValueOn(point1,res1);
  //
  f->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testGetValueOn2()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f->setMesh(m);
  DataArrayDouble *arr=DataArrayDouble::New();
  int nbOfCells=m->getNumberOfCells();
  arr->alloc(nbOfCells,3);
  f->setArray(arr);
  arr->decrRef();
  const double values1[15]={7.,107.,10007.,8.,108.,10008.,9.,109.,10009.,10.,110.,10010.,11.,111.,10011.};
  std::copy(values1,values1+15,arr->getPointer());
  const double loc[10]={-0.05,-0.05, 0.55,-0.25, 0.55,0.15, -0.05,0.45, 0.45,0.45};
  f->checkCoherency();
  DataArrayDouble *locs=f->getValueOnMulti(loc,5);
  CPPUNIT_ASSERT_EQUAL(5,locs->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,locs->getNumberOfComponents());
  for(int j=0;j<15;j++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(values1[j],locs->getIJ(0,j),1e-12);
  locs->decrRef();
  f->decrRef();
  // Testing ON_NODES
  f=MEDCouplingFieldDouble::New(ON_NODES,NO_TIME);
  f->setMesh(m);
  arr=DataArrayDouble::New();
  int nbOfNodes=m->getNumberOfNodes();
  arr->alloc(nbOfNodes,3);
  f->setArray(arr);
  arr->decrRef();
  const double values2[27]={7.,107.,10007.,8.,108.,10008.,9.,109.,10009.,10.,110.,10010.,11.,111.,10011.,12.,112.,10012.,13.,113.,10013.,14.,114.,10014.,15.,115.,10015.};
  std::copy(values2,values2+27,arr->getPointer());
  const double loc2[8]={0.5432,-0.2432, 0.5478,0.1528, 0.5432,-0.2432, 0.5432,-0.2432};
  const double expected2[12]={9.0272, 109.0272, 10009.0272, 11.4124,111.4124,10011.4124, 9.0272, 109.0272, 10009.0272, 9.0272, 109.0272, 10009.0272};
  f->checkCoherency();
  locs=f->getValueOnMulti(loc2,4);
  CPPUNIT_ASSERT_EQUAL(4,locs->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,locs->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],locs->getIJ(0,i),1e-12);
  f->decrRef();
  locs->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest::testDAIGetIdsNotEqual1()
{
  DataArrayInt *d=DataArrayInt::New();
  const int vals1[10]={2,3,5,6,8,5,5,6,1,-5};
  d->alloc(10,1);
  std::copy(vals1,vals1+10,d->getPointer());
  DataArrayInt *d2=d->getIdsNotEqual(5);
  CPPUNIT_ASSERT_EQUAL(7,d2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d2->getNumberOfComponents());
  const int expected1[7]={0,1,3,4,7,8,9};
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],d2->getIJ(0,i));
  d->rearrange(2);
  CPPUNIT_ASSERT_THROW(d->getIdsNotEqual(5),INTERP_KERNEL::Exception);
  const int vals2[3]={-4,5,6};
  std::vector<int> vals3(vals2,vals2+3);
  d->rearrange(1);
  DataArrayInt *d3=d->getIdsNotEqualList(vals3);
  CPPUNIT_ASSERT_EQUAL(5,d3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d3->getNumberOfComponents());
  const int expected2[5]={0,1,4,8,9};
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d3->getIJ(0,i));
  d3->decrRef();
  d->decrRef();
  d2->decrRef();
}

void MEDCouplingBasicsTest::testDAIComputeOffsets1()
{
  DataArrayInt *d=DataArrayInt::New();
  const int vals1[6]={3,5,1,2,0,8};
  const int expected1[6]={0,3,8,9,11,11};
  d->alloc(6,1);
  std::copy(vals1,vals1+6,d->getPointer());
  d->computeOffsets();
  CPPUNIT_ASSERT_EQUAL(6,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d->getNumberOfComponents());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],d->getIJ(0,i));
  d->decrRef();
}
