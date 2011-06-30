// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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
 * Not activated test ! To be implemented !
 */
void MEDCouplingBasicsTest::testQ1Localization1()
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

void MEDCouplingBasicsTest::testP2Localization1()
{
  MEDCouplingUMesh *m=MEDCouplingUMesh::New("testP2",2);
  const double coords[12]={0.,2.,3.5,0.,4.5,1.5,1.2,0.32,3.4,1.,2.1,2.4};
  const int conn[6]={0,1,2,3,4,5};
  DataArrayDouble *coo=DataArrayDouble::New();
  coo->alloc(6,2);
  std::copy(coords,coords+12,coo->getPointer());
  m->setCoords(coo);
  coo->decrRef();
  m->allocateCells(1);
  m->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,conn);
  m->finishInsertingCells();
  //
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  f->setMesh(m);
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(6,3);
  const double vals1[18]={1.2,2.3,3.4, 2.2,3.3,4.4, 3.2,4.3,5.4, 4.2,5.3,6.4, 5.2,6.3,7.4, 6.2,7.3,8.4};
  std::copy(vals1,vals1+18,da->getPointer());
  f->setArray(da);
  da->decrRef();
  //
  const double loc[2]={2.27,1.3};
  DataArrayDouble *locs=f->getValueOnMulti(loc,1);
  const double expected1[3]={6.0921164547752236, 7.1921164547752232, 8.2921164547752255};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],locs->getIJ(0,i),1e-12);
  locs->decrRef();
  //
  m->decrRef();
  f->decrRef();
}

void MEDCouplingBasicsTest::testP2Localization2()
{
  MEDCouplingUMesh *m=MEDCouplingUMesh::New("testP2_2",3);
  const double coords[30]={0.33312787792955395, -0.35155740179580952, -0.03567564825034563, 1.307146326477638, -0.57234557776250305, -0.08608044208272235, 0.5551834466499993, 0.62324964668794192, -0.014638951108536295, 0.37761817224442129, -0.38324019806913578, 0.96283164472856886, 0.79494856035658679, -0.40628057809270046, 0.0021004190225864614, 1.023740446371799, 0.07665912970471335, -0.072889657161871096, 0.54564584619517376, 0.11132872093429744, 0.039647326652013051, 0.27164784387819052, -0.42018012100866675, 0.46563376500745146, 0.89501965094896418, -0.56148455362735061, 0.43337469695473035, 0.49118025152924394, 0.093884938060727313, 0.47216346905220891};
  const int conn[10]={0,1,2,3,4,5,6,7,8,9};
  DataArrayDouble *coo=DataArrayDouble::New();
  coo->alloc(10,3);
  std::copy(coords,coords+30,coo->getPointer());
  m->setCoords(coo);
  coo->decrRef();
  m->allocateCells(1);
  m->insertNextCell(INTERP_KERNEL::NORM_TETRA10,10,conn);
  m->finishInsertingCells();
  //
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  f->setMesh(m);
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(10,1);
  const double vals1[10]={1.1,2.1,3.1,4.1,5.2,6.2,7.2,8.2,9.2,10.2};
  std::copy(vals1,vals1+10,da->getPointer());
  f->setArray(da);
  da->decrRef();
  //
  const double loc[3]={0.64637931739890486, -0.16185896817550552, 0.22678966365273748};
  DataArrayDouble *locs=f->getValueOnMulti(loc,1);
  const double expected1[1]={10.0844021968047};
  for(int i=0;i<1;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],locs->getIJ(0,i),1e-12);
  locs->decrRef();
  //
  m->decrRef();
  f->decrRef();
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

void MEDCouplingBasicsTest::testUMeshHexagonPrism1()
{
  const double coords[36]={
    0.8660254037844386, 0.5, 0.0, 0.0, 1.0, 0.0, -0.8660254037844386, 0.5, 0.0, -0.8660254037844386, -0.5, 0.0, 0.0, -1.0, 0.0, 0.8660254037844386, -0.5, 0.0,
    0.8660254037844386, 0.5, 2.0, 0.0, 1.0, 2.0, -0.8660254037844386, 0.5, 2.0, -0.8660254037844386, -0.5, 2.0, 0.0, -1.0, 2.0, 0.8660254037844386, -0.5, 2.0
  };
  const int conn[12]={1,2,3,4,5,0,7,8,9,10,11,6};
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New("MyFirstHexagonalPrism",3);
  DataArrayDouble *coo=DataArrayDouble::New();
  coo->alloc(12,3);
  std::copy(coords,coords+36,coo->getPointer());
  mesh->setCoords(coo);
  mesh->allocateCells(1);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXGP12,12,conn);
  mesh->finishInsertingCells();
  coo->decrRef();
  //
  mesh->checkCoherency();
  MEDCouplingFieldDouble *vols=mesh->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(1,vols->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,vols->getNumberOfComponents());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-5.196152422706632,vols->getIJ(0,0),1e-12);
  DataArrayDouble *bary=mesh->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(1,bary->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,bary->getNumberOfComponents());
  const double expected1[3]={0.,0.,1.};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],bary->getIJ(0,i),1e-12);
  DataArrayInt *d1=DataArrayInt::New();
  DataArrayInt *d2=DataArrayInt::New();
  DataArrayInt *d3=DataArrayInt::New();
  DataArrayInt *d4=DataArrayInt::New();
  MEDCouplingUMesh *m2=mesh->buildDescendingConnectivity(d1,d2,d3,d4);
  CPPUNIT_ASSERT_EQUAL(8,m2->getNumberOfCells());
  const int expected4[8][6]={{1,2,3,4,5,0},{7,6,11,10,9,8},{1,7,8,2},{2,8,9,3},{3,9,10,4},{4,10,11,5},{5,11,6,0},{0,6,7,1}};
  const INTERP_KERNEL::NormalizedCellType expected2[8]={INTERP_KERNEL::NORM_POLYGON, INTERP_KERNEL::NORM_POLYGON, INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_QUAD4};
  const int expected3[8]={6,6,4,4,4,4,4,4};
  for(int i=0;i<8;i++)
    {
      CPPUNIT_ASSERT(m2->getTypeOfCell(i)==expected2[i]);
      std::vector<int> v;
      m2->getNodeIdsOfCell(i,v);
      CPPUNIT_ASSERT((int)v.size()==expected3[i]);
      CPPUNIT_ASSERT(std::equal(expected4[i],expected4[i]+expected3[i],v.begin()));
    }
  d1->decrRef();
  d2->decrRef();
  d3->decrRef();
  d4->decrRef();
  m2->decrRef();
  //
  mesh->convertAllToPoly();
  CPPUNIT_ASSERT(INTERP_KERNEL::NORM_POLYHED==mesh->getTypeOfCell(0));
  mesh->unPolyze();
  CPPUNIT_ASSERT(INTERP_KERNEL::NORM_HEXGP12==mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(13,mesh->getMeshLength());
  //
  vols->decrRef();
  bary->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testDADCheckIsMonotonic()
{
  DataArrayDouble *da=DataArrayDouble::New();
  const double vals[4]={-1.,1.01,2.03,6.};
  da->alloc(2,2);
  std::copy(vals,vals+4,da->getPointer());
  CPPUNIT_ASSERT_THROW(da->isMonotonic(1e-12),INTERP_KERNEL::Exception);
  da->rearrange(1);
  CPPUNIT_ASSERT(da->isMonotonic(1e-12));
  da->checkMonotonic(1e-12);
  da->setIJ(2,0,6.1);
  CPPUNIT_ASSERT(!da->isMonotonic(1e-12));
  CPPUNIT_ASSERT_THROW(da->checkMonotonic(1e-12),INTERP_KERNEL::Exception);
  da->setIJ(2,0,5.99);
  CPPUNIT_ASSERT(da->isMonotonic(1e-12));
  CPPUNIT_ASSERT(!da->isMonotonic(1e-1));
  da->decrRef();
}

void MEDCouplingBasicsTest::testCheckCoherencyDeeper1()
{
  MEDCouplingUMesh *m=build3DSourceMesh_1();
  m->checkCoherency();
  m->checkCoherency1();
  m->getNodalConnectivity()->setIJ(8,0,-1);
  m->checkCoherency();
  CPPUNIT_ASSERT_THROW(m->checkCoherency1(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(8,0,-6);
  m->checkCoherency();
  CPPUNIT_ASSERT_THROW(m->checkCoherency1(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(8,0,9);//9>=NbOfNodes
  m->checkCoherency();
  CPPUNIT_ASSERT_THROW(m->checkCoherency1(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(8,0,8);//OK
  m->checkCoherency();
  m->checkCoherency1();
  const int elts[2]={1,5};
  std::vector<int> eltsV(elts,elts+2);
  m->convertToPolyTypes(eltsV);
  m->checkCoherency();
  m->checkCoherency1();
  m->getNodalConnectivity()->setIJ(2,0,9);//9>=NbOfNodes
  m->checkCoherency();
  CPPUNIT_ASSERT_THROW(m->checkCoherency1(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(2,0,-3);
  m->checkCoherency();
  CPPUNIT_ASSERT_THROW(m->checkCoherency1(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(2,0,-1);
  m->checkCoherency();
  CPPUNIT_ASSERT_THROW(m->checkCoherency1(),INTERP_KERNEL::Exception);//Throw because cell#0 is not a polyhedron
  m->getNodalConnectivity()->setIJ(2,0,4);
  m->checkCoherency();
  m->checkCoherency1();
  m->getNodalConnectivity()->setIJ(7,0,-1);
  m->checkCoherency();
  m->checkCoherency1();//OK because we are in polyhedron connec
  m->getNodalConnectivity()->setIJ(36,0,14);
  m->checkCoherency();
  CPPUNIT_ASSERT_THROW(m->checkCoherency1(),INTERP_KERNEL::Exception);//Throw beacause now cell 5 is a TETRA4 (14) so mimatch of number index and static type.
  m->decrRef();
}

void MEDCouplingBasicsTest::testUnPolyze2()
{
  MEDCouplingUMesh *m=MEDCouplingUMesh::New("jjj",3);
  DataArrayDouble *coo=DataArrayDouble::New();
  coo->alloc(4,3);
  coo->rearrange(1);
  coo->iota(0);
  coo->rearrange(3);
  m->setCoords(coo);
  coo->decrRef();
  m->allocateCells(2);
  const int conn[4]={0,1,2,3};
  m->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,conn);
  m->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,conn);
  m->finishInsertingCells();
  std::vector<const MEDCouplingUMesh *> ms(4,m);
  MEDCouplingUMesh *m2=MEDCouplingUMesh::MergeUMeshesOnSameCoords(ms);
  std::vector<int> temp(1,2);
  m2->convertToPolyTypes(temp);
  m2->unPolyze();
  CPPUNIT_ASSERT(INTERP_KERNEL::NORM_TETRA4==m2->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(40,m2->getMeshLength());
  std::vector<int> temp2;
  m2->getNodeIdsOfCell(2,temp2);
  CPPUNIT_ASSERT(4==(int)temp2.size());
  CPPUNIT_ASSERT(std::equal(conn,conn+4,temp2.begin()));
  m2->checkCoherency1();
  MEDCouplingMesh *m3=m2->deepCpy();
  m2->unPolyze();
  CPPUNIT_ASSERT(m3->isEqual(m2,1e-12));
  m3->decrRef();
  m->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest::testDACpyFrom1()
{
  DataArrayDouble *d=DataArrayDouble::New();
  d->alloc(12,1);
  d->iota(14.);
  d->rearrange(3);
  d->setName("Toto");
  d->setInfoOnComponent(0,"X [m]");
  d->setInfoOnComponent(1,"Y [m]");
  d->setInfoOnComponent(2,"Z [m]");
  //
  DataArrayDouble *d1=DataArrayDouble::New();
  CPPUNIT_ASSERT(!d->isEqual(*d1,1e-12));
  d1->cpyFrom(*d);
  CPPUNIT_ASSERT(d->isEqual(*d1,1e-12));
  d1->cpyFrom(*d);
  CPPUNIT_ASSERT(d->isEqual(*d1,1e-12));
  d1->rearrange(2);
  CPPUNIT_ASSERT(!d->isEqual(*d1,1e-12));
  d1->cpyFrom(*d);
  CPPUNIT_ASSERT(d->isEqual(*d1,1e-12));
  //
  DataArrayInt *d2=d->convertToIntArr();
  DataArrayInt *d4=DataArrayInt::New();
  CPPUNIT_ASSERT(!d2->isEqual(*d4));
  d4->cpyFrom(*d2);
  CPPUNIT_ASSERT(d2->isEqual(*d4));
  d4->cpyFrom(*d2);
  CPPUNIT_ASSERT(d2->isEqual(*d4));
  d4->rearrange(2);
  CPPUNIT_ASSERT(!d2->isEqual(*d4));
  d4->cpyFrom(*d2);
  CPPUNIT_ASSERT(d2->isEqual(*d4));
  //
  d->decrRef();
  d1->decrRef();
  d2->decrRef();
  d4->decrRef();
}

void MEDCouplingBasicsTest::testDAITransformWithIndArr1()
{
  const int tab1[4]={17,18,22,19};
  const int tab2[12]={0,1,1,3,3,0,1,3,2,2,3,0};
  const int expected[12]={17,18,18,19,19,17,18,19,22,22,19,17};
  DataArrayInt *d=DataArrayInt::New();
  d->alloc(4,1);
  std::copy(tab1,tab1+4,d->getPointer());
  DataArrayInt *d1=DataArrayInt::New();
  d1->alloc(12,1);
  std::copy(tab2,tab2+12,d1->getPointer());
  //
  d1->transformWithIndArr(d->getConstPointer(),d->getConstPointer()+d->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(12,d1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d1->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected[i],d1->getIJ(i,0));
  //
  d->decrRef();
  d1->decrRef();
}

void MEDCouplingBasicsTest::testDAIBuildPermArrPerLevel1()
{
  const int arr[12]={2,0,1,1,0,1,2,0,1,1,0,0};
  const int expected1[12]={10,0,5,6,1,7,11,2,8,9,3,4};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(12,1);
  std::copy(arr,arr+12,da->getPointer());
  DataArrayInt *da2=da->buildPermArrPerLevel();
  CPPUNIT_ASSERT_EQUAL(12,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(i,0));
  da->decrRef();
  da2->decrRef();
}

void MEDCouplingBasicsTest::testDAIOperations1()
{
  const int arr1[12]={-1,-2,4,7,3,2,6,6,4,3,0,1};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(4,3);
  std::copy(arr1,arr1+12,da->getPointer());
  DataArrayInt *da1=DataArrayInt::New();
  da1->alloc(12,1);
  da1->iota(2);
  CPPUNIT_ASSERT_THROW(DataArrayInt::Add(da,da1),INTERP_KERNEL::Exception);//not same number of tuples/Components
  da1->rearrange(3);
  DataArrayInt *da2=DataArrayInt::Add(da,da1);
  CPPUNIT_ASSERT_EQUAL(4,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da2->getNumberOfComponents());
  const int expected1[12]={1,1,8,12,9,9,14,15,14,14,12,14};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(0,i));
  da2->decrRef();
  da1->substractEqual(da);
  const int expected2[12]={3,5,0,-2,3,5,2,3,6,8,12,12};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],da1->getIJ(0,i));
  da1->rearrange(1); da1->iota(2); da1->rearrange(3);
  da1->addEqual(da);
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da1->getIJ(0,i));
  da1->rearrange(1); da1->iota(2); da1->rearrange(3);
  da2=DataArrayInt::Multiply(da,da1);
  CPPUNIT_ASSERT_EQUAL(4,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da2->getNumberOfComponents());
  const int expected3[12]={-2,-6,16,35,18,14,48,54,40,33,0,13};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],da2->getIJ(0,i));
  da2->decrRef();
  da->divideEqual(da1);
  CPPUNIT_ASSERT_EQUAL(4,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da->getNumberOfComponents());
  const int expected4[12]={0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],da->getIJ(0,i));
  std::copy(arr1,arr1+12,da->getPointer());
  da1->multiplyEqual(da);
  CPPUNIT_ASSERT_EQUAL(4,da1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da1->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],da1->getIJ(0,i));
  da1->rearrange(1); da1->iota(2); da1->rearrange(3);
  da2=DataArrayInt::Divide(da,da1);
  CPPUNIT_ASSERT_EQUAL(4,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da2->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],da2->getIJ(0,i));
  da2->decrRef();
  da1->applyInv(321);
  CPPUNIT_ASSERT_EQUAL(4,da1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da1->getNumberOfComponents());
  const int expected5[12]={160,107,80,64,53,45,40,35,32,29,26,24};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected5[i],da1->getIJ(0,i));
  da1->applyDivideBy(2);
  CPPUNIT_ASSERT_EQUAL(4,da1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da1->getNumberOfComponents());
  const int expected6[12]={80,53,40,32,26,22,20,17,16,14,13,12};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected6[i],da1->getIJ(0,i));
  const int expected7[12]={3,4,5,4,5,1,6,3,2,0,6,5};
  da1->applyModulus(7);
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected7[i],da1->getIJ(0,i));
  da1->applyLin(1,1);
  const int expected8[12]={3,3,3,3,3,1,3,3,0,0,3,3};
  da1->applyRModulus(3);
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected8[i],da1->getIJ(0,i));
  //
  da1->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest::testEmulateMEDMEMBDC1()
{
  MEDCouplingUMesh *m1=0;
  MEDCouplingUMesh *m=buildPointe_1(m1);
  DataArrayInt *da1=DataArrayInt::New();
  DataArrayInt *da2=DataArrayInt::New();
  DataArrayInt *da3=0;
  DataArrayInt *da4=0;
  DataArrayInt *da5=0;
  DataArrayInt *da0=0;
  MEDCouplingUMesh *m2=m->emulateMEDMEMBDC(m1,da1,da2,da3,da4,da5,da0);
  const int expected0[47]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,36,37,32,33,34,35,38,39,40,41,42,43,44,45,46};
  const int expected1[6]={1,32,29,23,41,36};
  CPPUNIT_ASSERT_EQUAL(47,da0->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da0->getNumberOfComponents());
  for(int i=0;i<47;i++)
    CPPUNIT_ASSERT_EQUAL(expected0[i],da0->getIJ(0,i));
  CPPUNIT_ASSERT_EQUAL(6,da5->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da5->getNumberOfComponents());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da5->getIJ(0,i));
  const int expected2[70]={0,1,2,3,4,0,5,6,7,4,8,9,1,7,10,11,12,13,14,5,15,16,17,8,18,19,20,10,21,22,23,2,13,24,25,21,16,26,27,12,19,28,29,15,22,30,31,18,36,26,28,30,24,37,32,33,34,35,38,36,39,40,41,42,37,38,43,44,45,46};
  CPPUNIT_ASSERT_EQUAL(70,da1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da1->getNumberOfComponents());
  for(int i=0;i<70;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],da1->getIJ(0,i));
  const int expected3[17]={0,4,8,12,16,20,24,28,32,36,40,44,48,53,58,64,70};
  CPPUNIT_ASSERT_EQUAL(17,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  for(int i=0;i<17;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],da2->getIJ(0,i));
  const int expected4[48]={0,2,4,6,7,9,11,12,14,16,17,19,20,22,24,25,27,29,30,32,34,35,37,39,40,42,43,45,46,48,49,51,52,53,54,55,56,58,60,62,63,64,65,66,67,68,69,70};
  //const int expected4[48]={0,2,4,6,7,9,11,12,14,16,17,19,20,22,24,25,27,29,30,32,34,35,37,39,40,42,43,45,46,48,49,51,52,54,56,57,58,59,60,62,63,64,65,66,67,68,69,70};
  CPPUNIT_ASSERT_EQUAL(48,da4->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da4->getNumberOfComponents());
  for(int i=0;i<48;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],da4->getIJ(0,i));
  const int expected5[70]={0,1,0,3,0,7,0,1,2,1,4,1,2,3,2,5,2,3,6,3,4,9,4,8,4,5,10,5,9,5,6,11,6,10,6,7,8,7,11,7,8,12,8,9,12,9,10,12,10,11,12,11,13,13,13,13,12,14,13,15,14,15,14,14,14,14,15,15,15,15};
  CPPUNIT_ASSERT_EQUAL(70,da3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da3->getNumberOfComponents());
  for(int i=0;i<70;i++)
    CPPUNIT_ASSERT_EQUAL(expected5[i],da3->getIJ(0,i));
  //
  da0->decrRef();
  da1->decrRef();
  da2->decrRef();
  da3->decrRef();
  da4->decrRef();
  da5->decrRef();
  //
  m2->decrRef();
  m1->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testGetLevArrPerCellTypes1()
{
  MEDCouplingUMesh *m1=0;
  MEDCouplingUMesh *m=buildPointe_1(m1);
  m1->decrRef();
  DataArrayInt *d0=DataArrayInt::New();
  DataArrayInt *d1=DataArrayInt::New();
  DataArrayInt *d2=DataArrayInt::New();
  DataArrayInt *d3=DataArrayInt::New();
  m1=m->buildDescendingConnectivity(d0,d1,d2,d3);
  d0->decrRef(); d1->decrRef(); d2->decrRef(); d3->decrRef();
  INTERP_KERNEL::NormalizedCellType order[2]={INTERP_KERNEL::NORM_TRI3,INTERP_KERNEL::NORM_QUAD4};
  DataArrayInt *da1=0;
  DataArrayInt *da0=m1->getLevArrPerCellTypes(order,order+2,da1);
  const int expected0[47]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,1,1,1,1,1,1,1,1};
  const int expected1[47]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,36,37,32,33,34,35,38,39,40,41,42,43,44,45,46};
  CPPUNIT_ASSERT_EQUAL(47,da0->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da0->getNumberOfComponents());
  for(int i=0;i<47;i++)
    CPPUNIT_ASSERT_EQUAL(expected0[i],da0->getIJ(0,i));
  CPPUNIT_ASSERT_EQUAL(2,da1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(36,da1->getIJ(0,0));//36 TRI3
  CPPUNIT_ASSERT_EQUAL(11,da1->getIJ(1,0));//11 QUAD4
  //
  DataArrayInt *da2=da0->buildPermArrPerLevel();
  //
  CPPUNIT_ASSERT_EQUAL(47,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  for(int i=0;i<47;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(0,i));
  da2->decrRef();
  da0->decrRef();
  da1->decrRef();
  //
  m->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest::testSortCellsInMEDFileFrmt1()
{
  MEDCouplingUMesh *m1=0;
  MEDCouplingUMesh *m=buildPointe_1(m1);
  MEDCouplingUMesh *m2=(MEDCouplingUMesh *)m->deepCpy();
  m->setCoords(0);
  const int vals[16]={0,1,2,14,3,12,4,5,15,6,7,8,9,10,11,13};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(16,1);
  std::copy(vals,vals+16,da->getPointer());
  DataArrayInt *daa=da->invertArrayN2O2O2N(16);
  m->renumberCells(daa->getConstPointer(),false);
  daa->decrRef();
  DataArrayInt *da2=m->sortCellsInMEDFileFrmt();
  CPPUNIT_ASSERT(m2->isEqual(m2,1e-12));
  CPPUNIT_ASSERT(da->isEqual(*da2));
  m2->decrRef();
  da2->decrRef();
  da->decrRef();
  m1->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testBuildPartAndReduceNodes1()
{
  MEDCouplingMesh *m=build2DTargetMesh_1();
  const int arr[2]={1,0};
  DataArrayInt *da;
  MEDCouplingMesh *m2=m->buildPartAndReduceNodes(arr,arr+2,da);
  CPPUNIT_ASSERT_EQUAL(5,m2->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,m2->getNumberOfCells());
  MEDCouplingFieldDouble *f=m2->getMeasureField(true);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,f->getArray()->getIJ(0,0),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,f->getArray()->getIJ(1,0),1e-12);
  f->decrRef();
  da->decrRef();
  m2->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testDAITransformWithIndArrR1()
{
  const int tab1[6]={2,4,5,3,6,7};
  const int tab2[12]={-1,-1,0,1,2,3,4,5,-1,-1,-1,-1};
  const int expected[6]={0,3,1,2,4,5};
  DataArrayInt *d=DataArrayInt::New();
  d->alloc(6,1);
  std::copy(tab1,tab1+6,d->getPointer());
  DataArrayInt *d1=DataArrayInt::New();
  d1->alloc(12,1);
  std::copy(tab2,tab2+12,d1->getPointer());
  //
  DataArrayInt *d3=d->transformWithIndArrR(d1->getConstPointer(),d1->getConstPointer()+d1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(6,d3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,d3->getNumberOfComponents());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected[i],d3->getIJ(i,0));
  d3->decrRef();
  //
  d->decrRef();
  d1->decrRef();
}

void MEDCouplingBasicsTest::testDAISplitByValueRange1()
{
  const int val1[9]={6,5,0,3,2,7,8,1,4};
  const int val2[3]={0,4,9};
  DataArrayInt *d=DataArrayInt::New();
  d->alloc(9,1);
  std::copy(val1,val1+9,d->getPointer());
  DataArrayInt *e=0,*f=0,*g=0;
  d->splitByValueRange(val2,val2+3,e,f,g);
  CPPUNIT_ASSERT_EQUAL(9,e->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,e->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(2,g->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,g->getNumberOfComponents());
  //
  const int expected1[9]={1,1,0,0,0,1,1,0,1};
  const int expected2[9]={2,1,0,3,2,3,4,1,0};
  for(int i=0;i<9;i++)
    {
      CPPUNIT_ASSERT_EQUAL(expected1[i],e->getIJ(i,0));
      CPPUNIT_ASSERT_EQUAL(expected2[i],f->getIJ(i,0));
    }
  CPPUNIT_ASSERT_EQUAL(0,g->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(1,g->getIJ(1,0));
  //
  e->decrRef();
  f->decrRef();
  g->decrRef();
  //
  d->setIJ(6,0,9);
  CPPUNIT_ASSERT_THROW(d->splitByValueRange(val2,val2+3,e,f,g),INTERP_KERNEL::Exception);
  //
  d->decrRef();
}

void MEDCouplingBasicsTest::testUMeshSplitProfilePerType1()
{
  const int val0[5]={2,0,1,3,4};
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  m->renumberCells(val0,false);
  std::vector<int> code;
  std::vector<DataArrayInt *> globIdsPerType;
  std::vector<DataArrayInt *> pfls;
  //
  const int val1[3]={0,2,3};
  DataArrayInt *d=DataArrayInt::New();
  d->alloc(3,1);
  d->setName("sup");
  std::copy(val1,val1+3,d->getPointer());
  m->splitProfilePerType(d,code,globIdsPerType,pfls);
  CPPUNIT_ASSERT_EQUAL(6,(int)code.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)globIdsPerType.size());
  const int expected1[6]={3,1,0, 4,2,1};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],code[i]);
  CPPUNIT_ASSERT_EQUAL(1,globIdsPerType[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,globIdsPerType[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(2,globIdsPerType[1]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,globIdsPerType[1]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(2,globIdsPerType[1]->getIJ(1,0));
  globIdsPerType[0]->decrRef();
  globIdsPerType[1]->decrRef();
  CPPUNIT_ASSERT_EQUAL(2,(int)pfls.size());
  CPPUNIT_ASSERT(std::string("sup")==pfls[0]->getName());
  CPPUNIT_ASSERT_EQUAL(1,pfls[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,pfls[0]->getIJ(0,0));
  CPPUNIT_ASSERT(std::string("sup")==pfls[1]->getName());
  CPPUNIT_ASSERT_EQUAL(2,pfls[1]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,pfls[1]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(1,pfls[1]->getIJ(1,0));
  pfls[0]->decrRef();
  pfls[1]->decrRef();
  d->decrRef();
  globIdsPerType.clear();
  pfls.clear();
  code.clear();
  //
  const int val2[4]={0,2,3,4};// all quad4 are selected here ! So no profile for Quads
  d=DataArrayInt::New();
  d->alloc(4,1);
  std::copy(val2,val2+4,d->getPointer());
  m->splitProfilePerType(d,code,globIdsPerType,pfls);
  CPPUNIT_ASSERT_EQUAL(6,(int)code.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)globIdsPerType.size());
  const int expected2[6]={3,1,0, 4,3,-1};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],code[i]);
  CPPUNIT_ASSERT_EQUAL(1,globIdsPerType[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,globIdsPerType[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(3,globIdsPerType[1]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,globIdsPerType[1]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(2,globIdsPerType[1]->getIJ(1,0));
  CPPUNIT_ASSERT_EQUAL(3,globIdsPerType[1]->getIJ(2,0));
  globIdsPerType[0]->decrRef();
  globIdsPerType[1]->decrRef();
  CPPUNIT_ASSERT_EQUAL(1,(int)pfls.size());
  CPPUNIT_ASSERT_EQUAL(1,pfls[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,pfls[0]->getIJ(0,0));
  pfls[0]->decrRef();
  d->decrRef();
  globIdsPerType.clear();
  pfls.clear();
  code.clear();
  //
  const int val3[3]={1,0,2};// all tri3 are selected here but not in the same order ! Profile requested for Tri3
  d=DataArrayInt::New();
  d->alloc(3,1);
  std::copy(val3,val3+3,d->getPointer());
  m->splitProfilePerType(d,code,globIdsPerType,pfls);
  CPPUNIT_ASSERT_EQUAL(6,(int)code.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)globIdsPerType.size());
  const int expected3[6]={3,2,0, 4,1,1};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],code[i]);
  CPPUNIT_ASSERT_EQUAL(2,globIdsPerType[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,globIdsPerType[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(1,globIdsPerType[0]->getIJ(1,0));
  CPPUNIT_ASSERT_EQUAL(1,globIdsPerType[1]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,globIdsPerType[1]->getIJ(0,0));
  globIdsPerType[0]->decrRef();
  globIdsPerType[1]->decrRef();
  CPPUNIT_ASSERT_EQUAL(2,(int)pfls.size());
  CPPUNIT_ASSERT_EQUAL(2,pfls[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,pfls[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(0,pfls[0]->getIJ(1,0));
  CPPUNIT_ASSERT_EQUAL(0,pfls[1]->getIJ(0,0));
  pfls[0]->decrRef();
  pfls[1]->decrRef();
  d->decrRef();
  globIdsPerType.clear();
  pfls.clear();
  code.clear();
  //
  const int val4[2]={3,4};// all tri3 are selected here but not in the same order ! Profile requested for Tri3
  d=DataArrayInt::New();
  d->alloc(2,1);
  std::copy(val4,val4+2,d->getPointer());
  m->splitProfilePerType(d,code,globIdsPerType,pfls);
  CPPUNIT_ASSERT_EQUAL(3,(int)code.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)globIdsPerType.size());
  const int expected4[3]={4,2,0};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],code[i]);
  CPPUNIT_ASSERT_EQUAL(2,globIdsPerType[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,globIdsPerType[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(1,globIdsPerType[0]->getIJ(1,0));
  globIdsPerType[0]->decrRef();
  CPPUNIT_ASSERT_EQUAL(1,(int)pfls.size());
  CPPUNIT_ASSERT_EQUAL(2,pfls[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,pfls[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(2,pfls[0]->getIJ(1,0));
  pfls[0]->decrRef();
  d->decrRef();
  globIdsPerType.clear();
  pfls.clear();
  code.clear();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest::testDAIBuildExplicitArrByRanges1()
{
  DataArrayInt *d=DataArrayInt::New();
  d->alloc(3,1);
  const int vals1[3]={0,2,3};
  std::copy(vals1,vals1+3,d->getPointer());
  DataArrayInt *e=DataArrayInt::New();
  e->alloc(6,1);
  const int vals2[6]={0,3,6,10,14,20};
  std::copy(vals2,vals2+6,e->getPointer());
  //
  DataArrayInt *f=d->buildExplicitArrByRanges(e);
  CPPUNIT_ASSERT_EQUAL(11,f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f->getNumberOfComponents());
  const int expected1[11]={0,1,2,6,7,8,9,10,11,12,13};
  for(int i=0;i<11;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],f->getIJ(i,0));
  //
  f->decrRef();
  e->decrRef();
  d->decrRef();
}
