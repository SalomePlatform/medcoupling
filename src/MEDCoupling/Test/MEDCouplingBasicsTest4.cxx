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

#include "MEDCouplingBasicsTest4.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingGaussLocalization.hxx"
#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingFieldOverTime.hxx"

#include <cmath>
#include <functional>
#include <iterator>

using namespace MEDCoupling;

void MEDCouplingBasicsTest4::testDescriptionInMeshTimeUnit1()
{
  static const char text1[]="totoTTEDD";
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  m->setDescription(text1);
  CPPUNIT_ASSERT(std::string(m->getDescription())==text1);
  MEDCouplingUMesh *m2=(MEDCouplingUMesh *)m->deepCopy();
  CPPUNIT_ASSERT(m->isEqual(m2,1e-12));
  CPPUNIT_ASSERT(std::string(m2->getDescription())==text1);
  m2->setDescription("ggg");
  CPPUNIT_ASSERT(!m->isEqual(m2,1e-12));
  m2->decrRef();
  //
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f->setTimeUnit(text1);
  CPPUNIT_ASSERT(std::string(f->getTimeUnit())==text1);
  MEDCouplingFieldDouble *f2=f->deepCopy();
  CPPUNIT_ASSERT(std::string(f2->getTimeUnit())==text1);
  f2->decrRef();
  //
  f->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest4::testMultiFields1()
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
  MEDCouplingMultiFields *mfs2=mfs->deepCopy();
  CPPUNIT_ASSERT(mfs->isEqual(mfs2,1e-12,1e-12));
  mfs2->decrRef();
  //
  mfs->decrRef();
}

void MEDCouplingBasicsTest4::testFieldOverTime1()
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

void MEDCouplingBasicsTest4::testDAICheckAndPreparePermutation1()
{
  const int vals1[]={9,10,0,6,4,11,3,7};
  const int expect1[]={5,6,0,3,2,7,1,4};
  const int vals2[]={9,10,0,6,10,11,3,7};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(8,1);
  std::copy(vals1,vals1+8,da->getPointer());
  DataArrayInt *da2=da->checkAndPreparePermutation();
 CPPUNIT_ASSERT_EQUAL(8,(int)da2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  for(int i=0;i<8;i++)
    CPPUNIT_ASSERT_EQUAL(expect1[i],da2->getIJ(i,0));
  da2->decrRef();
  da->decrRef();
  //
  da=DataArrayInt::New();
  da->alloc(8,1);
  da->iota(0);
  da2=da->checkAndPreparePermutation();
 CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  CPPUNIT_ASSERT(da2->isIota(8));
  da2->decrRef();
  da->decrRef();
  //
  da=DataArrayInt::New();
  da->alloc(8,1);
  std::copy(vals2,vals2+8,da->getPointer());
  CPPUNIT_ASSERT_THROW(da->checkAndPreparePermutation(),INTERP_KERNEL::Exception);
  da->decrRef();
}

void MEDCouplingBasicsTest4::testDAIChangeSurjectiveFormat1()
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
 CPPUNIT_ASSERT_EQUAL(5,(int)da2I->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(8,(int)da2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+5,da2I->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+8,da2->getConstPointer()));
  da2->decrRef();
  da2I->decrRef();
  //
  CPPUNIT_ASSERT_THROW(da->changeSurjectiveFormat(3,da2,da2I),INTERP_KERNEL::Exception);
  //
  da->decrRef();
}

void MEDCouplingBasicsTest4::testUMeshGetCellIdsLyingOnNodes1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  const int nodeIds1[5]={1,2,3,4,6};
  const int nodeIds2[2]={6,7};
  DataArrayInt *da=m->getCellIdsLyingOnNodes(nodeIds1,nodeIds1+5,true);
 CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(1,da->getIJ(0,0));
  da->decrRef();
  da=m->getCellIdsLyingOnNodes(nodeIds2,nodeIds2+2,false);
 CPPUNIT_ASSERT_EQUAL(2,(int)da->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,da->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(4,da->getIJ(1,0));
  da->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest4::testUMeshFindCellIdsOnBoundary1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  DataArrayInt *da5=m->findCellIdsOnBoundary();
  CPPUNIT_ASSERT(da5->isIota(5));
  //
  da5->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest4::testMeshSetTime1()
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
  MEDCouplingUMesh *m3=(MEDCouplingUMesh *)m1->deepCopy();
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
  MEDCouplingCMesh *c=(MEDCouplingCMesh *)b->deepCopy();
  CPPUNIT_ASSERT(c->isEqual(b,1e-12));
  tmp3=c->getTime(tmp1,tmp2);
  CPPUNIT_ASSERT_EQUAL(8,tmp1);
  CPPUNIT_ASSERT_EQUAL(100,tmp2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.67,tmp3,1e-12);
  c->decrRef();
  b->decrRef();
}

void MEDCouplingBasicsTest4::testApplyFuncTwo1()
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
  CPPUNIT_ASSERT_THROW(da->applyFuncCompo(1,"y+z"),INTERP_KERNEL::Exception);
  da->setInfoOnComponent(0,"x [m]");
  da->setInfoOnComponent(1,"y [mm]");
  da->setInfoOnComponent(2,"z [km]");
  
  CPPUNIT_ASSERT_THROW(da->applyFuncCompo(1,"x+y+zz+zzz"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->applyFuncCompo(1,"toto(x+y)"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->applyFuncCompo(1,"x/0"),INTERP_KERNEL::Exception);
  
  DataArrayDouble *da2=da->applyFuncCompo(1,"y+z");
 CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(5,(int)da2->getNumberOfTuples());
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
 CPPUNIT_ASSERT_EQUAL(3,(int)f1->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(5,(int)f1->getNumberOfTuples());
  f1->applyFuncCompo(1,"y+z");
 CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(5,(int)f1->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getArray()->getIJ(0,i),1e-12);
  //
  da->decrRef();
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest4::testApplyFuncThree1()
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
  CPPUNIT_ASSERT_THROW(da->applyFuncNamedCompo(1,vs,"y+z"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->applyFuncNamedCompo(1,vs,"x+Y+z+zz+zzz"),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->applyFuncNamedCompo(1,vs,"x/0."),INTERP_KERNEL::Exception);
  vs[1]="y";
  DataArrayDouble *da2=da->applyFuncNamedCompo(1,vs,"y+z");
  const double expected1[5]={32.,34.,36.,38.,40.};
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da2->getIJ(0,i),1e-12);
  da2->decrRef();
  std::vector<std::string> vs2(4); vs2[0]="x"; vs2[1]="y"; vs2[2]="z"; vs2[3]="a";
  CPPUNIT_ASSERT_THROW(da->applyFuncNamedCompo(1,vs2,"x+a"),INTERP_KERNEL::Exception);
  f1->setArray(da);
 CPPUNIT_ASSERT_EQUAL(3,(int)f1->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(5,(int)f1->getNumberOfTuples());
  f1->applyFuncNamedCompo(1,vs,"y+z");
 CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(5,(int)f1->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getArray()->getIJ(0,i),1e-12);
  //
  da->decrRef();
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest4::testFillFromAnalyticTwo1()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  m1->setTime(3.4,5,6); m1->setTimeUnit("us");
  int a,b;
  CPPUNIT_ASSERT_THROW(m1->fillFromAnalyticCompo(ON_NODES,1,"y+z"),INTERP_KERNEL::Exception);
  m1->getCoords()->setInfoOnComponent(0,"x [m]");
  m1->getCoords()->setInfoOnComponent(1,"y");
  m1->getCoords()->setInfoOnComponent(2,"z");
  MEDCouplingFieldDouble *f1=m1->fillFromAnalyticCompo(ON_NODES,1,"y+z");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.4,f1->getTime(a,b),1.e-14);
  CPPUNIT_ASSERT_EQUAL(5,a); CPPUNIT_ASSERT_EQUAL(6,b);
  CPPUNIT_ASSERT_EQUAL(std::string(f1->getTimeUnit()),std::string("us"));
 CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  const double expected1[9]={0.2, 0.7, 1.2, 0.7, 1.2, 1.7, 1.2, 1.7, 2.2};
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getArray()->getIJ(0,i),1e-12);
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest4::testFillFromAnalyticThree1()
{
  MEDCouplingUMesh *m1=build3DSurfTargetMesh_1();
  m1->setTime(3.4,5,6); m1->setTimeUnit("us");
  int a,b;
  std::vector<std::string> vs(3);
  vs[0]="x"; vs[1]="Y"; vs[2]="z";
  CPPUNIT_ASSERT_THROW(m1->fillFromAnalyticNamedCompo(ON_NODES,1,vs,"y+z"),INTERP_KERNEL::Exception);
  vs[1]="y";
  MEDCouplingFieldDouble *f1=m1->fillFromAnalyticNamedCompo(ON_NODES,1,vs,"y+z");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.4,f1->getTime(a,b),1.e-14);
  CPPUNIT_ASSERT_EQUAL(5,a); CPPUNIT_ASSERT_EQUAL(6,b);
  CPPUNIT_ASSERT_EQUAL(std::string(f1->getTimeUnit()),std::string("us"));
 CPPUNIT_ASSERT_EQUAL(1,(int)f1->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(9,(int)f1->getNumberOfTuples());
  const double expected1[9]={0.2, 0.7, 1.2, 0.7, 1.2, 1.7, 1.2, 1.7, 2.2};
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getArray()->getIJ(0,i),1e-12);
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest4::testDAUnitVar1()
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

void MEDCouplingBasicsTest4::testGaussCoordinates1()
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
 CPPUNIT_ASSERT_EQUAL(3,(int)resToTest->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(2,(int)resToTest->getNumberOfTuples());
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
  CPPUNIT_ASSERT_EQUAL(3,(int)resToTest->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(13,(int)resToTest->getNumberOfTuples());//2+3+4+4 gauss points for resp TRI3,TRI6,QUAD4,QUAD8
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
  CPPUNIT_ASSERT_EQUAL(3,(int)resToTest->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(8,(int)resToTest->getNumberOfTuples());//2+3+4+4 gauss points for resp TRI3,TRI6,QUAD4,QUAD8
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
void MEDCouplingBasicsTest4::testQ1Localization1()
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

void MEDCouplingBasicsTest4::testP2Localization1()
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

void MEDCouplingBasicsTest4::testP2Localization2()
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

void MEDCouplingBasicsTest4::testGetValueOn2()
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
  f->checkConsistencyLight();
  DataArrayDouble *locs=f->getValueOnMulti(loc,5);
 CPPUNIT_ASSERT_EQUAL(5,(int)locs->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)locs->getNumberOfComponents());
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
  f->checkConsistencyLight();
  locs=f->getValueOnMulti(loc2,4);
 CPPUNIT_ASSERT_EQUAL(4,(int)locs->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)locs->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],locs->getIJ(0,i),1e-12);
  f->decrRef();
  locs->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest4::testDAIGetIdsNotEqual1()
{
  DataArrayInt *d=DataArrayInt::New();
  const int vals1[10]={2,3,5,6,8,5,5,6,1,-5};
  d->alloc(10,1);
  std::copy(vals1,vals1+10,d->getPointer());
  DataArrayInt *d2=d->findIdsNotEqual(5);
 CPPUNIT_ASSERT_EQUAL(7,(int)d2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d2->getNumberOfComponents());
  const int expected1[7]={0,1,3,4,7,8,9};
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],d2->getIJ(0,i));
  d->rearrange(2);
  CPPUNIT_ASSERT_THROW(d->findIdsNotEqual(5),INTERP_KERNEL::Exception);
  const int vals2[3]={-4,5,6};
  std::vector<int> vals3(vals2,vals2+3);
  d->rearrange(1);
  DataArrayInt *d3=d->findIdsNotEqualList(&vals3[0],&vals3[0]+vals3.size());
 CPPUNIT_ASSERT_EQUAL(5,(int)d3->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d3->getNumberOfComponents());
  const int expected2[5]={0,1,4,8,9};
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d3->getIJ(0,i));
  d3->decrRef();
  d->decrRef();
  d2->decrRef();
}

void MEDCouplingBasicsTest4::testDAIComputeOffsets1()
{
  DataArrayInt *d=DataArrayInt::New();
  const int vals1[6]={3,5,1,2,0,8};
  const int expected1[6]={0,3,8,9,11,11};
  d->alloc(6,1);
  std::copy(vals1,vals1+6,d->getPointer());
  d->computeOffsets();
 CPPUNIT_ASSERT_EQUAL(6,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],d->getIJ(0,i));
  d->decrRef();
}

void MEDCouplingBasicsTest4::testUMeshHexagonPrism1()
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
  mesh->checkConsistencyLight();
  MEDCouplingFieldDouble *vols=mesh->getMeasureField(false);
 CPPUNIT_ASSERT_EQUAL(1,(int)vols->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)vols->getNumberOfComponents());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-5.196152422706632,vols->getIJ(0,0),1e-12);
  DataArrayDouble *bary=mesh->computeCellCenterOfMass();
 CPPUNIT_ASSERT_EQUAL(1,(int)bary->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)bary->getNumberOfComponents());
  const double expected1[3]={0.,0.,1.};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],bary->getIJ(0,i),1e-12);
  DataArrayInt *d1=DataArrayInt::New();
  DataArrayInt *d2=DataArrayInt::New();
  DataArrayInt *d3=DataArrayInt::New();
  DataArrayInt *d4=DataArrayInt::New();
  MEDCouplingUMesh *m2=mesh->buildDescendingConnectivity(d1,d2,d3,d4);
  CPPUNIT_ASSERT_EQUAL(8,(int)m2->getNumberOfCells());
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
  CPPUNIT_ASSERT_EQUAL(13,mesh->getNodalConnectivityArrayLen());
  //
  vols->decrRef();
  bary->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest4::testDADCheckIsMonotonic()
{
  DataArrayDouble *da=DataArrayDouble::New();
  const double vals[4]={-1.,1.01,2.03,6.};
  da->alloc(2,2);
  std::copy(vals,vals+4,da->getPointer());
  CPPUNIT_ASSERT_THROW(da->isMonotonic(true, 1e-12),INTERP_KERNEL::Exception);
  da->rearrange(1);
  CPPUNIT_ASSERT(da->isMonotonic(true, 1e-12));
  da->checkMonotonic(true, 1e-12);
  da->setIJ(2,0,6.1);
  CPPUNIT_ASSERT(!da->isMonotonic(true, 1e-12));
  CPPUNIT_ASSERT_THROW(da->checkMonotonic(true, 1e-12),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(da->checkMonotonic(false, 1e-12),INTERP_KERNEL::Exception);
  da->setIJ(2,0,5.99);
  CPPUNIT_ASSERT(da->isMonotonic(true, 1e-12));
  CPPUNIT_ASSERT(!da->isMonotonic(true, 1e-1));
  da->decrRef();
}

void MEDCouplingBasicsTest4::testCheckCoherencyDeeper1()
{
  MEDCouplingUMesh *m=build3DSourceMesh_1();
  m->checkConsistencyLight();
  m->checkConsistency();
  m->getNodalConnectivity()->setIJ(8,0,-1);
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_THROW(m->checkConsistency(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(8,0,-6);
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_THROW(m->checkConsistency(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(8,0,9);//9>=NbOfNodes
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_THROW(m->checkConsistency(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(8,0,8);//OK
  m->checkConsistencyLight();
  m->checkConsistency();
  const int elts[2]={1,5};
  std::vector<int> eltsV(elts,elts+2);
  m->convertToPolyTypes(&eltsV[0],&eltsV[0]+eltsV.size());
  m->checkConsistencyLight();
  m->checkConsistency();
  m->getNodalConnectivity()->setIJ(2,0,9);//9>=NbOfNodes
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_THROW(m->checkConsistency(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(2,0,-3);
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_THROW(m->checkConsistency(),INTERP_KERNEL::Exception);
  m->getNodalConnectivity()->setIJ(2,0,-1);
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_THROW(m->checkConsistency(),INTERP_KERNEL::Exception);//Throw because cell#0 is not a polyhedron
  m->getNodalConnectivity()->setIJ(2,0,4);
  m->checkConsistencyLight();
  m->checkConsistency();
  m->getNodalConnectivity()->setIJ(7,0,-1);
  m->checkConsistencyLight();
  m->checkConsistency();//OK because we are in polyhedron connec
  m->getNodalConnectivity()->setIJ(36,0,14);
  m->checkConsistencyLight();
  CPPUNIT_ASSERT_THROW(m->checkConsistency(),INTERP_KERNEL::Exception);//Throw beacause now cell 5 is a TETRA4 (14) so mimatch of number index and static type.
  m->decrRef();
}

void MEDCouplingBasicsTest4::testUnPolyze2()
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
  m2->convertToPolyTypes(&temp[0],&temp[0]+temp.size());
  m2->unPolyze();
  CPPUNIT_ASSERT(INTERP_KERNEL::NORM_TETRA4==m2->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(40,m2->getNodalConnectivityArrayLen());
  std::vector<int> temp2;
  m2->getNodeIdsOfCell(2,temp2);
  CPPUNIT_ASSERT(4==(int)temp2.size());
  CPPUNIT_ASSERT(std::equal(conn,conn+4,temp2.begin()));
  m2->checkConsistency();
  MEDCouplingMesh *m3=m2->deepCopy();
  m2->unPolyze();
  CPPUNIT_ASSERT(m3->isEqual(m2,1e-12));
  m3->decrRef();
  m->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest4::testDACpyFrom1()
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
  d1->deepCopyFrom(*d);
  CPPUNIT_ASSERT(d->isEqual(*d1,1e-12));
  d1->deepCopyFrom(*d);
  CPPUNIT_ASSERT(d->isEqual(*d1,1e-12));
  d1->rearrange(2);
  CPPUNIT_ASSERT(!d->isEqual(*d1,1e-12));
  d1->deepCopyFrom(*d);
  CPPUNIT_ASSERT(d->isEqual(*d1,1e-12));
  //
  MCAuto<DataArrayInt> d2=d->convertToIntArr();
  DataArrayInt *d4=DataArrayInt::New();
  CPPUNIT_ASSERT(!d2->isEqual(*d4));
  d4->deepCopyFrom(*d2);
  CPPUNIT_ASSERT(d2->isEqual(*d4));
  d4->deepCopyFrom(*d2);
  CPPUNIT_ASSERT(d2->isEqual(*d4));
  d4->rearrange(2);
  CPPUNIT_ASSERT(!d2->isEqual(*d4));
  d4->deepCopyFrom(*d2);
  CPPUNIT_ASSERT(d2->isEqual(*d4));
  //
  d->decrRef();
  d1->decrRef();
  d4->decrRef();
}

void MEDCouplingBasicsTest4::testDAITransformWithIndArr1()
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
 CPPUNIT_ASSERT_EQUAL(12,(int)d1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d1->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected[i],d1->getIJ(i,0));
  //
  d->decrRef();
  d1->decrRef();
}

void MEDCouplingBasicsTest4::testDAIBuildPermArrPerLevel1()
{
  const int arr[12]={2,0,1,1,0,1,2,0,1,1,0,0};
  const int expected1[12]={10,0,5,6,1,7,11,2,8,9,3,4};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(12,1);
  std::copy(arr,arr+12,da->getPointer());
  DataArrayInt *da2=da->buildPermArrPerLevel();
 CPPUNIT_ASSERT_EQUAL(12,(int)da2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(i,0));
  da->decrRef();
  da2->decrRef();
}

void MEDCouplingBasicsTest4::testDAIOperations1()
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
 CPPUNIT_ASSERT_EQUAL(4,(int)da2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)da2->getNumberOfComponents());
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
 CPPUNIT_ASSERT_EQUAL(4,(int)da2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)da2->getNumberOfComponents());
  const int expected3[12]={-2,-6,16,35,18,14,48,54,40,33,0,13};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],da2->getIJ(0,i));
  da2->decrRef();
  da->divideEqual(da1);
 CPPUNIT_ASSERT_EQUAL(4,(int)da->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)da->getNumberOfComponents());
  const int expected4[12]={0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],da->getIJ(0,i));
  std::copy(arr1,arr1+12,da->getPointer());
  da1->multiplyEqual(da);
 CPPUNIT_ASSERT_EQUAL(4,(int)da1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)da1->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],da1->getIJ(0,i));
  da1->rearrange(1); da1->iota(2); da1->rearrange(3);
  da2=DataArrayInt::Divide(da,da1);
 CPPUNIT_ASSERT_EQUAL(4,(int)da2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)da2->getNumberOfComponents());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],da2->getIJ(0,i));
  da2->decrRef();
  da1->applyInv(321);
 CPPUNIT_ASSERT_EQUAL(4,(int)da1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)da1->getNumberOfComponents());
  const int expected5[12]={160,107,80,64,53,45,40,35,32,29,26,24};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(expected5[i],da1->getIJ(0,i));
  da1->applyDivideBy(2);
 CPPUNIT_ASSERT_EQUAL(4,(int)da1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(3,(int)da1->getNumberOfComponents());
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

void MEDCouplingBasicsTest4::testEmulateMEDMEMBDC1()
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
 CPPUNIT_ASSERT_EQUAL(47,(int)da0->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da0->getNumberOfComponents());
  for(int i=0;i<47;i++)
    CPPUNIT_ASSERT_EQUAL(expected0[i],da0->getIJ(0,i));
 CPPUNIT_ASSERT_EQUAL(6,(int)da5->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da5->getNumberOfComponents());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da5->getIJ(0,i));
  const int expected2[70]={0,1,2,3,4,0,5,6,7,4,8,9,1,7,10,11,12,13,14,5,15,16,17,8,18,19,20,10,21,22,23,2,13,24,25,21,16,26,27,12,19,28,29,15,22,30,31,18,36,26,28,30,24,37,32,33,34,35,38,36,39,40,41,42,37,38,43,44,45,46};
 CPPUNIT_ASSERT_EQUAL(70,(int)da1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da1->getNumberOfComponents());
  for(int i=0;i<70;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],da1->getIJ(0,i));
  const int expected3[17]={0,4,8,12,16,20,24,28,32,36,40,44,48,53,58,64,70};
 CPPUNIT_ASSERT_EQUAL(17,(int)da2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  for(int i=0;i<17;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],da2->getIJ(0,i));
  const int expected4[48]={0,2,4,6,7,9,11,12,14,16,17,19,20,22,24,25,27,29,30,32,34,35,37,39,40,42,43,45,46,48,49,51,52,53,54,55,56,58,60,62,63,64,65,66,67,68,69,70};
  //const int expected4[48]={0,2,4,6,7,9,11,12,14,16,17,19,20,22,24,25,27,29,30,32,34,35,37,39,40,42,43,45,46,48,49,51,52,54,56,57,58,59,60,62,63,64,65,66,67,68,69,70};
 CPPUNIT_ASSERT_EQUAL(48,(int)da4->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da4->getNumberOfComponents());
  for(int i=0;i<48;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],da4->getIJ(0,i));
  const int expected5[70]={0,1,0,3,0,7,0,1,2,1,4,1,2,3,2,5,2,3,6,3,4,9,4,8,4,5,10,5,9,5,6,11,6,10,6,7,8,7,11,7,8,12,8,9,12,9,10,12,10,11,12,11,13,13,13,13,12,14,13,15,14,15,14,14,14,14,15,15,15,15};
 CPPUNIT_ASSERT_EQUAL(70,(int)da3->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da3->getNumberOfComponents());
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

void MEDCouplingBasicsTest4::testGetLevArrPerCellTypes1()
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
 CPPUNIT_ASSERT_EQUAL(47,(int)da0->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da0->getNumberOfComponents());
  for(int i=0;i<47;i++)
    CPPUNIT_ASSERT_EQUAL(expected0[i],da0->getIJ(0,i));
 CPPUNIT_ASSERT_EQUAL(2,(int)da1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(36,da1->getIJ(0,0));//36 TRI3
  CPPUNIT_ASSERT_EQUAL(11,da1->getIJ(1,0));//11 QUAD4
  //
  DataArrayInt *da2=da0->buildPermArrPerLevel();
  //
 CPPUNIT_ASSERT_EQUAL(47,(int)da2->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)da2->getNumberOfComponents());
  for(int i=0;i<47;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(0,i));
  da2->decrRef();
  da0->decrRef();
  da1->decrRef();
  //
  m->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest4::testSortCellsInMEDFileFrmt1()
{
  MEDCouplingUMesh *m1=0;
  MEDCouplingUMesh *m=buildPointe_1(m1);
  MEDCouplingUMesh *m2=(MEDCouplingUMesh *)m->deepCopy();
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

void MEDCouplingBasicsTest4::testBuildPartAndReduceNodes1()
{
  MEDCouplingMesh *m=build2DTargetMesh_1();
  const int arr[2]={1,0};
  DataArrayInt *da;
  MEDCouplingMesh *m2=m->buildPartAndReduceNodes(arr,arr+2,da);
  CPPUNIT_ASSERT_EQUAL(5,m2->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)m2->getNumberOfCells());
  MEDCouplingFieldDouble *f=m2->getMeasureField(true);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.125,f->getArray()->getIJ(0,0),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.25,f->getArray()->getIJ(1,0),1e-12);
  f->decrRef();
  da->decrRef();
  m2->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest4::testDAITransformWithIndArrR1()
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
 CPPUNIT_ASSERT_EQUAL(6,(int)d3->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d3->getNumberOfComponents());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected[i],d3->getIJ(i,0));
  d3->decrRef();
  //
  d->decrRef();
  d1->decrRef();
}

void MEDCouplingBasicsTest4::testDAISplitByValueRange1()
{
  const int val1[9]={6,5,0,3,2,7,8,1,4};
  const int val2[3]={0,4,9};
  DataArrayInt *d=DataArrayInt::New();
  d->alloc(9,1);
  std::copy(val1,val1+9,d->getPointer());
  DataArrayInt *ee=0,*f=0,*g=0;
  d->splitByValueRange(val2,val2+3,ee,f,g);
 CPPUNIT_ASSERT_EQUAL(9,(int)ee->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)ee->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(9,(int)f->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)f->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(2,(int)g->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)g->getNumberOfComponents());
  //
  const int expected1[9]={1,1,0,0,0,1,1,0,1};
  const int expected2[9]={2,1,0,3,2,3,4,1,0};
  for(int i=0;i<9;i++)
    {
      CPPUNIT_ASSERT_EQUAL(expected1[i],ee->getIJ(i,0));
      CPPUNIT_ASSERT_EQUAL(expected2[i],f->getIJ(i,0));
    }
  CPPUNIT_ASSERT_EQUAL(0,g->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(1,g->getIJ(1,0));
  //
  ee->decrRef();
  f->decrRef();
  g->decrRef();
  //
  d->setIJ(6,0,9);
  CPPUNIT_ASSERT_THROW(d->splitByValueRange(val2,val2+3,ee,f,g),INTERP_KERNEL::Exception);
  //
  d->decrRef();
}

void MEDCouplingBasicsTest4::testUMeshSplitProfilePerType1()
{
  const int val0[5]={2,0,1,3,4};
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  m->renumberCells(val0,false);
  std::vector<int> code;
  std::vector<DataArrayInt *> idsInPflPerType;
  std::vector<DataArrayInt *> pfls;
  //
  const int val1[3]={0,2,3};
  DataArrayInt *d=DataArrayInt::New();
  d->alloc(3,1);
  d->setName("sup");
  std::copy(val1,val1+3,d->getPointer());
  m->splitProfilePerType(d,code,idsInPflPerType,pfls);
  CPPUNIT_ASSERT_EQUAL(6,(int)code.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)idsInPflPerType.size());
  const int expected1[6]={3,1,0, 4,2,1};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],code[i]);
 CPPUNIT_ASSERT_EQUAL(1,(int)idsInPflPerType[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,idsInPflPerType[0]->getIJ(0,0));
 CPPUNIT_ASSERT_EQUAL(2,(int)idsInPflPerType[1]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,idsInPflPerType[1]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(2,idsInPflPerType[1]->getIJ(1,0));
  idsInPflPerType[0]->decrRef();
  idsInPflPerType[1]->decrRef();
  CPPUNIT_ASSERT_EQUAL(2,(int)pfls.size());
  CPPUNIT_ASSERT(std::string("sup")==pfls[0]->getName());
 CPPUNIT_ASSERT_EQUAL(1,(int)pfls[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,pfls[0]->getIJ(0,0));
  CPPUNIT_ASSERT(std::string("sup")==pfls[1]->getName());
 CPPUNIT_ASSERT_EQUAL(2,(int)pfls[1]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,pfls[1]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(1,pfls[1]->getIJ(1,0));
  pfls[0]->decrRef();
  pfls[1]->decrRef();
  d->decrRef();
  idsInPflPerType.clear();
  pfls.clear();
  code.clear();
  //
  const int val2[4]={0,2,3,4};// all quad4 are selected here ! So no profile for Quads
  d=DataArrayInt::New();
  d->alloc(4,1);
  std::copy(val2,val2+4,d->getPointer());
  m->splitProfilePerType(d,code,idsInPflPerType,pfls);
  CPPUNIT_ASSERT_EQUAL(6,(int)code.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)idsInPflPerType.size());
  const int expected2[6]={3,1,0, 4,3,-1};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],code[i]);
 CPPUNIT_ASSERT_EQUAL(1,(int)idsInPflPerType[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,idsInPflPerType[0]->getIJ(0,0));
 CPPUNIT_ASSERT_EQUAL(3,(int)idsInPflPerType[1]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,idsInPflPerType[1]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(2,idsInPflPerType[1]->getIJ(1,0));
  CPPUNIT_ASSERT_EQUAL(3,idsInPflPerType[1]->getIJ(2,0));
  idsInPflPerType[0]->decrRef();
  idsInPflPerType[1]->decrRef();
  CPPUNIT_ASSERT_EQUAL(1,(int)pfls.size());
 CPPUNIT_ASSERT_EQUAL(1,(int)pfls[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,pfls[0]->getIJ(0,0));
  pfls[0]->decrRef();
  d->decrRef();
  idsInPflPerType.clear();
  pfls.clear();
  code.clear();
  //
  const int val3[3]={1,0,2};// all tri3 are selected here but not in the same order ! Profile requested for Tri3
  d=DataArrayInt::New();
  d->alloc(3,1);
  std::copy(val3,val3+3,d->getPointer());
  m->splitProfilePerType(d,code,idsInPflPerType,pfls);
  CPPUNIT_ASSERT_EQUAL(6,(int)code.size());
  CPPUNIT_ASSERT_EQUAL(2,(int)idsInPflPerType.size());
  const int expected3[6]={3,2,0, 4,1,1};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected3[i],code[i]);
 CPPUNIT_ASSERT_EQUAL(2,(int)idsInPflPerType[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,idsInPflPerType[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(1,idsInPflPerType[0]->getIJ(1,0));
 CPPUNIT_ASSERT_EQUAL(1,(int)idsInPflPerType[1]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,idsInPflPerType[1]->getIJ(0,0));
  idsInPflPerType[0]->decrRef();
  idsInPflPerType[1]->decrRef();
  CPPUNIT_ASSERT_EQUAL(2,(int)pfls.size());
 CPPUNIT_ASSERT_EQUAL(2,(int)pfls[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,pfls[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(0,pfls[0]->getIJ(1,0));
  CPPUNIT_ASSERT_EQUAL(0,pfls[1]->getIJ(0,0));
  pfls[0]->decrRef();
  pfls[1]->decrRef();
  d->decrRef();
  idsInPflPerType.clear();
  pfls.clear();
  code.clear();
  //
  const int val4[2]={3,4};// all tri3 are selected here but not in the same order ! Profile requested for Tri3
  d=DataArrayInt::New();
  d->alloc(2,1);
  std::copy(val4,val4+2,d->getPointer());
  m->splitProfilePerType(d,code,idsInPflPerType,pfls);
  CPPUNIT_ASSERT_EQUAL(3,(int)code.size());
  CPPUNIT_ASSERT_EQUAL(1,(int)idsInPflPerType.size());
  const int expected4[3]={4,2,0};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_EQUAL(expected4[i],code[i]);
 CPPUNIT_ASSERT_EQUAL(2,(int)idsInPflPerType[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(0,idsInPflPerType[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(1,idsInPflPerType[0]->getIJ(1,0));
  idsInPflPerType[0]->decrRef();
  CPPUNIT_ASSERT_EQUAL(1,(int)pfls.size());
 CPPUNIT_ASSERT_EQUAL(2,(int)pfls[0]->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,pfls[0]->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(2,pfls[0]->getIJ(1,0));
  pfls[0]->decrRef();
  d->decrRef();
  idsInPflPerType.clear();
  pfls.clear();
  code.clear();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest4::testDAIBuildExplicitArrByRanges1()
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
 CPPUNIT_ASSERT_EQUAL(11,(int)f->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)f->getNumberOfComponents());
  const int expected1[11]={0,1,2,6,7,8,9,10,11,12,13};
  for(int i=0;i<11;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],f->getIJ(i,0));
  //
  f->decrRef();
  e->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest4::testDAIComputeOffsets2()
{
  DataArrayInt *d=DataArrayInt::New();
  const int vals1[6]={3,5,1,2,0,8};
  const int expected1[7]={0,3,8,9,11,11,19};
  d->alloc(6,1);
  std::copy(vals1,vals1+6,d->getPointer());
  d->computeOffsetsFull();
 CPPUNIT_ASSERT_EQUAL(7,(int)d->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(1,(int)d->getNumberOfComponents());
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],d->getIJ(0,i));
  d->decrRef();
}

void MEDCouplingBasicsTest4::testMergeField3()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  m->getCoords()->setInfoOnComponent(0,"x [m]");
  m->getCoords()->setInfoOnComponent(1,"z [km]");
  m->setName("m");
  m->setDescription("desc");
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setName("f1");
  f1->setMesh(m);
  DataArrayDouble *arr=DataArrayDouble::New();
  arr->alloc(5,2);
  arr->setInfoOnComponent(0,"X [m]");
  arr->setInfoOnComponent(1,"YY [mm]");
  arr->fillWithValue(2.);
  f1->setArray(arr);
  arr->decrRef();
  m->decrRef();
  //
  std::vector<const MEDCouplingFieldDouble *> tmp(1);
  tmp[0]=f1;
  MEDCouplingFieldDouble *f2=MEDCouplingFieldDouble::MergeFields(tmp);
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  //
  f1->decrRef();
  f2->decrRef();
}

void MEDCouplingBasicsTest4::testGetDistributionOfTypes1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  const int tab1[5]={2,0,1,3,4};
  CPPUNIT_ASSERT_THROW(m->getDistributionOfTypes(),INTERP_KERNEL::Exception);
  m->renumberCells(tab1,false);
  std::vector<int> code=m->getDistributionOfTypes();
  CPPUNIT_ASSERT_EQUAL(6,(int)code.size());
  CPPUNIT_ASSERT_EQUAL(3,code[0]);
  CPPUNIT_ASSERT_EQUAL(2,code[1]);
  CPPUNIT_ASSERT_EQUAL(-1,code[2]);
  CPPUNIT_ASSERT_EQUAL(4,code[3]);
  CPPUNIT_ASSERT_EQUAL(3,code[4]);
  CPPUNIT_ASSERT_EQUAL(-1,code[5]);
  m->decrRef();
}

void MEDCouplingBasicsTest4::testNorm2_1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f->setMesh(m);
  m->decrRef();
  //
  DataArrayDouble *d=DataArrayDouble::New();
  const double tab[10]={1.2,1.3,2.2,2.3,3.2,3.3,4.2,4.3,5.2,5.3};
  d->alloc(5,2);
  std::copy(tab,tab+10,d->getPointer());
  f->setArray(d);
  d->decrRef();
  f->checkConsistencyLight();
  //
  CPPUNIT_ASSERT_DOUBLES_EQUAL(11.209371079592289,f->norm2(),1e-14);
  //
  f->decrRef();
}

void MEDCouplingBasicsTest4::testNormMax1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f->setMesh(m);
  m->decrRef();
  //
  DataArrayDouble *d=DataArrayDouble::New();
  const double tab[10]={2.3,-1.2,6.3,-7.8,2.9,7.7,2.1,0.,3.6,-7.6};
  d->alloc(5,2);
  std::copy(tab,tab+10,d->getPointer());
  f->setArray(d);
  d->decrRef();
  f->checkConsistencyLight();
  //
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.8,f->normMax(),1e-14);
  //
  f->decrRef();
}

void MEDCouplingBasicsTest4::testFindAndCorrectBadOriented3DExtrudedCells1()
{
  const double coords[38*3]={0.0011180339887498999, -0.0011755705045849499, 0.0, -0.0012331070204200001, -0.0011755705045849499, 0.0, -0.00067557050458494599, -0.00145964954842536, 0.0, -0.00050000000000000001, -0.00086602540378443902, 0.0, 0.00140211303259031, -0.00061803398874989504, 0.0, 0.00086602540378443902, -0.00050000000000000001, 0.0, 0.001, 0.0, 0.0, 0.00034561537182258202, 0.000269164072574575, 0.0, 0.0, 0.001, 0.0, -0.00050000000000000001, 0.00086602540378443902, 0.0, -0.000269164072574575, 0.00034561537182258202, 0.0, -0.001, 0.0, 0.0, -0.00086602540378443902, -0.00050000000000000001, 0.0, -0.00034561537182258202, -0.000269164072574575, 0.0, 0.0, -0.001, 0.0, 0.00050000000000000001, -0.00086602540378443902, 0.0, 0.000269164072574575, -0.00034561537182258202, 0.0, 0.0015, -6.01853107621011e-36, 0.0, 0.00056049747291484397, -0.00145964954842536, 0.0, 0.0011180339887498999, -0.0011755705045849499, 0.00050000000000000001, -0.0012331070204200001, -0.0011755705045849499, 0.00050000000000000001, -0.00067557050458494599, -0.00145964954842536, 0.00050000000000000001, -0.00050000000000000001, -0.00086602540378443902, 0.00050000000000000001, 0.00140211303259031, -0.00061803398874989504, 0.00050000000000000001, 0.00086602540378443902, -0.00050000000000000001, 0.00050000000000000001, 0.001, 0.0, 0.00050000000000000001, 0.00034561537182258202, 0.000269164072574575, 0.00050000000000000001, 0.0, 0.001, 0.00050000000000000001, -0.00050000000000000001, 0.00086602540378443902, 0.00050000000000000001, -0.000269164072574575, 0.00034561537182258202, 0.00050000000000000001, -0.001, 0.0, 0.00050000000000000001, -0.00086602540378443902, -0.00050000000000000001, 0.00050000000000000001, -0.00034561537182258202, -0.000269164072574575, 0.00050000000000000001, 0.0, -0.001, 0.00050000000000000001, 0.00050000000000000001, -0.00086602540378443902, 0.00050000000000000001, 0.000269164072574575, -0.00034561537182258202, 0.00050000000000000001, 0.0015, -6.01853107621011e-36, 0.00050000000000000001, 0.00056049747291484397, -0.00145964954842536, 0.00050000000000000001};
  const int conn[56]={2, 1, 3, 21, 20, 22, 4, 0, 5, 23, 19, 24, 8, 9, 10, 27, 28, 29, 11, 12, 13, 30, 31, 32, 0, 18, 15, 5, 19, 37, 34, 24, 6, 17, 4, 5, 25, 36, 23, 24, 3, 14, 16, 13, 22, 33, 35, 32, 13, 16, 7, 10, 32, 35, 26, 29};
  const int connExp[64]={16, 2, 1, 3, 21, 20, 22, 16, 4, 0, 5, 23, 19, 24, 16, 8, 10, 9, 27, 29, 28, 16, 11, 13, 12, 30, 32, 31, 18, 0, 18, 15, 5, 19, 37, 34, 24,18, 6, 17, 4, 5, 25, 36, 23, 24, 18, 3, 13, 16, 14, 22, 32, 35, 33, 18, 13, 10, 7, 16, 32, 29, 26, 35};
  const int invalidCells[4]={2,3,6,7};
  MEDCouplingUMesh *m=MEDCouplingUMesh::New("Example",3);
  DataArrayDouble *coo=DataArrayDouble::New();
  coo->alloc(38,3);
  std::copy(coords,coords+114,coo->getPointer());
  m->setCoords(coo);
  coo->decrRef();
  m->allocateCells(8);
  m->insertNextCell(INTERP_KERNEL::NORM_PENTA6,6,conn);
  m->insertNextCell(INTERP_KERNEL::NORM_PENTA6,6,conn+6);
  m->insertNextCell(INTERP_KERNEL::NORM_PENTA6,6,conn+12);
  m->insertNextCell(INTERP_KERNEL::NORM_PENTA6,6,conn+18);
  m->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+24);
  m->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+32);
  m->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+40);
  m->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+48);
  m->finishInsertingCells();
  //
  DataArrayInt *v=m->findAndCorrectBadOriented3DExtrudedCells();
 CPPUNIT_ASSERT_EQUAL(4,(int)v->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(v->begin(),v->end(),invalidCells));
  CPPUNIT_ASSERT(std::equal(connExp,connExp+64,m->getNodalConnectivity()->getConstPointer()));
  v->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest4::testConvertExtrudedPolyhedra1()
{
  const int conn[72]={1,2,3,4, 5,6,7,8,9,10,11,12, 13,14,15,16, 17,18,19,20,21,22, 23,24,25,26,27,28, 29,30,31,32,33,34,35,36,37,38, 39,40,41,42,43,44,45,46, 47,48,49,50,51,52,53,54,55,56,57,58, 59,60,61,62,63,64,65,66,67,68,69,70,71,72};
  MEDCouplingUMesh *m=MEDCouplingUMesh::New("Example",3);
  DataArrayDouble *coo=DataArrayDouble::New();
  coo->alloc(73,3);
  coo->rearrange(1); coo->iota(0); coo->rearrange(3);
  m->setCoords(coo);
  coo->decrRef();
  m->allocateCells(9);
  m->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,conn);
  m->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+4);
  m->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,conn+12);
  m->insertNextCell(INTERP_KERNEL::NORM_POLYHED,6,conn+16);
  m->insertNextCell(INTERP_KERNEL::NORM_PENTA6,6,conn+22);
  m->insertNextCell(INTERP_KERNEL::NORM_POLYHED,10,conn+28);
  m->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+38);
  m->insertNextCell(INTERP_KERNEL::NORM_HEXGP12,12,conn+46);
  m->insertNextCell(INTERP_KERNEL::NORM_POLYHED,14,conn+58);
  m->finishInsertingCells();
  //
  m->convertExtrudedPolyhedra();
  DataArrayInt *da=m->getNodalConnectivity();
  DataArrayInt *dai=m->getNodalConnectivityIndex();
  CPPUNIT_ASSERT_EQUAL((std::size_t)10,dai->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)159,da->getNbOfElems());
  //
  const int expected1[159]={14,1,2,3,4,
                            18,5,6,7,8,9,10,11,12,
                            14,13,14,15,16,
                            31,17,18,19,-1,20,22,21,-1,17,20,21,18,-1,18,21,22,19,-1,19,22,20,17,
                            16,23,24,25,26,27,28,
                            31,29,30,31,32,33,-1,34,38,37,36,35,-1,29,34,35,30,-1,30,35,36,31,-1,31,36,37,32,-1,32,37,38,33,-1,33,38,34,29,
                            18,39,40,41,42,43,44,45,46,
                            22,47,48,49,50,51,52,53,54,55,56,57,58,
                            31,59,60,61,62,63,64,65,-1,66,72,71,70,69,68,67,-1,59,66,67,60,-1,60,67,68,61,-1,61,68,69,62,-1,62,69,70,63,-1,63,70,71,64,-1,64,71,72,65,-1,65,72,66,59};
  const int expected2[10]={0,5,14,19,42,49,86,95,108,159};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+159,da->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+10,dai->getConstPointer()));
  m->checkConsistency();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest4::testNonRegressionCopyTinyStrings()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m->getMeasureField(true);
  f1->getArray()->setInfoOnComponent(0,"P [N/m^2]");
  DataArrayDouble *bary=m->computeCellCenterOfMass();
  MEDCouplingFieldDouble *f2=f1->buildNewTimeReprFromThis(NO_TIME,false);
  f2->setArray(bary);
  CPPUNIT_ASSERT_THROW(f1->copyTinyAttrFrom(f2),INTERP_KERNEL::Exception);
  m->decrRef();
  f1->decrRef();
  bary->decrRef();
  f2->decrRef();
}

void MEDCouplingBasicsTest4::testDaDSetPartOfValuesAdv1()
{
  const double tab1[18]={3.,4.,5., 13.,14.,15., 23.,24.,25., 33.,34.,35., 43.,44.,45., 53.,54.,55.};
  const double tab2[9]={6.,7.,8., 16.,17.,18., 26.,27.,28.};
  const int tab3[6]={4,1, 2,2, 3,0};
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(6,3);
  std::copy(tab1,tab1+18,a->getPointer());
  DataArrayDouble *b=DataArrayDouble::New();
  b->alloc(3,3);
  std::copy(tab2,tab2+9,b->getPointer());
  DataArrayInt *c=DataArrayInt::New();
  c->alloc(3,2);
  std::copy(tab3,tab3+6,c->getPointer());
  //
  a->setPartOfValuesAdv(b,c);
  const double expected1[18]={3.,4.,5., 13.,14.,15., 26.,27.,28., 6.,7.,8., 16.,17.,18., 53.,54.,55.};
  std::equal(expected1,expected1+18,a->getConstPointer());
  //
  a->decrRef();
  b->decrRef();
  c->decrRef();
}

void MEDCouplingBasicsTest4::testUMeshBuildSetInstanceFromThis1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  MEDCouplingUMesh *m2=m->buildSetInstanceFromThis(3);
  CPPUNIT_ASSERT_EQUAL(m->getNodalConnectivity(),m2->getNodalConnectivity());
  CPPUNIT_ASSERT_EQUAL(m->getNodalConnectivityIndex(),m2->getNodalConnectivityIndex());
  CPPUNIT_ASSERT_EQUAL(m->getCoords(),m2->getCoords());
  m2->decrRef();
  m->decrRef();
  //
  m=MEDCouplingUMesh::New("toto",2);
  m2=m->buildSetInstanceFromThis(3);
  CPPUNIT_ASSERT_EQUAL(0,m2->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(0,(int)m2->getNumberOfCells());
  m->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest4::testUMeshMergeMeshesCVW1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  MEDCouplingUMesh *m2=MEDCouplingUMesh::New("toto",2);
  MEDCouplingUMesh *m3=MEDCouplingUMesh::MergeUMeshes(m,m2);
  m3->setName(m->getName().c_str());
  CPPUNIT_ASSERT(m->isEqual(m3,1e-12));
  m3->decrRef();
  m->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest4::testChangeUnderlyingMeshWithCMesh1()
{
  MEDCouplingCMesh* mesh=MEDCouplingCMesh::New();
  DataArrayDouble* coordsX=DataArrayDouble::New();
  double arrX[4] = { -1., 1., 2., 4. };
  coordsX->useArray(arrX,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsY=DataArrayDouble::New();
  double arrY[4] = { -2., 2., 4., 8. };
  coordsY->useArray(arrY,false, CPP_DEALLOC,4,1);
  DataArrayDouble* coordsZ=DataArrayDouble::New();
  double arrZ[4] = { -3., 3., 6., 12. };
  coordsZ->useArray(arrZ,false, CPP_DEALLOC,4,1);
  mesh->setCoords(coordsX,coordsY,coordsZ);
  coordsX->decrRef();
  coordsY->decrRef();
  coordsZ->decrRef();
  MEDCouplingMesh *mesh2=mesh->deepCopy();
  //
  static const int ids1[9]={0,1,2,10,11,12,20,21,22};
  for(const int *myId=ids1;myId!=ids1+9;myId++)
    {
      MEDCouplingFieldDouble *f=mesh->getMeasureField(true);
      f->changeUnderlyingMesh(mesh2,*myId,1e-12);
      f->decrRef();
    }
  mesh2->setName("uuuu");
  for(const int *myId=ids1+1;myId!=ids1+9;myId++)
    {
      MEDCouplingFieldDouble *f=mesh->getMeasureField(true);
      f->changeUnderlyingMesh(mesh2,*myId,1e-12);
      f->decrRef();
    }
  //
  mesh2->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest4::testDADFindCommonTuples1()
{
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(6,1);
  const double array1[6]={2.3,1.2,1.3,2.3,2.301,0.8};
  std::copy(array1,array1+6,da->getPointer());
  DataArrayInt *c=0,*cI=0;
  // nbOftuples=1
  da->findCommonTuples(1e-2,-1,c,cI);
  const int expected1[3]={0,3,4};
  const int expected2[2]={0,3};
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,c->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)2,cI->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,c->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+2,cI->getConstPointer()));
  c->decrRef();
  cI->decrRef();
  //
  da->findCommonTuples(2e-1,-1,c,cI);
  const int expected3[5]={0,3,4,1,2};
  const int expected4[3]={0,3,5};
  CPPUNIT_ASSERT_EQUAL((std::size_t)5,c->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,cI->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+5,c->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+3,cI->getConstPointer()));
  c->decrRef();
  cI->decrRef();
  // nbOftuples=2
  da->alloc(6,2);
  const double array2[12]={2.3,2.3,1.2,1.2,1.3,1.3,2.3,2.3,2.301,2.301,0.8,0.8};
  std::copy(array2,array2+12,da->getPointer());
  da->findCommonTuples(1e-2,-1,c,cI);
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,c->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)2,cI->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,c->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+2,cI->getConstPointer()));
  c->decrRef();
  cI->decrRef();
  //
  da->findCommonTuples(2e-1,-1,c,cI);
  CPPUNIT_ASSERT_EQUAL((std::size_t)5,c->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,cI->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+5,c->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+3,cI->getConstPointer()));
  c->decrRef();
  cI->decrRef();
  // nbOftuples=3
  da->alloc(6,3);
  const double array3[18]={2.3,2.3,2.3,1.2,1.2,1.2,1.3,1.3,1.3,2.3,2.3,2.3,2.301,2.301,2.301,0.8,0.8,0.8};
  std::copy(array3,array3+18,da->getPointer());
  da->findCommonTuples(1e-2,-1,c,cI);
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,c->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)2,cI->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,c->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+2,cI->getConstPointer()));
  c->decrRef();
  cI->decrRef();
  //
  da->findCommonTuples(2e-1,-1,c,cI);
  CPPUNIT_ASSERT_EQUAL((std::size_t)5,c->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,cI->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+5,c->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+3,cI->getConstPointer()));
  c->decrRef();
  cI->decrRef();
  //
  const double array11[6]={2.3,1.2,1.3,2.4,2.5,0.8};
  da->alloc(6,1);
  std::copy(array11,array11+6,da->getPointer());
  // nbOftuples=1, no common groups
  da->findCommonTuples(1e-2,-1,c,cI);
  CPPUNIT_ASSERT_EQUAL((std::size_t)0,c->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL((std::size_t)1,cI->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(0,cI->getIJ(0,0));
  
  da->alloc(6,5);  //bad NumberOfComponents
  CPPUNIT_ASSERT_THROW(da->findCommonTuples(1e-2,-1,c,cI),INTERP_KERNEL::Exception);
  
  c->decrRef();
  cI->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest4::testDABack1()
{
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(6,1);
  const double array1[6]={2.3,1.2,1.3,2.3,2.301,0.8};
  std::copy(array1,array1+6,da->getPointer());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.8,da->back(),1e-14);
  da->rearrange(2);
  CPPUNIT_ASSERT_THROW(da->back(),INTERP_KERNEL::Exception);
  da->alloc(0,1);
  CPPUNIT_ASSERT_THROW(da->back(),INTERP_KERNEL::Exception);
  da->decrRef();
  //
  DataArrayInt *da2=DataArrayInt::New();
  da2->alloc(4,1);
  const int array2[4]={4,7,8,2};
  std::copy(array2,array2+4,da2->getPointer());
  CPPUNIT_ASSERT_EQUAL(2,da2->back());
  da2->rearrange(2);
  CPPUNIT_ASSERT_THROW(da2->back(),INTERP_KERNEL::Exception);
  da2->alloc(0,1);
  CPPUNIT_ASSERT_THROW(da2->back(),INTERP_KERNEL::Exception);
  da2->decrRef();
}

void MEDCouplingBasicsTest4::testDADGetDifferentValues1()
{
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(6,1);
  const double array1[6]={2.3,1.2,1.3,2.3,2.301,0.8};
  std::copy(array1,array1+6,da->getPointer());
  //
  const double expected1[4]={2.301,1.2,1.3,0.8};
  DataArrayDouble *dv=da->getDifferentValues(1e-2);
  CPPUNIT_ASSERT_EQUAL((std::size_t)4,dv->getNbOfElems());
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],dv->getIJ(i,0),1e-14);
  dv->decrRef();
  //
  dv=da->getDifferentValues(2e-1);
  const double expected2[3]={2.301,1.3,0.8};
  CPPUNIT_ASSERT_EQUAL((std::size_t)3,dv->getNbOfElems());
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],dv->getIJ(i,0),1e-14);
  dv->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest4::testDAIBuildOld2NewArrayFromSurjectiveFormat2()
{
  const int arr[5]={0,3, 5,7,9};
  const int arrI[3]={0,2,5};
  DataArrayInt *a=DataArrayInt::New();
  a->alloc(5,1);
  std::copy(arr,arr+5,a->getPointer());
  DataArrayInt *b=DataArrayInt::New();
  b->alloc(3,1);
  std::copy(arrI,arrI+3,b->getPointer());
  int newNbTuple=-1;
  DataArrayInt *ret=DataArrayInt::ConvertIndexArrayToO2N(10,a->begin(),b->begin(),b->end(),newNbTuple);
  const int expected[10]={0,1,2,0,3,4,5,4,6,4};
  CPPUNIT_ASSERT_EQUAL((std::size_t)10,ret->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(7,newNbTuple);
 CPPUNIT_ASSERT_EQUAL(1,(int)ret->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expected,expected+10,ret->getConstPointer()));
  CPPUNIT_ASSERT_THROW(DataArrayInt::ConvertIndexArrayToO2N(9,a->begin(),b->begin(),b->end(),newNbTuple),INTERP_KERNEL::Exception);
  ret->decrRef();
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest4::testDADIReverse1()
{
  const int arr[6]={0,3,5,7,9,2};
  DataArrayInt *a=DataArrayInt::New();
  a->alloc(6,1);
  std::copy(arr,arr+6,a->getPointer());
  CPPUNIT_ASSERT_EQUAL(2,a->back());
  a->reverse();
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(arr[5-i],a->getIJ(i,0));
  a->alloc(5,1);
  std::copy(arr,arr+5,a->getPointer());
  a->reverse();
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_EQUAL(arr[4-i],a->getIJ(i,0));
  a->decrRef();
  //
  const double arr2[6]={0.,3.,5.,7.,9.,2.};
   DataArrayDouble *b=DataArrayDouble::New();
   b->alloc(6,1);
   std::copy(arr2,arr2+6,b->getPointer());
   b->reverse();
   for(int i=0;i<6;i++)
     CPPUNIT_ASSERT_DOUBLES_EQUAL(arr2[5-i],b->getIJ(i,0),1e-14);
   b->alloc(5,1);
   std::copy(arr,arr+5,b->getPointer());
   CPPUNIT_ASSERT_DOUBLES_EQUAL(9.,b->back(),1e-14);
   b->reverse();
   for(int i=0;i<5;i++)
     CPPUNIT_ASSERT_DOUBLES_EQUAL(arr2[4-i],b->getIJ(i,0),1e-14);
   b->decrRef();
}

void MEDCouplingBasicsTest4::testGetNodeIdsInUse1()
{
  MEDCouplingUMesh *m0=build2DTargetMesh_1();
  const int CellIds[2]={1,2};
  MEDCouplingUMesh *m1=static_cast<MEDCouplingUMesh *>(m0->buildPartOfMySelf(CellIds,CellIds+2,true));
  int newNbOfNodes=-1;
  DataArrayInt *arr=m1->getNodeIdsInUse(newNbOfNodes);
  const int expected[9]={-1,0,1,-1,2,3,-1,-1,-1};
  CPPUNIT_ASSERT_EQUAL(4,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL((std::size_t)9,arr->getNbOfElems());
  CPPUNIT_ASSERT(std::equal(expected,expected+9,arr->getConstPointer()));
  DataArrayInt *arr2=arr->invertArrayO2N2N2O(newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL((std::size_t)4,arr2->getNbOfElems());
  const int expected2[4]={1,2,4,5};
  CPPUNIT_ASSERT(std::equal(expected2,expected2+4,arr2->getConstPointer()));
  arr2->decrRef();
  arr->decrRef();
  m1->decrRef();
  m0->decrRef();
}

void MEDCouplingBasicsTest4::testBuildDescendingConnec2()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  DataArrayInt *desc=DataArrayInt::New();
  DataArrayInt *descIndx=DataArrayInt::New();
  DataArrayInt *revDesc=DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  //
  MEDCouplingUMesh *mesh2=mesh->buildDescendingConnectivity2(desc,descIndx,revDesc,revDescIndx);
  mesh2->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(1,mesh2->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(13,(int)mesh2->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL((std::size_t)14,revDescIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(14,(int)revDescIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)6,descIndx->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(6,(int)descIndx->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)18,desc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,(int)desc->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL((std::size_t)18,revDesc->getNbOfElems()); CPPUNIT_ASSERT_EQUAL(18,(int)revDesc->getNumberOfTuples());
  const int expected1[18]={1,2,3,4,-3,5,6, 7,8,-5,9,10,-2,11, 12,13,-7,-10};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+18,desc->getConstPointer()));
  const int expected2[6]={0,4,7,10,14,18};
  CPPUNIT_ASSERT(std::equal(expected2,expected2+6,descIndx->getConstPointer()));
  const int expected3[14]={0,1,3,5,6,8,9,11,12,13,15,16,17,18};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+14,revDescIndx->getConstPointer()));
  const int expected4[18]={0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+18,revDesc->getConstPointer()));
  DataArrayInt *conn=mesh2->getNodalConnectivity();
  DataArrayInt *connIndex=mesh2->getNodalConnectivityIndex();
  const int expected5[14]={0,3,6,9,12,15,18,21,24,27,30,33,36,39};
  CPPUNIT_ASSERT(std::equal(expected5,expected5+14,connIndex->getConstPointer()));
  const int expected6[39]={1, 0, 3, 1, 3, 4, 1, 4, 1, 1, 1, 0, 1, 4, 2, 1, 2, 1, 1, 4, 5, 1, 5, 2, 1, 6, 7, 1, 7, 4, 1, 3, 6, 1, 7, 8, 1, 8, 5};
  CPPUNIT_ASSERT(std::equal(expected6,expected6+39,conn->getConstPointer()));
  //
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
  mesh2->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest4::testIntersect2DMeshesTmp1()
{
  MEDCouplingCMesh *m1c=MEDCouplingCMesh::New();
  DataArrayDouble *coordX=DataArrayDouble::New();
  const double arrX[4]={-1., 1., 2., 4.};
  coordX->alloc(4,1);
  std::copy(arrX,arrX+4,coordX->getPointer());
  m1c->setCoordsAt(0,coordX);
  DataArrayDouble *coordY=DataArrayDouble::New();
  const double arrY[4]={-2., 2., 4., 8.};
  coordY->alloc(4,1);
  std::copy(arrY,arrY+4,coordY->getPointer());
  m1c->setCoordsAt(1,coordY);
  MEDCouplingUMesh *m1=m1c->buildUnstructured();
  const int subPart1[3]={3,4,5};
  MEDCouplingUMesh *m1bis=static_cast<MEDCouplingUMesh *>(m1->buildPartOfMySelf(subPart1,subPart1+3,false));
  MEDCouplingUMesh *m2tmp=static_cast<MEDCouplingUMesh *>(m1->deepCopy());
  const int subPart2[3]={0,1,2};
  MEDCouplingUMesh *m2=static_cast<MEDCouplingUMesh *>(m2tmp->buildPartOfMySelf(subPart2,subPart2+3,false));
  const double vec[2]={0.5,0.5};
  m2->translate(vec);
  // End of construction of input meshes m1bis and m2 -> start of specific part of the test
  DataArrayInt *d1=0,*d2=0;
  MEDCouplingUMesh *m3=MEDCouplingUMesh::Intersect2DMeshes(m1bis,m2,1e-10,d1,d2);
  const int expected1[8]={0,0,1,1,1,2,2,2};
  const int expected2[8]={0,-1,0,1,-1,1,2,-1};
 CPPUNIT_ASSERT_EQUAL(8,(int)d1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(8,(int)d2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(8,(int)m3->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(22,m3->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,m3->getSpaceDimension());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+8,d1->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+8,d2->getConstPointer()));
  const int expected3[44]={5,17,1,16,12,5,16,0,4,5,17,12,5,18,1,17,13,5,19,2,18,13,5,17,5,6,19,13,5,20,2,19,14,5,21,3,20,14,5,19,6,7,21,14};
  const int expected4[9]={0,5,12,17,22,28,33,38,44};
  const double expected5[44]={-1.0,2.0,1.0,2.0,2.0,2.0,4.0,2.0,-1.0,4.0,1.0,4.0,2.0,4.0,4.0,4.0,-0.5,-1.5,1.5,-1.5,2.5,-1.5,4.5,-1.5,-0.5,2.5,1.5,2.5,2.5,2.5,4.5,2.5,-0.5,2.0,1.0,2.5,1.5,2.0,2.0,2.5,2.5,2.0,4.0,2.5};
  CPPUNIT_ASSERT_EQUAL(44,(int)m3->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,(int)m3->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+44,m3->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+9,m3->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<44;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],m3->getCoords()->getIJ(0,i),1e-12);
  d1->decrRef();
  d2->decrRef();
  m3->decrRef();
  //
  m2->decrRef();
  m2tmp->decrRef();
  m1bis->decrRef();
  m1->decrRef();
  coordX->decrRef();
  coordY->decrRef();
  m1c->decrRef();
}

void MEDCouplingBasicsTest4::testFindNodesOnLine1()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  const double pt[2]={-0.3,-0.3};
  const double pt2[3]={0.,0.,0.};
  const double pt3[3]={-0.3,0.,0.};
  const double vec[2]={0.,1.};
  const double vec2[3]={1.,0.,0.};
  const double vec3[3]={0.,1.,1.};
  const int expected1[3]={0,3,6};
  std::vector<int> res;
  mesh->findNodesOnLine(pt,vec,1e-12,res);
  CPPUNIT_ASSERT_EQUAL(3,(int)res.size());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,res.begin()));
  res.clear();
  //
  mesh->changeSpaceDimension(3);
  mesh->rotate(pt2,vec2,M_PI/4.);
  mesh->findNodesOnLine(pt3,vec3,1e-12,res);
  CPPUNIT_ASSERT_EQUAL(3,(int)res.size());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,res.begin()));
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest4::testIntersect2DMeshesTmp2()
{
  MEDCouplingCMesh *m1c=MEDCouplingCMesh::New();
  DataArrayDouble *coordsX1=DataArrayDouble::New();
  const double arrX1[4]={ 0., 1., 1.5, 2. };
  coordsX1->alloc(4,1);
  std::copy(arrX1,arrX1+4,coordsX1->getPointer());
  m1c->setCoordsAt(0,coordsX1);
  DataArrayDouble *coordsY1=DataArrayDouble::New();
  const double arrY1[3]={ 0., 1.5, 3.};
  coordsY1->alloc(3,1);
  std::copy(arrY1,arrY1+3,coordsY1->getPointer());
  m1c->setCoordsAt(1,coordsY1);
  MEDCouplingUMesh *m1=m1c->buildUnstructured();
  //
  MEDCouplingCMesh *m2c=MEDCouplingCMesh::New();
  DataArrayDouble *coordsX2=DataArrayDouble::New();
  const double arrX2[3]={ 0., 1., 2. };
  coordsX2->alloc(3,1);
  std::copy(arrX2,arrX2+3,coordsX2->getPointer());
  m2c->setCoordsAt(0,coordsX2);
  DataArrayDouble *coordsY2=DataArrayDouble::New();
  coordsY2->alloc(3,1);
  const double arrY2[3]={ 0., 1., 3.};
  std::copy(arrY2,arrY2+3,coordsY2->getPointer());
  m2c->setCoordsAt(1,coordsY2);
  MEDCouplingUMesh *m2=m2c->buildUnstructured();
  //
  DataArrayInt *d1=0,*d2=0;
  MEDCouplingUMesh *m3=MEDCouplingUMesh::Intersect2DMeshes(m1,m2,1e-10,d1,d2);
  const int expected1[9]={0,0,1,1,2,2,3,4,5};
  const int expected2[9]={0,2,1,3,1,3,2,3,3};
 CPPUNIT_ASSERT_EQUAL(9,(int)d1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(9,(int)d2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(9,(int)m3->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(22,m3->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,m3->getSpaceDimension());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+9,d1->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+9,d2->getConstPointer()));
  const int expected3[45]={5,16,13,12,15,5,15,4,5,16,5,21,2,13,16,5,16,5,6,21,5,17,14,2,21,5,21,6,7,17,5,4,18,19,5,5,5,19,10,6,5,6,10,20,7};
  const int expected4[10]={0,5,10,15,20,25,30,35,40,45};
  const double expected5[44]={0.0,0.0,1.0,0.0,1.5,0.0,2.0,0.0,0.0,1.5,1.0,1.5,1.5,1.5,2.0,1.5,0.0,3.0,1.0,3.0,1.5,3.0,2.0,3.0,0.0,0.0,1.0,0.0,2.0,0.0,0.0,1.0,1.0,1.0,2.0,1.0,0.0,3.0,1.0,3.0,2.0,3.0,1.5,1.0};
  CPPUNIT_ASSERT_EQUAL(45,(int)m3->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(10,(int)m3->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+45,m3->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+10,m3->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<44;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],m3->getCoords()->getIJ(0,i),1e-12);
  d1->decrRef();
  d2->decrRef();
  m3->decrRef();
  //
  m1c->decrRef();
  coordsX1->decrRef();
  coordsY1->decrRef();
  m1->decrRef();
  m2c->decrRef();
  coordsX2->decrRef();
  coordsY2->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest4::testBuildPartOfMySelfSafe1()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  const int input1[4]={0,-1,4,2};
  const int input2[4]={0,4,5,4};
  CPPUNIT_ASSERT_THROW(mesh->buildPartOfMySelf(input1,input1+4,true),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_THROW(mesh->buildPartOfMySelf(input2,input2+4,true),INTERP_KERNEL::Exception);
  mesh->decrRef();
}

void MEDCouplingBasicsTest4::testIntersect2DMeshesTmp3()
{
  double m1Coords[50]={0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1,0.,-1.5,0.5,0.,1.25,0.,0.70710678118654757,0.70710678118654757,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.70710678118654757,0.70710678118654757,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.70710678118654757,-0.70710678118654757,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.70710678118654757,-0.70710678118654757,1.0606601717798214,-1.0606601717798214};
  int m1Conn[56]={0,3,1,13,11,9, 3,4,2,1,14,12,10,11, 5,3,0,15,13,17, 6,4,3,5,16,14,15,18, 5,0,7,17,21,19, 6,5,7,8,18,19,22,20, 0,1,7,9,23,21, 1,2,8,7,10,24,22,23};
  MEDCouplingUMesh *m1=MEDCouplingUMesh::New();
  m1->setMeshDimension(2);
  m1->allocateCells(8);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+6);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn+14);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+20);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn+28);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+34);
  m1->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,m1Conn+42);
  m1->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8,m1Conn+48);
  m1->finishInsertingCells();
  DataArrayDouble *myCoords1=DataArrayDouble::New();
  myCoords1->alloc(25,2);
  std::copy(m1Coords,m1Coords+50,myCoords1->getPointer());
  m1->setCoords(myCoords1);
  myCoords1->decrRef();
  //
  double m2Coords[30]={0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1,-1.1,-1.,0.,-1.,1.1,-1,1.7,-1.};
  int m2Conn[32]={0,3,2,1, 1,2,5,4, 7,6,3,0, 8,9,6,7, 7,0,12,11, 8,7,11,10, 0,1,13,12, 1,4,14,13};
  MEDCouplingUMesh *m2=MEDCouplingUMesh::New();
  m2->setMeshDimension(2);
  m2->allocateCells(8);
  for(int i=0;i<8;i++)
    m2->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,m2Conn+4*i);
  m2->finishInsertingCells();
  DataArrayDouble *myCoords2=DataArrayDouble::New();
  myCoords2->alloc(15,2);
  std::copy(m2Coords,m2Coords+30,myCoords2->getPointer());
  m2->setCoords(myCoords2);
  myCoords2->decrRef();
  //
  DataArrayInt *d1=0,*d2=0;
  MEDCouplingUMesh *m3=MEDCouplingUMesh::Intersect2DMeshes(m1,m2,1e-10,d1,d2);
  m3->unPolyze();
  const int expected1[16]={0,1,1,1,2,3,3,3,4,5,5,5,6,7,7,7};
  const int expected2[16]={0,0,1,-1,2,2,3,-1,4,4,5,-1,6,6,7,-1};
 CPPUNIT_ASSERT_EQUAL(16,(int)d1->getNumberOfTuples());
 CPPUNIT_ASSERT_EQUAL(16,(int)d2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(16,(int)m3->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(104,m3->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,m3->getSpaceDimension());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+16,d1->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected2,expected2+16,d2->getConstPointer()));
  const int expected3[136]={6,28,1,25,44,45,46,8,26,1,28,27,47,48,49,50,8,40,2,26,27,51,52,53,54,8,28,4,40,27,55,56,57,58,6,28,25,5,59,60,61,8,28,5,32,31,62,63,64,65,8,32,6,41,31,66,67,68,69,8,41,4,28,31,70,71,72,73,6,25,37,5,74,75,76,8,32,5,37,36,77,78,79,80,8,42,6,32,36,81,82,83,84,8,37,8,42,36,85,86,87,88,6,1,37,25,89,90,91,8,37,1,26,38,92,93,94,95,8,26,2,43,38,96,97,98,99,8,43,8,37,38,100,101,102,103};
  const int expected4[17]={0,7,16,25,34,41,50,59,68,75,84,93,102,109,118,127,136};
  const double expected5[208]={0.,0.,1.,0.,1.5,0.,0.,1.,0.,1.5,-1.,0.,-1.5,0.,0.,-1.,0.,-1.5,0.5,0.,1.25,0.,0.7071067811865476,0.7071067811865476,1.0606601717798214,1.0606601717798214,0.,0.5,0.,1.25,-0.7071067811865476,0.7071067811865476,-1.0606601717798214,1.0606601717798214,-0.5,0.,-1.25,0.,-0.7071067811865476,-0.7071067811865476,-1.0606601717798214,-1.0606601717798214,0.,-0.5,0.,-1.25,0.7071067811865476,-0.7071067811865476,1.0606601717798214,-1.0606601717798214,0.,0.,1.1,0.,1.1,1.,0.,1.,1.7,0.,1.7,1.,-1.1,1.,-1.1,0.,-1.7,0.,-1.7,1.,-1.7,-1.,-1.1,-1.,0.,-1.,1.1,-1.,1.7,-1.,1.118033988749895,1.,-1.118033988749895,1.,-1.118033988749895,-1.,1.118033988749895,-1.,0.7071067811865477,0.7071067811865476,0.5,0.,0.,0.5,1.05,0.,0.7071067811865475,0.7071067811865477,0.55,1.,1.1,0.5,1.4012585384440737,0.535233134659635,1.3,0.,1.1,0.5,1.1090169943749475,1.,0.,1.25,0.6123724356957946,1.369306393762915,1.1090169943749475,1.,0.55,1.,0.,0.5,-0.5,0.,-0.7071067811865477,0.7071067811865476,-0.7071067811865475,0.7071067811865477,-1.05,0.,-1.1,0.5,-0.55,1.,-1.3,0.,-1.4012585384440737,0.5352331346596344,-1.1090169943749475,1.,-1.1,0.5,-0.6123724356957941,1.3693063937629155,0.,1.25,-0.55,1.,-1.1090169943749475,1.,0.,-0.5,-0.7071067811865475,-0.7071067811865477,-0.5,0.,-1.05,0.,-0.7071067811865478,-0.7071067811865475,-0.55,-1.,-1.1,-0.5,-1.4012585384440734,-0.5352331346596354,-1.3,0.,-1.1,-0.5,-1.1090169943749475,-1.,0.,-1.25,-0.6123724356957945,-1.369306393762915,-1.1090169943749475,-1.,-0.55,-1.,0.7071067811865475,-0.7071067811865477,0.,-0.5,0.5,0.,0.7071067811865477,-0.7071067811865475,1.05,0.,1.1,-0.5,0.55,-1.,1.3,0.,1.4012585384440737,-0.535233134659635,1.1090169943749475,-1.,1.1,-0.5,0.6123724356957946,-1.369306393762915,0.,-1.25,0.55,-1.,1.1090169943749475,-1.0};
  CPPUNIT_ASSERT_EQUAL(136,(int)m3->getNodalConnectivity()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(17,(int)m3->getNodalConnectivityIndex()->getNumberOfTuples());
  CPPUNIT_ASSERT(std::equal(expected3,expected3+136,m3->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+17,m3->getNodalConnectivityIndex()->getConstPointer()));
  for(int i=0;i<208;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],m3->getCoords()->getIJ(0,i),1e-12);
  d1->decrRef();
  d2->decrRef();
  m3->decrRef();
  //
  m1->decrRef();
  m2->decrRef();
}
