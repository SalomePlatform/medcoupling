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

