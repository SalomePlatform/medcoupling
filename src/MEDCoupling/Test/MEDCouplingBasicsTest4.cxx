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
  dt.appendRepr(std::cout);
  //
  for(std::vector<MEDCouplingFieldDouble *>::iterator it=fs.begin();it!=fs.end();it++)
    (*it)->decrRef();
  fot->decrRef();
}
