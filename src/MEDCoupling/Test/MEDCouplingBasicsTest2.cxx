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

#include <cmath>
#include <functional>

using namespace ParaMEDMEM;

void MEDCouplingBasicsTest::testGaussPointField1()
{
  const double _a=0.446948490915965;
  const double _b=0.091576213509771;
  const double _p1=0.11169079483905;
  const double _p2=0.0549758718227661;
  const double refCoo1[6]={ 0.,0., 1.,0., 0.,1. };
  const double gsCoo1[12]={ 2*_b-1, 1-4*_b, 2*_b-1, 2.07*_b-1, 1-4*_b,
                            2*_b-1, 1-4*_a, 2*_a-1, 2*_a-1, 1-4*_a, 2*_a-1, 2*_a-1 };
  const double wg1[6]={ 4*_p2, 4*_p2, 4*_p2, 4*_p1, 4*_p1, 4*_p1 };
  std::vector<double> _refCoo1(refCoo1,refCoo1+6);
  std::vector<double> _gsCoo1(gsCoo1,gsCoo1+12);
  std::vector<double> _wg1(wg1,wg1+6);
  //
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_GAUSS_PT,NO_TIME);
  f->setMesh(m);
  CPPUNIT_ASSERT_EQUAL(0,f->getNbOfGaussLocalization());
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_TRI3,_refCoo1,_gsCoo1,_wg1);
  CPPUNIT_ASSERT_THROW(f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_QUAD4,_refCoo1,_gsCoo1,_wg1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_EQUAL(1,f->getNbOfGaussLocalization());
  const double refCoo2[8]={ 0.,0., 1.,0., 1.,1., 0.,1. };
  std::vector<double> _refCoo2(refCoo2,refCoo2+8);
  _gsCoo1.resize(4); _wg1.resize(2);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_QUAD4,_refCoo2,_gsCoo1,_wg1);
  CPPUNIT_ASSERT_EQUAL(2,f->getNbOfGaussLocalization());
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(18,2);
  double *ptr=array->getPointer();
  for(int i=0;i<18*2;i++)
    ptr[i]=(double)(i+1);
  f->setArray(array);
  f->setName("MyFirstFieldOnGaussPoint");
  array->decrRef();
  f->checkCoherency();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(27.,f->getIJK(2,5,0),1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(16.,f->getIJK(1,5,1),1e-14);
  //
  f->clearGaussLocalizations();
  CPPUNIT_ASSERT_EQUAL(0,f->getNbOfGaussLocalization());
  CPPUNIT_ASSERT_THROW(f->checkCoherency(),INTERP_KERNEL::Exception);
  int ids1[4]={0,1,3,4};
  CPPUNIT_ASSERT_THROW(f->setGaussLocalizationOnCells(ids1,ids1+4,_refCoo2,_gsCoo1,_wg1),INTERP_KERNEL::Exception);
  CPPUNIT_ASSERT_EQUAL(0,f->getNbOfGaussLocalization());
  int ids2[2]={0,4};
  f->setGaussLocalizationOnCells(ids2,ids2+2,_refCoo2,_gsCoo1,_wg1);
  CPPUNIT_ASSERT_EQUAL(1,f->getNbOfGaussLocalization());
  CPPUNIT_ASSERT_EQUAL(0,f->getGaussLocalizationIdOfOneCell(0));
  CPPUNIT_ASSERT_THROW(f->getGaussLocalizationIdOfOneCell(1),INTERP_KERNEL::Exception);
  int ids3[2]={1,2};
  f->setGaussLocalizationOnCells(ids3,ids3+2,_refCoo1,_gsCoo1,_wg1);
  CPPUNIT_ASSERT_EQUAL(2,f->getNbOfGaussLocalization());
  CPPUNIT_ASSERT_EQUAL(0,f->getGaussLocalizationIdOfOneCell(0));
  CPPUNIT_ASSERT_EQUAL(1,f->getGaussLocalizationIdOfOneCell(1));
  CPPUNIT_ASSERT_EQUAL(1,f->getGaussLocalizationIdOfOneCell(2));
  CPPUNIT_ASSERT_THROW(f->checkCoherency(),INTERP_KERNEL::Exception);//<- cell 3 has no localization
  int ids4[1]={3};
  std::vector<double> _gsCoo2(_gsCoo1);
  std::vector<double> _wg2(_wg1);
  _gsCoo2[0]=0.8888777776666; _wg2[0]=0.1234567892377;
  f->setGaussLocalizationOnCells(ids4,ids4+1,_refCoo2,_gsCoo2,_wg2);
  CPPUNIT_ASSERT_EQUAL(3,f->getNbOfGaussLocalization());
  std::vector<int> tmpIds;
  f->getCellIdsHavingGaussLocalization(0,tmpIds);
  CPPUNIT_ASSERT_EQUAL(2,(int)tmpIds.size());
  CPPUNIT_ASSERT(std::equal(ids2,ids2+2,tmpIds.begin()));
  CPPUNIT_ASSERT_THROW(f->checkCoherency(),INTERP_KERNEL::Exception);//<- it's always not ok because undelying array not with the good size.
  DataArrayDouble *array2=f->getArray()->substr(0,10);
  f->setArray(array2);
  array2->decrRef();
  f->checkCoherency();//<- here it is OK
  MEDCouplingFieldDouble *f2=f->clone(true);
  CPPUNIT_ASSERT(f->isEqual(f2,1e-14,1e-14));
  MEDCouplingGaussLocalization& gl1=f2->getGaussLocalization(0);
  double tmp=gl1.getGaussCoord(1,1);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.07*_b-1,tmp,1e-14);
  gl1.setGaussCoord(1,1,0.07);
  CPPUNIT_ASSERT(!f->isEqual(f2,1e-14,1e-14));
  gl1.setGaussCoord(1,1,tmp);
  CPPUNIT_ASSERT(f->isEqual(f2,1e-14,1e-14));
  f->decrRef();
  f2->checkCoherency();
  //
  f2->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testGaussPointNEField1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_GAUSS_NE,NO_TIME);
  f->setMesh(m);
  f->setName("MyFirstFieldOnNE");
  f->setDescription("MyDescriptionNE");
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(18,2);
  double *ptr=array->getPointer();
  for(int i=0;i<18*2;i++)
    ptr[i]=(double)(i+7);
  f->setArray(array);
  array->decrRef();
  //
  f->checkCoherency();
  MEDCouplingFieldDouble *f2=f->clone(true);
  CPPUNIT_ASSERT(f->isEqual(f2,1e-14,1e-14));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(21.,f->getIJK(2,0,0),1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(18.,f->getIJK(1,1,1),1e-14);
  f2->decrRef();
  //
  f->decrRef();
  m->decrRef();
}
