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
#include <iterator>

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

void MEDCouplingBasicsTest::testCellOrientation1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  double vec[3]={0.,0.,1.};
  std::vector<int> res1;
  CPPUNIT_ASSERT_THROW(m->are2DCellsNotCorrectlyOriented(vec,false,res1),INTERP_KERNEL::Exception);
  m->changeSpaceDimension(3);
  res1.clear();
  m->are2DCellsNotCorrectlyOriented(vec,false,res1);
  CPPUNIT_ASSERT(res1.empty());
  vec[2]=-1;
  m->are2DCellsNotCorrectlyOriented(vec,false,res1);
  CPPUNIT_ASSERT_EQUAL(5,(int)res1.size());
  res1.clear();
  //
  vec[2]=1.;
  // connectivity inversion
  int *conn=m->getNodalConnectivity()->getPointer();
  int tmp=conn[11];
  conn[11]=conn[12];
  conn[12]=tmp;
  m->are2DCellsNotCorrectlyOriented(vec,false,res1);
  CPPUNIT_ASSERT_EQUAL(1,(int)res1.size());
  CPPUNIT_ASSERT_EQUAL(2,res1[0]);
  res1.clear();
  m->orientCorrectly2DCells(vec,false);
  m->are2DCellsNotCorrectlyOriented(vec,false,res1);
  CPPUNIT_ASSERT(res1.empty());
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  m2->changeSpaceDimension(3);
  CPPUNIT_ASSERT(m->isEqual(m2,1e-12));
  m2->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest::testCellOrientation2()
{
  MEDCouplingUMesh *m1=0;
  MEDCouplingUMesh *m2=build3DExtrudedUMesh_1(m1);
  m1->decrRef();
  std::vector<int> res1;
  m2->arePolyhedronsNotCorrectlyOriented(res1);
  CPPUNIT_ASSERT_EQUAL(6,(int)res1.size());
  m2->orientCorrectlyPolyhedrons();
  res1.clear();
  m2->arePolyhedronsNotCorrectlyOriented(res1);
  CPPUNIT_ASSERT(res1.empty());
  m2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(18,m2->getNumberOfCells());
  int cellIds[3]={0,6,12};
  std::vector<int> cellIds2(cellIds,cellIds+3);
  m2->convertToPolyTypes(cellIds2);
  m2->orientCorrectlyPolyhedrons();
  res1.clear();
  m2->arePolyhedronsNotCorrectlyOriented(res1);
  CPPUNIT_ASSERT(res1.empty());
  MEDCouplingFieldDouble *f2=m2->getMeasureField(false);
  //Test to check global reverse in MEDCouplingUMesh::tryToCorrectPolyhedronOrientation
  MEDCouplingUMesh *m3=build2DTargetMesh_1();
  double vec[3]={0.,0.,-1.};//<- important for the test
  m3->changeSpaceDimension(3);
  const int ids1[5]={0,1,2,3,4};
  std::vector<int> ids2(ids1,ids1+5);
  m3->convertToPolyTypes(ids2);
  m3->orientCorrectly2DCells(vec,false);
  MEDCouplingUMesh *m4=buildCU1DMesh_U();
  m4->changeSpaceDimension(3);
  double center[3]={0.,0.,0.};
  double vector[3]={0.,1.,0.};
  m4->rotate(center,vector,-M_PI/2.);
  MEDCouplingUMesh *m5=m3->buildExtrudedMesh(m4,0);
  res1.clear();
  m5->arePolyhedronsNotCorrectlyOriented(res1);
  CPPUNIT_ASSERT_EQUAL(15,(int)res1.size());
  m5->orientCorrectlyPolyhedrons();
  res1.clear();
  m5->arePolyhedronsNotCorrectlyOriented(res1);
  CPPUNIT_ASSERT(res1.empty());
  MEDCouplingFieldDouble *f3=m5->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(15,f3->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f3->getNumberOfComponents());
  const double *f3Ptr=f3->getArray()->getConstPointer();
  const double expected1[15]={
    0.075,0.0375,0.0375,0.075,0.075,
    0.1125,0.05625,0.05625,0.1125,0.1125,
    0.0625,0.03125,0.03125,0.0625,0.0625
  };
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::abs(expected1[i]),f3Ptr[i],1e-12);
  f3->decrRef();
  DataArrayDouble *f4=m5->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(15,f4->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,f4->getNumberOfComponents());
  const double *f4Ptr=f4->getConstPointer();
  const double expected2[45]={
    -0.05,-0.05,0.15, 0.3666666666666667,-0.13333333333333333,0.15, 0.53333333333333333,0.033333333333333333,0.15, -0.05,0.45,0.15, 0.45,0.45,0.15,
    -0.05,-0.05,0.525, 0.3666666666666667,-0.13333333333333333,0.525, 0.53333333333333333,0.033333333333333333,0.525, -0.05,0.45,0.525, 0.45,0.45,0.525,
    -0.05,-0.05,0.875, 0.3666666666666667,-0.13333333333333333,0.875, 0.53333333333333333,0.033333333333333333,0.875, -0.05,0.45,0.875, 0.45,0.45,0.875
  };
  for(int i=0;i<45;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f4Ptr[i],1e-12);
  f4->decrRef();
  m5->decrRef();
  m3->decrRef();
  m4->decrRef();
  //
  f2->decrRef();
  m2->decrRef();
}

/*!
 * This test check polyhedron true barycenter computation. 
 */
void MEDCouplingBasicsTest::testPolyhedronBarycenter()
{
  int connN[]={0,3,2,1, -1, 4,5,6,7, -1, 0,4,7,3, -1, 3,7,6,2, -1, 2,6,5,1, -1, 1,5,4,0};
  double coords[]={0.,0.,0., 1.,0.,0., 1.,1.,0., 0.,1.,0., 0.,0.,1., 1.,0.,1., 1.,1.,1., 0.,1.,1., 0.5, 0.5, 0.5};
  MEDCouplingUMesh *meshN=MEDCouplingUMesh::New();
  meshN->setName("ForBary");
  meshN->setMeshDimension(3);
  meshN->allocateCells(4);
  meshN->insertNextCell(INTERP_KERNEL::NORM_POLYHED,29,connN);
  meshN->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(9,3);
  std::copy(coords,coords+27,myCoords->getPointer());
  meshN->setCoords(myCoords);
  myCoords->decrRef();
  meshN->checkCoherency();
  //
  std::vector<int> res1;
  meshN->arePolyhedronsNotCorrectlyOriented(res1);
  meshN->orientCorrectlyPolyhedrons();
  CPPUNIT_ASSERT(res1.empty());
  const double *ref,*daPtr;
  DataArrayDouble *da=meshN->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(1,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da->getNumberOfComponents());
  daPtr=da->getConstPointer();
  ref=meshN->getCoords()->getConstPointer()+24;
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ref[i],daPtr[i],1e-12);
  da->decrRef();
  //
  const double center[]={0.,0.,0.};
  const double vec[]={0.,2.78,0.};
  da=meshN->getBarycenterAndOwner();
  daPtr=da->getConstPointer();
  ref=meshN->getCoords()->getConstPointer()+24;
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ref[i],daPtr[i],1e-12);
  da->decrRef();
  //
  meshN->rotate(center,vec,M_PI/7.);
  meshN->translate(vec);
  da=meshN->getBarycenterAndOwner();
  daPtr=da->getConstPointer();
  ref=meshN->getCoords()->getConstPointer()+24;
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ref[i],daPtr[i],1e-12);
  da->decrRef();
  //
  const double center2[]={1.12,3.45,6.78};
  const double vec2[]={4.5,9.3,2.8};
  meshN->rotate(center2,vec2,M_E);
  meshN->translate(vec2);
  da=meshN->getBarycenterAndOwner();
  daPtr=da->getConstPointer();
  ref=meshN->getCoords()->getConstPointer()+24;
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(ref[i],daPtr[i],1e-10);
  da->decrRef();
  //
  meshN->decrRef();
}

void MEDCouplingBasicsTest::testNormL12Integ1D()
{
  MEDCouplingUMesh *m1=build1DTargetMesh_3();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(m1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(m1->getNumberOfCells(),3);
  const double arr[12]={-5.23,15.45,-25.56,6.67,-16.78,26.89,-7.91,17.23,-27.43,8.21,-18.63,28.72};
  std::copy(arr,arr+12,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //
  const double *ptr;
  DataArrayDouble *f3=m1->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(4,f3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f3->getNumberOfComponents());
  double expected9[4]={0.75,5.105,0.8,5.155};
  ptr=f3->getConstPointer();
   for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected9[i],ptr[i],1e-12);
  f3->decrRef();
  //
  MEDCouplingFieldDouble *f2=m1->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(4,f2->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  double expected1[4]={0.5,0.21,-0.6,-0.31};
  ptr=f2->getArray()->getConstPointer();
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],ptr[i],1e-12);
  f2->decrRef();
  double expected2[4]={0.5,0.21,0.6,0.31};
  f2=m1->getMeasureField(true);
  ptr=f2->getArray()->getConstPointer();
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],ptr[i],1e-12);
  f2->decrRef();
  //integral
  double res[3];
  f1->integral(false,res);
  double expected3[3]={0.9866,-0.3615,0.4217};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],res[i],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[0],f1->integral(0,false),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[1],f1->integral(1,false),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[2],f1->integral(2,false),1e-12);
  f1->integral(true,res);
  double expected4[3]={-3.4152,8.7639,-14.6879};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected4[i],res[i],1e-12);
  //normL1
  f1->normL1(res);
  double expected5[3]={11.3068,27.3621,43.7881};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[i],res[i],1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[0],f1->normL1(0),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[1],f1->normL1(1),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected5[2],f1->normL1(2),1e-12);
  //normL2
  f1->normL2(res);
  double expected7[3]={9.0252562290496776, 21.545259176904789, 34.433193070059595};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected7[i],res[i],1e-9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected7[0],f1->normL2(0),1e-9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected7[1],f1->normL2(1),1e-9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected7[2],f1->normL2(2),1e-9);
  //buildMeasureField
  MEDCouplingFieldDouble *f4=f1->buildMeasureField(false);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.2,f4->accumulate(0),1e-12);
  f4->decrRef();
  f4=f1->buildMeasureField(true);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.62,f4->accumulate(0),1e-12);
  f4->decrRef();
  //
  f1->decrRef();
  m1->decrRef();
  // Testing with 2D Curve
  m1=build2DCurveTargetMesh_3();
  f2=m1->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(4,f2->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  ptr=f2->getArray()->getConstPointer();
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(2.)*expected2[i],ptr[i],1e-12);
  f2->decrRef();
  f2=m1->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(4,f2->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  ptr=f2->getArray()->getConstPointer();
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i]*sqrt(2.),ptr[i],1e-12);
  f2->decrRef();
  //bary
  f3=m1->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(4,f3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f3->getNumberOfComponents());
  double expected10[8]={0.75,0.75,5.105,5.105,0.8,0.8,5.155,5.155};
  ptr=f3->getConstPointer();
   for(int i=0;i<8;i++)
     CPPUNIT_ASSERT_DOUBLES_EQUAL(expected10[i],ptr[i],1e-12);
  f3->decrRef();
  //
  f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(m1);
  array=DataArrayDouble::New();
  array->alloc(m1->getNumberOfCells(),3);
  std::copy(arr,arr+12,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->integral(false,res);
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(2.)*expected4[i],res[i],1e-12);
  f1->integral(true,res);
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(2.)*expected4[i],res[i],1e-12);
  f1->normL1(res);
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(2.)*expected5[i],res[i],1e-12);
  f1->normL2(res);
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(sqrt(sqrt(2.))*expected7[i],res[i],1e-12);
  //
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest::testAreaBary2D()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_3();
  MEDCouplingFieldDouble *f1=m1->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(10,f1->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  double expected1[10]={-0.5,-1,-1.5,-0.5,-1,  0.5,1,1.5,0.5,1};
  const double *ptr=f1->getArray()->getConstPointer();
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],ptr[i],1e-12);
  f1->decrRef();
  f1=m1->getMeasureField(true);
  ptr=f1->getArray()->getConstPointer();
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::abs(expected1[i]),ptr[i],1e-12);
  f1->decrRef();
  DataArrayDouble *f2=m1->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(10,f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f2->getNumberOfComponents());
  double expected2[20]={
    0.5,0.3333333333333333,0.5,0.5,0.5,0.77777777777777777,0.5,0.3333333333333333,0.5,0.5,
    0.5,0.3333333333333333,0.5,0.5,0.5,0.77777777777777777,0.5,0.3333333333333333,0.5,0.5,
  };
  ptr=f2->getConstPointer();
  for(int i=0;i<20;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],ptr[i],1e-12);
  f2->decrRef();
  m1->changeSpaceDimension(3);
  f1=m1->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(10,f1->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  ptr=f1->getArray()->getConstPointer();
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(std::abs(expected1[i]),ptr[i],1e-12);
  f1->decrRef();
  f2=m1->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(10,f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,f2->getNumberOfComponents());
  ptr=f2->getConstPointer();
  double expected3[30]={
    0.5,0.3333333333333333,0.,0.5,0.5,0.,0.5,0.77777777777777777,0.,0.5,0.3333333333333333,0.,0.5,0.5,0.,
    0.5,0.3333333333333333,0.,0.5,0.5,0.,0.5,0.77777777777777777,0.,0.5,0.3333333333333333,0.,0.5,0.5,0.
  };
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],ptr[i],1e-12);
  f2->decrRef();
  m1->decrRef();
}

/*!
 * This test check polyhedron true barycenter computation 2. 
 */
void MEDCouplingBasicsTest::testAreaBary3D()
{
  double coords [] = { 0.241310763507 , 0.0504777305619 , 0.0682283524903 , 0.252501053866 , -0.0625176732937 , 0.137272639894 ,
                       0.152262663601 , 0.241816569527 , 0.133812556197 , 0.18047750211 , -0.0789949051358 , 0.339098173401 ,
                       0.151741971857 , 0.238885278571 , 0.137715037333 , 0.242532155481 , -0.0928169086456 , 0.0678043417367 ,
                       0.240941965335 , -0.015461491464 , 0.0617186345825 , 0.24127650112 , 0.0499427876717 , 0.0679634099148 ,
                       -0.145828917428 , 0.206291632565 , 0.0310071927543 , 0.0125651775307 , 0.266262085828 , 0.105228430543 ,
                       -0.0994066533286 , 0.233224271238 , 0.0572213839567 , -0.0951345338317 , 0.234819509426 , 0.0592126284538 ,
                       0.136580574205 , -0.205486212579 , 0.0572866072014 , 0.0637270784978 , -0.168886355238 , 0.446614057077 ,
                       0.041337157151 , -0.213402568198 , 0.372407095999 , 0.0411601970268 , -0.202387875756 , 0.411334979491 ,
                       -0.108355701857 , 0.193636239335 , 0.204886756738 , 0.00639779029829 , 0.155296981517 , 0.252585892979 ,
                       0.0262473111702 , -0.112919732543 , 0.424286639249 ,-0.224103052733 , -0.139430015438 , -0.0122352295701 ,
                       -0.0312760589481 , -0.274272003594 , 0.0323959636568 , -0.166663422532 , -0.217754445175 , 0.00392109070364 ,
                       -0.30586619777 , -0.0475168041091 , -0.0144585228182 , -0.280881480586 , 0.135571293538 , 0.00623923647986 ,
                       -0.25548538234 , 0.156819217766 , 0.0645277879769 , -0.131567009284 , 0.184133752309 , 0.206021802753 ,
                       -0.196204010965 , 0.151602971681 , 0.212974777736 , -0.183713879463 , 0.0802946639531 , 0.260115662599 ,
                       -0.244241178767 , -0.0738873389604 , 0.144590565817 , -0.155804057829 , -0.164892720025 , 0.210613950558 ,
                       -0.170950800428 , -0.215099334026 , 0.00610122860092 , -0.30552634869 , -0.0490020791904 , -0.0132786533145 ,
                       0.271831011884 , 0.15105657296 , 0.0230534827908 , 0.281919192283 , 0.0898544306288 , -0.0625201489143 ,
                       0.260240727276 , -0.0120688706637 , -0.0532316588626 , 0.244947737722 , 0.0197984684293 , 0.0309341209233 ,
                       0.23439631578 , 0.229825279875 , 0.0508520585381 , 0.160921316875 , 0.265078502128 , 0.121716560626 ,
                       -0.315088694175 , 0.0747700471918 , -0.245836615071 , -0.327728781776 , 0.0857114674649 , -0.239431905957 ,
                       -0.308385460634 , 0.145142997084 , -0.149886828433 , 0.0488236045164 , 0.309462801914 , 0.0849169148265 ,
                       -0.0244964803395 , 0.33145611751 , -0.0476415818061 , 0.0060567994229 , 0.32418412014 , 0.0367779543812 ,
                       -0.0950221448063 , 0.236675326003 , 0.0572594453983 , 0.248723023186 , 0.0886648784791 , -0.176629430538 ,
                       0.116796984 , 0.256596599567 , -0.292863523603 , 0.118024552914 , 0.229154257843 , -0.34233232501 ,
                       0.217507892549 , -0.0417822335742 , -0.176771782888 , -0.224429321304 , 0.0125595300114 , -0.362064725588 ,
                       0.0937301100955 , -0.0500824832657 , -0.299713548444 , -0.244162220397 , 0.0383853931293 , -0.389856984411 ,
                       -0.0281989366102 , 0.097392811563 , -0.458244577284 , -0.385010847162 , 0.10122766194 , -0.140052859922 ,
                       -0.377936358012 , 0.110875172128 , -0.176207095463 , 0.244483045556 , -0.0991073977045 , 0.0575134372934 ,
                       0.262605120167 , -0.100243191645 , -0.0495620806935 , 0.240306880972 , -0.136153701579 , -0.114745281696 ,
                       0.215763176129 , -0.0836766059189 , -0.183249640616 , 0.237870396603 , -0.132449578286 , -0.121598854639 ,
                       -0.0637683083097 , -0.27921020214 , -0.149112321992 , -0.0856211014977 , -0.2973233473 , -0.0446878139589 ,
                       0.104675342288 , -0.0625908305324 , -0.290346256534 , 0.0248264249186 , -0.247797708548 , -0.165830884019 ,
                       0.0719302438309 , -0.178468260473 , -0.211432157345 , 0.142871843159 , -0.208769948542 , 0.0454101128246 ,
                       0.167803379307 , -0.207851396623 , -0.088802726124 , 0.12868717152 , -0.230920439715 , 0.00760508389036 ,
                       -0.0372812069535 , -0.286740286332 , 0.00963701291166 };

  int connN [] = { /*polyhedron 0*/
    0 , 1 , 3 , 4 , 2 , -1 , 1 , 5 , 6 , 7 , 0 , -1 , 0 , 7 , 8 , 10 , 11 , 9 , 2 , -1 , 1 , 5 , 12 , 14 , 15 , 13 , 3 , -1 , 16 , 9 , 2 , 4 , 17 , -1
    , 4 , 3 , 13 , 18 , 17 , -1 , 5 , 6 , 19 , 21 , 20 , 12 , -1 , 6 , 7 , 8 , 23 , 22 , 19 , -1 , 23 , 24 , 10 , 8 , -1 , 25 , 11 , 9 , 16 , -1
    , 24 , 26 , 25 , 11 , 10 , -1 , 12 , 14 , 20 , -1 , 27 , 28 , 29 , 15 , 13 , 18 , -1 , 14 , 15 , 29 , 30 , 21 , 20 , -1 , 26 , 27 , 18 , 17 , 16 , 25 , -1
    , 22 , 19 , 21 , 30 , 31 , -1 , 22 , 31 , 28 , 27 , 26 , 24 , 23 , -1 , 31 , 30 , 29 , 28,
    /* polyhedron 1*/
    0 , 7 , 8 , 10 , 11 , 9 , 2 , -1 , 32 , 0 , 7 , 35 , 34 , 33 , -1 , 32 , 0 , 2 , 37 , 36 , -1 , 35 , 7 , 8 , 40 , 39 , 38 , -1
    , 2 , 37 , 41 , 9 , -1 , 40 , 8 , 10 , 44 , 43 , 42 , -1 , 41 , 9 , 11 , 44 , 43 , -1 , 44 , 11 , 10 , -1 , 32 , 33 , 45 , 47 , 46 , 36 , -1
    , 33 , 34 , 48 , 45 , -1 , 35 , 34 , 48 , 50 , 49 , 38 , -1 , 41 , 43 , 42 , 46 , 36 , 37 , -1 , 38 , 39 , 51 , 49 , -1
    , 39 , 40 , 42 , 46 , 47 , 52 , 51 , -1 , 45 , 47 , 52 , 50 , 48 , -1 , 52 , 51 , 49 , 50,
    /* polyhedron 2*/
    6 , 7 , 8 , 23 , 22 , 19 , -1 , 6 , 35 , 7 , -1 , 6 , 35 , 38 , 19 , -1 , 35 , 7 , 8 , 40 , 39 , 38 , -1 , 53 , 22 , 19 , 38 , 39 , 54 , -1
    , 23 , 53 , 54 , 40 , 8 , -1 , 53 , 22 , 23 , -1 , 39 , 54 , 40,
    /*polyhedron 3*/
    35 , 34 , 48 , 50 , 49 , 38 , -1 , 6 , 35 , 34 , 56 , 55 , 5 , -1 , 6 , 35 , 38 , 19 , -1 , 34 , 56 , 57 , 59 , 58 , 48 , -1
    , 60 , 61 , 21 , 19 , 38 , 49 , -1 , 62 , 50 , 48 , 58 , -1 , 60 , 63 , 64 , 62 , 50 , 49 , -1 , 5 , 6 , 19 , 21 , 20 , 12 , -1
    , 55 , 5 , 12 , 65 , -1 , 66 , 67 , 65 , 55 , 56 , 57 , -1 , 63 , 66 , 57 , 59 , 64 , -1 , 64 , 62 , 58 , 59 , -1
    , 60 , 63 , 66 , 67 , 68 , 61 , -1 , 61 , 68 , 20 , 21 , -1 , 67 , 68 , 20 , 12 , 65};

  double barys[]={ -0.0165220465527 , -0.0190922868195 , 0.158882733414 ,
                   0.0287618656076 , 0.135874379934 , -0.14601588119 ,
                   -0.147128055553 , 0.0465995097041 , -0.049391174453 ,
                   -0.00142506732317 , -0.0996953090351 , -0.115159183132 };
  MEDCouplingUMesh *meshN=MEDCouplingUMesh::New();
  meshN->setName("ForBary");
  meshN->setMeshDimension(3);
  meshN->allocateCells(4);
  meshN->insertNextCell(INTERP_KERNEL::NORM_POLYHED,113,connN);
  meshN->insertNextCell(INTERP_KERNEL::NORM_POLYHED,99,connN+113);
  meshN->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,connN+212);
  meshN->insertNextCell(INTERP_KERNEL::NORM_POLYHED,92,connN+255);
  meshN->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(69,3);
  std::copy(coords,coords+207,myCoords->getPointer());
  meshN->setCoords(myCoords);
  myCoords->decrRef();
  meshN->checkCoherency();
  std::vector<int> res1;
  meshN->arePolyhedronsNotCorrectlyOriented(res1);
  meshN->orientCorrectlyPolyhedrons();
  res1.clear();
  meshN->arePolyhedronsNotCorrectlyOriented(res1);
  CPPUNIT_ASSERT(res1.empty());
  //
  DataArrayDouble *da=meshN->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(4,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da->getNumberOfComponents());
  const double *daPtr=da->getConstPointer();
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(barys[i],daPtr[i],1e-12);
  da->decrRef();
  //
  meshN->decrRef();
}

void MEDCouplingBasicsTest::testRenumberCellsForFields()
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
  const int renumber1[5]={3,1,0,4,2};
  double res[3];
  const double loc[]={-0.05,-0.05, 0.55,-0.25, 0.55,0.15, -0.05,0.45, 0.45,0.45};
  for(int j=0;j<5;j++)
    {
      f->getValueOn(loc+2*j,res);
      for(int i=0;i<3;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(values1[i+3*j],res[i],1e-12);
    }
  f->renumberCells(renumber1,false);
  const double *ptr=f->getArray()->getConstPointer();
  const double expected1[15]={9.,109.,10009.,8.,108.,10008.,11.,111.,10011.,7.,107.,10007.,10.,110.,10010.};
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],ptr[i],1e-12);
  //check that fields remains the same geometrically
  for(int j=0;j<5;j++)
    {
      f->getValueOn(loc+2*j,res);
      for(int i=0;i<3;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(values1[i+3*j],res[i],1e-12);
    }
  f->decrRef();
  //On gauss
  f=MEDCouplingFieldDouble::New(ON_GAUSS_PT,NO_TIME);
  f->setMesh(m);
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
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_TRI3,_refCoo1,_gsCoo1,_wg1);
  const double refCoo2[8]={ 0.,0., 1.,0., 1.,1., 0.,1. };
  std::vector<double> _refCoo2(refCoo2,refCoo2+8);
  _gsCoo1.resize(4); _wg1.resize(2);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_QUAD4,_refCoo2,_gsCoo1,_wg1);
  arr=DataArrayDouble::New();
  arr->alloc(18,2);
  const double values2[36]={1.,1001.,2.,1002., 11.,1011.,12.,1012.,13.,1013.,14.,1014.,15.,1015.,16.,1016., 21.,1021.,22.,1022.,23.,1023.,24.,1024.,25.,1025.,26.,1026., 31.,1031.,32.,1032., 41.,1041.,42.,1042.};
  std::copy(values2,values2+36,arr->getPointer());
  f->setArray(arr);
  arr->decrRef();
  f->checkCoherency();
  MEDCouplingFieldDouble *fCpy=f->clone(true);
  CPPUNIT_ASSERT(f->isEqual(fCpy,1e-12,1e-12));
  f->renumberCells(renumber1,false);
  CPPUNIT_ASSERT(!f->isEqual(fCpy,1e-12,1e-12));
  double expected2[36]={21.,1021.,22.,1022.,23.,1023.,24.,1024.,25.,1025.,26.,1026., 11.,1011.,12.,1012.,13.,1013.,14.,1014.,15.,1015.,16.,1016., 41.,1041.,42.,1042., 1.,1001.,2.,1002., 31.,1031.,32.,1032.};
  ptr=f->getArray()->getConstPointer();
  for(int i=0;i<36;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],ptr[i],1e-12);
  const int renumber2[5]={2,1,4,0,3};//reverse renumber1
  f->renumberCells(renumber2,false);
  CPPUNIT_ASSERT(f->isEqual(fCpy,1e-12,1e-12));
  fCpy->decrRef();
  f->decrRef();
  //GaussNE
  f=MEDCouplingFieldDouble::New(ON_GAUSS_NE,NO_TIME);
  f->setMesh(m);
  arr=DataArrayDouble::New();
  arr->alloc(18,2);
  const double values3[36]={1.,1001.,2.,1002.,3.,1003.,4.,1004., 11.,1011.,12.,1012.,13.,1013., 21.,1021.,22.,1022.,23.,1023., 31.,1031.,32.,1032.,33.,1033.,34.,1034., 41.,1041.,42.,1042.,43.,1043.,44.,1044.};
  std::copy(values3,values3+36,arr->getPointer());
  f->setArray(arr);
  arr->decrRef();
  f->checkCoherency();
  fCpy=f->clone(true);
  CPPUNIT_ASSERT(f->isEqual(fCpy,1e-12,1e-12));
  f->renumberCells(renumber1,false);
  CPPUNIT_ASSERT(!f->isEqual(fCpy,1e-12,1e-12));
  double expected3[36]={21.,1021.,22.,1022.,23.,1023.,11.,1011.,12.,1012.,13.,1013.,41.,1041.,42.,1042.,43.,1043.,44.,1044.,1.,1001.,2.,1002.,3.,1003.,4.,1004.,31.,1031.,32.,1032.,33.,1033.,34.,1034.};
  ptr=f->getArray()->getConstPointer();
  for(int i=0;i<36;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],ptr[i],1e-12);
  f->renumberCells(renumber2,false);//perform reverse operation of renumbering to check that the resulting field is equal.
  CPPUNIT_ASSERT(f->isEqual(fCpy,1e-12,1e-12));
  fCpy->decrRef();
  f->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest::testRenumberNodesForFields()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_NODES,NO_TIME);
  f->setMesh(m);
  DataArrayDouble *arr=DataArrayDouble::New();
  int nbOfNodes=m->getNumberOfNodes();
  arr->alloc(nbOfNodes,3);
  f->setArray(arr);
  arr->decrRef();
  const double values1[27]={7.,107.,10007.,8.,108.,10008.,9.,109.,10009.,10.,110.,10010.,11.,111.,10011.,12.,112.,10012.,13.,113.,10013.,14.,114.,10014.,15.,115.,10015.};
  std::copy(values1,values1+27,arr->getPointer());
  f->checkCoherency();
  const int renumber1[9]={0,4,1,3,5,2,6,7,8};
  double res[3];
  const double loc[]={0.5432,-0.2432, 0.5478,0.1528};
  const double expected1[6]={9.0272, 109.0272, 10009.0272, 11.4124,111.4124,10011.4124};
  for(int j=0;j<2;j++)
    {
      f->getValueOn(loc+2*j,res);
      for(int i=0;i<3;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i+3*j],res[i],1e-12);
    }
  MEDCouplingFieldDouble *fCpy=f->clone(true);
  CPPUNIT_ASSERT(f->isEqual(fCpy,1e-12,1e-12));
  f->renumberNodes(renumber1);
  CPPUNIT_ASSERT(!f->isEqual(fCpy,1e-12,1e-12));
  for(int j=0;j<2;j++)
    {
      f->getValueOn(loc+2*j,res);
      for(int i=0;i<3;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i+3*j],res[i],1e-12);
    }
  const double expected2[27]={7.,107.,10007.,9.,109.,10009.,12.,112.,10012.,10.,110.,10010.,8.,108.,10008.,11.,111.,10011.,13.,113.,10013.,14.,114.,10014.,15.,115.,10015.};
  for(int i=0;i<27;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f->getArray()->getConstPointer()[i],1e-12);
  const int renumber2[9]={0,2,5,3,1,4,6,7,8};//reverse of renumber2
  f->renumberNodes(renumber2);
  CPPUNIT_ASSERT(f->isEqual(fCpy,1e-12,1e-12));
  fCpy->decrRef();
  //
  m->decrRef();
  f->decrRef();
}

void MEDCouplingBasicsTest::testConvertQuadraticCellsToLinear()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_3();
  mesh->checkCoherency();
  const std::set<INTERP_KERNEL::NormalizedCellType>& types=mesh->getAllTypes();
  CPPUNIT_ASSERT_EQUAL(5,(int)types.size());
  INTERP_KERNEL::NormalizedCellType expected1[5]={INTERP_KERNEL::NORM_POLYGON, INTERP_KERNEL::NORM_TRI3, INTERP_KERNEL::NORM_QUAD4, INTERP_KERNEL::NORM_TRI6, INTERP_KERNEL::NORM_QUAD8};
  std::set<INTERP_KERNEL::NormalizedCellType> expected1Bis(expected1,expected1+5);
  CPPUNIT_ASSERT(expected1Bis==types);
  CPPUNIT_ASSERT(mesh->isPresenceOfQuadratic());
  CPPUNIT_ASSERT_EQUAL(62,mesh->getMeshLength());
  MEDCouplingFieldDouble *f1=mesh->getMeasureField(false);
  //
  mesh->convertQuadraticCellsToLinear();
  CPPUNIT_ASSERT(!mesh->isPresenceOfQuadratic());
  //
  mesh->checkCoherency();
  MEDCouplingFieldDouble *f2=mesh->getMeasureField(false);
  CPPUNIT_ASSERT(f1->getArray()->isEqual(*f2->getArray(),1e-12));
  CPPUNIT_ASSERT_EQUAL(48,mesh->getMeshLength());
  const std::set<INTERP_KERNEL::NormalizedCellType>& types2=mesh->getAllTypes();
  CPPUNIT_ASSERT_EQUAL(3,(int)types.size());
  INTERP_KERNEL::NormalizedCellType expected2[3]={INTERP_KERNEL::NORM_POLYGON, INTERP_KERNEL::NORM_TRI3, INTERP_KERNEL::NORM_QUAD4};
  std::set<INTERP_KERNEL::NormalizedCellType> expected2Bis(expected2,expected2+3);
  CPPUNIT_ASSERT(expected2Bis==types2);
  //
  f1->decrRef();
  f2->decrRef();
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testCheckGeoEquivalWith()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_3();
  MEDCouplingUMesh *mesh2=build2DTargetMesh_3();
  DataArrayInt *cellCor,*nodeCor;
  //First test mesh1
  mesh1->checkGeoEquivalWith(mesh1,0,1e-12,cellCor,nodeCor);//deepEqual
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  mesh1->checkGeoEquivalWith(mesh1,1,1e-12,cellCor,nodeCor);//fastEqual
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  mesh1->checkGeoEquivalWith(mesh1,10,1e-12,cellCor,nodeCor);//deepEqual with geo permutations
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  //Second test mesh1 and mesh2 are 2 different meshes instance
  mesh1->checkGeoEquivalWith(mesh2,0,1e-12,cellCor,nodeCor);//deepEqual
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  mesh1->checkGeoEquivalWith(mesh2,1,1e-12,cellCor,nodeCor);//fastEqual
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  mesh1->checkGeoEquivalWith(mesh2,10,1e-12,cellCor,nodeCor);//deepEqual with geo permutations
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  //Third test : cell permutation by keeping the first the middle and the last as it is.
  const int renum[]={0,2,1,3,4,5,6,8,7,9};
  mesh2->renumberCells(renum,false);
  CPPUNIT_ASSERT_THROW(mesh1->checkGeoEquivalWith(mesh2,0,1e-12,cellCor,nodeCor),INTERP_KERNEL::Exception);//deepEqual fails
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  mesh1->checkGeoEquivalWith(mesh2,1,1e-12,cellCor,nodeCor);//fastEqual do not see anything
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  mesh1->checkGeoEquivalWith(mesh2,10,1e-12,cellCor,nodeCor);//deepEqual with geo permutations
  CPPUNIT_ASSERT(cellCor);
  CPPUNIT_ASSERT_EQUAL(10,cellCor->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,cellCor->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(renum,renum+10,cellCor->getConstPointer()));
  CPPUNIT_ASSERT(nodeCor==0);
  cellCor->decrRef();
  cellCor=0;
  CPPUNIT_ASSERT(nodeCor==0);
  //4th test : cell and node permutation by keeping the first the middle and the last as it is.
  mesh2->decrRef();
  mesh2=build2DTargetMesh_3();
  const int renum2[]={0,2,1,3,4,5,6,8,7,9,10};
  mesh2->renumberCells(renum,false);
  mesh2->renumberNodes(renum2,11);
  CPPUNIT_ASSERT_THROW(mesh1->checkGeoEquivalWith(mesh2,0,1e-12,cellCor,nodeCor),INTERP_KERNEL::Exception);//deepEqual fails
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  mesh1->checkGeoEquivalWith(mesh2,1,1e-12,cellCor,nodeCor);//fastEqual do not see anything
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  mesh1->checkGeoEquivalWith(mesh2,10,1e-12,cellCor,nodeCor);//deepEqual with geo permutations
  CPPUNIT_ASSERT(cellCor);
  CPPUNIT_ASSERT_EQUAL(10,cellCor->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,cellCor->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(renum,renum+10,cellCor->getConstPointer()));
  CPPUNIT_ASSERT(nodeCor);
  CPPUNIT_ASSERT_EQUAL(11,nodeCor->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,nodeCor->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(renum2,renum2+11,nodeCor->getConstPointer()));
  cellCor->decrRef();
  cellCor=0;
  nodeCor->decrRef();
  nodeCor=0;
  //5th test : modification of the last cell to check fastCheck detection.
  mesh2->decrRef();
  mesh2=build2DTargetMesh_3();
  const int renum3[]={0,2,1,3,4,5,6,8,9,7};
  mesh2->renumberCells(renum3,false);
  mesh2->renumberNodes(renum2,11);
  bool isExcep=false;
  try { mesh1->checkGeoEquivalWith(mesh2,0,1e-12,cellCor,nodeCor);//deepEqual fails
  }
  catch(INTERP_KERNEL::Exception& e) { isExcep=true; }
  CPPUNIT_ASSERT(isExcep); isExcep=false;
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  try { mesh1->checkGeoEquivalWith(mesh2,1,1e-12,cellCor,nodeCor);//fastEqual has detected something
  }
  catch(INTERP_KERNEL::Exception& e) { isExcep=true; }
  CPPUNIT_ASSERT(isExcep); isExcep=false;
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor==0);
  mesh2->checkGeoEquivalWith(mesh1,10,1e-12,cellCor,nodeCor);//deepEqual with geo permutations
  CPPUNIT_ASSERT(cellCor);
  CPPUNIT_ASSERT_EQUAL(10,cellCor->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,cellCor->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(renum3,renum3+10,cellCor->getConstPointer()));
  CPPUNIT_ASSERT(nodeCor);
  CPPUNIT_ASSERT_EQUAL(11,nodeCor->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,nodeCor->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(renum2,renum2+11,nodeCor->getConstPointer()));
  cellCor->decrRef();
  cellCor=0;
  nodeCor->decrRef();
  nodeCor=0;
  //
  mesh1->decrRef();
  mesh2->decrRef();
}

void MEDCouplingBasicsTest::testCheckGeoEquivalWith2()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_4();
  MEDCouplingUMesh *mesh2=build2DTargetMesh_1();
  DataArrayInt *cellCor,*nodeCor;
  mesh1->checkGeoEquivalWith(mesh2,10,1e-12,cellCor,nodeCor);
  CPPUNIT_ASSERT(cellCor==0);
  CPPUNIT_ASSERT(nodeCor!=0);
  const int expected1[9]={0, 1, 3, 4, 5, 6, 7, 8, 9};
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],nodeCor->getIJ(i,0));
  nodeCor->decrRef();
  //
  mesh1->decrRef();
  mesh2->decrRef();
}

void MEDCouplingBasicsTest::testCopyTinyStringsFromOnFields()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  int nbOfCells=m->getNumberOfCells();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,LINEAR_TIME);
  f->setMesh(m);
  f->setName("a");
  f->setDescription("b");
  DataArrayDouble *a1=DataArrayDouble::New();
  a1->alloc(nbOfCells,2);
  a1->fillWithZero();
  a1->setInfoOnComponent(0,"c");
  a1->setInfoOnComponent(1,"d");
  DataArrayDouble *a2=a1->deepCpy();
  a2->setInfoOnComponent(0,"e");
  a2->setInfoOnComponent(1,"f");
  f->setArray(a1);
  f->setEndArray(a2);
  f->setEndTime(3.,3,4);
  a2->decrRef();
  a1->decrRef();
  m->setName("g");
  m->getCoords()->setInfoOnComponent(0,"h");
  m->getCoords()->setInfoOnComponent(1,"i");
  m->getCoords()->setInfoOnComponent(2,"j");
  //
  f->checkCoherency();
  MEDCouplingFieldDouble *f2=f->clone(true);
  CPPUNIT_ASSERT(f2->isEqual(f,1e-12,1e-12));
  f2->setName("smth");
  CPPUNIT_ASSERT(!f2->isEqual(f,1e-12,1e-12));
  f2->copyTinyStringsFrom(f);
  CPPUNIT_ASSERT(f2->isEqual(f,1e-12,1e-12));
  f2->setDescription("GGG");
  CPPUNIT_ASSERT(!f2->isEqual(f,1e-12,1e-12));
  f2->copyTinyStringsFrom(f);
  CPPUNIT_ASSERT(f2->isEqual(f,1e-12,1e-12));
  f2->getArray()->setInfoOnComponent(0,"mmmm");
  CPPUNIT_ASSERT(!f2->isEqual(f,1e-12,1e-12));
  f2->copyTinyStringsFrom(f);
  CPPUNIT_ASSERT(f2->isEqual(f,1e-12,1e-12));
  f2->getEndArray()->setInfoOnComponent(1,"mmmm");
  CPPUNIT_ASSERT(!f2->isEqual(f,1e-12,1e-12));
  f2->copyTinyStringsFrom(f);
  CPPUNIT_ASSERT(f2->isEqual(f,1e-12,1e-12));
  f2->decrRef();
  MEDCouplingUMesh *m2=m->clone(true);
  CPPUNIT_ASSERT(m2->isEqual(m,1e-12));
  m2->setName("123");
  CPPUNIT_ASSERT(!m2->isEqual(m,1e-12));
  m2->copyTinyStringsFrom(m);
  CPPUNIT_ASSERT(m2->isEqual(m,1e-12));
  m2->getCoords()->setInfoOnComponent(1,"eee");
  CPPUNIT_ASSERT(!m2->isEqual(m,1e-12));
  m2->copyTinyStringsFrom(m);
  CPPUNIT_ASSERT(m2->isEqual(m,1e-12));
  m2->decrRef();
  //
  f->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testTryToShareSameCoordsPermute()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  MEDCouplingUMesh *m2=build3DSurfTargetMesh_1();
  CPPUNIT_ASSERT(m->getCoords()!=m2->getCoords());
  m->tryToShareSameCoordsPermute(*m2,1e-12);
  CPPUNIT_ASSERT(m->getCoords()==m2->getCoords());
  CPPUNIT_ASSERT(m2->isEqual(m,1e-12));
  const int renum1[9]={1,2,0,5,8,7,4,3,6};
  m->renumberNodes(renum1,9);
  CPPUNIT_ASSERT(m->getCoords()!=m2->getCoords());
  CPPUNIT_ASSERT(!m2->isEqual(m,1e-12));
  m->tryToShareSameCoordsPermute(*m2,1e-12);
  CPPUNIT_ASSERT(m->getCoords()==m2->getCoords());
  CPPUNIT_ASSERT(m2->isEqual(m,1e-12));
  m2->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testTryToShareSameCoordsPermute2()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_4();
  double targetCoords[8]={-0.3,-0.3, 0.2,-0.3, -0.3,0.2, 0.2,0.2 };
  int targetConn[4]={0,2,3,1};
  MEDCouplingUMesh *m2=MEDCouplingUMesh::New();
  m2->setMeshDimension(2);
  m2->allocateCells(1);
  m2->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
  m2->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(4,2);
  std::copy(targetCoords,targetCoords+8,myCoords->getPointer());
  m2->setCoords(myCoords);
  myCoords->decrRef();
  m2->checkCoherency();
  m1->checkCoherency();
  //
  const double expected1[5]={0.25,0.125,0.125,0.25,0.25};
  MEDCouplingFieldDouble *f1=m1->getMeasureField(false);
  MEDCouplingFieldDouble *f2=m2->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(5,f1->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f2->getArray()->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getIJ(i,0),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[0],f2->getIJ(0,0),1e-12);
  f2->decrRef();
  f1->decrRef();
  CPPUNIT_ASSERT_THROW(m1->tryToShareSameCoordsPermute(*m2,1e-12),INTERP_KERNEL::Exception);// <- here in this order the sharing is impossible.
  // Let's go for deeper test of tryToShareSameCoordsPermute
  m2->tryToShareSameCoordsPermute(*m1,1e-12);
  f1=m1->getMeasureField(false);
  f2=m2->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(5,f1->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f2->getArray()->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getIJ(i,0),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[0],f2->getIJ(0,0),1e-12);
  //
  f2->decrRef();
  f1->decrRef();
  //
  m1->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest::testChangeUnderlyingMesh1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_3();
  MEDCouplingUMesh *mesh2=build2DTargetMesh_3();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),2);
  const double arr[20]={7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.};
  std::copy(arr,arr+20,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //
  const int renum[]={0,2,1,3,4,5,6,8,7,9};
  mesh2->renumberCells(renum,false);
  CPPUNIT_ASSERT(f1->getMesh()==mesh1);
  f1->changeUnderlyingMesh(mesh1,10,1e-12);// nothing done only to check that nothing done.
  CPPUNIT_ASSERT(f1->getMesh()==mesh1);
  f1->changeUnderlyingMesh(mesh2,10,1e-12);
  CPPUNIT_ASSERT(f1->getMesh()==mesh2);
  const double expected1[20]={7.,107.,9.,109.,8.,108.,10.,110.,11.,111.,12.,112.,13.,113.,15.,115.,14.,114.,16.,116.};
  for(int i=0;i<20;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getArray()->getIJ(0,i),1e-12);
  f1->decrRef();
  //
  f1=MEDCouplingFieldDouble::New(ON_NODES,NO_TIME);
  f1->setMesh(mesh1);
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfNodes(),2);
  const double arr2[22]={7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.,17.,117.};
  std::copy(arr2,arr2+22,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //
  const int renum2[]={0,2,10,3,4,5,6,8,7,9,1};
  mesh2->renumberNodes(renum2,11);
  CPPUNIT_ASSERT(f1->getMesh()==mesh1);
  f1->changeUnderlyingMesh(mesh2,10,1e-12);
  CPPUNIT_ASSERT(f1->getMesh()==mesh2);
  const double expected2[22]={7.,107.,9.,109.,17.,117.,10.,110.,11.,111.,12.,112.,13.,113.,15.,115.,14.,114.,16.,116.,8.,108.};
  for(int i=0;i<22;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f1->getArray()->getIJ(0,i),1e-12);
  f1->decrRef();
  //
  mesh1->decrRef();
  mesh2->decrRef();
}

void MEDCouplingBasicsTest::testGetMaxValue1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  int nbOfCells=m->getNumberOfCells();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,LINEAR_TIME);
  f->setMesh(m);
  DataArrayDouble *a1=DataArrayDouble::New();
  a1->alloc(nbOfCells,1);
  const double val1[5]={3.,4.,5.,6.,7.};
  std::copy(val1,val1+5,a1->getPointer());
  DataArrayDouble *a2=DataArrayDouble::New();
  a2->alloc(nbOfCells,1);
  const double val2[5]={0.,1.,2.,8.,7.};
  std::copy(val2,val2+5,a2->getPointer());
  f->setArray(a1);
  f->setEndArray(a2);
  f->setEndTime(3.,3,4);
  f->checkCoherency();
  //
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.,f->getMaxValue(),1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,f->getMinValue(),1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.,f->getAverageValue(),1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(5.125,f->getWeightedAverageValue(),1e-14);
  a1->setIJ(0,2,9.5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9.5,f->getMaxValue(),1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,f->getMinValue(),1e-14);
  a2->setIJ(0,0,9.);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9.5,f->getMaxValue(),1e-14);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,f->getMinValue(),1e-14);
  //
  a2->decrRef();
  a1->decrRef();
  m->decrRef();
  f->decrRef();
}

void MEDCouplingBasicsTest::testSubstractInPlaceDM1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_3();
  MEDCouplingUMesh *mesh2=build2DTargetMesh_3();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),2);
  const double arr[20]={7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.};
  std::copy(arr,arr+20,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //
  CPPUNIT_ASSERT_EQUAL(10,f1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(20,f1->getNumberOfValues());
  //
  const int renum[]={0,2,1,3,4,5,6,8,7,9};
  mesh2->renumberCells(renum,false);
  //
  MEDCouplingFieldDouble *f2=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f2->setMesh(mesh2);
  array=DataArrayDouble::New();
  array->alloc(mesh2->getNumberOfCells(),2);
  const double arr2[20]={7.1,107.1,9.1,109.1,8.1,108.1,10.1,110.1,11.1,111.1,12.1,112.1,13.1,113.1,15.1,115.1,14.1,114.1,16.1,116.1};
  std::copy(arr2,arr2+20,array->getPointer());
  f2->setArray(array);
  array->decrRef();
  //
  f1->substractInPlaceDM(f2,10,1e-12);
  f1->applyFunc(1,"abs(x+y+0.2)");
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,f1->getMaxValue(),1e-14);
  //
  f1->decrRef();
  f2->decrRef();
  mesh1->decrRef();
  mesh2->decrRef();
}

void MEDCouplingBasicsTest::testDotCrossProduct1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_3();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setTime(2.3,5,6);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),3);
  const double arr1[30]={7.,107.,207.,8.,108.,208.,9.,109.,209.,10.,110.,210.,11.,111.,211.,12.,112.,212.,13.,113.,213.,14.,114.,214.,15.,115.,215.,16.,116.,216.};
  std::copy(arr1,arr1+30,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  MEDCouplingFieldDouble *f2=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f2->setTime(7.8,4,5);
  f2->setMesh(mesh1);
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),3);
  const double arr2[30]={1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.};
  std::copy(arr2,arr2+30,array->getPointer());
  f2->setArray(array);
  array->decrRef();
  //
  MEDCouplingFieldDouble *f3=f1->dot(*f2);
  const double expected1[10]={842.,1820.,2816.,3830.,4862.,5912.,6980.,8066.,9170.,10292.};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f3->getIJ(i,0),1e-9);
  f3->decrRef();
  //
  MEDCouplingFieldDouble *f4=f1->crossProduct(*f2);
  const double expected2[30]={-93., 186., -93., -392., 784., -392., -691., 1382., -691., -990., 1980., -990., -1289., 2578., -1289., -1588., 3176., -1588., -1887., 3774., -1887., -2186., 4372., -2186., -2485., 4970., -2485., -2784., 5568., -2784.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f4->getIJ(0,i),1e-9);
  f4->decrRef();
  //
  f2->decrRef();
  f1->decrRef();
  mesh1->decrRef();
}

void MEDCouplingBasicsTest::testMinMaxFields1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_3();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setTime(2.3,5,6);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),3);
  const double arr1[30]={7.,107.,207.,8.,108.,208.,9.,109.,209.,10.,110.,210.,11.,111.,211.,12.,112.,212.,13.,113.,213.,14.,114.,214.,15.,115.,215.,16.,116.,216.};
  std::copy(arr1,arr1+30,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  MEDCouplingFieldDouble *f2=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f2->setTime(7.8,4,5);
  f2->setMesh(mesh1);
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),3);
  const double arr2[30]={6.,108.,206.,9.,107.,209.,8.,110.,208.,11.,109.,211.,10.,112.,210.,13.,111.,213.,12.,114.,212.,15.,113.,215.,14.,116.,214.,17.,115.,217.};
  std::copy(arr2,arr2+30,array->getPointer());
  f2->setArray(array);
  array->decrRef();
  //
  MEDCouplingFieldDouble *f3=f1->max(*f2);
  const double expected1[30]={7.,108.,207.,9.,108.,209.,9.,110.,209.,11.,110.,211.,11.,112.,211.,13.,112.,213.,13.,114.,213.,15.,114.,215.,15.,116.,215.,17.,116.,217.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f3->getIJ(0,i),1e-9);
  f3->decrRef();
  //
  MEDCouplingFieldDouble *f4=f1->min(*f2);
  const double expected2[30]={6.,107.,206.,8.,107.,208.,8.,109.,208.,10.,109.,210.,10.,111.,210.,12.,111.,212.,12.,113.,212.,14.,113.,214.,14.,115.,214.,16.,115.,216.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f4->getIJ(0,i),1e-9);
  f4->decrRef();
  //
  f2->decrRef();
  f1->decrRef();
  mesh1->decrRef();
}

void MEDCouplingBasicsTest::testApplyLin1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_3();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,LINEAR_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),2);
  const double arr[20]={7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.,13.,113.,14.,114.,15.,115.,16.,116.};
  std::copy(arr,arr+20,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //
  f1->applyLin(2.,3.,0);
  const double expected1[20]={17.,107.,19.,108.,21.,109.,23.,110.,25.,111.,27.,112.,29.,113.,31.,114.,33.,115.,35.,116.};
  for(int i=0;i<20;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getIJ(0,i),1e-9);
  //
  const double arr2[20]={2.,102.,3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.};
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),2);
  std::copy(arr2,arr2+20,array->getPointer());
  f1->setEndArray(array);
  array->decrRef();
  //
  f1->applyLin(4.,5.,1);
  //
  const double expected2[20]={17.,433.,19.,437.,21.,441.,23.,445.,25.,449.,27.,453.,29.,457.,31.,461.,33.,465.,35.,469.};
  for(int i=0;i<20;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f1->getIJ(0,i),1e-9);
  const double expected3[20]={2.,413.,3.,417.,4.,421.,5.,425.,6.,429.,7.,433.,8.,437.,9.,441.,10.,445.,11.,449.};
  for(int i=0;i<20;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],f1->getEndArray()->getIJ(0,i),1e-9);
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testGetIdsInRange1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_3();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setTime(2.3,5,6);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),1);
  const double arr1[10]={2.,8.,6.,5.,11.,7.,9.,3.,10.,4.};
  std::copy(arr1,arr1+10,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //
  f1->checkCoherency();
  DataArrayInt *da=f1->getIdsInRange(2.9,7.1);
  CPPUNIT_ASSERT_EQUAL(5,da->getNbOfElems());
  const int expected1[5]={2,3,5,7,9};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+5,da->getConstPointer()));
  da->decrRef();
  da=f1->getIdsInRange(8.,12.);
  CPPUNIT_ASSERT_EQUAL(4,da->getNbOfElems());
  const int expected2[4]={1,4,6,8};
  CPPUNIT_ASSERT(std::equal(expected2,expected2+4,da->getConstPointer()));
  da->decrRef();
  //
  f1->decrRef();
  mesh1->decrRef();
}

void MEDCouplingBasicsTest::testBuildSubPart1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setTime(2.3,5,6);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),2);
  const double arr1[10]={3.,103.,4.,104.,5.,105.,6.,106.,7.,107.};
  std::copy(arr1,arr1+10,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //
  const int part1[3]={2,1,4};
  MEDCouplingFieldDouble *f2=f1->buildSubPart(part1,part1+3);
  CPPUNIT_ASSERT_EQUAL(3,f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f2->getNumberOfComponents());
  const double expected1[6]={5.,105.,4.,104.,7.,107.};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f2->getIJ(0,i),expected1[i],1e-12);
  CPPUNIT_ASSERT_EQUAL(3,f2->getMesh()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(6,f2->getMesh()->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getMeshDimension());
  MEDCouplingUMesh *m2C=dynamic_cast<MEDCouplingUMesh *>(const_cast<MEDCouplingMesh *>(f2->getMesh()));
  CPPUNIT_ASSERT_EQUAL(13,m2C->getMeshLength());
  const double expected2[12]={0.2, -0.3, 0.7, -0.3, 0.2, 0.2, 0.7, 0.2, 0.2, 0.7, 0.7, 0.7};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],m2C->getCoords()->getIJ(0,i),1.e-12);
  const double expected3[13]={3,2,3,1,3,0,2,1,4,4,5,3,2};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+13,m2C->getNodalConnectivity()->getConstPointer()));
  const double expected4[4]={0,4,8,13};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+4,m2C->getNodalConnectivityIndex()->getConstPointer()));
  f2->decrRef();
  f1->decrRef();
  // Test with field on nodes.
  f1=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  f1->setTime(2.3,5,6);
  f1->setMesh(mesh1);
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfNodes(),2);
  const double arr2[18]={3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.};
  std::copy(arr2,arr2+18,array->getPointer());  
  f1->setArray(array);
  array->decrRef();
  const int part2[4]={1,4,2,5};
  f2=f1->buildSubPart(part2,part2+4);
  CPPUNIT_ASSERT_EQUAL(4,f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f2->getNumberOfComponents());
  const double expected5[8]={4.,104.,5.,105.,7.,107.,8.,108.};
  for(int i=0;i<8;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f2->getIJ(0,i),expected5[i],1e-12);
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(4,f2->getMesh()->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getMeshDimension());
  m2C=dynamic_cast<MEDCouplingUMesh *>(const_cast<MEDCouplingMesh *>(f2->getMesh()));
  CPPUNIT_ASSERT_EQUAL(8,m2C->getMeshLength());
  for(int i=0;i<8;i++)//8 is not an error
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],m2C->getCoords()->getIJ(0,i),1.e-12);
  CPPUNIT_ASSERT(std::equal(expected3,expected3+4,m2C->getNodalConnectivity()->getConstPointer()+4));
  CPPUNIT_ASSERT(std::equal(expected3+4,expected3+8,m2C->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+3,m2C->getNodalConnectivityIndex()->getConstPointer()));
  f2->decrRef();
  //idem previous because nodes of cell#4 are not fully present in part3 
  const int part3[5]={1,4,2,5,7};
  DataArrayInt *arrr=DataArrayInt::New();
  arrr->alloc(5,1);
  std::copy(part3,part3+5,arrr->getPointer());
  f2=f1->buildSubPart(arrr);
  arrr->decrRef();
  CPPUNIT_ASSERT_EQUAL(4,f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f2->getNumberOfComponents());
  for(int i=0;i<8;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f2->getIJ(0,i),expected5[i],1e-12);
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(4,f2->getMesh()->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getMeshDimension());
  m2C=dynamic_cast<MEDCouplingUMesh *>(const_cast<MEDCouplingMesh *>(f2->getMesh()));
  CPPUNIT_ASSERT_EQUAL(8,m2C->getMeshLength());
  for(int i=0;i<8;i++)//8 is not an error
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],m2C->getCoords()->getIJ(0,i),1.e-12);
  CPPUNIT_ASSERT(std::equal(expected3,expected3+4,m2C->getNodalConnectivity()->getConstPointer()+4));
  CPPUNIT_ASSERT(std::equal(expected3+4,expected3+8,m2C->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+3,m2C->getNodalConnectivityIndex()->getConstPointer()));
  f2->decrRef();
  //
  const int part4[6]={1,4,2,5,7,8};
  f2=f1->buildSubPart(part4,part4+6);
  CPPUNIT_ASSERT_EQUAL(6,f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f2->getNumberOfComponents());
  const double expected6[12]={4.,104.,5.,105.,7.,107.,8.,108.,10.,110.,11.,111.};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(f2->getIJ(0,i),expected6[i],1e-12);
  CPPUNIT_ASSERT_EQUAL(3,f2->getMesh()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(6,f2->getMesh()->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getMeshDimension());
  m2C=dynamic_cast<MEDCouplingUMesh *>(const_cast<MEDCouplingMesh *>(f2->getMesh()));
  CPPUNIT_ASSERT_EQUAL(13,m2C->getMeshLength());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],m2C->getCoords()->getIJ(0,i),1.e-12);
  CPPUNIT_ASSERT(std::equal(expected3,expected3+4,m2C->getNodalConnectivity()->getConstPointer()+4));
  CPPUNIT_ASSERT(std::equal(expected3+4,expected3+8,m2C->getNodalConnectivity()->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected3+8,expected3+13,m2C->getNodalConnectivity()->getConstPointer()+8));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+4,m2C->getNodalConnectivityIndex()->getConstPointer()));
  f2->decrRef();
  //
  f1->decrRef();
  mesh1->decrRef();
}

void MEDCouplingBasicsTest::testDoublyContractedProduct1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),6);
  const double arr1[30]={7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5};
  std::copy(arr1,arr1+30,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  MEDCouplingFieldDouble *f2=f1->doublyContractedProduct();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3906.56,f2->getIJ(i,0),1e-9);
  f2->decrRef();
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testDeterminant1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,CONST_ON_TIME_INTERVAL);
  f1->setTime(2.3,5,6);
  f1->setEndTime(3.8,7,3);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),4);
  const double arr1[20]={1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5};
  std::copy(arr1,arr1+20,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //4 components
  f1->checkCoherency();
  MEDCouplingFieldDouble *f2=f1->determinant();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(CONST_ON_TIME_INTERVAL,f2->getTimeDiscretization());
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfValues());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(-2.42,f2->getIJ(i,0),1e-13);
  f2->decrRef();
  f1->decrRef();
  //6 components multi arrays with end array not defined
  f1=MEDCouplingFieldDouble::New(ON_NODES,LINEAR_TIME);
  f1->setTime(2.3,5,6);
  f1->setEndTime(3.8,7,3);
  f1->setMesh(mesh1);
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfNodes(),6);
  const double arr2[54]={1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7,
                         1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7};
  std::copy(arr2,arr2+54,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  CPPUNIT_ASSERT_THROW(f1->checkCoherency(),INTERP_KERNEL::Exception);//no end array specified !
  //
  f2=f1->determinant();
  CPPUNIT_ASSERT_EQUAL(LINEAR_TIME,f2->getTimeDiscretization());
  CPPUNIT_ASSERT_EQUAL(1,f2->getArray()->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f2->getNumberOfTuples());
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(137.335,f2->getIJ(i,0),1e-10);
  f2->decrRef();
  //6 components multi arrays with end array defined
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfNodes(),6);
  const double arr3[54]={7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5,
                         7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5};
  std::copy(arr3,arr3+54,array->getPointer());
  f1->setEndArray(array);
  array->decrRef();
  f1->checkCoherency();
  f2=f1->determinant();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(LINEAR_TIME,f2->getTimeDiscretization());
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f2->getNumberOfTuples());
  int it,order;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.3,f2->getTime(it,order),1e-12);
  CPPUNIT_ASSERT_EQUAL(5,it); CPPUNIT_ASSERT_EQUAL(6,order);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.8,f2->getEndTime(it,order),1e-12);
  CPPUNIT_ASSERT_EQUAL(7,it); CPPUNIT_ASSERT_EQUAL(3,order);
  for(int i=0;i<9;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(137.335,f2->getIJ(i,0),1e-10);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1289.685,f2->getEndArray()->getIJ(i,0),1e-9);
    }
  f2->decrRef();
  f1->decrRef();
  //9 components
  f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setTime(7.8,10,2);
  f1->setMesh(mesh1);
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),9);
  const double arr4[45]={1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1};
  std::copy(arr4,arr4+45,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //
  f1->checkCoherency();
  f2=f1->determinant();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(ONE_TIME,f2->getTimeDiscretization());
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.8,f2->getTime(it,order),1e-12);
  CPPUNIT_ASSERT_EQUAL(10,it); CPPUNIT_ASSERT_EQUAL(2,order);
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(3.267,f2->getIJ(i,0),1e-13);
  f2->decrRef();
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testEigenValues1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),6);
  const double arr1[30]={1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7};
  std::copy(arr1,arr1+30,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  MEDCouplingFieldDouble *f2=f1->eigenValues();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(3,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  const double expected1[3]={13.638813677891717,-4.502313844635971,-2.2364998332557486};
  for(int i=0;i<5;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[0],f2->getIJ(i,0),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[1],f2->getIJ(i,1),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[2],f2->getIJ(i,2),1e-13);
    }
  f2->decrRef();
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testEigenVectors1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),6);
  const double arr1[30]={1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7};
  std::copy(arr1,arr1+30,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  MEDCouplingFieldDouble *f2=f1->eigenVectors();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(9,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  const double expected1[9]={
    0.5424262364180696, 0.5351201064614425, 0.6476266283176001,//eigenvect 0
    0.7381111277307373, 0.06458838384003074, -0.6715804522117897,//eigenvect 1
    -0.4012053603397987, 0.8423032781211455, -0.3599436712889738//eigenvect 2
  };
  for(int i=0;i<5;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[0],f2->getIJ(i,0),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[1],f2->getIJ(i,1),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[2],f2->getIJ(i,2),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[3],f2->getIJ(i,3),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[4],f2->getIJ(i,4),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[5],f2->getIJ(i,5),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[6],f2->getIJ(i,6),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[7],f2->getIJ(i,7),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[8],f2->getIJ(i,8),1e-13);
    }
  f2->decrRef();
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testInverse1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),9);
  const double arr1[45]={1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1};
  std::copy(arr1,arr1+45,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  MEDCouplingFieldDouble *f2=f1->inverse();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(9,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  const double expected1[9]={-2.6538108356290113, 2.855831037649208, -1.1111111111111067, 3.461891643709813, -4.775022956841121, 2.2222222222222143, -1.1111111111111054, 2.222222222222214, -1.1111111111111072};
  for(int i=0;i<5;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[0],f2->getIJ(i,0),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[1],f2->getIJ(i,1),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[2],f2->getIJ(i,2),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[3],f2->getIJ(i,3),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[4],f2->getIJ(i,4),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[5],f2->getIJ(i,5),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[6],f2->getIJ(i,6),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[7],f2->getIJ(i,7),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[8],f2->getIJ(i,8),1e-13);
    }
  f2->decrRef();
  //
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),6);
  const double arr3[30]={7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5};
  std::copy(arr3,arr3+30,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  f2=f1->inverse();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(6,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  const double expected3[6]={-0.3617705098531818, -0.8678630828458127, -0.026843764174972983, 0.5539957431465833, 0.13133439560823013, -0.05301294502145887};
  for(int i=0;i<5;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[0],f2->getIJ(i,0),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[1],f2->getIJ(i,1),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[2],f2->getIJ(i,2),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[3],f2->getIJ(i,3),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[4],f2->getIJ(i,4),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[5],f2->getIJ(i,5),1e-13);
    }
  f2->decrRef();
  //
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),4);
  const double arr2[20]={1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5};
  std::copy(arr2,arr2+20,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  f2=f1->inverse();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(4,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  const double expected2[4]={-1.8595041322314059, 0.9504132231404963, 1.404958677685951, -0.49586776859504156};
  for(int i=0;i<5;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[0],f2->getIJ(i,0),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[1],f2->getIJ(i,1),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[2],f2->getIJ(i,2),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[3],f2->getIJ(i,3),1e-13);
    }
  f2->decrRef();
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testTrace1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),9);
  const double arr1[45]={1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1, 1.2,2.3,3.4,4.5,5.6,6.7,7.8,8.9,9.1};
  std::copy(arr1,arr1+45,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  MEDCouplingFieldDouble *f2=f1->trace();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(15.9,f2->getIJ(i,0),1e-13);
  f2->decrRef();
  //
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),6);
  const double arr3[30]={7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5, 7.8,8.9,9.1,10.2,23.4,34.5};
  std::copy(arr3,arr3+30,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  f2=f1->trace();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(25.8,f2->getIJ(i,0),1e-13);
  f2->decrRef();
  //
  array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),4);
  const double arr2[20]={1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5, 1.2,2.3,3.4,4.5};
  std::copy(arr2,arr2+20,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  f2=f1->trace();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.7,f2->getIJ(i,0),1e-13);
  f2->decrRef();
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testDeviator1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),6);
  const double arr1[30]={1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7, 1.2,2.3,3.4,4.5,5.6,6.7};
  std::copy(arr1,arr1+30,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  MEDCouplingFieldDouble *f2=f1->deviator();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(6,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  const double expected1[6]={-1.1,0.,1.1,4.5,5.6,6.7};
  for(int i=0;i<5;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[0],f2->getIJ(i,0),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[1],f2->getIJ(i,1),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[2],f2->getIJ(i,2),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[3],f2->getIJ(i,3),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[4],f2->getIJ(i,4),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[5],f2->getIJ(i,5),1e-13);
    }
  f2->decrRef();
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testMagnitude1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),5);
  const double arr1[25]={1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6, 1.2,2.3,3.4,4.5,5.6};
  std::copy(arr1,arr1+25,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  MEDCouplingFieldDouble *f2=f1->magnitude();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(8.3606219864313918,f2->getIJ(i,0),1e-13);
  f2->decrRef();
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testMaxPerTuple1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),5);
  const double arr1[25]={1.2,2.3,3.4,4.5,5.6, 1.2,3.4,4.5,5.6,2.3, 3.4,4.5,5.6,1.2,2.3, 5.6,1.2,2.3,3.4,4.5, 4.5,5.6,1.2,2.3,3.4};
  std::copy(arr1,arr1+25,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  MEDCouplingFieldDouble *f2=f1->maxPerTuple();
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(5.6,f2->getIJ(i,0),1e-13);
  f2->decrRef();
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testChangeNbOfComponents()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),5);
  const double arr1[25]={1.2,2.3,3.4,4.5,5.6, 1.2,3.4,4.5,5.6,2.3, 3.4,4.5,5.6,1.2,2.3, 5.6,1.2,2.3,3.4,4.5, 4.5,5.6,1.2,2.3,3.4};
  std::copy(arr1,arr1+25,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  f1->changeNbOfComponents(3,7.77);
  f1->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(3,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  const double expected1[15]={1.2,2.3,3.4, 1.2,3.4,4.5, 3.4,4.5,5.6, 5.6,1.2,2.3, 4.5,5.6,1.2};
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getIJ(0,i),1e-13);
  f1->changeNbOfComponents(4,7.77);
  f1->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(4,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  const double expected2[20]={1.2,2.3,3.4,7.77, 1.2,3.4,4.5,7.77, 3.4,4.5,5.6,7.77, 5.6,1.2,2.3,7.77, 4.5,5.6,1.2,7.77};
  for(int i=0;i<20;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f1->getIJ(0,i),1e-13);
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testSortPerTuple1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  f1->setMesh(mesh1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),5);
  const double arr1[25]={1.2,2.3,3.4,4.5,5.6, 1.2,3.4,4.5,5.6,2.3, 3.4,4.5,5.6,1.2,2.3, 5.6,1.2,2.3,3.4,4.5, 4.5,5.6,1.2,2.3,3.4};
  std::copy(arr1,arr1+25,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  f1->checkCoherency();
  //
  f1->sortPerTuple(true);
  f1->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  for(int i=0;i<5;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[0],f1->getIJ(i,0),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[1],f1->getIJ(i,1),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[2],f1->getIJ(i,2),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[3],f1->getIJ(i,3),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[4],f1->getIJ(i,4),1e-13);
    }
  //
  f1->sortPerTuple(false);
  f1->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  for(int i=0;i<5;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[4],f1->getIJ(i,0),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[3],f1->getIJ(i,1),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[2],f1->getIJ(i,2),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[1],f1->getIJ(i,3),1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(arr1[0],f1->getIJ(i,4),1e-13);
    }
  //
  mesh1->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testIsEqualWithoutConsideringStr1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  MEDCouplingUMesh *mesh2=build2DTargetMesh_1();
  DataArrayInt *da1,*da2;
  //
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  mesh2->setName("rr");
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  mesh1->checkDeepEquivalWith(mesh2,2,1e-12,da1,da2);
  CPPUNIT_ASSERT_THROW(mesh1->checkGeoEquivalWith(mesh2,0,1e-12,da1,da2),INTERP_KERNEL::Exception);
  mesh2->setName("");
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  mesh2->getCoords()->setInfoOnComponent(0,"tty");
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  mesh2->getCoords()->setInfoOnComponent(0,"");
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  mesh2->getCoords()->setInfoOnComponent(1,"tty");
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  mesh2->getCoords()->setInfoOnComponent(1,"");
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  double tmp=mesh2->getCoords()->getIJ(0,3);
  mesh2->getCoords()->setIJ(0,3,9999.);
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(!mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  mesh2->getCoords()->setIJ(0,3,tmp);
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  int tmp2=mesh2->getNodalConnectivity()->getIJ(0,4);
  mesh2->getNodalConnectivity()->setIJ(0,4,0);
  CPPUNIT_ASSERT(!mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(!mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  mesh2->getNodalConnectivity()->setIJ(0,4,tmp2);
  CPPUNIT_ASSERT(mesh1->isEqual(mesh2,1e-12));
  CPPUNIT_ASSERT(mesh1->isEqualWithoutConsideringStr(mesh2,1e-12));
  //
  MEDCouplingFieldDouble *f1=mesh1->getMeasureField(true);
  MEDCouplingFieldDouble *f2=mesh2->getMeasureField(true);
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  CPPUNIT_ASSERT(f1->isEqualWithoutConsideringStr(f2,1e-12,1e-12));
  f2->setName("ftest");
  CPPUNIT_ASSERT(!f1->isEqual(f2,1e-12,1e-12));
  CPPUNIT_ASSERT(f1->isEqualWithoutConsideringStr(f2,1e-12,1e-12));
  f1->setName("ftest");
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  CPPUNIT_ASSERT(f1->isEqualWithoutConsideringStr(f2,1e-12,1e-12));
  //
  f2->getArray()->setInfoOnComponent(0,"eee");
  CPPUNIT_ASSERT(!f1->isEqual(f2,1e-12,1e-12));
  CPPUNIT_ASSERT(f1->isEqualWithoutConsideringStr(f2,1e-12,1e-12));
  f2->getArray()->setInfoOnComponent(0,"");
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  CPPUNIT_ASSERT(f1->isEqualWithoutConsideringStr(f2,1e-12,1e-12));
  //
  f2->getArray()->setIJ(1,0,0.123);
  CPPUNIT_ASSERT(!f1->isEqual(f2,1e-12,1e-12));
  CPPUNIT_ASSERT(!f1->isEqualWithoutConsideringStr(f2,1e-12,1e-12));
  f2->getArray()->setIJ(1,0,0.125);
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  CPPUNIT_ASSERT(f1->isEqualWithoutConsideringStr(f2,1e-12,1e-12));
  //
  f1->decrRef();
  f2->decrRef();
  //
  mesh1->decrRef();
  mesh2->decrRef();
}

void MEDCouplingBasicsTest::testGetNodeIdsOfCell1()
{
  MEDCouplingUMesh *mesh1=build2DTargetMesh_1();
  std::vector<int> nodeIds;
  mesh1->getNodeIdsOfCell(1,nodeIds);
  CPPUNIT_ASSERT_EQUAL(3,(int)nodeIds.size());
  CPPUNIT_ASSERT_EQUAL(1,nodeIds[0]);
  CPPUNIT_ASSERT_EQUAL(4,nodeIds[1]);
  CPPUNIT_ASSERT_EQUAL(2,nodeIds[2]);
  std::vector<double> coords;
  mesh1->getCoordinatesOfNode(4,coords);
  CPPUNIT_ASSERT_EQUAL(2,(int)coords.size());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2,coords[0],1e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.2,coords[1],1e-13);
  mesh1->decrRef();
}

void MEDCouplingBasicsTest::testGetEdgeRatioField1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=m1->getEdgeRatioField();
  CPPUNIT_ASSERT_EQUAL(m1->getNumberOfCells(),f1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  const double expected1[5]={1.,1.4142135623730951, 1.4142135623730951,1.,1.};
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getIJ(i,0),1e-14);
  f1->decrRef();
  m1->decrRef();
  //
  m1=build3DSurfTargetMesh_1();
  f1=m1->getEdgeRatioField();
  CPPUNIT_ASSERT_EQUAL(m1->getNumberOfCells(),f1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  const double expected2[5]={1.4142135623730951, 1.7320508075688772, 1.7320508075688772, 1.4142135623730951, 1.4142135623730951};
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f1->getIJ(i,0),1e-14);
  f1->decrRef();
  m1->decrRef();
}

void MEDCouplingBasicsTest::testFillFromAnalytic3()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  CPPUNIT_ASSERT_THROW(f1->fillFromAnalytic(1,"y+x"),INTERP_KERNEL::Exception);
  f1->setMesh(m);
  f1->setName("myField");
  f1->fillFromAnalytic(1,"y+x");
  f1->checkCoherency();
  CPPUNIT_ASSERT(std::string(f1->getName())=="myField");
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_CELLS);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==ONE_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  double values1[5]={-0.1,0.23333333333333336,0.56666666666666665,0.4,0.9};
  const double *tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+5,values1,values1,std::minus<double>());
  std::transform(values1,values1+5,values1,std::ptr_fun<double,double>(fabs));
  double max=*std::max_element(values1,values1+5);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=MEDCouplingFieldDouble::New(ON_NODES,CONST_ON_TIME_INTERVAL);
  f1->setMesh(m);
  f1->setEndTime(1.2,3,4);
  f1->fillFromAnalytic(1,"y+2*x");
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==CONST_ON_TIME_INTERVAL);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  double values2[9]={-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+9,values2,values2,std::minus<double>());
  std::transform(values2,values2+9,values2,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2,values2+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  f1=MEDCouplingFieldDouble::New(ON_NODES,LINEAR_TIME);
  f1->setMesh(m);
  f1->setEndTime(1.2,3,4);
  f1->fillFromAnalytic(1,"2.*x+y");
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==LINEAR_TIME);
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  tmp=f1->getArray()->getConstPointer();
  double values2Bis[9]={-0.9,0.1,1.1,-0.4,0.6,1.6,0.1,1.1,2.1};
  double values2BisBis[9];
  std::transform(tmp,tmp+9,values2Bis,values2BisBis,std::minus<double>());
  std::transform(values2,values2+9,values2BisBis,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2BisBis,values2BisBis+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  tmp=f1->getEndArray()->getConstPointer();
  std::transform(tmp,tmp+9,values2Bis,values2BisBis,std::minus<double>());
  std::transform(values2,values2+9,values2BisBis,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values2BisBis,values2BisBis+9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  f1->decrRef();
  //
  f1=MEDCouplingFieldDouble::New(ON_NODES,NO_TIME);
  f1->setMesh(m);
  f1->fillFromAnalytic(2,"(x+y)*IVec+2*(x+y)*JVec");
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->getTypeOfField()==ON_NODES);
  CPPUNIT_ASSERT(f1->getTimeDiscretization()==NO_TIME);
  CPPUNIT_ASSERT_EQUAL(2,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  double values3[18]={-0.6,-1.2,-0.1,-0.2,0.4,0.8,-0.1,-0.2,0.4,0.8,0.9,1.8,0.4,0.8,0.9,1.8,1.4,2.8};
  tmp=f1->getArray()->getConstPointer();
  std::transform(tmp,tmp+18,values3,values3,std::minus<double>());
  std::transform(values3,values3+18,values3,std::ptr_fun<double,double>(fabs));
  max=*std::max_element(values3,values3+18);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,max,1.e-12);
  double values4[2];
  f1->accumulate(values4);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3.6,values4[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(7.2,values4[1],1.e-12);
  f1->integral(true,values4);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.5,values4[0],1.e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(1.,values4[1],1.e-12);
  f1->decrRef();
  //
  f1=MEDCouplingFieldDouble::New(ON_NODES,NO_TIME);
  f1->setMesh(m);
  CPPUNIT_ASSERT_THROW(f1->fillFromAnalytic(1,"1./(x-0.2)"),INTERP_KERNEL::Exception);
  //
  m->decrRef();
  f1->decrRef();
}

void MEDCouplingBasicsTest::testFieldDoubleOpEqual1()
{
  MEDCouplingUMesh *m=build2DTargetMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  CPPUNIT_ASSERT_THROW((*f1)=0.07,INTERP_KERNEL::Exception);
  f1->setMesh(m);
  (*f1)=0.07;
  f1->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.07,f1->getIJ(i,0),1e-16);
  (*f1)=0.09;
  f1->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f1->getNumberOfTuples());
  for(int i=0;i<5;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.09,f1->getIJ(i,0),1e-16);
  f1->decrRef();
  //
  f1=MEDCouplingFieldDouble::New(ON_NODES,LINEAR_TIME);
  f1->setEndTime(4.5,2,3);
  f1->setMesh(m);
  (*f1)=0.08;
  f1->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(1,f1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getNumberOfTuples());
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08,f1->getIJ(i,0),1e-16);
  CPPUNIT_ASSERT_EQUAL(1,f1->getEndArray()->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(9,f1->getEndArray()->getNumberOfTuples());
  for(int i=0;i<9;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(0.08,f1->getEndArray()->getIJ(i,0),1e-16);
  f1->decrRef();
  //
  m->decrRef();
}

void MEDCouplingBasicsTest::testAreaBary3D2()
{
  const double coordsForHexa8[24]={
    -75.45749305371, 180.95495078401, 39.515472018008,
    -9.755591679144, 23.394927935279, 5.108794294848,
    14.337630157832, 61.705351002702, 160.42422501908,
    -27.273893776752, 167.567731083961, 192.830034145464,
    //
    99.857193154796,264.499264735586,-8.287335493412,
    144.939882761126,156.38626563134,-31.896173894226,
    161.34096835726,182.4654895809,73.832387065572,
    132.680430393685,255.37973247196,96.15235602819
  };
  const double volHexa8=3258520.29637466;
  const double baryHexa8[3]={43.925705821778, 155.31893955289, 65.874418109644};

  const double coordsForPenta6[18]={
    -68.199829618726,178.938498373416,62.608505919588,
    8.461744647847,76.653979804423,165.00018874933,
    -27.273893776752,167.567731083961,192.830034145464,
    //
    106.586501038965,262.629609408327,13.124533008813,
    155.465082847275,197.414118382622,78.408350795821,
    132.680430393685,255.37973247196,96.15235602819
  };
  const double volPenta6=944849.868507338;
  const double baryPenta6[3]={39.631002313543,182.692711783428,106.98540473964};
  
  const double coordsForPyra5[15]={
    132.680430393685,255.37973247196,96.15235602819,
    -27.273893776752,167.567731083961,192.830034145464,
    8.461744647847,76.653979804423,165.00018874933,
    155.465082847275,197.414118382622,78.408350795821,
    //
    -68.199829618726,178.938498373416,62.608505919588
  };
  const double volPyra5=756943.92980254;
  const double baryPyra5[3]={29.204294116618,172.540129749156,118.01035951483};
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New("Bary3D2",3);
  DataArrayDouble *coo=DataArrayDouble::New();
  coo->alloc(19,3);
  double *tmp=std::copy(coordsForHexa8,coordsForHexa8+24,coo->getPointer());
  tmp=std::copy(coordsForPenta6,coordsForPenta6+18,tmp);
  std::copy(coordsForPyra5,coordsForPyra5+15,tmp);
  mesh->setCoords(coo);
  coo->decrRef();
  //
  int tmpConn[8]={0,1,2,3,4,5,6,7};
  mesh->allocateCells(3);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,tmpConn);
  std::transform(tmpConn,tmpConn+8,tmpConn,std::bind2nd(std::plus<int>(),8));
  mesh->insertNextCell(INTERP_KERNEL::NORM_PENTA6,6,tmpConn);
  std::transform(tmpConn,tmpConn+8,tmpConn,std::bind2nd(std::plus<int>(),6));
  mesh->insertNextCell(INTERP_KERNEL::NORM_PYRA5,5,tmpConn);
  mesh->finishInsertingCells();
  mesh->checkCoherency();
  bool isMerged;
  int newNebOfNodes;
  DataArrayInt *da=mesh->mergeNodes(1e-7,isMerged,newNebOfNodes);
  da->decrRef();
  CPPUNIT_ASSERT_EQUAL(12,newNebOfNodes);
  MEDCouplingFieldDouble *vols=mesh->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(3,vols->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,vols->getNumberOfComponents());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(volHexa8,vols->getIJ(0,0),1e-6);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(volPenta6,vols->getIJ(1,0),1e-7);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(volPyra5,vols->getIJ(2,0),1e-7);
  vols->decrRef();
  DataArrayDouble *bary=mesh->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(3,bary->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,bary->getNumberOfComponents());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(baryHexa8[0],bary->getIJ(0,0),1e-11);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(baryHexa8[1],bary->getIJ(0,1),1e-11);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(baryHexa8[2],bary->getIJ(0,2),1e-11);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(baryPenta6[0],bary->getIJ(1,0),1e-11);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(baryPenta6[1],bary->getIJ(1,1),1e-11);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(baryPenta6[2],bary->getIJ(1,2),1e-11);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(baryPyra5[0],bary->getIJ(2,0),1e-11);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(baryPyra5[1],bary->getIJ(2,1),1e-11);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(baryPyra5[2],bary->getIJ(2,2),1e-11);
  bary->decrRef();
  //
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testGetMeasureFieldCMesh1()
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
  m->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(4,m->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(1,m->getSpaceDimension());
  MEDCouplingFieldDouble *f=m->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(3,f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f->getNumberOfComponents());
  const double expected1[3]={1.1,2.4,4.4};
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f->getIJ(i,0),1e-12);
  f->decrRef();
  DataArrayDouble *coords=m->getCoordinatesAndOwner();
  CPPUNIT_ASSERT_EQUAL(4,coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,coords->getNumberOfComponents());
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(discX[i],coords->getIJ(i,0),1e-12);
  coords->decrRef();
  coords=m->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(3,coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,coords->getNumberOfComponents());
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
  m->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(12,m->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(6,m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(2,m->getSpaceDimension());
  f=m->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(6,f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f->getNumberOfComponents());
  const double expected2[6]={12.21,26.64,48.84,24.64,53.76,98.56};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f->getIJ(i,0),1e-12);
  f->decrRef();
  coords=m->getCoordinatesAndOwner();
  CPPUNIT_ASSERT_EQUAL(12,coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,coords->getNumberOfComponents());
  const double expected2_2[24]={2.3,12.3,3.4,12.3,5.8,12.3,10.2,12.3, 2.3,23.4,3.4,23.4,5.8,23.4,10.2,23.4, 2.3,45.8,3.4,45.8,5.8,45.8,10.2,45.8};
  for(int i=0;i<24;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2_2[i],coords->getIJ(0,i),1e-12);
  coords->decrRef();
  coords=m->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(6,coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,coords->getNumberOfComponents());
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
  m->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(60,m->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(24,m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(3,m->getSpaceDimension());
  f=m->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(24,f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,f->getNumberOfComponents());
  const double expected3[24]={23.199, 50.616, 92.796, 46.816, 102.144, 187.264, 0.6105, 1.332, 2.442, 1.232, 2.688, 4.928, 10.7448, 23.4432, 42.9792, 21.6832, 47.3088, 86.7328, 6.5934, 14.3856, 26.3736, 13.3056, 29.0304, 53.2224};
  for(int i=0;i<24;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],f->getIJ(i,0),1e-12);
  f->decrRef();
  coords=m->getCoordinatesAndOwner();
  CPPUNIT_ASSERT_EQUAL(60,coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,coords->getNumberOfComponents());
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
  coords=m->getBarycenterAndOwner();
  CPPUNIT_ASSERT_EQUAL(24,coords->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,coords->getNumberOfComponents());
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

void MEDCouplingBasicsTest::testFieldDoubleZipCoords1()
{
  MEDCouplingUMesh *m=build2DTargetMeshMergeNode_1();
  MEDCouplingFieldDouble *f=m->fillFromAnalytic(ON_NODES,2,"x*2.");
  f->getArray()->setInfoOnComponent(0,"titi");
  f->getArray()->setInfoOnComponent(1,"tutu");
  f->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(18,f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f->getNumberOfComponents());
  const double expected1[36]={-0.6, -0.6, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 0.4, 0.4};
  for(int i=0;i<36;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(f->zipCoords());
  f->checkCoherency();
  const double expected2[30]={-0.6, -0.6, 1.4, 1.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 1.4, 1.4, -0.6, -0.6, 0.4, 0.4, 1.4, 1.4, 0.4, 0.4};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(!f->zipCoords());
  f->checkCoherency();
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(std::string(f->getArray()->getInfoOnComponent(0))=="titi");
  CPPUNIT_ASSERT(std::string(f->getArray()->getInfoOnComponent(1))=="tutu");
  f->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testFieldDoubleZipConnectivity1()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  const int cells1[3]={2,3,4};
  MEDCouplingPointSet *m3_1=m2->buildPartOfMySelf(cells1,cells1+3,true);
  MEDCouplingUMesh *m3=dynamic_cast<MEDCouplingUMesh *>(m3_1);
  CPPUNIT_ASSERT(m3);
  m2->decrRef();
  MEDCouplingUMesh *m4=build2DSourceMesh_1();
  MEDCouplingUMesh *m5=MEDCouplingUMesh::mergeUMeshes(m1,m3);
  m1->decrRef();
  m3->decrRef();
  MEDCouplingUMesh *m6=MEDCouplingUMesh::mergeUMeshes(m5,m4);
  m4->decrRef();
  m5->decrRef();
  //
  CPPUNIT_ASSERT_EQUAL(10,m6->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(22,m6->getNumberOfNodes());
  bool areNodesMerged;
  int newNbOfNodes;
  DataArrayInt *arr=m6->mergeNodes(1e-13,areNodesMerged,newNbOfNodes);
  CPPUNIT_ASSERT_EQUAL(9,m6->getNumberOfNodes());
  arr->decrRef();
  MEDCouplingFieldDouble *f=m6->fillFromAnalytic(ON_CELLS,2,"x");
  MEDCouplingFieldDouble *f2=m6->fillFromAnalytic(ON_NODES,2,"x");
  CPPUNIT_ASSERT_EQUAL(10,f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f->getNumberOfComponents());
  const double expected1[20]={-0.05, -0.05, 0.3666666666666667, 0.3666666666666667, 0.53333333333333321, 0.53333333333333321,
                              -0.05, -0.05, 0.45, 0.45, 0.53333333333333321, 0.53333333333333321, -0.05, -0.05, 0.45, 0.45,
                              0.36666666666666659, 0.36666666666666659, 0.033333333333333326, 0.033333333333333326};
  for(int i=0;i<20;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f->getIJ(0,i),1e-12);
  f->getArray()->setInfoOnComponent(0,"titi");
  f->getArray()->setInfoOnComponent(1,"tutu");
  f->checkCoherency();
  CPPUNIT_ASSERT(f->zipConnectivity(0));
  const double expected2[14]={-0.05, -0.05, 0.3666666666666667, 0.3666666666666667, 0.53333333333333321, 0.53333333333333321,
                              -0.05, -0.05, 0.45, 0.45, 0.36666666666666659, 0.36666666666666659, 0.033333333333333326, 0.033333333333333326};
  CPPUNIT_ASSERT_EQUAL(7,f->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f->getNumberOfComponents());
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],f->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(std::string(f->getArray()->getInfoOnComponent(0))=="titi");
  CPPUNIT_ASSERT(std::string(f->getArray()->getInfoOnComponent(1))=="tutu");
  CPPUNIT_ASSERT(!f->zipConnectivity(0));
  f->decrRef();
  //
  const double expected3[18]={-0.3, -0.3, 0.2, 0.2, 0.7, 0.7, -0.3, -0.3, 0.2, 0.2, 0.7, 0.7, 
                              -0.3, -0.3, 0.2, 0.2, 0.7, 0.7};
  CPPUNIT_ASSERT_EQUAL(9,f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f2->getNumberOfComponents());
  for(int i=0;i<18;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],f2->getIJ(0,i),1e-12);
  CPPUNIT_ASSERT(f2->zipConnectivity(0));
  CPPUNIT_ASSERT_EQUAL(9,f2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,f2->getNumberOfComponents());
  for(int i=0;i<18;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],f2->getIJ(0,i),1e-12);
  f2->decrRef();
  //
  m6->decrRef();
}

void MEDCouplingBasicsTest::testDaDoubleRenumber1()
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
  CPPUNIT_ASSERT_EQUAL(7,b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,b->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(7,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(1))=="tata");
  const int expected2[14]={3, 13, 2, 12, 7, 17, 1, 11, 6, 16, 5, 15, 4, 14};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d->getIJ(0,i));
  c->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest::testDaDoubleRenumberAndReduce1()
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
  CPPUNIT_ASSERT_EQUAL(5,b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,b->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(5,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(1))=="tata");
  const int expected2[10]={5,15,3,13,1,11,7,17,6,16};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d->getIJ(0,i));
  c->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest::testDaDoubleRenumberInPlace1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(7,2);
  const double arr1[14]={1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1};
  std::copy(arr1,arr1+14,a->getPointer());
  //
  const int arr2[7]={3,1,0,6,5,4,2};
  a->renumberInPlace(arr2);
  CPPUNIT_ASSERT_EQUAL(7,a->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,a->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(7,c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,c->getNumberOfComponents());
  const int expected2[14]={3, 13, 2, 12, 7, 17, 1, 11, 6, 16, 5, 15, 4, 14};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],c->getIJ(0,i));
  c->decrRef();
}

void MEDCouplingBasicsTest::testDaDoubleRenumberR1()
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
  CPPUNIT_ASSERT_EQUAL(7,b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,b->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(7,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(1))=="tata");
  const int expected2[14]={4, 14, 2, 12, 1, 11, 7, 17, 6, 16, 5, 15, 3, 13};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d->getIJ(0,i));
  c->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest::testDaDoubleRenumberInPlaceR1()
{
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(7,2);
  const double arr1[14]={1.1,11.1,2.1,12.1,3.1,13.1,4.1,14.1,5.1,15.1,6.1,16.1,7.1,17.1};
  std::copy(arr1,arr1+14,a->getPointer());
  //
  const int arr2[7]={3,1,0,6,5,4,2};
  a->renumberInPlaceR(arr2);
  CPPUNIT_ASSERT_EQUAL(7,a->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,a->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(7,c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,c->getNumberOfComponents());
  const int expected2[14]={4, 14, 2, 12, 1, 11, 7, 17, 6, 16, 5, 15, 3, 13};
  for(int i=0;i<14;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],c->getIJ(0,i));
  c->decrRef();
}

void MEDCouplingBasicsTest::testDaDoubleSelectByTupleId1()
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
  CPPUNIT_ASSERT_EQUAL(5,b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,b->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(5,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,d->getNumberOfComponents());
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(0))=="toto");
  CPPUNIT_ASSERT(std::string(d->getInfoOnComponent(1))=="tata");
  const int expected2[10]={5,15,3,13,1,11,7,17,6,16};
  for(int i=0;i<10;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],d->getIJ(0,i));
  c->decrRef();
  d->decrRef();
}

void MEDCouplingBasicsTest::testDaDoubleGetMinMaxValues1()
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
  CPPUNIT_ASSERT_EQUAL(3,ws->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,ws->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(3,ws->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,ws->getNumberOfComponents());
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],ws->getIJ(i,0));
  ws->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest::testFieldDoubleGetMinMaxValues2()
{
  MEDCouplingUMesh *m1=0;
  MEDCouplingUMesh *m2=build3DExtrudedUMesh_1(m1);
  m1->decrRef();
  CPPUNIT_ASSERT_EQUAL(18,m2->getNumberOfCells());
  const double arr1[18]={8.71,4.53,-12.41,8.71,-8.71,8.7099,4.55,8.71,5.55,6.77,-1e-200,4.55,8.7099,0.,1.23,0.,2.22,8.71};
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  DataArrayDouble *a=DataArrayDouble::New();
  a->alloc(18,1);
  std::copy(arr1,arr1+18,a->getPointer());
  f->setArray(a);
  a->decrRef();
  f->setMesh(m2);
  //
  f->checkCoherency();
  double m=f->getMaxValue();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.71,m,1e-12);
  DataArrayInt *ws;
  m=f->getMaxValue2(ws);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(8.71,m,1e-12);
  CPPUNIT_ASSERT_EQUAL(4,ws->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,ws->getNumberOfComponents());
  const int expected1[4]={0,3,7,17};
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],ws->getIJ(i,0));
  ws->decrRef();
  //
  const double arr2[18]={-8.71,-4.53,12.41,-8.71,8.71,-8.7099,-4.55,-8.71,-5.55,-6.77,1e-200,-4.55,-8.7099,0.,-1.23,0.,-2.22,-8.71};
  std::copy(arr2,arr2+18,a->getPointer());
  f->checkCoherency();
  m=f->getMinValue();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-8.71,m,1e-12);
  m=f->getMinValue2(ws);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-8.71,m,1e-12);
  CPPUNIT_ASSERT_EQUAL(4,ws->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,ws->getNumberOfComponents());
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],ws->getIJ(i,0));
  ws->decrRef();
  //
  f->decrRef();
  m2->decrRef();
}

void MEDCouplingBasicsTest::testBuildUnstructuredCMesh1()
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
  m->checkCoherency();
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
  m2->checkCoherency();
  MEDCouplingFieldDouble *f1=m->getMeasureField(false);
  MEDCouplingFieldDouble *f2=m2->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(f1->getNumberOfTuples(),3);
  CPPUNIT_ASSERT_EQUAL(f2->getNumberOfTuples(),3);
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
  m2->checkCoherency();
  f1=m->getMeasureField(false);
  f2=m2->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(f1->getNumberOfTuples(),6);
  CPPUNIT_ASSERT_EQUAL(f2->getNumberOfTuples(),6);
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
  m2->checkCoherency();
  f1=m->getMeasureField(false);
  f2=m2->getMeasureField(false);
  CPPUNIT_ASSERT_EQUAL(f1->getNumberOfTuples(),24);
  CPPUNIT_ASSERT_EQUAL(f2->getNumberOfTuples(),24);
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

void MEDCouplingBasicsTest::testDataArrayIntInvertO2NNO21()
{
  const int arr1[6]={2,0,4,1,5,3};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(6,1);
  std::copy(arr1,arr1+6,da->getPointer());
  DataArrayInt *da2=da->invertArrayO2N2N2O(6);
  CPPUNIT_ASSERT_EQUAL(6,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(6,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
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

void MEDCouplingBasicsTest::testKeepSetSelectedComponent1()
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
  DataArrayDouble *a2=a1->keepSelectedComponents(arr2V);
  CPPUNIT_ASSERT_EQUAL(6,a2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,a2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(0))=="bbbb");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(1))=="cccc");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(2))=="bbbb");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(3))=="cccc");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(4))=="aaaa");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(5))=="aaaa");
  const double expected1[30]={2.,3.,2.,3.,1.,1., 12.,13.,12.,13.,11.,11., 22.,23.,22.,23.,21.,21., 32.,33.,32.,33.,31.,31., 42.,43.,42.,43.,41.,41.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],a2->getIJ(0,i),1e-14);
  DataArrayInt *a3=a1->convertToIntArr();
  DataArrayInt *a4=a3->keepSelectedComponents(arr2V);
  CPPUNIT_ASSERT_EQUAL(6,a4->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,a4->getNumberOfTuples());
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
  DataArrayDouble *a5=a1->keepSelectedComponents(arr3V);
  a5->setInfoOnComponent(0,"eeee");
  a5->setInfoOnComponent(1,"ffff");
  const int arr4[2]={1,2};
  std::vector<int> arr4V(arr4,arr4+2);
  a2->setSelectedComponents(a5,arr4V);
  CPPUNIT_ASSERT_EQUAL(6,a2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,a2->getNumberOfTuples());
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(0))=="bbbb");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(1))=="eeee");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(2))=="ffff");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(3))=="cccc");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(4))=="aaaa");
  CPPUNIT_ASSERT(std::string(a2->getInfoOnComponent(5))=="aaaa");
  const double expected2[30]={2.,4.,3.,3.,1.,1., 12.,14.,13.,13.,11.,11., 22.,24.,23.,23.,21.,21., 32.,34.,33.,33.,31.,31., 42.,44.,43.,43.,41.,41.};
  for(int i=0;i<30;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],a2->getIJ(0,i),1e-14);
  DataArrayInt *a6=a5->convertToIntArr();
  a6->setInfoOnComponent(0,"eeee");
  a6->setInfoOnComponent(1,"ffff");
  a4->setSelectedComponents(a6,arr4V);
  CPPUNIT_ASSERT_EQUAL(6,a4->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,a4->getNumberOfTuples());
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
  a6->decrRef();
  a5->decrRef();
  a4->decrRef();
  a3->decrRef();
  a2->decrRef();
  a1->decrRef();
}

void MEDCouplingBasicsTest::testKeepSetSelectedComponent2()
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
  f1->checkCoherency();
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
  f2->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(6,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
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
  f5->checkCoherency();
  const int arr4[2]={1,2};
  std::vector<int> arr4V(arr4,arr4+2);
  f2->setSelectedComponents(f5,arr4V);
  CPPUNIT_ASSERT_EQUAL(6,f2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(5,f2->getNumberOfTuples());
  f2->checkCoherency();
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

void MEDCouplingBasicsTest::testDAIGetIdsEqual1()
{
  const int tab1[7]={5,-2,-4,-2,3,2,-2};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(7,1);
  std::copy(tab1,tab1+7,da->getPointer());
  DataArrayInt *da2=da->getIdsEqual(-2);
  CPPUNIT_ASSERT_EQUAL(3,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  const int expected1[3]={1,3,6};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,da2->getConstPointer()));
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest::testDAIGetIdsEqualList1()
{
  const int tab1[7]={5,-2,-4,-2,3,2,-2};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(7,1);
  std::copy(tab1,tab1+7,da->getPointer());
  const int tab2[3]={3,-2,0};
  std::vector<int> tab2V(tab2,tab2+3);
  DataArrayInt *da2=da->getIdsEqualList(tab2V);
  CPPUNIT_ASSERT_EQUAL(4,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  const int expected1[4]={1,3,4,6};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+4,da2->getConstPointer()));
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest::testDAFromNoInterlace1()
{
  const int tab1[15]={1,11,21,31,41,2,12,22,32,42,3,13,23,33,43};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(5,3);
  std::copy(tab1,tab1+15,da->getPointer());
  DataArrayInt *da2=da->fromNoInterlace();
  const int expected1[15]={1,2,3,11,12,13,21,22,23,31,32,33,41,42,43};
  CPPUNIT_ASSERT_EQUAL(5,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da2->getNumberOfComponents());// it's not a bug. Avoid to have 1 million components !
  CPPUNIT_ASSERT(std::equal(expected1,expected1+15,da2->getConstPointer()));
  DataArrayDouble *da3=da->convertToDblArr();
  DataArrayDouble *da4=da3->fromNoInterlace();
  CPPUNIT_ASSERT_EQUAL(5,da4->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da4->getNumberOfComponents());// it's not a bug. Avoid to have 1 million components !
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)expected1[i],da4->getIJ(0,i),1e-14);
  da4->decrRef();
  da3->decrRef();
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest::testDAToNoInterlace1()
{
  const int tab1[15]={1,2,3,11,12,13,21,22,23,31,32,33,41,42,43};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(5,3);
  std::copy(tab1,tab1+15,da->getPointer());
  DataArrayInt *da2=da->toNoInterlace();
  const int expected1[15]={1,11,21,31,41,2,12,22,32,42,3,13,23,33,43};
  CPPUNIT_ASSERT_EQUAL(5,da2->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da2->getNumberOfComponents());// it's not a bug. Avoid to have 1 million components !
  CPPUNIT_ASSERT(std::equal(expected1,expected1+15,da2->getConstPointer()));
  DataArrayDouble *da3=da->convertToDblArr();
  DataArrayDouble *da4=da3->toNoInterlace();
  CPPUNIT_ASSERT_EQUAL(5,da4->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,da4->getNumberOfComponents());// it's not a bug. Avoid to have 1 million components !
  for(int i=0;i<15;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)expected1[i],da4->getIJ(0,i),1e-14);
  da4->decrRef();
  da3->decrRef();
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest::testDAIsUniform1()
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
  DataArrayDouble *da2=da->convertToDblArr();
  CPPUNIT_ASSERT(da2->isUniform(1.,1e-12));
  da2->setIJ(1,0,1.+1.e-13);
  CPPUNIT_ASSERT(da2->isUniform(1.,1e-12));
  da2->setIJ(1,0,1.+1.e-11);
  CPPUNIT_ASSERT(!da2->isUniform(1.,1e-12));
  da2->decrRef();
  da->decrRef();
}

void MEDCouplingBasicsTest::testDADFromPolarToCart1()
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

void MEDCouplingBasicsTest::testDADFromCylToCart1()
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

void MEDCouplingBasicsTest::testDADFromSpherToCart1()
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

void MEDCouplingBasicsTest::testUnPolyze1()
{
  const int elts[8]={0,1,2,3,4,5,6,7};
  std::vector<int> eltsV(elts,elts+8);
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  mesh->convertToPolyTypes(eltsV);
  mesh->unPolyze();
  MEDCouplingUMesh *mesh2=build3DTargetMesh_1();
  mesh->checkCoherency();
  CPPUNIT_ASSERT(mesh->isEqual(mesh2,1e-12));
  mesh->convertToPolyTypes(eltsV);
  CPPUNIT_ASSERT(!mesh->isEqual(mesh2,1e-12));
  mesh->getNodalConnectivity()->setIJ(0,6,10);
  mesh->getNodalConnectivity()->setIJ(0,7,9);
  mesh->getNodalConnectivity()->setIJ(0,8,12);
  mesh->getNodalConnectivity()->setIJ(0,9,13);
  mesh->unPolyze();
  CPPUNIT_ASSERT(mesh->isEqual(mesh2,1e-12));
  mesh->convertToPolyTypes(eltsV);
  mesh->getNodalConnectivity()->setIJ(0,6,12);
  mesh->getNodalConnectivity()->setIJ(0,7,13);
  mesh->getNodalConnectivity()->setIJ(0,8,10);
  mesh->getNodalConnectivity()->setIJ(0,9,9);
  mesh->unPolyze();
  CPPUNIT_ASSERT(mesh->isEqual(mesh2,1e-12));
  mesh->convertToPolyTypes(eltsV);
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
  mesh->convertToPolyTypes(eltsV);
  CPPUNIT_ASSERT(!mesh->isEqual(mesh2,1e-12));
  mesh->unPolyze();
  CPPUNIT_ASSERT(mesh->isEqual(mesh2,1e-12));
  mesh->decrRef();
  mesh2->decrRef();
}

void MEDCouplingBasicsTest::testConvertDegeneratedCells1()
{
  MEDCouplingUMesh *mesh=build3DTargetMesh_1();
  int conn[32]={0,1,3,3,9,10,12,12, 0,1,3,4,9,9,9,9, 1,1,1,1,10,12,9,10, 10,11,12,9,1,1,1,1};
  mesh->allocateCells(4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+8);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+16);
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+24);
  mesh->finishInsertingCells();
  mesh->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(4,mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,mesh->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,mesh->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_HEXA8,mesh->getTypeOfCell(3));
  MEDCouplingFieldDouble *f1=mesh->getMeasureField(true);
  mesh->convertDegeneratedCells();
  mesh->checkCoherency();
  MEDCouplingFieldDouble *f2=mesh->getMeasureField(true);
  CPPUNIT_ASSERT_EQUAL(4,mesh->getNumberOfCells());
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

void MEDCouplingBasicsTest::testGetNodeIdsNearPoints1()
{
  MEDCouplingUMesh *mesh=build2DTargetMesh_1();
  DataArrayDouble *coords=mesh->getCoords();
  DataArrayDouble *tmp=DataArrayDouble::New();
  tmp->alloc(3,2);
  const double vals[6]={0.2,0.2,0.1,0.2,0.2,0.2};
  std::copy(vals,vals+6,tmp->getPointer());
  DataArrayDouble *tmp2=DataArrayDouble::aggregate(coords,tmp);
  tmp->decrRef();
  mesh->setCoords(tmp2);
  tmp2->decrRef();
  const double pts[6]={0.2,0.2,0.1,0.3,-0.3,0.7};
  std::vector<int> c=mesh->getNodeIdsNearPoint(pts,1e-7);
  CPPUNIT_ASSERT_EQUAL(3,(int)c.size());
  CPPUNIT_ASSERT_EQUAL(4,c[0]);
  CPPUNIT_ASSERT_EQUAL(9,c[1]);
  CPPUNIT_ASSERT_EQUAL(11,c[2]);
  c.clear();
  std::vector<int> cI;
  mesh->getNodeIdsNearPoints(pts,3,1e-7,c,cI);
  CPPUNIT_ASSERT_EQUAL(4,(int)cI.size());
  CPPUNIT_ASSERT_EQUAL(4,(int)c.size());
  CPPUNIT_ASSERT_EQUAL(4,c[0]);
  CPPUNIT_ASSERT_EQUAL(9,c[1]);
  CPPUNIT_ASSERT_EQUAL(11,c[2]);
  CPPUNIT_ASSERT_EQUAL(6,c[3]);
  CPPUNIT_ASSERT_EQUAL(0,cI[0]);
  CPPUNIT_ASSERT_EQUAL(3,cI[1]);
  CPPUNIT_ASSERT_EQUAL(3,cI[2]);
  CPPUNIT_ASSERT_EQUAL(4,cI[3]);
  mesh->decrRef();
}

void MEDCouplingBasicsTest::testFieldCopyTinyAttrFrom1()
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
void MEDCouplingBasicsTest::testExtrudedMesh5()
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
  i->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(36,i->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(37,i->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(12,i->getNumberOfCellsWithType(INTERP_KERNEL::NORM_TRI3));
  CPPUNIT_ASSERT_EQUAL(24,i->getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD4));
  const double expected1[3]={0.25,0.75,2.0625};
  MEDCouplingFieldDouble *j=i->getMeasureField(true);
  for(int i=0;i<12;i++)
    for(int k=0;k<3;k++)
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[k],j->getIJ(0,i*3+k),1e-10);
  const double expected2[72]={0.62200846792814113, 0.16666666666681595, 1.4513530918323276, 0.38888888888923495, 2.6293994326053212, 0.7045454545460802, 0.45534180126145435, 0.45534180126150181, 1.0624642029433926, 1.0624642029435025, 1.9248539780597826, 1.9248539780599816, 0.16666666666661334, 0.62200846792815856, 0.38888888888876294, 1.4513530918323678, 0.70454545454522521, 2.629399432605394, -0.16666666666674007, 0.62200846792812436, -0.38888888888906142, 1.4513530918322881, -0.70454545454576778, 2.6293994326052488, -0.45534180126154766, 0.45534180126140844, -1.0624642029436118, 1.0624642029432834, -1.9248539780601803, 1.9248539780595841, -0.62200846792817499, 0.1666666666665495, -1.451353091832408, 0.388888888888613, -2.6293994326054668, 0.70454545454495332, -0.62200846792810593, -0.16666666666680507, -1.451353091832247, -0.38888888888921297, -2.6293994326051746, -0.70454545454604123, -0.45534180126135926, -0.45534180126159562, -1.0624642029431723, -1.0624642029437235, -1.9248539780593836, -1.9248539780603811, -0.1666666666664828, -0.62200846792819242, -0.38888888888846079, -1.4513530918324489, -0.70454545454467987, -2.6293994326055397, 0.16666666666687083, -0.62200846792808862, 0.38888888888936374, -1.4513530918322073, 0.70454545454631357, -2.6293994326051022, 0.45534180126164348, -0.45534180126131207, 1.0624642029438327, -1.0624642029430627, 1.9248539780605791, -1.9248539780591853, 0.62200846792821063, -0.16666666666641802, 1.4513530918324888, -0.38888888888831086, 2.6293994326056125, -0.70454545454440853};
  DataArrayDouble *m=i->getBarycenterAndOwner();
  for(int i=0;i<72;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],m->getIJ(0,i),1e-10);
  //
  m->decrRef();
  j->decrRef();
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

/*!
 * 1D -> 2D extrusion without rotation
 */
void MEDCouplingBasicsTest::testExtrudedMesh6()
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
  g->checkCoherency();
  const double expected1[]={ 0.4330127018922193, 0.4330127018922193, 0.649519052838329, 1.2990381056766578, 1.299038105676658, 1.948557158514987, 2.1650635094610955, 2.1650635094610964, 3.2475952641916446, 3.031088913245533, 3.0310889132455352, 4.546633369868303 };
  MEDCouplingFieldDouble *f1=g->getMeasureField(true);
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],f1->getIJ(0,i),1e-12);
  
  const double expected2[]={0.625, 0.21650635094610962, 1.625, 0.21650635094610959, 2.8750000000000004, 0.21650635094610965, 1.1250000000000002, 1.0825317547305482, 2.125, 1.0825317547305482, 3.3750000000000004, 1.0825317547305484, 2.125, 2.8145825622994254, 3.125, 2.8145825622994254, 4.375, 2.8145825622994254, 3.6250000000000009, 5.4126587736527414, 4.625, 5.4126587736527414, 5.875, 5.4126587736527414};
  DataArrayDouble *f2=g->getBarycenterAndOwner();
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
void MEDCouplingBasicsTest::testExtrudedMesh7()
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

void MEDCouplingBasicsTest::testSimplexize1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  std::vector<int> v(1);
  v[0]=3;
  m->convertToPolyTypes(v);
  DataArrayInt *da=m->simplexize(0);
  CPPUNIT_ASSERT_EQUAL(7,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da->getNumberOfComponents());
  const int expected2[7]={0,0,1,2,3,4,4};
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],da->getIJ(i,0));
  m->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(7,m->getNumberOfCells());
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
  std::set<INTERP_KERNEL::NormalizedCellType> types=m->getAllTypes();
  CPPUNIT_ASSERT_EQUAL(2,(int)types.size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*(types.begin()));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_POLYGON,*(++(types.begin())));
  f->decrRef();
  da->decrRef();
  m->decrRef();
  //
  m=build3DSurfTargetMesh_1();
  v[0]=3;
  m->convertToPolyTypes(v);
  da=m->simplexize(1);
  CPPUNIT_ASSERT_EQUAL(7,da->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,da->getNumberOfComponents());
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expected2[i],da->getIJ(i,0));
  m->checkCoherency();
  types=m->getAllTypes();
  CPPUNIT_ASSERT_EQUAL(2,(int)types.size());
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_TRI3,*(types.begin()));
  CPPUNIT_ASSERT_EQUAL(INTERP_KERNEL::NORM_POLYGON,*(++(types.begin())));
  CPPUNIT_ASSERT_EQUAL(7,m->getNumberOfCells());
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

void MEDCouplingBasicsTest::testSimplexize2()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  std::vector<int> v(1);
  v[0]=3;
  m->convertToPolyTypes(v);
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setMesh(m);
  DataArrayDouble *arr=DataArrayDouble::New();
  const double arr1[10]={10.,110.,20.,120.,30.,130.,40.,140.,50.,150.};
  arr->alloc(5,2);
  std::copy(arr1,arr1+10,arr->getPointer());
  f1->setArray(arr);
  arr->decrRef();
  //
  f1->checkCoherency();
  CPPUNIT_ASSERT(f1->simplexize(0));
  f1->checkCoherency();
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

void MEDCouplingBasicsTest::testDAMeld1()
{
  DataArrayDouble *da1=DataArrayDouble::New();
  da1->alloc(7,2);
  DataArrayDouble *da2=DataArrayDouble::New();
  da2->alloc(7,1);
  //
  da1->fillWithValue(7.);
  da2->iota(0.);
  DataArrayDouble *da3=da2->applyFunc(3,"10*x*IVec+100*x*JVec+1000*x*KVec");
  //
  da1->setInfoOnComponent(0,"c0da1");
  da1->setInfoOnComponent(1,"c1da1");
  da3->setInfoOnComponent(0,"c0da3");
  da3->setInfoOnComponent(1,"c1da3");
  da3->setInfoOnComponent(2,"c2da3");
  //
  DataArrayDouble *da1C=da1->deepCpy();
  da1->meldWith(da3);
  CPPUNIT_ASSERT_EQUAL(5,da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,da1->getNumberOfTuples());
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
  DataArrayInt *dai1=da1C->convertToIntArr();
  DataArrayInt *dai3=da3->convertToIntArr();
  dai1->meldWith(dai3);
  CPPUNIT_ASSERT_EQUAL(5,dai1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,dai1->getNumberOfTuples());
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(0)=="c0da1");
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(1)=="c1da1");
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(2)=="c0da3");
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(3)=="c1da3");
  CPPUNIT_ASSERT(dai1->getInfoOnComponent(4)=="c2da3");
  for(int i=0;i<35;i++)
    CPPUNIT_ASSERT_EQUAL((int)expected1[i],dai1->getIJ(0,i));
  // test of static method DataArrayDouble::meld
  DataArrayDouble *da4=DataArrayDouble::meld(da1C,da3);
  CPPUNIT_ASSERT_EQUAL(5,da4->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,da4->getNumberOfTuples());
  CPPUNIT_ASSERT(da4->getInfoOnComponent(0)=="c0da1");
  CPPUNIT_ASSERT(da4->getInfoOnComponent(1)=="c1da1");
  CPPUNIT_ASSERT(da4->getInfoOnComponent(2)=="c0da3");
  CPPUNIT_ASSERT(da4->getInfoOnComponent(3)=="c1da3");
  CPPUNIT_ASSERT(da4->getInfoOnComponent(4)=="c2da3");
  for(int i=0;i<35;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],da4->getIJ(0,i),1e-10);
  // test of static method DataArrayInt::meld
  dai1->decrRef();
  dai1=da1C->convertToIntArr();
  DataArrayInt *dai4=DataArrayInt::meld(dai1,dai3);
  CPPUNIT_ASSERT_EQUAL(5,dai4->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(7,dai4->getNumberOfTuples());
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
  dai3->decrRef();
  dai1->decrRef();
  da1C->decrRef();
  da1->decrRef();
  da2->decrRef();
  da3->decrRef();
}

void MEDCouplingBasicsTest::testFieldMeld1()
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
  f1->checkCoherency();
  //
  MEDCouplingFieldDouble *f2=f1->deepCpy();
  f2->setMesh(f1->getMesh());
  f2->checkCoherency();
  f2->changeNbOfComponents(2,5.);
  (*f2)=5.;
  f2->getArray()->setInfoOnComponent(0,"bbb");
  f2->getArray()->setInfoOnComponent(1,"ccc");
  f2->checkCoherency();
  //
  MEDCouplingFieldDouble *f3=MEDCouplingFieldDouble::meldFields(f2,f1);
  f3->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(5,f3->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,f3->getNumberOfComponents());
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
  MEDCouplingFieldDouble *f6=MEDCouplingFieldDouble::meldFields(f4,f5);
  f6->checkCoherency();
  CPPUNIT_ASSERT_EQUAL(5,f6->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(3,f6->getNumberOfComponents());
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

void MEDCouplingBasicsTest::testMergeNodes2()
{
  MEDCouplingUMesh *m1=build2DTargetMesh_1();
  MEDCouplingUMesh *m2=build2DTargetMesh_1();
  const double vec[2]={0.002,0.};
  m2->translate(vec);
  //
  std::vector<const MEDCouplingUMesh *> tmp(2);
  tmp[0]=m1;
  tmp[1]=m2;
  MEDCouplingUMesh *m3=MEDCouplingUMesh::mergeUMeshes(tmp);
  bool b;
  int newNbOfNodes;
  DataArrayInt *da=m3->mergeNodes2(0.01,b,newNbOfNodes);
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

void MEDCouplingBasicsTest::testMergeField2()
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
  MEDCouplingFieldDouble *f4=MEDCouplingFieldDouble::mergeFields(tmp);
  CPPUNIT_ASSERT_EQUAL(15,f4->getMesh()->getNumberOfCells());
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

void MEDCouplingBasicsTest::testDAIBuildComplement1()
{
  DataArrayInt *a=DataArrayInt::New();
  const int tab[4]={3,1,7,8};
  a->alloc(4,1);
  std::copy(tab,tab+4,a->getPointer());
  DataArrayInt *b=a->buildComplement(12);
  CPPUNIT_ASSERT_EQUAL(8,b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,b->getNumberOfComponents());
  const int expected1[8]={0,2,4,5,6,9,10,11};
  for(int i=0;i<8;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],b->getIJ(0,i));
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest::testDAIBuildUnion1()
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
  CPPUNIT_ASSERT_EQUAL(7,b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,b->getNumberOfComponents());
  const int expected1[7]={0,1,3,5,7,8,18};
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],b->getIJ(0,i));
  c->decrRef();
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest::testDAIBuildIntersection1()
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
  CPPUNIT_ASSERT_EQUAL(2,b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,b->getNumberOfComponents());
  const int expected1[2]={3,8};
  for(int i=0;i<2;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],b->getIJ(0,i));
  c->decrRef();
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest::testDAIDeltaShiftIndex1()
{
  DataArrayInt *a=DataArrayInt::New();
  const int tab[7]={1,3,6,7,7,9,15};
  a->alloc(7,1);
  std::copy(tab,tab+7,a->getPointer());
  DataArrayInt *b=a->deltaShiftIndex();
  CPPUNIT_ASSERT_EQUAL(6,b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,b->getNumberOfComponents());
  const int expected1[6]={2,3,1,0,2,6};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],b->getIJ(0,i));
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest::testDaDoubleSelectByTupleIdSafe1()
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
  CPPUNIT_ASSERT_EQUAL(5,b->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,b->getNumberOfComponents());
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
  CPPUNIT_ASSERT_EQUAL(5,d->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(2,d->getNumberOfComponents());
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

void MEDCouplingBasicsTest::testAreCellsIncludedIn1()
{
  MEDCouplingUMesh *m=build3DSurfTargetMesh_1();
  const int pt[2]={1,3};
  MEDCouplingUMesh *m2=(MEDCouplingUMesh *)m->buildPartOfMySelf(pt,pt+2,true);
  DataArrayInt *tmp;
  CPPUNIT_ASSERT(m->areCellsIncludedIn(m2,0,tmp));
  CPPUNIT_ASSERT_EQUAL(2,tmp->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,tmp->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(pt[0],tmp->getIJ(0,0));
  CPPUNIT_ASSERT_EQUAL(pt[1],tmp->getIJ(0,1));
  tmp->decrRef();
  CPPUNIT_ASSERT(!m2->areCellsIncludedIn(m,0,tmp));
  tmp->decrRef();
  m2->decrRef();
  m->decrRef();
}

void MEDCouplingBasicsTest::testDAIBuildSubstraction1()
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
  CPPUNIT_ASSERT_EQUAL(3,c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,c->getNumberOfComponents());
  const int expected1[3]={2,6,8};
  CPPUNIT_ASSERT(std::equal(expected1,expected1+3,c->getConstPointer()));
  //
  c->decrRef();
  b->decrRef();
  a->decrRef();
}

void MEDCouplingBasicsTest::testBuildOrthogonalField2()
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
  CPPUNIT_ASSERT_EQUAL(2,da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(13,da1->getNumberOfTuples());
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

void MEDCouplingBasicsTest::testUMInsertNextCell1()
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
  CPPUNIT_ASSERT_THROW(targetMesh->insertNextCell(INTERP_KERNEL::NORM_POINT0,1,targetConn),INTERP_KERNEL::Exception);
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
  targetMesh->checkCoherency();
  targetMesh->decrRef();
}

void MEDCouplingBasicsTest::testFieldOperatorDivDiffComp1()
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
  f2->checkCoherency();
  //
  MEDCouplingFieldDouble *f3=(*f1)/(*f2);
  CPPUNIT_ASSERT_THROW((*f2)/(*f1),INTERP_KERNEL::Exception);
  f3->checkCoherency();
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

void MEDCouplingBasicsTest::testDARearrange1()
{
  DataArrayInt *da1=DataArrayInt::New();
  da1->alloc(12,1);
  da1->iota(0);
  const int *ptr=da1->getConstPointer();
  //
  CPPUNIT_ASSERT_EQUAL(12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(1,da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(12,da1->getNumberOfTuples());
  da1->rearrange(4);
  CPPUNIT_ASSERT(ptr==da1->getConstPointer());
  CPPUNIT_ASSERT_EQUAL(12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(4,da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,da1->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(i,da1->getIJ(0,i));
  //
  da1->rearrange(6);
  CPPUNIT_ASSERT(ptr==da1->getConstPointer());
  CPPUNIT_ASSERT_EQUAL(12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(6,da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(2,da1->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(i,da1->getIJ(0,i));
  //
  CPPUNIT_ASSERT_THROW(da1->rearrange(7),INTERP_KERNEL::Exception);
  //
  da1->rearrange(12);
  CPPUNIT_ASSERT(ptr==da1->getConstPointer());
  CPPUNIT_ASSERT_EQUAL(12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(12,da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(1,da1->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(i,da1->getIJ(0,i));
  //
  da1->rearrange(3);
  CPPUNIT_ASSERT(ptr==da1->getConstPointer());
  CPPUNIT_ASSERT_EQUAL(12,da1->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(3,da1->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(4,da1->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(i,da1->getIJ(0,i));
  //double
  DataArrayDouble *da2=da1->convertToDblArr();
  da1->decrRef();
  const double *ptr2=da2->getConstPointer();
  //
  CPPUNIT_ASSERT_EQUAL(12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(3,da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(4,da2->getNumberOfTuples());
  da2->rearrange(4);
  CPPUNIT_ASSERT(ptr2==da2->getConstPointer());
  CPPUNIT_ASSERT_EQUAL(12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(4,da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(3,da2->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)i,da2->getIJ(0,i),1e-14);
  //
  da2->rearrange(6);
  CPPUNIT_ASSERT(ptr2==da2->getConstPointer());
  CPPUNIT_ASSERT_EQUAL(12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(6,da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(2,da2->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)i,da2->getIJ(0,i),1e-14);
  //
  CPPUNIT_ASSERT_THROW(da2->rearrange(7),INTERP_KERNEL::Exception);
  //
  da2->rearrange(1);
  CPPUNIT_ASSERT(ptr2==da2->getConstPointer());
  CPPUNIT_ASSERT_EQUAL(12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(1,da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(12,da2->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)i,da2->getIJ(0,i),1e-14);
  //
  da2->rearrange(3);
  CPPUNIT_ASSERT(ptr2==da2->getConstPointer());
  CPPUNIT_ASSERT_EQUAL(12,da2->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(3,da2->getNumberOfComponents());
  CPPUNIT_ASSERT_EQUAL(4,da2->getNumberOfTuples());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL((double)i,da2->getIJ(0,i),1e-14);
  da2->decrRef();
}

void MEDCouplingBasicsTest::testGetDifferentValues1()
{
  DataArrayInt *da1=DataArrayInt::New();
  const int arr[12]={1,2,3,2,2,3,5,1,5,5,2,2};
  da1->alloc(4,3);
  std::copy(arr,arr+12,da1->getPointer());
  std::set<int> s=da1->getDifferentValues();
  const int expected1[4]={1,2,3,5};
  CPPUNIT_ASSERT_EQUAL(4,(int)s.size());
  CPPUNIT_ASSERT(std::equal(expected1,expected1+4,s.begin()));
  da1->decrRef();
}

