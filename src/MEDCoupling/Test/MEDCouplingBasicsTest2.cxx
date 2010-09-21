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
  MEDCouplingUMesh *m5=m3->buildExtrudedMeshFromThis(m4,0);
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
  //buildWeightingField
  MEDCouplingFieldDouble *f4=f1->buildWeightingField(false);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(-0.2,f4->accumulate(0),1e-12);
  f4->decrRef();
  f4=f1->buildWeightingField(true);
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
  DataArrayDouble *a2=a1->deepCopy();
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
