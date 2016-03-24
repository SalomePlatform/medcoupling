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

#include "ThreeDSurfProjectionTest.hxx"
#include "PlanarIntersector.txx"

class MyMeshType
{
public:
  static const int MY_SPACEDIM=3;
  static const int MY_MESHDIM=3;
  static const INTERP_KERNEL::NumberingPolicy My_numPol=INTERP_KERNEL::ALL_C_MODE;
  typedef int MyConnType;
};

class MyMatrixType
{
};

void INTERP_TEST::ThreeDSurfProjectionTest::test1()
{
  // Two triangles coo and coo2 are perfectly // each others with a distance equal to 1e-6.
  // A little rotation to make it more funny.
  //coo=DataArrayDouble([0.,0.,0.,1.,0.,0.,0.,1.,0.],3,3)
  //eps=1e-6
  //coo2=DataArrayDouble([0.,0.,eps,1.,0.,eps,0.,1.,eps],3,3)
  //MEDCouplingPointSet.Rotate3DAlg([0.,0.,0.],[2.,1.,3.],0.3,coo)
  //MEDCouplingPointSet.Rotate3DAlg([0.,0.,0.],[2.,1.,3.],0.3,coo2)
  const double coo[9]={0.,0.,0.,0.96809749223257568,0.24332379388106262,-0.059839592782071335,-0.23056279077409292,0.95852673990234838,0.16753294721527912};
  const double coo2[9]={9.8122602102980502e-08,-1.4839144255482456e-7,9.8404874611628791e-7,0.96809759035517784,0.24332364548962007,-0.059838608733325221,-0.23056269265149082,0.9585265915109058,0.16753393126402524};
  double *tmp0(new double[9]),*tmp1(new double[9]);
  int ret;
  //eps=1e-2. eps is a tolerance to detect that two points are the same or not in a same polygon.
  // here the max 3D distance is 1e-5 > 1e-6 so 1 is expected
  std::copy(coo,coo+9,tmp0);
  std::copy(coo2,coo2+9,tmp1);
  ret=INTERP_KERNEL::PlanarIntersector<MyMeshType,MyMatrixType>::Projection(tmp0,tmp1,3,3,1e-2,1e-5/* <- */,-1.,0.5,true);
  CPPUNIT_ASSERT_EQUAL(1,ret);
  const double expected0[9]={0.,0.,0.,1.,0.,0.,0.,1.,0.};
  for(int i=0;i<9;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected0[i],tmp0[i],1e-15);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected0[i],tmp1[i],1e-15);
    }
  // here the max 3D distance is 1e-8 < 1e-6 so 0 is expected
  std::copy(coo,coo+9,tmp0);
  std::copy(coo2,coo2+9,tmp1);
  ret=INTERP_KERNEL::PlanarIntersector<MyMeshType,MyMatrixType>::Projection(tmp0,tmp1,3,3,1e-2,1e-8/* <- */,-1.,0.5,true);
  CPPUNIT_ASSERT_EQUAL(0,ret);
  // here testing when max 3D distance is 1e-5 > 1e-6 with inverted cells
  std::copy(coo,coo+9,tmp0);
  std::copy(coo2,coo2+3,tmp1+6); std::copy(coo2+3,coo2+6,tmp1+3); std::copy(coo2+6,coo2+9,tmp1);
  ret=INTERP_KERNEL::PlanarIntersector<MyMeshType,MyMatrixType>::Projection(tmp0,tmp1,3,3,1e-2,1e-5/* <- */,-1.,0.5,true);
  CPPUNIT_ASSERT_EQUAL(-1,ret);
  const double expected1[9]={-0.7071067811865476,-0.7071067811865476,0.,0.,-1.4142135623730951,0.,-1.4142135623730951,-1.4142135623730951,0.};
  const double expected2[9]={-1.4142135623730951,-1.4142135623730951,0.,0.,-1.4142135623730951,0.,-0.7071067811865476,-0.7071067811865476,0.};
  for(int i=0;i<9;i++)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],tmp0[i],1e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],tmp1[i],1e-14);
    }
  //
  delete [] tmp0;
  delete [] tmp1;
}

void INTERP_TEST::ThreeDSurfProjectionTest::test2()
{// here the two triangles have their center of inertia very close (eps) but the angle between the two planes is "big"
  //coo=DataArrayDouble([0.,0.,0.,1.,0.,0.,0.,1.,0.],3,3)
  //coocpy=coo.deepCopy()
  //MEDCouplingPointSet.Rotate3DAlg([0.,0.,0.],[-1,-1.,0.],pi/3,coocpy)
  //coocpy+=[eps*sqrt(3)/2,eps/2,eps*0.]
  //
  const double coo[9]={0.,0.,0.,0.96809749223257568,0.24332379388106262,-0.059839592782071335,-0.23056279077409292,0.95852673990234838,0.16753294721527912};
  const double coo2[9]={7.2311562622637225e-07,6.8998795679738294e-07,3.1943866106249849e-08,0.72852072144314628,0.33125439126063028,0.5996079016637561,0.0090154262465889021,0.87059752249869415,-0.49191448334281612};
  double *tmp0(new double[9]),*tmp1(new double[9]);
  int ret;
  //eps=1e-2. eps is a tolerance to detect that two points are the same or not in a same polygon.
  // here the max 3D distance is 1e-5 > 1e-6 so 1 is expected
  std::copy(coo,coo+9,tmp0);
  std::copy(coo2,coo2+9,tmp1);
  ret=INTERP_KERNEL::PlanarIntersector<MyMeshType,MyMatrixType>::Projection(tmp0,tmp1,3,3,1e-2,1e-5/* <- */,-1.,0.5,true);
  CPPUNIT_ASSERT_EQUAL(1,ret);
  // here the max 3D distance is 1e-8 < 1e-6 so 0 is expected
  std::copy(coo,coo+9,tmp0);
  std::copy(coo2,coo2+9,tmp1);
  ret=INTERP_KERNEL::PlanarIntersector<MyMeshType,MyMatrixType>::Projection(tmp0,tmp1,3,3,1e-2,1e-8/* <- */,-1.,0.5,true);
  CPPUNIT_ASSERT_EQUAL(0,ret);
  // again max 3D distance is 1e-5 > 1e-6 so 1 is expected
  std::copy(coo,coo+9,tmp0);
  std::copy(coo2,coo2+9,tmp1);
  ret=INTERP_KERNEL::PlanarIntersector<MyMeshType,MyMatrixType>::Projection(tmp0,tmp1,3,3,1e-2,1e-5/* <- */,-1.,0.5,true);
  CPPUNIT_ASSERT_EQUAL(1,ret);
  // again max 3D distance is 1e-5 > 1e-6 but minDot set to 0.8. 0 expected. because the angle is pi/4 so cos(pi/3) > 0.8
  std::copy(coo,coo+9,tmp0);
  std::copy(coo2,coo2+9,tmp1);
  ret=INTERP_KERNEL::PlanarIntersector<MyMeshType,MyMatrixType>::Projection(tmp0,tmp1,3,3,1e-2,1e-5/* <- */,0.8/* <- */,0.5,true);
  CPPUNIT_ASSERT_EQUAL(0,ret);
  // again max 3D distance is 1e-5 > 1e-6 but minDot set to 0.7. 1 expected. because the angle is pi/4 so cos(pi/3) < 0.49
  std::copy(coo,coo+9,tmp0);
  std::copy(coo2,coo2+9,tmp1);
  ret=INTERP_KERNEL::PlanarIntersector<MyMeshType,MyMatrixType>::Projection(tmp0,tmp1,3,3,1e-2,1e-5/* <- */,0.49/* <- */,0.5,true);
  CPPUNIT_ASSERT_EQUAL(1,ret);
  // again max 3D distance is 1e-5 > 1e-6 but minDot set to 0.7. 0 expected. because the angle is pi/4 so cos(pi/3) > 0.51
  std::copy(coo,coo+9,tmp0);
  std::copy(coo2,coo2+9,tmp1);
  ret=INTERP_KERNEL::PlanarIntersector<MyMeshType,MyMatrixType>::Projection(tmp0,tmp1,3,3,1e-2,1e-5/* <- */,0.51/* <- */,0.5,true);
  CPPUNIT_ASSERT_EQUAL(0,ret);
  //
  delete [] tmp0;
  delete [] tmp1;
}
