//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#include "MEDCouplingRemapperTest.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingRemapper.hxx"

#include "MEDCouplingBasicsTest.hxx"

#include <cmath>

using namespace ParaMEDMEM;

void MEDCouplingRemapperTest::test2DInterpP0P0_1()
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingBasicsTest::build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=MEDCouplingBasicsTest::build2DTargetMesh_1();
  //
  MEDCouplingRemapper remapper;
  remapper.setPrecision(1e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  
  MEDCouplingFieldDouble *srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  double *ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  MEDCouplingFieldDouble *trgfield=remapper.transferField(srcField,4.57);
  const double *values=trgfield->getArray()->getConstPointer();
  const double valuesExpected[5]={7.5 ,7. ,7.,8.,7.5};
  CPPUNIT_ASSERT_EQUAL(5,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(IntegralGlobConstraint);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected2[5]={3.75 ,1.75 ,1.75,4.,3.75};
  CPPUNIT_ASSERT_EQUAL(5,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(ConservativeVolumic);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(5,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(IntegralGlobConstraint);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(5,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(Integral);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(5,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(RevIntegral);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(5,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->decrRef();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingRemapperTest::test2DInterpP0P0R_1()
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingBasicsTest::build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=MEDCouplingBasicsTest::build2DTargetMesh_1();
  //
  MEDCouplingRemapper remapper;
  remapper.setPrecision(1e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  
  MEDCouplingFieldDouble *targetField=MEDCouplingFieldDouble::New(ON_CELLS);
  targetField->setNature(ConservativeVolumic);
  targetField->setMesh(sourceMesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(targetMesh->getNumberOfCells(),1);
  targetField->setArray(array);
  double *ptr=array->getPointer();
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  //
  MEDCouplingFieldDouble *srcfield=remapper.reverseTransferField(targetField,4.57);
  const double *values=srcfield->getArray()->getConstPointer();
  const double valuesExpected[2]={8.75 ,9.5};
  CPPUNIT_ASSERT_EQUAL(2,srcfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,srcfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  srcfield->decrRef();
  //
  targetField->setNature(IntegralGlobConstraint);
  srcfield=remapper.reverseTransferField(targetField,4.57);
  values=srcfield->getArray()->getConstPointer();
  const double valuesExpected2[2]={26., 19.};
  CPPUNIT_ASSERT_EQUAL(2,srcfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,srcfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  srcfield->decrRef();
  //
  targetField->decrRef();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingRemapperTest::test2DInterpMultiMethods()
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingBasicsTest::build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=MEDCouplingBasicsTest::build2DTargetMesh_1();
  //
  MEDCouplingRemapper remapper;
  remapper.setPrecision(1e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  
  MEDCouplingFieldDouble *srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  double *ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  MEDCouplingFieldDouble *trgfield=remapper.transferField(srcField,4.57);
  const double *values=trgfield->getArray()->getConstPointer();
  const double valuesExpected[5]={7.5 ,7. ,7.,8.,7.5};
  CPPUNIT_ASSERT_EQUAL(5,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  srcField->decrRef();
  //
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P1P0"));
  srcField=MEDCouplingFieldDouble::New(ON_NODES);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfNodes(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfNodes();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected2[5]={7.,7.666666666666667,8.6666666666666661,8.8333333333333339,10.};
  CPPUNIT_ASSERT_EQUAL(5,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  srcField->decrRef();
  //
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(targetMesh,sourceMesh,"P0P1"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(targetMesh);
  array=DataArrayDouble::New();
  array->alloc(targetMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected3[4]={7.5,8.5,10.,10.625};
  CPPUNIT_ASSERT_EQUAL(4,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected3[i0],values[i0],1e-12);
  trgfield->decrRef();
  srcField->decrRef();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
  //
  sourceMesh=MEDCouplingBasicsTest::build2DSourceMesh_1();
  targetMesh=MEDCouplingBasicsTest::build2DTargetMesh_2();
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P1P1"));
  srcField=MEDCouplingFieldDouble::New(ON_NODES);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfNodes(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfNodes();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected4[9]={ 7.,7.35,8.,7.7,8.2857142857142865,
                                    9.5333333333333332,9.,9.7666666666666657,10.};
  CPPUNIT_ASSERT_EQUAL(9,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<9;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected4[i0],values[i0],1e-12);
  trgfield->decrRef();
  srcField->decrRef();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingRemapperTest::testMultiDimCombi()
{
  // ------------- 2D
  MEDCouplingUMesh *sourceMesh=MEDCouplingBasicsTest::build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=MEDCouplingBasicsTest::build2DTargetMesh_1();
  //
  MEDCouplingRemapper remapper;
  remapper.setPrecision(1e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  MEDCouplingFieldDouble *srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  double *ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  MEDCouplingFieldDouble *trgField=remapper.transferField(srcField,4.57);
  const double *values=trgField->getArray()->getConstPointer();
  const double valuesExpected[5]={7.5 ,7. ,7.,8.,7.5};
  CPPUNIT_ASSERT_EQUAL(5,trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgField->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgField->decrRef();
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  // ------------- 3D Surf
  sourceMesh=MEDCouplingBasicsTest::build3DSurfSourceMesh_1();
  targetMesh=MEDCouplingBasicsTest::build3DSurfTargetMesh_1();
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+8);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  CPPUNIT_ASSERT_EQUAL(5,trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgField->getArray()->getNumberOfComponents());
  const double valuesExpected2[5]={8.5,8.,8.,9.,8.5};
  values=trgField->getArray()->getConstPointer();
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgField->decrRef();
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  // ------------- 3D
  sourceMesh=MEDCouplingBasicsTest::build3DSourceMesh_1();
  targetMesh=MEDCouplingBasicsTest::build3DTargetMesh_1();
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  CPPUNIT_ASSERT_EQUAL(8,trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgField->getArray()->getNumberOfComponents());
  const double valuesExpected3[8]={13.166666666666668, 13.888888888888888, 10.722222222222223, 10.870370370370372,
                                   14.555555555555555, 13.888888888888889, 14.444444444444443, 11.72222222222222};
  values=trgField->getArray()->getConstPointer();
  for(int i0=0;i0<8;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected3[i0],values[i0],1e-12);
  trgField->decrRef();
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  // ------------- 3D -> 1D
  sourceMesh=MEDCouplingBasicsTest::build3DTargetMesh_1();
  targetMesh=MEDCouplingBasicsTest::build1DTargetMesh_1();
  remapper.setIntersectionType(INTERP_KERNEL::PointLocator);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  CPPUNIT_ASSERT_EQUAL(8,trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgField->getArray()->getNumberOfComponents());
  const double valuesExpected4[8]={7.,11.,8.,12.,9.,13.,10.,14.};
  values=trgField->getArray()->getConstPointer();
  for(int i0=0;i0<8;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected4[i0],values[i0],1e-12);
  trgField->decrRef();
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  // ------------- 1D -> 3D
  sourceMesh=MEDCouplingBasicsTest::build1DTargetMesh_1();
  targetMesh=MEDCouplingBasicsTest::build3DTargetMesh_1();
  remapper.setIntersectionType(INTERP_KERNEL::PointLocator);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  CPPUNIT_ASSERT_EQUAL(8,trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgField->getArray()->getNumberOfComponents());
  const double valuesExpected5[8]={7.,9.,11.,13.,8.,10.,12.,14.};
  values=trgField->getArray()->getConstPointer();
  for(int i0=0;i0<8;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected5[i0],values[i0],1e-12);
  trgField->decrRef();
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  // ------------- 2D -> 1D
  sourceMesh=MEDCouplingBasicsTest::build2DTargetMesh_1();
  targetMesh=build1DTargetMesh_2();
  remapper.setIntersectionType(INTERP_KERNEL::PointLocator);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  const double valuesExpected8[5]={9.,8.,11.,7.,11.};
  values=trgField->getArray()->getConstPointer();
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected8[i0],values[i0],1e-12);
  trgField->decrRef();
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  // ------------- 1D -> 2D
  sourceMesh=build1DTargetMesh_2();
  targetMesh=MEDCouplingBasicsTest::build2DTargetMesh_1();
  remapper.setIntersectionType(INTERP_KERNEL::PointLocator);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  const double valuesExpected9[5]={10.,8.,7.,4.57,10.};
  values=trgField->getArray()->getConstPointer();
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected9[i0],values[i0],1e-12);
  trgField->decrRef();
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  // ------------- 2D -> -1D
  sourceMesh=MEDCouplingBasicsTest::build2DTargetMesh_1();
  targetMesh=MEDCouplingUMesh::New("an example of -1 D mesh",-1);
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  trgField=remapper.transferField(srcField,4.57);
  values=trgField->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(1,trgField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgField->getNumberOfComponents());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,values[0],1e-14);
  srcField->decrRef();
  srcField=remapper.reverseTransferField(trgField,4.220173);
  CPPUNIT_ASSERT_EQUAL(5,srcField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,srcField->getNumberOfComponents());
  values=srcField->getArray()->getConstPointer();
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,values[i0],1e-14);
  srcField->decrRef();
  trgField->setNature(Integral);
  srcField=remapper.reverseTransferField(trgField,4.220173);
  CPPUNIT_ASSERT_EQUAL(5,srcField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,srcField->getNumberOfComponents());
  values=srcField->getArray()->getConstPointer();
  const double valuesExpected6[5]={2.28125,1.140625,1.140625,2.28125,2.28125};
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected6[i0],values[i0],1e-14);
  srcField->decrRef();
  trgField->decrRef();
  // ------------- -1D -> 2D
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(targetMesh,sourceMesh,"P0P0"));
  trgField=MEDCouplingFieldDouble::New(ON_CELLS);
  trgField->setNature(ConservativeVolumic);
  trgField->setMesh(targetMesh);
  array=DataArrayDouble::New();
  array->alloc(targetMesh->getNumberOfCells(),1);
  trgField->setArray(array);
  ptr=array->getPointer();
  ptr[0]=7.;
  array->decrRef();
  srcField=remapper.transferField(trgField,4.221073);
  values=srcField->getArray()->getConstPointer();
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(7.,values[i0],1e-14);
  srcField->decrRef();
  trgField->setNature(IntegralGlobConstraint);
  srcField=remapper.transferField(trgField,4.221073);
  values=srcField->getArray()->getConstPointer();
  const double valuesExpected7[5]={1.75,0.875,0.875,1.75,1.75};
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected7[i0],values[i0],1e-14);
  srcField->decrRef();
  trgField->setNature(Integral);
  srcField=remapper.transferField(trgField,4.221073);
  values=srcField->getArray()->getConstPointer();
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected7[i0],values[i0],1e-14);
  //
  srcField->decrRef();
  trgField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingRemapperTest::testNatureOfField()
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingBasicsTest::build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_3();
  //
  MEDCouplingRemapper remapper;
  remapper.setPrecision(1e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  MEDCouplingFieldDouble *srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  double *ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  MEDCouplingFieldDouble *trgfield=remapper.transferField(srcField,4.220173);
  const double *values=trgfield->getArray()->getConstPointer();
  const double valuesExpected[4]={7.75, 7.0625, 4.220173,8.0};
  CPPUNIT_ASSERT_EQUAL(4,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(IntegralGlobConstraint);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected2[4]={2.8374999999999999, 7.3624999999999998, 4.220173, 4.7999999999999998};
  CPPUNIT_ASSERT_EQUAL(4,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(Integral);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected3[4]={1.24, 4.5199999999999996, 4.220173, 1.9199999999999999};
  CPPUNIT_ASSERT_EQUAL(4,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected3[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(RevIntegral);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected9[4]={2.48, 3.766666666666666, 4.220173, 1.9199999999999999};
  CPPUNIT_ASSERT_EQUAL(4,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected9[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->decrRef();
  // REVERSE ***********
  trgfield=MEDCouplingFieldDouble::New(ON_CELLS);
  trgfield->setNature(ConservativeVolumic);
  trgfield->setMesh(targetMesh);
  array=DataArrayDouble::New();
  array->alloc(targetMesh->getNumberOfCells(),1);
  trgfield->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  srcField=remapper.reverseTransferField(trgfield,4.220173);
  values=srcField->getArray()->getConstPointer();
  const double valuesExpected4[2]={7.9375, 8.9};
  CPPUNIT_ASSERT_EQUAL(2,srcField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,srcField->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected4[i0],values[i0],1e-12);
  srcField->decrRef();
  //
  trgfield->decrRef();
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
  // REVERSE ALL
  sourceMesh=build2DTargetMesh_3();
  targetMesh=MEDCouplingBasicsTest::build2DSourceMesh_1();
  //
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ConservativeVolumic);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected5[2]={7.9375, 8.9};
  CPPUNIT_ASSERT_EQUAL(2,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected5[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(IntegralGlobConstraint);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected6[4]={9.25, 15.75};
  CPPUNIT_ASSERT_EQUAL(2,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected6[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(Integral);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected7[2]={4.56, 4.3466666666666667};
  CPPUNIT_ASSERT_EQUAL(2,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected7[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(RevIntegral);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected10[2]={5.08, 3.56};
  CPPUNIT_ASSERT_EQUAL(2,trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected10[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->decrRef();
  // REVERSE ***********
  trgfield=MEDCouplingFieldDouble::New(ON_CELLS);
  trgfield->setNature(ConservativeVolumic);
  trgfield->setMesh(targetMesh);
  array=DataArrayDouble::New();
  array->alloc(targetMesh->getNumberOfCells(),1);
  trgfield->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<targetMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  srcField=remapper.reverseTransferField(trgfield,4.220173);
  values=srcField->getArray()->getConstPointer();
  const double valuesExpected8[4]={7.75, 7.0625,4.220173, 8.0};
  CPPUNIT_ASSERT_EQUAL(4,srcField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,srcField->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected8[i0],values[i0],1e-12);
  srcField->decrRef();
  //
  trgfield->decrRef();
  //
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingRemapperTest::testExtruded()
{
  MEDCouplingUMesh *mesh2DS=0;
  MEDCouplingUMesh *mesh3DS=build3DExtrudedUMesh_1(mesh2DS);
  MEDCouplingExtrudedMesh *extS=MEDCouplingExtrudedMesh::New(mesh3DS,mesh2DS,1);
  mesh3DS->decrRef();
  mesh2DS->decrRef();
  MEDCouplingUMesh *mesh2DT=0;
  MEDCouplingUMesh *mesh3DT=build3DExtrudedUMesh_1(mesh2DT);
  MEDCouplingExtrudedMesh *extT=MEDCouplingExtrudedMesh::New(mesh3DT,mesh2DT,1);
  //
  //
  mesh3DT->decrRef();
  mesh2DT->decrRef();
  //
  extS->decrRef();
  extT->decrRef();
}

MEDCouplingUMesh *MEDCouplingRemapperTest::build1DTargetMesh_2()
{
  double targetCoords[20]={
    0.59,0.09, 0.69,0.19, 0.21,-0.29,0.31,-0.19, 0.45,0.25,0.65,0.45,
    -0.2,-0.2,0.11,0.11, 0.25,0.25, 0.45,0.45
  };
  int targetConn[10]={0,1, 2,3, 4,5, 6,7, 8,9};

  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New("my name of mesh 1D 2",1);
  targetMesh->allocateCells(5);
  for(int i=0;i<5;i++)
    targetMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,targetConn+2*i);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(10,2);
  std::copy(targetCoords,targetCoords+20,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingRemapperTest::build2DTargetMesh_3()
{
  double targetCoords[20]={-0.6,-0.4, -0.1,-0.4, 1.1,-0.4, 2.1,-0.4,
                           -0.6,0.1,  -0.1,0.1,  1.1,0.1,  2.1,0.1,
                           -0.6,1.1,  -0.1,1.1};
  int targetConn[16]={0,4,5,1, 1,5,6,2, 2,6,7,3, 4,8,9,5};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(4);
  for(int i=0;i<4;i++)
    targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+4*i);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(10,2);
  std::copy(targetCoords,targetCoords+20,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDCouplingRemapperTest::build3DExtrudedUMesh_1(MEDCouplingUMesh *&mesh2D)
{
  double coords[180]={
    0.,0.,0., 1.,1.,0., 1.,1.25,0., 0.,1.,0., 1.,1.5,0., 2.,0.,0., 2.,1.,0., 1.,2.,0., 0.,2.,0., 3.,1.,0.,
    3.,2.,0., 0.,1.,0., 1.,3.,0., 2.,2.,0., 2.,3.,0.,
    0.,0.,1., 1.,1.,1., 1.,1.25,1., 0.,1.,1., 1.,1.5,1., 2.,0.,1., 2.,1.,1., 1.,2.,1., 0.,2.,1., 3.,1.,1.,
    3.,2.,1., 0.,1.,1., 1.,3.,1., 2.,2.,1., 2.,3.,1.,
    0.,0.,2., 1.,1.,2., 1.,1.25,2., 0.,1.,2., 1.,1.5,2., 2.,0.,2., 2.,1.,2., 1.,2.,2., 0.,2.,2., 3.,1.,2.,
    3.,2.,2., 0.,1.,2., 1.,3.,2., 2.,2.,2., 2.,3.,2.,
    0.,0.,3., 1.,1.,3., 1.,1.25,3., 0.,1.,3., 1.,1.5,3., 2.,0.,3., 2.,1.,3., 1.,2.,3., 0.,2.,3., 3.,1.,3.,
    3.,2.,3., 0.,1.,3., 1.,3.,3., 2.,2.,3., 2.,3.,3.};

  int conn[354]={
    // 0
    0,11,1,3,15,26,16,18,   1,2,4,7,13,6,-1,1,16,21,6,-1,6,21,28,13,-1,13,7,22,28,-1,7,4,19,22,-1,4,2,17,19,-1,2,1,16,17,-1,16,21,28,22,19,17,
    1,6,5,3,16,21,20,18,   13,10,9,6,28,25,24,21,
    11,8,7,4,2,1,-1,11,26,16,1,-1,1,16,17,2,-1,2,17,19,4,-1,4,19,22,7,-1,7,8,23,22,-1,8,11,26,23,-1,26,16,17,19,22,23,
    7,12,14,13,22,27,29,28,
    // 1
    15,26,16,18,30,41,31,33,   16,17,19,22,28,21,-1,16,31,36,21,-1,21,36,43,28,-1,28,22,37,43,-1,22,19,34,37,-1,19,17,32,34,-1,17,16,31,32,-1,31,36,43,37,34,32,
    16,21,20,18,31,36,35,33,   28,25,24,21,43,40,39,36,
    26,23,22,19,17,16,-1,26,41,31,16,-1,16,31,32,17,-1,17,32,34,19,-1,19,34,37,22,-1,22,23,38,37,-1,23,26,41,38,-1,41,31,32,34,37,38,
    22,27,29,28,37,42,44,43,
    // 2
    30,41,31,33,45,56,46,48,  31,32,34,37,43,36,-1,31,46,51,36,-1,36,51,58,43,-1,43,37,52,58,-1,37,34,49,52,-1,34,32,47,49,-1,32,31,46,47,-1,46,51,58,52,49,47,
    31,36,35,33,46,51,50,48,  43,40,39,36,58,55,54,51,
    41,38,37,34,32,31,-1,41,56,46,31,-1,31,46,47,32,-1,32,47,49,34,-1,34,49,52,37,-1,37,38,53,52,-1,38,41,56,53,-1,56,46,47,49,52,53,
    37,42,44,43,52,57,59,58
  };
  int conn2[28]={7,12,14,13, 11,8,7,4,2,1, 13,10,9,6, 1,6,5,3, 1,2,4,7,13,6, 0,11,1,3};
  //
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  ret->setMeshDimension(3);
  ret->allocateCells(18);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+8);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+51);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+59);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+67);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+110);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+118);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+126);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+169);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+177);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+185);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+228);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+236);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+244);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+287);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+295);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+303);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+346);
  //
  ret->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(60,3);
  std::copy(coords,coords+180,myCoords->getPointer());
  ret->setCoords(myCoords);
  //
  mesh2D=MEDCouplingUMesh::New();
  mesh2D->setMeshDimension(2);
  mesh2D->allocateCells(6);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_POLYGON,6,conn2+4);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+10);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+14);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_POLYGON,6,conn2+18);
  mesh2D->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+24);
  mesh2D->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}
