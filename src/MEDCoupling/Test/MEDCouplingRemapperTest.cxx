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

#include "MEDCouplingRemapperTest.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingRemapper.hxx"

#include "MEDCouplingBasicsTest.hxx"

#include <cmath>
#include <numeric>

using namespace MEDCoupling;

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
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(5,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(ExtensiveConservation);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected2[5]={3.75 ,1.75 ,1.75,4.,3.75};
  CPPUNIT_ASSERT_EQUAL(5,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(IntensiveMaximum);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(5,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(ExtensiveConservation);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(5,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(ExtensiveMaximum);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(5,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(IntensiveConservation);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(5,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
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
  targetField->setNature(IntensiveMaximum);
  targetField->setMesh(targetMesh);
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
  CPPUNIT_ASSERT_EQUAL(2,(int)srcfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)srcfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  srcfield->decrRef();
  //
  targetField->setNature(ExtensiveConservation);
  srcfield=remapper.reverseTransferField(targetField,4.57);
  values=srcfield->getArray()->getConstPointer();
  const double valuesExpected2[2]={26., 19.};
  CPPUNIT_ASSERT_EQUAL(2,(int)srcfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)srcfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  srcfield->decrRef();
  //
  targetField->decrRef();
  //clean up
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingRemapperTest::test1DInterp_1()
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingBasicsTest::build1DSourceMesh_2();
  MEDCouplingUMesh *targetMesh=MEDCouplingBasicsTest::build1DTargetMesh_2();
  //
  MEDCouplingRemapper remapper;
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  MEDCouplingFieldDouble *srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(sourceMesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  double *ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  //
  MEDCouplingFieldDouble *trgfield=remapper.transferField(srcField,4.57);
  const double *values=trgfield->getArray()->getConstPointer();
  const double valuesExpected1[2]={9.0540540540540526,7.4};
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected1[i0],values[i0],1e-12);
  trgfield->decrRef();
  const double valuesExpected2[2]={24.75,5.75};
  srcField->setNature(ExtensiveMaximum);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  const double valuesExpected3[2]={24.75,9.25};
  srcField->setNature(ExtensiveConservation);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected3[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  const double valuesExpected4[2]={7.4444444444444446,7.4};
  srcField->setNature(IntensiveConservation);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected4[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  //2D Curve
  sourceMesh=MEDCouplingBasicsTest::build2DCurveSourceMesh_2();
  targetMesh=MEDCouplingBasicsTest::build2DCurveTargetMesh_2();
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  //
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected1[i0],values[i0],1e-12);
  trgfield->decrRef();
  srcField->setNature(ExtensiveMaximum);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(ExtensiveConservation);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected3[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(IntensiveConservation);
  trgfield=remapper.transferField(srcField,4.57);
  values=trgfield->getArray()->getConstPointer();
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected4[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->decrRef();
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
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(5,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  srcField->decrRef();
  //
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P1P0"));
  srcField=MEDCouplingFieldDouble::New(ON_NODES);
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(5,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  srcField->decrRef();
  //
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(targetMesh,sourceMesh,"P0P1"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(4,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
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
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(9,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
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
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(5,(int)trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getArray()->getNumberOfComponents());
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
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+8);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  CPPUNIT_ASSERT_EQUAL(5,(int)trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getArray()->getNumberOfComponents());
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
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  CPPUNIT_ASSERT_EQUAL(8,(int)trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getArray()->getNumberOfComponents());
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
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  CPPUNIT_ASSERT_EQUAL(8,(int)trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getArray()->getNumberOfComponents());
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
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  CPPUNIT_ASSERT_EQUAL(8,(int)trgField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getArray()->getNumberOfComponents());
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
  srcField->setNature(IntensiveMaximum);
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
  srcField->setNature(IntensiveMaximum);
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
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getNumberOfComponents());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,values[0],1e-14);
  srcField->decrRef();
  srcField=remapper.reverseTransferField(trgField,4.220173);
  CPPUNIT_ASSERT_EQUAL(5,(int)srcField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)srcField->getNumberOfComponents());
  values=srcField->getArray()->getConstPointer();
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,values[i0],1e-14);
  srcField->decrRef();
  trgField->setNature(ExtensiveMaximum);
  srcField=remapper.reverseTransferField(trgField,4.220173);
  CPPUNIT_ASSERT_EQUAL(5,(int)srcField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)srcField->getNumberOfComponents());
  values=srcField->getArray()->getConstPointer();
  const double valuesExpected6[5]={2.28125,1.140625,1.140625,2.28125,2.28125};
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected6[i0],values[i0],1e-14);
  srcField->decrRef();
  trgField->decrRef();
  // ------------- -1D -> 2D
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(targetMesh,sourceMesh,"P0P0"));
  trgField=MEDCouplingFieldDouble::New(ON_CELLS);
  trgField->setNature(IntensiveMaximum);
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
  trgField->setNature(ExtensiveConservation);
  srcField=remapper.transferField(trgField,4.221073);
  values=srcField->getArray()->getConstPointer();
  const double valuesExpected7[5]={1.75,0.875,0.875,1.75,1.75};
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected7[i0],values[i0],1e-14);
  srcField->decrRef();
  trgField->setNature(ExtensiveMaximum);
  srcField=remapper.transferField(trgField,4.221073);
  values=srcField->getArray()->getConstPointer();
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected7[i0],values[i0],1e-14);
  srcField->decrRef();
  trgField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  //------------- 1D -> 2D
  const int conn[8]={0,1,1,2,2,3,3,0};
  const int conn2[12]={6,7,5,4,2,7,6,3,0,4,5,1};
  const double coords1[]={0.17,0.93,0.56,0.93,0.56,0.25,0.17,0.52};
  const double coords2[]={0.,0.,1.,0.,1.,1.,0.,1.,0.,0.5,1.,0.5,0.,0.8,1.,0.8};
  sourceMesh=MEDCouplingUMesh::New("src1D",1);
  sourceMesh->allocateCells(4);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+4);
  sourceMesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+6);
  sourceMesh->finishInsertingCells();
  array=DataArrayDouble::New(); array->alloc(4,2);
  std::copy(coords1,coords1+8,array->getPointer());
  sourceMesh->setCoords(array); array->decrRef();
  targetMesh=MEDCouplingUMesh::New("trg2D",2);
  targetMesh->allocateCells(3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+4);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn2+8);
  targetMesh->finishInsertingCells();
  array=DataArrayDouble::New(); array->alloc(8,2);
  std::copy(coords2,coords2+16,array->getPointer());
  targetMesh->setCoords(array); array->decrRef();
  remapper.setPrecision(1e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Geometric2D);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(4,1); array->iota(2.);
  srcField->setArray(array); array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  const double valuesExpected10[3]={3.9674868868103834, 2.8, 3.6372633449255796};
  CPPUNIT_ASSERT_EQUAL(3,(int)trgField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getNumberOfComponents());
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected10[i],trgField->getIJ(i,0),1e-13);
  srcField->decrRef();
  trgField->decrRef();
  //------------- 2D -> 1D
  remapper.setPrecision(1e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Geometric2D);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(targetMesh,sourceMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(targetMesh);
  array=DataArrayDouble::New();
  array->alloc(3,1); array->iota(2.);
  srcField->setArray(array); array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  const double valuesExpected11[4]={3., 2.9264705882352944, 3.8518518518518516, 2.3170731707317076};
  CPPUNIT_ASSERT_EQUAL(4,(int)trgField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getNumberOfComponents());
  for(int i=0;i<4;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected11[i],trgField->getIJ(i,0),1e-13);
  srcField->decrRef();
  trgField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
  //------------- 2D -> 3D
  sourceMesh=MEDCouplingBasicsTest::build3D2DSourceMesh();
  targetMesh=MEDCouplingBasicsTest::build3D2DTargetMesh();
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(sourceMesh);
  array=DataArrayDouble::New();
  array->alloc(7,1); array->iota(2.);
  srcField->setArray(array); array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  const double valuesExpected12[3]={5.70909090909091, 6.08362715128042, 6.92857142857143};
  CPPUNIT_ASSERT_EQUAL(3,(int)trgField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getNumberOfComponents());
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected12[i],trgField->getIJ(i,0),1e-13);
  srcField->decrRef();
  trgField->decrRef();
  //------------- 3D -> 2D
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(targetMesh,sourceMesh,"P0P0"));
  srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(targetMesh);
  array=DataArrayDouble::New();
  array->alloc(3,1); array->iota(2.);
  srcField->setArray(array); array->decrRef();
  trgField=remapper.transferField(srcField,4.57);
  const double valuesExpected13[7]={3., 4., 2.5, 2.909090909090909, 2., 3.5, 3.3571428571428572};
  CPPUNIT_ASSERT_EQUAL(7,(int)trgField->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgField->getNumberOfComponents());
  for(int i=0;i<7;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected13[i],trgField->getIJ(i,0),1e-13);
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
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(4,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(ExtensiveConservation);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected2[4]={2.8374999999999999, 7.3624999999999998, 4.220173, 4.7999999999999998};
  CPPUNIT_ASSERT_EQUAL(4,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(ExtensiveMaximum);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected3[4]={1.24, 4.5199999999999996, 4.220173, 1.9199999999999999};
  CPPUNIT_ASSERT_EQUAL(4,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected3[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(IntensiveConservation);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected9[4]={2.48, 3.766666666666666, 4.220173, 1.9199999999999999};
  CPPUNIT_ASSERT_EQUAL(4,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected9[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->decrRef();
  // REVERSE ***********
  trgfield=MEDCouplingFieldDouble::New(ON_CELLS);
  trgfield->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(2,(int)srcField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)srcField->getArray()->getNumberOfComponents());
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
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected5[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(ExtensiveConservation);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected6[4]={9.25, 15.75};
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected6[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(ExtensiveMaximum);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected7[2]={4.56, 4.3466666666666667};
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected7[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->setNature(IntensiveConservation);
  trgfield=remapper.transferField(srcField,4.220173);
  values=trgfield->getArray()->getConstPointer();
  const double valuesExpected10[2]={5.08, 3.56};
  CPPUNIT_ASSERT_EQUAL(2,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<2;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected10[i0],values[i0],1e-12);
  trgfield->decrRef();
  //
  srcField->decrRef();
  // REVERSE ***********
  trgfield=MEDCouplingFieldDouble::New(ON_CELLS);
  trgfield->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(4,(int)srcField->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)srcField->getArray()->getNumberOfComponents());
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
  MEDCouplingMappedExtrudedMesh *extS=MEDCouplingMappedExtrudedMesh::New(mesh3DS,mesh2DS,1);
  mesh3DS->decrRef();
  mesh2DS->decrRef();
  MEDCouplingUMesh *mesh2DT=0;
  MEDCouplingUMesh *mesh3DT=build3DExtrudedUMesh_1(mesh2DT);
  MEDCouplingMappedExtrudedMesh *extT=MEDCouplingMappedExtrudedMesh::New(mesh3DT,mesh2DT,1);
  //
  //
  mesh3DT->decrRef();
  mesh2DT->decrRef();
  //
  extS->decrRef();
  extT->decrRef();
}

void MEDCouplingRemapperTest::testExtruded2()
{
  MEDCouplingUMesh *meshN,*meshTT,*meshTF;
  MEDCouplingBasicsTest::build3DExtrudedUMesh_2(meshN,meshTT,meshTF);
  std::vector<int> n;
  double pt[3]={300.,300.,0.};
  double v[3]={0.,0.,2.};
  meshN->findNodesOnPlane(pt,v,1e-12,n);
  MEDCouplingUMesh *meshN2D=(MEDCouplingUMesh *)meshN->buildFacePartOfMySelfNode(&n[0],&n[0]+n.size(),true);
  n.clear();
  bool b=false;
  int newNbOfNodes;
  DataArrayInt *da=meshTT->mergeNodes(1e-12,b,newNbOfNodes);
  CPPUNIT_ASSERT(b);
  da->decrRef();
  meshTT->findNodesOnPlane(pt,v,1e-12,n);
  MEDCouplingUMesh *meshTT2D=(MEDCouplingUMesh *)meshTT->buildFacePartOfMySelfNode(&n[0],&n[0]+n.size(),true);
  n.clear();
  meshTF->findNodesOnPlane(pt,v,1e-12,n);
  MEDCouplingUMesh *meshTF2D=(MEDCouplingUMesh *)meshTF->buildFacePartOfMySelfNode(&n[0],&n[0]+n.size(),true);
  n.clear();
  //
  MEDCouplingMappedExtrudedMesh *meshNE=MEDCouplingMappedExtrudedMesh::New(meshN,meshN2D,0);
  MEDCouplingMappedExtrudedMesh *meshTTE=MEDCouplingMappedExtrudedMesh::New(meshTT,meshTT2D,0);
  MEDCouplingMappedExtrudedMesh *meshTFE=MEDCouplingMappedExtrudedMesh::New(meshTF,meshTF2D,0);
  //
  MEDCouplingRemapper remapper;
  remapper.setPrecision(1e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(meshNE,meshTTE,"P0P0"));
  MEDCouplingFieldDouble *srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(ExtensiveConservation);
  srcField->setMesh(meshNE);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(meshNE->getNumberOfCells(),1);
  srcField->setArray(array);
  double vals1[40]={
    1000.,1000.,1020.,1030.,1040.,1000.,1000.,1070.,1080.,1090.,1000.,1000.,1120.,1130.,1140.,1000.,1000.,1170.,1180.,1190.,
    2000.,2000.,2020.,2030.,2040.,2000.,2000.,2070.,2080.,2090.,2000.,2000.,2120.,2130.,2140.,2000.,2000.,2170.,2180.,2190.,
  };
  CPPUNIT_ASSERT_EQUAL((int)(sizeof(vals1)/sizeof(double)),(int)meshNE->getNumberOfCells());
  std::copy(vals1,vals1+meshNE->getNumberOfCells(),array->getPointer());
  array->decrRef();
  MEDCouplingFieldDouble *trgField=remapper.transferField(srcField,4.220173);
  double expected1[200]={
    800.,800.,800.,800.,800.,800.,800.,800.,800.,800.,1600.,1600.,1600.,1600.,1600.,1600.,1600.,1600.,1600.,1600.,
    102.,102.,102.,102.,102.,102.,102.,102.,102.,102.,202.,202.,202.,202.,202.,202.,202.,202.,202.,202.,
    103.,103.,103.,103.,103.,103.,103.,103.,103.,103.,203.,203.,203.,203.,203.,203.,203.,203.,203.,203.,
    104.,104.,104.,104.,104.,104.,104.,104.,104.,104.,204.,204.,204.,204.,204.,204.,204.,204.,204.,204.,
    219.,219.,219.,219.,219.,219.,219.,219.,219.,219.,419.,419.,419.,419.,419.,419.,419.,419.,419.,419.,
    221.,221.,221.,221.,221.,221.,221.,221.,221.,221.,421.,421.,421.,421.,421.,421.,421.,421.,421.,421.,
    223.,223.,223.,223.,223.,223.,223.,223.,223.,223.,423.,423.,423.,423.,423.,423.,423.,423.,423.,423.,
    117.,117.,117.,117.,117.,117.,117.,117.,117.,117.,217.,217.,217.,217.,217.,217.,217.,217.,217.,217.,
    118.,118.,118.,118.,118.,118.,118.,118.,118.,118.,218.,218.,218.,218.,218.,218.,218.,218.,218.,218.,
    119.,119.,119.,119.,119.,119.,119.,119.,119.,119.,219.,219.,219.,219.,219.,219.,219.,219.,219.,219.
    };
  for(int i=0;i<200;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],trgField->getArray()->getConstPointer()[i],1e-3);//1e-3 precision due to non coincidence in 1D mesh
  CPPUNIT_ASSERT_DOUBLES_EQUAL(std::accumulate(expected1,expected1+200,0.),std::accumulate(vals1,vals1+40,0.),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(std::accumulate(expected1,expected1+200,0.),std::accumulate(trgField->getArray()->getConstPointer(),trgField->getArray()->getConstPointer()+200,0.),1e-10);
  trgField->decrRef();
  //
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(meshNE,meshTFE,"P0P0"));
  trgField=remapper.transferField(srcField,4.220173);
  double expected2[340]={25.5, 51.25, 51.75, 26., 79., 158.75, 160.25, 80.5, 85.25, 171.25, 172.75, 86.75, 29.25, 58.75, 59.25, 29.75, 25.5, 51.25, 51.75, 26., 79., 158.75,
                         160.25, 80.5, 85.25, 171.25, 172.75, 86.75, 29.25, 58.75, 59.25, 29.75, 25.5, 51.25, 51.75, 26., 79., 158.75, 160.25, 80.5, 85.25, 171.25, 172.75, 86.75,
                         29.25, 58.75, 59.25, 29.75, 25.5, 51.25, 51.75, 26., 79., 158.75, 160.25, 80.5, 85.25, 171.25, 172.75, 86.75, 29.25, 58.75, 59.25, 29.75, 25.5, 51.25, 51.75,
                         26., 79., 158.75, 160.25, 80.5, 85.25, 171.25, 172.75, 86.75, 29.25, 58.75, 59.25, 29.75, 25.5, 51.25, 51.75, 26., 79., 158.75, 160.25, 80.5, 85.25, 171.25,
                         172.75, 86.75, 29.25, 58.75, 59.25, 29.75, 25.5, 51.25, 51.75, 26., 79., 158.75, 160.25, 80.5, 85.25, 171.25, 172.75, 86.75, 29.25, 58.75, 59.25, 29.75, 25.5,
                         51.25, 51.75, 26., 79., 158.75, 160.25, 80.5, 85.25, 171.25, 172.75, 86.75, 29.25, 58.75, 59.25, 29.75, 25.5, 51.25, 51.75, 26., 79., 158.75, 160.25, 80.5,
                         85.25, 171.25, 172.75, 86.75, 29.25, 58.75, 59.25, 29.75, 25.5, 51.25, 51.75, 26., 79., 158.75, 160.25, 80.5, 85.25, 171.25, 172.75, 86.75, 29.25, 58.75, 59.25,
                         29.75, 50.5, 101.25, 101.75, 51., 154., 308.75, 310.25, 155.5, 160.25, 321.25, 322.75, 161.75, 54.25, 108.75, 109.25, 54.75, 50.5, 101.25, 101.75, 51., 154.,
                         308.75, 310.25, 155.5, 160.25, 321.25, 322.75, 161.75, 54.25, 108.75, 109.25, 54.75, 50.5, 101.25, 101.75, 51., 154., 308.75, 310.25, 155.5, 160.25, 321.25, 322.75,
                         161.75, 54.25, 108.75, 109.25, 54.75, 50.5, 101.25, 101.75, 51., 154., 308.75, 310.25, 155.5, 160.25, 321.25, 322.75, 161.75, 54.25, 108.75, 109.25, 54.75, 50.5,
                         101.25, 101.75, 51., 154., 308.75, 310.25, 155.5, 160.25, 321.25, 322.75, 161.75, 54.25, 108.75, 109.25, 54.75, 50.5, 101.25, 101.75, 51., 154., 308.75, 310.25,
                         155.5, 160.25, 321.25, 322.75, 161.75, 54.25, 108.75, 109.25, 54.75, 50.5, 101.25, 101.75, 51., 154., 308.75, 310.25, 155.5, 160.25, 321.25, 322.75, 161.75, 54.25,
                         108.75, 109.25, 54.75, 50.5, 101.25, 101.75, 51., 154., 308.75, 310.25, 155.5, 160.25, 321.25, 322.75, 161.75, 54.25, 108.75, 109.25, 54.75, 50.5, 101.25, 101.75,
                         51., 154., 308.75, 310.25, 155.5, 160.25, 321.25, 322.75, 161.75, 54.25, 108.75, 109.25, 54.75, 50.5, 101.25, 101.75, 51., 154., 308.75, 310.25, 155.5, 160.25, 321.25,
                         322.75, 161.75, 54.25, 108.75, 109.25, 54.75, 800., 800., 800., 800., 800., 800., 800., 800., 800., 800., 1600., 1600., 1600., 1600., 1600., 1600., 1600.,
                         1600., 1600., 1600.};
  for(int i=0;i<340;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],trgField->getArray()->getConstPointer()[i],1e-3);//1e-3 precision due to non coincidence in 1D mesh
  CPPUNIT_ASSERT_DOUBLES_EQUAL(std::accumulate(expected2,expected2+340,0.),std::accumulate(vals1,vals1+40,0.),1e-10);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(std::accumulate(expected2,expected2+340,0.),std::accumulate(trgField->getArray()->getConstPointer(),trgField->getArray()->getConstPointer()+340,0.),1e-10);
  trgField->decrRef();
  srcField->decrRef();
  //
  double vals2[200]={
    100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1200., 1300., 1400., 1500., 1600., 1700., 1800., 1900., 2000,
    101., 201., 301., 401., 501., 601., 701., 801., 901., 1001., 1101., 1201., 1301., 1401., 1501., 1601., 1701., 1801., 1901., 2001,
    102., 202., 302., 402., 502., 602., 702., 802., 902., 1002., 1102., 1202., 1302., 1402., 1502., 1602., 1702., 1802., 1902., 2002,
    103., 203., 303., 403., 503., 603., 703., 803., 903., 1003., 1103., 1203., 1303., 1403., 1503., 1603., 1703., 1803., 1903., 2003,
    104., 204., 304., 404., 504., 604., 704., 804., 904., 1004., 1104., 1204., 1304., 1404., 1504., 1604., 1704., 1804., 1904., 2004,
    105., 205., 305., 405., 505., 605., 705., 805., 905., 1005., 1105., 1205., 1305., 1405., 1505., 1605., 1705., 1805., 1905., 2005,
    106., 206., 306., 406., 506., 606., 706., 806., 906., 1006., 1106., 1206., 1306., 1406., 1506., 1606., 1706., 1806., 1906., 2006,
    107., 207., 307., 407., 507., 607., 707., 807., 907., 1007., 1107., 1207., 1307., 1407., 1507., 1607., 1707., 1807., 1907., 2007,
    108., 208., 308., 408., 508., 608., 708., 808., 908., 1008., 1108., 1208., 1308., 1408., 1508., 1608., 1708., 1808., 1908., 2008,
    109., 209., 309., 409., 509., 609., 709., 809., 909., 1009., 1109., 1209., 1309., 1409., 1509., 1609., 1709., 1809., 1909., 2009.
  };
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(meshNE,meshTTE,"P0P0"));
  trgField=MEDCouplingFieldDouble::New(ON_CELLS);
  trgField->setNature(IntensiveMaximum);
  trgField->setMesh(meshTTE);
  array=DataArrayDouble::New();
  array->alloc(meshTTE->getNumberOfCells(),1);
  trgField->setArray(array);
  std::copy(vals2,vals2+meshTTE->getNumberOfCells(),array->getPointer());
  array->decrRef();
  srcField=remapper.reverseTransferField(trgField,4.220173);
  double expected3[40]={
    550.,550.,551.,552.,553.,550.,550.,554.,555.,556.,550.,550.,554.,555.,556.,550.,550.,557.,558.,559.,
    1550.,1550.,1551.,1552.,1553.,1550.,1550.,1554.,1555.,1556.,1550.,1550.,1554.,1555.,1556.,1550.,1550.,1557.,1558.,1559.
  };
  for(int i=0;i<40;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected3[i],srcField->getArray()->getConstPointer()[i],1e-3);//1e-3 precision due to non coincidence in 1D mesh
  srcField->decrRef();
  trgField->decrRef();
  //
  double vals3[340]={
    100., 101., 102., 103., 104., 105., 106., 107., 108., 109., 110., 111., 112., 113., 114., 115.,
    200., 201., 202., 203., 204., 205., 206., 207., 208., 209., 210., 211., 212., 213., 214., 215.,
    300., 301., 302., 303., 304., 305., 306., 307., 308., 309., 310., 311., 312., 313., 314., 315.,
    400., 401., 402., 403., 404., 405., 406., 407., 408., 409., 410., 411., 412., 413., 414., 415.,
    500., 501., 502., 503., 504., 505., 506., 507., 508., 509., 510., 511., 512., 513., 514., 515.,
    600., 601., 602., 603., 604., 605., 606., 607., 608., 609., 610., 611., 612., 613., 614., 615.,
    700., 701., 702., 703., 704., 705., 706., 707., 708., 709., 710., 711., 712., 713., 714., 715.,
    800., 801., 802., 803., 804., 805., 806., 807., 808., 809., 810., 811., 812., 813., 814., 815.,
    900., 901., 902., 903., 904., 905., 906., 907., 908., 909., 910., 911., 912., 913., 914., 915.,
    1000., 1001., 1002., 1003., 1004., 1005., 1006., 1007., 1008., 1009., 1010., 1011., 1012., 1013., 1014., 1015.,
    1100., 1101., 1102., 1103., 1104., 1105., 1106., 1107., 1108., 1109., 1110., 1111., 1112., 1113., 1114., 1115.,
    1200., 1201., 1202., 1203., 1204., 1205., 1206., 1207., 1208., 1209., 1210., 1211., 1212., 1213., 1214., 1215.,
    1300., 1301., 1302., 1303., 1304., 1305., 1306., 1307., 1308., 1309., 1310., 1311., 1312., 1313., 1314., 1315.,
    1400., 1401., 1402., 1403., 1404., 1405., 1406., 1407., 1408., 1409., 1410., 1411., 1412., 1413., 1414., 1415.,
    1500., 1501., 1502., 1503., 1504., 1505., 1506., 1507., 1508., 1509., 1510., 1511., 1512., 1513., 1514., 1515.,
    1600., 1601., 1602., 1603., 1604., 1605., 1606., 1607., 1608., 1609., 1610., 1611., 1612., 1613., 1614., 1615.,
    1700., 1701., 1702., 1703., 1704., 1705., 1706., 1707., 1708., 1709., 1710., 1711., 1712., 1713., 1714., 1715.,
    1800., 1801., 1802., 1803., 1804., 1805., 1806., 1807., 1808., 1809., 1810., 1811., 1812., 1813., 1814., 1815.,
    1900., 1901., 1902., 1903., 1904., 1905., 1906., 1907., 1908., 1909., 1910., 1911., 1912., 1913., 1914., 1915.,
    2000., 2001., 2002., 2003., 2004., 2005., 2006., 2007., 2008., 2009., 2010., 2011., 2012., 2013., 2014., 2015.,
    116.,216.,316.,416.,516.,616.,716.,816.,916.,1016.,1116.,1216.,1316.,1416.,1516.,1616.,1716.,1816.,1916.,2016.
  };
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(meshNE,meshTFE,"P0P0"));
  trgField=MEDCouplingFieldDouble::New(ON_CELLS);
  trgField->setNature(IntensiveMaximum);
  trgField->setMesh(meshTFE);
  array=DataArrayDouble::New();
  array->alloc(meshTFE->getNumberOfCells(),1);
  trgField->setArray(array);
  std::copy(vals3,vals3+meshTFE->getNumberOfCells(),array->getPointer());
  array->decrRef();
  srcField=remapper.reverseTransferField(trgField,4.220173);
  double expected4[40]={
    566.,566.,552.5,553.5,554.5,566.,566.,554.5,555.5,556.5,566.,566.,558.5,559.5,560.5,566.,566.,560.5,561.5,562.5,
    1566.,1566.,1552.5,1553.5,1554.5,1566.,1566.,1554.5,1555.5,1556.5,1566.,1566.,1558.5,1559.5,1560.5,1566.,1566.,1560.5,1561.5,1562.5
  };
  for(int i=0;i<40;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected4[i],srcField->getArray()->getConstPointer()[i],1e-3);//1e-3 precision due to non coincidence in 1D mesh
  srcField->decrRef();
  trgField->decrRef();
  //
  meshN2D->decrRef();
  meshTT2D->decrRef();
  meshTF2D->decrRef();
  meshNE->decrRef();
  meshTTE->decrRef();
  meshTFE->decrRef();
  meshN->decrRef();
  meshTT->decrRef();
  meshTF->decrRef();
}

void MEDCouplingRemapperTest::testPrepareEx1()
{
  MEDCouplingUMesh *sourceMesh=MEDCouplingBasicsTest::build2DSourceMesh_1();
  MEDCouplingUMesh *targetMesh=build2DTargetMesh_3();
  //
  MEDCouplingRemapper remapper;
  remapper.setPrecision(1e-12);
  remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
  MEDCouplingFieldTemplate *srcFt=MEDCouplingFieldTemplate::New(ON_CELLS);
  MEDCouplingFieldTemplate *trgFt=MEDCouplingFieldTemplate::New(ON_CELLS);
  srcFt->setMesh(sourceMesh);
  trgFt->setMesh(targetMesh);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepareEx(srcFt,trgFt));
  srcFt->decrRef();
  trgFt->decrRef();
  MEDCouplingFieldDouble *srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(IntensiveMaximum);
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
  CPPUNIT_ASSERT_EQUAL(4,(int)trgfield->getArray()->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)trgfield->getArray()->getNumberOfComponents());
  for(int i0=0;i0<4;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected[i0],values[i0],1e-12);
  trgfield->decrRef();
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
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

void MEDCouplingRemapperTest::testPartialTransfer1()
{
  MEDCouplingRemapper remapper;
  MEDCouplingUMesh *sourceMesh=build1DTargetMesh_2();
  MEDCouplingUMesh *targetMesh=MEDCouplingBasicsTest::build2DTargetMesh_1();
  remapper.setIntersectionType(INTERP_KERNEL::PointLocator);
  CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(sourceMesh,targetMesh,"P0P0"));
  MEDCouplingFieldDouble *srcField=MEDCouplingFieldDouble::New(ON_CELLS);
  srcField->setNature(IntensiveMaximum);
  srcField->setMesh(sourceMesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(sourceMesh->getNumberOfCells(),1);
  srcField->setArray(array);
  double *ptr=array->getPointer();
  for(int i=0;i<sourceMesh->getNumberOfCells();i++)
    ptr[i]=(double)(i+7);
  array->decrRef();
  MEDCouplingFieldDouble *trgField=MEDCouplingFieldDouble::New(ON_CELLS);
  trgField->setNature(IntensiveMaximum);
  trgField->setMesh(targetMesh);
  array=DataArrayDouble::New();
  array->alloc(targetMesh->getNumberOfCells(),1);
  ptr=array->getPointer();
  std::fill(ptr,ptr+targetMesh->getNumberOfCells(),96.3);
  trgField->setArray(array);
  array->decrRef();
  remapper.partialTransfer(srcField,trgField);
  const double valuesExpected9[5]={10.,8.,7.,96.3,10.};
  const double *values=trgField->getArray()->getConstPointer();
  for(int i0=0;i0<5;i0++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected9[i0],values[i0],1e-12);
  trgField->decrRef();
  srcField->decrRef();
  sourceMesh->decrRef();
  targetMesh->decrRef();
}

void MEDCouplingRemapperTest::testBugNonRegression1()
{
  // source
  DataArrayDouble *coordsSrc(DataArrayDouble::New());
  const double coordsSrcData[18]={-6.25,3.6084391824351605,264.85199999999998,-6.25,3.6084391824351605,289.05200000000002,-6.2499999999999991,-3.6084391824351618,264.85199999999998,-6.2499999999999991,-3.6084391824351618,289.05200000000002,-1.7763568394002505e-15,4.4408920985006262e-15,264.85199999999998,-1.7763568394002505e-15,4.4408920985006262e-15,289.05200000000002};
  coordsSrc->useArray(coordsSrcData,false,CPP_DEALLOC,6,3);
  DataArrayInt *connSrc(DataArrayInt::New()),*connISrc(DataArrayInt::New());
  const int connSrcData[7]={16,2,0,4,3,1,5};
  connSrc->useArray(connSrcData,false,CPP_DEALLOC,7,1);
  const int connISrcData[2]={0,7};
  connISrc->useArray(connISrcData,false,CPP_DEALLOC,2,1);
  MEDCouplingUMesh *srcMesh(MEDCouplingUMesh::New("source",3));
  srcMesh->setCoords(coordsSrc);
  srcMesh->setConnectivity(connSrc,connISrc,true);
  coordsSrc->decrRef(); connSrc->decrRef(); connISrc->decrRef();
  // target
  DataArrayDouble *coordsTrg(DataArrayDouble::New());
const double coordsTrgData[36]={-2,1.1547005383792521,264.85199999999998,-2,0.57735026918962618,264.85199999999998,-2.5,0.2886751345948132,264.85199999999998,-2.5,1.443375672974065,264.85199999999998,-3.0000000000000004,1.1547005383792526,264.85199999999998,-3.0000000000000004,0.57735026918962662,264.85199999999998,-2,1.1547005383792521,289.05200000000002,-2,0.57735026918962618,289.05200000000002,-2.5,0.2886751345948132,289.05200000000002,-2.5,1.443375672974065,289.05200000000002,-3.0000000000000004,1.1547005383792526,289.05200000000002,-3.0000000000000004,0.57735026918962662,289.05200000000002};
 coordsTrg->useArray(coordsTrgData,false,CPP_DEALLOC,12,3);
 DataArrayInt *connTrg=DataArrayInt::New();
 const int connTrgData[44]={31,0,1,2,5,4,3,-1,7,6,9,10,11,8,-1,3,9,6,0,-1,4,10,9,3,-1,5,11,10,4,-1,2,8,11,5,-1,1,7,8,2,-1,0,6,7,1};
 connTrg->useArray(connTrgData,false,CPP_DEALLOC,44,1);
 DataArrayInt *connITrg=DataArrayInt::New();
 const int connITrgData[2]={0,44};
 connITrg->useArray(connITrgData,false,CPP_DEALLOC,2,1);
 MEDCouplingUMesh *trgMesh=MEDCouplingUMesh::New("target",3);
 trgMesh->setCoords(coordsTrg);
 trgMesh->setConnectivity(connTrg,connITrg,true);
 coordsTrg->decrRef(); connTrg->decrRef(); connITrg->decrRef();
 // Go !
 const double valExpected(20.957814771583468);
 MEDCouplingRemapper remapper;
 remapper.setPrecision(1e-12);
 remapper.setIntersectionType(INTERP_KERNEL::Triangulation);
 CPPUNIT_ASSERT_EQUAL(1,remapper.prepare(srcMesh,trgMesh,"P0P0"));
 std::vector<std::map<int,double> > matrx(remapper.getCrudeMatrix());
 CPPUNIT_ASSERT_EQUAL(1,(int)matrx.size());
 CPPUNIT_ASSERT_EQUAL(1,(int)matrx[0].size());
 CPPUNIT_ASSERT_DOUBLES_EQUAL(valExpected,matrx[0][0],1e-13);
 //
 srcMesh->decrRef(); trgMesh->decrRef();
}

