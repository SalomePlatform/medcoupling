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
// Author : Anthony Geay (EDF R&D)

#include "MEDCouplingBasicsTest.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingMultiFields.hxx"


void CppExample_MEDCouplingFieldDouble_WriteVTK()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_WriteVTK_1]
  // mesh1
  const double coords[3] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh1 = MEDCouplingCMesh::New();
  mesh1->setCoords(coordsArr,coordsArr); // mesh becomes a 2D one

  // 3 fields (lying on the same mesh!)
  MCAuto<MEDCouplingFieldDouble> field1 =
    mesh1->getMeasureField( true );
  MCAuto<MEDCouplingFieldDouble> field2 =
    mesh1->buildOrthogonalField();
  MCAuto<MEDCouplingFieldDouble> field3 =
    mesh1->fillFromAnalytic( ON_CELLS, 1, "x");
  field2->setName( "Normal" ); //  name is necessary!
  field3->setName( "Barycenter" ); //  name is necessary!

  // WriteVTK
  const char fileName[] = "testExample_MEDCouplingFieldDouble_WriteVTK.vtk";
  std::vector<const MEDCouplingFieldDouble *> fs( 3 ); // field series
  fs[0] = field1;
  fs[1] = field2;
  fs[2] = field3;
  MEDCouplingFieldDouble::WriteVTK( fileName, fs );
  //! [CppSnippet_MEDCouplingFieldDouble_WriteVTK_1]
  remove(fileName);
}

void CppExample_MEDCouplingFieldDouble_MaxFields()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_MaxFields_1]
  const double vals1[4]   = {0.,2., 4.,6.}; // for field 1
  const double vals2[4]   = {2.,0., 6.,4.}; // for field 2
  const double valsMax[4] = {2.,2., 6.,6.}; // expected max field
  const double valsMin[4] = {0.,0., 4.,4.}; // expected min field
  // field 1
  MCAuto<DataArrayDouble> valsArr1 = DataArrayDouble::New();
  valsArr1->useExternalArrayWithRWAccess( vals1, 2,2 ); // 2 tuples per 2 components
  MCAuto<MEDCouplingFieldDouble> field1 = MEDCouplingFieldDouble::New( ON_NODES );
  field1->setArray( valsArr1 );
  // field 2
  MCAuto<DataArrayDouble> valsArr2 = DataArrayDouble::New();
  valsArr2->useExternalArrayWithRWAccess( vals2, 2,2 ); // 2 tuples per 2 components
  MCAuto<MEDCouplingFieldDouble> field2 = MEDCouplingFieldDouble::New( ON_NODES );
  field2->setArray( valsArr2 );
  // max field 
  MCAuto<MEDCouplingFieldDouble> fieldMax = MEDCouplingFieldDouble::MaxFields( field1, field2 );
  CPPUNIT_ASSERT( std::equal( valsMax, valsMax+4, fieldMax->getArray()->getConstPointer() )); // fieldMax == valsMax
  // min field 
  MCAuto<MEDCouplingFieldDouble> fieldMin = MEDCouplingFieldDouble::MinFields( field1, field2 );
  CPPUNIT_ASSERT( std::equal( valsMin, valsMin+4, fieldMin->getArray()->getConstPointer() )); // fieldMin == valsMin
  //! [CppSnippet_MEDCouplingFieldDouble_MaxFields_1]
}

void CppExample_MEDCouplingFieldDouble_MergeFields()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_MergeFields_1]
  // mesh1
  const double coords[3] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh1 = MEDCouplingCMesh::New();
  mesh1->setCoords(coordsArr); // mesh becomes a 1D
  // field1
  MCAuto<MEDCouplingFieldDouble> field1 =
    mesh1->fillFromAnalytic( ON_CELLS, 1, "x");

  // mesh2 and field2
  MCAuto<MEDCouplingFieldDouble> field2 =
    field1->cloneWithMesh( true );
  double vec[1] = { 5. };
  (const_cast<MEDCoupling::MEDCouplingMesh *>(field2->getMesh()))->translate(vec); // translate mesh2
  field2->applyFunc("x + 5"); // "translate" field2

  // concatenate field1 and field2
  MCAuto<MEDCouplingFieldDouble> field3 =
    MEDCouplingFieldDouble::MergeFields( field1, field2 );
  std::vector<const MEDCouplingFieldDouble *> fields( 2 );
  fields[0] = field1;
  fields[1] = field2;
  MCAuto<MEDCouplingFieldDouble> field4 =
    MEDCouplingFieldDouble::MergeFields( fields );
  //! [CppSnippet_MEDCouplingFieldDouble_MergeFields_1]
}

void CppExample_MEDCouplingFieldDouble_substractInPlaceDM()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_substractInPlaceDM_1]
  const double coords1[4] = {0.,1.,2.,3.};
  const double coords2[4] = {2.,1.,0.,3.}; //  #0 <==> #2
  // mesh 1
  MCAuto<MEDCouplingUMesh> mesh1 = MEDCouplingUMesh::New();
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords1, 4, 1 );
  mesh1->setCoords(coordsArr);
  mesh1->setMeshDimension(0);
  mesh1->allocateCells(0);
  mesh1->finishInsertingCells();
  // mesh 2
  MCAuto<MEDCouplingUMesh> mesh2 =
    (MEDCouplingUMesh*) mesh1->deepCopy();
  mesh2->getCoords()->useExternalArrayWithRWAccess( coords2, 4, 1 );
  //! [CppSnippet_MEDCouplingFieldDouble_substractInPlaceDM_1]
  //! [CppSnippet_MEDCouplingFieldDouble_substractInPlaceDM_2]
  MCAuto<MEDCouplingFieldDouble> field1 =
    mesh1->fillFromAnalytic( MEDCoupling::ON_NODES,1,"x"); // field1 values == coords1
  MCAuto<MEDCouplingFieldDouble> field2 =
    mesh2->fillFromAnalytic( MEDCoupling::ON_NODES,1,"x"); // field2 values == coords2
  const double levOfCheck = 10; // nodes can be permuted
  field1->substractInPlaceDM( field2, levOfCheck, 1e-13, 0 ); // values #0 and #2 must swap
  //! [CppSnippet_MEDCouplingFieldDouble_substractInPlaceDM_2]
  //! [CppSnippet_MEDCouplingFieldDouble_substractInPlaceDM_3]
  field2->applyFunc( 1, 0.0 ); // all field2 values == 0.0
  CPPUNIT_ASSERT( field1->isEqual( field2, 1e-13, 1e-13 )); // field1 == field2 == 0.0
  //! [CppSnippet_MEDCouplingFieldDouble_substractInPlaceDM_3]
}

void CppExample_MEDCouplingFieldDouble_changeUnderlyingMesh()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_changeUnderlyingMesh_1]
  const double coords1[4] = {0.,1.,2.,3.};
  const double coords2[4] = {2.,1.,0.,3.}; //  #0 <==> #2
  // mesh 1
  MCAuto<MEDCouplingUMesh> mesh1 = MEDCouplingUMesh::New();
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords1, 4, 1 );
  mesh1->setCoords(coordsArr);
  mesh1->setMeshDimension(0);
  mesh1->allocateCells(0);
  mesh1->finishInsertingCells();
  // mesh 2
  MCAuto<MEDCouplingUMesh> mesh2 =
    (MEDCouplingUMesh*) mesh1->deepCopy();
  mesh2->getCoords()->useExternalArrayWithRWAccess( coords2, 4, 1 );
  //! [CppSnippet_MEDCouplingFieldDouble_changeUnderlyingMesh_1]
  //! [CppSnippet_MEDCouplingFieldDouble_changeUnderlyingMesh_2]
  MCAuto<MEDCouplingFieldDouble> field =
    mesh1->fillFromAnalytic( MEDCoupling::ON_NODES,1,"x"); // field values == coords1
  const double levOfCheck = 10; // nodes can be permuted
  field->changeUnderlyingMesh( mesh2, levOfCheck, 1e-13, 0 ); // values #0 and #2 must swap
  CPPUNIT_ASSERT( std::equal( coords2, coords2+4, field->getArray()->getConstPointer() ));
  //! [CppSnippet_MEDCouplingFieldDouble_changeUnderlyingMesh_2]
}

void CppExample_MEDCouplingFieldDouble_applyFunc_same_nb_comp()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc_same_nb_comp_1]
  const double v[4] = {1.,2., 3.,4.};
  MCAuto<DataArrayDouble> array = DataArrayDouble::New();
  array->useExternalArrayWithRWAccess( v, 2, 2 ); // 2 tuples per 2 components
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setArray( array );
  const char func[] = "IVec * v + JVec * w*w + 10";
  field->applyFunc( 2, func );
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 2 ); // 2 components remains
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc_same_nb_comp_1]
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc_same_nb_comp_2]
  const double* v2 = field->getArray()->getConstPointer();
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v2[0], 10 + v[0], 13 );      // "10 + IVec * v"  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v2[1], 10 + v[1]*v[1], 13 ); // "10 + JVec * v*v"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v2[2], 10 + v[2], 13 );      // "10 + IVec * v"  
  CPPUNIT_ASSERT_DOUBLES_EQUAL( v2[3], 10 + v[3]*v[3], 13 ); // "10 + JVec * v*v"
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc_same_nb_comp_2]
}

void CppExample_MEDCouplingFieldDouble_applyFunc3()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc3_1]
  // create a 2D vector field
  const double values[4] = {1.,1., 2.,1.};
  MCAuto<DataArrayDouble> array = DataArrayDouble::New();
  array->useExternalArrayWithRWAccess( values, 2, 2 ); // 2 tuples per 2 components
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setArray( array );
  // transform the field to a 3D vector field
  const char func[] = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10";
  const char* varNames[2] = { "a", "b" }; // names used to refer to X and Y components
  std::vector<std::string> varNamesVec( varNames, varNames+2 );
  field->applyFuncNamedCompo( 3, varNamesVec, func ); // require 3 components 
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 ); // 3 components as required
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc3_1]
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc3_2]
  double vec1[3]; // vector #1
  field->getArray()->getTuple( 1, vec1 );
  const double a = values[2], b = values[3]; // initial components of the vector #1
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vec1[0], 10 + b, 13 ); // "10 + IVec * b"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vec1[1], 10 + a, 13 ); // "10 + JVec * a"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vec1[2], 10 + sqrt(a*a+b*b), 13 ); // "10 + KVec * sqrt( a*a + b*b )"
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc3_2]
}

void CppExample_MEDCouplingFieldDouble_applyFunc2()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc2_1]
  // create a 2D vector field
  const double values[4] = {1.,1., 2.,1.};
  MCAuto<DataArrayDouble> array = DataArrayDouble::New();
  array->useExternalArrayWithRWAccess( values, 2, 2 ); // 2 tuples per 2 components
  array->setInfoOnComponent(0,"a"); // name used to refer to X component within a function
  array->setInfoOnComponent(1,"b"); // name used to refer to Y component within a function
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setArray( array );
  // transform the field to a 3D vector field
  const char func[] = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10";
  field->applyFuncCompo( 3, func ); // require 3 components 
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 ); // 3 components as required
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc2_1]
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc2_2]
  double vec1[3]; // vector #1
  field->getArray()->getTuple( 1, vec1 );
  const double a = values[2], b = values[3]; // initial components of the vector #1
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vec1[0], 10 + b, 13 ); // "10 + IVec * b"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vec1[1], 10 + a, 13 ); // "10 + JVec * a"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vec1[2], 10 + sqrt(a*a+b*b), 13 ); // "10 + KVec * sqrt( a*a + b*b )"
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc2_2]
}

void CppExample_MEDCouplingFieldDouble_applyFunc()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc_1]
  // create a 2D vector field
  const double values[4] = {1.,1., 2.,1.};
  MCAuto<DataArrayDouble> array = DataArrayDouble::New();
  array->useExternalArrayWithRWAccess( values, 2, 2 ); // 2 tuples per 2 components
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setArray( array );
  // transform the field to a 3D vector field
  const char func[] = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10";
  field->applyFunc( 3, func ); // require 3 components 
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 ); // 3 components as required
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc_1]
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc_2]
  double vec1[3]; // vector #1
  field->getArray()->getTuple( 1, vec1 );
  const double a = values[2], b = values[3]; // initial components of the vector #1
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vec1[0], 10 + b, 13 ); // "10 + IVec * b"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vec1[1], 10 + a, 13 ); // "10 + JVec * a"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vec1[2], 10 + sqrt(a*a+b*b), 13 ); // "10 + KVec * sqrt( a*a + b*b )"
  //! [CppSnippet_MEDCouplingFieldDouble_applyFunc_2]
}

void CppExample_MEDCouplingFieldDouble_applyFunc_val()
{
  using namespace MEDCoupling;
  //! [Snippet_MEDCouplingFieldDouble_applyFunc_val_1]
  // mesh
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh = MEDCouplingCMesh::New();
  mesh->setCoords(coordsArr,coordsArr); // mesh becomes a 2D structured mesh
  // field
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setMesh( mesh );
  field->fillFromAnalytic(2,"IVec * x + JVec * y"); // 2 components
  //! [Snippet_MEDCouplingFieldDouble_applyFunc_val_1]
  //! [Snippet_MEDCouplingFieldDouble_applyFunc_val_2]
  const double newValue = 7.;
  field->applyFunc( 3, newValue ); // # 3 components are required
  CPPUNIT_ASSERT( field->getIJ(1,0) == newValue ); // a value is as expected
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 );
  CPPUNIT_ASSERT( field->getNumberOfTuples() == mesh->getNumberOfCells() );
  //! [Snippet_MEDCouplingFieldDouble_applyFunc_val_2]
}

void CppExample_MEDCouplingFieldDouble_fillFromAnalytic3()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic3_1]
  const double coords[4] = {0.,2.,4.,6.}; // 6. is not used
  MCAuto<DataArrayDouble> x = DataArrayDouble::New();
  x->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<DataArrayDouble> y = DataArrayDouble::New();
  y->useExternalArrayWithRWAccess( coords, 2, 1 );
  MCAuto<MEDCouplingCMesh> mesh=MEDCouplingCMesh::New();
  mesh->setCoords(x,y);
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic3_1]
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic3_2]
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setMesh( mesh );
  const char func[] = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10";
  const char* varNames[2] = { "a", "b" }; // names used to refer to X and Y coord components
  std::vector<std::string> varNamesVec( varNames, varNames+2 );
  field->fillFromAnalyticNamedCompo( 3, varNamesVec, func );
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic3_2]
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic3_3]
  double val1[3]; // a value (vector) of the cell #1
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 ); // 3 components in the field
  field->getArray()->getTuple( 1, val1 );
  //
  MCAuto<DataArrayDouble> bc =
    mesh->computeCellCenterOfMass(); // func is applied to barycenters of cells
  double bc1[2]; // coordinates of the second point
  bc->getTuple( 1, bc1 );
  //
  double dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] );  // "sqrt( a*a + b*b )"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[0], 10 + bc1[1], 13 ); // "10 + IVec * b"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[1], 10 + bc1[0], 13 ); // "10 + JVec * a"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[2], 10 + dist  , 13 ); // "10 + KVec * sqrt( a*a + b*b )"
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic3_3]
}

void CppExample_MEDCouplingFieldDouble_fillFromAnalytic2()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic2_1]
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> x = DataArrayDouble::New();
  x->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<DataArrayDouble> y = DataArrayDouble::New();
  y->useExternalArrayWithRWAccess( coords, 2, 1 );
  x->setInfoOnComponent(0,"a"); //  name used to refer to X coordinate within a function
  y->setInfoOnComponent(0,"b"); //  name used to refer to Y coordinate within a function
  MCAuto<MEDCouplingCMesh> mesh=MEDCouplingCMesh::New();
  mesh->setCoords(x,y);
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic2_1]
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic2_2]
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setMesh( mesh );
  const char func[] = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10";
  field->fillFromAnalytic( 3, func );
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic2_2]
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic2_3]
  double val1[3]; // a value (vector) of the cell #1
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 ); // 3 components in the field
  field->getArray()->getTuple( 1, val1 );
  //
  MCAuto<DataArrayDouble> bc =
    mesh->computeCellCenterOfMass(); // func is applied to barycenters of cells
  double bc1[2]; // coordinates of the second point
  bc->getTuple( 1, bc1 );
  //
  double dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] );  // "sqrt( a*a + b*b )"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[0], 10 + bc1[1], 13 ); // "10 + IVec * b"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[1], 10 + bc1[0], 13 ); // "10 + JVec * a"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[2], 10 + dist  , 13 ); // "10 + KVec * sqrt( a*a + b*b )"
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic2_3]
}

void CppExample_MEDCouplingFieldDouble_fillFromAnalytic()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic_1]
  const double coords[3] = {0.,2.,4};
  MCAuto<DataArrayDouble> x = DataArrayDouble::New();
  x->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<DataArrayDouble> y = DataArrayDouble::New();
  y->useExternalArrayWithRWAccess( coords, 2, 1 );
  MCAuto<MEDCouplingCMesh> mesh=MEDCouplingCMesh::New();
  mesh->setCoords(x,y);
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic_1]
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic_2]
  const char func[] = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10";
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setMesh( mesh );
  field->fillFromAnalytic( 3, func );
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic_2]
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic_3]
  double val1[3]; // a value (vector) of the cell #1
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 ); // 3 components in the field
  field->getArray()->getTuple( 1, val1 );
  //
  MCAuto<DataArrayDouble> bc =
    mesh->computeCellCenterOfMass(); // func is applied to barycenters of cells
  double bc1[2]; // coordinates of the second point
  bc->getTuple( 1, bc1 );
  //
  double dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] );  // "sqrt( a*a + b*b )"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[0], 10 + bc1[1], 13 ); // "10 + IVec * b"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[1], 10 + bc1[0], 13 ); // "10 + JVec * a"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[2], 10 + dist  , 13 ); // "10 + KVec * sqrt( a*a + b*b )"
  //! [CppSnippet_MEDCouplingFieldDouble_fillFromAnalytic_3]
}

//! [Snippet_MEDCouplingFieldDouble_fillFromAnalytic_c_func_0]
bool getNewValue(const double *pos, double *res)
{
  res[0] = pos[0];
  res[1] = pos[1];
  res[2] = sqrt( pos[0]*pos[0] + pos[1]*pos[1] );
  return true;
}
//! [Snippet_MEDCouplingFieldDouble_fillFromAnalytic_c_func_0]

void CppExample_MEDCouplingFieldDouble_fillFromAnalytic_c_func()
{
  using namespace MEDCoupling;
  //! [Snippet_MEDCouplingFieldDouble_fillFromAnalytic_c_func_1]
  // mesh
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh = MEDCouplingCMesh::New();
  mesh->setCoords(coordsArr,coordsArr); // mesh becomes a 2D structured mesh
  // field
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setMesh( mesh );
  field->fillFromAnalytic( 3, &getNewValue ); // 3 components are required
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 );
  CPPUNIT_ASSERT( field->getNumberOfTuples() == mesh->getNumberOfCells() );
  //! [Snippet_MEDCouplingFieldDouble_fillFromAnalytic_c_func_1]
}

void CppExample_MEDCouplingFieldDouble_applyFunc_c_func()
{
  using namespace MEDCoupling;
  //! [Snippet_MEDCouplingFieldDouble_applyFunc_c_func_1]
  // mesh
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh = MEDCouplingCMesh::New();
  mesh->setCoords(coordsArr,coordsArr); // mesh becomes a 2D structured mesh
  // field
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS );
  field->setMesh( mesh );
  MCAuto<DataArrayDouble> bc = mesh->computeCellCenterOfMass();
  field->setArray( bc ); // 2 components here as the mesh is 2D
  //! [Snippet_MEDCouplingFieldDouble_applyFunc_c_func_1]
  //! [Snippet_MEDCouplingFieldDouble_applyFunc_c_func_2]
  field->applyFunc( 3, &getNewValue ); // 3 components are required
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 );
  CPPUNIT_ASSERT( field->getNumberOfTuples() == mesh->getNumberOfCells() );
  //! [Snippet_MEDCouplingFieldDouble_applyFunc_c_func_2]
}

void CppExample_MEDCouplingFieldDouble_getValueOn_time()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_time_1]
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh = MEDCouplingCMesh::New();
  mesh->setCoords(coordsArr,coordsArr);
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_time_1]
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_time_2]
  MCAuto<MEDCouplingFieldDouble> field =
    MEDCouplingFieldDouble::New( MEDCoupling::ON_CELLS, MEDCoupling::LINEAR_TIME );
  field->setMesh( mesh );
  field->fillFromAnalytic( 1,"10"); // all values == 10.
  MCAuto<DataArrayDouble> array2 =
    DataArrayDouble::Add( field->getArray(), field->getArray() ); // == 2 * field->getArray()
  field->setEndArray( array2 ); // all values == 20.
  const double time1 = 1.1, time2 = 22.;
  field->setStartTime( time1, 0, 0 );
  field->setEndTime  ( time2, 0, 0 );
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_time_2]
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_time_3]
  const double pos[2] = { 1., 1. }; // we are in 2D space
  double value[1]; // the field is scalar <-> 1 component
  field->getValueOn( pos, 0.5*( time1 + time2 ), value );
  CPPUNIT_ASSERT( fabs( value[0] - 0.5*( 10. + 20. )) < 1e-13 ); 
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_time_3]
}

void CppExample_MEDCouplingFieldDouble_getValueOnMulti()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOnMulti_1]
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh = MEDCouplingCMesh::New();
  mesh->setCoords(coordsArr,coordsArr);
  MCAuto<MEDCouplingFieldDouble> field =
    mesh->fillFromAnalytic( MEDCoupling::ON_CELLS,1,"x+y");
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOnMulti_1]
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOnMulti_2]
  // field values are located at cell barycenters
  MCAuto<DataArrayDouble> bc = mesh->computeCellCenterOfMass();
  MCAuto<DataArrayDouble> valArray =
    field->getValueOnMulti( bc->getConstPointer(), bc->getNumberOfTuples() );
  CPPUNIT_ASSERT( valArray->isEqual( * field->getArray(), 1e-13 ));
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOnMulti_2]
}

void CppExample_MEDCouplingFieldDouble_getValueOn()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_1]
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh = MEDCouplingCMesh::New();
  mesh->setCoords(coordsArr,coordsArr);
  MCAuto<MEDCouplingFieldDouble> field =
    mesh->fillFromAnalytic( MEDCoupling::ON_CELLS,1,"x+y");
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_1]
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_2]
  // field values are located at cell barycenters
  MCAuto<DataArrayDouble> bc = mesh->computeCellCenterOfMass();
  std::vector<double> vals( field->getNumberOfTuples() ); // array to collect values returned by getValueOn()
  double cellBC[2]; // we are in 2D space
  for ( int i = 0; i < bc->getNumberOfTuples(); ++i )
  {
    bc->getTuple( i, cellBC );
    field->getValueOn( cellBC, & vals[i] );
  }
  CPPUNIT_ASSERT( std::equal( vals.begin(), vals.end(), field->getArray()->getConstPointer() ));
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOn_2]
}

void CppExample_MEDCouplingFieldDouble_getValueOnPos()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOnPos_1]
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh = MEDCouplingCMesh::New();
  mesh->setCoords(coordsArr,coordsArr);
  MCAuto<MEDCouplingFieldDouble> field =
    mesh->fillFromAnalytic( MEDCoupling::ON_CELLS,1,"x+y");
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOnPos_1]
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOnPos_2]
  double val11[1]; // 1 == field->getNumberOfComponents()
  field->getValueOnPos( 1,1,-1, val11 );
  // field values are located at cell barycenters
  MCAuto<DataArrayDouble> bc = mesh->computeCellCenterOfMass();
  CPPUNIT_ASSERT( val11[0] == bc->getIJ(3,0) + bc->getIJ(3,1) );
  //! [CppSnippet_MEDCouplingFieldDouble_getValueOnPos_2]
}

void CppExample_MEDCouplingFieldDouble_renumberNodes()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_renumberNodes_1]
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> cmesh = MEDCouplingCMesh::New();
  cmesh->setCoords(coordsArr,coordsArr);
  MCAuto<MEDCouplingUMesh> mesh = cmesh->buildUnstructured();
  //! [CppSnippet_MEDCouplingFieldDouble_renumberNodes_1]
  //! [CppSnippet_MEDCouplingFieldDouble_renumberNodes_2]
  MCAuto<MEDCouplingFieldDouble> field =
    mesh->fillFromAnalytic( MEDCoupling::ON_NODES,2,"IVec*x+JVec*y");
  const DataArrayDouble* values = field->getArray();
  const DataArrayDouble* nodeCoords = mesh->getCoords();
  CPPUNIT_ASSERT( values->isEqualWithoutConsideringStr( *nodeCoords, 1e-13 ));
  //! [CppSnippet_MEDCouplingFieldDouble_renumberNodes_2]
  //! [CppSnippet_MEDCouplingFieldDouble_renumberNodes_3]
  const int renumber[9] = { 8, 7, 6, 5, 4, 3, 2, 1, 0 };
  field->renumberNodes(renumber,false);
  const MEDCouplingMesh* mesh2 = field->getMesh(); // field now refers to another mesh
  values = field->getArray();
  nodeCoords = (static_cast<const MEDCouplingUMesh*>(mesh2))->getCoords();
  CPPUNIT_ASSERT( values->isEqualWithoutConsideringStr( *nodeCoords, 1e-13 ));
  //! [CppSnippet_MEDCouplingFieldDouble_renumberNodes_3]
}

void CppExample_MEDCouplingFieldDouble_renumberCells()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_renumberCells_1]
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> cmesh = MEDCouplingCMesh::New();
  cmesh->setCoords(coordsArr,coordsArr);
  MCAuto<MEDCouplingUMesh> mesh = cmesh->buildUnstructured();
  //! [CppSnippet_MEDCouplingFieldDouble_renumberCells_1]
  //! [CppSnippet_MEDCouplingFieldDouble_renumberCells_2]
  MCAuto<MEDCouplingFieldDouble> field =
    mesh->fillFromAnalytic( MEDCoupling::ON_CELLS,2,"IVec*x+JVec*y");
  const DataArrayDouble* values = field->getArray();
  MCAuto<DataArrayDouble> bc = mesh->computeCellCenterOfMass();
  CPPUNIT_ASSERT( values->isEqualWithoutConsideringStr( *bc, 1e-13 ));
  //! [CppSnippet_MEDCouplingFieldDouble_renumberCells_2]
  //! [CppSnippet_MEDCouplingFieldDouble_renumberCells_3]
  const int renumber[4] = { 3, 2, 1, 0 };
  field->renumberCells(renumber,false);
  const MEDCouplingMesh* mesh2 = field->getMesh(); // field now refers to another mesh
  values = field->getArray();
  bc = mesh2->computeCellCenterOfMass();
  CPPUNIT_ASSERT( values->isEqualWithoutConsideringStr( *bc, 1e-13 ));
  //! [CppSnippet_MEDCouplingFieldDouble_renumberCells_3]
}

void CppExample_MEDCouplingFieldDouble_buildNewTimeReprFromThis()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingFieldDouble_buildNewTimeReprFromThis_1]
  const double coords[4] = {0.,2.,4.};
  MCAuto<DataArrayDouble> coordsArr = DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh = MEDCouplingCMesh::New();
  mesh->setCoords(coordsArr,coordsArr);
  MCAuto<MEDCouplingFieldDouble> field1 =
    mesh->fillFromAnalytic( MEDCoupling::ON_NODES,1,"x+y");
  CPPUNIT_ASSERT( field1->getTimeDiscretization() == MEDCoupling::ONE_TIME );
  //! [CppSnippet_MEDCouplingFieldDouble_buildNewTimeReprFromThis_1]
  //! [CppSnippet_MEDCouplingFieldDouble_buildNewTimeReprFromThis_2]
  MCAuto<MEDCouplingFieldDouble> field2 =
    field1->buildNewTimeReprFromThis( MEDCoupling::NO_TIME, false );
  CPPUNIT_ASSERT( field2->getTimeDiscretization() == MEDCoupling::NO_TIME );
  //! [CppSnippet_MEDCouplingFieldDouble_buildNewTimeReprFromThis_2]
}

void CppExample_MEDCouplingMesh_fillFromAnalytic3()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic3_1]
  const double coords[4] = {0.,2.,4.,6.}; // 6. is not used
  MCAuto<DataArrayDouble> x = DataArrayDouble::New();
  x->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<DataArrayDouble> y = DataArrayDouble::New();
  y->useExternalArrayWithRWAccess( coords, 2, 1 );
  MCAuto<MEDCouplingCMesh> mesh=MEDCouplingCMesh::New();
  mesh->setCoords(x,y);
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic3_1]
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic3_2]
  const char func[] = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10";
  const char* varNames[2] = { "a", "b" }; // names used to refer to X and Y coord components
  std::vector<std::string> varNamesVec( varNames, varNames+2 );
  MCAuto<MEDCouplingFieldDouble> field =
    mesh->fillFromAnalyticNamedCompo( MEDCoupling::ON_CELLS, 3, varNamesVec, func );
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic3_2]
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic3_3]
  double val1[3]; // a value (vector) of the cell #1
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 ); // 3 components in the field
  field->getArray()->getTuple( 1, val1 );
  //
  MCAuto<DataArrayDouble> bc =
    mesh->computeCellCenterOfMass(); // func is applied to barycenters of cells
  double bc1[2]; // coordinates of the second point
  bc->getTuple( 1, bc1 );
  //
  double dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] );  // "sqrt( a*a + b*b )"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[0], 10 + bc1[1], 13 ); // "10 + IVec * b"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[1], 10 + bc1[0], 13 ); // "10 + JVec * a"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[2], 10 + dist  , 13 ); // "10 + KVec * sqrt( a*a + b*b )"
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic3_3]
}

void CppExample_MEDCouplingMesh_fillFromAnalytic2()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic2_1]
  const double coords[4] = {0.,2.,4.,6.}; // 6. is not used
  MCAuto<DataArrayDouble> x = DataArrayDouble::New();
  x->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<DataArrayDouble> y = DataArrayDouble::New();
  y->useExternalArrayWithRWAccess( coords, 2, 1 );
  x->setInfoOnComponent(0,"a"); //  name used to refer to X coordinate within a function
  y->setInfoOnComponent(0,"b"); //  name used to refer to Y coordinate within a function
  MCAuto<MEDCouplingCMesh> mesh=MEDCouplingCMesh::New();
  mesh->setCoords(x,y);
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic2_1]
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic2_2]
  const char func[] = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10";
  MCAuto<MEDCouplingFieldDouble> field =
    mesh->fillFromAnalyticCompo( MEDCoupling::ON_CELLS, 3, func );
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic2_2]
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic2_3]
  double val1[3]; // a value (vector) of the cell #1
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 ); // 3 components in the field
  field->getArray()->getTuple( 1, val1 );
  //
  MCAuto<DataArrayDouble> bc =
    mesh->computeCellCenterOfMass(); // func is applied to barycenters of cells
  double bc1[2]; // coordinates of the second point
  bc->getTuple( 1, bc1 );
  //
  double dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] );  // "sqrt( a*a + b*b )"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[0], 10 + bc1[1], 13 ); // "10 + IVec * b"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[1], 10 + bc1[0], 13 ); // "10 + JVec * a"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[2], 10 + dist  , 13 ); // "10 + KVec * sqrt( a*a + b*b )"
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic2_3]
}

void CppExample_MEDCouplingMesh_fillFromAnalytic()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic_1]
  const double coords[4] = {0.,2.,4.,6.}; // 6. is not used
  MCAuto<DataArrayDouble> x = DataArrayDouble::New();
  x->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<DataArrayDouble> y = DataArrayDouble::New();
  y->useExternalArrayWithRWAccess( coords, 2, 1 );
  MCAuto<MEDCouplingCMesh> mesh=MEDCouplingCMesh::New();
  mesh->setCoords(x,y);
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic_1]
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic_2]
  const char func[] = "IVec * b + JVec * a + KVec * sqrt( a*a + b*b ) + 10";
  MCAuto<MEDCouplingFieldDouble> field =
    mesh->fillFromAnalytic( MEDCoupling::ON_CELLS, 3, func );
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic_2]
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic_3]
  double val1[3]; // a value (vector) of the cell #1
  CPPUNIT_ASSERT( field->getNumberOfComponents() == 3 ); // 3 components in the field
  field->getArray()->getTuple( 1, val1 );
  //
  MCAuto<DataArrayDouble> bc =
    mesh->computeCellCenterOfMass(); // func is applied to barycenters of cells
  double bc1[2]; // coordinates of the second point
  bc->getTuple( 1, bc1 );
  //
  double dist = sqrt( bc1[0]*bc1[0] + bc1[1]*bc1[1] );  // "sqrt( a*a + b*b )"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[0], 10 + bc1[1], 13 ); // "10 + IVec * b"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[1], 10 + bc1[0], 13 ); // "10 + JVec * a"
  CPPUNIT_ASSERT_DOUBLES_EQUAL( val1[2], 10 + dist  , 13 ); // "10 + KVec * sqrt( a*a + b*b )"
  //! [CppSnippet_MEDCouplingMesh_fillFromAnalytic_3]
}

void CppExample_MEDCouplingCMesh_getCoordsAt()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingCMesh_getCoordsAt_1]
  const double coords[3] = {1.,2.,4.};
  MCAuto<DataArrayDouble> x = DataArrayDouble::New();
  x->useExternalArrayWithRWAccess( coords, 3, 1 );
  MCAuto<MEDCouplingCMesh> mesh=MEDCouplingCMesh::New();
  mesh->setCoordsAt(0,x);
  const DataArrayDouble* x2=mesh->getCoordsAt(0);
  CPPUNIT_ASSERT( x2->isEqual( *x, 1e-13 ));
  //! [CppSnippet_MEDCouplingCMesh_getCoordsAt_1]
}

void CppExample_MEDCouplingUMesh_areCellsIncludedIn()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_areCellsIncludedIn_1]
  MCAuto<MEDCouplingUMesh> mesh1=MEDCouplingUMesh::New();
  mesh1->setMeshDimension(2);
  mesh1->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh1->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // #0
  mesh1->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // #1
  mesh1->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // #2
  mesh1->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // #3
  mesh1->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // #4
  mesh1->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh1->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_areCellsIncludedIn_1]
  //! [CppSnippet_MEDCouplingUMesh_areCellsIncludedIn_2]
  const int cells2[3] = { 4,2,0 }; // even cells selected
  MCAuto<MEDCouplingUMesh> mesh2 =
    (MEDCouplingUMesh*) mesh1->buildPartOfMySelf( cells2, cells2+3, true );
  //! [CppSnippet_MEDCouplingUMesh_areCellsIncludedIn_2]
  //! [CppSnippet_MEDCouplingUMesh_areCellsIncludedIn_3]
  int compType = 0; // the strongest policy
  DataArrayInt *corr2to1, *corr1to2;
  // a larger mesh1 includes a smaller mesh2
  CPPUNIT_ASSERT( mesh1->areCellsIncludedIn( mesh2, compType, corr2to1 ));
  CPPUNIT_ASSERT( std::equal( cells2, cells2+3, corr2to1->getConstPointer() ));
  //! [CppSnippet_MEDCouplingUMesh_areCellsIncludedIn_3]
  //! [CppSnippet_MEDCouplingUMesh_areCellsIncludedIn_4]
  // the smaller mesh2 does NOT include the larger mesh1
  CPPUNIT_ASSERT( ! mesh2->areCellsIncludedIn( mesh1, compType, corr1to2 ));
  const int corr1to2Expected[5] = {2, 3, 1, 4, 0};
  CPPUNIT_ASSERT(std::equal( corr1to2Expected, corr1to2Expected+5, corr1to2->getConstPointer() ));
  //! [CppSnippet_MEDCouplingUMesh_areCellsIncludedIn_4]
  corr2to1->decrRef();
  corr1to2->decrRef();
}

void CppExample_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells_1]
  // 2D coordinates of 5 base nodes
  const double coords[5*2]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2 };
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 5, 2 );
  // coordinates of 5 top nodes
  MCAuto<DataArrayDouble> coordsArr2 = coordsArr->deepCopy();
  // 3D coordinates of base + top nodes
  coordsArr  = coordsArr-> changeNbOfComponents( 3, 0 );
  coordsArr2 = coordsArr2->changeNbOfComponents( 3, 1 );
  coordsArr = DataArrayDouble::Aggregate( coordsArr, coordsArr2 );
  // mesh
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  mesh->setMeshDimension(3);
  mesh->allocateCells(2);
  // connectivity of reversed HEXA8 and PENTA6
  const int conn[8+6]={0,1,4,3, 5,6,9,8, 1,2,4, 6,7,9};
  mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8, 8,conn+0);
  mesh->insertNextCell(INTERP_KERNEL::NORM_PENTA6,6,conn+8);
  mesh->finishInsertingCells();
  //! [CppSnippet_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells_1]
  //! [CppSnippet_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells_2]
  MCAuto<DataArrayInt> fixedCells =
    mesh->findAndCorrectBadOriented3DExtrudedCells();
  CPPUNIT_ASSERT( fixedCells->getNumberOfTuples() == 2 ); // 2 cells fixed
  fixedCells = mesh->findAndCorrectBadOriented3DExtrudedCells();
  CPPUNIT_ASSERT( fixedCells->getNumberOfTuples() == 0 ); // no bad cells
  //! [CppSnippet_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells_2]
}

void CppExample_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented_1]
  // 2D coordinates of 5 base nodes
  const double coords[5*2]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2 };
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess( coords, 5, 2 );
  // coordinates of 5 top nodes
  MCAuto<DataArrayDouble> coordsArr2 = coordsArr->deepCopy();
  // 3D coordinates of base + top nodes
  coordsArr  = coordsArr-> changeNbOfComponents( 3, 0 );
  coordsArr2 = coordsArr2->changeNbOfComponents( 3, 1 );
  coordsArr = DataArrayDouble::Aggregate( coordsArr, coordsArr2 );
  // mesh
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  mesh->setMeshDimension(3);
  mesh->allocateCells(2);
  // connectivity of a HEXA8 + a reversed PENTA6
  const int conn[8+6]={0,3,4,1, 5,8,9,6, 1,2,4, 6,7,9};
  mesh->insertNextCell(INTERP_KERNEL::NORM_POLYHED,8,conn); //  "extruded" polyhedron
  mesh->insertNextCell(INTERP_KERNEL::NORM_POLYHED,6,conn+8);
  mesh->finishInsertingCells();
  // fix connectivity of NORM_POLYHED's
  mesh->convertExtrudedPolyhedra();
  //! [CppSnippet_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented_1]
  //! [CppSnippet_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented_2]
  std::vector<int> badCellIds;
  mesh->arePolyhedronsNotCorrectlyOriented( badCellIds );
  CPPUNIT_ASSERT( badCellIds.size() == 1 ); //  one polyhedron is KO
  // fix invalid rolyherdons
  mesh->orientCorrectlyPolyhedrons();
  // re-check orientation
  badCellIds.clear(); // as badCellIds is not cleared by arePolyhedronsNotCorrectlyOriented()
  mesh->arePolyhedronsNotCorrectlyOriented( badCellIds );
  CPPUNIT_ASSERT( badCellIds.size() == 0 ); // connectivity is OK
  //! [CppSnippet_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented_2]
}

void CppExample_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,2,4, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  mesh->changeSpaceDimension(3);
  //! [CppSnippet_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented_1]
  //! [CppSnippet_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented_2]
  const double vec[3] = {0.,0.,-1.};
  std::vector<int> badCellIds;
  mesh->are2DCellsNotCorrectlyOriented( vec, false, badCellIds );
  CPPUNIT_ASSERT( badCellIds.size() == 1 ); //  one cell is reversed
  // fix orientation
  mesh->orientCorrectly2DCells( vec, false );
  // re-check orientation
  badCellIds.clear(); // as badCellIds is not cleared by are2DCellsNotCorrectlyOriented()
  mesh->are2DCellsNotCorrectlyOriented( vec, false, badCellIds );
  CPPUNIT_ASSERT( badCellIds.size() == 0 ); // the orientation is OK
  //! [CppSnippet_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented_2]
}

void CppExample_MEDCouplingUMesh_getCellsContainingPoints()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_getCellsContainingPoints_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);   
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14);
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_getCellsContainingPoints_1]
  //! [CppSnippet_MEDCouplingUMesh_getCellsContainingPoints_2]
  const double pos[3*2] = { 10., 10,               // point out of the mesh
                            0.3, 0.3,              // point located somewhere inside the mesh
                            coords[2], coords[3]}; // point at the node #1
  const double eps = 1e-4; // ball radius
  MCAuto<DataArrayInt> cells, cellsIndex;
  mesh->getCellsContainingPoints( pos, 3, eps, cells, cellsIndex );
  const int cellsExpected[3]={4, 0, 1};
  const int cellsIndexExpected[4]={0, 0, 1, 3};
  CPPUNIT_ASSERT(std::equal( cellsExpected,      cellsExpected+3,      cells->begin()));
  CPPUNIT_ASSERT(std::equal( cellsIndexExpected, cellsIndexExpected+4, cellsIndex->begin()));
  //! [CppSnippet_MEDCouplingUMesh_getCellsContainingPoints_2]
}

void CppExample_MEDCouplingUMesh_getCellsContainingPoint()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_getCellsContainingPoint_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_getCellsContainingPoint_1]
  //! [CppSnippet_MEDCouplingUMesh_getCellsContainingPoint_2]
  const double* coords4  = coords + 4*2; // coordinates of the node #4
  const double eps = 1e-4; // ball radius
  const double pos[2] = { coords4[0] + eps, coords4[1] - eps }; // ball center
  std::vector<int> cellIds;
  mesh->getCellsContainingPoint( pos, eps, cellIds );
  CPPUNIT_ASSERT ( (int)cellIds.size() == mesh->getNumberOfCells() );
  //! [CppSnippet_MEDCouplingUMesh_getCellsContainingPoint_2]
}

void CppExample_MEDCouplingUMesh_buildPartOrthogonalField()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_buildPartOrthogonalField_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_buildPartOrthogonalField_1]
  //! [CppSnippet_MEDCouplingUMesh_buildPartOrthogonalField_2]
  const int part[4] = {1,2,3,4}; // cell #0 is omitted
  MCAuto<MEDCouplingFieldDouble> vecField=
    mesh->buildPartOrthogonalField( part, part+4 );
  CPPUNIT_ASSERT ( vecField->getArray()->getNumberOfTuples() == 4 );
  CPPUNIT_ASSERT ( vecField->getArray()->getNumberOfComponents() == 3 );
  //! [CppSnippet_MEDCouplingUMesh_buildPartOrthogonalField_2]
}

void CppExample_MEDCouplingUMesh_getPartMeasureField()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_getPartMeasureField_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,2,4, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_getPartMeasureField_1]
  //! [CppSnippet_MEDCouplingUMesh_getPartMeasureField_2]
  const bool isAbs = true;
  const int part[4] = {1,2,3,4}; // cell #0 is omitted
  MCAuto<DataArrayDouble> areaArr=
    mesh->getPartMeasureField( isAbs, part, part+4 );
  CPPUNIT_ASSERT( areaArr->getIJ(0,0) > 0 ); // orientation ignored
  areaArr=mesh->getPartMeasureField( !isAbs, part, part+4 );
  CPPUNIT_ASSERT( areaArr->getIJ(0,0) < 0 ); // orientation considered
  CPPUNIT_ASSERT ( areaArr->getNumberOfTuples() == 4 );
  //! [CppSnippet_MEDCouplingUMesh_getPartMeasureField_2]
  //! [CppSnippet_MEDCouplingUMesh_getPartMeasureField_3]
  const int cellIds[4] = {1,2,3,4}; // cell #0 is omitted
  MCAuto<DataArrayDouble> baryCenters=
    mesh->getPartBarycenterAndOwner( cellIds, cellIds+4 );
  CPPUNIT_ASSERT( baryCenters->getNumberOfTuples() == 4 );
  CPPUNIT_ASSERT( baryCenters->getNumberOfComponents() == mesh->getSpaceDimension() );
  //! [CppSnippet_MEDCouplingUMesh_getPartMeasureField_3]
}

void CppExample_MEDCouplingUMesh_getCellsInBoundingBox()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_getCellsInBoundingBox_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(1);
  const double coords[3*2]={0.,0., 0.,1., 1.,1};
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 3,2);
  mesh->setCoords(coordsArr);
  mesh->allocateCells(1);
  const int conn[3]={0,1,2};
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn);
  mesh->finishInsertingCells();
  //! [CppSnippet_MEDCouplingUMesh_getCellsInBoundingBox_1]
  //! [CppSnippet_MEDCouplingUMesh_getCellsInBoundingBox_2]
  const double bbox[] = {1., 1., 1.001,1.001}; // xMin, xMax, yMin, yMax
  MCAuto<DataArrayInt> cellIdsArr =
    mesh->getCellsInBoundingBox( bbox, 0.0 );
  CPPUNIT_ASSERT( cellIdsArr->getNumberOfTuples() == 0 );
  cellIdsArr = mesh->getCellsInBoundingBox( bbox, 0.1 );
  CPPUNIT_ASSERT( cellIdsArr->getNumberOfTuples() == 1 );
  //! [CppSnippet_MEDCouplingUMesh_getCellsInBoundingBox_2]
}

void CppExample_MEDCouplingUMesh_renumberNodesInConn()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_renumberNodesInConn_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(1);
  const int conn[4]={4,3,2,1};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
  mesh->finishInsertingCells();
  //! [CppSnippet_MEDCouplingUMesh_renumberNodesInConn_1]
  //! [CppSnippet_MEDCouplingUMesh_renumberNodesInConn_2]
  const int old2newIds[] = {-1,3,2,1,0};
  mesh->renumberNodesInConn( old2newIds );
  const int nodes0Expected[] = {0,1,2,3};
  std::vector<int> nodes0;
  mesh->getNodeIdsOfCell( 0, nodes0 );
  CPPUNIT_ASSERT(std::equal( nodes0Expected, nodes0Expected+4, &nodes0[0] ));
  //! [CppSnippet_MEDCouplingUMesh_renumberNodesInConn_2]
}

void CppExample_MEDCouplingUMesh_renumberNodes()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_renumberNodes_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  const double coords[4*2]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.3};
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 4,2);
  mesh->setCoords(coordsArr);
  mesh->allocateCells(0);
  mesh->finishInsertingCells();
  //! [CppSnippet_MEDCouplingUMesh_renumberNodes_1]
  //! [CppSnippet_MEDCouplingUMesh_renumberNodes_2]
  const int newIds[] = { 2,1,0,-1 };
  mesh->renumberNodes(newIds, 3);
  coordsArr = mesh->getCoordinatesAndOwner(); // get a shorten array
  const double coordsExpected[3*2]={0.7,-0.3, 0.2,-0.3, -0.3,-0.3};
  MCAuto<DataArrayDouble> coordsExpectedArr=DataArrayDouble::New();
  coordsExpectedArr->useExternalArrayWithRWAccess(coordsExpected, 3,2);
  CPPUNIT_ASSERT( coordsExpectedArr->isEqual( *coordsArr, 1e-13 ));
  //! [CppSnippet_MEDCouplingUMesh_renumberNodes_2]
  //! [CppSnippet_MEDCouplingUMesh_renumberNodes_3]
  coordsArr->useExternalArrayWithRWAccess(coords, 4,2); // restore old nodes
  const int newIds2[] = { 2,1,0,2 };
  mesh->renumberNodesCenter(newIds2, 3);
  coordsArr = mesh->getCoordinatesAndOwner(); // get a shorten array
  const double coordsExpected2[3*2]={0.7,-0.3, 0.2,-0.3, -0.3, 0.0};
  MCAuto<DataArrayDouble> coordsExpectedArr2=DataArrayDouble::New();
  coordsExpectedArr2->useExternalArrayWithRWAccess(coordsExpected2, 3,2);
  CPPUNIT_ASSERT( coordsExpectedArr2->isEqual( *coordsArr, 1e-13 ));
  //! [CppSnippet_MEDCouplingUMesh_renumberNodes_3]
}

void CppExample_MEDCouplingUMesh_findBoundaryNodes()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_findBoundaryNodes_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);   
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14);
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_findBoundaryNodes_1]
  //! [CppSnippet_MEDCouplingUMesh_findBoundaryNodes_2]
  MCAuto<DataArrayInt> nodeIdsArr=mesh->findBoundaryNodes();
  CPPUNIT_ASSERT( nodeIdsArr->getNumberOfTuples() == mesh->getNumberOfNodes() - 1 );
  //! [CppSnippet_MEDCouplingUMesh_findBoundaryNodes_2]
}

void CppExample_MEDCouplingUMesh_buildBoundaryMesh()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_buildBoundaryMesh_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);   
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14);
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_buildBoundaryMesh_1]
  //! [CppSnippet_MEDCouplingUMesh_buildBoundaryMesh_2]
  MCAuto<MEDCouplingPointSet> mesh1=mesh->buildBoundaryMesh(true);
  MCAuto<MEDCouplingPointSet> mesh2=mesh->buildBoundaryMesh(false);
  CPPUNIT_ASSERT(  coordsArr->isEqual( *mesh1->getCoords(), 1e-13 )); // same nodes
  CPPUNIT_ASSERT( !coordsArr->isEqual( *mesh2->getCoords(), 1e-13 )); // different nodes
  //! [CppSnippet_MEDCouplingUMesh_buildBoundaryMesh_2]
}

void CppExample_MEDCouplingUMesh_buildFacePartOfMySelfNode()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_buildFacePartOfMySelfNode_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_buildFacePartOfMySelfNode_1]
  //! [CppSnippet_MEDCouplingUMesh_buildFacePartOfMySelfNode_2]
  std::vector<int> nodes;
  mesh->getNodeIdsOfCell( 0, nodes );
  const bool allNodes = true;
  MCAuto<MEDCouplingUMesh> mesh1 =
    (MEDCouplingUMesh*)mesh->buildFacePartOfMySelfNode( &nodes[0],&nodes[0]+nodes.size(),allNodes);
  CPPUNIT_ASSERT( mesh1->getNumberOfCells() == 4 ); // 4 segments bounding QUAD4 #0 only
  MCAuto<MEDCouplingUMesh> mesh2 =
    (MEDCouplingUMesh*)mesh->buildFacePartOfMySelfNode( &nodes[0],&nodes[0]+nodes.size(),!allNodes);
  CPPUNIT_ASSERT( mesh2->getNumberOfCells() == 9 ); // more segments added
  //! [CppSnippet_MEDCouplingUMesh_buildFacePartOfMySelfNode_2]
}

void CppExample_MEDCouplingUMesh_buildPartOfMySelfNode()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_buildPartOfMySelfNode_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_buildPartOfMySelfNode_1]
  //! [CppSnippet_MEDCouplingUMesh_buildPartOfMySelfNode_2]
  std::vector<int> nodes;
  mesh->getNodeIdsOfCell( 0, nodes );
  const bool allNodes = true;
  MCAuto<MEDCouplingUMesh> mesh1 =
    (MEDCouplingUMesh*)mesh->buildPartOfMySelfNode( &nodes[0], &nodes[0]+nodes.size(), allNodes);
  MCAuto<MEDCouplingUMesh> mesh2 =
    (MEDCouplingUMesh*)mesh->buildPartOfMySelfNode( &nodes[0], &nodes[0]+nodes.size(),!allNodes);
  CPPUNIT_ASSERT_EQUAL( (int)mesh1->getNumberOfCells(), 1 );
  CPPUNIT_ASSERT_EQUAL( mesh2->getNumberOfCells(), mesh->getNumberOfCells() );
  //! [CppSnippet_MEDCouplingUMesh_buildPartOfMySelfNode_2]
}

void CppExample_MEDCouplingUMesh_getCellIdsLyingOnNodes()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_getCellIdsLyingOnNodes_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_getCellIdsLyingOnNodes_1]
  //! [CppSnippet_MEDCouplingUMesh_getCellIdsLyingOnNodes_2]
  std::vector<int> nodes;
  mesh->getNodeIdsOfCell( 0, nodes );
  const bool allNodes = true;
  DataArrayInt* cellIdsArr1 = mesh->getCellIdsLyingOnNodes( &nodes[0], &nodes[0]+nodes.size(), allNodes);
  DataArrayInt* cellIdsArr2 = mesh->getCellIdsLyingOnNodes( &nodes[0], &nodes[0]+nodes.size(),!allNodes);
  CPPUNIT_ASSERT_EQUAL( (int)cellIdsArr1->getNumberOfTuples(), 1 );
  CPPUNIT_ASSERT_EQUAL( (int)cellIdsArr2->getNumberOfTuples(), (int)mesh->getNumberOfCells() );
  //! [CppSnippet_MEDCouplingUMesh_getCellIdsLyingOnNodes_2]
  cellIdsArr1->decrRef();
  cellIdsArr2->decrRef();
}

void CppExample_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);   
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14);
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds_1]
  //! [CppSnippet_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds_2]
  const int cellIds[2]={1,2};
  std::vector<int> nodes;
  mesh->getNodeIdsOfCell( cellIds[0], nodes );
  mesh->getNodeIdsOfCell( cellIds[1], nodes );
  DataArrayInt* cellIdsArr = mesh->getCellIdsFullyIncludedInNodeIds( &nodes[0], &nodes[0]+nodes.size());
  CPPUNIT_ASSERT(std::equal( cellIds, cellIds+2, cellIdsArr->getPointer() ));
  //! [CppSnippet_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds_2]
  cellIdsArr->decrRef();
}

void CppExample_MEDCouplingUMesh_buildPartOfMySelf()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_buildPartOfMySelf_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_buildPartOfMySelf_1]
  //! [CppSnippet_MEDCouplingUMesh_buildPartOfMySelf_2]
  const int cellIds[2]={1,2};
  MEDCouplingUMesh* mesh2=(MEDCouplingUMesh*)mesh->buildPartOfMySelf(cellIds,cellIds+2,true);
  MEDCouplingUMesh* mesh3=(MEDCouplingUMesh*)mesh->buildPartOfMySelf(cellIds,cellIds+2,false);
  CPPUNIT_ASSERT(  coordsArr->isEqual( *mesh2->getCoords(), 1e-13 )); // same nodes
  CPPUNIT_ASSERT( !coordsArr->isEqual( *mesh3->getCoords(), 1e-13 )); // different nodes
  for ( int i = 0; i < 2; ++i )
    {
      std::vector<int> nodes1, nodes2;
      mesh ->getNodeIdsOfCell(cellIds[i], nodes1);
      mesh2->getNodeIdsOfCell(i, nodes2);
      CPPUNIT_ASSERT( nodes1 == nodes2 ); // cell #cellIds[i] was copied
    }
  //! [CppSnippet_MEDCouplingUMesh_buildPartOfMySelf_2]
  mesh2->decrRef();
  mesh3->decrRef();
}

void CppExample_MEDCouplingUMesh_mergeNodes()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_mergeNodes_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);
  mesh->finishInsertingCells();
  const double coords[6*2]={0.3,-0.301,  // #0
                            0.2,-0.3,    // #1
                            0.3,-0.302,  // #2 ~~ #0
                            1.1,0.0,     // #3
                            1.1,0.0,     // #4 == #3
                            0.3,-0.303}; // #5 ~~ #0
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(6,2);
  std::copy(coords,coords+6*2,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_mergeNodes_1]
  //! [CppSnippet_MEDCouplingUMesh_mergeNodes_2]
  bool areNodesMerged; int newNbOfNodes;
  MCAuto<DataArrayInt> arr=
    mesh->mergeNodes(0.004,areNodesMerged,newNbOfNodes);
  const int idsExpected[6] = {0, 1, 0, 2, 2, 0};
  CPPUNIT_ASSERT(std::equal(idsExpected,idsExpected+6,arr->getPointer()));
  CPPUNIT_ASSERT( areNodesMerged );
  CPPUNIT_ASSERT_EQUAL( 3, newNbOfNodes );
  //! [CppSnippet_MEDCouplingUMesh_mergeNodes_2]
  //! [CppSnippet_MEDCouplingUMesh_mergeNodes_3]
  const double* baryCoords2 = coords + 2*2; // initial coordinates of node #2
  coordsArr=mesh->getCoordinatesAndOwner(); // retrieve a new shorten coord array
  CPPUNIT_ASSERT( fabs( baryCoords2[1] - coordsArr->getIJ(0,1)) > 1e-4 ); // Y of node #0 differs from that of baryCoords2
  // restore coordinates
  coordsArr->alloc(6,2);
  std::copy(coords,coords+6*2,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  // call mergeNodesCenter()
  arr = mesh->mergeNodesCenter(0.004,areNodesMerged,newNbOfNodes);
  coordsArr=mesh->getCoordinatesAndOwner(); // retrieve a new shorten coord array
  CPPUNIT_ASSERT_DOUBLES_EQUAL( baryCoords2[1], coordsArr->getIJ(0,1), 13 ); // Y of node #0 equals to that of baryCoords2
  //! [CppSnippet_MEDCouplingUMesh_mergeNodes_3]
}

void CppExample_MEDCouplingUMesh_zipConnectivityTraducer()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_zipConnectivityTraducer_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[11]={0,3,4,1, 1,4,2, 4,1,0,3};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+0); // 0     
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4); // 1     
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4); // 2 == 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+0); // 3 == 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+7); // 4 ~~ 0
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_zipConnectivityTraducer_1]
  //! [CppSnippet_MEDCouplingUMesh_zipConnectivityTraducer_2]
  const int oldNbCells = mesh->getNumberOfCells();
  DataArrayInt *arr = mesh->zipConnectivityTraducer(0);
  CPPUNIT_ASSERT_EQUAL( oldNbCells-2, (int)mesh->getNumberOfCells() );
  const int idsExpected[5] = {0, 1, 1, 0, 2};
  CPPUNIT_ASSERT(std::equal(idsExpected,idsExpected+5,arr->getPointer()));
  //! [CppSnippet_MEDCouplingUMesh_zipConnectivityTraducer_2]
  arr->decrRef();
}

void CppExample_MEDCouplingUMesh_zipCoordsTraducer()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_zipCoordsTraducer_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);   
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14);
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_zipCoordsTraducer_1]
  //! [CppSnippet_MEDCouplingUMesh_zipCoordsTraducer_2]
  const int cellIds[2]={1,2};
  MEDCouplingUMesh* mesh2=(MEDCouplingUMesh*)mesh->buildPartOfMySelf(cellIds,cellIds+2,true);
  DataArrayInt *arr=mesh2->zipCoordsTraducer();
  CPPUNIT_ASSERT_EQUAL( 4, mesh2->getNumberOfNodes() ); // nb of nodes decreased
  CPPUNIT_ASSERT_EQUAL( (int)mesh->getNumberOfNodes(), (int)arr->getNumberOfTuples() );
  const int idsExpected[9] = {-1,0,1,-1,2,3,-1,-1,-1}; // -1 for unused nodes
  CPPUNIT_ASSERT(std::equal(idsExpected,idsExpected+9,arr->getPointer()));
  //! [CppSnippet_MEDCouplingUMesh_zipCoordsTraducer_2]
  mesh2->decrRef();
  arr->decrRef();
}

void CppExample_MEDCouplingUMesh_getNodeIdsInUse()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_getNodeIdsInUse_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);   
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7); 
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14);
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_getNodeIdsInUse_1]
  //! [CppSnippet_MEDCouplingUMesh_getNodeIdsInUse_2]
  const int cellIds[2]={1,2};
  MEDCouplingUMesh* mesh2=(MEDCouplingUMesh*)mesh->buildPartOfMySelf(cellIds,cellIds+2,true);
  int newNbOfNodes = 0;
  DataArrayInt *arr=mesh2->getNodeIdsInUse( newNbOfNodes );
  const int idsExpected[9] = {-1,0,1,-1,2,3,-1,-1,-1};
  CPPUNIT_ASSERT(std::equal(idsExpected,idsExpected+9,arr->getPointer()));
  //! [CppSnippet_MEDCouplingUMesh_getNodeIdsInUse_2]
  //! [CppSnippet_MEDCouplingUMesh_getNodeIdsInUse_3]
  DataArrayInt *arr2=arr->invertArrayO2N2N2O(newNbOfNodes);
  const int idsExpected2[4] = {1,2,4,5};
  CPPUNIT_ASSERT(std::equal(idsExpected2,idsExpected2+4,arr2->getPointer()));
  //! [CppSnippet_MEDCouplingUMesh_getNodeIdsInUse_3]
  mesh2->decrRef();
  arr->decrRef();
  arr2->decrRef();
}

void CppExample_MEDCouplingUMesh_convertToPolyTypes()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_convertToPolyTypes_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_convertToPolyTypes_1]
  //! [CppSnippet_MEDCouplingUMesh_convertToPolyTypes_2]
  const int cells[2]={1,3};
  mesh->convertToPolyTypes(cells, cells+2);
  CPPUNIT_ASSERT( mesh->getTypeOfCell(0) == INTERP_KERNEL::NORM_QUAD4 );
  CPPUNIT_ASSERT( mesh->getTypeOfCell(1) == INTERP_KERNEL::NORM_POLYGON );
  CPPUNIT_ASSERT( mesh->getTypeOfCell(2) == INTERP_KERNEL::NORM_TRI3 );
  CPPUNIT_ASSERT( mesh->getTypeOfCell(3) == INTERP_KERNEL::NORM_POLYGON );
  //! [CppSnippet_MEDCouplingUMesh_convertToPolyTypes_2]
}

void CppExample_MEDCouplingUMesh_buildDescendingConnectivity2()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity2_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity2_1]
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity2_2]
  DataArrayInt *desc       =DataArrayInt::New();
  DataArrayInt *descIndx   =DataArrayInt::New();
  DataArrayInt *revDesc    =DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  MEDCouplingUMesh * mesh2 = mesh->buildDescendingConnectivity2(desc,descIndx,revDesc,revDescIndx);
  const int descExpected[]        = {1,2,3,4,-3,5,6,7,8,-5,9,10,-2,11,12,13,-7,-10};
  const int descIndxExpected[]    = {0,4,7,10,14,18};
  const int revDescExpected[]     = {0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4};
  const int revDescIndxExpected[] = {0,1,3,5,6,8,9,11,12,13,15,16,17,18};
  CPPUNIT_ASSERT(std::equal(descExpected,descExpected+18,desc->getPointer()));
  CPPUNIT_ASSERT(std::equal(descIndxExpected,descIndxExpected+6,descIndx->getPointer()));
  CPPUNIT_ASSERT(std::equal(revDescExpected,revDescExpected+18,revDesc->getPointer()));
  CPPUNIT_ASSERT(std::equal(revDescIndxExpected,revDescIndxExpected+14,revDescIndx->getPointer()));
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity2_2]
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity2_3]
  const int cell2ConnExpect[] = {4,1};
  std::vector<int> cell2Conn;
  mesh2->getNodeIdsOfCell( 3-1, cell2Conn ); // cell #3 in FORTRAN mode
  CPPUNIT_ASSERT(std::equal(cell2ConnExpect,cell2ConnExpect+2,&cell2Conn[0]));
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity2_3]
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
  mesh2->decrRef();
}

void CppExample_MEDCouplingUMesh_buildDescendingConnectivity()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity_1]
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity_2]
  DataArrayInt *desc       =DataArrayInt::New();
  DataArrayInt *descIndx   =DataArrayInt::New();
  DataArrayInt *revDesc    =DataArrayInt::New();
  DataArrayInt *revDescIndx=DataArrayInt::New();
  MEDCouplingUMesh * mesh2 = mesh->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  const int descExpected[]        = {0,1,2,3, 2,4,5, 6,7,4, 8,9,1,10, 11,12,6,9};
  const int descIndxExpected[]    = {0,4,7,10,14,18};
  const int revDescExpected[]     = {0, 0,3, 0,1, 0, 1,2, 1, 2,4, 2, 3, 3,4, 3, 4, 4};
  const int revDescIndxExpected[] = {0,1,3,5,6,8,9,11,12,13,15,16,17,18};
  CPPUNIT_ASSERT(std::equal(descExpected,descExpected+18,desc->getPointer()));
  CPPUNIT_ASSERT(std::equal(descIndxExpected,descIndxExpected+6,descIndx->getPointer()));
  CPPUNIT_ASSERT(std::equal(revDescExpected,revDescExpected+18,revDesc->getPointer()));
  CPPUNIT_ASSERT(std::equal(revDescIndxExpected,revDescIndxExpected+14,revDescIndx->getPointer()));
  //! [CppSnippet_MEDCouplingUMesh_buildDescendingConnectivity_2]
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
  mesh2->decrRef();
}

void CppExample_MEDCouplingUMesh_getReverseNodalConnectivity()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_getReverseNodalConnectivity_1]
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  mesh->allocateCells(5);
  const int conn[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);    // 0
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+4);  // 1
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+7);  // 2
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+10); // 3
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+14); // 4
  mesh->finishInsertingCells();
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->alloc(9,2);
  const double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  std::copy(coords,coords+18,coordsArr->getPointer());
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingUMesh_getReverseNodalConnectivity_1]
  //! [CppSnippet_MEDCouplingUMesh_getReverseNodalConnectivity_2]
  DataArrayInt *revNodal=DataArrayInt::New();
  DataArrayInt *revNodalIndx=DataArrayInt::New();
  mesh->getReverseNodalConnectivity(revNodal,revNodalIndx);
  const int revNodalExpected[18]={0,0,1,1,2,0,3,0,1,2,3,4,2,4,3,3,4,4};
  const int revNodalIndexExpected[10]={0,1,3,5,7,12,14,15,17,18};
  CPPUNIT_ASSERT(std::equal(revNodalExpected,revNodalExpected+18,revNodal->getPointer()));
  CPPUNIT_ASSERT(std::equal(revNodalIndexExpected,revNodalIndexExpected+10,revNodalIndx->getPointer()));
  //! [CppSnippet_MEDCouplingUMesh_getReverseNodalConnectivity_2]
  revNodal->decrRef();
  revNodalIndx->decrRef();
}

void CppExample_MEDCouplingUMesh_checkDeepEquivalWith()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingUMesh_checkDeepEquivalWith_1]
  // mesh 1
  MEDCouplingUMesh *mesh1=MEDCouplingUMesh::New();
  const double coords[4*2]={0.0,0.0,  // #0
                            1.0,0.0,  // #1
                            1.0,1.0,  // #2
                            0.0,1.0}; // #3
  {
    mesh1->setMeshDimension(2);
    MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
    coordsArr->useExternalArrayWithRWAccess( coords, 4, 2 );
    mesh1->setCoords(coordsArr);
    mesh1->allocateCells(2);
    const int conn[6]={0,1,2, 1,2,3};
    mesh1->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+0);  // #0
    mesh1->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+3);  // #1
    mesh1->finishInsertingCells();
  }
  // mesh 2
  MEDCouplingUMesh *mesh2=MEDCouplingUMesh::New();
  const double coords2[4*2]={0.0,1.0,    // #0 = #3
                             0.0,0.0,    // #1 = #0
                             1.0,0.0,    // #2 = #1
                             1.0,1.001}; // #3 ~ #2
  {
    mesh2->setMeshDimension(2);
    MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
    coordsArr->useExternalArrayWithRWAccess( coords2, 4, 2 );
    mesh2->setCoords(coordsArr);
    mesh2->allocateCells(2);
    const int conn[6]={2,3,0, 3,1,2};
    mesh2->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+0);  // #0 = #1
    mesh2->insertNextCell(INTERP_KERNEL::NORM_TRI3,3, conn+3);  // #1 ~ #0
    mesh2->finishInsertingCells();
  }
  //! [CppSnippet_MEDCouplingUMesh_checkDeepEquivalWith_1]
  //! [CppSnippet_MEDCouplingUMesh_checkDeepEquivalWith_2]
  int cellCompPol = 1; // "permuted same orientation" - policy of medium severity
  DataArrayInt *nOld2New, *cOld2New;
  mesh1->checkDeepEquivalWith( mesh2, cellCompPol, 0.002, cOld2New, nOld2New );
  const int nOld2NewExpected[4] = { 3, 0, 1, 2 };
  const int cOld2NewExpected[2] = { 1, 0 };
  CPPUNIT_ASSERT(std::equal(nOld2NewExpected,nOld2NewExpected+4,nOld2New->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(cOld2NewExpected,cOld2NewExpected+2,cOld2New->getConstPointer()));
  //! [CppSnippet_MEDCouplingUMesh_checkDeepEquivalWith_2]
  //! [CppSnippet_MEDCouplingUMesh_checkDeepEquivalWith_3]
  cOld2New->decrRef(); // else memory leaks
  CPPUNIT_ASSERT_THROW ( mesh1->checkDeepEquivalOnSameNodesWith( mesh2, cellCompPol, 0.002, cOld2New ), INTERP_KERNEL::Exception );
  mesh2->setCoords( mesh1->getCoords() ); // make meshes share the same coordinates array
  mesh2->allocateCells(2);
  const int conn[6]={1,2,3, 1,0,2};
  mesh2->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+0); // #0 = #1
  mesh2->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+3); // #1 ~ #0
  mesh2->finishInsertingCells();
  cellCompPol = 2; // the weakest policy
  mesh1->checkDeepEquivalOnSameNodesWith( mesh2, cellCompPol, 0, cOld2New );
  //! [CppSnippet_MEDCouplingUMesh_checkDeepEquivalWith_3]
  nOld2New->decrRef();
  cOld2New->decrRef();
  mesh1->decrRef();
  mesh2->decrRef();
}

void CppExample_MEDCouplingPointSet_scale()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingPointSet_scale_1]
  double coords[4*2]={0.0,0.0, 1.0,0.0, 1.0,1.0, 0.0,1.0}; // 2D coordinates of 4 nodes
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 4,2);
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  DataArrayDouble *initCoords = coordsArr->deepCopy();
  //! [CppSnippet_MEDCouplingPointSet_scale_1]
  //! [CppSnippet_MEDCouplingPointSet_scale_2]
  const double center[2] = {0.,0.};
  const double factor = 2.;
  mesh->scale( center, factor );
  //! [CppSnippet_MEDCouplingPointSet_scale_2]
  //! [CppSnippet_MEDCouplingPointSet_scale_3]
  const DataArrayDouble * coordsArr2 = mesh->getCoords();
  CPPUNIT_ASSERT( coordsArr2->isEqualWithoutConsideringStr( *initCoords, 1.0 ));
  CPPUNIT_ASSERT( !coordsArr2->isEqualWithoutConsideringStr( *initCoords, 0.9 ));
  // release data
  initCoords->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_scale_3]
}

void CppExample_MEDCouplingPointSet_translate()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingPointSet_translate_1]
  double coords[4*2]={0.0,0.0, 1.0,0.0, 1.0,1.0, 0.0,1.0}; // 2D coordinates of 4 nodes
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 4,2);
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  DataArrayDouble *initCoords = coordsArr->deepCopy();
  //! [CppSnippet_MEDCouplingPointSet_translate_1]
  //! [CppSnippet_MEDCouplingPointSet_translate_2]
  double vector[2] = {1.,1.};
  mesh->translate( vector );
  //! [CppSnippet_MEDCouplingPointSet_translate_2]
  //! [CppSnippet_MEDCouplingPointSet_translate_3]
  const DataArrayDouble * coordsArr2 = mesh->getCoords();
  CPPUNIT_ASSERT( coordsArr2->isEqualWithoutConsideringStr( *initCoords, 1.0 ));
  CPPUNIT_ASSERT( !coordsArr2->isEqualWithoutConsideringStr( *initCoords, 0.9 ));
  // release data
  initCoords->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_translate_3]
}

void CppExample_MEDCouplingPointSet_rotate()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingPointSet_rotate_1]
  double coords[4*2]={0.0,0.0, 0.1,0.0, 0.1,0.1, 0.0,0.1}; // 2D coordinates of 4 nodes
  double coordsOrig[4*2];
  std::copy(coords,coords+sizeof(coords)/sizeof(double),coordsOrig);//keep tracks of initial values
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 4,2);
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingPointSet_rotate_1]
  //! [CppSnippet_MEDCouplingPointSet_rotate_2]
  double center[3] = {0.,0.,0.}; // it suits for 2D as well
  double vector[3] = {0.,0.,1.}; // it is not used in 2D
  mesh->rotate( center, vector, -M_PI/2); // warning here C++ 'coords' array (defined above) has been modified !
  //! [CppSnippet_MEDCouplingPointSet_rotate_2]
  //! [CppSnippet_MEDCouplingPointSet_rotate_3]
  mesh->changeSpaceDimension(3);
  mesh->rotate( center, vector, +M_PI/2);
  //! [CppSnippet_MEDCouplingPointSet_rotate_3]
  //! [CppSnippet_MEDCouplingPointSet_rotate_4]
  mesh->changeSpaceDimension(2);
  const DataArrayDouble * coordsArr2 = mesh->getCoords();
  coordsArr->useExternalArrayWithRWAccess(coordsOrig, 4,2);
  CPPUNIT_ASSERT( coordsArr2->isEqualWithoutConsideringStr( *coordsArr, 1e-13 ));
  //! [CppSnippet_MEDCouplingPointSet_rotate_4]
}

void CppExample_MEDCouplingPointSet_getBoundingBox()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingPointSet_getBoundingBox_1]
  double cc[2*3]={0.0, 0.1, 0.2, // 3D coordinates of 2 nodes
                  2.0, 2.1, 2.2};
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(cc, 2,3);
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingPointSet_getBoundingBox_1]
  //! [CppSnippet_MEDCouplingPointSet_getBoundingBox_2]
  double bbox[3][2];
  mesh->getBoundingBox( (double*) bbox );

  // check the returned coordinates of extremum points of the bounding box
  for ( int i = 0; i < 2; ++i )   // point id
    for ( int j = 0; j < 3; ++j ) // component
      CPPUNIT_ASSERT_DOUBLES_EQUAL( cc[ i*3 + j ], bbox[j][i], 1e-13);
  //! [CppSnippet_MEDCouplingPointSet_getBoundingBox_2]
}

void CppExample_MEDCouplingPointSet_getNodeIdsNearPoint()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoint_1]
  // 2D coordinates of 5 nodes
  double coords[5*2]={0.3,-0.30001, // #0
                      0.2,-0.3,   // #1
                      0.3,-0.30002, // #2
                      1.1,0.0,    // #3
                      0.3,-0.30003};// #4
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 5,2);
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoint_1]
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoint_2]
  double point [2]={0.3, -0.3}; // point close to nodes #0, #2 and #4
  DataArrayInt *ids = mesh->getNodeIdsNearPoint(point, 1e-2);

  // check found ids
  const int expectedIDs[3] = {0,2,4};
  DataArrayInt * okIDs = ids->findIdsEqualList ( expectedIDs, expectedIDs+3 );
  CPPUNIT_ASSERT_EQUAL(3, (int)okIDs->getNumberOfTuples());

  // release data
  ids->decrRef();
  okIDs->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoint_2]
}
void CppExample_MEDCouplingPointSet_getNodeIdsNearPoints()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoints_1]
  // 2D coordinates of 7 nodes
  double coords[7*2]={0.3,-0.301, // #0
                      0.2,-0.3,   // #1
                      0.3,-0.302, // #2
                      1.1,0.0,    // #3
                      1.1,0.0,    // #4
                      1.1,0.002,  // #5
                      0.3,-0.303};// #6
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 7,2);
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoints_1]
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoints_2]
  const int nbOfPoints = 3;
  double points [nbOfPoints*2]={0.2,-0.30001,  // ~ node #1
                                0.0, 0.0,
                                1.1, 0.002}; // ~ nodes #3, #4 and #5
  DataArrayInt *ids, *idsIndex;
  mesh->getNodeIdsNearPoints(points, nbOfPoints, 1e-1,ids,idsIndex);

  // check found ids (i.e. contents of 'ids' array)
  const int expectedIDs[4] = {1, 3, 4, 5};
  DataArrayInt * okIDs = ids->findIdsEqualList ( expectedIDs, expectedIDs+4 );
  CPPUNIT_ASSERT_EQUAL(4, (int)okIDs->getNumberOfTuples());

  // release data
  ids->decrRef();
  idsIndex->decrRef();
  okIDs->decrRef();
  //! [CppSnippet_MEDCouplingPointSet_getNodeIdsNearPoints_2]
}

void CppExample_MEDCouplingPointSet_findCommonNodes()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingPointSet_findCommonNodes_1]
  double coords[6*2]={0.3,-0.301, // 0
                      0.2,-0.3,   // 1
                      0.3,-0.302, // 2
                      1.1,0.0,    // 3
                      1.1,0.0,    // 4
                      0.3,-0.303};// 5
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 6,2);
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingPointSet_findCommonNodes_1]
  //! [CppSnippet_MEDCouplingPointSet_findCommonNodes_2]
  DataArrayInt *com, *comI;
  mesh->findCommonNodes(1e-13,-1,com,comI);
  CPPUNIT_ASSERT_EQUAL(2, (int)com->getNumberOfTuples());
  com->decrRef(); comI->decrRef();
  mesh->findCommonNodes(0.004,-1,com,comI);
  CPPUNIT_ASSERT_EQUAL(5, (int)com->getNumberOfTuples());
  //! [CppSnippet_MEDCouplingPointSet_findCommonNodes_2]
  com->decrRef(); comI->decrRef();
}

void CppExample_MEDCouplingPointSet_getCoordinatesOfNode()
{
  using namespace MEDCoupling;
  //! [CppSnippet_MEDCouplingPointSet_getCoordinatesOfNode_1]
  double coords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3};
  MCAuto<DataArrayDouble> coordsArr=DataArrayDouble::New();
  coordsArr->useExternalArrayWithRWAccess(coords, 3,2);
  MCAuto<MEDCouplingUMesh> mesh=MEDCouplingUMesh::New();
  mesh->setCoords(coordsArr);
  //! [CppSnippet_MEDCouplingPointSet_getCoordinatesOfNode_1]
  //! [CppSnippet_MEDCouplingPointSet_getCoordinatesOfNode_2]
  std::vector<double> coords2;
  mesh->getCoordinatesOfNode(1,coords2);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(coords[2],coords2[0],1e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(coords[3],coords2[1],1e-13);
  //! [CppSnippet_MEDCouplingPointSet_getCoordinatesOfNode_2]
}

void CppExample_DataArrayInt_buildPermutationArr()
{
  using namespace MEDCoupling;
  //! [CppSnippet_DataArrayInt_buildPermutationArr_1]
  DataArrayInt *a=DataArrayInt::New();
  const int vala[5]={4,5,6,7,8};
  a->alloc(5,1);
  std::copy(vala,vala+5,a->getPointer());
  DataArrayInt *b=DataArrayInt::New();
  const int valb[5]={5,4,8,6,7};
  b->alloc(5,1);
  std::copy(valb,valb+5,b->getPointer());
  DataArrayInt *c=a->buildPermutationArr(*b);
  //! [CppSnippet_DataArrayInt_buildPermutationArr_1]
  const int expect1[5]={1,0,4,2,3};
  CPPUNIT_ASSERT_EQUAL(5,(int)c->getNumberOfTuples());
  CPPUNIT_ASSERT_EQUAL(1,(int)c->getNumberOfComponents());
  CPPUNIT_ASSERT(std::equal(expect1,expect1+5,c->getConstPointer()));
  CPPUNIT_ASSERT(a->isEqualWithoutConsideringStrAndOrder(*b));
  a->decrRef();
  b->decrRef();
  c->decrRef();
}

void CppExample_DataArrayInt_invertArrayO2N2N2O()
{
  using namespace MEDCoupling;
  //! [CppSnippet_DataArrayInt_invertArrayO2N2N2O_1]
  const int arr1[6]={2,0,4,1,5,3};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(6,1);
  std::copy(arr1,arr1+6,da->getPointer());
  DataArrayInt *da2=da->invertArrayO2N2N2O(6);
  const int expected1[6]={1,3,0,5,2,4};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(i,0));
  //! [CppSnippet_DataArrayInt_invertArrayO2N2N2O_1]
  da->decrRef();
  da2->decrRef();
}

void CppExample_DataArrayInt_invertArrayN2O2O2N()
{
  using namespace MEDCoupling;
  //! [CppSnippet_DataArrayInt_invertArrayN2O2O2N_1]
  const int arr1[6]={2,0,4,1,5,3};
  DataArrayInt *da=DataArrayInt::New();
  da->alloc(6,1);
  std::copy(arr1,arr1+6,da->getPointer());
  DataArrayInt *da2=da->invertArrayN2O2O2N(6);
  const int expected1[6]={1,3,0,5,2,4};
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(expected1[i],da2->getIJ(i,0));
  //! [CppSnippet_DataArrayInt_invertArrayN2O2O2N_1]
  da->decrRef();
  da2->decrRef();
}

void CppExample_DataArrayDouble_getIdsInRange()
{
  using namespace MEDCoupling;
  //! [CppSnippet_DataArrayDouble_getIdsInRange_1]
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(10,1);
  da->iota();

  DataArrayInt* da2 = da->findIdsInRange( 2.5, 6 );
  //! [CppSnippet_DataArrayDouble_getIdsInRange_1]
  da->decrRef();
  da2->decrRef();
}

void CppExample_DataArrayDouble_findCommonTuples()
{
  using namespace MEDCoupling;
  //! [CppSnippet_DataArrayDouble_findCommonTuples1]
  DataArrayDouble *da=DataArrayDouble::New();
  da->alloc(6,2);
  const double array2[12]={2.3,2.3, // 0
                           1.2,1.2, // 1
                           1.3,1.3, // 2
                           2.3,2.3, // 3
                           2.301,   // 4
                           2.301,   // 5
                           0.8,0.8};// 6
  std::copy(array2,array2+12,da->getPointer());
  //! [CppSnippet_DataArrayDouble_findCommonTuples1]
  //! [CppSnippet_DataArrayDouble_findCommonTuples2]
  DataArrayInt *c=0,*cI=0;
  da->findCommonTuples(1.01e-1,-1,c,cI);

  const int expected3[5]={0,3,4,1,2};
  const int expected4[3]={0,3,5};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+5,c->getConstPointer()));
  CPPUNIT_ASSERT(std::equal(expected4,expected4+3,cI->getConstPointer()));
  c->decrRef();
  cI->decrRef();
  da->decrRef();
  //! [CppSnippet_DataArrayDouble_findCommonTuples2]
}

void CppExample_DataArrayDouble_Meld1()
{
  using namespace MEDCoupling;
  //! [CppSnippet_DataArrayDouble_Meld1_1]
  const int sameNbTuples = 7;

  DataArrayDouble *da1=DataArrayDouble::New();
  da1->alloc(sameNbTuples,2);
  da1->fillWithValue(7.);
  da1->setInfoOnComponent(0,"c0da1");
  da1->setInfoOnComponent(1,"c1da1");

  DataArrayDouble *da2=DataArrayDouble::New();
  da2->alloc(sameNbTuples,1);
  da2->iota(0.);
  da2->setInfoOnComponent(0,"c0da2");

  da1->meldWith(da2);
  //! [CppSnippet_DataArrayDouble_Meld1_1]
  //! [CppSnippet_DataArrayDouble_Meld1_2]
  da1->decrRef();
  da2->decrRef();
  //! [CppSnippet_DataArrayDouble_Meld1_2]
}

void CppExample_DataArrayInt_Meld1()
{
  using namespace MEDCoupling;
  //! [CppSnippet_DataArrayInt_Meld1_1]
  const int sameNbTuples = 7;

  DataArrayInt *da1=DataArrayInt::New();
  da1->alloc(sameNbTuples,2);
  da1->fillWithValue(7);
  da1->setInfoOnComponent(0,"c0da1");
  da1->setInfoOnComponent(1,"c1da1");

  DataArrayInt *da2=DataArrayInt::New();
  da2->alloc(sameNbTuples,1);
  da2->iota(0);
  da2->setInfoOnComponent(0,"c0da2");

  da1->meldWith(da2);
  //! [CppSnippet_DataArrayInt_Meld1_1]
  //! [CppSnippet_DataArrayInt_Meld1_2]
  da1->decrRef();
  da2->decrRef();
  //! [CppSnippet_DataArrayInt_Meld1_2]
}

void CppExampleFieldDoubleBuildSubPart1()
{
  //! [CppSnippetFieldDoubleBuildSubPart1_1]
  MEDCoupling::MEDCouplingUMesh *mesh1=MEDCoupling::MEDCouplingBasicsTest::build2DTargetMesh_1();
  MEDCoupling::MEDCouplingFieldDouble *f1=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::ONE_TIME);
  f1->setTime(2.3,5,6);
  f1->setMesh(mesh1);
  MEDCoupling::DataArrayDouble *array=MEDCoupling::DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfCells(),2);
  const double arr1[10]={3.,103.,4.,104.,5.,105.,6.,106.,7.,107.};
  std::copy(arr1,arr1+10,array->getPointer());
  f1->setArray(array);
  array->decrRef();
  //! [CppSnippetFieldDoubleBuildSubPart1_1]
  //! [CppSnippetFieldDoubleBuildSubPart1_2]
  const int part1[3]={2,1,4};
  MEDCoupling::MEDCouplingFieldDouble *f2=f1->buildSubPart(part1,part1+3);
  //! [CppSnippetFieldDoubleBuildSubPart1_2]
  f2->zipCoords();
  CPPUNIT_ASSERT_EQUAL(3,(int)f2->getMesh()->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(6,f2->getMesh()->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,f2->getMesh()->getMeshDimension());
  MEDCoupling::MEDCouplingUMesh *m2C=dynamic_cast<MEDCoupling::MEDCouplingUMesh *>(const_cast<MEDCoupling::MEDCouplingMesh *>(f2->getMesh()));
  CPPUNIT_ASSERT_EQUAL(13,m2C->getNodalConnectivityArrayLen());
  const double expected2[12]={0.2, -0.3, 0.7, -0.3, 0.2, 0.2, 0.7, 0.2, 0.2, 0.7, 0.7, 0.7};
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],m2C->getCoords()->getIJ(0,i),1.e-12);
  const double expected3[13]={3,2,3,1,3,0,2,1,4,4,5,3,2};
  CPPUNIT_ASSERT(std::equal(expected3,expected3+13,m2C->getNodalConnectivity()->getConstPointer()));
  const double expected4[4]={0,4,8,13};
  CPPUNIT_ASSERT(std::equal(expected4,expected4+4,m2C->getNodalConnectivityIndex()->getConstPointer()));
  f2->decrRef();
  f1->decrRef();
  //! [CppSnippetFieldDoubleBuildSubPart1_3]
  f1=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES,MEDCoupling::ONE_TIME);
  f1->setTime(2.3,5,6);
  f1->setMesh(mesh1);
  array=MEDCoupling::DataArrayDouble::New();
  array->alloc(mesh1->getNumberOfNodes(),2);
  const double arr2[18]={3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.};
  std::copy(arr2,arr2+18,array->getPointer());  
  f1->setArray(array);
  array->decrRef();
  //! [CppSnippetFieldDoubleBuildSubPart1_3]
  //! [CppSnippetFieldDoubleBuildSubPart1_4]
  const int part2[2]={1,2};
  f2=f1->buildSubPart(part2,part2+2);
  //! [CppSnippetFieldDoubleBuildSubPart1_4]
  f2->decrRef();
  //idem previous because nodes of cell#4 are not fully present in part3 
  const int part3[2]={1,2};
  MEDCoupling::DataArrayInt *arrr=MEDCoupling::DataArrayInt::New();
  arrr->alloc(2,1);
  std::copy(part3,part3+2,arrr->getPointer());
  f2=f1->buildSubPart(arrr);
  arrr->decrRef();
  f2->decrRef();
  //
  const int part4[3]={1,2,4};
  f2=f1->buildSubPart(part4,part4+3);
  f2->decrRef();
  //
  f1->decrRef();
  mesh1->decrRef();
  return;
}

void CppSnippetUMeshStdBuild1()
{
  //! [CppSnippetUMeshStdBuild1_1]
  double coords[27]={-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                     0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. };
  int nodalConnPerCell[18]={0,3,4,1, 1,4,2, 4,5,2, 6,7,4,3, 7,8,5,4};
  //! [CppSnippetUMeshStdBuild1_1]
  //! [CppSnippetUMeshStdBuild1_2]
  MEDCoupling::MEDCouplingUMesh *mesh=MEDCoupling::MEDCouplingUMesh::New("My2DMesh",2);
  //! [CppSnippetUMeshStdBuild1_2]
  //! [CppSnippetUMeshStdBuild1_3]
  mesh->allocateCells(5);//You can put more than 5 if you want but not less.
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,nodalConnPerCell);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,nodalConnPerCell+4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,nodalConnPerCell+7);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,nodalConnPerCell+10);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,nodalConnPerCell+14);
  mesh->finishInsertingCells();
  //! [CppSnippetUMeshStdBuild1_3]
  //! [CppSnippetUMeshStdBuild1_4]
  MEDCoupling::DataArrayDouble *coordsArr=MEDCoupling::DataArrayDouble::New();
  coordsArr->alloc(9,3);//here coordsArr are declared to have 3 components, mesh will deduce that its spaceDim==3. 
  std::copy(coords,coords+27,coordsArr->getPointer());
  mesh->setCoords(coordsArr);//coordsArr contains 9 tuples, that is to say mesh contains 9 nodes.
  coordsArr->decrRef();
  //! [CppSnippetUMeshStdBuild1_4]
  mesh->checkConsistencyLight();
  //! [CppSnippetUMeshStdBuild1_5]
  mesh->decrRef();
  //! [CppSnippetUMeshStdBuild1_5]
}

void CppSnippetCMeshStdBuild1()
{
  //! [CppSnippetCMeshStdBuild1_1]
  double XCoords[9]={-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.,0.1,0.37,0.45,0.47,0.49,1.007};
  MEDCoupling::DataArrayDouble *arrX=MEDCoupling::DataArrayDouble::New();
  arrX->alloc(9,1);
  std::copy(XCoords,XCoords+9,arrX->getPointer());
  arrX->setInfoOnComponent(0,"X [m]");
  MEDCoupling::DataArrayDouble *arrY=MEDCoupling::DataArrayDouble::New();
  arrY->alloc(7,1);
  std::copy(YCoords,YCoords+7,arrY->getPointer());
  arrY->setInfoOnComponent(0,"Y [m]");
  //! [CppSnippetCMeshStdBuild1_1]
  //! [CppSnippetCMeshStdBuild1_2]
  MEDCoupling::MEDCouplingCMesh *mesh=MEDCoupling::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY);
  arrX->decrRef();
  arrY->decrRef();
  //! [CppSnippetCMeshStdBuild1_2]
  //! [CppSnippetCMeshStdBuild1_3]
  CPPUNIT_ASSERT_EQUAL(8*6,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9*7,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  //! [CppSnippetCMeshStdBuild1_3]
  mesh->decrRef();
  mesh=MEDCoupling::MEDCouplingCMesh::New("My2D_CMesh");
  arrX=MEDCoupling::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  arrY=MEDCoupling::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]");
  //! [CppSnippetCMeshStdBuild1_2bis]
  mesh->setCoordsAt(0,arrX);
  arrX->decrRef();
  mesh->setCoordsAt(1,arrY);
  arrY->decrRef();
  //! [CppSnippetCMeshStdBuild1_2bis]
  CPPUNIT_ASSERT_EQUAL(8*6,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(9*7,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  //! [CppSnippetCMeshStdBuild1_4]
  mesh->decrRef();
  //! [CppSnippetCMeshStdBuild1_4]
}

void CppSnippetUMeshAdvBuild1()
{
  //! [CppSnippetUMeshAdvBuild1_1]
  double coords[27]={-0.3,-0.3,0.,   0.2,-0.3,0.,   0.7,-0.3,0.,   -0.3,0.2,0.,   0.2,0.2,0., 
                     0.7,0.2,0.,    -0.3,0.7,0.,    0.2,0.7,0.,     0.7,0.7,0. };
  int nodalConnPerCell[23]={4,0,3,4,1, 3,1,4,2, 3,4,5,2, 4,6,7,4,3, 4,7,8,5,4};
  int nodalConnPerCellIndex[6]={0,5,9,13,18,23};
  //! [CppSnippetUMeshAdvBuild1_1]
  //! [CppSnippetUMeshAdvBuild1_2]
  MEDCoupling::MEDCouplingUMesh *mesh=MEDCoupling::MEDCouplingUMesh::New("My2DMesh",2);
  //! [CppSnippetUMeshAdvBuild1_2]
  //! [CppSnippetUMeshAdvBuild1_3]
  MEDCoupling::DataArrayInt *nodalConn=MEDCoupling::DataArrayInt::New();
  nodalConn->alloc(23,1);
  std::copy(nodalConnPerCell,nodalConnPerCell+23,nodalConn->getPointer());
  MEDCoupling::DataArrayInt *nodalConnI=MEDCoupling::DataArrayInt::New();
  nodalConnI->alloc(6,1);
  std::copy(nodalConnPerCellIndex,nodalConnPerCellIndex+6,nodalConnI->getPointer());
  mesh->setConnectivity(nodalConn,nodalConnI,true);
  nodalConn->decrRef();// nodalConn DataArrayInt instance is owned by mesh after call to setConnectivity method. No more need here -> decrRef()
  nodalConnI->decrRef();// nodalConnI DataArrayInt instance is owned by mesh after call to setConnectivity method. No more need here -> decrRef()
  //! [CppSnippetUMeshAdvBuild1_3]
  //! [CppSnippetUMeshAdvBuild1_4]
  MEDCoupling::DataArrayDouble *coordsArr=MEDCoupling::DataArrayDouble::New();
  coordsArr->alloc(9,3);//here coordsArr are declared to have 3 components, mesh will deduce that its spaceDim==3. 
  std::copy(coords,coords+27,coordsArr->getPointer());
  mesh->setCoords(coordsArr);//coordsArr contains 9 tuples, that is to say mesh contains 9 nodes.
  coordsArr->decrRef();
  //! [CppSnippetUMeshAdvBuild1_4]
  mesh->checkConsistencyLight();
  //! [CppSnippetUMeshAdvBuild1_5]
  mesh->decrRef();
  //! [CppSnippetUMeshAdvBuild1_5]
}

void CppSnippetDataArrayBuild1()
{
  //! [CppSnippetDataArrayBuild1_0]
  const int nbOfNodes=12;
  double coords[3*nbOfNodes]={2.,3.,4.,3.,4.,5.,4.,5.,6.,5.,6.,7.,6.,7.,8.,7.,8.,9.,8.,9.,10.,9.,10.,11.,10.,11.,12.,11.,12.,13.,12.,13.,14.,13.,14.,15.};
  //
  MEDCoupling::DataArrayDouble *coordsArr=0;
  double *tmp=0;
  //! [CppSnippetDataArrayBuild1_0]
  //
  //! [CppSnippetDataArrayBuild1_1]
  coordsArr=MEDCoupling::DataArrayDouble::New();
  coordsArr->useArray(coords,false,MEDCoupling::CPP_DEALLOC,nbOfNodes,3);
  //now use coordsArr as you need
  //...
  //coordsArr is no more useful here : release it
  coordsArr->decrRef();
  //! [CppSnippetDataArrayBuild1_1]
  //! [CppSnippetDataArrayBuild1_2]
  coordsArr=MEDCoupling::DataArrayDouble::New();
  tmp=new double[3*nbOfNodes];
  std::copy(coords,coords+3*nbOfNodes,tmp);
  coordsArr->useArray(tmp,true,MEDCoupling::CPP_DEALLOC,nbOfNodes,3);
  //now use coordsArr as you need
  //...
  //coordsArr is no more useful, release it
  coordsArr->decrRef();
  //! [CppSnippetDataArrayBuild1_2]
  //! [CppSnippetDataArrayBuild1_3]
  coordsArr=MEDCoupling::DataArrayDouble::New();
  tmp=(double *)malloc(3*nbOfNodes*sizeof(double));
  std::copy(coords,coords+3*nbOfNodes,tmp);
  coordsArr->useArray(tmp,true,MEDCoupling::C_DEALLOC,nbOfNodes,3);
  //now use coordsArr as you need
  //...
  //coordsArr is no more useful here : release it
  coordsArr->decrRef();
  //! [CppSnippetDataArrayBuild1_3]
  //! [CppSnippetDataArrayBuild1_4]
  coordsArr=MEDCoupling::DataArrayDouble::New();
  coordsArr->alloc(nbOfNodes,3);
  tmp=coordsArr->getPointer();
  std::copy(coords,coords+3*nbOfNodes,tmp);
  coordsArr->declareAsNew();//you have modified data pointed by internal pointer notify object
  //now use coordsArr as you need
  //...
  //coordsArr is no more useful here : release it
  coordsArr->decrRef();
  //! [CppSnippetDataArrayBuild1_4]
  coordsArr=MEDCoupling::DataArrayDouble::New();
  coordsArr->alloc(nbOfNodes,3);
  tmp=coordsArr->getPointer();
  std::copy(coords,coords+3*nbOfNodes,tmp);
  MEDCoupling::DataArrayDouble *coordsArrCpy=0;
  //! [CppSnippetDataArrayBuild1_5]
  coordsArrCpy=coordsArr->deepCopy();
  //! [CppSnippetDataArrayBuild1_5]
  //! [CppSnippetDataArrayBuild1_6]
  CPPUNIT_ASSERT(coordsArrCpy->isEqual(*coordsArr,1e-12));
  coordsArrCpy->setIJ(0,0,1000.);
  CPPUNIT_ASSERT(!coordsArrCpy->isEqual(*coordsArr,1e-12));//coordsArrCpy only has been modified
  //! [CppSnippetDataArrayBuild1_6]
  //! [CppSnippetDataArrayBuild1_7]
  coordsArrCpy->decrRef();
  //! [CppSnippetDataArrayBuild1_7]
  //! [CppSnippetDataArrayBuild1_5bis]
  coordsArrCpy=coordsArr->performCopyOrIncrRef(true);
  //! [CppSnippetDataArrayBuild1_5bis]
  CPPUNIT_ASSERT(coordsArrCpy->isEqual(*coordsArr,1e-12));
  coordsArrCpy->setIJ(0,0,1000.);
  CPPUNIT_ASSERT(!coordsArrCpy->isEqual(*coordsArr,1e-12));//coordsArrCpy only has been modified
  coordsArrCpy->decrRef();
  //! [CppSnippetDataArrayBuild1_8]
  coordsArrCpy=coordsArr->performCopyOrIncrRef(false);
  //! [CppSnippetDataArrayBuild1_8]
  //! [CppSnippetDataArrayBuild1_9]
  CPPUNIT_ASSERT(coordsArrCpy->isEqual(*coordsArr,1e-12));
  coordsArrCpy->setIJ(0,0,1000.);
  CPPUNIT_ASSERT(coordsArrCpy->isEqual(*coordsArr,1e-12));//coordsArr and coordsArrCpy have been modified simultaneously
  //! [CppSnippetDataArrayBuild1_9]
  //! [CppSnippetDataArrayBuild1_10]
  coordsArrCpy->decrRef();
  //! [CppSnippetDataArrayBuild1_10]
  //! [CppSnippetDataArrayBuild1_11]
  coordsArrCpy=MEDCoupling::DataArrayDouble::New();
  //! [CppSnippetDataArrayBuild1_11]
  //! [CppSnippetDataArrayBuild1_12]
  coordsArrCpy->deepCopyFrom(*coordsArr);
  //! [CppSnippetDataArrayBuild1_12]
  //! [CppSnippetDataArrayBuild1_13]
  CPPUNIT_ASSERT(coordsArrCpy->isEqual(*coordsArr,1e-12));
  coordsArrCpy->setIJ(0,0,2000.);
  CPPUNIT_ASSERT(!coordsArrCpy->isEqual(*coordsArr,1e-12));//coordsArrCpy only has been modified
  //! [CppSnippetDataArrayBuild1_13]
  //! [CppSnippetDataArrayBuild1_14]
  coordsArrCpy->decrRef();
  //! [CppSnippetDataArrayBuild1_14]
  coordsArr->decrRef();
  //! [CppSnippetDataArrayBuild1_14]
}

void CppSnippetFieldDoubleBuild1()
{
  double XCoords[9]={-0.3,0.07,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.07,0.1,0.37,0.45,0.47,0.49,1.007};
  MEDCoupling::DataArrayDouble *arrX=MEDCoupling::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  MEDCoupling::DataArrayDouble *arrY=MEDCoupling::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  MEDCoupling::MEDCouplingCMesh *mesh=MEDCoupling::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild1_1]
  MEDCoupling::MEDCouplingFieldDouble* fieldOnCells=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::NO_TIME);
  fieldOnCells->setName("MyTensorFieldOnCellNoTime");
  fieldOnCells->setMesh(mesh);
  mesh->decrRef(); // no more need of mesh because mesh has been attached to fieldOnCells
  MEDCoupling::DataArrayDouble *array=MEDCoupling::DataArrayDouble::New();
  array->alloc(fieldOnCells->getMesh()->getNumberOfCells(),9);//Implicitely fieldOnCells will be a 9 components field.
  array->fillWithValue(7.);
  fieldOnCells->setArray(array);
  array->decrRef();
  // fieldOnCells is now usable
  // ...
  // fieldOnCells is no more useful here : release it
  fieldOnCells->decrRef();
  //! [CppSnippetFieldDoubleBuild1_1]
  arrX=MEDCoupling::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  arrY=MEDCoupling::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  mesh=MEDCoupling::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild1_2]
  MEDCoupling::MEDCouplingFieldDouble *f1=mesh->fillFromAnalytic(MEDCoupling::ON_CELLS,1,"x*x+y*y*3+2.*x");//f1 is scalar
  MEDCoupling::MEDCouplingFieldDouble *f2=mesh->fillFromAnalytic(MEDCoupling::ON_CELLS,1,"cos(x+y/x)");//f2 is scalar too
  MEDCoupling::MEDCouplingFieldDouble *f2bis=mesh->fillFromAnalytic(MEDCoupling::ON_CELLS,2,"x*x*IVec+3*y*JVec");//f2bis is a vectors field
  MEDCoupling::MEDCouplingFieldDouble *f3=(*f1)+(*f2);//f3 scalar
  MEDCoupling::MEDCouplingFieldDouble *f4=(*f3)/(*f2);//f4 scalar
  f2bis->applyFunc(1,"sqrt(x*x+y*y)");//f2bis becomes scalar
  MEDCoupling::MEDCouplingFieldDouble *f5=(*f2bis)*(*f4);//f5 scalar
  const double pos1[2]={0.48,0.38};
  double res;
  f4->getValueOn(pos1,&res);//f4 is scalar so the returned value is of size 1.
  // ...
  //! [CppSnippetFieldDoubleBuild1_2]
  mesh->decrRef();
  //! [CppSnippetFieldDoubleBuild1_3]
  // f1, f2, f2bis, f3, f4, f5 are no more useful here : release them
  f1->decrRef();
  f2->decrRef();
  f2bis->decrRef();
  f3->decrRef();
  f4->decrRef();
  f5->decrRef();
  //! [CppSnippetFieldDoubleBuild1_3]
}

void CppSnippetFieldDoubleBuild2()
{
  double XCoords[9]={-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.,0.1,0.37,0.45,0.47,0.49,1.007};
  MEDCoupling::DataArrayDouble *arrX=MEDCoupling::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  MEDCoupling::DataArrayDouble *arrY=MEDCoupling::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  MEDCoupling::MEDCouplingCMesh *mesh=MEDCoupling::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild2_1]
  MEDCoupling::MEDCouplingFieldDouble* fieldOnNodes=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES,MEDCoupling::NO_TIME);
  fieldOnNodes->setName("MyScalarFieldOnNodeNoTime");
  fieldOnNodes->setMesh(mesh);
  mesh->decrRef(); // no more need of mesh because mesh has been attached to fieldOnNodes
  MEDCoupling::DataArrayDouble *array=MEDCoupling::DataArrayDouble::New();
  array->alloc(fieldOnNodes->getMesh()->getNumberOfNodes(),1);//Implicitely fieldOnNodes will be a 1 component field.
  array->fillWithValue(8.);
  fieldOnNodes->setArray(array);
  array->decrRef();
  // fieldOnNodes is now usable
  // ...
  // fieldOnNodes is no more useful here : release it
  fieldOnNodes->decrRef();
  //! [CppSnippetFieldDoubleBuild2_1]
}

void CppSnippetFieldDoubleBuild3()
{
  double XCoords[9]={-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.,0.1,0.37,0.45,0.47,0.49,1.007};
  MEDCoupling::DataArrayDouble *arrX=MEDCoupling::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  MEDCoupling::DataArrayDouble *arrY=MEDCoupling::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  MEDCoupling::MEDCouplingCMesh *mesh=MEDCoupling::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild3_1]
  MEDCoupling::MEDCouplingFieldDouble* fieldOnCells=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_CELLS,MEDCoupling::ONE_TIME);
  fieldOnCells->setName("MyTensorFieldOnCellNoTime");
  fieldOnCells->setTimeUnit("ms"); // Time unit is ms.
  fieldOnCells->setTime(4.22,2,-1); // Time attached is 4.22 ms, iteration id is 2 and order id (or sub iteration id) is -1
  fieldOnCells->setMesh(mesh);
  mesh->decrRef(); // no more need of mesh because mesh has been attached to fieldOnCells
  MEDCoupling::DataArrayDouble *array=MEDCoupling::DataArrayDouble::New();
  array->alloc(fieldOnCells->getMesh()->getNumberOfCells(),2);//Implicitely fieldOnCells will be a 2 components field.
  array->fillWithValue(7.);
  fieldOnCells->setArray(array);
  array->decrRef();
  // fieldOnCells is now usable
  // ...
  // fieldOnCells is no more useful here : release it
  fieldOnCells->decrRef();
  //! [CppSnippetFieldDoubleBuild3_1]
}

void CppSnippetFieldDoubleBuild4()
{
  double XCoords[9]={-0.3,0.,0.1,0.3,0.45,0.47,0.49,1.,1.22};
  double YCoords[7]={0.,0.1,0.37,0.45,0.47,0.49,1.007};
  MEDCoupling::DataArrayDouble *arrX=MEDCoupling::DataArrayDouble::New(); arrX->alloc(9,1); std::copy(XCoords,XCoords+9,arrX->getPointer()); arrX->setInfoOnComponent(0,"X [m]");
  MEDCoupling::DataArrayDouble *arrY=MEDCoupling::DataArrayDouble::New(); arrY->alloc(7,1); std::copy(YCoords,YCoords+7,arrY->getPointer()); arrY->setInfoOnComponent(0,"Y [m]"); 
  MEDCoupling::MEDCouplingCMesh *mesh=MEDCoupling::MEDCouplingCMesh::New("My2D_CMesh");
  mesh->setCoords(arrX,arrY); arrX->decrRef(); arrY->decrRef();
  //! [CppSnippetFieldDoubleBuild4_1]
  MEDCoupling::MEDCouplingFieldDouble* fieldOnNodes=MEDCoupling::MEDCouplingFieldDouble::New(MEDCoupling::ON_NODES,MEDCoupling::CONST_ON_TIME_INTERVAL);
  fieldOnNodes->setName("MyVecFieldOnNodeWithConstTime");
  fieldOnNodes->setTimeUnit("ms"); // Time unit is ms.
  fieldOnNodes->setStartTime(4.22,2,-1);
  fieldOnNodes->setEndTime(6.44,4,-1); // fieldOnNodes is defined in interval [4.22 ms,6.44 ms] 
  fieldOnNodes->setMesh(mesh);
  mesh->decrRef(); // no more need of mesh because mesh has been attached to fieldOnNodes
  MEDCoupling::DataArrayDouble *array=MEDCoupling::DataArrayDouble::New();
  array->alloc(fieldOnNodes->getMesh()->getNumberOfNodes(),3);//Implicitely fieldOnNodes will be a 3 components field.
  array->fillWithValue(8.);
  fieldOnNodes->setArray(array);
  array->decrRef();
  // fieldOnNodes is now usable
  // ...
  // fieldOnNodes is no more useful here : release it
  fieldOnNodes->decrRef();
  //! [CppSnippetFieldDoubleBuild4_1]
}

int main(int argc, char *argv[])
{
  CppExample_MEDCouplingFieldDouble_WriteVTK();
  CppExample_MEDCouplingFieldDouble_MaxFields();
  CppExample_MEDCouplingFieldDouble_MergeFields();
  CppExample_MEDCouplingFieldDouble_substractInPlaceDM();
  CppExample_MEDCouplingFieldDouble_changeUnderlyingMesh();
  CppExample_MEDCouplingFieldDouble_applyFunc_same_nb_comp();
  CppExample_MEDCouplingFieldDouble_applyFunc3();
  CppExample_MEDCouplingFieldDouble_applyFunc2();
  CppExample_MEDCouplingFieldDouble_applyFunc();
  CppExample_MEDCouplingFieldDouble_applyFunc_val();
  CppExample_MEDCouplingFieldDouble_fillFromAnalytic3();
  CppExample_MEDCouplingFieldDouble_fillFromAnalytic2();
  CppExample_MEDCouplingFieldDouble_fillFromAnalytic();
  CppExample_MEDCouplingFieldDouble_fillFromAnalytic_c_func();
  CppExample_MEDCouplingFieldDouble_applyFunc_c_func();
  CppExample_MEDCouplingFieldDouble_getValueOn_time();
  CppExample_MEDCouplingFieldDouble_getValueOnMulti();
  CppExample_MEDCouplingFieldDouble_getValueOn();
  CppExample_MEDCouplingFieldDouble_getValueOnPos();
  CppExample_MEDCouplingFieldDouble_renumberNodes();
  CppExample_MEDCouplingFieldDouble_renumberCells();
  CppExample_MEDCouplingFieldDouble_buildNewTimeReprFromThis();
  CppExample_MEDCouplingMesh_fillFromAnalytic3();
  CppExample_MEDCouplingMesh_fillFromAnalytic2();
  CppExample_MEDCouplingMesh_fillFromAnalytic();
  CppExample_MEDCouplingCMesh_getCoordsAt();
  CppExample_MEDCouplingUMesh_areCellsIncludedIn();
  CppExample_MEDCouplingUMesh_findAndCorrectBadOriented3DExtrudedCells();
  CppExample_MEDCouplingUMesh_arePolyhedronsNotCorrectlyOriented();
  CppExample_MEDCouplingUMesh_are2DCellsNotCorrectlyOriented();
  CppExample_MEDCouplingUMesh_getCellsContainingPoints();
  CppExample_MEDCouplingUMesh_getCellsContainingPoint();
  CppExample_MEDCouplingUMesh_buildPartOrthogonalField();
  CppExample_MEDCouplingUMesh_getPartMeasureField();
  CppExample_MEDCouplingUMesh_getCellsInBoundingBox();
  CppExample_MEDCouplingUMesh_renumberNodesInConn();
  CppExample_MEDCouplingUMesh_renumberNodes();
  CppExample_MEDCouplingUMesh_findBoundaryNodes();
  CppExample_MEDCouplingUMesh_buildBoundaryMesh();
  CppExample_MEDCouplingUMesh_buildFacePartOfMySelfNode();
  CppExample_MEDCouplingUMesh_buildPartOfMySelfNode();
  CppExample_MEDCouplingUMesh_getCellIdsLyingOnNodes();
  CppExample_MEDCouplingUMesh_getCellIdsFullyIncludedInNodeIds();
  CppExample_MEDCouplingUMesh_buildPartOfMySelf();
  CppExample_MEDCouplingUMesh_mergeNodes();
  CppExample_MEDCouplingUMesh_zipConnectivityTraducer();
  CppExample_MEDCouplingUMesh_zipCoordsTraducer();
  CppExample_MEDCouplingUMesh_getNodeIdsInUse();
  CppExample_MEDCouplingUMesh_convertToPolyTypes();
  CppExample_MEDCouplingUMesh_buildDescendingConnectivity2();
  CppExample_MEDCouplingUMesh_buildDescendingConnectivity();
  CppExample_MEDCouplingUMesh_getReverseNodalConnectivity();
  CppExample_MEDCouplingUMesh_checkDeepEquivalWith();
  CppExample_MEDCouplingPointSet_scale();
  CppExample_MEDCouplingPointSet_translate();
  CppExample_MEDCouplingPointSet_rotate();
  CppExample_MEDCouplingPointSet_getBoundingBox();
  CppExample_MEDCouplingPointSet_getNodeIdsNearPoint();
  CppExample_MEDCouplingPointSet_getNodeIdsNearPoints();
  CppExample_MEDCouplingPointSet_findCommonNodes();
  CppExample_MEDCouplingPointSet_getCoordinatesOfNode();
  CppExample_DataArrayInt_buildPermutationArr();
  CppExample_DataArrayInt_invertArrayO2N2N2O();
  CppExample_DataArrayInt_invertArrayN2O2O2N();
  CppExample_DataArrayDouble_getIdsInRange();
  CppExample_DataArrayDouble_findCommonTuples();
  CppExample_DataArrayDouble_Meld1();
  CppExampleFieldDoubleBuildSubPart1();
  CppSnippetUMeshStdBuild1();
  CppSnippetUMeshAdvBuild1();
  CppSnippetDataArrayBuild1();
  CppSnippetCMeshStdBuild1();
  CppSnippetFieldDoubleBuild1();
  CppSnippetFieldDoubleBuild2();
  CppSnippetFieldDoubleBuild3();
  CppSnippetFieldDoubleBuild4();

  return 0;
}
