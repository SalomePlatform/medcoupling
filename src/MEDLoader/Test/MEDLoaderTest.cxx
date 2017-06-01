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

#include "MEDLoaderTest.hxx"
#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "TestInterpKernelUtils.hxx"  // getResourceFile()

#include <cmath>
#include <numeric>

using namespace MEDCoupling;

void MEDLoaderTest::testMesh1DRW()
{
  MEDCouplingUMesh *mesh=build1DMesh_1();
  mesh->checkConsistencyLight();
  WriteUMesh("file1.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=ReadUMeshFromFile("file1.med",mesh->getName().c_str(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh2DCurveRW()
{
  MEDCouplingUMesh *mesh=build2DCurveMesh_1();
  mesh->checkConsistencyLight();
  WriteUMesh("file2.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=ReadUMeshFromFile("file2.med",mesh->getName().c_str(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh2DRW()
{
  MEDCouplingUMesh *mesh=build2DMesh_1();
  mesh->checkConsistencyLight();
  WriteUMesh("file3.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=ReadUMeshFromFile("file3.med",mesh->getName().c_str(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh3DSurfRW()
{
  MEDCouplingUMesh *mesh=build3DSurfMesh_1();
  mesh->checkConsistencyLight();
  WriteUMesh("file4.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=ReadUMeshFromFile("file4.med",mesh->getName().c_str(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh3DRW()
{
  MEDCouplingUMesh *mesh=build3DMesh_1();
  mesh->checkConsistencyLight();
  WriteUMesh("file5.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=ReadUMeshFromFile("file5.med",mesh->getName().c_str(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

/*!
 * Most basic test : one and only one MEDCoupling field in a new file.
 */
void MEDLoaderTest::testFieldRW1()
{
  MEDCouplingFieldDouble *f1=buildVecFieldOnCells_1();
  WriteField("file6.med",f1,true);
  MEDCouplingFieldDouble *f2=ReadFieldCell("file6.med",f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),0,1);
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f1->decrRef();
  f2->decrRef();
  //
  f1=buildVecFieldOnNodes_1();
  WriteField("file7.med",f1,true);
  f2=ReadFieldNode("file7.med",f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),2,3);
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  // testing kind message on error of field type.
  CPPUNIT_ASSERT_THROW(ReadFieldCell("file7.med",f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),2,3),INTERP_KERNEL::Exception);
  //
  f1->decrRef();
  f2->decrRef();
}

/*!
 * Multi field writing in a same file.
 */
void MEDLoaderTest::testFieldRW2()
{
  const char fileName[]="file8.med";
  static const double VAL1=12345.67890314;
  static const double VAL2=-1111111111111.;
  MEDCouplingFieldDouble *f1=buildVecFieldOnCells_1();
  WriteField(fileName,f1,true);
  f1->setTime(10.,8,9);
  double *tmp=f1->getArray()->getPointer();
  tmp[0]=VAL1;
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(10.14,18,19);
  tmp[0]=VAL2;
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  //retrieving time steps...
  MEDCouplingFieldDouble *f2=ReadFieldCell(fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),8,9);
  f1->setTime(10.,8,9);
  tmp[0]=VAL1;
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f2->decrRef();
  f2=ReadFieldCell(fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),0,1);
  MEDCouplingFieldDouble *f3=buildVecFieldOnCells_1();
  CPPUNIT_ASSERT(f3->isEqual(f2,1e-12,1e-12));
  f3->decrRef();
  f2->decrRef();
  f2=ReadFieldCell(fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),18,19);
  f1->setTime(10.14,18,19);
  tmp[0]=VAL2;
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  //test of throw on invalid (dt,it)
  CPPUNIT_ASSERT_THROW(ReadFieldCell(fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),28,19),INTERP_KERNEL::Exception);
  f2->decrRef();
  f1->decrRef();
  //ON NODES
  f1=buildVecFieldOnNodes_1();
  const char fileName2[]="file9.med";
  WriteField(fileName2,f1,true);
  f1->setTime(110.,108,109);
  tmp=f1->getArray()->getPointer();
  tmp[3]=VAL1;
  WriteFieldUsingAlreadyWrittenMesh(fileName2,f1);
  f1->setTime(210.,208,209);
  tmp[3]=VAL2;
  WriteFieldUsingAlreadyWrittenMesh(fileName2,f1);
  f2=ReadFieldNode(fileName2,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),108,109);
  f1->setTime(110.,108,109);
  tmp[3]=VAL1;
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f2->decrRef();
  f2=ReadFieldNode(fileName2,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),2,3);
  f3=buildVecFieldOnNodes_1();
  CPPUNIT_ASSERT(f3->isEqual(f2,1e-12,1e-12));
  f3->decrRef();
  f2->decrRef();
  f2=ReadFieldNode(fileName2,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),208,209);
  f1->setTime(210.,208,209);
  tmp[3]=VAL2;
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f2->decrRef();
  f1->decrRef();
}

/*!
 * Multi field in a same file, but this field has several
 */
void MEDLoaderTest::testFieldRW3()
{
  const char fileName[]="file11.med";
  static const double VAL1=12345.67890314;
  static const double VAL2=-1111111111111.;
  const char name1[]="AField";
  const char name3[]="AMesh1";
  MEDCouplingFieldDouble *f1=buildVecFieldOnCells_1();
  (const_cast<MEDCouplingMesh *>(f1->getMesh()))->setName(name3);
  f1->setName(name1);
  f1->setTime(10.,8,9);
  double *tmp=f1->getArray()->getPointer();
  tmp[0]=VAL1;
  WriteField(fileName,f1,true);
  f1->setTime(10.14,18,19);
  tmp[0]=VAL2;
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(10.55,28,29);
  tmp[0]=3*VAL1;
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(10.66,38,39);
  tmp[0]=3*VAL2;
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(10.77,48,49);
  tmp[0]=4*VAL2;
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  //ON NODES
  f1->decrRef();
  f1=buildVecFieldOnNodes_1();
  f1->setName(name1);
  (const_cast<MEDCouplingMesh *>(f1->getMesh()))->setName(name3);
  f1->setTime(110.,8,9);
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(110.,108,109);
  tmp=f1->getArray()->getPointer();
  tmp[3]=VAL1;
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(210.,208,209);
  tmp[3]=VAL2;
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  //
  std::vector< std::pair<int,int> > it1=GetCellFieldIterations(fileName,name3,name1);
  CPPUNIT_ASSERT_EQUAL(5,(int)it1.size());
  CPPUNIT_ASSERT_EQUAL(8,it1[0].first); CPPUNIT_ASSERT_EQUAL(9,it1[0].second);
  CPPUNIT_ASSERT_EQUAL(18,it1[1].first); CPPUNIT_ASSERT_EQUAL(19,it1[1].second);
  CPPUNIT_ASSERT_EQUAL(28,it1[2].first); CPPUNIT_ASSERT_EQUAL(29,it1[2].second);
  CPPUNIT_ASSERT_EQUAL(38,it1[3].first); CPPUNIT_ASSERT_EQUAL(39,it1[3].second);
  CPPUNIT_ASSERT_EQUAL(48,it1[4].first); CPPUNIT_ASSERT_EQUAL(49,it1[4].second);
  std::vector< std::pair<int,int> > it3=GetNodeFieldIterations(fileName,name3,name1);
  CPPUNIT_ASSERT_EQUAL(3,(int)it3.size());
  CPPUNIT_ASSERT_EQUAL(8,it3[0].first); CPPUNIT_ASSERT_EQUAL(9,it3[0].second);
  CPPUNIT_ASSERT_EQUAL(108,it3[1].first); CPPUNIT_ASSERT_EQUAL(109,it3[1].second);
  CPPUNIT_ASSERT_EQUAL(208,it3[2].first); CPPUNIT_ASSERT_EQUAL(209,it3[2].second);
  //
  f1->decrRef();
  //
  f1=ReadFieldCell(fileName,name3,0,name1,8,9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(VAL1,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  f1=ReadFieldCell(fileName,name3,0,name1,18,19);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(VAL2,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  f1=ReadFieldCell(fileName,name3,0,name1,28,29);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3*VAL1,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  f1=ReadFieldCell(fileName,name3,0,name1,38,39);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3*VAL2,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  f1=ReadFieldCell(fileName,name3,0,name1,48,49);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4*VAL2,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  //
  f1=ReadFieldNode(fileName,name3,0,name1,8,9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(71.,f1->getArray()->getConstPointer()[3],1e-13);
  f1->decrRef();
  f1=ReadFieldNode(fileName,name3,0,name1,108,109);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(VAL1,f1->getArray()->getConstPointer()[3],1e-13);
  f1->decrRef();
  f1=ReadFieldNode(fileName,name3,0,name1,208,209);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(VAL2,f1->getArray()->getConstPointer()[3],1e-13);
  f1->decrRef();
}

void MEDLoaderTest::testMultiMeshRW1()
{
  const char fileName[]="file10.med";
  MEDCouplingUMesh *mesh1=build3DMesh_1();
  const int part1[5]={1,2,4,13,15};
  MEDCouplingUMesh *mesh2=(MEDCouplingUMesh *)mesh1->buildPartOfMySelf(part1,part1+5,true);
  mesh2->setName("mesh2");
  const int part2[4]={3,4,13,14};
  MEDCouplingUMesh *mesh3=(MEDCouplingUMesh *)mesh1->buildPartOfMySelf(part2,part2+4,true);
  mesh3->setName("mesh3");
  MEDCouplingUMesh *mesh4=MEDCouplingUMesh::New();
  mesh4->setName("mesh4");
  mesh4->setMeshDimension(3);
  mesh4->allocateCells(1);
  int conn[4]={0,11,1,3};
  mesh4->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4,conn);
  mesh4->finishInsertingCells();
  mesh4->setCoords(mesh1->getCoords());
  std::vector<const MEDCouplingUMesh *> meshes;
  meshes.push_back(mesh1);
  meshes.push_back(mesh2);
  meshes.push_back(mesh3);
  meshes.push_back(mesh4);
  const char mnane[]="3DToto";
  WriteUMeshesPartition(fileName,mnane,meshes,true);
  //
  MEDCouplingUMesh *mesh5=ReadUMeshFromFile(fileName,mnane);
  mesh1->setName(mnane);
  const int part3[18]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17};
  MEDCouplingUMesh *mesh6=(MEDCouplingUMesh *)mesh5->buildPartOfMySelf(part3,part3+18,true);
  mesh6->setName(mnane);
  mesh5->decrRef();
  CPPUNIT_ASSERT(mesh6->isEqual(mesh1,1e-12));
  mesh6->decrRef();
  std::vector<std::string> grps=GetMeshGroupsNames(fileName,mnane);
  CPPUNIT_ASSERT_EQUAL(4,(int)grps.size());
  CPPUNIT_ASSERT(std::find(grps.begin(),grps.end(),std::string("mesh2"))!=grps.end());
  CPPUNIT_ASSERT(std::find(grps.begin(),grps.end(),std::string("mesh3"))!=grps.end());
  CPPUNIT_ASSERT(std::find(grps.begin(),grps.end(),std::string("mesh4"))!=grps.end());
  CPPUNIT_ASSERT(std::find(grps.begin(),grps.end(),std::string("3DMesh_1"))!=grps.end());
  //
  std::vector<std::string> vec;
  vec.push_back(std::string("mesh2"));
  MEDCouplingUMesh *mesh2_2=ReadUMeshFromGroups(fileName,mnane,0,vec);
  CPPUNIT_ASSERT(mesh2_2->isEqual(mesh2,1e-12));
  mesh2_2->decrRef();
  vec.clear(); vec.push_back(std::string("mesh3"));
  MEDCouplingUMesh *mesh3_2=ReadUMeshFromGroups(fileName,mnane,0,vec);
  CPPUNIT_ASSERT(mesh3_2->isEqual(mesh3,1e-12));
  mesh3_2->decrRef();
  vec.clear(); vec.push_back(std::string("mesh4"));
  MEDCouplingUMesh *mesh4_2=ReadUMeshFromGroups(fileName,mnane,0,vec);
  CPPUNIT_ASSERT(mesh4_2->isEqual(mesh4,1e-12));
  mesh4_2->decrRef();
  vec.clear(); vec.push_back(std::string("3DMesh_1"));
  MEDCouplingUMesh *mesh1_2=ReadUMeshFromGroups(fileName,mnane,0,vec);
  mesh1->setName("3DMesh_1");
  CPPUNIT_ASSERT(mesh1_2->isEqual(mesh1,1e-12));
  mesh1_2->decrRef();
  //
  vec.clear(); vec.push_back(std::string("Family_-3")); vec.push_back(std::string("Family_-5"));
  mesh2_2=ReadUMeshFromFamilies(fileName,mnane,0,vec);
  mesh2_2->setName("mesh2");
  CPPUNIT_ASSERT(mesh2_2->isEqual(mesh2,1e-12));
  mesh2_2->decrRef();
  //
  std::vector<std::string> ret(GetMeshFamiliesNamesOnGroup(fileName,"3DToto","3DMesh_1"));
  std::set<std::string> s(ret.begin(),ret.end());
  std::set<std::string> ref_s;
  ref_s.insert("Family_-2");
  ref_s.insert("Family_-3");
  ref_s.insert("Family_-4");
  ref_s.insert("Family_-5");
  CPPUNIT_ASSERT_EQUAL(4,(int)ret.size());
  CPPUNIT_ASSERT(s==ref_s);
  //
  std::vector<std::string> ret1=GetMeshGroupsNamesOnFamily(fileName,"3DToto","Family_-3");
  CPPUNIT_ASSERT_EQUAL(2,(int)ret1.size());
  CPPUNIT_ASSERT(ret1[0]=="3DMesh_1");
  CPPUNIT_ASSERT(ret1[1]=="mesh2");
  //
  mesh4->decrRef();
  mesh3->decrRef();
  mesh2->decrRef();
  mesh1->decrRef();
}

void MEDLoaderTest::testFieldProfilRW1()
{
  const char fileName[]="file12.med";
  MEDCouplingUMesh *mesh1=build3DMesh_1();
  bool b;
  int newNbOfNodes;
  DataArrayInt *da=mesh1->mergeNodes(1e-12,b,newNbOfNodes);
  da->decrRef();
  WriteUMesh(fileName,mesh1,true);
  const int part1[5]={1,2,4,13,15};
  MEDCouplingUMesh *mesh2=(MEDCouplingUMesh *)mesh1->buildPartOfMySelf(part1,part1+5,true);
  mesh2->setName(mesh1->getName().c_str());//<- important for the test
  //
  int nbOfCells=mesh2->getNumberOfCells();
  CPPUNIT_ASSERT_EQUAL(5,nbOfCells);
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setName("VectorFieldOnCells");
  f1->setMesh(mesh2);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(nbOfCells,2);
  f1->setArray(array);
  array->decrRef();
  double *tmp=array->getPointer();
  const double arr1[10]={71.,171.,10.,110.,20.,120.,30.,130.,40.,140.};
  std::copy(arr1,arr1+10,tmp);
  f1->setTime(3.14,2,7);
  f1->checkConsistencyLight();
  //
  WriteField(fileName,f1,false);//<- false important for the test
  //
  MEDCouplingFieldDouble *f2=ReadFieldCell(fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),2,7);
  std::vector<MEDCoupling::TypeOfField> types=GetTypesOfField(fileName,f1->getMesh()->getName().c_str(),f1->getName().c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)types.size());
  CPPUNIT_ASSERT(types[0]==ON_CELLS);
  f2->checkConsistencyLight();
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  //
  f2->decrRef();
  f1->decrRef();
  mesh1->decrRef();
  mesh2->decrRef();
}

/*!
 * Test MED file profiles.
 */
void MEDLoaderTest::testFieldNodeProfilRW1()
{
  const char fileName[]="file19.med";
  const char fileName2[]="file20.med";
  MEDCouplingUMesh *m=build2DMesh_1();
  int nbOfNodes=m->getNumberOfNodes();
  WriteUMesh(fileName,m,true);
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  f1->setName("VFieldOnNodes");
  f1->setMesh(m);
  DataArrayDouble *array=DataArrayDouble::New();
  const double arr1[24]={1.,101.,2.,102.,3.,103.,4.,104.,5.,105.,6.,106.,7.,107.,8.,108.,9.,109.,10.,110.,11.,111.,12.,112.};
  array->alloc(nbOfNodes,2);
  std::copy(arr1,arr1+24,array->getPointer());
  f1->setArray(array);
  array->setInfoOnComponent(0,"tyty [mm]");
  array->setInfoOnComponent(1,"uiop [MW]");
  array->decrRef();
  f1->setTime(3.14,2,7);
  f1->checkConsistencyLight();
  const int arr2[2]={1,4};//node ids are 2,4,5,3,6,7
  MEDCouplingFieldDouble *f2=f1->buildSubPart(arr2,arr2+2);
  (const_cast<MEDCouplingMesh *>(f2->getMesh()))->setName(f1->getMesh()->getName().c_str());
  WriteField(fileName,f2,false);//<- false important for the test
  //
  MEDCouplingFieldDouble *f3=ReadFieldNode(fileName,f2->getMesh()->getName().c_str(),0,f2->getName().c_str(),2,7);
  f3->checkConsistencyLight();
  CPPUNIT_ASSERT(f3->isEqual(f2,1e-12,1e-12));
  f3->decrRef();
  //
  const int arr3[6]={1,3,0,5,2,4};
  f2->renumberNodes(arr3);
  WriteUMesh(fileName2,m,true);
  WriteField(fileName2,f2,false);//<- false important for the test
  f3=ReadFieldNode(fileName2,f2->getMesh()->getName().c_str(),0,f2->getName().c_str(),2,7);
  f3->checkConsistencyLight();
  CPPUNIT_ASSERT(f3->isEqual(f2,1e-12,1e-12));
  f3->decrRef();
  f2->decrRef();
  //
  f1->decrRef();
  m->decrRef();
}

void MEDLoaderTest::testFieldNodeProfilRW2()
{
  const char fileName[]="file23.med";
  MEDCouplingUMesh *mesh=build3DSurfMesh_1();
  WriteUMesh(fileName,mesh,true);
  //
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  f1->setName("FieldMix");
  f1->setMesh(mesh);
  const double arr2[24]={
    1071.,1171.,1010.,1110.,1020.,1120.,1030.,1130.,1040.,1140.,1050.,1150.,
    1060.,1160.,1070.,1170.,1080.,1180.,1090.,1190.,1091.,1191.,1092.,1192.
  };
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(12,2);
  f1->setArray(array);
  array->setInfoOnComponent(0,"plkj [mm]");
  array->setInfoOnComponent(1,"pqqqss [mm]");
  array->decrRef();
  double *tmp=array->getPointer();
  std::copy(arr2,arr2+24,tmp);
  f1->setTime(3.17,2,7);
  //
  const int renumArr[12]={3,7,2,1,5,11,10,0,9,6,8,4};
  f1->renumberNodes(renumArr);
  f1->checkConsistencyLight();
  WriteField(fileName,f1,false);//<- false important for the test
  MEDCouplingFieldDouble *f2=ReadFieldNode(fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),2,7);
  CPPUNIT_ASSERT(f2->isEqual(f1,1e-12,1e-12));
  //
  f2->decrRef();
  mesh->decrRef();
  f1->decrRef();
}

void MEDLoaderTest::testFieldGaussRW1()
{
  const char fileName[]="file13.med";
  MEDCouplingFieldDouble *f1=buildVecFieldOnGauss_1();
  WriteField(fileName,f1,true);
  MCAuto<MEDCouplingFieldDouble> f2(ReadField(ON_GAUSS_PT,fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),1,5));
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f1->decrRef();
}

void MEDLoaderTest::testFieldGaussNERW1()
{
  const char fileName[]="file14.med";
  MEDCouplingFieldDouble *f1=buildVecFieldOnGaussNE_1();
  WriteField(fileName,f1,true);
  std::vector<MEDCoupling::TypeOfField> tof(GetTypesOfField(fileName,"2DMesh_2","MyFieldOnGaussNE"));
  CPPUNIT_ASSERT_EQUAL(1,(int)tof.size());
  CPPUNIT_ASSERT(ON_GAUSS_NE==tof[0]);
  MCAuto<MEDCouplingFieldDouble> f2(ReadField(ON_GAUSS_NE,fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),1,5));
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f1->decrRef();
}

void MEDLoaderTest::testLittleStrings1()
{
  std::string s("azeeeerrrtty");
  MEDLoaderBase::zipEqualConsChar(s,3);
  CPPUNIT_ASSERT(s=="azertty");
}

void MEDLoaderTest::testSplitIntoNameAndUnit1()
{
  std::string s(" []");
  std::string c,u;
  MEDLoaderBase::splitIntoNameAndUnit(s,c,u);
  CPPUNIT_ASSERT(c.empty());
  CPPUNIT_ASSERT(u.empty());
  s="   lmmm  kki jjj      ";
  MEDLoaderBase::strip(s);
  CPPUNIT_ASSERT(s=="lmmm  kki jjj");
  s=" ";
  MEDLoaderBase::strip(s);
  CPPUNIT_ASSERT(s.empty());
  s="";
  MEDLoaderBase::strip(s);
  CPPUNIT_ASSERT(s.empty());
  s="      ";
  MEDLoaderBase::strip(s);
  CPPUNIT_ASSERT(s.empty());
  s="     pp";
  MEDLoaderBase::strip(s);
  CPPUNIT_ASSERT(s=="pp");
}

void MEDLoaderTest::testMesh3DSurfShuffleRW()
{
  const char fileName[]="file15.med";
  MEDCouplingUMesh *mesh=build3DSurfMesh_1();
  const int renumber1[6]={2,5,1,0,3,4};
  mesh->renumberCells(renumber1,false);
  mesh->checkConsistencyLight();
  WriteUMesh(fileName,mesh,true);
  MEDCouplingUMesh *mesh_rw=ReadUMeshFromFile(fileName,mesh->getName().c_str(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testFieldShuffleRW1()
{
  const char fileName[]="file16.med";
  MEDCouplingUMesh *mesh=build3DSurfMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setName("FieldOnCellsShuffle");
  f1->setMesh(mesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(6,2);
  f1->setArray(array);
  array->decrRef();
  double *tmp=array->getPointer();
  const double arr1[12]={71.,171.,10.,110.,20.,120.,30.,130.,40.,140.,50.,150.};
  std::copy(arr1,arr1+12,tmp);
  f1->setTime(3.14,2,7);
  f1->checkConsistencyLight();
  //
  const int renumber1[6]={2,1,5,0,3,4};
  f1->renumberCells(renumber1,false);
  WriteField(fileName,f1,true);
  MEDCouplingFieldDouble *f2=ReadFieldCell(fileName,mesh->getName().c_str(),0,f1->getName().c_str(),2,7);
  CPPUNIT_ASSERT(f2->isEqual(f1,1e-12,1e-12));
  f2->decrRef();
  //
  mesh->decrRef();
  f1->decrRef();
}

/*!
 * Shuffle de cells but no profile. Like pointe.med
 */
void MEDLoaderTest::testMultiFieldShuffleRW1()
{
  const char fileName[]="file17.med";
  MEDCouplingUMesh *m=build3DMesh_2();
  CPPUNIT_ASSERT_EQUAL(20,(int)m->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(45,m->getNumberOfNodes());
  const int polys[3]={1,4,6};
  std::vector<int> poly2(polys,polys+3);
  m->convertToPolyTypes(&poly2[0],&poly2[0]+poly2.size());
  const int renum[20]={1,3,2,8,9,12,13,16,19,0,4,7,5,15,14,17,10,18,6,11};
  m->renumberCells(renum,false);
  m->orientCorrectlyPolyhedrons();
  // Writing
  WriteUMesh(fileName,m,true);
  MEDCouplingFieldDouble *f1Tmp=m->getMeasureField(false);
  MEDCouplingFieldDouble *f1=f1Tmp->buildNewTimeReprFromThis(ONE_TIME,false);
  f1Tmp->decrRef();
  f1->setTime(0.,1,2);
  MEDCouplingFieldDouble *f_1=f1->cloneWithMesh(true);
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->applyFunc("2*x");
  f1->setTime(0.01,3,4);
  MEDCouplingFieldDouble *f_2=f1->cloneWithMesh(true);
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->applyFunc("2*x/3");
  f1->setTime(0.02,5,6);
  MEDCouplingFieldDouble *f_3=f1->cloneWithMesh(true);
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->decrRef();
  // Reading
  std::vector<std::pair<int,int> > its;
  its.push_back(std::pair<int,int>(1,2));
  its.push_back(std::pair<int,int>(3,4));
  its.push_back(std::pair<int,int>(5,6));
  std::vector<MEDCouplingFieldDouble *> fs=ReadFieldsOnSameMesh(ON_CELLS,fileName,f_1->getMesh()->getName().c_str(),0,f_1->getName().c_str(),its);
  CPPUNIT_ASSERT_EQUAL(3,(int)fs.size());
  const MEDCouplingMesh *mm=fs[0]->getMesh();
  CPPUNIT_ASSERT(fs[0]->isEqual(f_1,1e-12,1e-12));
  CPPUNIT_ASSERT(fs[1]->isEqual(f_2,1e-12,1e-12));
  CPPUNIT_ASSERT(fs[2]->isEqual(f_3,1e-12,1e-12));
  CPPUNIT_ASSERT(mm==fs[1]->getMesh());// <- important for the test
  CPPUNIT_ASSERT(mm==fs[2]->getMesh());// <- important for the test
  for(std::vector<MEDCouplingFieldDouble *>::iterator iter=fs.begin();iter!=fs.end();iter++)
    (*iter)->decrRef();
  //
  f_1->decrRef();
  f_2->decrRef();
  f_3->decrRef();
  //
  m->decrRef();
}

void MEDLoaderTest::testWriteUMeshesRW1()
{
  const char fileName[]="file18.med";
  MEDCouplingUMesh *m3d=build3DMesh_2();
  const double pt[3]={0.,0.,-0.3};
  const double vec[3]={0.,0.,1.};
  std::vector<int> nodes;
  m3d->findNodesOnPlane(pt,vec,1e-12,nodes);
  MEDCouplingUMesh *m2d=(MEDCouplingUMesh *)m3d->buildFacePartOfMySelfNode(&nodes[0],&nodes[0]+nodes.size(),true);
  const int renumber[5]={1,2,0,4,3};
  m2d->renumberCells(renumber,false);
  m2d->setName("ExampleOfMultiDimW");
  std::vector<const MEDCouplingUMesh *> meshes;
  meshes.push_back(m2d);
  meshes.push_back(m3d);
  WriteUMeshes(fileName,meshes,true);
  MEDCouplingUMesh *m3d_bis=ReadUMeshFromFile(fileName,m2d->getName().c_str(),0);
  CPPUNIT_ASSERT(!m3d_bis->isEqual(m3d,1e-12));
  m3d_bis->setName(m3d->getName().c_str());
  CPPUNIT_ASSERT(m3d_bis->isEqual(m3d,1e-12));
  MEDCouplingUMesh *m2d_bis=ReadUMeshFromFile(fileName,m2d->getName().c_str(),-1);//-1 for faces
  CPPUNIT_ASSERT(m2d_bis->isEqual(m2d,1e-12));
  // Creation of a field on faces.
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setName("FieldOnFacesShuffle");
  f1->setMesh(m2d);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(m2d->getNumberOfCells(),2);
  array->setInfoOnComponent(0,"plkj [mm]");
  array->setInfoOnComponent(1,"pqqqss [mm]");
  f1->setArray(array);
  array->decrRef();
  double *tmp=array->getPointer();
  const double arr1[10]={71.,171.,10.,110.,20.,120.,30.,130.,40.,140.};
  std::copy(arr1,arr1+10,tmp);
  f1->setTime(3.14,2,7);
  f1->checkConsistencyLight();
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  MEDCouplingFieldDouble *f2=ReadFieldCell(fileName,f1->getMesh()->getName().c_str(),-1,f1->getName().c_str(),2,7);
  CPPUNIT_ASSERT(f2->isEqual(f1,1e-12,1e-12));
  f1->decrRef();
  f2->decrRef();
  //
  m2d_bis->decrRef();
  m3d_bis->decrRef();
  m2d->decrRef();
  m3d->decrRef();
}

void MEDLoaderTest::testMixCellAndNodesFieldRW1()
{
  const char fileName[]="file21.med";
  MEDCouplingUMesh *mesh=build3DSurfMesh_1();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setName("FieldMix");
  f1->setMesh(mesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(6,2);
  f1->setArray(array);
  array->setInfoOnComponent(0,"plkj [mm]");
  array->setInfoOnComponent(1,"pqqqss [mm]");
  array->decrRef();
  double *tmp=array->getPointer();
  const double arr1[12]={71.,171.,10.,110.,20.,120.,30.,130.,40.,140.,50.,150.};
  std::copy(arr1,arr1+12,tmp);
  f1->setTime(3.14,2,7);
  f1->checkConsistencyLight();
  //
  MEDCouplingFieldDouble *f2=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  f2->setName("FieldMix");
  f2->setMesh(mesh);
  array=DataArrayDouble::New();
  array->alloc(12,2);
  f2->setArray(array);
  array->setInfoOnComponent(0,"plkj [mm]");
  array->setInfoOnComponent(1,"pqqqss [mm]");
  array->decrRef();
  tmp=array->getPointer();
  const double arr2[24]={
    1071.,1171.,1010.,1110.,1020.,1120.,1030.,1130.,1040.,1140.,1050.,1150.,
    1060.,1160.,1070.,1170.,1080.,1180.,1090.,1190.,1091.,1191.,1092.,1192.
  };
  std::copy(arr2,arr2+24,tmp);
  f2->setTime(3.14,2,7);
  f2->checkConsistencyLight();
  //
  WriteField(fileName,f1,true);
  std::vector<MEDCoupling::TypeOfField> ts=GetTypesOfField(fileName,f1->getMesh()->getName().c_str(),f1->getName().c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)ts.size());
  CPPUNIT_ASSERT_EQUAL(ON_CELLS,ts[0]);
  std::vector<std::string> fs=GetAllFieldNamesOnMesh(fileName,f1->getMesh()->getName().c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)fs.size());
  CPPUNIT_ASSERT(fs[0]=="FieldMix");
  WriteFieldUsingAlreadyWrittenMesh(fileName,f2);
  fs=GetAllFieldNamesOnMesh(fileName,f1->getMesh()->getName().c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)fs.size());
  CPPUNIT_ASSERT(fs[0]=="FieldMix");
  //
  ts=GetTypesOfField(fileName,f1->getMesh()->getName().c_str(),f1->getName().c_str());
  CPPUNIT_ASSERT_EQUAL(2,(int)ts.size());
  CPPUNIT_ASSERT_EQUAL(ON_NODES,ts[0]);
  CPPUNIT_ASSERT_EQUAL(ON_CELLS,ts[1]);
  //
  MEDCouplingFieldDouble *f3=ReadFieldNode(fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),2,7);
  CPPUNIT_ASSERT(f3->isEqual(f2,1e-12,1e-12));
  f3->decrRef();
  f3=ReadFieldCell(fileName,f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),2,7);
  CPPUNIT_ASSERT(f3->isEqual(f1,1e-12,1e-12));
  f3->decrRef();
  //
  f1->decrRef();
  f2->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testGetAllFieldNamesRW1()
{
  const char fileName[]="file22.med";
  MEDCouplingUMesh *mesh=build2DMesh_2();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  f1->setName("Field1");
  f1->setTime(3.44,5,6);
  f1->setMesh(mesh);
  f1->fillFromAnalytic(2,"x+y");
  WriteField(fileName,f1,true);
  f1->setTime(1002.3,7,8);
  f1->fillFromAnalytic(2,"x+77.*y");
  WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setName("Field2");
  WriteField(fileName,f1,false);
  f1->setName("Field3");
  mesh->setName("2DMesh_2Bis");
  WriteField(fileName,f1,false);
  f1->decrRef();
  f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setName("Field8");
  f1->setTime(8.99,7,9);
  f1->setMesh(mesh);
  f1->fillFromAnalytic(3,"3*x+y");
  WriteField(fileName,f1,false);
  f1->decrRef();
  std::vector<std::string> fs=GetAllFieldNames(fileName);
  CPPUNIT_ASSERT_EQUAL(4,(int)fs.size());
  CPPUNIT_ASSERT(fs[0]=="Field1");
  CPPUNIT_ASSERT(fs[1]=="Field2");
  CPPUNIT_ASSERT(fs[2]=="Field3");
  CPPUNIT_ASSERT(fs[3]=="Field8");
  mesh->decrRef();
}


void MEDLoaderTest::testMEDLoaderRead1()
{
  using namespace std;
  using namespace INTERP_KERNEL;

  string fileName= INTERP_TEST::getResourceFile("pointe.med", 3);
  vector<string> meshNames=GetMeshNames(fileName.c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)meshNames.size());
  MEDCouplingUMesh *mesh=ReadUMeshFromFile(fileName.c_str(),meshNames[0].c_str(),0);
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)mesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,mesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,mesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,mesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,mesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,mesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  mesh->decrRef();
  //
  vector<string> families=GetMeshFamiliesNames(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(8,(int)families.size());
  CPPUNIT_ASSERT(families[2]=="FAMILLE_ELEMENT_3");
  //
  vector<string> families2;
  families2.push_back(families[2]);
  mesh=ReadUMeshFromFamilies(fileName.c_str(),meshNames[0].c_str(),0,families2);
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,mesh->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL((std::size_t)11,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(132,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+11,0));
  CPPUNIT_ASSERT_EQUAL(16,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+3,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  mesh->decrRef();
  //
  vector<string> groups=GetMeshGroupsNames(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(5,(int)groups.size());
  CPPUNIT_ASSERT(groups[0]=="groupe1");
  CPPUNIT_ASSERT(groups[1]=="groupe2");
  CPPUNIT_ASSERT(groups[2]=="groupe3");
  CPPUNIT_ASSERT(groups[3]=="groupe4");
  CPPUNIT_ASSERT(groups[4]=="groupe5");
  vector<string> groups2;
  groups2.push_back(groups[0]);
  mesh=ReadUMeshFromGroups(fileName.c_str(),meshNames[0].c_str(),0,groups2);
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(7,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  for(int i=0;i<6;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,mesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,mesh->getTypeOfCell(6));
  CPPUNIT_ASSERT_EQUAL((std::size_t)36,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(254,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+36,0));
  CPPUNIT_ASSERT_EQUAL(141,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+8,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  mesh->decrRef();
  //
  std::vector<std::string> fieldsName=GetCellFieldNamesOnMesh(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(2,(int)fieldsName.size());
  CPPUNIT_ASSERT(fieldsName[0]=="fieldcelldoublescalar");
  CPPUNIT_ASSERT(fieldsName[1]=="fieldcelldoublevector");
  std::vector<std::pair<int,int> > its0=GetCellFieldIterations(fileName.c_str(),meshNames[0].c_str(),fieldsName[0].c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)its0.size());
  CPPUNIT_ASSERT_EQUAL(-1,its0[0].first);
  CPPUNIT_ASSERT_EQUAL(-1,its0[0].second);
  std::vector<std::pair<int,int> > its1=GetCellFieldIterations(fileName.c_str(),meshNames[0].c_str(),fieldsName[1].c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)its1.size());
  CPPUNIT_ASSERT_EQUAL(-1,its1[0].first);
  CPPUNIT_ASSERT_EQUAL(-1,its1[0].second);
  //
  MEDCouplingFieldDouble *field0=ReadFieldCell(fileName.c_str(),meshNames[0].c_str(),0,fieldsName[0].c_str(),its0[0].first,its0[0].second);
  field0->checkConsistencyLight();
  CPPUNIT_ASSERT(field0->getName()==fieldsName[0]);
 CPPUNIT_ASSERT_EQUAL(1,(int)field0->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(16,(int)field0->getNumberOfTuples());
  const double expectedValues[16]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,2.,3.,3.,2.};
  double diffValue[16];
  std::transform(field0->getArray()->getPointer(),field0->getArray()->getPointer()+16,expectedValues,diffValue,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue,diffValue+16),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue,diffValue+16),1e-12);
  const MEDCouplingUMesh *constMesh=dynamic_cast<const MEDCouplingUMesh *>(field0->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,(int)constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+57,0),1e-12);
  field0->decrRef();
  //
  MEDCouplingFieldDouble *field1=ReadFieldCell(fileName.c_str(),meshNames[0].c_str(),0,fieldsName[1].c_str(),its1[0].first,its1[0].second);
  field1->checkConsistencyLight();
  CPPUNIT_ASSERT(field1->getName()==fieldsName[1]);
 CPPUNIT_ASSERT_EQUAL(3,(int)field1->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(16,(int)field1->getNumberOfTuples());
  const double expectedValues2[48]={1.,0.,1.,1.,0.,1.,1.,0.,1.,2.,1.,0.,2.,1.,0.,2.,1.,0.,3.,0.,1.,3.,0.,1.,3.,0.,1.,4.,1.,0.,4.,1.,0.,4.,1.,0.,5.,0.,0.,6.,1.,1.,6.,0.,0.,5.,1.,1.};
  double diffValue2[48];
  std::transform(field1->getArray()->getPointer(),field1->getArray()->getPointer()+48,expectedValues2,diffValue2,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue2,diffValue2+48),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue2,diffValue2+48),1e-12);
  constMesh=dynamic_cast<const MEDCouplingUMesh *>(field1->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,(int)constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+57,0),1e-12);
  field1->decrRef();
  //fields on nodes
  std::vector<std::string> fieldsNameNode=GetNodeFieldNamesOnMesh(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(2,(int)fieldsNameNode.size());
  CPPUNIT_ASSERT(fieldsNameNode[0]=="fieldnodedouble");
  CPPUNIT_ASSERT(fieldsNameNode[1]=="fieldnodeint");
  std::vector<std::pair<int,int> > its0Node=GetNodeFieldIterations(fileName.c_str(),meshNames[0].c_str(),fieldsNameNode[0].c_str());
  CPPUNIT_ASSERT_EQUAL(3,(int)its0Node.size());
  CPPUNIT_ASSERT_EQUAL(-1,its0Node[0].first);
  CPPUNIT_ASSERT_EQUAL(-1,its0Node[0].second);
  CPPUNIT_ASSERT_EQUAL(1,its0Node[1].first);
  CPPUNIT_ASSERT_EQUAL(-1,its0Node[1].second);
  CPPUNIT_ASSERT_EQUAL(2,its0Node[2].first);
  CPPUNIT_ASSERT_EQUAL(-1,its0Node[2].second);
  MEDCouplingFieldDouble *field0Nodes=ReadFieldNode(fileName.c_str(),meshNames[0].c_str(),0,fieldsNameNode[0].c_str(),its0Node[0].first,its0Node[0].second);
  field0Nodes->checkConsistencyLight();
  CPPUNIT_ASSERT(field0Nodes->getName()==fieldsNameNode[0]);
 CPPUNIT_ASSERT_EQUAL(1,(int)field0Nodes->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(19,(int)field0Nodes->getNumberOfTuples());
  const double expectedValues3[19]={1.,1.,1.,2.,2.,2.,3.,3.,3.,4.,4.,4.,5.,5.,5.,6.,6.,6.,7.};
  double diffValue3[19];
  std::transform(field0Nodes->getArray()->getPointer(),field0Nodes->getArray()->getPointer()+19,expectedValues3,diffValue3,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue3,diffValue3+19),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue3,diffValue3+19),1e-12);
  constMesh=dynamic_cast<const MEDCouplingUMesh *>(field0Nodes->getMesh());
  CPPUNIT_ASSERT(constMesh);
  field0Nodes->decrRef();
  //
  field0Nodes=ReadFieldNode(fileName.c_str(),meshNames[0].c_str(),0,fieldsNameNode[0].c_str(),its0Node[2].first,its0Node[2].second);
  field0Nodes->checkConsistencyLight();
  CPPUNIT_ASSERT(field0Nodes->getName()==fieldsNameNode[0]);
 CPPUNIT_ASSERT_EQUAL(1,(int)field0Nodes->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(19,(int)field0Nodes->getNumberOfTuples());
  const double expectedValues4[19]={1.,2.,2.,2.,3.,3.,3.,4.,4.,4.,5.,5.,5.,6.,6.,6.,7.,7.,7.};
  std::transform(field0Nodes->getArray()->getPointer(),field0Nodes->getArray()->getPointer()+19,expectedValues4,diffValue3,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue3,diffValue3+19),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue3,diffValue3+19),1e-12);
  constMesh=dynamic_cast<const MEDCouplingUMesh *>(field0Nodes->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,(int)constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+57,0),1e-12);
  field0Nodes->decrRef();
  //
  field0Nodes=ReadFieldNode(fileName.c_str(),meshNames[0].c_str(),0,fieldsNameNode[0].c_str(),its0Node[0].first,its0Node[0].second);
  field0Nodes->checkConsistencyLight();
  CPPUNIT_ASSERT(field0Nodes->getName()==fieldsNameNode[0]);
 CPPUNIT_ASSERT_EQUAL(1,(int)field0Nodes->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(19,(int)field0Nodes->getNumberOfTuples());
  const double expectedValues5[19]={1.,1.,1.,2.,2.,2.,3.,3.,3.,4.,4.,4.,5.,5.,5.,6.,6.,6.,7.};
  std::transform(field0Nodes->getArray()->getPointer(),field0Nodes->getArray()->getPointer()+19,expectedValues5,diffValue3,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue3,diffValue3+19),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue3,diffValue3+19),1e-12);
  constMesh=dynamic_cast<const MEDCouplingUMesh *>(field0Nodes->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(16,(int)constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<12;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_HEXA8,constMesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_PYRA5,constMesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL((std::size_t)90,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(701,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+90,0));
  CPPUNIT_ASSERT_EQUAL(705,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+17,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(46.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+57,0),1e-12);
  field0Nodes->decrRef();
}

void MEDLoaderTest::testMEDLoaderPolygonRead()
{
  using namespace std;
  using namespace INTERP_KERNEL;

  string fileName=INTERP_TEST::getResourceFile("polygones.med", 3);
  vector<string> meshNames=GetMeshNames(fileName.c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)meshNames.size());
  CPPUNIT_ASSERT(meshNames[0]=="Bord");
  MEDCouplingUMesh *mesh=ReadUMeshFromFile(fileName.c_str(),meshNames[0].c_str(),0);
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(538,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(579,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  for(int i=0;i<514;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(i));
  for(int i=514;i<538;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+1737,0),1e-12);
  const double expectedVals1[12]={1.4851585216522212,-0.5,0.,1.4851585216522212,-0.4,0.,1.4851585216522212,-0.3,0., 1.5741585216522211, -0.5, 0. };
  double diffValue1[12];
  std::transform(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+12,expectedVals1,diffValue1,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue1,diffValue1+12),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue1,diffValue1+12),1e-12);
  CPPUNIT_ASSERT_EQUAL((std::size_t)2768,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(651050,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+2768,0));
  CPPUNIT_ASSERT_EQUAL(725943,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+539,0));
  mesh->decrRef();
  //
  std::vector<std::string> fieldsName=GetCellFieldNamesOnMesh(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(3,(int)fieldsName.size());
  CPPUNIT_ASSERT(fieldsName[0]=="bord_:_distorsion");
  CPPUNIT_ASSERT(fieldsName[1]=="bord_:_familles");
  CPPUNIT_ASSERT(fieldsName[2]=="bord_:_non-ortho");
  std::vector<std::pair<int,int> > its0=GetCellFieldIterations(fileName.c_str(),meshNames[0].c_str(),fieldsName[0].c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)its0.size());
  MEDCouplingFieldDouble *field=ReadFieldCell(fileName.c_str(),meshNames[0].c_str(),0,fieldsName[0].c_str(),its0[0].first,its0[0].second);
  field->checkConsistencyLight();
  CPPUNIT_ASSERT(field->getName()==fieldsName[0]);
 CPPUNIT_ASSERT_EQUAL(1,(int)field->getNumberOfComponents());
 CPPUNIT_ASSERT_EQUAL(538,(int)field->getNumberOfTuples());
  const MEDCouplingUMesh *constMesh=dynamic_cast<const MEDCouplingUMesh *>(field->getMesh());
  CPPUNIT_ASSERT(constMesh);
  CPPUNIT_ASSERT_EQUAL(3,constMesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,constMesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(538,(int)constMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(579,constMesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)constMesh->getAllGeoTypes().size());
  for(int i=0;i<514;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,constMesh->getTypeOfCell(i));
  for(int i=514;i<538;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,constMesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,std::accumulate(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+1737,0),1e-12);
  std::transform(constMesh->getCoords()->getConstPointer(),constMesh->getCoords()->getConstPointer()+12,expectedVals1,diffValue1,std::minus<double>());
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::max_element(diffValue1,diffValue1+12),1e-12);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(0.,*std::min_element(diffValue1,diffValue1+12),1e-12);
  CPPUNIT_ASSERT_EQUAL((std::size_t)2768,constMesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(651050,std::accumulate(constMesh->getNodalConnectivity()->getConstPointer(),constMesh->getNodalConnectivity()->getConstPointer()+2768,0));
  CPPUNIT_ASSERT_EQUAL(725943,std::accumulate(constMesh->getNodalConnectivityIndex()->getConstPointer(),constMesh->getNodalConnectivityIndex()->getConstPointer()+539,0));
  const double *values=field->getArray()->getPointer();
  CPPUNIT_ASSERT_DOUBLES_EQUAL(2.87214203182918,std::accumulate(values,values+538,0.),1e-12);
  field->decrRef();
}

void MEDLoaderTest::testMEDLoaderPolyhedronRead()
{
  using namespace std;
  using namespace INTERP_KERNEL;

  string fileName=INTERP_TEST::getResourceFile("poly3D.med", 3);
  vector<string> meshNames=GetMeshNames(fileName.c_str());
  CPPUNIT_ASSERT_EQUAL(1,(int)meshNames.size());
  CPPUNIT_ASSERT(meshNames[0]=="poly3D");
  MEDCouplingUMesh *mesh=ReadUMeshFromFile(fileName.c_str(),meshNames[0].c_str(),0);
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(2,(int)mesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(NORM_TETRA4,mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(NORM_POLYHED,mesh->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(NORM_POLYHED,mesh->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL((std::size_t)98,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(725,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+98,0));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(110.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  CPPUNIT_ASSERT_EQUAL(155,std::accumulate(mesh->getNodalConnectivityIndex()->getPointer(),mesh->getNodalConnectivityIndex()->getPointer()+4,0));
  mesh->decrRef();
  //
  mesh=ReadUMeshFromFile(fileName.c_str(),meshNames[0].c_str(),-1);
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(17,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,(int)mesh->getAllGeoTypes().size());
  CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(0));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(1));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(2));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(3));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(4));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(5));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(6));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(7));
  CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(8));
  CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(9));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(10));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(11));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(12));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(13));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(14));
  CPPUNIT_ASSERT_EQUAL(NORM_QUAD4,mesh->getTypeOfCell(15));
  CPPUNIT_ASSERT_EQUAL(NORM_TRI3,mesh->getTypeOfCell(16));
  CPPUNIT_ASSERT_DOUBLES_EQUAL(110.,std::accumulate(mesh->getCoords()->getPointer(),mesh->getCoords()->getPointer()+57,0),1e-12);
  CPPUNIT_ASSERT_EQUAL((std::size_t)83,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(619,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+83,0));
  mesh->decrRef();
  //
  vector<string> families=GetMeshFamiliesNames(fileName.c_str(),meshNames[0].c_str());
  CPPUNIT_ASSERT_EQUAL(4,(int)families.size());
  CPPUNIT_ASSERT(families[0]=="FAMILLE_FACE_POLYGONS3");
  CPPUNIT_ASSERT(families[1]=="FAMILLE_FACE_QUAD41");
  CPPUNIT_ASSERT(families[2]=="FAMILLE_FACE_TRIA32");
  CPPUNIT_ASSERT(families[3]=="FAMILLE_ZERO");
  vector<string> families2;
  families2.push_back(families[0]);
  mesh=ReadUMeshFromFamilies(fileName.c_str(),meshNames[0].c_str(),-1,families2);
  mesh->checkConsistencyLight();
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(2,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(3,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(1,(int)mesh->getAllGeoTypes().size());
  for(int i=0;i<3;i++)
    CPPUNIT_ASSERT_EQUAL(NORM_POLYGON,mesh->getTypeOfCell(i));
  CPPUNIT_ASSERT_EQUAL((std::size_t)19,mesh->getNodalConnectivity()->getNbOfElems());
  CPPUNIT_ASSERT_EQUAL(117,std::accumulate(mesh->getNodalConnectivity()->getPointer(),mesh->getNodalConnectivity()->getPointer()+19,0));
  mesh->decrRef();
  //
  mesh=ReadUMeshFromFamilies(fileName.c_str(),meshNames[0].c_str(),0,families2);
  CPPUNIT_ASSERT_EQUAL(3,mesh->getSpaceDimension());
  CPPUNIT_ASSERT_EQUAL(0,(int)mesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(19,mesh->getNumberOfNodes());
  CPPUNIT_ASSERT_EQUAL(3,mesh->getMeshDimension());
  CPPUNIT_ASSERT_EQUAL(0,(int)mesh->getAllGeoTypes().size());
  mesh->decrRef();
}

MEDCouplingUMesh *MEDLoaderTest::build1DMesh_1()
{
  double coords[6]={ 0.0, 0.3, 0.75, 1.0, 1.4, 1.3 };
  int conn[9]={ 0,1, 1,2, 2,3 , 3,4,5};
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setName("1DMesh_1");
  mesh->setMeshDimension(1);
  mesh->allocateCells(4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG3,3,conn+6);
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(6,1);
  myCoords->setInfoOnComponent(0,"tototototototot [m*m*m*m*m*m*m*m]");
  std::copy(coords,coords+6,myCoords->getPointer());
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  return mesh;
}

MEDCouplingUMesh *MEDLoaderTest::build2DCurveMesh_1()
{
  double coords[12]={ 0.0,0.0, 0.3,0.3, 0.75,0.75, 1.0,1.0, 1.4,1.4, 1.3,1.3 };
  int conn[9]={ 0,1, 1,2, 2,3 , 3,4,5};
  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setName("2DCurveMesh_1");
  mesh->setMeshDimension(1);
  mesh->allocateCells(4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+4);
  mesh->insertNextCell(INTERP_KERNEL::NORM_SEG3,3,conn+6);
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(6,2);
  std::copy(coords,coords+12,myCoords->getPointer());
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  return mesh;
}

MEDCouplingUMesh *MEDLoaderTest::build2DMesh_1()
{
  double targetCoords[24]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7, -0.05,0.95, 0.2,1.2, 0.45,0.95 };
  int targetConn[24]={1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(6);
  targetMesh->setName("2DMesh_1");
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,targetConn+6);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+16);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POLYGON,4,targetConn+20);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(12,2);
  myCoords->setInfoOnComponent(0,"tototototototot [m]");
  myCoords->setInfoOnComponent(1,"energie [kW]");
  std::copy(targetCoords,targetCoords+24,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDLoaderTest::build2DMesh_2()
{
  double targetCoords[24]={
    -0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7,
    -0.05,0.95, 0.2,1.2, 0.45,0.95
  };
  int targetConn[24]={1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(5);
  targetMesh->setName("2DMesh_2");
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+16);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,targetConn+6);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(12,2);
  myCoords->setInfoOnComponent(0,"toto [m]");
  myCoords->setInfoOnComponent(1,"energie [kW]");
  std::copy(targetCoords,targetCoords+24,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDLoaderTest::build3DSurfMesh_1()
{
  double targetCoords[36]={
    -0.3,-0.3,-0.3, 0.2,-0.3,-0.3, 0.7,-0.3,-0.3, -0.3,0.2,-0.3, 0.2,0.2,-0.3, 0.7,0.2,-0.3, -0.3,0.7,-0.3, 0.2,0.7,-0.3, 0.7,0.7,-0.3
    ,-0.05,0.95,-0.3, 0.2,1.2,-0.3, 0.45,0.95,-0.3
  };
  int targetConn[24]={1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(6);
  targetMesh->setName("3DSurfMesh_1");
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+16);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,targetConn+6);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POLYGON,4,targetConn+20);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(12,3);
  myCoords->setInfoOnComponent(0,"toto [m]");
  myCoords->setInfoOnComponent(2,"ff [km]");//component 1 is not set for test
  std::copy(targetCoords,targetCoords+36,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDLoaderTest::build3DMesh_1()
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
  //
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  ret->setName("3DMesh_1");
  ret->setMeshDimension(3);
  ret->allocateCells(18);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+51);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+59);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+110);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+118);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+169);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+177);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+228);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+236);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+287);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+295);
  ret->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,conn+346);
  //
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+8);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+67);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+126);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+185);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+244);
  ret->insertNextCell(INTERP_KERNEL::NORM_POLYHED,43,conn+303);
  //
  ret->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(60,3);
  myCoords->setInfoOnComponent(0,"titi [m]");
  myCoords->setInfoOnComponent(1,"density power [MW/m^3]");
  myCoords->setInfoOnComponent(2,"t [kW]");
  std::copy(coords,coords+180,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
  return ret;
}

MEDCouplingUMesh *MEDLoaderTest::build3DMesh_2()
{
  MEDCouplingUMesh *m3dsurfBase=build3DSurfMesh_1();
  int numbers[5]={0,1,2,3,5};
  MEDCouplingUMesh *m3dsurf=(MEDCouplingUMesh *)m3dsurfBase->buildPartOfMySelf(numbers,numbers+5,false);
  m3dsurfBase->decrRef();
  MEDCouplingUMesh *m1dBase=build1DMesh_1();
  int numbers2[4]={0,1,2,3};
  MEDCouplingUMesh *m1d=(MEDCouplingUMesh *)m1dBase->buildPartOfMySelf(numbers2,numbers2+4,false);
  m1dBase->decrRef();
  m1d->changeSpaceDimension(3);
  const double vec[3]={0.,1.,0.};
  const double pt[3]={0.,0.,0.};
  m1d->rotate(pt,vec,-M_PI/2.);
  MEDCouplingUMesh *ret=m3dsurf->buildExtrudedMesh(m1d,0);
  m1d->decrRef();
  m3dsurf->decrRef();
  return ret;
}

MEDCouplingFieldDouble *MEDLoaderTest::buildVecFieldOnCells_1()
{
  MEDCouplingUMesh *mesh=build3DSurfMesh_1();
  int nbOfCells=mesh->getNumberOfCells();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setName("VectorFieldOnCells");
  f1->setMesh(mesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(nbOfCells,3);
  array->setInfoOnComponent(0,"power [MW/m^3]");
  array->setInfoOnComponent(1,"density [g/cm^3]");
  array->setInfoOnComponent(2,"temperature [K]");
  f1->setArray(array);
  array->decrRef();
  double *tmp=array->getPointer();
  const double arr1[18]={0.,10.,20.,1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.};
  std::copy(arr1,arr1+18,tmp);
  f1->setTime(2.,0,1);
  f1->checkConsistencyLight();
  mesh->decrRef();
  return f1;
}

MEDCouplingFieldDouble *MEDLoaderTest::buildVecFieldOnNodes_1()
{
  MEDCouplingUMesh *mesh=build3DSurfMesh_1();
  int nbOfNodes=mesh->getNumberOfNodes();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  f1->setName("VectorFieldOnNodes");
  f1->setMesh(mesh);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(nbOfNodes,3);
  f1->setArray(array);
  array->setInfoOnComponent(0,"power [MW/m^3]");
  array->setInfoOnComponent(1,"density [g/cm^3]");
  array->setInfoOnComponent(2,"temperature [K]");
  array->decrRef();
  double *tmp=array->getPointer();
  const double arr1[36]={
    70.,80.,90.,71.,81.,91.,72.,82.,92.,73.,83.,93.,74.,84.,94.,75.,85.,95.,
    1000.,10010.,10020.,1001.,10011.,10021.,1002.,10012.,10022.,1003.,10013.,10023.,1004.,10014.,10024.,1005.,10015.,10025.,
  };
  std::copy(arr1,arr1+36,tmp);
  f1->setTime(2.12,2,3);
  f1->checkConsistencyLight();
  mesh->decrRef();
  return f1;
}

MEDCouplingFieldDouble *MEDLoaderTest::buildVecFieldOnGauss_1()
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
  MEDCouplingUMesh *m=build2DMesh_2();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_GAUSS_PT,ONE_TIME);
  f->setTime(3.14,1,5);
  f->setMesh(m);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_TRI3,_refCoo1,_gsCoo1,_wg1);
  const double refCoo2[12]={-1.0,1.0, -1.0,-1.0, 1.0,-1.0, -1.0,0.0, 0.0,-1.0, 0.0,0.0 };
  std::vector<double> _refCoo2(refCoo2,refCoo2+12);
  std::vector<double> _gsCoo2(_gsCoo1);
  std::vector<double> _wg2(_wg1);
  _gsCoo2.resize(6); _wg2.resize(3);
  const double refCoo3[8]={ 0.,0., 1.,0., 1.,1., 0.,1. };
  std::vector<double> _refCoo3(refCoo3,refCoo3+8);
  _gsCoo1.resize(4); _wg1.resize(2);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_QUAD4,_refCoo3,_gsCoo1,_wg1);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_TRI6,_refCoo2,_gsCoo2,_wg2);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(19,2);
  double *ptr=array->getPointer();
  for(int i=0;i<19*2;i++)
    ptr[i]=(double)(i+7);
  f->setArray(array);
  f->setName("MyFirstFieldOnGaussPoint");
  array->setInfoOnComponent(0,"power [MW/m^3]");
  array->setInfoOnComponent(1,"density");
  array->decrRef();
  f->checkConsistencyLight();
  m->decrRef();
  return f;
}

MEDCouplingFieldDouble *MEDLoaderTest::buildVecFieldOnGaussNE_1()
{
  MEDCouplingUMesh *m=build2DMesh_2();
  MEDCouplingFieldDouble *f=MEDCouplingFieldDouble::New(ON_GAUSS_NE,ONE_TIME);
  f->setTime(3.14,1,5);
  f->setMesh(m);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(20,2);
  double *ptr=array->getPointer();
  for(int i=0;i<20*2;i++)
    ptr[i]=(double)(i+8);
  f->setArray(array);
  array->setInfoOnComponent(0,"power [W]");
  array->setInfoOnComponent(1,"temperature");
  f->setName("MyFieldOnGaussNE");
  array->decrRef();
  f->checkConsistencyLight();
  m->decrRef();
  return f;
}



