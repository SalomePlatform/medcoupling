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

#include "MEDLoaderTest.hxx"
#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"

using namespace ParaMEDMEM;

void MEDLoaderTest::testMesh1DRW()
{
  MEDCouplingUMesh *mesh=build1DMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file1.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file1.med",mesh->getName(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh2DCurveRW()
{
  MEDCouplingUMesh *mesh=build2DCurveMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file2.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file2.med",mesh->getName(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh2DRW()
{
  MEDCouplingUMesh *mesh=build2DMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file3.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file3.med",mesh->getName(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh3DSurfRW()
{
  MEDCouplingUMesh *mesh=build3DSurfMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file4.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file4.med",mesh->getName(),0);
  CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
  mesh_rw->decrRef();
  mesh->decrRef();
}

void MEDLoaderTest::testMesh3DRW()
{
  MEDCouplingUMesh *mesh=build3DMesh_1();
  mesh->checkCoherency();
  MEDLoader::WriteUMesh("file5.med",mesh,true);
  MEDCouplingUMesh *mesh_rw=MEDLoader::ReadUMeshFromFile("file5.med",mesh->getName(),0);
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
  MEDLoader::WriteField("file6.med",f1,true);
  MEDCouplingFieldDouble *f2=MEDLoader::ReadFieldDoubleCell("file6.med",f1->getMesh()->getName(),0,f1->getName(),0,1);
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f1->decrRef();
  f2->decrRef();
  //
  f1=buildVecFieldOnNodes_1();
  MEDLoader::WriteField("file7.med",f1,true);
  f2=MEDLoader::ReadFieldDoubleNode("file7.med",f1->getMesh()->getName(),0,f1->getName(),2,3);
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
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
  MEDLoader::WriteField(fileName,f1,true);
  f1->setTime(10.,8,9);
  double *tmp=f1->getArray()->getPointer();
  tmp[0]=VAL1;
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(10.14,18,19);
  tmp[0]=VAL2;
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  //retrieving time steps...
  MEDCouplingFieldDouble *f2=MEDLoader::ReadFieldDoubleCell(fileName,f1->getMesh()->getName(),0,f1->getName(),8,9);
  f1->setTime(10.,8,9);
  tmp[0]=VAL1;
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f2->decrRef();
  f2=MEDLoader::ReadFieldDoubleCell(fileName,f1->getMesh()->getName(),0,f1->getName(),0,1);
  MEDCouplingFieldDouble *f3=buildVecFieldOnCells_1();
  CPPUNIT_ASSERT(f3->isEqual(f2,1e-12,1e-12));
  f3->decrRef();
  f2->decrRef();
  f2=MEDLoader::ReadFieldDoubleCell(fileName,f1->getMesh()->getName(),0,f1->getName(),18,19);
  f1->setTime(10.14,18,19);
  tmp[0]=VAL2;
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f2->decrRef();
  f1->decrRef();
  //ON NODES
  f1=buildVecFieldOnNodes_1();
  const char fileName2[]="file9.med";
  MEDLoader::WriteField(fileName2,f1,true);
  f1->setTime(110.,108,109);
  tmp=f1->getArray()->getPointer();
  tmp[3]=VAL1;
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName2,f1);
  f1->setTime(210.,208,209);
  tmp[3]=VAL2;
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName2,f1);
  f2=MEDLoader::ReadFieldDoubleNode(fileName2,f1->getMesh()->getName(),0,f1->getName(),108,109);
  f1->setTime(110.,108,109);
  tmp[3]=VAL1;
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f2->decrRef();
  f2=MEDLoader::ReadFieldDoubleNode(fileName2,f1->getMesh()->getName(),0,f1->getName(),2,3);
  f3=buildVecFieldOnNodes_1();
  CPPUNIT_ASSERT(f3->isEqual(f2,1e-12,1e-12));
  f3->decrRef();
  f2->decrRef();
  f2=MEDLoader::ReadFieldDoubleNode(fileName2,f1->getMesh()->getName(),0,f1->getName(),208,209);
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
  const char name2[]="AMesh2";
  MEDCouplingFieldDouble *f1=buildVecFieldOnCells_1();
  ((MEDCouplingMesh *)f1->getMesh())->setName(name3);
  f1->setName(name1);
  f1->setTime(10.,8,9);
  double *tmp=f1->getArray()->getPointer();
  tmp[0]=VAL1;
  MEDLoader::WriteField(fileName,f1,true);
  f1->setTime(10.14,18,19);
  tmp[0]=VAL2;
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  ((MEDCouplingMesh *)f1->getMesh())->setName(name2);
  f1->setTime(10.55,28,29);
  tmp[0]=3*VAL1;
  MEDLoader::WriteField(fileName,f1,false);
  f1->setTime(10.66,38,39);
  tmp[0]=3*VAL2;
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(10.77,48,49);
  tmp[0]=4*VAL2;
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  //ON NODES
  f1->decrRef();
  f1=buildVecFieldOnNodes_1();
  f1->setName(name1);
  ((MEDCouplingMesh *)f1->getMesh())->setName(name2);
  f1->setTime(110.,8,9);
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(110.,108,109);
  tmp=f1->getArray()->getPointer();
  tmp[3]=VAL1;
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  f1->setTime(210.,208,209);
  tmp[3]=VAL2;
  MEDLoader::WriteFieldUsingAlreadyWrittenMesh(fileName,f1);
  //
  std::vector< std::pair<int,int> > it1=MEDLoader::GetCellFieldIterations(fileName,name3,name1);
  CPPUNIT_ASSERT_EQUAL(2,(int)it1.size());
  CPPUNIT_ASSERT_EQUAL(8,it1[0].first); CPPUNIT_ASSERT_EQUAL(9,it1[0].second);
  CPPUNIT_ASSERT_EQUAL(18,it1[1].first); CPPUNIT_ASSERT_EQUAL(19,it1[1].second);
  std::vector< std::pair<int,int> > it2=MEDLoader::GetCellFieldIterations(fileName,name2,name1);
  CPPUNIT_ASSERT_EQUAL(3,(int)it2.size());
  CPPUNIT_ASSERT_EQUAL(28,it2[0].first); CPPUNIT_ASSERT_EQUAL(29,it2[0].second);
  CPPUNIT_ASSERT_EQUAL(38,it2[1].first); CPPUNIT_ASSERT_EQUAL(39,it2[1].second);
  CPPUNIT_ASSERT_EQUAL(48,it2[2].first); CPPUNIT_ASSERT_EQUAL(49,it2[2].second);
  std::vector< std::pair<int,int> > it3=MEDLoader::GetNodeFieldIterations(fileName,name2,name1);
  CPPUNIT_ASSERT_EQUAL(3,(int)it3.size());
  CPPUNIT_ASSERT_EQUAL(8,it3[0].first); CPPUNIT_ASSERT_EQUAL(9,it3[0].second);
  CPPUNIT_ASSERT_EQUAL(108,it3[1].first); CPPUNIT_ASSERT_EQUAL(109,it3[1].second);
  CPPUNIT_ASSERT_EQUAL(208,it3[2].first); CPPUNIT_ASSERT_EQUAL(209,it3[2].second);
  std::vector< std::pair<int,int> > it4=MEDLoader::GetNodeFieldIterations(fileName,name3,name1);
  CPPUNIT_ASSERT(it4.empty());
  //
  f1->decrRef();
  //
  f1=MEDLoader::ReadFieldDoubleCell(fileName,name3,0,name1,8,9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(VAL1,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  f1=MEDLoader::ReadFieldDoubleCell(fileName,name3,0,name1,18,19);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(VAL2,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  f1=MEDLoader::ReadFieldDoubleCell(fileName,name2,0,name1,28,29);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3*VAL1,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  f1=MEDLoader::ReadFieldDoubleCell(fileName,name2,0,name1,38,39);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3*VAL2,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  f1=MEDLoader::ReadFieldDoubleCell(fileName,name2,0,name1,48,49);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(4*VAL2,f1->getArray()->getConstPointer()[0],1e-13);
  f1->decrRef();
  //
  f1=MEDLoader::ReadFieldDoubleNode(fileName,name2,0,name1,8,9);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(71.,f1->getArray()->getConstPointer()[3],1e-13);
  f1->decrRef();
  f1=MEDLoader::ReadFieldDoubleNode(fileName,name2,0,name1,108,109);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(VAL1,f1->getArray()->getConstPointer()[3],1e-13);
  f1->decrRef();
  f1=MEDLoader::ReadFieldDoubleNode(fileName,name2,0,name1,208,209);
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
  std::vector<MEDCouplingUMesh *> meshes;
  meshes.push_back(mesh1);
  meshes.push_back(mesh2);
  meshes.push_back(mesh3);
  meshes.push_back(mesh4);
  const char mnane[]="3DToto";
  MEDLoader::WriteUMeshes(fileName,mnane,meshes,true);
  //
  MEDCouplingUMesh *mesh5=MEDLoader::ReadUMeshFromFile(fileName,mnane);
  mesh1->setName(mnane);
  const int part3[18]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
  MEDCouplingUMesh *mesh6=(MEDCouplingUMesh *)mesh5->buildPartOfMySelf(part3,part3+18,true);
  mesh6->setName(mnane);
  mesh5->decrRef();
  CPPUNIT_ASSERT(mesh6->isEqual(mesh1,1e-12));
  mesh6->decrRef();
  std::vector<std::string> grps=MEDLoader::GetMeshGroupsNames(fileName,mnane);
  CPPUNIT_ASSERT_EQUAL(4,(int)grps.size());
  CPPUNIT_ASSERT(std::find(grps.begin(),grps.end(),std::string("mesh2"))!=grps.end());
  CPPUNIT_ASSERT(std::find(grps.begin(),grps.end(),std::string("mesh3"))!=grps.end());
  CPPUNIT_ASSERT(std::find(grps.begin(),grps.end(),std::string("mesh4"))!=grps.end());
  CPPUNIT_ASSERT(std::find(grps.begin(),grps.end(),std::string("3DMesh_1"))!=grps.end());
  //
  std::vector<std::string> vec;
  vec.push_back(std::string("mesh2"));
  MEDCouplingUMesh *mesh2_2=MEDLoader::ReadUMeshFromGroups(fileName,mnane,0,vec);
  CPPUNIT_ASSERT(mesh2_2->isEqual(mesh2,1e-12));
  mesh2_2->decrRef();
  vec.clear(); vec.push_back(std::string("mesh3"));
  MEDCouplingUMesh *mesh3_2=MEDLoader::ReadUMeshFromGroups(fileName,mnane,0,vec);
  CPPUNIT_ASSERT(mesh3_2->isEqual(mesh3,1e-12));
  mesh3_2->decrRef();
  vec.clear(); vec.push_back(std::string("mesh4"));
  MEDCouplingUMesh *mesh4_2=MEDLoader::ReadUMeshFromGroups(fileName,mnane,0,vec);
  CPPUNIT_ASSERT(mesh4_2->isEqual(mesh4,1e-12));
  mesh4_2->decrRef();
  vec.clear(); vec.push_back(std::string("3DMesh_1"));
  MEDCouplingUMesh *mesh1_2=MEDLoader::ReadUMeshFromGroups(fileName,mnane,0,vec);
  mesh1->setName("3DMesh_1");
  CPPUNIT_ASSERT(mesh1_2->isEqual(mesh1,1e-12));
  mesh1_2->decrRef();
  //
  vec.clear(); vec.push_back(std::string("Family_4")); vec.push_back(std::string("Family_2"));
  mesh2_2=MEDLoader::ReadUMeshFromFamilies(fileName,mnane,0,vec);
  mesh2_2->setName("mesh2");
  CPPUNIT_ASSERT(mesh2_2->isEqual(mesh2,1e-12));
  mesh2_2->decrRef();
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
  DataArrayInt *da=mesh1->mergeNodes(1e-12,b);
  da->decrRef();
  MEDLoader::WriteUMesh(fileName,mesh1,true);
  const int part1[5]={1,2,4,13,15};
  MEDCouplingUMesh *mesh2=(MEDCouplingUMesh *)mesh1->buildPartOfMySelf(part1,part1+5,true);
  mesh2->setName(mesh1->getName());//<- important for the test
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
  f1->checkCoherency();
  //
  MEDLoader::WriteField(fileName,f1,false);//<- false important for the test
  //
  MEDCouplingFieldDouble *f2=MEDLoader::ReadFieldDoubleCell(fileName,f1->getMesh()->getName(),0,f1->getName(),2,7);
  f2->checkCoherency();
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  //
  f2->decrRef();
  f1->decrRef();
  mesh1->decrRef();
  mesh2->decrRef();
}

void MEDLoaderTest::testFieldGaussRW1()
{
  const char fileName[]="file13.med";
  MEDCouplingFieldDouble *f1=buildVecFieldOnGauss_1();
  MEDLoader::WriteField(fileName,f1,true);
  MEDCouplingFieldDouble *f2=MEDLoader::ReadFieldDouble(ON_GAUSS_PT,fileName,f1->getMesh()->getName(),0,f1->getName(),1,5);
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f2->decrRef();
  f1->decrRef();
}

void MEDLoaderTest::testFieldGaussNERW1()
{
  const char fileName[]="file14.med";
  MEDCouplingFieldDouble *f1=buildVecFieldOnGaussNE_1();
  MEDLoader::WriteField(fileName,f1,true);
  MEDCouplingFieldDouble *f2=MEDLoader::ReadFieldDouble(ON_GAUSS_NE,fileName,f1->getMesh()->getName(),0,f1->getName(),1,5);
  CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
  f2->decrRef();
  f1->decrRef();
}

void MEDLoaderTest::testLittleStrings1()
{
  std::string s("azeeeerrrtty");
  MEDLoaderBase::zipEqualConsChar(s,3);
  CPPUNIT_ASSERT(s=="azertty");
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
  myCoords->setInfoOnComponent(0,"tototototototot (m*m*m*m*m*m*m*m)");
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
  double targetCoords[24]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
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
  myCoords->setInfoOnComponent(0,"tototototototot (m)");
  myCoords->setInfoOnComponent(1,"energie (kW)");
  std::copy(targetCoords,targetCoords+24,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDLoaderTest::build2DMesh_2()
{
  double targetCoords[24]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
  int targetConn[24]={1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(5);
  targetMesh->setName("2DMesh_2");
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,targetConn+6);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+16);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(12,2);
  myCoords->setInfoOnComponent(0,"toto (m)");
  myCoords->setInfoOnComponent(1,"energie (kW)");
  std::copy(targetCoords,targetCoords+24,myCoords->getPointer());
  targetMesh->setCoords(myCoords);
  myCoords->decrRef();
  return targetMesh;
}

MEDCouplingUMesh *MEDLoaderTest::build3DSurfMesh_1()
{
  double targetCoords[36]={-0.3,-0.3,-0.3, 0.2,-0.3,-0.3, 0.7,-0.3,-0.3, -0.3,0.2,-0.3, 0.2,0.2,-0.3, 0.7,0.2,-0.3, -0.3,0.7,-0.3, 0.2,0.7,-0.3, 0.7,0.7,-0.3 };
  int targetConn[24]={1,4,2, 4,5,2, 6,10,8,9,11,7, 0,3,4,1, 6,7,4,3, 7,8,5,4};
  MEDCouplingUMesh *targetMesh=MEDCouplingUMesh::New();
  targetMesh->setMeshDimension(2);
  targetMesh->allocateCells(6);
  targetMesh->setName("3DSurfMesh_1");
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_TRI6,6,targetConn+6);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+12);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+16);
  targetMesh->insertNextCell(INTERP_KERNEL::NORM_POLYGON,4,targetConn+20);
  targetMesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(12,3);
  myCoords->setInfoOnComponent(0,"toto (m)");
  myCoords->setInfoOnComponent(2,"ff (km)");//component 1 is not set for test
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
  myCoords->setInfoOnComponent(0,"titi (m)");
  myCoords->setInfoOnComponent(1,"density power (MW/m^3)");
  myCoords->setInfoOnComponent(2,"t (kW)");
  std::copy(coords,coords+180,myCoords->getPointer());
  ret->setCoords(myCoords);
  myCoords->decrRef();
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
  array->setInfoOnComponent(0,"power (MW/m^3)");
  array->setInfoOnComponent(1,"density (g/cm^3)");
  array->setInfoOnComponent(2,"temperature (K)");
  f1->setArray(array);
  array->decrRef();
  double *tmp=array->getPointer();
  const double arr1[18]={0.,10.,20.,1.,11.,21.,2.,12.,22.,3.,13.,23.,4.,14.,24.,5.,15.,25.};
  std::copy(arr1,arr1+18,tmp);
  f1->setTime(2.,0,1);
  f1->checkCoherency();
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
  array->setInfoOnComponent(0,"power (MW/m^3)");
  array->setInfoOnComponent(1,"density (g/cm^3)");
  array->setInfoOnComponent(2,"temperature (K)");
  array->decrRef();
  double *tmp=array->getPointer();
  const double arr1[36]={
    70.,80.,90.,71.,81.,91.,72.,82.,92.,73.,83.,93.,74.,84.,94.,75.,85.,95.,
    1000.,10010.,10020.,1001.,10011.,10021.,1002.,10012.,10022.,1003.,10013.,10023.,1004.,10014.,10024.,1005.,10015.,10025.,
  };
  std::copy(arr1,arr1+36,tmp);
  f1->setTime(2.12,2,3);
  f1->checkCoherency();
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
  _gsCoo1.resize(6); _wg1.resize(3);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_TRI6,_refCoo2,_gsCoo1,_wg1);
  const double refCoo3[8]={ 0.,0., 1.,0., 1.,1., 0.,1. };
  std::vector<double> _refCoo3(refCoo3,refCoo3+8);
  _gsCoo1.resize(4); _wg1.resize(2);
  f->setGaussLocalizationOnType(INTERP_KERNEL::NORM_QUAD4,_refCoo3,_gsCoo1,_wg1);
  DataArrayDouble *array=DataArrayDouble::New();
  array->alloc(19,2);
  double *ptr=array->getPointer();
  for(int i=0;i<19*2;i++)
    ptr[i]=(double)(i+7);
  f->setArray(array);
  f->setName("MyFirstFieldOnGaussPoint");
  array->setInfoOnComponent(0,"power (MW/m^3)");
  array->setInfoOnComponent(1,"density");
  array->decrRef();
  f->checkCoherency();
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
  array->setInfoOnComponent(0,"power (W)");
  array->setInfoOnComponent(1,"temperature");
  f->setName("MyFieldOnGaussNE");
  array->decrRef();
  f->checkCoherency();
  m->decrRef();
  return f;
}
