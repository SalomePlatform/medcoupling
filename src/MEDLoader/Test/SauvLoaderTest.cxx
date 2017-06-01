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

#include "SauvLoaderTest.hxx"

#include "SauvReader.hxx"
#include "SauvWriter.hxx"
#include "MEDFileData.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "TestInterpKernelUtils.hxx"  // getResourceFile()

#ifdef WIN32
#include <windows.h>
#else
# include <unistd.h>
#endif

#include <vector>
#include <string>

using namespace MEDCoupling;

void SauvLoaderTest::testSauv2Med()
{
  // read a file containing all types of readable piles
  std::string file = INTERP_TEST::getResourceFile("allPillesTest.sauv", 3);
  MCAuto<SauvReader> sr=SauvReader::New(file.c_str());
  MCAuto<MEDFileData> d2=sr->loadInMEDFileDS();
  // write MED
  d2->write("allPillesTest.med",0);
  // check
  CPPUNIT_ASSERT_EQUAL(1,d2->getNumberOfMeshes());
  CPPUNIT_ASSERT_EQUAL(8+97,d2->getNumberOfFields());
  MEDFileMesh * m = d2->getMeshes()->getMeshAtPos(0);
  CPPUNIT_ASSERT_EQUAL(17,int(m->getGroupsNames().size()));
}

void SauvLoaderTest::testMed2SauvOnAMeshWithVoidFamily()
{
  // Create a mesh with 2 quads.
  const int spaceDim = 2;
  const int nbOfNodes = 6;
  double coords[nbOfNodes*spaceDim] = {0,0, 1,0, 1,1, 0,1, 2,0, 2,1};
  int conn[8]={0,1,2,3, 1,4,5,2};
  MCAuto<MEDCouplingUMesh> mesh2d=MEDCouplingUMesh::New("Mesh",spaceDim);
  mesh2d->allocateCells(2);
  mesh2d->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
  mesh2d->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+4);
  mesh2d->finishInsertingCells();
  MCAuto<DataArrayDouble> myCoords=DataArrayDouble::New();
  myCoords->alloc(nbOfNodes,spaceDim);
  std::copy(coords,coords+nbOfNodes*spaceDim,myCoords->getPointer());
  mesh2d->setCoords(myCoords);

  // create a MedFileUMesh
  MCAuto<MEDFileUMesh> m= MEDFileUMesh::New();
  m->setMeshAtLevel(0,mesh2d);

  // Create families and groups

  MCAuto<DataArrayInt> fam = DataArrayInt::New();
  fam->alloc(2,1);
  int elemsFams[2] = {-2,-3};
  std::copy(elemsFams,elemsFams+2,fam->getPointer());

  m->setFamilyFieldArr(0,fam);

  std::map<std::string,int> theFamilies;
  theFamilies["FAM_-1"]=-1;
  theFamilies["FAM_-2"]=-2;
  theFamilies["FAM_-3"]=-3;

  std::map<std::string, std::vector<std::string> > theGroups;
  theGroups["Group1"].push_back("FAM_-2");
  theGroups["Group2"].push_back("FAM_-3");
  theGroups["Grouptot"].push_back("FAM_-1");
  theGroups["Grouptot"].push_back("FAM_-2");
  theGroups["Grouptot"].push_back("FAM_-3");
  m->setFamilyInfo(theFamilies);
  m->setGroupInfo(theGroups);

  // write to MED for visual check
  //const char* medFile = "mesh_with_void_family.med";
  //m->write(medFile, 2);

  // write to SAUV
  const char* sauvFile = "mesh_with_void_family.sauv";
  MCAuto<MEDFileData> medData = MEDFileData::New();
  MCAuto<MEDFileMeshes> medMeshes = MEDFileMeshes::New();
  MCAuto<SauvWriter> sw=SauvWriter::New();
  medMeshes->setMeshAtPos(0, m);
  medData->setMeshes(medMeshes);
  sw->setMEDFileDS(medData);
  sw->write(sauvFile);

  // read SAUV and check groups
  MCAuto<SauvReader> sr=SauvReader::New(sauvFile);
  MCAuto<MEDFileData> d2=sr->loadInMEDFileDS();
  MEDFileUMesh* m2 = static_cast<MEDFileUMesh*>( d2->getMeshes()->getMeshAtPos(0) );
  MCAuto<MEDCouplingUMesh> group1 = m2->getGroup(0, "Group1");
  CPPUNIT_ASSERT_EQUAL(1,(int)group1->getNumberOfCells());
  MCAuto<MEDCouplingUMesh> group2 = m2->getGroup(0, "Group2");
  CPPUNIT_ASSERT_EQUAL(1,(int)group2->getNumberOfCells());
  MCAuto<MEDCouplingUMesh> grptot = m2->getGroup(0, "Grouptot");
  CPPUNIT_ASSERT_EQUAL(2,(int)grptot->getNumberOfCells());
}

void SauvLoaderTest::testSauv2MedOnA3SubsField()
{
  // read SAUV
  std::string sauvFile = INTERP_TEST::getResourceFile("portico_3subs.sauv", 3);
  MCAuto<SauvReader> sr=SauvReader::New(sauvFile.c_str());
  MCAuto<MEDFileData> d2=sr->loadInMEDFileDS();
  // check mesh
  MEDFileUMesh* m2 = static_cast<MEDFileUMesh*>(d2->getMeshes()->getMeshAtPos(0));
  MCAuto<MEDCouplingUMesh> mesh1d = m2->getMeshAtLevel(0);
  MCAuto<MEDCouplingFieldDouble> length1dField = mesh1d->getMeasureField(0);
  std::cout << "Length of 1d elements: " << length1dField->accumulate(0) << std::endl;
  CPPUNIT_ASSERT_DOUBLES_EQUAL(3, length1dField->accumulate(0), 1e-12);
  // check field
  MCAuto<MEDFileFieldMultiTS> field =
    dynamic_cast<MEDFileFieldMultiTS *>(d2->getFields()->getFieldWithName("CHAM1D"));
  std::cout << "Number of components in field: " << field->getInfo().size() << std::endl;
  CPPUNIT_ASSERT_EQUAL(6,(int)field->getInfo().size());
  std::vector< std::pair<int,int> > timesteps = field->getIterations();

  MCAuto<MEDCouplingFieldDouble> field1d =
    field->getFieldOnMeshAtLevel(ON_GAUSS_NE, timesteps[0].first, timesteps[0].second, 0, m2);

  // Check first component of the field
  // 2 gauss points per element => 12 values
  double values[12] = {
      -7.687500000000e-03,
      -7.687500000000e-03,
      -4.562500000000e-03,
      -4.562500000000e-03,
      -8.208333333333e-03,
      -8.208333333333e-03,
      -6.125000000000e-03,
      -6.125000000000e-03,
      -4.041666666666e-03,
      -4.041666666666e-03,
      -6.111413346910e-07,
      -6.111413346910e-07};

  for (int i=0; i < field1d->getNumberOfTuples(); i++)
  {
    CPPUNIT_ASSERT_DOUBLES_EQUAL( values[i], field1d->getIJ(i, 0), 1e-12 );
  }
}

void SauvLoaderTest::testMed2Sauv()
{
  // read pointe.med
  std::string file = INTERP_TEST::getResourceFile("pointe.med", 3);
  MCAuto<MEDFileData> pointeMed=MEDFileData::New(file.c_str());

  // add 3 faces to pointeMed
  MEDFileUMesh* pointeMedMesh = static_cast<MEDFileUMesh*>(pointeMed->getMeshes()->getMeshAtPos(0));
  MCAuto<MEDCouplingUMesh> pointeM1D = MEDCouplingUMesh::New();
  DataArrayDouble     *coords = pointeMedMesh->getCoords();
  pointeM1D->setCoords( coords );
  pointeM1D->setMeshDimension( 2 );
  pointeM1D->allocateCells( 3 );
  int conn[]=
    {
      0,1,2, 0,1,3, 10,11,12,13
    };
  pointeM1D->insertNextCell( INTERP_KERNEL::NORM_TRI3, 3, conn);
  pointeM1D->insertNextCell( INTERP_KERNEL::NORM_TRI3, 3, conn+3);
  pointeM1D->insertNextCell( INTERP_KERNEL::NORM_QUAD4, 4, conn+6);
  pointeM1D->finishInsertingCells();
  pointeMedMesh->setMeshAtLevel( -1, pointeM1D );
  pointeMed->getMeshes()->setMeshAtPos( 0, pointeMedMesh );

  // add a field on 2 faces to pointeMed
  MCAuto<MEDFileFieldMultiTS> ff1=MEDFileFieldMultiTS::New();
  MCAuto<MEDCouplingFieldDouble> f1=MEDCouplingFieldDouble::New(ON_GAUSS_NE,ONE_TIME);
  f1->setMesh( pointeM1D );
  f1->setName("Field on 2 faces");
  MCAuto<DataArrayDouble> d=DataArrayDouble::New();
  d->alloc(3+4,2);
  d->setInfoOnComponent(0,"sigX [MPa]");
  d->setInfoOnComponent(1,"sigY [GPa]");
  double vals[2*(3+4)] =
    {
      311,312,321,322,331,332,411,412,421,422,431,432,441,442
    };
  std::copy(vals,vals+d->getNbOfElems(),d->getPointer());
  f1->setArray(d);
  MCAuto<DataArrayInt> da=DataArrayInt::New();
  int ids[] =
    {
      0,2
    };
  da->alloc(2,1);
  std::copy(ids,ids+da->getNbOfElems(),da->getPointer());
  da->setName("sup2");
  ff1->appendFieldProfile(f1,pointeMedMesh,-1,da);
  pointeMed->getFields()->pushField( ff1 );

  // remove "fieldnodeint"
  MEDFileFields* pointeFields = pointeMed->getFields();
  for ( int i = 0; i < pointeFields->getNumberOfFields(); ++i )
    {
      MCAuto<MEDFileAnyTypeFieldMultiTS> ts = pointeFields->getFieldAtPos(i);
      if ( std::string("fieldnodeint") == ts->getName())
        {
          pointeFields->destroyFieldAtPos( i );
          break;
        }
    }
  // write pointeMed to SAUV
  const char* sauvFile = "pointe.sauv";
  MCAuto<SauvWriter> sw=SauvWriter::New();
  sw->setMEDFileDS(pointeMed);
  sw->write(sauvFile);

  // read SAUV and check
  MCAuto<SauvReader> sr=SauvReader::New(sauvFile);
  MCAuto<MEDFileData> d2=sr->loadInMEDFileDS();
  CPPUNIT_ASSERT_EQUAL(1,d2->getNumberOfMeshes());
  CPPUNIT_ASSERT_EQUAL(4,d2->getNumberOfFields());
  MEDFileUMesh * m = static_cast<MEDFileUMesh*>( d2->getMeshes()->getMeshAtPos(0) );
  CPPUNIT_ASSERT_EQUAL(std::string("maa1"),std::string(m->getName() ));
  CPPUNIT_ASSERT_EQUAL(3,m->getMeshDimension());
  std::vector<std::string > groups = m->getGroupsNames();
  CPPUNIT_ASSERT_EQUAL(6,(int)groups.size());
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe1") != groups.end() );
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe2") != groups.end() );
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe3") != groups.end() );
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe4") != groups.end() );
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe5") != groups.end() );
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"maa1") != groups.end() );
  CPPUNIT_ASSERT_EQUAL(16,m->getSizeAtLevel(0));
  MCAuto<MEDCouplingMesh> um0 = m->getMeshAtLevel(0);
  CPPUNIT_ASSERT_EQUAL(12, (int)um0->getNumberOfCellsWithType( INTERP_KERNEL::NORM_TETRA4 ));
  CPPUNIT_ASSERT_EQUAL(2,  (int)um0->getNumberOfCellsWithType( INTERP_KERNEL::NORM_PYRA5 ));
  CPPUNIT_ASSERT_EQUAL(2,  (int)um0->getNumberOfCellsWithType( INTERP_KERNEL::NORM_HEXA8 ));
  MCAuto<MEDCouplingMesh> um1 = m->getMeshAtLevel(-1);
  CPPUNIT_ASSERT_EQUAL(1, (int)um1->getNumberOfCellsWithType( INTERP_KERNEL::NORM_TRI3 ));
  //CPPUNIT_ASSERT_EQUAL(2, um1->getNumberOfCellsWithType( INTERP_KERNEL::NORM_TRI3 ));
  MCAuto<MEDCouplingUMesh> pointeUM0 =
    static_cast<MEDCouplingUMesh*>( pointeMedMesh->getMeshAtLevel(0));
  DataArrayDouble *coo = m->getCoords();
  DataArrayDouble *pointeCoo = pointeMedMesh->getCoords();
  CPPUNIT_ASSERT(coo->isEqualWithoutConsideringStr(*pointeCoo,1e-12));
  MCAuto<MEDCouplingFieldDouble> vol = um0->getMeasureField(0);
  MCAuto<MEDCouplingFieldDouble> pointeVol = pointeUM0->getMeasureField(0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vol->accumulate(0), pointeVol->accumulate(0),1e-12);
  // check fields
  // fieldnodedouble
  MCAuto<MEDFileFieldMultiTS> fieldnodedoubleTS1 =
    dynamic_cast<MEDFileFieldMultiTS *>(pointeMed->getFields()->getFieldWithName("fieldnodedouble"));
  MCAuto<MEDFileFieldMultiTS> fieldnodedoubleTS2 =
    dynamic_cast<MEDFileFieldMultiTS *>(d2->getFields()->getFieldWithName("fieldnodedouble"));
  CPPUNIT_ASSERT_EQUAL( fieldnodedoubleTS1->getInfo().size(), fieldnodedoubleTS2->getInfo().size());
  for ( size_t i = 0; i < fieldnodedoubleTS1->getInfo().size(); ++i )
    CPPUNIT_ASSERT_EQUAL( fieldnodedoubleTS1->getInfo()[i], fieldnodedoubleTS2->getInfo()[i]);
  CPPUNIT_ASSERT_EQUAL( fieldnodedoubleTS1->getNumberOfTS(), fieldnodedoubleTS2->getNumberOfTS());
  std::vector< std::pair<int,int> > io1 = fieldnodedoubleTS1->getIterations();
  std::vector< std::pair<int,int> > io2 = fieldnodedoubleTS2->getIterations();
  for ( int i =0; i < fieldnodedoubleTS1->getNumberOfTS(); ++i )
    {
      MCAuto<MEDCouplingFieldDouble> fnd1 =
        fieldnodedoubleTS1->getFieldOnMeshAtLevel(ON_NODES, io1[i].first,io1[i].second,pointeUM0);
      MCAuto<MEDCouplingFieldDouble> fnd2 =
        fieldnodedoubleTS2->getFieldOnMeshAtLevel(ON_NODES, io2[i].first,io2[i].second,um0);
      CPPUNIT_ASSERT( fnd1->getArray()->isEqual( *fnd2->getArray(), 1e-12 ));
    }
  // fieldcelldoublevector
  MCAuto<MEDFileFieldMultiTS> fieldcelldoublevectorTS1 =
    dynamic_cast<MEDFileFieldMultiTS *>(pointeMed->getFields()->getFieldWithName("fieldcelldoublevector"));
  MCAuto<MEDFileFieldMultiTS> fieldcelldoublevectorTS2 =
    dynamic_cast<MEDFileFieldMultiTS *>(d2->getFields()->getFieldWithName("fieldcelldoublevector"));
  CPPUNIT_ASSERT_EQUAL( fieldcelldoublevectorTS1->getInfo().size(), fieldcelldoublevectorTS2->getInfo().size());
  for ( size_t i = 0; i < fieldcelldoublevectorTS1->getInfo().size(); ++i )
    CPPUNIT_ASSERT_EQUAL( fieldcelldoublevectorTS1->getInfo()[i], fieldcelldoublevectorTS2->getInfo()[i]);
  CPPUNIT_ASSERT_EQUAL( fieldcelldoublevectorTS1->getNumberOfTS(), fieldcelldoublevectorTS2->getNumberOfTS());
  io1 = fieldcelldoublevectorTS1->getIterations();
  io2 = fieldcelldoublevectorTS2->getIterations();
  for ( int i =0; i < fieldcelldoublevectorTS1->getNumberOfTS(); ++i )
    {
      MCAuto<MEDCouplingFieldDouble> fnd1 =
        fieldcelldoublevectorTS1->getFieldOnMeshAtLevel(ON_CELLS, io1[i].first,io1[i].second,pointeUM0);
      MCAuto<MEDCouplingFieldDouble> fnd2 =
        fieldcelldoublevectorTS2->getFieldOnMeshAtLevel(ON_CELLS, io2[i].first,io2[i].second,um0);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( fnd1->accumulate(0), fnd2->accumulate(0), 1e-12 );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( fnd1->accumulate(1), fnd2->accumulate(1), 1e-12 );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( fnd1->accumulate(2), fnd2->accumulate(2), 1e-12 );
    }
  // "Field on 2 faces"
  MCAuto<MEDFileFieldMultiTS> fieldOnFaces =
    dynamic_cast<MEDFileFieldMultiTS *>(d2->getFields()->getFieldWithName(f1->getName().c_str()));
  io1 = fieldOnFaces->getIterations();
  MCAuto<MEDCouplingFieldDouble> fof =
    fieldOnFaces->getFieldOnMeshAtLevel(f1->getTypeOfField(),io1[0].first,io1[0].second,um1);
  CPPUNIT_ASSERT( d->isEqual( *fof->getArray(), 1e-12 ));
}

void SauvLoaderTest::testCellsWithLingNames()
{
  // test IMP 3285: [CEA 1778] SauvReader: only keep the meshes named in the table MED_MAIL
  std::string file = INTERP_TEST::getResourceFile("test_MED_MAIL.sauv", 3);
  MCAuto<SauvReader> sr=SauvReader::New(file.c_str());
  MCAuto<MEDFileData> d2=sr->loadInMEDFileDS();
  // check that the mesh contains
  // - Nombre de noeuds : 74
  // - Nombre de mailles de type MED_TRIA3 : 6
  // - Nombre de mailles de type MED_QUAD4 : 43
  // - Nombre de mailles de type MED_HEXA8 : 24
  // - Nombre de mailles de type MED_PENTA6 : 3
  MEDFileUMesh* m = static_cast<MEDFileUMesh*>( d2->getMeshes()->getMeshAtPos(0));
  CPPUNIT_ASSERT_EQUAL(6,  m->getNumberOfCellsWithType( INTERP_KERNEL::NORM_TRI3 ));
  CPPUNIT_ASSERT_EQUAL(43, m->getNumberOfCellsWithType( INTERP_KERNEL::NORM_QUAD4 ));
  CPPUNIT_ASSERT_EQUAL(24, m->getNumberOfCellsWithType( INTERP_KERNEL::NORM_HEXA8 ));
  CPPUNIT_ASSERT_EQUAL(3,  m->getNumberOfCellsWithType( INTERP_KERNEL::NORM_PENTA6 ));
}

void SauvLoaderTest::tearDown()
{
  const int nbFilesToRemove = 3;
  const char* fileToRemove[nbFilesToRemove] = { "allPillesTest.med", "pointe.sauv", "mesh_with_void_family.sauv" };
  for ( int i = 0; i < nbFilesToRemove; ++i )
  {
#ifdef WIN32
    if (GetFileAttributes(fileToRemove[i]) != INVALID_FILE_ATTRIBUTES)
#else
      if (access(fileToRemove[i], F_OK) == 0)
#endif
        remove(fileToRemove[i]);
  }
}
