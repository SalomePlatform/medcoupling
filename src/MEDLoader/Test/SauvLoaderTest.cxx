// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifdef WNT
#include <windows.h>
#endif

#include <vector>
#include <string>

using namespace ParaMEDMEM;

void SauvLoaderTest::testSauv2Med()
{
  // read a file containing all types of readable piles
  std::string file = getResourceFile("allPillesTest.sauv");
  MEDCouplingAutoRefCountObjectPtr<SauvReader> sr=SauvReader::New(file.c_str());
  MEDCouplingAutoRefCountObjectPtr<MEDFileData> d2=sr->loadInMEDFileDS();
  // write MED
  d2->write("allPillesTest.med",0);
  // check 
  CPPUNIT_ASSERT_EQUAL(1,d2->getNumberOfMeshes());
  CPPUNIT_ASSERT_EQUAL(8+97,d2->getNumberOfFields());
  MEDFileMesh * m = d2->getMeshes()->getMeshAtPos(0);
  CPPUNIT_ASSERT_EQUAL(17,int(m->getGroupsNames().size()));
}

void SauvLoaderTest::testMed2Sauv()
{
  // read pointe.med
  std::string file = getResourceFile("pointe.med");
  MEDCouplingAutoRefCountObjectPtr<MEDFileData> pointeMed=MEDFileData::New(file.c_str());

  // add 3 faces to pointeMed
  MEDFileUMesh* pointeMedMesh = static_cast<MEDFileUMesh*>(pointeMed->getMeshes()->getMeshAtPos(0));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> pointeM1D = MEDCouplingUMesh::New();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble>     coords = pointeMedMesh->getCoords();
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
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> ff1=MEDFileFieldMultiTS::New();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> f1=MEDCouplingFieldDouble::New(ON_GAUSS_NE,ONE_TIME);
  f1->setMesh( pointeM1D );
  f1->setName("Field on 2 faces");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> d=DataArrayDouble::New();
  d->alloc(3+4,2);
  d->setInfoOnComponent(0,"sigX [MPa]");
  d->setInfoOnComponent(1,"sigY [GPa]");
  double vals[2*(3+4)] =
    {
      311,312,321,322,331,332,411,412,421,422,431,432,441,442
    };
  std::copy(vals,vals+d->getNbOfElems(),d->getPointer());
  f1->setArray(d);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> da=DataArrayInt::New();
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
      MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> ts = pointeFields->getFieldAtPos(i);
      if ( std::string("fieldnodeint") == ts->getName())
        {
          pointeFields->destroyFieldAtPos( i );
          break;
        }
    }
  // write pointeMed to SAUV
  const char* sauvFile = "pointe.sauv";
  MEDCouplingAutoRefCountObjectPtr<SauvWriter> sw=SauvWriter::New();
  sw->setMEDFileDS(pointeMed);
  sw->write(sauvFile);

  // read SAUV and check
  MEDCouplingAutoRefCountObjectPtr<SauvReader> sr=SauvReader::New(sauvFile);
  MEDCouplingAutoRefCountObjectPtr<MEDFileData> d2=sr->loadInMEDFileDS();
  CPPUNIT_ASSERT_EQUAL(1,d2->getNumberOfMeshes());
  CPPUNIT_ASSERT_EQUAL(4,d2->getNumberOfFields());
  MEDFileUMesh * m = static_cast<MEDFileUMesh*>( d2->getMeshes()->getMeshAtPos(0) );
  CPPUNIT_ASSERT_EQUAL(std::string("maa1"),std::string(m->getName() ));
  CPPUNIT_ASSERT_EQUAL(3,m->getMeshDimension());
  std::vector<std::string > groups = m->getGroupsNames();
  CPPUNIT_ASSERT_EQUAL(5,(int)groups.size());
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe1") != groups.end() );
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe2") != groups.end() );
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe3") != groups.end() );
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe4") != groups.end() );
  CPPUNIT_ASSERT( std::find(groups.begin(),groups.end(),"groupe5") != groups.end() );
  CPPUNIT_ASSERT_EQUAL(16,m->getSizeAtLevel(0));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> um0 = m->getGenMeshAtLevel(0);
  CPPUNIT_ASSERT_EQUAL(12, um0->getNumberOfCellsWithType( INTERP_KERNEL::NORM_TETRA4 ));
  CPPUNIT_ASSERT_EQUAL(2,  um0->getNumberOfCellsWithType( INTERP_KERNEL::NORM_PYRA5 ));
  CPPUNIT_ASSERT_EQUAL(2,  um0->getNumberOfCellsWithType( INTERP_KERNEL::NORM_HEXA8 ));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> um1 = m->getGenMeshAtLevel(-1);
  CPPUNIT_ASSERT_EQUAL(2, um1->getNumberOfCellsWithType( INTERP_KERNEL::NORM_TRI3 ));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> pointeUM0 =
    static_cast<MEDCouplingUMesh*>( pointeMedMesh->getGenMeshAtLevel(0));
  MEDCouplingAutoRefCountObjectPtr< DataArrayDouble > coo = m->getCoords();
  MEDCouplingAutoRefCountObjectPtr< DataArrayDouble > pointeCoo = pointeMedMesh->getCoords();
  CPPUNIT_ASSERT(coo->isEqualWithoutConsideringStr(*pointeCoo,1e-12));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> vol = um0->getMeasureField(0);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> pointeVol = pointeUM0->getMeasureField(0);
  CPPUNIT_ASSERT_DOUBLES_EQUAL( vol->accumulate(0), pointeVol->accumulate(0),1e-12);
  // check fields
  // fieldnodedouble
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> fieldnodedoubleTS1 =
    pointeMed->getFields()->getFieldWithName("fieldnodedouble");
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> fieldnodedoubleTS2 =
    d2->getFields()->getFieldWithName("fieldnodedouble");
  CPPUNIT_ASSERT_EQUAL( fieldnodedoubleTS1->getInfo().size(), fieldnodedoubleTS2->getInfo().size());
  for ( size_t i = 0; i < fieldnodedoubleTS1->getInfo().size(); ++i )
    CPPUNIT_ASSERT_EQUAL( fieldnodedoubleTS1->getInfo()[i], fieldnodedoubleTS2->getInfo()[i]);
  CPPUNIT_ASSERT_EQUAL( fieldnodedoubleTS1->getNumberOfTS(), fieldnodedoubleTS2->getNumberOfTS());
  std::vector< std::pair<int,int> > io1 = fieldnodedoubleTS1->getIterations();
  std::vector< std::pair<int,int> > io2 = fieldnodedoubleTS2->getIterations();
  for ( int i =0; i < fieldnodedoubleTS1->getNumberOfTS(); ++i )
    {
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> fnd1 =
        fieldnodedoubleTS1->getFieldOnMeshAtLevel(ON_NODES, io1[i].first,io1[i].second,pointeUM0);
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> fnd2 =
        fieldnodedoubleTS2->getFieldOnMeshAtLevel(ON_NODES, io2[i].first,io2[i].second,um0);
      CPPUNIT_ASSERT( fnd1->getArray()->isEqual( *fnd2->getArray(), 1e-12 ));
    }
  // fieldcelldoublevector
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> fieldcelldoublevectorTS1 =
    pointeMed->getFields()->getFieldWithName("fieldcelldoublevector");
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> fieldcelldoublevectorTS2 =
    d2->getFields()->getFieldWithName("fieldcelldoublevector");
  CPPUNIT_ASSERT_EQUAL( fieldcelldoublevectorTS1->getInfo().size(), fieldcelldoublevectorTS2->getInfo().size());
  for ( size_t i = 0; i < fieldcelldoublevectorTS1->getInfo().size(); ++i )
    CPPUNIT_ASSERT_EQUAL( fieldcelldoublevectorTS1->getInfo()[i], fieldcelldoublevectorTS2->getInfo()[i]);
  CPPUNIT_ASSERT_EQUAL( fieldcelldoublevectorTS1->getNumberOfTS(), fieldcelldoublevectorTS2->getNumberOfTS());
  io1 = fieldcelldoublevectorTS1->getIterations();
  io2 = fieldcelldoublevectorTS2->getIterations();
  for ( int i =0; i < fieldcelldoublevectorTS1->getNumberOfTS(); ++i )
    {
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> fnd1 =
        fieldcelldoublevectorTS1->getFieldOnMeshAtLevel(ON_CELLS, io1[i].first,io1[i].second,pointeUM0);
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> fnd2 =
        fieldcelldoublevectorTS2->getFieldOnMeshAtLevel(ON_CELLS, io2[i].first,io2[i].second,um0);
      CPPUNIT_ASSERT_DOUBLES_EQUAL( fnd1->accumulate(0), fnd2->accumulate(0), 1e-12 );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( fnd1->accumulate(1), fnd2->accumulate(1), 1e-12 );
      CPPUNIT_ASSERT_DOUBLES_EQUAL( fnd1->accumulate(2), fnd2->accumulate(2), 1e-12 );
    }
  // "Field on 2 faces"
  MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTS> fieldOnFaces =
    d2->getFields()->getFieldWithName(f1->getName());
  io1 = fieldOnFaces->getIterations();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> fof =
    fieldOnFaces->getFieldOnMeshAtLevel(f1->getTypeOfField(),io1[0].first,io1[0].second,um1);
  CPPUNIT_ASSERT( d->isEqual( *fof->getArray(), 1e-12 ));
}

void SauvLoaderTest::tearDown()
{
  const int nbFilesToRemove = 2;
  const char* fileToRemove[nbFilesToRemove] = { "allPillesTest.med", "pointe.sauv" };
  for ( int i = 0; i < nbFilesToRemove; ++i )
    {
#ifdef WNT
      if (GetFileAttributes(fileToRemove[i]) != INVALID_FILE_ATTRIBUTES)
#else
        if (access(fileToRemove[i], F_OK) == 0)
#endif
      remove(fileToRemove[i]);
  }
}

std::string SauvLoaderTest::getResourceFile( const std::string& filename )
{
  std::string resourceFile = "";

  if ( getenv("top_srcdir") ) {
    // we are in 'make check' step
    resourceFile = getenv("top_srcdir");
    resourceFile += "/resources/";
  }
  else if ( getenv("MED_ROOT_DIR") ) {
    // use MED_ROOT_DIR env.var
    resourceFile = getenv("MED_ROOT_DIR");
    resourceFile += "/share/salome/resources/med/";
  }
  resourceFile += filename;
#ifdef WNT
  std::string fixedpath = resourceFile;
  for ( int i=0; i < fixedpath.length(); ++i )
    if (fixedpath[i] == '/')
      fixedpath[i] = '\\';
  return fixedpath;
#endif
  return resourceFile;
}
