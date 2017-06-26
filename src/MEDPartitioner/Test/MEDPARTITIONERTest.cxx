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

#include "MEDPARTITIONERTest.hxx"

#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "CellModel.hxx"
#include "MEDFileMesh.hxx"
#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMappedExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingMultiFields.hxx"

#include <cppunit/TestAssert.h>

#include <sstream>
#include <fstream>
#include <cmath>
#include <list>
#include <stdexcept>
#include <cstdlib>
#include <vector>
#ifdef WIN32
#include<direct.h>
#define getcwd _getcwd
#else
#include <unistd.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;
using namespace MEDCoupling;
using namespace MEDPARTITIONER;

void MEDPARTITIONERTest::setSize(int ni, int nj, int nk)
{
  this->_ni=ni;  //nb of hexa9
  this->_nj=nj;
  this->_nk=nk;
  this->_ntot=_ni*_nj*_nk;
  string ijk=IntToStr(ni)+"x"+IntToStr(nj)+"x"+IntToStr(nk);
  this->_file_name="tmp_testMesh_"+ijk+".med";
  this->_file_name_with_faces="tmp_testMeshWithFaces_"+ijk+".med";
  string ij=IntToStr(ni)+"x"+IntToStr(nj);
  this->_file_name2="tmp_testMesh_"+ij+".med";
  this->_mesh_name="testMesh";
}

void MEDPARTITIONERTest::setSmallSize()
{
  setSize(2,3,5); //nb of hexa9
}

void MEDPARTITIONERTest::setMedianSize()
{
  setSize(20,30,50); //nb of hexa9
}

void MEDPARTITIONERTest::setbigSize()
{
  setSize(200,300,500); //nb of hexa9
}

std::string MEDPARTITIONERTest::getPartitionerExe() const
{
  std::string execName;
  if ( getenv("MEDCOUPLING_ROOT_DIR") )
    {
      execName=getenv("MEDCOUPLING_ROOT_DIR");  //.../INSTALL/MED
      execName+="/bin/medpartitioner";
      std::ifstream my_file(execName.c_str());
      if (my_file.good())
        return execName;
    }
  execName = getcwd(NULL, 0);
  execName += "/../../MEDPartitioner/medpartitioner";
  if (! std::ifstream(execName.c_str()))
    CPPUNIT_FAIL("Can't find medpartitioner, please set MEDCOUPLING_ROOT_DIR");
  return execName;
}

// ============================================================================
/*!
 *  Set up the environment called at every CPPUNIT_TEST ()
 */
// ============================================================================
void MEDPARTITIONERTest::setUp()
{
  this->_verbose=0;
#if defined(HAVE_MPI)
  if (MyGlobals::_Rank==-1)  //do once only
    {
      MPI_Init(0,0);
      MPI_Comm_size(MPI_COMM_WORLD, &MyGlobals::_World_Size);
      MPI_Comm_rank(MPI_COMM_WORLD, &MyGlobals::_Rank);
    }
#else
  //sequential : no MPI
  MyGlobals::_World_Size=1;
  MyGlobals::_Rank=0;
#endif

  if (_verbose>10)
    {
#if defined(HAVE_MPI)
      cout<<"\ndefined(HAVE_MPI)"<<endl;
#else
      cout<<"\nNOT defined(HAVE_MPI)"<<endl;
#endif
#if defined(MED_ENABLE_PARMETIS)
      cout<<"defined(MED_ENABLE_PARMETIS)"<<endl;
#else
      cout<<"NOT defined(MED_ENABLE_PARMETIS)"<<endl;
#endif
#if defined(MED_ENABLE_METIS)
      cout<<"defined(MED_ENABLE_METIS)"<<endl;
#else
      cout<<"NOT defined(MED_ENABLE_METIS)"<<endl;
#endif
#if defined(MED_ENABLE_SCOTCH)
      cout<<"defined(MED_ENABLE_SCOTCH)"<<endl;
#else
      cout<<"NOT defined(MED_ENABLE_SCOTCH)"<<endl;
#endif
    }
}

// ============================================================================
/*!
 *  - delete
 */
// ============================================================================
void MEDPARTITIONERTest::tearDown()
{
}

MEDCoupling::MEDCouplingUMesh * MEDPARTITIONERTest::buildCUBE3DMesh()
//only hexa8
{
  vector<int> conn;
  vector<double> coor;
  for (int k=0; k<=_nk; k++)
    for (int j=0; j<=_nj; j++)
      for (int i=0; i<=_ni; i++)
        {
          coor.push_back(i+.1);
          coor.push_back(j+.2);
          coor.push_back(k+.3);
        }
  int ii;
  for (int k=0; k<_nk; k++)
    for (int j=0; j<_nj; j++)
      for (int i=0; i<_ni; i++)
        {
          ii=i + j*(_ni+1) + k*(_ni+1)*(_nj+1);
          conn.push_back(ii);
          conn.push_back(ii+1);
          ii=ii + _ni + 2 ;
          conn.push_back(ii);
          conn.push_back(ii-1);

          ii=i + j*(_ni+1) + (k+1)*(_ni+1)*(_nj+1);
          conn.push_back(ii);
          conn.push_back(ii+1);
          ii=ii + _ni + 2 ;
          conn.push_back(ii);
          conn.push_back(ii-1);
        }

  /*
  if (_verbose)  //only for debug
    {
      cout<< "\nnb coor " << (_ni+1)*(_nj+1)*(_nk+1)*3 << " " << coor.size() << endl;
      for (int i=0; i<(int)coor.size(); i++)
        cout << coor[i] << " ";
      cout << endl;
      cout << "\nnb conn " << (_ni)*(_nj)*(_nk)*8 << " " << conn.size() << endl;
      for (int i=0; i<(int)conn.size(); i=i+8)
        {
          for (int j=0; j<8; j++)
            cout << conn[i+j] << " ";
          cout << endl;
        }
      cout << endl;
    }
  */

  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(3);
  int nbc=conn.size()/8; //nb of cells
  int nbv=coor.size()/3; //nb of vertices
  mesh->allocateCells(nbc);
  for(int i=0; i<nbc; i++)
    {
      int onehexa[8];
      std::copy(conn.begin()+i*8,conn.begin()+(i+1)*8,onehexa);
      if (false) //(_verbose)
        {
          for (int j=0; j<8; j++) cout<<onehexa[j]<<" ";
          cout<<endl;
        }
      mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,onehexa);
    }
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbv,3);
  std::copy(coor.begin(),coor.end(),myCoords->getPointer());
  mesh->setCoords(myCoords);
  mesh->setName(_mesh_name.c_str());
  myCoords->decrRef();
  mesh->checkConsistencyLight();
  return mesh;
}

MEDCoupling::MEDCouplingUMesh * MEDPARTITIONERTest::buildCARRE3DMesh()
//only quad4 in oblique (k=j)
{
  vector<int> conn;
  vector<double> coor;
  for (int j=0; j<=_nj; j++)
    for (int i=0; i<=_ni; i++)
      {
        int k=j;
        coor.push_back(i+.1);
        coor.push_back(j+.2);
        coor.push_back(k+.3);
      }
  int ii;
  int k=0;
  for (int j=0; j<_nj; j++)
    for (int i=0; i<_ni; i++)
      {
        ii=i + j*(_ni+1) + k*(_ni+1)*(_nj+1);
        conn.push_back(ii);
        conn.push_back(ii+1);
        ii=ii + _ni + 2 ;
        conn.push_back(ii);
        conn.push_back(ii-1);
      }

  if (false) //(_verbose)
    {
      cout<<"\nnb coor "<<(_ni+1)*(_nj+1)*3<<" "<<coor.size()<<endl;
      for (int i=0; i<(int)coor.size(); i++)
        cout << coor[i] << " ";
      cout<<endl;
      cout<<"\nnb conn "<<(_ni)*(_nj)*4<<" "<<conn.size()<<endl;
      for (int i=0; i<(int)conn.size(); i=i+4)
        {
          for (int j=0; j<4; j++) cout<<conn[i+j]<<" ";
          cout<<endl;
        }
      cout<<endl;
    }

  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  int nbc=conn.size()/4; //nb of cells
  int nbv=coor.size()/3; //nb of vertices
  mesh->allocateCells(nbc);
  for(int i=0; i<nbc; i++)
    {
      int onequa[4];
      std::copy(conn.begin()+i*4,conn.begin()+(i+1)*4,onequa);
      if (false) //(_verbose)
        {
          for (int j=0; j<4; j++) cout<<onequa[j]<<" ";
          cout<<endl;
        }
      mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,onequa);
    }
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbv,3);
  std::copy(coor.begin(),coor.end(),myCoords->getPointer());
  mesh->setCoords(myCoords);
  mesh->setName(_mesh_name.c_str());
  myCoords->decrRef();
  mesh->checkConsistencyLight();
  return mesh;
}

MEDCoupling::MEDCouplingUMesh * MEDPARTITIONERTest::buildFACE3DMesh()
//only quad4 on a global face of the CUBE3D (k=0)
{
  vector<int> conn;
  vector<double> coor;
  for (int j=0; j<=_nj; j++)
    for (int i=0; i<=_ni; i++)
      {
        int k=0;
        coor.push_back(i+.1);
        coor.push_back(j+.2);
        coor.push_back(k+.3);
      }
  int ii;
  int k=0;
  for (int j=0; j<_nj; j++)
    for (int i=0; i<_ni; i++)
      {
        ii=i + j*(_ni+1) + k*(_ni+1)*(_nj+1);
        conn.push_back(ii);
        conn.push_back(ii+1);
        ii=ii + _ni + 2 ;
        conn.push_back(ii);
        conn.push_back(ii-1);
      }

  if (false) //(_verbose)
    {
      cout<<"\nnb coor "<<(_ni+1)*(_nj+1)*3<<" "<<coor.size()<<endl;
      for (int i=0; i<(int)coor.size(); i++)
        cout << coor[i] << " ";
      cout<<endl;
      cout<<"\nnb conn "<<(_ni)*(_nj)*4<<" "<<conn.size()<<endl;
      for (int i=0; i<(int)conn.size(); i=i+4)
        {
          for (int j=0; j<4; j++)
            cout << conn[i+j] << " ";
          cout << endl;
        }
      cout << endl;
    }

  MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
  mesh->setMeshDimension(2);
  int nbc=conn.size()/4; //nb of cells
  int nbv=coor.size()/3; //nb of vertices
  mesh->allocateCells(nbc);
  for(int i=0; i<nbc; i++)
    {
      int onequa[4];
      std::copy(conn.begin()+i*4,conn.begin()+(i+1)*4,onequa);
      if (false) //(_verbose)
        {
          for (int j=0; j<4; j++) cout<<onequa[j]<<" ";
          cout<<endl;
        }
      mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,onequa);
    }
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbv,3);
  std::copy(coor.begin(),coor.end(),myCoords->getPointer());
  mesh->setCoords(myCoords);
  mesh->setName(_mesh_name.c_str());
  myCoords->decrRef();
  mesh->checkConsistencyLight();
  return mesh;
}

MEDCouplingFieldDouble * MEDPARTITIONERTest::buildVecFieldOnCells(string myfileName)
{
  //int ni=2,nj=3,nk=5; //nb of hexa9
  vector<double> field;
  for (int k=0; k<_nk; k++)
    for (int j=0; j<_nj; j++)
      for (int i=0; i<_ni; i++)
        {
          field.push_back(i+.1);
          field.push_back(j+.2);
          field.push_back(k+.3);
        }

  MEDCouplingUMesh *mesh=ReadUMeshFromFile(myfileName.c_str(),_mesh_name.c_str(),0);
  int nbOfCells=mesh->getNumberOfCells();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME);
  f1->setName("VectorFieldOnCells");
  f1->setDescription("DescriptionOfFieldOnCells"); //not saved in file?
  f1->setMesh(mesh);
  DataArrayDouble *myField=DataArrayDouble::New();
  myField->alloc(nbOfCells,3);
  std::copy(field.begin(),field.end(),myField->getPointer());
  f1->setArray(myField);
  myField->setInfoOnComponent(0,"vx");
  myField->setInfoOnComponent(1,"vy");
  myField->setInfoOnComponent(2,"vz");
  myField->decrRef();
  f1->setTime(2.,0,1);
  f1->checkConsistencyLight();
  mesh->decrRef();
  return f1;
}

MEDCouplingFieldDouble * MEDPARTITIONERTest::buildVecFieldOnNodes()
{
  //int ni=2,nj=3,nk=5; //nb of hexa9
  vector<double> field;
  for (int k=0; k<=_nk; k++)
    for (int j=0; j<=_nj; j++)
      for (int i=0; i<=_ni; i++)
        {
          field.push_back(i+.1);
          field.push_back(j+.2);
          field.push_back(k+.3);
        }

  MEDCouplingUMesh *mesh=ReadUMeshFromFile(_file_name.c_str(),_mesh_name.c_str(),0);
  int nbOfNodes=mesh->getNumberOfNodes();
  MEDCouplingFieldDouble *f1=MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME);
  f1->setName("VectorFieldOnNodes");
  f1->setDescription("DescriptionOfFieldOnNodes"); //not saved in file?
  f1->setMesh(mesh);
  DataArrayDouble *myField=DataArrayDouble::New();
  myField->alloc(nbOfNodes,3);
  std::copy(field.begin(),field.end(),myField->getPointer());
  f1->setArray(myField);
  myField->setInfoOnComponent(0,"vx");
  myField->setInfoOnComponent(1,"vy");
  myField->setInfoOnComponent(2,"vz");
  myField->decrRef();
  f1->setTime(2.,0,1);
  f1->checkConsistencyLight();
  mesh->decrRef();
  return f1;
}


void MEDPARTITIONERTest::createTestMeshWithoutField()
{
  {
    MEDCouplingUMesh * mesh = buildCUBE3DMesh();
    WriteUMesh(_file_name.c_str(),mesh,true);
    if (_verbose) cout<<endl<<_file_name<<" created"<<endl;
    if (_ntot<1000000) //too long
      {
        MEDCouplingUMesh *mesh_rw=ReadUMeshFromFile(_file_name.c_str(),mesh->getName().c_str(),0);
        if (_verbose) cout<<_file_name<<" reread"<<endl;
        CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
        mesh_rw->decrRef();
      }
    mesh->decrRef();
  }

  {
    vector<const MEDCoupling::MEDCouplingUMesh*> meshes;
    MEDCouplingUMesh * mesh1 = buildCUBE3DMesh();
    MEDCouplingUMesh * mesh2 = buildFACE3DMesh();
    mesh1->setName("testMesh");
    mesh2->setName("theFaces");
    mesh2->tryToShareSameCoordsPermute(*mesh1, 1e-9);
    mesh2->checkConsistencyLight();
    mesh1->checkConsistencyLight();
    meshes.push_back(mesh1);
    meshes.push_back(mesh2);
    WriteUMeshes(_file_name_with_faces.c_str(), meshes, true);

    MEDCoupling::MEDFileUMesh* mfm=MEDCoupling::MEDFileUMesh::New(_file_name_with_faces.c_str(), mesh1->getName().c_str());
    DataArrayInt* FacesFam=DataArrayInt::New();
    FacesFam->alloc(mfm->getSizeAtLevel(-1),1);
    FacesFam->fillWithValue(-1);
    DataArrayInt* CellsFam=DataArrayInt::New();
    CellsFam->alloc(mfm->getSizeAtLevel(0),1);
    CellsFam->fillWithValue(1);
    mfm->setFamilyFieldArr(-1,FacesFam);
    mfm->setFamilyFieldArr(0,CellsFam);
    map<string,int> theFamilies;
    theFamilies["FAMILLE_ZERO"]=0;
    theFamilies["FamilyFaces"]=-1;
    theFamilies["FamilyCells"]=1;
    map<string, vector<string> > theGroups;
    theGroups["GroupFaces"].push_back("FamilyFaces");
    theGroups["GroupCells"].push_back("FamilyCells");
    mfm->setFamilyInfo(theFamilies);
    mfm->setGroupInfo(theGroups);
    mfm->write(_file_name_with_faces.c_str(),0);
    FacesFam->decrRef();
    CellsFam->decrRef();

    /*ce truc marche pas!
      MEDCoupling::MEDFileUMesh* mfm=MEDCoupling::MEDFileUMesh::New(_file_name_with_faces.c_str(), mesh1->getName());
      vector<const MEDCoupling::MEDCouplingUMesh*> ms;
      ms.push_back(mesh2);
      mfm->setGroupsFromScratch(-1, ms);
      mfm->write(_file_name_with_faces.c_str(),0);
    */

    if (_verbose) cout<<endl<<_file_name_with_faces<<" created"<<endl;
    if (_ntot<1000000) //too long
      {
        MEDCouplingUMesh *mesh_rw=ReadUMeshFromFile(_file_name_with_faces.c_str(),mesh1->getName().c_str(),0);
        if (_verbose) cout<<_file_name_with_faces<<" reread"<<endl;
        CPPUNIT_ASSERT(mesh1->isEqual(mesh_rw,1e-12));
        mesh_rw->decrRef();
      }
    mesh1->decrRef();
    mesh2->decrRef();
    mfm->decrRef();
  }

  {
    MEDCouplingUMesh * mesh = buildCARRE3DMesh();
    WriteUMesh(_file_name2.c_str(),mesh,true);
    if (_verbose) cout<<endl<<_file_name2<<" created"<<endl;
    MEDCouplingUMesh *mesh_rw=ReadUMeshFromFile(_file_name2.c_str(),mesh->getName().c_str(),0);
    if (_verbose) cout<<_file_name2<<" reread"<<endl;
    CPPUNIT_ASSERT(mesh->isEqual(mesh_rw,1e-12));
    mesh_rw->decrRef();
    mesh->decrRef();
  }
}

/*
create a set of nbx*nby*nbz files mesh of ni*ny*nz cells
*/
void MEDPARTITIONERTest::createHugeTestMesh(int ni, int nj, int nk, int nbx, int nby, int nbz, int nbTarget)
{
  setSize(ni,nj,nk);
  _nb_target_huge=nbTarget;
  MEDCouplingUMesh * mesh = buildCUBE3DMesh();
  //int nbx=1, nby=1, nbz=2;
  std::vector< double > cooDep,cooFin;
  mesh->getCoordinatesOfNode(0, cooDep);
  mesh->getCoordinatesOfNode(mesh->getNumberOfNodes()-1, cooFin);
  //cout<<endl<<cooDep[0]<<" "<<cooDep[1]<<" "<<cooDep[2]<<endl;
  //cout<<cooFin[0]<<" "<<cooFin[1]<<" "<<cooFin[2]<<endl;

  string tagXml="\
<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n \
<root>\n \
  <version maj=\"2\" min=\"3\" ver=\"1\"/>\n \
  <description what=\"\" when=\"YYMMDDHHmm\"/>\n \
  <content>\n \
    <mesh name=\"testMesh\"/>\n \
  </content>\n \
  <splitting>\n \
    <subdomain number=\"$subdomainNumber\"/>\n \
    <global_numbering present=\"no\"/>\n \
  </splitting>\n \
  <files>\n$tagSubfile \
  </files>\n \
  <mapping>\n$tagMesh \
  </mapping>\n \
</root>\n";

  string tagSubfiles, tagSubfile="\
    <subfile id=\"$xyz\">\n \
      <name>$fileName</name>\n \
      <machine>localhost</machine>\n \
    </subfile>\n";
  string tagMeshes, tagMesh="\
    <mesh name=\"testMesh\">\n \
      <chunk subdomain=\"$xyz\">\n \
        <name>testMesh</name>\n \
      </chunk>\n \
    </mesh>\n";

  int xyz=1;
  string sxyz;
  DataArrayDouble* coordsInit=mesh->getCoords()->deepCopy();
  double* ptrInit=coordsInit->getPointer();
  double deltax=cooFin[0]-cooDep[0];
  double deltay=cooFin[1]-cooDep[1];
  double deltaz=cooFin[2]-cooDep[2];

  double dz=0.;
  for (int z=0; z<nbz; z++)
    {
      double dy=0.;
      for (int y=0; y<nby; y++)
        {
          double dx=0.;
          for (int x=0; x<nbx; x++)
            {
              string fileName;
              sxyz=IntToStr(xyz);
              fileName="tmp_testMeshHuge_"+IntToStr(_ni)+"x"+IntToStr(_nj)+"x"+IntToStr(_nk)+"_"+sxyz+".med";

              DataArrayDouble* coords=mesh->getCoords();
              //int nbOfComp=coords->getNumberOfComponents();  //be 3D
              int nbOfTuple=coords->getNumberOfTuples();
              double* ptr=coords->getPointer();
              double* ptrini=ptrInit;
              for (int i=0; i<nbOfTuple; i++)
                {
                  *ptr=(*ptrini)+dx; ptr++; ptrini++; //be 3D
                  *ptr=(*ptrini)+dy; ptr++; ptrini++;
                  *ptr=(*ptrini)+dz; ptr++; ptrini++;
                }

              WriteUMesh(fileName.c_str(),mesh,true);

              tagSubfiles+=tagSubfile;
              tagSubfiles.replace(tagSubfiles.find("$xyz"),4,sxyz);
              tagSubfiles.replace(tagSubfiles.find("$fileName"),9,fileName);

              tagMeshes+=tagMesh;
              tagMeshes.replace(tagMeshes.find("$xyz"),4,sxyz);
              xyz++;
              dx+=deltax;
            }
          dy+=deltay;
        }
      dz+=deltaz;
    }
  coordsInit->decrRef();

  tagXml.replace(tagXml.find("$subdomainNumber"),16,sxyz);
  tagXml.replace(tagXml.find("$tagSubfile"),11,tagSubfiles);
  tagXml.replace(tagXml.find("$tagMesh"),8,tagMeshes);

  string nameFileXml;
  _file_name_huge_xml="tmp_testMeshHuge_"+IntToStr(_ni)+"x"+IntToStr(_nj)+"x"+IntToStr(_nk)+"_"+sxyz+".xml";
  std::ofstream f(_file_name_huge_xml.c_str());
  f<<tagXml;
  f.close();
  //cout<<"\n"<<tagXml<<endl;
  if (_verbose)
    cout<<endl<<nameFileXml<<" created"<<endl;
  mesh->decrRef();
}

void MEDPARTITIONERTest::createTestMeshWithVecFieldOnCells()
{
  {
    string name=_file_name;
    MEDCouplingFieldDouble *f1=buildVecFieldOnCells(name);
    name.replace(name.find(".med"),4,"_WithVecFieldOnCells.med");
    WriteField(name.c_str(),f1,true);
    f1->setTime(3.,1,1);  //time,it,order
    f1->applyFunc("x/2.");
    WriteField(name.c_str(),f1,false);
    if (_verbose) cout<<endl<<name<<" created"<<endl;
    if (_ntot<1000000) //too long
      {
        MCAuto<MEDCouplingFieldDouble> f2=ReadFieldCell(name.c_str(),f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),0,1);
        //DataArrayDouble *res=f2->getArray();
        if (_verbose) cout<<name<<" reread"<<endl;
        //CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12));
      }
    f1->decrRef();
  }
  {
    string name=_file_name;
    MEDCouplingFieldDouble *f1=buildVecFieldOnCells(name);
    name.replace(name.find(".med"),4,"_WithVecFieldOnGaussNe.med");
    MEDCouplingFieldDouble *f3=MEDCouplingFieldDouble::New(ON_GAUSS_NE,ONE_TIME);
    f3->setMesh(f1->getMesh());
    //cout<<"\nNumberOfMeshPlacesExpected "<<f3->getNumberOfMeshPlacesExpected()<<" "
    //                                     /*<<getNumberOfTuples(f1->getMesh())<<" "*/
    //                                     <<f3->getMesh()->getNumberOfNodes()<<" "
    //                                     <<f3->getMesh()->getNumberOfCells()<<endl;
    f3->setName("MyFieldOnGaussNE");
    f3->setDescription("MyDescriptionNE");
    DataArrayDouble *array=DataArrayDouble::New();
    //int nb=f1->getMesh()->getNumberOfNodes();

    /*8 pt de gauss by cell
      int nb=f3->getMesh()->getNumberOfCells()*8;
      array->alloc(nb,2);
      double *ptr=array->getPointer();
      for (int i=0; i<nb*2; i=i+2) {ptr[i]=(double)(i/8) ; ptr[i]=2.*(double)(i/8);}
    */

    //more nbptgauss=8 by default needs set MEDCouplingFieldDiscretizationPerCell
    //theory: (may be) http://www.code-aster.org/V2/doc/v9/fr/man_r/r3/r3.06.03.pdf
    int nbptgauss=8; //nb pt de gauss by cell
    int nbcell=f3->getMesh()->getNumberOfCells();
    int nb=nbcell*nbptgauss;
    int nbcomp=2;
    array->alloc(nb,nbcomp);
    double *ptr=array->getPointer();
    int ii=0;
    for (int i=0; i<nbcell; i++)
      for (int j=0; j<nbptgauss; j++)
        for (int k=0; k<nbcomp; k++)
          {
            //123.4 for 12th cell,3rd component, 4th gausspoint
            ptr[ii]=(double)((i+1)*10+(k+1))+((double)(j+1))/10.;
            ii++;
          }
    array->setInfoOnComponent(0,"vGx");
    array->setInfoOnComponent(1,"vGy");
    f3->setTime(4.,5,6);
    f3->setArray(array);
    array->decrRef();
    WriteField(name.c_str(),f3,true);
    if (_verbose) cout<<endl<<name<<" created"<<endl;
    f3->checkConsistencyLight();
    f1->decrRef();
    if (_ntot<1000000) //too long
      {
        MCAuto<MEDCouplingFieldDouble> f4=ReadField(ON_GAUSS_NE, name.c_str(), f3->getMesh()->getName().c_str(), 0, "MyFieldOnGaussNE", 5, 6);
        if (_verbose) cout<<"MyFieldOnGaussNE reread"<<endl;
      }
    f3->decrRef();
  }
  {
    string name=_file_name_with_faces;
    MEDCouplingFieldDouble *f1=buildVecFieldOnCells(name);
    name.replace(name.find(".med"),4,"_WithVecFieldOnCells.med");
    WriteField(name.c_str(),f1,true);
    if (_verbose) cout<<endl<<name<<" created"<<endl;
    if (_ntot<1000000) //too long
      {
        MCAuto<MEDCouplingFieldDouble> f2=ReadFieldCell(name.c_str(),f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),0,1);
        if (_verbose) cout<<name<<" reread"<<endl;
        //CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12)); assertion failed!!
      }
    f1->decrRef();
  }
}

void MEDPARTITIONERTest::createTestMeshWithVecFieldOnNodes()
{
  MEDCouplingFieldDouble *f1=buildVecFieldOnNodes();
  string name=_file_name;
  name.replace(name.find(".med"),4,"_WithVecFieldOnNodes.med");
  WriteField(name.c_str(),f1,true);
  if (_verbose) cout<<endl<<name<<" created"<<endl;
  if (_ntot<1000000) //too long
    {
      MCAuto<MEDCouplingFieldDouble> f2=ReadFieldNode(name.c_str(),f1->getMesh()->getName().c_str(),0,f1->getName().c_str(),0,1);
      if (_verbose) cout<<name<<" reread"<<endl;
      //CPPUNIT_ASSERT(f1->isEqual(f2,1e-12,1e-12)); assertion failed!!
    }
  f1->decrRef();
}

void MEDPARTITIONERTest::verifyTestMeshWithVecFieldOnNodes()
{
  string name=_file_name;
  name.replace(name.find(".med"),4,"_WithVecFieldOnNodes.med");
  MEDCouplingUMesh * m=ReadUMeshFromFile(name.c_str(),_mesh_name.c_str(),0);
  std::set<INTERP_KERNEL::NormalizedCellType> types(m->getAllGeoTypes());
  if (_verbose)
    {
      cout<<"\n types in "<<name<<" : ";
      //for (std::set<INTERP_KERNEL::NormalizedCellType>::iterator t=types.begin(); t!=types.end(); ++t) cout<<" "<<*t;
      for (std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator t=types.begin(); t!=types.end(); ++t)
        {
          //INTERP_KERNEL::CellModel essai=INTERP_KERNEL::CellModel::GetCellModel(*t);
          cout<<" "<<(INTERP_KERNEL::CellModel::GetCellModel(*t)).getRepr();
        }
      cout<<endl;
    }
  m->decrRef();

  MEDFileUMesh * mf = MEDFileUMesh::New(_file_name.c_str(),_mesh_name.c_str(),-1,-1);
  vector<int> lev;
  lev=mf->getNonEmptyLevels();
  if (_verbose)
    {
      cout<<" levels in "<<name<<" : ";
      for (vector<int>::iterator l=lev.begin(); l!=lev.end(); ++l) cout<<" "<<*l;
      cout<<endl;
    }
  mf->decrRef();
}

void MEDPARTITIONERTest::createTestMeshes()
{
  createTestMeshWithoutField();
  createTestMeshWithVecFieldOnCells();
  createTestMeshWithVecFieldOnNodes();
}

void MEDPARTITIONERTest::deleteTestMeshes()
{
  string cmd="rm *tmp_testMesh*";
  if (_verbose) cout<<endl<<cmd<<endl;
  system(cmd.c_str());  //may be not if debug
}

void MEDPARTITIONERTest::testMeshCollectionSingle()
{
  setSmallSize();
  createTestMeshes();
  MyGlobals::_World_Size=1;
  MyGlobals::_Rank=0;
  string fileName=_file_name_with_faces;
  MEDPARTITIONER::ParaDomainSelector parallelizer(false);
  MEDPARTITIONER::MeshCollection collection(fileName,parallelizer);
  CPPUNIT_ASSERT(collection.isParallelMode());
  CPPUNIT_ASSERT_EQUAL(3, collection.getMeshDimension());
  CPPUNIT_ASSERT(collection.getName()=="testMesh");
  CPPUNIT_ASSERT_EQUAL(1,collection.getNbOfLocalMeshes());
  CPPUNIT_ASSERT_EQUAL(1,collection.getNbOfGlobalMeshes());
  CPPUNIT_ASSERT_EQUAL(_ni*_nj*_nk,collection.getNbOfLocalCells());
  CPPUNIT_ASSERT_EQUAL(_ni*_nj,collection.getNbOfLocalFaces());
}

void MEDPARTITIONERTest::testMeshCollectionXml()
{
  setSmallSize();
  createHugeTestMesh(_ni, _nj, _nk, 2, 2, 2, 32); //xml but not so huge
  string fileName=_file_name_huge_xml;
  MEDPARTITIONER::ParaDomainSelector parallelizer(false);
  MEDPARTITIONER::MeshCollection collection(fileName,parallelizer);
  CPPUNIT_ASSERT(collection.isParallelMode());
  CPPUNIT_ASSERT_EQUAL(3, collection.getMeshDimension());
  CPPUNIT_ASSERT(collection.getName()=="testMesh");
  CPPUNIT_ASSERT_EQUAL(8,collection.getNbOfLocalMeshes());
  CPPUNIT_ASSERT_EQUAL(8,collection.getNbOfGlobalMeshes());
  CPPUNIT_ASSERT_EQUAL(_ni*_nj*_nk*8,collection.getNbOfLocalCells());
  CPPUNIT_ASSERT_EQUAL(0,collection.getNbOfLocalFaces());
}


//#################for metis



#if defined(MED_ENABLE_METIS)
void MEDPARTITIONERTest::testMeshCollectionSinglePartitionMetis()
{
  setSmallSize();
  createTestMeshes();
  //MyGlobals::_Verbose=500;
  string fileName=_file_name_with_faces;
  int ndomains=2;
  bool split_family=false;
  bool empty_groups=false;
  MEDPARTITIONER::ParaDomainSelector parallelizer(false);
  MEDPARTITIONER::MeshCollection collection(fileName,parallelizer);

  MEDPARTITIONER::ParallelTopology* aPT = (MEDPARTITIONER::ParallelTopology*) collection.getTopology();
  aPT->setGlobalNumerotationDefault(collection.getParaDomainSelector());
  //Creating the graph and partitioning it
  auto_ptr< MEDPARTITIONER::Topology > new_topo;
  new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::METIS) );
  //Creating a new mesh collection from the partitioning
  MEDPARTITIONER::MeshCollection new_collection(collection,new_topo.get(),split_family,empty_groups);

  //example to create files
  //MyGlobals::_General_Informations.clear();
  //MyGlobals::_General_Informations.push_back(SerializeFromString("finalMeshName=Merge"));
  //if (MyGlobals::_Verbose>100) cout << "generalInformations : \n"<<ReprVectorOfString(MyGlobals::_General_Informations);
  //new_collection.write("ttmp")

  CPPUNIT_ASSERT(new_collection.isParallelMode());
  CPPUNIT_ASSERT_EQUAL(3, new_collection.getMeshDimension());
  CPPUNIT_ASSERT(new_collection.getName()==collection.getName());
  CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfLocalMeshes());
  CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfGlobalMeshes());
  CPPUNIT_ASSERT_EQUAL(collection.getNbOfLocalCells(),new_collection.getNbOfLocalCells());
  CPPUNIT_ASSERT_EQUAL(collection.getNbOfLocalFaces(),new_collection.getNbOfLocalFaces());
}

void MEDPARTITIONERTest::testMeshCollectionComplexPartitionMetis()
{
  setSmallSize();
  createHugeTestMesh(_ni, _nj, _nk, 2, 2, 2, 32); //xml on 2*2*2 meshes but not so huge
  string fileName=_file_name_huge_xml;
  bool split_family=false;
  bool empty_groups=false;
  MEDPARTITIONER::ParaDomainSelector parallelizer(false);
  MEDPARTITIONER::MeshCollection collection(fileName,parallelizer);

  MEDPARTITIONER::ParallelTopology* aPT = (MEDPARTITIONER::ParallelTopology*) collection.getTopology();
  aPT->setGlobalNumerotationDefault(collection.getParaDomainSelector());

  for (int ndomains=2 ; ndomains<=16 ; ndomains++)
    {
      //Creating the graph and partitioning it
      auto_ptr< MEDPARTITIONER::Topology > new_topo;
      new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::METIS) );
      //Creating a new mesh collection from the partitioning
      MEDPARTITIONER::MeshCollection new_collection(collection,new_topo.get(),split_family,empty_groups);

      CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfLocalMeshes());
      CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfGlobalMeshes());
      CPPUNIT_ASSERT_EQUAL(collection.getNbOfLocalCells(),new_collection.getNbOfLocalCells());
      CPPUNIT_ASSERT_EQUAL(collection.getNbOfLocalFaces(),new_collection.getNbOfLocalFaces());
    }
}

void MEDPARTITIONERTest::testMetisSmallSize()
{
  //#if !defined(HAVE_MPI)
  setSmallSize();
  createTestMeshes();
  std::string MetisOrScotch("metis");
  launchMetisOrScotchMedpartitionerOnTestMeshes(MetisOrScotch);
  verifyMetisOrScotchMedpartitionerOnSmallSizeForMesh(MetisOrScotch);
  verifyMetisOrScotchMedpartitionerOnSmallSizeForFieldOnCells(MetisOrScotch);
  verifyMetisOrScotchMedpartitionerOnSmallSizeForFieldOnGaussNe(MetisOrScotch);
  //#endif
}
#endif


//#################for scotch


#if defined(MED_ENABLE_SCOTCH)
void MEDPARTITIONERTest::testMeshCollectionSinglePartitionScotch()
{
  setSmallSize();
  createTestMeshes();
  //MyGlobals::_Verbose=500;
  string fileName=_file_name_with_faces;
  int ndomains=2;
  bool split_family=false;
  bool empty_groups=false;
  MEDPARTITIONER::ParaDomainSelector parallelizer(false);
  MEDPARTITIONER::MeshCollection collection(fileName,parallelizer);

  MEDPARTITIONER::ParallelTopology* aPT = (MEDPARTITIONER::ParallelTopology*) collection.getTopology();
  aPT->setGlobalNumerotationDefault(collection.getParaDomainSelector());
  //Creating the graph and partitioning it
  auto_ptr< MEDPARTITIONER::Topology > new_topo;
  new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::SCOTCH) );
  //Creating a new mesh collection from the partitioning
  MEDPARTITIONER::MeshCollection new_collection(collection,new_topo.get(),split_family,empty_groups);

  //example to create files
  //MyGlobals::_General_Informations.clear();
  //MyGlobals::_General_Informations.push_back(SerializeFromString("finalMeshName=Merge"));
  //if (MyGlobals::_Verbose>100) cout << "generalInformations : \n"<<ReprVectorOfString(MyGlobals::_General_Informations);
  //new_collection.write("ttmp")

  CPPUNIT_ASSERT(new_collection.isParallelMode());
  CPPUNIT_ASSERT_EQUAL(3, new_collection.getMeshDimension());
  CPPUNIT_ASSERT(new_collection.getName()==collection.getName());
  CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfLocalMeshes());
  CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfGlobalMeshes());
  CPPUNIT_ASSERT_EQUAL(collection.getNbOfLocalCells(),new_collection.getNbOfLocalCells());
  CPPUNIT_ASSERT_EQUAL(collection.getNbOfLocalFaces(),new_collection.getNbOfLocalFaces());
}

void MEDPARTITIONERTest::testMeshCollectionComplexPartitionScotch()
{
  setSmallSize();
  createHugeTestMesh(_ni, _nj, _nk, 2, 2, 2, 32); //xml on 2*2*2 meshes but not so huge
  string fileName=_file_name_huge_xml;
  bool split_family=false;
  bool empty_groups=false;
  MEDPARTITIONER::ParaDomainSelector parallelizer(false);
  MEDPARTITIONER::MeshCollection collection(fileName,parallelizer);

  MEDPARTITIONER::ParallelTopology* aPT = (MEDPARTITIONER::ParallelTopology*) collection.getTopology();
  aPT->setGlobalNumerotationDefault(collection.getParaDomainSelector());

  for (int ndomains=2 ; ndomains<=16 ; ndomains++)
    {
      //Creating the graph and partitioning it
      auto_ptr< MEDPARTITIONER::Topology > new_topo;
      new_topo.reset( collection.createPartition(ndomains,MEDPARTITIONER::Graph::SCOTCH) );
      //Creating a new mesh collection from the partitioning
      MEDPARTITIONER::MeshCollection new_collection(collection,new_topo.get(),split_family,empty_groups);

      CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfLocalMeshes());
      CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfGlobalMeshes());
      CPPUNIT_ASSERT_EQUAL(collection.getNbOfLocalCells(),new_collection.getNbOfLocalCells());
      CPPUNIT_ASSERT_EQUAL(collection.getNbOfLocalFaces(),new_collection.getNbOfLocalFaces());
    }
}

void MEDPARTITIONERTest::testScotchSmallSize()
{
  //#if !defined(HAVE_MPI)
  setSmallSize();
  createTestMeshes();
  std::string MetisOrScotch("scotch");
  launchMetisOrScotchMedpartitionerOnTestMeshes(MetisOrScotch);
  verifyMetisOrScotchMedpartitionerOnSmallSizeForMesh(MetisOrScotch);
  verifyMetisOrScotchMedpartitionerOnSmallSizeForFieldOnCells(MetisOrScotch);
  verifyMetisOrScotchMedpartitionerOnSmallSizeForFieldOnGaussNe(MetisOrScotch);
  //#endif
}
#endif

void MEDPARTITIONERTest::launchMetisOrScotchMedpartitionerOnTestMeshes(std::string MetisOrScotch)
{
  int res;
  string cmd,execName,sourceName,targetName;

  execName=getPartitionerExe();

  cmd="which "+execName+" 2>/dev/null 1>/dev/null";  //no trace
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL_MESSAGE(execName + " - INVALID PATH TO medpartitioner", 0, res);

  cmd=execName+" --ndomains=2 --split-method="+MetisOrScotch;  //on same proc
  sourceName=_file_name;
  targetName=_file_name;
  targetName.replace(targetName.find(".med"),4,"_partitionedTo2_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);

  cmd=execName+" --ndomains=5 --split-method="+MetisOrScotch; //on less proc
  sourceName=_file_name;
  targetName=_file_name;
  targetName.replace(targetName.find(".med"),4,"_partitionedTo5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);

  cmd=execName+" --ndomains=1 --split-method="+MetisOrScotch;  //on 1 proc
  sourceName=targetName+".xml";
  targetName=_file_name;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);

  cmd=execName+" --ndomains=1 --split-method="+MetisOrScotch;  //on more proc
  //sourceName=targetName+".xml";
  targetName=_file_name;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
}

void MEDPARTITIONERTest::verifyMetisOrScotchMedpartitionerOnSmallSizeForMesh(std::string MetisOrScotch)
{
  int res;
  string fileName,cmd,execName,sourceName,targetName,input;
  execName=getPartitionerExe();
  fileName=_file_name_with_faces;

  MEDCoupling::MEDFileUMesh* initialMesh=MEDCoupling::MEDFileUMesh::New(fileName.c_str(),_mesh_name.c_str());
  MEDCoupling::MEDCouplingUMesh* cellMesh=initialMesh->getLevel0Mesh(false);
  MEDCoupling::MEDCouplingUMesh* faceMesh=initialMesh->getLevelM1Mesh(false);

  cmd=execName+" --ndomains=5 --split-method="+MetisOrScotch;  //on same proc
  sourceName=fileName;
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_partitionedTo5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  input=targetName+".xml";

  MEDPARTITIONER::ParaDomainSelector parallelizer(false);
  MEDPARTITIONER::MeshCollection collection(input,parallelizer);
  CPPUNIT_ASSERT_EQUAL(3, collection.getMeshDimension());
  std::vector<MEDCoupling::MEDCouplingUMesh*>cellMeshes=collection.getMesh();
  CPPUNIT_ASSERT_EQUAL(5, (int) cellMeshes.size());
  int nbcells=0;
  for (std::size_t i = 0; i < cellMeshes.size(); i++)
    nbcells+=cellMeshes[i]->getNumberOfCells();
  CPPUNIT_ASSERT_EQUAL((int)cellMesh->getNumberOfCells(), nbcells);

  std::vector<MEDCoupling::MEDCouplingUMesh*>faceMeshes=collection.getFaceMesh();
  CPPUNIT_ASSERT_EQUAL(5, (int) faceMeshes.size());
  int nbfaces=0;
  for (std::size_t i=0; i < faceMeshes.size(); i++)
    nbfaces+=faceMeshes[i]->getNumberOfCells();
  CPPUNIT_ASSERT_EQUAL((int)faceMesh->getNumberOfCells(), nbfaces);

  //merge split meshes and test equality
  cmd=execName+" --ndomains=1 --split-method="+MetisOrScotch;  //on same proc
  sourceName=targetName+".xml";
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);

  string refusedName=targetName+"1.med";
  MEDCoupling::MEDFileUMesh* refusedMesh=MEDCoupling::MEDFileUMesh::New(refusedName.c_str(),_mesh_name.c_str());
  MEDCoupling::MEDCouplingUMesh* refusedCellMesh=refusedMesh->getLevel0Mesh(false);
  MEDCoupling::MEDCouplingUMesh* refusedFaceMesh=refusedMesh->getLevelM1Mesh(false);

  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), refusedCellMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(faceMesh->getNumberOfCells(), refusedFaceMesh->getNumberOfCells());

  /*not the good job
    MEDCoupling::MEDCouplingMesh* mergeCell=cellMesh->mergeMyselfWith(refusedCellMesh);
    CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), mergeCell->getNumberOfCells());

    MEDCoupling::MEDCouplingMesh* mergeFace=faceMesh->mergeMyselfWith(refusedFaceMesh);
    CPPUNIT_ASSERT_EQUAL(faceMesh->getNumberOfCells(), mergeFace->getNumberOfCells());

    CPPUNIT_ASSERT(faceMesh->isEqual(refusedFaceMesh,1e-12));
  */

  std::vector<const MEDCouplingUMesh *> meshes;
  std::vector<DataArrayInt *> corr;
  meshes.push_back(cellMesh);
  refusedCellMesh->tryToShareSameCoordsPermute(*cellMesh, 1e-9);
  meshes.push_back(refusedCellMesh);
  MEDCouplingUMesh* fusedCell=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,0,corr);
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), fusedCell->getNumberOfCells());

  meshes.resize(0);
  for (std::size_t i = 0; i < corr.size(); i++)
    corr[i]->decrRef();
  corr.resize(0);
  meshes.push_back(faceMesh);
  refusedFaceMesh->tryToShareSameCoordsPermute(*faceMesh, 1e-9);
  meshes.push_back(refusedFaceMesh);
  MEDCouplingUMesh* fusedFace=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,0,corr);
  CPPUNIT_ASSERT_EQUAL(faceMesh->getNumberOfCells(), fusedFace->getNumberOfCells());

  for (std::size_t i = 0; i < corr.size(); i++)
    corr[i]->decrRef();
  fusedFace->decrRef();
  refusedFaceMesh->decrRef();
  faceMesh->decrRef();
  fusedCell->decrRef();
  refusedCellMesh->decrRef();
  refusedMesh->decrRef();
  cellMesh->decrRef();
  initialMesh->decrRef();
  //done in ~collection
  //for (int i = 0; i < faceMeshes.size(); i++) faceMeshes[i]->decrRef();
  //for (int i = 0; i < cellMeshes.size(); i++) cellMeshes[i]->decrRef();
}

void MEDPARTITIONERTest::verifyMetisOrScotchMedpartitionerOnSmallSizeForFieldOnCells(std::string MetisOrScotch)
{
  int res;
  string fileName,cmd,execName,sourceName,targetName,input;
  execName=getPartitionerExe();
  fileName=_file_name;
  fileName.replace(fileName.find(".med"),4,"_WithVecFieldOnCells.med");

  MEDCoupling::MEDFileUMesh* initialMesh=MEDCoupling::MEDFileUMesh::New(fileName.c_str(),_mesh_name.c_str());
  MEDCoupling::MEDCouplingUMesh* cellMesh=initialMesh->getLevel0Mesh(false);

  cmd=execName+" --ndomains=5 --split-method="+MetisOrScotch;  //on same proc
  sourceName=fileName;
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_partitionedTo5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  input=targetName+".xml";

  //merge split meshes and test equality
  cmd=execName+" --ndomains=1 --split-method="+MetisOrScotch;  //on same proc
  sourceName=targetName+".xml";
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);

  string refusedName=targetName+"1.med";
  MEDCoupling::MEDFileUMesh* refusedMesh=MEDCoupling::MEDFileUMesh::New(refusedName.c_str(),_mesh_name.c_str());
  MEDCoupling::MEDCouplingUMesh* refusedCellMesh=refusedMesh->getLevel0Mesh(false);

  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), refusedCellMesh->getNumberOfCells());

  std::vector<const MEDCouplingUMesh *> meshes;
  std::vector<DataArrayInt *> corr;
  meshes.push_back(cellMesh);
  refusedCellMesh->tryToShareSameCoordsPermute(*cellMesh, 1e-9);
  meshes.push_back(refusedCellMesh);
  MEDCouplingUMesh* fusedCell=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,0,corr);
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), fusedCell->getNumberOfCells());

  MCAuto<MEDCouplingFieldDouble> field1=ReadFieldCell(fileName.c_str(),initialMesh->getName().c_str(),0,"VectorFieldOnCells",0,1);
  MCAuto<MEDCouplingFieldDouble> field2=ReadFieldCell(refusedName.c_str(),refusedCellMesh->getName().c_str(),0,"VectorFieldOnCells",0,1);

  int nbcells=corr[1]->getNumberOfTuples();
  CPPUNIT_ASSERT_EQUAL((int)cellMesh->getNumberOfCells(), nbcells);
  //use corr to test equality of field
  DataArrayDouble* f1=field1->getArray();
  DataArrayDouble* f2=field2->getArray();
  if (_verbose>300)
    {
      cout<<"\nf1 : "<<f1->reprZip();
      cout<<"\nf2 : "<<f2->reprZip(); //field2->advancedRepradvancedRepr();
      for (std::size_t i = 0; i < corr.size(); i++)
        cout << "\ncorr " << i << " : " << corr[i]->reprZip();

    }
  int nbequal=0;
  int nbcomp=field1->getNumberOfComponents();
  double* p1=f1->getPointer();
  double* p2=f2->getPointer();
  int* pc=corr[1]->getPointer();
  for (int i = 0; i < nbcells; i++)
    {
      int i1=pc[i]*nbcomp;
      int i2=i*nbcomp;
      for (int j = 0; j < nbcomp; j++)
        {
          if (p1[i1+j]==p2[i2+j]) nbequal++;
          //cout<<" "<<p1[i1+j]<<"="<<p2[i2+j];
        }
    }
  CPPUNIT_ASSERT_EQUAL(nbcells*nbcomp, nbequal);

  for (std::size_t i = 0; i < corr.size(); i++)
    corr[i]->decrRef();
  fusedCell->decrRef();
  refusedMesh->decrRef();
  refusedCellMesh->decrRef();
  cellMesh->decrRef();
  initialMesh->decrRef();
}

void MEDPARTITIONERTest::verifyMetisOrScotchMedpartitionerOnSmallSizeForFieldOnGaussNe(std::string MetisOrScotch)
{
  int res;
  string fileName,cmd,execName,sourceName,targetName,input;
  execName=getPartitionerExe();
  fileName=_file_name;
  fileName.replace(fileName.find(".med"),4,"_WithVecFieldOnGaussNe.med");

  MEDCoupling::MEDFileUMesh* initialMesh=MEDCoupling::MEDFileUMesh::New(fileName.c_str(),_mesh_name.c_str());
  MEDCoupling::MEDCouplingUMesh* cellMesh=initialMesh->getLevel0Mesh(false);

  cmd=execName+" --ndomains=5 --split-method="+MetisOrScotch;  //on same proc
  sourceName=fileName;
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_partitionedTo5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  input=targetName+".xml";

  //merge split meshes and test equality
  cmd=execName+" --ndomains=1 --split-method="+MetisOrScotch;  //on same proc
  sourceName=targetName+".xml";
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);

  string refusedName=targetName+"1.med";
  MEDCoupling::MEDFileUMesh* refusedMesh=MEDCoupling::MEDFileUMesh::New(refusedName.c_str(),_mesh_name.c_str());
  MEDCoupling::MEDCouplingUMesh* refusedCellMesh=refusedMesh->getLevel0Mesh(false);

  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), refusedCellMesh->getNumberOfCells());

  std::vector<const MEDCouplingUMesh *> meshes;
  std::vector<DataArrayInt *> corr;
  meshes.push_back(cellMesh);
  refusedCellMesh->tryToShareSameCoordsPermute(*cellMesh, 1e-9);
  meshes.push_back(refusedCellMesh);
  MEDCouplingUMesh* fusedCell=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,0,corr);
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), fusedCell->getNumberOfCells());

  MCAuto<MEDCouplingFieldDouble> field1=ReadField(ON_GAUSS_NE,fileName.c_str(),initialMesh->getName().c_str(),0,"MyFieldOnGaussNE",5,6);
  MCAuto<MEDCouplingFieldDouble> field2=ReadField(ON_GAUSS_NE,refusedName.c_str(),refusedCellMesh->getName().c_str(),0,"MyFieldOnGaussNE",5,6);

  int nbcells=corr[1]->getNumberOfTuples();
  CPPUNIT_ASSERT_EQUAL((int)cellMesh->getNumberOfCells(), nbcells);
  //use corr to test equality of field
  DataArrayDouble* f1=field1->getArray();
  DataArrayDouble* f2=field2->getArray();
  if (_verbose>300)
    {
      cout << "\nf1 : " << f1->reprZip(); //123.4 for 12th cell,3rd component, 4th gausspoint
      cout << "\nf2 : " << f2->reprZip(); //field2->advancedRepradvancedRepr();
      for (std::size_t i = 0; i < corr.size(); i++)
        cout << "\ncorr " << i << " : " << corr[i]->reprZip();

    }
  int nbequal=0;
  int nbptgauss=8;
  int nbcomp=field1->getNumberOfComponents();
  double* p1=f1->getPointer();
  double* p2=f2->getPointer();
  int* pc=corr[1]->getPointer();
  for (int i = 0; i < nbcells; i++)
    {
      int i1=pc[i]*nbcomp*nbptgauss;
      int i2=i*nbcomp*nbptgauss;
      for (int j = 0; j < nbcomp*nbptgauss; j++)
        {
          if (p1[i1+j]==p2[i2+j]) nbequal++;
          //cout<<" "<<p1[i1+j]<<"="<<p2[i2+j];
        }
    }
  CPPUNIT_ASSERT_EQUAL(nbcells*nbcomp*nbptgauss, nbequal);

  for (std::size_t i = 0; i < corr.size(); i++)
    corr[i]->decrRef();
  fusedCell->decrRef();
  refusedMesh->decrRef();
  refusedCellMesh->decrRef();
  cellMesh->decrRef();
  initialMesh->decrRef();
}

//================================================================================
/*!
 * \brief Test for 0021756: [CEA 602] MEDPartitioner improvements
 */
//================================================================================

void MEDPARTITIONERTest::testCreateBoundaryFaces2D()
{
  // Fixed complains are:
  // - 2D is not available
  // - groups and family handling is bugged (probably due to bug in the handling
  //   of arrayTo in castIntField())
  // - creates boundary faces option is not handled

  // Create a 2D mesh in a file

  const char fileName[] = "tmp_testCreateBoundaryFaces2D.med";

  const int idFam1 = 3, idFam2 = 2;
  int nbFam1, nbFam2, nbc;
  {
    const int nbX = 20, nbY = 15;
    vector<int> conn;
    vector<double> coor;
    for (int j=0; j<=nbY; j++)
      for (int i=0; i<=nbX; i++)
        {
          coor.push_back(i+.1);
          coor.push_back(j+.2);
        }
    int ii;
    for (int j=0; j<nbY; j++)
      for (int i=0; i<nbX; i++)
        {
          ii=i + j*(nbX+1);
          conn.push_back(ii);
          conn.push_back(ii+1);
          ii=ii + nbX + 2 ;
          conn.push_back(ii);
          conn.push_back(ii-1);
        }
    MEDCouplingUMesh *mesh=MEDCouplingUMesh::New();
    mesh->setMeshDimension(2);

    nbc=conn.size()/4; //nb of cells
    mesh->allocateCells(nbc);
    int* pConn = &conn[0];
    for(int i=0; i<nbc; i++, pConn+=4)
      mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,pConn);
    mesh->finishInsertingCells();

    int nbv=coor.size()/2; //nb of vertices
    DataArrayDouble *myCoords=DataArrayDouble::New();
    myCoords->useArray( &coor[0], /*ownership=*/false, CPP_DEALLOC, nbv, 2 );
    mesh->setCoords(myCoords);
    mesh->setName("FacesIn2D");
    myCoords->decrRef();
    mesh->checkConsistencyLight();

    // groups of cells
    DataArrayInt* cellsFam=DataArrayInt::New();
    cellsFam->alloc(nbc,1);
    nbFam1 = nbc/3, nbFam2 = nbc/2;
    int iE = 0;
    for ( int i = 0; i < nbFam1; ++i ) cellsFam->getPointer()[ iE++ ] = idFam1;
    for ( int i = 0; i < nbFam2; ++i ) cellsFam->getPointer()[ iE++ ] = idFam2;
    for (            ; iE < nbc;     ) cellsFam->getPointer()[ iE++ ] = 0;
    map<string,int> theFamilies;
    theFamilies["FAMILLE_ZERO"]=0;
    theFamilies["Family1"     ]=idFam1;
    theFamilies["Family2"     ]=idFam2;
    map<string, vector<string> > theGroups;
    theGroups["Group1"].push_back("Family1");
    theGroups["Group2"].push_back("Family2");

    // write mesh
    MEDFileUMesh * fileMesh = MEDFileUMesh::New();
    fileMesh->setMeshAtLevel(0, mesh);
    fileMesh->setFamilyInfo(theFamilies);
    fileMesh->setGroupInfo(theGroups);
    fileMesh->setFamilyFieldArr(0, cellsFam);
    fileMesh->write(fileName,2);

    cellsFam->decrRef();
    mesh    ->decrRef();
    fileMesh->decrRef();

  } // mesh creation

  // Partition the mesh into 4 parts

  const int ndomains = 4;
  ParaDomainSelector parallelizer(false);
  MeshCollection collection(fileName,parallelizer);
  ParallelTopology* aPT = (ParallelTopology*) collection.getTopology();
  aPT->setGlobalNumerotationDefault(collection.getParaDomainSelector());

  std::auto_ptr< Topology > new_topo;
#if defined(MED_ENABLE_METIS) || defined(MED_ENABLE_PARMETIS)
  new_topo.reset( collection.createPartition(ndomains,Graph::METIS) );
#endif
#if defined(MED_ENABLE_SCOTCH)
  if ( !new_topo.get() )
    new_topo.reset( collection.createPartition(ndomains,Graph::SCOTCH) );
#endif
  if ( !new_topo.get() )
    return;

  // Check that "2D is available"

  const char xmlName[] = "tmp_testCreateBoundaryFaces2D";
  {
    MyGlobals::_Create_Boundary_Faces = true;
    MeshCollection new_collection(collection,new_topo.get());

    CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfLocalMeshes());
    CPPUNIT_ASSERT_EQUAL(ndomains,new_collection.getNbOfGlobalMeshes());
    CPPUNIT_ASSERT_EQUAL(collection.getNbOfLocalCells(),new_collection.getNbOfLocalCells());
    CPPUNIT_ASSERT_EQUAL(0,collection.getNbOfLocalFaces());
    CPPUNIT_ASSERT      (new_collection.getNbOfLocalFaces() > 0 );

    MyGlobals::_General_Informations.clear();
    MyGlobals::_General_Informations.push_back(SerializeFromString("finalMeshName=2D"));
    new_collection.write( xmlName );
  }

  // Check that "groups and family handling is NOT bugged"

  MeshCollection new_collection(std::string(xmlName)+".xml");
  std::map< int, int > famId2nb; // count total nb of cells in divided families
  std::map< int, int >::iterator id2nn;
  {
    const std::vector<MEDCoupling::DataArrayInt*>& famIdsVec = new_collection.getCellFamilyIds();
    for ( size_t i = 0; i < famIdsVec.size(); ++i )
      {
        MEDCoupling::DataArrayInt* famIdsArr = famIdsVec[i];
        for ( int j = famIdsArr->getNbOfElems()-1; j >= 0; --j )
          {
            id2nn = famId2nb.insert( make_pair( famIdsArr->getPointer()[j], 0 )).first;
            id2nn->second++;
          }
      }
  }
  CPPUNIT_ASSERT_EQUAL( 3, (int) famId2nb.size() ); // 3 fams/groups in all
  CPPUNIT_ASSERT_EQUAL( 1, (int) famId2nb.count( 0      ));
  CPPUNIT_ASSERT_EQUAL( 1, (int) famId2nb.count( idFam1 ));
  CPPUNIT_ASSERT_EQUAL( 1, (int) famId2nb.count( idFam2 ));
  CPPUNIT_ASSERT_EQUAL( nbFam1, famId2nb[ idFam1 ]);
  CPPUNIT_ASSERT_EQUAL( nbFam2, famId2nb[ idFam2 ]);
  CPPUNIT_ASSERT_EQUAL( nbc - nbFam1 - nbFam2, famId2nb[ 0 ]);

  // Check that "creates boundary faces option is handled"

  famId2nb.clear();
  const std::vector<MEDCoupling::DataArrayInt*>& famIdsVec = new_collection.getFaceFamilyIds();
  for ( size_t i = 0; i < famIdsVec.size(); ++i )
    {
      MEDCoupling::DataArrayInt* famIdsArr = famIdsVec[i];
      for ( int j = famIdsArr->getNbOfElems()-1; j >= 0; --j )
        {
          id2nn = famId2nb.insert( make_pair( famIdsArr->getPointer()[j], 0 )).first;
          id2nn->second++;
        }
    }

  CPPUNIT_ASSERT( !famId2nb.empty() );

  // for each "JOINT_n_p_..." group there must be "JOINT_p_n_..." group
  // of the same size
  std::map<std::string,int>& famName2id = new_collection.getFamilyInfo();
  std::map<std::string,int>::iterator na2id = famName2id.begin(), na2id2;
  std::set< int > okFamIds;
  okFamIds.insert(0);
  for ( ; na2id != famName2id.end(); ++na2id )
    {
      if ( okFamIds.count( na2id->second ) || na2id->first[0] != 'J')
        continue;
      na2id2 = na2id;
      bool groupOK = false;
      while ( !groupOK && ++na2id2 != famName2id.end() )
        groupOK = ( na2id2->first.find_first_not_of( na2id->first ) == std::string::npos );

      CPPUNIT_ASSERT( groupOK );
      CPPUNIT_ASSERT( na2id->second != na2id2->second);
      CPPUNIT_ASSERT_EQUAL( 1, (int) famId2nb.count( na2id2->second ));
      CPPUNIT_ASSERT_EQUAL( 1, (int) famId2nb.count( na2id->second ));
      CPPUNIT_ASSERT_EQUAL( (int) famId2nb[ na2id2->second ],
                            (int) famId2nb[ na2id->second ]);
      okFamIds.insert( na2id2->second );
    }
}
