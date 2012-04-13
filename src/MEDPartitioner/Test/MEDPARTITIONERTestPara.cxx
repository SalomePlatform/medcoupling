// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingMultiFields.hxx"

#include <cppunit/TestAssert.h>

#include <sstream>
#include <cmath>
#include <list>
#include <stdexcept>
#include <cstdlib>
#include <vector>

#include <mpi.h>

using namespace std;
using namespace ParaMEDMEM;
using namespace MEDPARTITIONER;

#if defined(HAVE_MPI2)
void MEDPARTITIONERTest::verifyMedpartitionerOnSmallSizeForMesh()
{
  int res;
  string fileName,cmd,execName,sourceName,targetName,input;
  execName=getenv("MED_ROOT_DIR");  //.../INSTALL/MED
  execName+="/bin/salome/medpartitioner_para";
  fileName=_file_name_with_faces;
  
  ParaMEDMEM::MEDFileUMesh* initialMesh=ParaMEDMEM::MEDFileUMesh::New(fileName.c_str(),_mesh_name.c_str());
  ParaMEDMEM::MEDCouplingUMesh* cellMesh=initialMesh->getLevel0Mesh(false);
  ParaMEDMEM::MEDCouplingUMesh* faceMesh=initialMesh->getLevelM1Mesh(false);
  
  cmd="mpirun -np 5 "+execName+" --ndomains=5 --split-method=metis";  //on same proc
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
  std::vector<ParaMEDMEM::MEDCouplingUMesh*>cellMeshes=collection.getMesh();
  CPPUNIT_ASSERT_EQUAL(5, (int) cellMeshes.size());
  int nbcells=0;
  for (std::size_t i = 0; i < cellMeshes.size(); i++)
    nbcells+=cellMeshes[i]->getNumberOfCells();
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), nbcells);
  
  std::vector<ParaMEDMEM::MEDCouplingUMesh*>faceMeshes=collection.getFaceMesh();
  CPPUNIT_ASSERT_EQUAL(5, (int) faceMeshes.size());
  int nbfaces=0;
  for (std::size_t i=0; i < faceMeshes.size(); i++)
    nbfaces+=faceMeshes[i]->getNumberOfCells();
  CPPUNIT_ASSERT_EQUAL(faceMesh->getNumberOfCells(), nbfaces);
  
  //merge split meshes and test equality
  cmd="mpirun -np 1 "+execName+" --ndomains=1 --split-method=metis";  //on same proc
  sourceName=targetName+".xml";
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  
  string refusedName=targetName+"1.med";
  ParaMEDMEM::MEDFileUMesh* refusedMesh=ParaMEDMEM::MEDFileUMesh::New(refusedName.c_str(),_mesh_name.c_str());
  ParaMEDMEM::MEDCouplingUMesh* refusedCellMesh=refusedMesh->getLevel0Mesh(false);
  ParaMEDMEM::MEDCouplingUMesh* refusedFaceMesh=refusedMesh->getLevelM1Mesh(false);
  
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), refusedCellMesh->getNumberOfCells());
  CPPUNIT_ASSERT_EQUAL(faceMesh->getNumberOfCells(), refusedFaceMesh->getNumberOfCells());
  
  /*not the good job
    ParaMEDMEM::MEDCouplingMesh* mergeCell=cellMesh->mergeMyselfWith(refusedCellMesh);
    CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), mergeCell->getNumberOfCells());
  
    ParaMEDMEM::MEDCouplingMesh* mergeFace=faceMesh->mergeMyselfWith(refusedFaceMesh);
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
  cellMesh->decrRef();
  //done in ~collection
  //for (int i = 0; i < faceMeshes.size(); i++) faceMeshes[i]->decrRef();
  //for (int i = 0; i < cellMeshes.size(); i++) cellMeshes[i]->decrRef();
}

void MEDPARTITIONERTest::verifyMedpartitionerOnSmallSizeForFieldOnCells()
{
  int res;
  string fileName,cmd,execName,sourceName,targetName,input;
  execName=getenv("MED_ROOT_DIR");  //.../INSTALL/MED
  execName+="/bin/salome/medpartitioner_para";
  fileName=_file_name;
  fileName.replace(fileName.find(".med"),4,"_WithVecFieldOnCells.med");
  
  ParaMEDMEM::MEDFileUMesh* initialMesh=ParaMEDMEM::MEDFileUMesh::New(fileName.c_str(),_mesh_name.c_str());
  ParaMEDMEM::MEDCouplingUMesh* cellMesh=initialMesh->getLevel0Mesh(false);
  
  cmd="mpirun -np 5 "+execName+" --ndomains=5 --split-method=metis";  //on same proc
  sourceName=fileName;
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_partitionedTo5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  input=targetName+".xml";
  
  //merge split meshes and test equality
  cmd="mpirun -np 1 "+execName+" --ndomains=1 --split-method=metis";  //on same proc
  sourceName=targetName+".xml";
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  
  string refusedName=targetName+"1.med";
  ParaMEDMEM::MEDFileUMesh* refusedMesh=ParaMEDMEM::MEDFileUMesh::New(refusedName.c_str(),_mesh_name.c_str());
  ParaMEDMEM::MEDCouplingUMesh* refusedCellMesh=refusedMesh->getLevel0Mesh(false);
  
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), refusedCellMesh->getNumberOfCells());
  
  std::vector<const MEDCouplingUMesh *> meshes;
  std::vector<DataArrayInt *> corr;
  meshes.push_back(cellMesh);
  refusedCellMesh->tryToShareSameCoordsPermute(*cellMesh, 1e-9);
  meshes.push_back(refusedCellMesh);
  MEDCouplingUMesh* fusedCell=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,0,corr);
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), fusedCell->getNumberOfCells());
  
  MEDCouplingFieldDouble* field1=MEDLoader::ReadFieldCell(fileName.c_str(),initialMesh->getName(),0,"VectorFieldOnCells",0,1);
  MEDCouplingFieldDouble* field2=MEDLoader::ReadFieldCell(refusedName.c_str(),refusedCellMesh->getName(),0,"VectorFieldOnCells",0,1);
  
  int nbcells=corr[1]->getNumberOfTuples();
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), nbcells);
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
  field1->decrRef();
  field2->decrRef();
  fusedCell->decrRef();
  refusedCellMesh->decrRef();
  cellMesh->decrRef();
}

void MEDPARTITIONERTest::verifyMedpartitionerOnSmallSizeForFieldOnGaussNe()
{
  int res;
  string fileName,cmd,execName,sourceName,targetName,input;
  execName=getenv("MED_ROOT_DIR");  //.../INSTALL/MED
  execName+="/bin/salome/medpartitioner_para";
  fileName=_file_name;
  fileName.replace(fileName.find(".med"),4,"_WithVecFieldOnGaussNe.med");
  
  ParaMEDMEM::MEDFileUMesh* initialMesh=ParaMEDMEM::MEDFileUMesh::New(fileName.c_str(),_mesh_name.c_str());
  ParaMEDMEM::MEDCouplingUMesh* cellMesh=initialMesh->getLevel0Mesh(false);
  
  cmd="mpirun -np 5 "+execName+" --ndomains=5 --split-method=metis";  //on same proc
  sourceName=fileName;
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_partitionedTo5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  input=targetName+".xml";
  
  //merge split meshes and test equality
  cmd="mpirun -np 1 "+execName+" --ndomains=1 --split-method=metis";  //on same proc
  sourceName=targetName+".xml";
  targetName=fileName;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  
  string refusedName=targetName+"1.med";
  ParaMEDMEM::MEDFileUMesh* refusedMesh=ParaMEDMEM::MEDFileUMesh::New(refusedName.c_str(),_mesh_name.c_str());
  ParaMEDMEM::MEDCouplingUMesh* refusedCellMesh=refusedMesh->getLevel0Mesh(false);
  
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), refusedCellMesh->getNumberOfCells());
  
  std::vector<const MEDCouplingUMesh *> meshes;
  std::vector<DataArrayInt *> corr;
  meshes.push_back(cellMesh);
  refusedCellMesh->tryToShareSameCoordsPermute(*cellMesh, 1e-9);
  meshes.push_back(refusedCellMesh);
  MEDCouplingUMesh* fusedCell=MEDCouplingUMesh::FuseUMeshesOnSameCoords(meshes,0,corr);
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), fusedCell->getNumberOfCells());
  
  MEDCouplingFieldDouble* field1=MEDLoader::ReadField(ON_GAUSS_NE,fileName.c_str(),initialMesh->getName(),0,"MyFieldOnGaussNE",5,6);
  MEDCouplingFieldDouble* field2=MEDLoader::ReadField(ON_GAUSS_NE,refusedName.c_str(),refusedCellMesh->getName(),0,"MyFieldOnGaussNE",5,6);
  
  int nbcells=corr[1]->getNumberOfTuples();
  CPPUNIT_ASSERT_EQUAL(cellMesh->getNumberOfCells(), nbcells);
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
  field1->decrRef();
  field2->decrRef();
  fusedCell->decrRef();
  refusedCellMesh->decrRef();
  cellMesh->decrRef();
}

void MEDPARTITIONERTest::launchMedpartitionerOnTestMeshes()
{
  
  /* examples 
     export INFI=/home/vb144235/resources/blade.med
     //no need export MESH=Fuse_1
     export INFI=tmp_testMeshxxx.med
     //no need export MESH=testMesh
     mpirun -np 2 medpartitioner_para --input-file=$INFI --output-file=ttmp1_ --ndomains=4
     mpirun -np 5 medpartitioner_para --input-file=ttmp1_.xml --output-file=ttmp2_ --ndomains=5
     mpirun -np 2 valgrind  medpartitioner_para --input-file=tmp_testMesh_20x30x50.med  --output-file=ttmp1petit_ --ndomains=4  --dump-cpu-memory --verbose=111
  */
  int res;
  string cmd,execName,sourceName,targetName;
  
  res=system("which mpirun 2>/dev/null 1>/dev/null"); //no trace
  CPPUNIT_ASSERT_EQUAL(0, res);
  
  execName=getenv("MED_ROOT_DIR");  //.../INSTALL/MED
  execName+="/bin/salome/medpartitioner_para";
  
  cmd="which "+execName+" 2>/dev/null 1>/dev/null";  //no trace
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  
  cmd="mpirun -np 2 "+execName+" --ndomains=2 --split-method=metis";  //on same proc
  sourceName=_file_name;
  targetName=_file_name;
  targetName.replace(targetName.find(".med"),4,"_partitionedTo2_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  
  cmd="mpirun -np 3 "+execName+" --ndomains=5 --split-method=metis"; //on less proc
  sourceName=_file_name;
  targetName=_file_name;
  targetName.replace(targetName.find(".med"),4,"_partitionedTo5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
  
  cmd="mpirun -np 1 "+execName+" --ndomains=1 --split-method=metis";  //on 1 proc
  sourceName=targetName+".xml";
  targetName=_file_name;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);

  cmd="mpirun -np 8 "+execName+" --ndomains=1 --split-method=metis";  //on more proc
  //sourceName=targetName+".xml";
  targetName=_file_name;
  targetName.replace(targetName.find(".med"),4,"_remergedFrom5_");
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
}  

void MEDPARTITIONERTest::launchMedpartitionerOnHugeTestMeshes()
{
  int res=0;
  string cmd,execName,sourceName,targetName;
  execName=getenv("MED_ROOT_DIR");  //.../INSTALL/MED
  execName+="/bin/salome/medpartitioner_para";

  string snbTarget=IntToStr(_nb_target_huge);
  cmd="mpirun -np "+snbTarget+" "+execName+" --ndomains="+snbTarget+" --split-method=metis";  //on same proc
  sourceName=_file_name_huge_xml;
  targetName=_file_name_huge_xml;
  string tmp="_partitionedTo"+snbTarget+"_";
  targetName.replace(targetName.find(".xml"),4,tmp);
  cmd+=" --input-file="+sourceName+" --output-file="+targetName+" --verbose="+IntToStr(_verbose);
  if (_verbose) cout<<endl<<cmd<<endl;
  res=system(cmd.c_str());
  CPPUNIT_ASSERT_EQUAL(0, res);
}  

void MEDPARTITIONERTest::testMpirunSmallSize()
{
  setSmallSize();
  createTestMeshes();
  launchMedpartitionerOnTestMeshes();
  verifyMedpartitionerOnSmallSizeForMesh();
  verifyMedpartitionerOnSmallSizeForFieldOnCells();
  verifyMedpartitionerOnSmallSizeForFieldOnGaussNe();
}

void MEDPARTITIONERTest::testMpirunMedianSize()
{
  setMedianSize();
  createTestMeshes();
  launchMedpartitionerOnTestMeshes();
}

void MEDPARTITIONERTest::testMpirunHugeSize()
{
  //setBigSize(); //may be a lot for now
  setMedianSize();
  //create a set of nbx*nby*nbz files mesh of ni*ny*nz cells
  //_verbose=1;
  createHugeTestMesh(_ni, _nj, _nk, 2, 2, 2, 32); //it is now to know how far we are going to test
  launchMedpartitionerOnHugeTestMeshes();
}
#endif
