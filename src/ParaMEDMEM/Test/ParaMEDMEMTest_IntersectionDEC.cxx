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
#include "ParaMEDMEMTest.hxx"
#include <cppunit/TestAssert.h>

#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "Topology.hxx"
#include "DEC.hxx"
#include "MxN_Mapping.hxx"
#include "IntersectionDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ComponentTopology.hxx"
#include "ICoCoMEDField.hxx"
#include "MEDLoader.hxx"
 
#include <string>

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES


using namespace std;
using namespace ParaMEDMEM;
 
void ParaMEDMEMTest::testIntersectionDEC_2D()
{
  testIntersectionDEC_2D_("P0","P0");
}

void ParaMEDMEMTest::testIntersectionDEC_3D()
{
  testIntersectionDEC_3D_("P0","P0");
}

void ParaMEDMEMTest::testIntersectionDEC_2DP0P1()
{
  //testIntersectionDEC_2D_("P0","P1");
}

/*
 * Check methods defined in IntersectionDEC.hxx
 *
 IntersectionDEC();
 IntersectionDEC(ProcessorGroup& local_group, ProcessorGroup& distant_group);
 virtual ~IntersectionDEC();
 void synchronize();
 void recvData();
 void sendData();
*/
 
void ParaMEDMEMTest::testIntersectionDEC_2D_(const char *srcMeth, const char *targetMeth)
{
  std::string srcM(srcMeth);
  std::string targetM(targetMeth);
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //the test is meant to run on five processors
  if (size !=5) return ;
   
  int nproc_source = 3;
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  
  ParaMEDMEM::CommInterface interface;
    
  ParaMEDMEM::ProcessorGroup* self_group = new ParaMEDMEM::MPIProcessorGroup(interface,self_procs);
  ParaMEDMEM::ProcessorGroup* target_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_target);
  ParaMEDMEM::ProcessorGroup* source_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_source);
  
  //loading the geometry for the source group

  ParaMEDMEM::IntersectionDEC dec (*source_group,*target_group);

  ParaMEDMEM::MEDCouplingUMesh* mesh;
  ParaMEDMEM::ParaMESH* paramesh;
  ParaMEDMEM::ParaFIELD* parafield;
  ICoCo::Field* icocofield ;
  
  string data_dir                   = getenv("MED_ROOT_DIR");
  string tmp_dir                    = getenv("TMP");
  if (tmp_dir == "")
    tmp_dir = "/tmp";
  string filename_xml1              = data_dir + "/share/salome/resources/med/square1_split";
  string filename_xml2              = data_dir + "/share/salome/resources/med/square2_split"; 
  string filename_seq_wr            = tmp_dir + "/";
  string filename_seq_med           = tmp_dir + "/myWrField_seq_pointe221.med";
  
  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (source_group->containsMyRank())
    {
      string master = filename_xml1;
      
      ostringstream strstream;
      strstream <<master<<rank+1<<".med";
      ostringstream meshname ;
      meshname<< "Mesh_2_"<< rank+1;
      
      mesh=MEDLoader::ReadUMeshFromFile(strstream.str().c_str(),meshname.str().c_str(),0);
      
    
      paramesh=new ParaMESH (mesh,*source_group,"source mesh");
    
      //      ParaMEDMEM::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT( support,*source_group);
      ParaMEDMEM::ComponentTopology comptopo;
      if(srcM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(ConservativeVolumic);
        }
      else
        parafield = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);
      int nb_local;
      if(srcM=="P0")
        nb_local=mesh->getNumberOfCells();
      else
        nb_local=mesh->getNumberOfNodes();
      //      double * value= new double[nb_local];
      double *value=parafield->getField()->getArray()->getPointer();
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=1.0;
    
      //      ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);
      icocofield=new ICoCo::MEDField(paramesh,parafield);
      dec.setMethod(srcMeth);
      dec.attachLocalField(icocofield);
    }
  
  //loading the geometry for the target group
  if (target_group->containsMyRank())
    {
      string master= filename_xml2;
      ostringstream strstream;
      strstream << master<<(rank-nproc_source+1)<<".med";
      ostringstream meshname ;
      meshname<< "Mesh_3_"<<rank-nproc_source+1;
      mesh = MEDLoader::ReadUMeshFromFile(strstream.str().c_str(),meshname.str().c_str(),0);
      
      paramesh=new ParaMESH (mesh,*target_group,"target mesh");
      //      ParaMEDMEM::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT(support,*target_group);
      ParaMEDMEM::ComponentTopology comptopo;
      if(targetM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(ConservativeVolumic);
        }
      else
        parafield = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);
      int nb_local;
      if(targetM=="P0")
        nb_local=mesh->getNumberOfCells();
      else
        nb_local=mesh->getNumberOfNodes();
      //      double * value= new double[nb_local];
      double *value=parafield->getField()->getArray()->getPointer();
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=0.0;
      //      ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);
      icocofield=new ICoCo::MEDField(paramesh,parafield);
      dec.setMethod(targetMeth);
      dec.attachLocalField(icocofield);
    }
    
  
  //attaching a DEC to the source group 
  double field_before_int;
  double field_after_int;
  
  if (source_group->containsMyRank())
    { 
      field_before_int = parafield->getVolumeIntegral(0);
      dec.synchronize();
      cout<<"DEC usage"<<endl;
      dec.setForcedRenormalization(false);

      dec.sendData();
      MEDLoader::writeParaMesh("./sourcesquareb",paramesh);
      if (source_group->myRank()==0)
        aRemover.Register("./sourcesquareb");
      ostringstream filename;
      filename<<"./sourcesquareb_"<<source_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      MEDLoader::writeParaField("./sourcesquareb","boundary",parafield);
   
      dec.recvData();
      cout <<"writing"<<endl;
      MEDLoader::writeParaMesh("./sourcesquare",paramesh);
      if (source_group->myRank()==0)
        aRemover.Register("./sourcesquare");
      MEDLoader::writeParaField("./sourcesquare","boundary",parafield);
      
     
      filename<<"./sourcesquare_"<<source_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      field_after_int = parafield->getVolumeIntegral(0);
      
      
      //      MPI_Bcast(&field_before_int,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      //       MPI_Bcast(&field_after_int,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(field_before_int, field_after_int, 1e-6);
    
    }
  
  //attaching a DEC to the target group
  if (target_group->containsMyRank())
    {
      dec.synchronize();
      dec.setForcedRenormalization(false);

      dec.recvData();
      MEDLoader::writeParaMesh("./targetsquareb",paramesh);
      MEDLoader::writeParaField("./targetsquareb", "boundary",parafield);
      if (target_group->myRank()==0)
        aRemover.Register("./targetsquareb");
      ostringstream filename;
      filename<<"./targetsquareb_"<<target_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      dec.sendData();
      MEDLoader::writeParaMesh("./targetsquare",paramesh);
      MEDLoader::writeParaField("./targetsquare", "boundary",parafield);
      
      if (target_group->myRank()==0)
        aRemover.Register("./targetsquareb");
      
      filename<<"./targetsquareb_"<<target_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      //    double field_before_int, field_after_int;
      //       MPI_Bcast(&field_before_int,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      //       MPI_Bcast(&field_after_int,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      
      //      CPPUNIT_ASSERT_DOUBLES_EQUAL(field_before_int, field_after_int, 1e-6);
    
    }
  
  delete source_group;
  delete target_group;
  delete self_group;
  delete parafield;
  delete paramesh;
  mesh->decrRef();

  delete icocofield;

  MPI_Barrier(MPI_COMM_WORLD);
  cout << "end of IntersectionDEC_2D test"<<endl;
}

void ParaMEDMEMTest::testIntersectionDEC_3D_(const char *srcMeth, const char *targetMeth)
{
  std::string srcM(srcMeth);
  std::string targetM(targetMeth);
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //the test is meant to run on five processors
  if (size !=3) return ;
   
  int nproc_source = 2;
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  
  ParaMEDMEM::CommInterface interface;
    
  ParaMEDMEM::ProcessorGroup* self_group = new ParaMEDMEM::MPIProcessorGroup(interface,self_procs);
  ParaMEDMEM::ProcessorGroup* target_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_target);
  ParaMEDMEM::ProcessorGroup* source_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_source);
  
  //loading the geometry for the source group

  ParaMEDMEM::IntersectionDEC dec (*source_group,*target_group);

  ParaMEDMEM::MEDCouplingUMesh* mesh;
  ParaMEDMEM::ParaMESH* paramesh;
  ParaMEDMEM::ParaFIELD* parafield;
  ICoCo::Field* icocofield ;
  
  string data_dir                   = getenv("MED_ROOT_DIR");
  string tmp_dir                    = getenv("TMP");
  if (tmp_dir == "")
    tmp_dir = "/tmp";
  string filename_xml1              = data_dir + "/share/salome/resources/med/Mesh3D_10_2d";
  string filename_xml2              = data_dir + "/share/salome/resources/med/Mesh3D_11"; 
  string filename_seq_wr            = tmp_dir + "/";
  string filename_seq_med           = tmp_dir + "/myWrField_seq_pointe221.med";
  
  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (source_group->containsMyRank())
    {
      string master = filename_xml1;
      
      ostringstream strstream;
      strstream <<master<<rank+1<<".med";
      ostringstream meshname ;
      meshname<< "Mesh_3_"<< rank+1;
      
      mesh=MEDLoader::ReadUMeshFromFile(strstream.str().c_str(),meshname.str().c_str(),0);
      
    
      paramesh=new ParaMESH (mesh,*source_group,"source mesh");
    
      //      ParaMEDMEM::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT( support,*source_group);
      ParaMEDMEM::ComponentTopology comptopo;
      if(srcM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(ConservativeVolumic);
        }
      else
        parafield = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);
      int nb_local;
      if(srcM=="P0")
        nb_local=mesh->getNumberOfCells();
      else
        nb_local=mesh->getNumberOfNodes();
      //      double * value= new double[nb_local];
      double *value=parafield->getField()->getArray()->getPointer();
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=1.0;
    
      //      ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);
      icocofield=new ICoCo::MEDField(paramesh,parafield);
      dec.setMethod(srcMeth);
      dec.attachLocalField(icocofield);
    }
  
  //loading the geometry for the target group
  if (target_group->containsMyRank())
    {
      string master= filename_xml2;
      ostringstream strstream;
      strstream << master << ".med";
      ostringstream meshname ;
      meshname<< "Mesh_6";
      mesh = MEDLoader::ReadUMeshFromFile(strstream.str().c_str(),meshname.str().c_str(),0);
      
      paramesh=new ParaMESH (mesh,*target_group,"target mesh");
      //      ParaMEDMEM::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT(support,*target_group);
      ParaMEDMEM::ComponentTopology comptopo;
      if(targetM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(ConservativeVolumic);
        }
      else
        parafield = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);
      int nb_local;
      if(targetM=="P0")
        nb_local=mesh->getNumberOfCells();
      else
        nb_local=mesh->getNumberOfNodes();
      //      double * value= new double[nb_local];
      double *value=parafield->getField()->getArray()->getPointer();
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=0.0;
      //      ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);
      icocofield=new ICoCo::MEDField(paramesh,parafield);
      dec.setMethod(targetMeth);
      dec.attachLocalField(icocofield);
    }  
  //attaching a DEC to the source group 
  double field_before_int;
  double field_after_int;
  
  if (source_group->containsMyRank())
    { 
      field_before_int = parafield->getVolumeIntegral(0);
      dec.synchronize();
      cout<<"DEC usage"<<endl;
      dec.setForcedRenormalization(false);

      dec.sendData();
      MEDLoader::writeParaMesh("./sourcesquareb",paramesh);
      if (source_group->myRank()==0)
        aRemover.Register("./sourcesquareb");
      ostringstream filename;
      filename<<"./sourcesquareb_"<<source_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      MEDLoader::writeParaField("./sourcesquareb","boundary",parafield);
   
      dec.recvData();
      cout <<"writing"<<endl;
      MEDLoader::writeParaMesh("./sourcesquare",paramesh);
      if (source_group->myRank()==0)
        aRemover.Register("./sourcesquare");
      MEDLoader::writeParaField("./sourcesquare","boundary",parafield);
      
     
      filename<<"./sourcesquare_"<<source_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      field_after_int = parafield->getVolumeIntegral(0);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(field_before_int, field_after_int, 1e-6);
    
    }
  
  //attaching a DEC to the target group
  if (target_group->containsMyRank())
    {
      dec.synchronize();
      dec.setForcedRenormalization(false);

      dec.recvData();
      MEDLoader::writeParaMesh("./targetsquareb",paramesh);
      MEDLoader::writeParaField("./targetsquareb", "boundary",parafield);
      if (target_group->myRank()==0)
        aRemover.Register("./targetsquareb");
      ostringstream filename;
      filename<<"./targetsquareb_"<<target_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      dec.sendData();
      MEDLoader::writeParaMesh("./targetsquare",paramesh);
      MEDLoader::writeParaField("./targetsquare", "boundary",parafield);
      
      if (target_group->myRank()==0)
        aRemover.Register("./targetsquareb");
      
      filename<<"./targetsquareb_"<<target_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
    }
  delete source_group;
  delete target_group;
  delete self_group;
  delete parafield;
  delete paramesh;
  mesh->decrRef();

  delete icocofield;

  MPI_Barrier(MPI_COMM_WORLD);
  cout << "end of IntersectionDEC_3D test"<<endl;
}

//Synchronous tests without interpolation with native mode (AllToAll(v) from lam/MPI:
void ParaMEDMEMTest::testSynchronousEqualIntersectionWithoutInterpNativeDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.1,1,0.1,1,false,false,false,"P0","P0");
}

//Synchronous tests without interpolation :
void ParaMEDMEMTest::testSynchronousEqualIntersectionWithoutInterpDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.1,1,0.1,1,true,false,false,"P0","P0");
}

//Synchronous tests with interpolation :
void ParaMEDMEMTest::testSynchronousEqualIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.1,1,0.1,1,true,false,true,"P0","P0");
}
void ParaMEDMEMTest::testSynchronousFasterSourceIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.09,1,0.1,1,true,false,true,"P0","P0");
}
void ParaMEDMEMTest::testSynchronousSlowerSourceIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.11,1,0.1,1,true,false,true,"P0","P0");
}
void ParaMEDMEMTest::testSynchronousSlowSourceIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.11,1,0.01,1,true,false,true,"P0","P0");
}
void ParaMEDMEMTest::testSynchronousFastSourceIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.01,1,0.11,1,true,false,true,"P0","P0");
}

//Asynchronous tests with interpolation :
void ParaMEDMEMTest::testAsynchronousEqualIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.1,1,0.1,1,true,true,true,"P0","P0");
}
void ParaMEDMEMTest::testAsynchronousFasterSourceIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.09,1,0.1,1,true,true,true,"P0","P0");
}
void ParaMEDMEMTest::testAsynchronousSlowerSourceIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.11,1,0.1,1,true,true,true,"P0","P0");
}
void ParaMEDMEMTest::testAsynchronousSlowSourceIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.11,1,0.01,1,true,true,true,"P0","P0");
}
void ParaMEDMEMTest::testAsynchronousFastSourceIntersectionDEC_2D()
{
  testAsynchronousIntersectionDEC_2D(0.01,1,0.11,1,true,true,true,"P0","P0");
}

void ParaMEDMEMTest::testIntersectionDECNonOverlapp_2D_P0P0()
{
  //
  const double sourceCoordsAll[2][8]={{0.4,0.5,0.4,1.5,1.6,1.5,1.6,0.5},
                                      {0.3,-0.5,1.6,-0.5,1.6,-1.5,0.3,-1.5}};
  const double targetCoordsAll[3][16]={{0.7,1.45,0.7,1.65,0.9,1.65,0.9,1.45,  1.1,1.4,1.1,1.6,1.3,1.6,1.3,1.4},
                                       {0.7,-0.6,0.7,0.7,0.9,0.7,0.9,-0.6,  1.1,-0.7,1.1,0.6,1.3,0.6,1.3,-0.7},
                                       {0.7,-1.55,0.7,-1.35,0.9,-1.35,0.9,-1.55,  1.1,-1.65,1.1,-1.45,1.3,-1.45,1.3,-1.65}};
  int conn4All[8]={0,1,2,3,4,5,6,7};
  double targetResults[3][2]={{34.,34.},{38.333333333333336,42.666666666666664},{47.,47.}};
  double targetResults2[3][2]={{0.28333333333333344,0.56666666666666687},{1.8564102564102569,2.0128205128205132},{1.0846153846153845,0.36153846153846159}};
  double targetResults3[3][2]={{3.7777777777777781,7.5555555555555562},{24.511111111111113,26.355555555555558},{14.1,4.7}};
  //
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=5)
    return ;
  int nproc_source = 2;
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  //
  ParaMEDMEM::MEDCouplingUMesh *mesh=0;
  ParaMEDMEM::ParaMESH *paramesh=0;
  ParaMEDMEM::ParaFIELD* parafield=0;
  ICoCo::Field* icocofield=0;
  //
  ParaMEDMEM::CommInterface interface;
  //
  ProcessorGroup* self_group = new ParaMEDMEM::MPIProcessorGroup(interface,self_procs);
  ProcessorGroup* target_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_source);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if(source_group->containsMyRank())
    {
      std::ostringstream stream; stream << "sourcemesh2D proc " << rank;
      mesh=MEDCouplingUMesh::New(stream.str().c_str(),2);
      mesh->allocateCells(2);
      mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn4All);
      mesh->finishInsertingCells();
      DataArrayDouble *myCoords=DataArrayDouble::New();
      myCoords->alloc(4,2);
      const double *sourceCoords=sourceCoordsAll[rank];
      std::copy(sourceCoords,sourceCoords+8,myCoords->getPointer());
      mesh->setCoords(myCoords);
      myCoords->decrRef();
      paramesh=new ParaMESH(mesh,*source_group,"source mesh");
      ParaMEDMEM::ComponentTopology comptopo;
      parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      double *value=parafield->getField()->getArray()->getPointer();
      value[0]=34+13*((double)rank);
    }
  else
    {
      std::ostringstream stream; stream << "targetmesh2D proc " << rank-nproc_source;
      mesh=MEDCouplingUMesh::New(stream.str().c_str(),2);
      mesh->allocateCells(2);
      mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn4All);
      mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn4All+4);
      mesh->finishInsertingCells();
      DataArrayDouble *myCoords=DataArrayDouble::New();
      myCoords->alloc(8,2);
      const double *targetCoords=targetCoordsAll[rank-nproc_source];
      std::copy(targetCoords,targetCoords+16,myCoords->getPointer());
      mesh->setCoords(myCoords);
      myCoords->decrRef();
      paramesh=new ParaMESH (mesh,*target_group,"target mesh");
      ParaMEDMEM::ComponentTopology comptopo;
      parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
    }
  //test 1 - Conservative volumic
  ParaMEDMEM::IntersectionDEC dec(*source_group,*target_group);
  parafield->getField()->setNature(ConservativeVolumic);
  if (source_group->containsMyRank())
    { 
      dec.setMethod("P0");
      dec.attachLocalField(parafield);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.sendData();
    }
  else
    {
      dec.setMethod("P0");
      dec.attachLocalField(parafield);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      const double *expected=targetResults[rank-nproc_source];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[0],res[0],1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[1],res[1],1e-13);
    }
  //test 2 - Integral
  ParaMEDMEM::IntersectionDEC dec2(*source_group,*target_group);
  parafield->getField()->setNature(Integral);
  if (source_group->containsMyRank())
    { 
      dec2.setMethod("P0");
      dec2.attachLocalField(parafield);
      dec2.synchronize();
      dec2.setForcedRenormalization(false);
      dec2.sendData();
    }
  else
    {
      dec2.setMethod("P0");
      dec2.attachLocalField(parafield);
      dec2.synchronize();
      dec2.setForcedRenormalization(false);
      dec2.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      const double *expected=targetResults2[rank-nproc_source];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[0],res[0],1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[1],res[1],1e-13);
    }
  //test 3 - Integral with global constraint
  ParaMEDMEM::IntersectionDEC dec3(*source_group,*target_group);
  parafield->getField()->setNature(IntegralGlobConstraint);
  if (source_group->containsMyRank())
    { 
      dec3.setMethod("P0");
      dec3.attachLocalField(parafield);
      dec3.synchronize();
      dec3.setForcedRenormalization(false);
      dec3.sendData();
    }
  else
    {
      dec3.setMethod("P0");
      dec3.attachLocalField(parafield);
      dec3.synchronize();
      dec3.setForcedRenormalization(false);
      dec3.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      const double *expected=targetResults3[rank-nproc_source];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[0],res[0],1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[1],res[1],1e-13);
    }
  //test 4 - Conservative volumic reversed
  ParaMEDMEM::IntersectionDEC dec4(*source_group,*target_group);
  parafield->getField()->setNature(ConservativeVolumic);
  if (source_group->containsMyRank())
    { 
      dec4.setMethod("P0");
      dec4.attachLocalField(parafield);
      dec4.synchronize();
      dec4.setForcedRenormalization(false);
      dec4.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_EQUAL(1,parafield->getField()->getNumberOfTuples());
      const double expected[]={37.8518518518519,43.5333333333333};
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[rank],res[0],1e-13);
    }
  else
    {
      dec4.setMethod("P0");
      dec4.attachLocalField(parafield);
      dec4.synchronize();
      dec4.setForcedRenormalization(false);
      double *res=parafield->getField()->getArray()->getPointer();
      const double *toSet=targetResults[rank-nproc_source];
      res[0]=toSet[0];
      res[1]=toSet[1];
      dec4.sendData();
    }
  //test 5 - Integral reversed
  ParaMEDMEM::IntersectionDEC dec5(*source_group,*target_group);
  parafield->getField()->setNature(Integral);
  if (source_group->containsMyRank())
    { 
      dec5.setMethod("P0");
      dec5.attachLocalField(parafield);
      dec5.synchronize();
      dec5.setForcedRenormalization(false);
      dec5.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_EQUAL(1,parafield->getField()->getNumberOfTuples());
      const double expected[]={0.794600591715977,1.35631163708087};
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[rank],res[0],1e-13);
    }
  else
    {
      dec5.setMethod("P0");
      dec5.attachLocalField(parafield);
      dec5.synchronize();
      dec5.setForcedRenormalization(false);
      double *res=parafield->getField()->getArray()->getPointer();
      const double *toSet=targetResults2[rank-nproc_source];
      res[0]=toSet[0];
      res[1]=toSet[1];
      dec5.sendData();
    }
  //test 6 - Integral with global constraint reversed
  ParaMEDMEM::IntersectionDEC dec6(*source_group,*target_group);
  parafield->getField()->setNature(IntegralGlobConstraint);
  if (source_group->containsMyRank())
    { 
      dec6.setMethod("P0");
      dec6.attachLocalField(parafield);
      dec6.synchronize();
      dec6.setForcedRenormalization(false);
      dec6.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_EQUAL(1,parafield->getField()->getNumberOfTuples());
      const double expected[]={36.4592592592593,44.5407407407407};
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[rank],res[0],1e-13);
    }
  else
    {
      dec6.setMethod("P0");
      dec6.attachLocalField(parafield);
      dec6.synchronize();
      dec6.setForcedRenormalization(false);
      double *res=parafield->getField()->getArray()->getPointer();
      const double *toSet=targetResults3[rank-nproc_source];
      res[0]=toSet[0];
      res[1]=toSet[1];
      dec6.sendData();
    }
  //
  delete parafield;
  mesh->decrRef();
  delete paramesh;
  delete self_group;
  delete target_group;
  delete source_group;
  //
  MPI_Barrier(MPI_COMM_WORLD);
}

void ParaMEDMEMTest::testIntersectionDECNonOverlapp_2D_P0P1P1P0()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=5)
    return ;
  int nproc_source = 2;
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  //
  ParaMEDMEM::MEDCouplingUMesh *mesh=0;
  ParaMEDMEM::ParaMESH *paramesh=0;
  ParaMEDMEM::ParaFIELD *parafieldP0=0,*parafieldP1=0;
  ICoCo::Field* icocofield=0;
  //
  ParaMEDMEM::CommInterface interface;
  //
  ProcessorGroup* self_group = new ParaMEDMEM::MPIProcessorGroup(interface,self_procs);
  ProcessorGroup* target_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_source);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if(source_group->containsMyRank())
    {
      if(rank==0)
        {
          double coords[6]={-0.3,-0.3, 0.7,0.7, 0.7,-0.3};
          int conn[3]={0,1,2};
          int globalNode[3]={1,2,0};
          mesh=MEDCouplingUMesh::New("Source mesh Proc0",2);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(3,2);
          std::copy(coords,coords+6,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      if(rank==1)
        {
          double coords[6]={-0.3,-0.3, -0.3,0.7, 0.7,0.7};
          int conn[3]={0,1,2};
          int globalNode[3]={1,3,2};
          mesh=MEDCouplingUMesh::New("Source mesh Proc1",2);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(3,2);
          std::copy(coords,coords+6,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      paramesh=new ParaMESH(mesh,*source_group,"source mesh");
      ParaMEDMEM::ComponentTopology comptopo;
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafieldP1 = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);
      double *valueP0=parafieldP0->getField()->getArray()->getPointer();
      double *valueP1=parafieldP1->getField()->getArray()->getPointer();
      parafieldP0->getField()->setNature(ConservativeVolumic);
      parafieldP1->getField()->setNature(ConservativeVolumic);
      if(rank==0)
        {
          valueP0[0]=31.;
          valueP1[0]=34.; valueP1[1]=77.; valueP1[2]=53.;
        }
      if(rank==1)
        {
          valueP0[0]=47.;
          valueP1[0]=34.; valueP1[1]=57.; valueP1[2]=77.;
        }
    }
  else
    {
      if(rank==2)
        {
          double coords[10]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2 };
          int conn[7]={0,3,4,1, 1,4,2};
          int globalNode[5]={4,3,0,2,1};
          mesh=MEDCouplingUMesh::New("Target mesh Proc2",2);
          mesh->allocateCells(2);
          mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
          mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+4);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(5,2);
          std::copy(coords,coords+10,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      if(rank==3)
        {
          double coords[6]={0.2,0.2, 0.7,-0.3, 0.7,0.2};
          int conn[3]={0,2,1};
          int globalNode[3]={1,0,5};
          mesh=MEDCouplingUMesh::New("Target mesh Proc3",2);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(3,2);
          std::copy(coords,coords+6,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      if(rank==4)
        {
          double coords[12]={-0.3,0.2, -0.3,0.7, 0.2,0.7, 0.2,0.2, 0.7,0.7, 0.7,0.2};
          int conn[8]={0,1,2,3, 3,2,4,5};
          int globalNode[6]={2,6,7,1,8,5};
          mesh=MEDCouplingUMesh::New("Target mesh Proc4",2);
          mesh->allocateCells(2);
          mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
          mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+4);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(6,2);
          std::copy(coords,coords+12,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      ParaMEDMEM::ComponentTopology comptopo;
      paramesh=new ParaMESH(mesh,*target_group,"target mesh");
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafieldP1 = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);
      parafieldP0->getField()->setNature(ConservativeVolumic);
      parafieldP1->getField()->setNature(ConservativeVolumic);
    }
  // test 1 - P0 P1
  ParaMEDMEM::IntersectionDEC dec(*source_group,*target_group);
  if (source_group->containsMyRank())
    { 
      dec.setMethod("P0");
      dec.attachLocalField(parafieldP0);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.sendData();
    }
  else
    {
      dec.setMethod("P1");
      dec.attachLocalField(parafieldP1);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.recvData();
      /*const double *res=parafield->getField()->getArray()->getConstPointer();
      const double *expected=targetResults[rank-nproc_source];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[0],res[0],1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[1],res[1],1e-13);*/
    }
  //
  delete parafieldP0;
  delete parafieldP1;
  mesh->decrRef();
  delete paramesh;
  delete self_group;
  delete target_group;
  delete source_group;
  //
  MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * Tests an asynchronous exchange between two codes
 * one sends data with dtA as an interval, the max time being tmaxA
 * the other one receives with dtB as an interval, the max time being tmaxB
 */
void ParaMEDMEMTest::testAsynchronousIntersectionDEC_2D(double dtA, double tmaxA, 
                                                        double dtB, double tmaxB, bool WithPointToPoint, bool Asynchronous,
                                                        bool WithInterp, const char *srcMeth, const char *targetMeth)
{
  std::string srcM(srcMeth);
  std::string targetM(targetMeth);
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 
  //the test is meant to run on five processors
  if (size !=5) return ;
   
  int nproc_source = 3;
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  
  ParaMEDMEM::CommInterface interface;
    
  ParaMEDMEM::ProcessorGroup* self_group = new ParaMEDMEM::MPIProcessorGroup(interface,self_procs);
  ParaMEDMEM::ProcessorGroup* target_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_target);
  ParaMEDMEM::ProcessorGroup* source_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_source);
    
  //loading the geometry for the source group

  ParaMEDMEM::IntersectionDEC dec (*source_group,*target_group);
  
  ParaMEDMEM::MEDCouplingUMesh* mesh;
  ParaMEDMEM::ParaMESH* paramesh;
  ParaMEDMEM::ParaFIELD* parafield;
  
  double * value ;
  ICoCo::Field* icocofield ;

  string data_dir                   = getenv("MED_ROOT_DIR");
  string tmp_dir                    = getenv("TMP");
  if (tmp_dir == "")
    tmp_dir = "/tmp";
  string filename_xml1              = data_dir + "/share/salome/resources/med/square1_split";
  string filename_xml2              = data_dir + "/share/salome/resources/med/square2_split"; 
  string filename_seq_wr            = tmp_dir + "/";
  string filename_seq_med           = tmp_dir + "/myWrField_seq_pointe221.med";
  
  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;
  
  MPI_Barrier(MPI_COMM_WORLD);

  if (source_group->containsMyRank())
    {
      string master = filename_xml1;
      
      ostringstream strstream;
      strstream <<master<<rank+1<<".med";
      ostringstream meshname ;
      meshname<< "Mesh_2_"<< rank+1;
      
      mesh=MEDLoader::ReadUMeshFromFile(strstream.str().c_str(),meshname.str().c_str(),0);

      paramesh=new ParaMESH (mesh,*source_group,"source mesh");
    
      //      ParaMEDMEM::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT( support,*source_group);
      ParaMEDMEM::ComponentTopology comptopo;
      if(srcM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(ConservativeVolumic);//InvertIntegral);//ConservativeVolumic);
        }
      else
        parafield = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);

      int nb_local;
      if(srcM=="P0")
        nb_local=mesh->getNumberOfCells();
      else
        nb_local=mesh->getNumberOfNodes();
      //      double * value= new double[nb_local];
      double *value=parafield->getField()->getArray()->getPointer();
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=0.0;
    
      //      ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);
      icocofield=new ICoCo::MEDField(paramesh,parafield);
     
      dec.attachLocalField(icocofield);


    }
  
  //loading the geometry for the target group
  if (target_group->containsMyRank())
    {
      string master= filename_xml2;
      ostringstream strstream;
      strstream << master<<(rank-nproc_source+1)<<".med";
      ostringstream meshname ;
      meshname<< "Mesh_3_"<<rank-nproc_source+1;
      
      mesh = MEDLoader::ReadUMeshFromFile(strstream.str().c_str(),meshname.str().c_str(),0);

      paramesh=new ParaMESH (mesh,*target_group,"target mesh");
      //      ParaMEDMEM::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT(support,*target_group);
      ParaMEDMEM::ComponentTopology comptopo;
      if(targetM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(ConservativeVolumic);//InvertIntegral);//ConservativeVolumic);
        }
      else
        parafield = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);
      
      int nb_local;
      if(targetM=="P0")
        nb_local=mesh->getNumberOfCells();
      else
        nb_local=mesh->getNumberOfNodes();
                        
      double *value=parafield->getField()->getArray()->getPointer();
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=0.0;
      //      ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);
      icocofield=new ICoCo::MEDField(paramesh,parafield);
      
      dec.attachLocalField(icocofield);
    }
    
  
  //attaching a DEC to the source group 
  
  if (source_group->containsMyRank())
    { 
      cout<<"DEC usage"<<endl;
      dec.setAsynchronous(Asynchronous);
      if ( WithInterp ) {
        dec.setTimeInterpolationMethod(LinearTimeInterp);
      }
      if ( WithPointToPoint ) {
        dec.setAllToAllMethod(PointToPoint);
      }
      else {
        dec.setAllToAllMethod(Native);
      }
      dec.synchronize();
      dec.setForcedRenormalization(false);
      for (double time=0; time<tmaxA+1e-10; time+=dtA)
        {
          cout << "testAsynchronousIntersectionDEC_2D" << rank << " time " << time
               << " dtA " << dtA << " tmaxA " << tmaxA << endl ;
          if ( time+dtA < tmaxA+1e-7 ) {
            dec.sendData( time , dtA );
          }
          else {
            dec.sendData( time , 0 );
          }
          double* value = parafield->getField()->getArray()->getPointer();
          int nb_local=parafield->getField()->getMesh()->getNumberOfCells();
          for (int i=0; i<nb_local;i++)
            value[i]= time+dtA;

       
        }
    }
  
  //attaching a DEC to the target group
  if (target_group->containsMyRank())
    {
      cout<<"DEC usage"<<endl;
      dec.setAsynchronous(Asynchronous);
      if ( WithInterp ) {
        dec.setTimeInterpolationMethod(LinearTimeInterp);
      }
      if ( WithPointToPoint ) {
        dec.setAllToAllMethod(PointToPoint);
      }
      else {
        dec.setAllToAllMethod(Native);
      }
      dec.synchronize();
      dec.setForcedRenormalization(false);
      vector<double> times;
      for (double time=0; time<tmaxB+1e-10; time+=dtB)
        {
          cout << "testAsynchronousIntersectionDEC_2D" << rank << " time " << time
               << " dtB " << dtB << " tmaxB " << tmaxB << endl ;
          dec.recvData( time );
          double vi = parafield->getVolumeIntegral(0);
          cout << "testAsynchronousIntersectionDEC_2D" << rank << " time " << time
               << " VolumeIntegral " << vi
               << " time*10000 " << time*10000 << endl ;
          
          CPPUNIT_ASSERT_DOUBLES_EQUAL(vi,time*10000,0.001);
        }
      
    }
  
  delete source_group;
  delete target_group;
  delete self_group;
  delete parafield ;
  delete paramesh ;
  mesh->decrRef() ;
  delete icocofield ;

  cout << "testAsynchronousIntersectionDEC_2D" << rank << " MPI_Barrier " << endl ;
 
  if (Asynchronous) MPI_Barrier(MPI_COMM_WORLD);
  cout << "end of IntersectionDEC_2D test"<<endl;
}
