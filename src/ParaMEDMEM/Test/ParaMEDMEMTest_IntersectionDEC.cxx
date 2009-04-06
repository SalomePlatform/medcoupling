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
        parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
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
        parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
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
        parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
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
        parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
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
        parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
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
        parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
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
