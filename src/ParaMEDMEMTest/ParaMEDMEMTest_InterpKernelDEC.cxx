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

#include "ParaMEDMEMTest.hxx"
#include <cppunit/TestAssert.h>

#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "Topology.hxx"
#include "DEC.hxx"
#include "MxN_Mapping.hxx"
#include "InterpKernelDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ComponentTopology.hxx"
#include "ICoCoMEDField.hxx"
#include "ParaMEDLoader.hxx"
#include "MEDLoader.hxx"
#include "TestInterpKernelUtils.hxx"

 
#include <string>
#include <iterator>

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES


using namespace std;
using namespace MEDCoupling;

void ParaMEDMEMTest::testInterpKernelDEC_2D()
{
  testInterpKernelDEC_2D_("P0","P0");
}

void ParaMEDMEMTest::testInterpKernelDEC2_2D()
{
  testInterpKernelDEC2_2D_("P0","P0");
}

void ParaMEDMEMTest::testInterpKernelDEC_3D()
{
  testInterpKernelDEC_3D_("P0","P0");
}

void ParaMEDMEMTest::testInterpKernelDEC_2DP0P1()
{
  //testInterpKernelDEC_2D_("P0","P1");
}

void ParaMEDMEMTest::testInterpKernelDEC_1D()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=5)
    return ;
  int nproc_source = 3;
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  //
  MEDCoupling::MEDCouplingUMesh *mesh=0;
  MEDCoupling::ParaMESH *paramesh=0;
  MEDCoupling::ParaFIELD *parafieldP0=0;
  //
  MEDCoupling::CommInterface interface;
  //
  ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if(source_group->containsMyRank())
    {
      if(rank==0)
        {
          double coords[4]={0.3,0.7, 0.9,1.0};
          int conn[4]={0,1,2,3};
          mesh=MEDCouplingUMesh::New("Source mesh Proc0",1);
          mesh->allocateCells(2);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(4,1);
          std::copy(coords,coords+4,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      if(rank==1)
        {
          double coords[2]={0.7,0.9};
          int conn[2]={0,1};
          mesh=MEDCouplingUMesh::New("Source mesh Proc1",1);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(2,1);
          std::copy(coords,coords+2,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      if(rank==2)
        {
          double coords[2]={1.,1.12};
          int conn[2]={0,1};
          mesh=MEDCouplingUMesh::New("Source mesh Proc2",1);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(2,1);
          std::copy(coords,coords+2,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      paramesh=new ParaMESH(mesh,*source_group,"source mesh");
      MEDCoupling::ComponentTopology comptopo;
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      double *valueP0=parafieldP0->getField()->getArray()->getPointer();
      parafieldP0->getField()->setNature(IntensiveMaximum);
      if(rank==0)
        {
          valueP0[0]=7.; valueP0[1]=8.;
        }
      if(rank==1)
        {
          valueP0[0]=9.;
        }
      if(rank==2)
        {
          valueP0[0]=10.;
        }
    }
  else
    {
      const char targetMeshName[]="target mesh";
      if(rank==3)
        {
          double coords[2]={0.5,0.75};
          int conn[2]={0,1};
          mesh=MEDCouplingUMesh::New("Target mesh Proc3",1);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(2,1);
          std::copy(coords,coords+2,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          paramesh=new ParaMESH(mesh,*target_group,targetMeshName);
        }
      if(rank==4)
        {
          double coords[2]={0.75,1.2};
          int conn[2]={0,1};
          mesh=MEDCouplingUMesh::New("Target mesh Proc4",1);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(2,1);
          std::copy(coords,coords+2,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          paramesh=new ParaMESH(mesh,*target_group,targetMeshName);
        }
      MEDCoupling::ComponentTopology comptopo;
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafieldP0->getField()->setNature(IntensiveMaximum);
    }
  // test 1
  MEDCoupling::InterpKernelDEC dec(*source_group,*target_group);
  if (source_group->containsMyRank())
    { 
      dec.setMethod("P0");
      dec.attachLocalField(parafieldP0);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.sendData();
      dec.recvData();
      const double *valueP0=parafieldP0->getField()->getArray()->getPointer();
      if(rank==0)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(7.4,valueP0[0],1e-7);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0540540540540544,valueP0[1],1e-7);
        }
      if(rank==1)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(8.64054054054054,valueP0[0],1e-7);
        }
      if(rank==2)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0540540540540544,valueP0[0],1e-7);
        }
    }
  else
    {
      dec.setMethod("P0");
      dec.attachLocalField(parafieldP0);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.recvData();
      const double *res=parafieldP0->getField()->getArray()->getConstPointer();
      if(rank==3)
        {
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP0->getField()->getNumberOfTuples());
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP0->getField()->getNumberOfComponents());
          CPPUNIT_ASSERT_DOUBLES_EQUAL(7.4,res[0],1e-12);
        }
      if(rank==4)
        {
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP0->getField()->getNumberOfTuples());
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP0->getField()->getNumberOfComponents());
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0540540540540526,res[0],1e-12);
        }
      dec.sendData();
    }
  //
  delete parafieldP0;
  mesh->decrRef();
  delete paramesh;
  delete self_group;
  delete target_group;
  delete source_group;
  //
  MPI_Barrier(MPI_COMM_WORLD);
}

void ParaMEDMEMTest::testInterpKernelDEC_2DCurve()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=5)
    return ;
  int nproc_source = 3;
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  //
  MEDCoupling::MEDCouplingUMesh *mesh=0;
  MEDCoupling::ParaMESH *paramesh=0;
  MEDCoupling::ParaFIELD *parafieldP0=0;
  //
  MEDCoupling::CommInterface interface;
  //
  ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if(source_group->containsMyRank())
    {
      if(rank==0)
        {
          double coords[8]={0.3,0.3,0.7,0.7, 0.9,0.9,1.0,1.0};
          int conn[4]={0,1,2,3};
          mesh=MEDCouplingUMesh::New("Source mesh Proc0",1);
          mesh->allocateCells(2);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn+2);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(4,2);
          std::copy(coords,coords+8,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      if(rank==1)
        {
          double coords[4]={0.7,0.7,0.9,0.9};
          int conn[2]={0,1};
          mesh=MEDCouplingUMesh::New("Source mesh Proc1",1);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(2,2);
          std::copy(coords,coords+4,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      if(rank==2)
        {
          double coords[4]={1.,1.,1.12,1.12};
          int conn[2]={0,1};
          mesh=MEDCouplingUMesh::New("Source mesh Proc2",1);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(2,2);
          std::copy(coords,coords+4,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
        }
      paramesh=new ParaMESH(mesh,*source_group,"source mesh");
      MEDCoupling::ComponentTopology comptopo;
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      double *valueP0=parafieldP0->getField()->getArray()->getPointer();
      parafieldP0->getField()->setNature(IntensiveMaximum);
      if(rank==0)
        {
          valueP0[0]=7.; valueP0[1]=8.;
        }
      if(rank==1)
        {
          valueP0[0]=9.;
        }
      if(rank==2)
        {
          valueP0[0]=10.;
        }
    }
  else
    {
      const char targetMeshName[]="target mesh";
      if(rank==3)
        {
          double coords[4]={0.5,0.5,0.75,0.75};
          int conn[2]={0,1};
          mesh=MEDCouplingUMesh::New("Target mesh Proc3",1);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(2,2);
          std::copy(coords,coords+4,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          paramesh=new ParaMESH(mesh,*target_group,targetMeshName);
        }
      if(rank==4)
        {
          double coords[4]={0.75,0.75,1.2,1.2};
          int conn[2]={0,1};
          mesh=MEDCouplingUMesh::New("Target mesh Proc4",1);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(2,2);
          std::copy(coords,coords+4,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          paramesh=new ParaMESH(mesh,*target_group,targetMeshName);
        }
      MEDCoupling::ComponentTopology comptopo;
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafieldP0->getField()->setNature(IntensiveMaximum);
    }
  // test 1
  MEDCoupling::InterpKernelDEC dec(*source_group,*target_group);
  if (source_group->containsMyRank())
    { 
      dec.setMethod("P0");
      dec.attachLocalField(parafieldP0);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.sendData();
      dec.recvData();
      const double *valueP0=parafieldP0->getField()->getArray()->getPointer();
      if(rank==0)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(7.4,valueP0[0],1e-7);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0540540540540544,valueP0[1],1e-7);
        }
      if(rank==1)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(8.64054054054054,valueP0[0],1e-7);
        }
      if(rank==2)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0540540540540544,valueP0[0],1e-7);
        }
    }
  else
    {
      dec.setMethod("P0");
      dec.attachLocalField(parafieldP0);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.recvData();
      const double *res=parafieldP0->getField()->getArray()->getConstPointer();
      if(rank==3)
        {
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP0->getField()->getNumberOfTuples());
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP0->getField()->getNumberOfComponents());
          CPPUNIT_ASSERT_DOUBLES_EQUAL(7.4,res[0],1e-12);
        }
      if(rank==4)
        {
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP0->getField()->getNumberOfTuples());
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP0->getField()->getNumberOfComponents());
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0540540540540526,res[0],1e-12);
        }
      dec.sendData();
    }
  //
  delete parafieldP0;
  mesh->decrRef();
  delete paramesh;
  delete self_group;
  delete target_group;
  delete source_group;
  //
  MPI_Barrier(MPI_COMM_WORLD);
}


/*
 * Check methods defined in InterpKernelDEC.hxx
 *
 InterpKernelDEC();
 InterpKernelDEC(ProcessorGroup& local_group, ProcessorGroup& distant_group);
 virtual ~InterpKernelDEC();
 void synchronize();
 void recvData();
 void sendData();
*/
 
void ParaMEDMEMTest::testInterpKernelDEC_2D_(const char *srcMeth, const char *targetMeth)
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
  
  MEDCoupling::CommInterface interface;
    
  MEDCoupling::ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  MEDCoupling::ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  MEDCoupling::ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  
  //loading the geometry for the source group

  MEDCoupling::InterpKernelDEC dec (*source_group,*target_group);

  MEDCoupling::MEDCouplingUMesh* mesh;
  MEDCoupling::ParaMESH* paramesh;
  MEDCoupling::ParaFIELD* parafield;
  ICoCo::MEDField* icocofield ;
  
  string filename_xml1              = "square1_split";
  string filename_xml2              = "square2_split";
  //string filename_seq_wr            = makeTmpFile("");
  //string filename_seq_med           = makeTmpFile("myWrField_seq_pointe221.med");
  
  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (source_group->containsMyRank())
    {
      string master = filename_xml1;
      
      ostringstream strstream;
      strstream <<master<<rank+1<<".med";
      string fName = INTERP_TEST::getResourceFile(strstream.str());
      ostringstream meshname ;
      meshname<< "Mesh_2_"<< rank+1;
      
      mesh=ReadUMeshFromFile(fName.c_str(),meshname.str().c_str(),0);
      
    
      paramesh=new ParaMESH (mesh,*source_group,"source mesh");
    
      //      MEDCoupling::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT( support,*source_group);
      MEDCoupling::ComponentTopology comptopo;
      if(srcM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(IntensiveMaximum);
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
      icocofield=new ICoCo::MEDField(parafield->getField());
      dec.setMethod(srcMeth);
      dec.attachLocalField(icocofield);
    }
  
  //loading the geometry for the target group
  if (target_group->containsMyRank())
    {
      string master= filename_xml2;
      ostringstream strstream;
      strstream << master<<(rank-nproc_source+1)<<".med";
      string fName = INTERP_TEST::getResourceFile(strstream.str());
      ostringstream meshname ;
      meshname<< "Mesh_3_"<<rank-nproc_source+1;
      mesh = ReadUMeshFromFile(fName.c_str(),meshname.str().c_str(),0);
      
      paramesh=new ParaMESH (mesh,*target_group,"target mesh");
      //      MEDCoupling::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT(support,*target_group);
      MEDCoupling::ComponentTopology comptopo;
      if(targetM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(IntensiveMaximum);
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
      icocofield=new ICoCo::MEDField(parafield->getField());
      dec.setMethod(targetMeth);
      dec.attachLocalField(icocofield);
    }
    
  
  //attaching a DEC to the source group 
  double field_before_int;
  double field_after_int;
  
  if (source_group->containsMyRank())
    { 
      field_before_int = parafield->getVolumeIntegral(0,true);
      dec.synchronize();
      cout<<"DEC usage"<<endl;
      dec.setForcedRenormalization(false);

      dec.sendData();
      ParaMEDLoader::WriteParaMesh("./sourcesquareb",paramesh);
      if (source_group->myRank()==0)
        aRemover.Register("./sourcesquareb");
      ostringstream filename;
      filename<<"./sourcesquareb_"<<source_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      //WriteField("./sourcesquareb",parafield->getField());
   
      dec.recvData();
      cout <<"writing"<<endl;
      ParaMEDLoader::WriteParaMesh("./sourcesquare",paramesh);
      if (source_group->myRank()==0)
        aRemover.Register("./sourcesquare");
      //WriteField("./sourcesquare",parafield->getField());
      
     
      filename<<"./sourcesquare_"<<source_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      field_after_int = parafield->getVolumeIntegral(0,true);
      
      
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
      ParaMEDLoader::WriteParaMesh("./targetsquareb",paramesh);
      //WriteField("./targetsquareb",parafield->getField());
      if (target_group->myRank()==0)
        aRemover.Register("./targetsquareb");
      ostringstream filename;
      filename<<"./targetsquareb_"<<target_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      dec.sendData();
      ParaMEDLoader::WriteParaMesh("./targetsquare",paramesh);
      //WriteField("./targetsquare",parafield->getField());
      
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
  cout << "end of InterpKernelDEC_2D test"<<endl;
}

void ParaMEDMEMTest::testInterpKernelDEC2_2D_(const char *srcMeth, const char *targetMeth)
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
  
  MEDCoupling::CommInterface interface;
    
  MEDCoupling::ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  MEDCoupling::ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  MEDCoupling::ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  
  //loading the geometry for the source group

  MEDCoupling::InterpKernelDEC dec (*source_group,*target_group);

  MEDCoupling::MEDCouplingUMesh* mesh;
  MEDCoupling::MEDCouplingFieldDouble* mcfield;
  
  string filename_xml1              = "square1_split";
  string filename_xml2              = "square2_split";
  
  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (source_group->containsMyRank())
    {
      string master = filename_xml1;
      
      ostringstream strstream;
      strstream <<master<<rank+1<<".med";
      string fName = INTERP_TEST::getResourceFile(strstream.str());
      ostringstream meshname ;
      meshname<< "Mesh_2_"<< rank+1;
      
      mesh=ReadUMeshFromFile(fName.c_str(),meshname.str().c_str(),0);
      MEDCoupling::ComponentTopology comptopo;
      if(srcM=="P0")
        {
          mcfield = MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
          mcfield->setMesh(mesh);
          DataArrayDouble *array=DataArrayDouble::New();
          array->alloc(mcfield->getNumberOfTuples(),1);
          mcfield->setArray(array);
          array->decrRef();
          mcfield->setNature(IntensiveMaximum);
        }
      else
        {
          mcfield = MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
          mcfield->setMesh(mesh);
          DataArrayDouble *array=DataArrayDouble::New();
          array->alloc(mcfield->getNumberOfTuples(),1);
          mcfield->setArray(array);
          array->decrRef();
        }
      int nb_local;
      if(srcM=="P0")
        nb_local=mesh->getNumberOfCells();
      else
        nb_local=mesh->getNumberOfNodes();
      double *value=mcfield->getArray()->getPointer();
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=1.0;
      dec.setMethod(srcMeth);
      dec.attachLocalField(mcfield);
      dec.attachLocalField(mcfield);
    }
  
  //loading the geometry for the target group
  if (target_group->containsMyRank())
    {
      string master= filename_xml2;
      ostringstream strstream;
      strstream << master<<(rank-nproc_source+1)<<".med";
      string fName = INTERP_TEST::getResourceFile(strstream.str());
      ostringstream meshname ;
      meshname<< "Mesh_3_"<<rank-nproc_source+1;
      mesh = ReadUMeshFromFile(fName.c_str(),meshname.str().c_str(),0);
      MEDCoupling::ComponentTopology comptopo;
      if(targetM=="P0")
        {
          mcfield = MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
          mcfield->setMesh(mesh);
          DataArrayDouble *array=DataArrayDouble::New();
          array->alloc(mcfield->getNumberOfTuples(),1);
          mcfield->setArray(array);
          array->decrRef();
          mcfield->setNature(IntensiveMaximum);
        }
      else
        {
          mcfield = MEDCouplingFieldDouble::New(ON_NODES,NO_TIME);
          mcfield->setMesh(mesh);
          DataArrayDouble *array=DataArrayDouble::New();
          array->alloc(mcfield->getNumberOfTuples(),1);
          mcfield->setArray(array);
          array->decrRef();
        }
      int nb_local;
      if(targetM=="P0")
        nb_local=mesh->getNumberOfCells();
      else
        nb_local=mesh->getNumberOfNodes();
      double *value=mcfield->getArray()->getPointer();
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=0.0;
      dec.setMethod(targetMeth);
      dec.attachLocalField(mcfield);
      dec.attachLocalField(mcfield);
    }
    
  
  //attaching a DEC to the source group 

  if (source_group->containsMyRank())
    { 
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.sendData();
      dec.recvData();
    }
  
  //attaching a DEC to the target group
  if (target_group->containsMyRank())
    {
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.recvData();
      dec.sendData();
    }
  delete source_group;
  delete target_group;
  delete self_group;
  mcfield->decrRef();
  mesh->decrRef();

  MPI_Barrier(MPI_COMM_WORLD);
  cout << "end of InterpKernelDEC2_2D test"<<endl;
}

void ParaMEDMEMTest::testInterpKernelDEC_3D_(const char *srcMeth, const char *targetMeth)
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
  
  MEDCoupling::CommInterface interface;
    
  MEDCoupling::ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  MEDCoupling::ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  MEDCoupling::ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  
  //loading the geometry for the source group

  MEDCoupling::InterpKernelDEC dec (*source_group,*target_group);

  MEDCoupling::MEDCouplingUMesh* mesh;
  MEDCoupling::ParaMESH* paramesh;
  MEDCoupling::ParaFIELD* parafield;
  ICoCo::MEDField* icocofield ;
  
  char * tmp_dir_c                    = getenv("TMP");
  string tmp_dir;
  if (tmp_dir_c != NULL)
    tmp_dir = string(tmp_dir_c);
  else
    tmp_dir = "/tmp";
  string filename_xml1              = "Mesh3D_10_2d";
  string filename_xml2              = "Mesh3D_11";
  //string filename_seq_wr            = makeTmpFile("");
  //string filename_seq_med           = makeTmpFile("myWrField_seq_pointe221.med");
  
  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (source_group->containsMyRank())
    {
      string master = filename_xml1;
      
      ostringstream strstream;
      strstream <<master<<rank+1<<".med";
      std::string fName = INTERP_TEST::getResourceFile(strstream.str());
      ostringstream meshname ;
      meshname<< "Mesh_3_"<< rank+1;
      
      mesh=ReadUMeshFromFile(fName.c_str(),meshname.str().c_str(),0);
      
    
      paramesh=new ParaMESH (mesh,*source_group,"source mesh");
    
      //      MEDCoupling::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT( support,*source_group);
      MEDCoupling::ComponentTopology comptopo;
      if(srcM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(IntensiveMaximum);
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
      icocofield=new ICoCo::MEDField(parafield->getField());
      dec.setMethod(srcMeth);
      dec.attachLocalField(icocofield);
    }
  
  //loading the geometry for the target group
  if (target_group->containsMyRank())
    {
      string master= filename_xml2;
      ostringstream strstream;
      strstream << master << ".med";
      std::string fName = INTERP_TEST::getResourceFile(strstream.str());
      ostringstream meshname ;
      meshname<< "Mesh_6";
      mesh = ReadUMeshFromFile(fName.c_str(),meshname.str().c_str(),0);
      
      paramesh=new ParaMESH (mesh,*target_group,"target mesh");
      //      MEDCoupling::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT(support,*target_group);
      MEDCoupling::ComponentTopology comptopo;
      if(targetM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(IntensiveMaximum);
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
      icocofield=new ICoCo::MEDField(parafield->getField());
      dec.setMethod(targetMeth);
      dec.attachLocalField(icocofield);
    }  
  //attaching a DEC to the source group 
  double field_before_int;
  double field_after_int;
  
  if (source_group->containsMyRank())
    { 
      field_before_int = parafield->getVolumeIntegral(0,true);
      dec.synchronize();
      cout<<"DEC usage"<<endl;
      dec.setForcedRenormalization(false);

      dec.sendData();
      ParaMEDLoader::WriteParaMesh("./sourcesquareb",paramesh);
      if (source_group->myRank()==0)
        aRemover.Register("./sourcesquareb");
      ostringstream filename;
      filename<<"./sourcesquareb_"<<source_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      //WriteField("./sourcesquareb",parafield->getField());
   
      dec.recvData();
      cout <<"writing"<<endl;
      ParaMEDLoader::WriteParaMesh("./sourcesquare",paramesh);
      if (source_group->myRank()==0)
        aRemover.Register("./sourcesquare");
      //WriteField("./sourcesquare",parafield->getField());
      
     
      filename<<"./sourcesquare_"<<source_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      field_after_int = parafield->getVolumeIntegral(0,true);

      CPPUNIT_ASSERT_DOUBLES_EQUAL(field_before_int, field_after_int, 1e-6);
    
    }
  
  //attaching a DEC to the target group
  if (target_group->containsMyRank())
    {
      dec.synchronize();
      dec.setForcedRenormalization(false);

      dec.recvData();
      ParaMEDLoader::WriteParaMesh("./targetsquareb",paramesh);
      //WriteField("./targetsquareb",parafield->getField());
      if (target_group->myRank()==0)
        aRemover.Register("./targetsquareb");
      ostringstream filename;
      filename<<"./targetsquareb_"<<target_group->myRank()+1;
      aRemover.Register(filename.str().c_str());
      dec.sendData();
      ParaMEDLoader::WriteParaMesh("./targetsquare",paramesh);
      //WriteField("./targetsquare",parafield->getField());
      
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
  cout << "end of InterpKernelDEC_3D test"<<endl;
}

//Synchronous tests without interpolation with native mode (AllToAll(v) from lam/MPI:
void ParaMEDMEMTest::testSynchronousEqualInterpKernelWithoutInterpNativeDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.1,1,0.1,1,false,false,false,"P0","P0");
}

//Synchronous tests without interpolation :
void ParaMEDMEMTest::testSynchronousEqualInterpKernelWithoutInterpDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.1,1,0.1,1,true,false,false,"P0","P0");
}

//Synchronous tests with interpolation :
void ParaMEDMEMTest::testSynchronousEqualInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.1,1,0.1,1,true,false,true,"P0","P0");
}
void ParaMEDMEMTest::testSynchronousFasterSourceInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.09,1,0.1,1,true,false,true,"P0","P0");
}
void ParaMEDMEMTest::testSynchronousSlowerSourceInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.11,1,0.1,1,true,false,true,"P0","P0");
}
void ParaMEDMEMTest::testSynchronousSlowSourceInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.11,1,0.01,1,true,false,true,"P0","P0");
}
void ParaMEDMEMTest::testSynchronousFastSourceInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.01,1,0.11,1,true,false,true,"P0","P0");
}

//Asynchronous tests with interpolation :
void ParaMEDMEMTest::testAsynchronousEqualInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.1,1,0.1,1,true,true,true,"P0","P0");
}
void ParaMEDMEMTest::testAsynchronousFasterSourceInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.09,1,0.1,1,true,true,true,"P0","P0");
}
void ParaMEDMEMTest::testAsynchronousSlowerSourceInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.11,1,0.1,1,true,true,true,"P0","P0");
}
void ParaMEDMEMTest::testAsynchronousSlowSourceInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.11,1,0.01,1,true,true,true,"P0","P0");
}
void ParaMEDMEMTest::testAsynchronousFastSourceInterpKernelDEC_2D()
{
  testAsynchronousInterpKernelDEC_2D(0.01,1,0.11,1,true,true,true,"P0","P0");
}

void ParaMEDMEMTest::testInterpKernelDECNonOverlapp_2D_P0P0()
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
  double targetResults4[3][2]={{8.5,17},{8.8461538461538431, 9.8461538461538449},{35.25,11.75}};
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
  MEDCoupling::MEDCouplingUMesh *mesh=0;
  MEDCoupling::ParaMESH *paramesh=0;
  MEDCoupling::ParaFIELD* parafield=0;
  //
  MEDCoupling::CommInterface interface;
  //
  ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
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
      MEDCoupling::ComponentTopology comptopo;
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
      MEDCoupling::ComponentTopology comptopo;
      parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
    }
  //test 1 - Conservative volumic
  MEDCoupling::InterpKernelDEC dec(*source_group,*target_group);
  parafield->getField()->setNature(IntensiveMaximum);
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
  //test 2 - ExtensiveMaximum
  MEDCoupling::InterpKernelDEC dec2(*source_group,*target_group);
  parafield->getField()->setNature(ExtensiveMaximum);
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
  //test 3 - ExtensiveMaximum with global constraint
  MEDCoupling::InterpKernelDEC dec3(*source_group,*target_group);
  parafield->getField()->setNature(ExtensiveConservation);
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
  //test 4 - IntensiveConservation
  MEDCoupling::InterpKernelDEC dec4(*source_group,*target_group);
  parafield->getField()->setNature(IntensiveConservation);
  if (source_group->containsMyRank())
    { 
      dec4.setMethod("P0");
      dec4.attachLocalField(parafield);
      dec4.synchronize();
      dec4.setForcedRenormalization(false);
      dec4.sendData();
    }
  else
    {
      dec4.setMethod("P0");
      dec4.attachLocalField(parafield);
      dec4.synchronize();
      dec4.setForcedRenormalization(false);
      dec4.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      const double *expected=targetResults4[rank-nproc_source];
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[0],res[0],1e-13);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[1],res[1],1e-13);
    }
  //test 5 - Conservative volumic reversed
  MEDCoupling::InterpKernelDEC dec5(*source_group,*target_group);
  parafield->getField()->setNature(IntensiveMaximum);
  if (source_group->containsMyRank())
    { 
      dec5.setMethod("P0");
      dec5.attachLocalField(parafield);
      dec5.synchronize();
      dec5.setForcedRenormalization(false);
      dec5.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_EQUAL(1,(int)parafield->getField()->getNumberOfTuples());
      const double expected[]={37.8518518518519,43.5333333333333};
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[rank],res[0],1e-13);
    }
  else
    {
      dec5.setMethod("P0");
      dec5.attachLocalField(parafield);
      dec5.synchronize();
      dec5.setForcedRenormalization(false);
      double *res=parafield->getField()->getArray()->getPointer();
      const double *toSet=targetResults[rank-nproc_source];
      res[0]=toSet[0];
      res[1]=toSet[1];
      dec5.sendData();
    }
  //test 6 - ExtensiveMaximum reversed
  MEDCoupling::InterpKernelDEC dec6(*source_group,*target_group);
  parafield->getField()->setNature(ExtensiveMaximum);
  if (source_group->containsMyRank())
    { 
      dec6.setMethod("P0");
      dec6.attachLocalField(parafield);
      dec6.synchronize();
      dec6.setForcedRenormalization(false);
      dec6.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_EQUAL(1,(int)parafield->getField()->getNumberOfTuples());
      const double expected[]={0.794600591715977,1.35631163708087};
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[rank],res[0],1e-13);
    }
  else
    {
      dec6.setMethod("P0");
      dec6.attachLocalField(parafield);
      dec6.synchronize();
      dec6.setForcedRenormalization(false);
      double *res=parafield->getField()->getArray()->getPointer();
      const double *toSet=targetResults2[rank-nproc_source];
      res[0]=toSet[0];
      res[1]=toSet[1];
      dec6.sendData();
    }
  //test 7 - ExtensiveMaximum with global constraint reversed
  MEDCoupling::InterpKernelDEC dec7(*source_group,*target_group);
  parafield->getField()->setNature(ExtensiveConservation);
  if (source_group->containsMyRank())
    { 
      dec7.setMethod("P0");
      dec7.attachLocalField(parafield);
      dec7.synchronize();
      dec7.setForcedRenormalization(false);
      dec7.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_EQUAL(1,(int)parafield->getField()->getNumberOfTuples());
      const double expected[]={36.4592592592593,44.5407407407407};
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[rank],res[0],1e-13);
    }
  else
    {
      dec7.setMethod("P0");
      dec7.attachLocalField(parafield);
      dec7.synchronize();
      dec7.setForcedRenormalization(false);
      double *res=parafield->getField()->getArray()->getPointer();
      const double *toSet=targetResults3[rank-nproc_source];
      res[0]=toSet[0];
      res[1]=toSet[1];
      dec7.sendData();
    }
  //test 8 - ExtensiveMaximum with IntensiveConservation reversed
  MEDCoupling::InterpKernelDEC dec8(*source_group,*target_group);
  parafield->getField()->setNature(IntensiveConservation);
  if (source_group->containsMyRank())
    { 
      dec8.setMethod("P0");
      dec8.attachLocalField(parafield);
      dec8.synchronize();
      dec8.setForcedRenormalization(false);
      dec8.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_EQUAL(1,(int)parafield->getField()->getNumberOfTuples());
      const double expected[]={0.81314102564102553,1.3428994082840233};
      CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[rank],res[0],1e-13);
    }
  else
    {
      dec8.setMethod("P0");
      dec8.attachLocalField(parafield);
      dec8.synchronize();
      dec8.setForcedRenormalization(false);
      double *res=parafield->getField()->getArray()->getPointer();
      const double *toSet=targetResults4[rank-nproc_source];
      res[0]=toSet[0];
      res[1]=toSet[1];
      dec8.sendData();
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

void ParaMEDMEMTest::testInterpKernelDECNonOverlapp_2D_P0P1P1P0()
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
  MEDCoupling::MEDCouplingUMesh *mesh=0;
  MEDCoupling::ParaMESH *paramesh=0;
  MEDCoupling::ParaFIELD *parafieldP0=0,*parafieldP1=0;
  //
  MEDCoupling::CommInterface interface;
  //
  ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if(source_group->containsMyRank())
    {
      if(rank==0)
        {
          double coords[6]={-0.3,-0.3, 0.7,0.7, 0.7,-0.3};
          int conn[3]={0,1,2};
          //int globalNode[3]={1,2,0};
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
          //int globalNode[3]={1,3,2};
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
      MEDCoupling::ComponentTopology comptopo;
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafieldP1 = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);
      double *valueP0=parafieldP0->getField()->getArray()->getPointer();
      double *valueP1=parafieldP1->getField()->getArray()->getPointer();
      parafieldP0->getField()->setNature(IntensiveMaximum);
      parafieldP1->getField()->setNature(IntensiveMaximum);
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
      const char targetMeshName[]="target mesh";
      if(rank==2)
        {
          double coords[10]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2 };
          int conn[7]={0,3,4,1, 1,4,2};
          //int globalNode[5]={4,3,0,2,1};
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
          paramesh=new ParaMESH(mesh,*target_group,targetMeshName);
          DataArrayInt *da=DataArrayInt::New();
          const int globalNumberingP2[5]={0,1,2,3,4};
          da->useArray(globalNumberingP2,false,CPP_DEALLOC,5,1);
          paramesh->setNodeGlobal(da);
          da->decrRef();
        }
      if(rank==3)
        {
          double coords[6]={0.2,0.2, 0.7,-0.3, 0.7,0.2};
          int conn[3]={0,2,1};
          //int globalNode[3]={1,0,5};
          mesh=MEDCouplingUMesh::New("Target mesh Proc3",2);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(3,2);
          std::copy(coords,coords+6,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          paramesh=new ParaMESH(mesh,*target_group,targetMeshName);
          DataArrayInt *da=DataArrayInt::New();
          const int globalNumberingP3[3]={4,2,5};
          da->useArray(globalNumberingP3,false,CPP_DEALLOC,3,1);
          paramesh->setNodeGlobal(da);
          da->decrRef();
        }
      if(rank==4)
        {
          double coords[12]={-0.3,0.2, -0.3,0.7, 0.2,0.7, 0.2,0.2, 0.7,0.7, 0.7,0.2};
          int conn[8]={0,1,2,3, 3,2,4,5};
          //int globalNode[6]={2,6,7,1,8,5};
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
          paramesh=new ParaMESH(mesh,*target_group,targetMeshName);
          DataArrayInt *da=DataArrayInt::New();
          const int globalNumberingP4[6]={3,6,7,4,8,5};
          da->useArray(globalNumberingP4,false,CPP_DEALLOC,6,1);
          paramesh->setNodeGlobal(da);
          da->decrRef();
        }
      MEDCoupling::ComponentTopology comptopo;
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafieldP1 = new ParaFIELD(ON_NODES,NO_TIME,paramesh, comptopo);
      parafieldP0->getField()->setNature(IntensiveMaximum);
      parafieldP1->getField()->setNature(IntensiveMaximum);
    }
  // test 1 - P0 P1
  MEDCoupling::InterpKernelDEC dec(*source_group,*target_group);
  if (source_group->containsMyRank())
    { 
      dec.setMethod("P0");
      dec.attachLocalField(parafieldP0);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.sendData();
      dec.recvData();
      const double *valueP0=parafieldP0->getField()->getArray()->getPointer();
      if(rank==0)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(34.42857143,valueP0[0],1e-7);
        }
      if(rank==1)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(44.,valueP0[0],1e-7);
        }
    }
  else
    {
      dec.setMethod("P1");
      dec.attachLocalField(parafieldP1);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.recvData();
      const double *res=parafieldP1->getField()->getArray()->getConstPointer();
      if(rank==2)
        {
          const double expectP2[5]={39.0, 31.0, 31.0, 47.0, 39.0};
          CPPUNIT_ASSERT_EQUAL(5,(int)parafieldP1->getField()->getNumberOfTuples());
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP1->getField()->getNumberOfComponents());
          for(int kk=0;kk<5;kk++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expectP2[kk],res[kk],1e-12);
        }
      if(rank==3)
        {
          const double expectP3[3]={39.0, 31.0, 31.0};
          CPPUNIT_ASSERT_EQUAL(3,(int)parafieldP1->getField()->getNumberOfTuples());
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP1->getField()->getNumberOfComponents());
          for(int kk=0;kk<3;kk++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expectP3[kk],res[kk],1e-12);
        }
      if(rank==4)
        {
          const double expectP4[6]={47.0, 47.0, 47.0, 39.0, 39.0, 31.0};
          CPPUNIT_ASSERT_EQUAL(6,(int)parafieldP1->getField()->getNumberOfTuples());
          CPPUNIT_ASSERT_EQUAL(1,(int)parafieldP1->getField()->getNumberOfComponents());
          for(int kk=0;kk<6;kk++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expectP4[kk],res[kk],1e-12);
        }
      dec.sendData();
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

void ParaMEDMEMTest::testInterpKernelDEC2DM1D_P0P0()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=3)
    return ;
  int nproc_source=2;
  set<int> procs_source;
  set<int> procs_target;
  //
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source;i<size; i++)
    procs_target.insert(i);
  //
  MEDCoupling::MEDCouplingUMesh *mesh=0;
  MEDCoupling::ParaMESH *paramesh=0;
  MEDCoupling::ParaFIELD *parafield=0;
  //
  MEDCoupling::CommInterface interface;
  //
  ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if(source_group->containsMyRank())
    {
      double targetCoords[18]={-0.3,-0.3, 0.2,-0.3, 0.7,-0.3, -0.3,0.2, 0.2,0.2, 0.7,0.2, -0.3,0.7, 0.2,0.7, 0.7,0.7 };
      mesh=MEDCouplingUMesh::New();
      mesh->setMeshDimension(2);
      DataArrayDouble *myCoords=DataArrayDouble::New();
      myCoords->alloc(9,2);
      std::copy(targetCoords,targetCoords+18,myCoords->getPointer());
      mesh->setCoords(myCoords);
      myCoords->decrRef();
      if(rank==0)
        {
          int targetConn[7]={0,3,4,1, 1,4,2};
          mesh->allocateCells(2);
          mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
          mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+4);
          mesh->finishInsertingCells();
        }
      else
        { 
          int targetConn[11]={4,5,2, 6,7,4,3, 7,8,5,4};
          mesh->allocateCells(3);
          mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
          mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+3);
          mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn+7);
          mesh->finishInsertingCells();
        }
      MEDCoupling::ComponentTopology comptopo;
      paramesh=new ParaMESH(mesh,*source_group,"source mesh");
      parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafield->getField()->setNature(IntensiveMaximum);
      double *vals=parafield->getField()->getArray()->getPointer();
      if(rank==0)
        { vals[0]=7.; vals[1]=8.; }
      else
        { vals[0]=9.; vals[1]=10.; vals[2]=11.; }
    }
  else
    {
      mesh=MEDCouplingUMesh::New("an example of -1 D mesh",-1);
      MEDCoupling::ComponentTopology comptopo;
      paramesh=new ParaMESH(mesh,*target_group,"target mesh");
      parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafield->getField()->setNature(IntensiveMaximum);
    }
  MEDCoupling::InterpKernelDEC dec(*source_group,*target_group);
  if(source_group->containsMyRank())
    {
      dec.setMethod("P0");
      dec.attachLocalField(parafield);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.sendData();
      dec.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      if(rank==0)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[0],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[1],1e-12);
        }
      else
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[0],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[1],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[2],1e-12);
        }
    }
  else
    {
      dec.setMethod("P0");
      dec.attachLocalField(parafield);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[0],1e-12);
      dec.sendData();
    }
  MEDCoupling::InterpKernelDEC dec2(*source_group,*target_group);
  dec2.setMethod("P0");
  parafield->getField()->setNature(ExtensiveConservation);
  if(source_group->containsMyRank())
    {
      double *vals=parafield->getField()->getArray()->getPointer();
      if(rank==0)
        { vals[0]=7.; vals[1]=8.; }
      else
        { vals[0]=9.; vals[1]=10.; vals[2]=11.; }
      dec2.attachLocalField(parafield);
      dec2.synchronize();
      dec2.sendData();
      dec2.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      if(rank==0)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(11.25,res[0],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(5.625,res[1],1e-12);
        }
      else
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(5.625,res[0],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(11.25,res[1],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(11.25,res[2],1e-12);
        }
    }
  else
    {
      dec2.attachLocalField(parafield);
      dec2.synchronize();
      dec2.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(45.,res[0],1e-12);
      dec2.sendData();
    }
  //
  MEDCoupling::InterpKernelDEC dec3(*source_group,*target_group);
  dec3.setMethod("P0");
  parafield->getField()->setNature(ExtensiveMaximum);
  if(source_group->containsMyRank())
    {
      double *vals=parafield->getField()->getArray()->getPointer();
      if(rank==0)
        { vals[0]=7.; vals[1]=8.; }
      else
        { vals[0]=9.; vals[1]=10.; vals[2]=11.; }
      dec3.attachLocalField(parafield);
      dec3.synchronize();
      dec3.sendData();
      dec3.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      if(rank==0)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(11.25,res[0],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(5.625,res[1],1e-12);
        }
      else
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(5.625,res[0],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(11.25,res[1],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(11.25,res[2],1e-12);
        }
    }
  else
    {
      dec3.attachLocalField(parafield);
      dec3.synchronize();
      dec3.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(45.,res[0],1e-12);
      dec3.sendData();
    }
  //
  MEDCoupling::InterpKernelDEC dec4(*source_group,*target_group);
  dec4.setMethod("P0");
  parafield->getField()->setNature(IntensiveConservation);
  if(source_group->containsMyRank())
    {
      double *vals=parafield->getField()->getArray()->getPointer();
      if(rank==0)
        { vals[0]=7.; vals[1]=8.; }
      else
        { vals[0]=9.; vals[1]=10.; vals[2]=11.; }
      dec4.attachLocalField(parafield);
      dec4.synchronize();
      dec4.sendData();
      dec4.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
       if(rank==0)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[0],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[1],1e-12);
        }
      else
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[0],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[1],1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[2],1e-12);
        }
    }
  else
    {
      dec4.attachLocalField(parafield);
      dec4.synchronize();
      dec4.recvData();
      const double *res=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(9.125,res[0],1e-12);
      dec4.sendData();
    }
  delete parafield;
  delete paramesh;
  mesh->decrRef();
  delete target_group;
  delete source_group;
  //
  MPI_Barrier(MPI_COMM_WORLD);
}

void ParaMEDMEMTest::testInterpKernelDECPartialProcs()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=3)
    return ;
  set<int> procs_source;
  set<int> procs_target;
  //
  procs_source.insert(0);
  procs_target.insert(1);
  //
  MEDCoupling::MEDCouplingUMesh *mesh=0;
  MEDCoupling::ParaMESH *paramesh=0;
  MEDCoupling::ParaFIELD *parafield=0;
  //
  MEDCoupling::CommInterface interface;
  //
  MPI_Barrier(MPI_COMM_WORLD);
  double targetCoords[8]={ 0.,0., 1., 0., 0., 1., 1., 1. };
  CommInterface comm;
  int grpIds[2]={0,1};
  MPI_Group grp,group_world;
  comm.commGroup(MPI_COMM_WORLD,&group_world);
  comm.groupIncl(group_world,2,grpIds,&grp);
  MPI_Comm partialComm;
  comm.commCreate(MPI_COMM_WORLD,grp,&partialComm);
  //
  ProcessorGroup* target_group=0;
  ProcessorGroup* source_group=0;
  //
  MEDCoupling::InterpKernelDEC *dec=0;
  if(rank==0 || rank==1)
    {
      target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target,partialComm);
      source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source,partialComm);
      if(source_group->containsMyRank())
        {    
          mesh=MEDCouplingUMesh::New();
          mesh->setMeshDimension(2);
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(4,2);
          std::copy(targetCoords,targetCoords+8,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          int targetConn[4]={0,2,3,1};
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,targetConn);
          mesh->finishInsertingCells();
          MEDCoupling::ComponentTopology comptopo;
          paramesh=new ParaMESH(mesh,*source_group,"source mesh");
          parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(IntensiveMaximum);
          double *vals=parafield->getField()->getArray()->getPointer();
          vals[0]=7.;
          dec=new MEDCoupling::InterpKernelDEC(*source_group,*target_group);
          dec->attachLocalField(parafield);
          dec->synchronize();
          dec->sendData();
          dec->recvData();
        }
      else
        {
          mesh=MEDCouplingUMesh::New();
          mesh->setMeshDimension(2);
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(4,2);
          std::copy(targetCoords,targetCoords+8,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          int targetConn[6]={0,2,1,2,3,1};
          mesh->allocateCells(2);
          mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn);
          mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,targetConn+3);
          mesh->finishInsertingCells();
          MEDCoupling::ComponentTopology comptopo;
          paramesh=new ParaMESH(mesh,*target_group,"target mesh");
          parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(IntensiveMaximum);
          dec=new MEDCoupling::InterpKernelDEC(*source_group,*target_group);
          dec->attachLocalField(parafield);
          dec->synchronize();
          dec->recvData();
          dec->sendData();
        }
    }
  delete parafield;
  delete paramesh;
  if(mesh)
    mesh->decrRef();
  delete target_group;
  delete source_group;
  delete dec;
  if(partialComm != MPI_COMM_NULL)
    comm.commFree(&partialComm);
  comm.groupFree(&grp);
  comm.groupFree(&group_world);
  MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * This test reproduces bug of Gauthier on 13/9/2010 concerning 3DSurf meshes.
 * It is possible to lead to dead lock in InterpKernelDEC when 3DSurfMeshes global bounding boxes intersects whereas cell bounding box intersecting only on one side.
 */
void ParaMEDMEMTest::testInterpKernelDEC3DSurfEmptyBBox()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=3)
    return ;
  int nproc_source = 1;
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  //
  MEDCoupling::MEDCouplingUMesh *mesh=0;
  MEDCoupling::ParaMESH *paramesh=0;
  MEDCoupling::ParaFIELD *parafieldP0=0;
  //
  MEDCoupling::CommInterface interface;
  //
  ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if(source_group->containsMyRank())
    {
      double coords[15]={1.,0.,0., 2.,0.,0., 2.,2.,0., 0.,2.,0., 0.5,0.5,1.};
      int conn[7]={0,1,2,3,0,3,4};
      mesh=MEDCouplingUMesh::New("Source mesh Proc0",2);
      mesh->allocateCells(2);
      mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
      mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+4);
      mesh->finishInsertingCells();
      DataArrayDouble *myCoords=DataArrayDouble::New();
      myCoords->alloc(5,3);
      std::copy(coords,coords+15,myCoords->getPointer());
      mesh->setCoords(myCoords);
      myCoords->decrRef();
      //
      paramesh=new ParaMESH(mesh,*source_group,"source mesh");
      MEDCoupling::ComponentTopology comptopo;
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      double *valueP0=parafieldP0->getField()->getArray()->getPointer();
      parafieldP0->getField()->setNature(IntensiveMaximum);
      valueP0[0]=7.; valueP0[1]=8.;
    }
  else
    {
      const char targetMeshName[]="target mesh";
      if(rank==1)
        {
          double coords[12]={0.25,0.25,0.5, 0.,0.25,0.5, 0.,0.,0.5, 0.25,0.,0.5};
          int conn[4]={0,1,2,3};
          mesh=MEDCouplingUMesh::New("Target mesh Proc1",2);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(4,3);
          std::copy(coords,coords+12,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          paramesh=new ParaMESH(mesh,*target_group,targetMeshName);
        }
      if(rank==2)
        {
          double coords[12]={0.,0.25,0.5, 0.,0.,0.5, -1.,0.,0.5, -1.,0.25,0.5};
          int conn[4]={0,1,2,3};
          mesh=MEDCouplingUMesh::New("Target mesh Proc2",2);
          mesh->allocateCells(1);
          mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(4,3);
          std::copy(coords,coords+12,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          paramesh=new ParaMESH(mesh,*target_group,targetMeshName);
        }
      MEDCoupling::ComponentTopology comptopo;
      parafieldP0 = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafieldP0->getField()->setNature(IntensiveMaximum);
    }
  // test 1
  MEDCoupling::InterpKernelDEC dec(*source_group,*target_group);
  if (source_group->containsMyRank())
    { 
      dec.setMethod("P0");
      dec.attachLocalField(parafieldP0);
      dec.synchronize();
      // dec.setForcedRenormalization(false);
      // dec.sendData();
      // dec.recvData();
      // const double *valueP0=parafieldP0->getField()->getArray()->getPointer();
      // if(rank==0)
      //   {
      //     CPPUNIT_ASSERT_DOUBLES_EQUAL(7.4,valueP0[0],1e-7);
      //     CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0540540540540544,valueP0[1],1e-7);
      //   }
      // if(rank==1)
      //   {
      //     CPPUNIT_ASSERT_DOUBLES_EQUAL(8.64054054054054,valueP0[0],1e-7);
      //   }
      // if(rank==2)
      //   {
      //     CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0540540540540544,valueP0[0],1e-7);
      //   }
    }
  else
    {
      dec.setMethod("P0");
      dec.attachLocalField(parafieldP0);
      dec.synchronize();
      // dec.setForcedRenormalization(false);
      // dec.recvData();
      // const double *res=parafieldP0->getField()->getArray()->getConstPointer();
      // if(rank==3)
      //   {
      //     CPPUNIT_ASSERT_EQUAL(1,parafieldP0->getField()->getNumberOfTuples());
      //     CPPUNIT_ASSERT_EQUAL(1,parafieldP0->getField()->getNumberOfComponents());
      //     CPPUNIT_ASSERT_DOUBLES_EQUAL(7.4,res[0],1e-12);
      //   }
      // if(rank==4)
      //   {
      //     CPPUNIT_ASSERT_EQUAL(1,parafieldP0->getField()->getNumberOfTuples());
      //     CPPUNIT_ASSERT_EQUAL(1,parafieldP0->getField()->getNumberOfComponents());
      //     CPPUNIT_ASSERT_DOUBLES_EQUAL(9.0540540540540526,res[0],1e-12);
      //   }
      // dec.sendData();
    }
  //
  delete parafieldP0;
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
void ParaMEDMEMTest::testAsynchronousInterpKernelDEC_2D(double dtA, double tmaxA, 
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
  
  MEDCoupling::CommInterface interface;
    
  MEDCoupling::ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  MEDCoupling::ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  MEDCoupling::ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
    
  //loading the geometry for the source group

  MEDCoupling::InterpKernelDEC dec (*source_group,*target_group);
  
  MEDCoupling::MEDCouplingUMesh* mesh;
  MEDCoupling::ParaMESH* paramesh;
  MEDCoupling::ParaFIELD* parafield;
  
  ICoCo::MEDField* icocofield ;

  char * tmp_dir_c                    = getenv("TMP");
  string tmp_dir;
  if (tmp_dir_c != NULL)
    tmp_dir = string(tmp_dir_c);
  else
    tmp_dir = "/tmp";
  string filename_xml1              = "square1_split";
  string filename_xml2              = "square2_split";
  //string filename_seq_wr            = makeTmpFile("");
  //string filename_seq_med           = makeTmpFile("myWrField_seq_pointe221.med");
  
  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;
  
  MPI_Barrier(MPI_COMM_WORLD);

  if (source_group->containsMyRank())
    {
      string master = filename_xml1;
      
      ostringstream strstream;
      strstream <<master<<rank+1<<".med";
      string fName = INTERP_TEST::getResourceFile(strstream.str());
      ostringstream meshname ;
      meshname<< "Mesh_2_"<< rank+1;
      
      mesh=ReadUMeshFromFile(fName.c_str(),meshname.str().c_str(),0);

      paramesh=new ParaMESH (mesh,*source_group,"source mesh");
    
      //      MEDCoupling::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT( support,*source_group);
      MEDCoupling::ComponentTopology comptopo;
      if(srcM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(IntensiveMaximum);//InvertIntegral);//IntensiveMaximum);
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
      icocofield=new ICoCo::MEDField(parafield->getField());
     
      dec.attachLocalField(icocofield);


    }
  
  //loading the geometry for the target group
  if (target_group->containsMyRank())
    {
      string master= filename_xml2;
      ostringstream strstream;
      strstream << master<<(rank-nproc_source+1)<<".med";
      string fName = INTERP_TEST::getResourceFile(strstream.str());
      ostringstream meshname ;
      meshname<< "Mesh_3_"<<rank-nproc_source+1;
      
      mesh = ReadUMeshFromFile(fName.c_str(),meshname.str().c_str(),0);

      paramesh=new ParaMESH (mesh,*target_group,"target mesh");
      //      MEDCoupling::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT(support,*target_group);
      MEDCoupling::ComponentTopology comptopo;
      if(targetM=="P0")
        {
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
          parafield->getField()->setNature(IntensiveMaximum);//InvertIntegral);//IntensiveMaximum);
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
      icocofield=new ICoCo::MEDField(parafield->getField());
      
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
          cout << "testAsynchronousInterpKernelDEC_2D" << rank << " time " << time
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
          cout << "testAsynchronousInterpKernelDEC_2D" << rank << " time " << time
               << " dtB " << dtB << " tmaxB " << tmaxB << endl ;
          dec.recvData( time );
          double vi = parafield->getVolumeIntegral(0,true);
          cout << "testAsynchronousInterpKernelDEC_2D" << rank << " time " << time
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

  cout << "testAsynchronousInterpKernelDEC_2D" << rank << " MPI_Barrier " << endl ;
 
  if (Asynchronous) MPI_Barrier(MPI_COMM_WORLD);
  cout << "end of InterpKernelDEC_2D test"<<endl;
}
