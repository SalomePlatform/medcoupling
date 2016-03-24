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

#ifdef MED_ENABLE_FVM

#include "ParaMEDMEMTest.hxx"
#include <cppunit/TestAssert.h>

#include "MEDMEM_Exception.hxx"
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "Topology.hxx"
#include "DEC.hxx"
#include "NonCoincidentDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "UnstructuredParaSUPPORT.hxx"
#include "ICoCoMEDField.hxx"

#include <string>

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES


using namespace std;
using namespace MEDCoupling;
using namespace MEDMEM;

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

void ParaMEDMEMTest::testNonCoincidentDEC_2D()
{

  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  //the test is meant to run on five processors
  if (size !=5) return ;

  testNonCoincidentDEC( "square1_split",
                        "Mesh_2",
                        "square2_split",
                        "Mesh_3",
                        3,
                        1e-6);
}

void ParaMEDMEMTest::testNonCoincidentDEC_3D()
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  //the test is meant to run on five processors
  if (size !=4) return ;

  testNonCoincidentDEC( "blade_12000_split2",
                        "Mesh_1",
                        "blade_3000_split2",
                        "Mesh_1",
                        2,
                        1e4);
}

void ParaMEDMEMTest::testNonCoincidentDEC(const string& filename1,
                                          const string& meshname1,
                                          const string& filename2,
                                          const string& meshname2,
                                          int nproc_source,
                                          double epsilon)
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

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

  MEDCoupling::ParaMESH* source_mesh=0;
  MEDCoupling::ParaMESH* target_mesh=0;
  MEDCoupling::ParaSUPPORT* parasupport=0;
  //loading the geometry for the source group

  MEDCoupling::NonCoincidentDEC dec (*source_group,*target_group);

  MEDMEM::MESH* mesh;
  MEDMEM::SUPPORT* support;
  MEDMEM::FIELD<double>* field;
  MEDCoupling::ParaMESH* paramesh;
  MEDCoupling::ParaFIELD* parafield;

  string filename_xml1              = getResourceFile(filename1);
  string filename_xml2              = getResourceFile(filename2);
  //string filename_seq_wr            = makeTmpFile("");
  //string filename_seq_med           = makeTmpFile("myWrField_seq_pointe221.med");

  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;
  //aRemover.Register(filename_seq_wr);
  //aRemover.Register(filename_seq_med);
  MPI_Barrier(MPI_COMM_WORLD);
  ICoCo::Field* icocofield;
  if (source_group->containsMyRank())
    {
      string master = filename_xml1;

      ostringstream strstream;
      strstream <<master<<rank+1<<".med";
      ostringstream meshname ;
      meshname<< meshname1<<"_"<< rank+1;

      CPPUNIT_ASSERT_NO_THROW(mesh = new MESH(MED_DRIVER,strstream.str(),meshname.str()));
      support=new MEDMEM::SUPPORT(mesh,"all elements",MED_EN::MED_CELL);

      paramesh=new ParaMESH (*mesh,*source_group,"source mesh");

      parasupport=new UnstructuredParaSUPPORT( support,*source_group);
      MEDCoupling::ComponentTopology comptopo;
      parafield = new ParaFIELD(parasupport, comptopo);


      int nb_local=support->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS);
      double * value= new double[nb_local];
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=1.0;
      parafield->getField()->setValue(value);

      icocofield=new ICoCo::MEDField(paramesh,parafield);

      dec.attachLocalField(icocofield);
      delete [] value;
    }

  //loading the geometry for the target group
  if (target_group->containsMyRank())
    {
      string master= filename_xml2;
      ostringstream strstream;
      strstream << master<<(rank-nproc_source+1)<<".med";
      ostringstream meshname ;
      meshname<< meshname2<<"_"<<rank-nproc_source+1;

      CPPUNIT_ASSERT_NO_THROW(mesh = new MESH(MED_DRIVER,strstream.str(),meshname.str()));
      support=new MEDMEM::SUPPORT(mesh,"all elements",MED_EN::MED_CELL);

      paramesh=new ParaMESH (*mesh,*target_group,"target mesh");
      parasupport=new UnstructuredParaSUPPORT(support,*target_group);
      MEDCoupling::ComponentTopology comptopo;
      parafield = new ParaFIELD(parasupport, comptopo);


      int nb_local=support->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS);
      double * value= new double[nb_local];
      for(int ielem=0; ielem<nb_local;ielem++)
        value[ielem]=0.0;
      parafield->getField()->setValue(value);
      icocofield=new ICoCo::MEDField(paramesh,parafield);

      dec.attachLocalField(icocofield);
      delete [] value;
    }


  //attaching a DEC to the source group
  double field_before_int;
  double field_after_int;

  if (source_group->containsMyRank())
    {
      field_before_int = parafield->getVolumeIntegral(1);
      MPI_Bcast(&field_before_int, 1,MPI_DOUBLE, 0,MPI_COMM_WORLD);
      dec.synchronize();
      cout<<"DEC usage"<<endl;
      dec.setOption("ForcedRenormalization",false);

      dec.sendData();
      //      paramesh->write(MED_DRIVER,"./sourcesquarenc");
      //parafield->write(MED_DRIVER,"./sourcesquarenc","boundary");


    }

  //attaching a DEC to the target group
  if (target_group->containsMyRank())
    {
      MPI_Bcast(&field_before_int, 1,MPI_DOUBLE, 0,MPI_COMM_WORLD);

      dec.synchronize();
      dec.setOption("ForcedRenormalization",false);
      dec.recvData();
      //paramesh->write(MED_DRIVER, "./targetsquarenc");
      //parafield->write(MED_DRIVER, "./targetsquarenc", "boundary");
      field_after_int = parafield->getVolumeIntegral(1);

    }
  MPI_Bcast(&field_before_int,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  MPI_Bcast(&field_after_int, 1,MPI_DOUBLE, size-1,MPI_COMM_WORLD);

  CPPUNIT_ASSERT_DOUBLES_EQUAL(field_before_int, field_after_int, epsilon);

  delete source_group;
  delete target_group;
  delete self_group;
  delete icocofield;
  delete paramesh;
  delete parafield;
  delete support;
  delete parasupport;
  delete mesh;
  MPI_Barrier(MPI_COMM_WORLD);

}
#endif
