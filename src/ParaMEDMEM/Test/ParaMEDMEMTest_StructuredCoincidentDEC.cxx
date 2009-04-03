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
#include "StructuredCoincidentDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ICoCoMEDField.hxx"
#include "MEDLoader.hxx"

#include <string>

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES

using namespace std;
using namespace ParaMEDMEM;

/*
 * Check methods defined in StructuredCoincidentDEC.hxx
 *
 StructuredCoincidentDEC();
 StructuredCoincidentDEC(ProcessorGroup& local_group, ProcessorGroup& distant_group);
 virtual ~StructuredCoincidentDEC();
 void synchronize();
 void recvData();
 void sendData();
*/

void ParaMEDMEMTest::testStructuredCoincidentDEC() {
  string testname="ParaMEDMEM - testStructured CoincidentDEC";
  //  MPI_Init(&argc, &argv); 
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (size<4) {
    return;
  }

  ParaMEDMEM::CommInterface interface;

  ParaMEDMEM::MPIProcessorGroup self_group (interface,rank,rank);
  ParaMEDMEM::MPIProcessorGroup target_group(interface,3,size-1);
  ParaMEDMEM::MPIProcessorGroup source_group (interface,0,2);

  ParaMEDMEM::MEDCouplingUMesh* mesh;
  ParaMEDMEM::ParaMESH* paramesh;
  ParaMEDMEM::ParaFIELD* parafield;

  string data_dir = getenv("MED_ROOT_DIR");
  string tmp_dir = getenv("TMP");
  if (tmp_dir == "")
    tmp_dir = "/tmp";
  string filename_xml1 = data_dir
    + "/share/salome/resources/med/square1_split";
  string filename_2 = data_dir + "/share/salome/resources/med/square1.med";
  string filename_seq_wr = tmp_dir + "/";
  string filename_seq_med = tmp_dir + "/myWrField_seq_pointe221.med";

  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;

  //loading the geometry for the source group

  ParaMEDMEM::StructuredCoincidentDEC dec(source_group, target_group);

  MPI_Barrier(MPI_COMM_WORLD);
  if (source_group.containsMyRank()) {
    string master = filename_xml1;

    ostringstream strstream;
    strstream <<master<<rank+1<<".med";
    ostringstream meshname;
    meshname<< "Mesh_2_"<< rank+1;

    mesh=MEDLoader::ReadUMeshFromFile(strstream.str().c_str(),meshname.str().c_str(),0);
    

    paramesh=new ParaMESH (mesh,source_group,"source mesh");

    ParaMEDMEM::ComponentTopology comptopo(6);
    parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);

    int nb_local=mesh->getNumberOfCells();
    const int* global_numbering = paramesh->getGlobalNumberingCell();
    
    double *value=parafield->getField()->getArray()->getPointer();
    for(int ielem=0; ielem<nb_local;ielem++)
      for (int icomp=0; icomp<6; icomp++)
        value[ielem*6+icomp]=global_numbering[ielem]*6+icomp;

    ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);

    dec.attachLocalField(icocofield);
    dec.synchronize();
    dec.sendData();
    delete icocofield;
  }

  //loading the geometry for the target group
  if (target_group.containsMyRank()) {

    string meshname2("Mesh_2");
    mesh = MEDLoader::ReadUMeshFromFile(filename_2.c_str(),meshname2.c_str(),0);
    
    paramesh=new ParaMESH (mesh,self_group,"target mesh");
    ParaMEDMEM::ComponentTopology comptopo(6, &target_group);

    parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);

    int nb_local=mesh->getNumberOfCells();
    double *value=parafield->getField()->getArray()->getPointer();
    for (int ielem=0; ielem<nb_local; ielem++)
      for (int icomp=0; icomp<comptopo.nbLocalComponents(); icomp++)
        value[ielem*comptopo.nbLocalComponents()+icomp]=0.0;
    ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);

    dec.attachLocalField(icocofield);
    dec.synchronize();
    dec.recvData();

    //checking validity of field
    const double* recv_value = parafield->getField()->getArray()->getPointer();
    for (int i=0; i< nb_local; i++) {
      int first = comptopo.firstLocalComponent();
      for (int icomp = 0; icomp < comptopo.nbLocalComponents(); icomp++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(recv_value[i*comptopo.nbLocalComponents()+icomp],(double)(i*6+icomp+first),1e-12);
    }
    delete icocofield;
  }
  delete parafield;
  delete paramesh;
  mesh->decrRef();

  //  MPI_Barrier(MPI_COMM_WORLD);

}
