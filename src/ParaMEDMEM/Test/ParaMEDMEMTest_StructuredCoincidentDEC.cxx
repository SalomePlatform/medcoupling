// Copyright (C) 2006  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
// CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
//
// This library is distributed in the hope that it will be useful
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com

#include "ParaMEDMEMTest.hxx"
#include <cppunit/TestAssert.h>

#include "MEDMEM_Exception.hxx"
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "Topology.hxx"
#include "DEC.hxx"
#include "StructuredCoincidentDEC.hxx"
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
using namespace ParaMEDMEM;
using namespace MEDMEM;

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
	//	MPI_Init(&argc, &argv); 
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

	MEDMEM::MESH* mesh;
	MEDMEM::SUPPORT* support;
	ParaMEDMEM::ParaMESH* paramesh;
	ParaMEDMEM::ParaFIELD* parafield;

	string data_dir = getenv("MED_ROOT_DIR");
	string tmp_dir = getenv("TMP");
	if (tmp_dir == "")
		tmp_dir = "/tmp";
	string filename_xml1 = data_dir
			+ "/share/salome/resources/MedFiles/square1_split";
	string filename_2 = data_dir + "/share/salome/resources/MedFiles/square1.med";
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

		CPPUNIT_ASSERT_NO_THROW(mesh = new MESH(MED_DRIVER,strstream.str(),meshname.str()));
		support=new MEDMEM::SUPPORT(mesh,"all elements",MED_EN::MED_CELL);

		paramesh=new ParaMESH (*mesh,source_group,"source mesh");

		ParaMEDMEM::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT( support,source_group);
		ParaMEDMEM::ComponentTopology comptopo(6);
		parafield = new ParaFIELD(parasupport, comptopo);

		int nb_local=support->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS);
		const int* global_numbering = parasupport->getGlobalNumbering();
		double * value= new double[nb_local*6];
		for (int ielem=0; ielem<nb_local; ielem++)
			for (int icomp=0; icomp<6; icomp++)
				value[ielem*6+icomp]=global_numbering[ielem]*6+icomp;

		parafield->getField()->setValue(value);

		ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);

		dec.attachLocalField(icocofield);
		dec.synchronize();
		dec.sendData();
		delete icocofield;
		delete[] value;
		delete parasupport;
	}

	//loading the geometry for the target group
	if (target_group.containsMyRank()) {

		string meshname2("Mesh_2");
		CPPUNIT_ASSERT_NO_THROW(mesh = new MESH(MED_DRIVER,filename_2,meshname2));
		support=new MEDMEM::SUPPORT(mesh,"all elements",MED_EN::MED_CELL);
		
		paramesh=new ParaMESH (*mesh,self_group,"target mesh");
		ParaMEDMEM::ParaSUPPORT* parasupport=new UnstructuredParaSUPPORT(support,self_group);
		ParaMEDMEM::ComponentTopology comptopo(6, &target_group);

		parafield = new ParaFIELD(parasupport, comptopo);

		int nb_local=support->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS);
		double * value= new double[nb_local*comptopo.nbLocalComponents()];
		for (int ielem=0; ielem<nb_local; ielem++)
			for (int icomp=0; icomp<comptopo.nbLocalComponents(); icomp++)
				value[ielem*comptopo.nbLocalComponents()+icomp]=0.0;

		parafield->getField()->setValue(value);
		ICoCo::Field* icocofield=new ICoCo::MEDField(paramesh,parafield);

		dec.attachLocalField(icocofield);
		dec.synchronize();
		dec.recvData();

		//checking validity of field
		const double* recv_value = parafield->getField()->getValue();
		for (int i=0; i< nb_local; i++) {
			int first = comptopo.firstLocalComponent();
			for (int icomp = 0; icomp < comptopo.nbLocalComponents(); icomp++)
				CPPUNIT_ASSERT_DOUBLES_EQUAL(recv_value[i*comptopo.nbLocalComponents()+icomp],(double)(i*6+icomp+first),1e-12);
		}
		delete icocofield;
		delete [] value;
		delete parasupport;
	}
	delete parafield;
	delete paramesh;
	delete support;
	delete mesh;

//	MPI_Barrier(MPI_COMM_WORLD);

}
