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

#include <cppunit/extensions/HelperMacros.h>

#include "MPI2Connector.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "InterpKernelDEC.hxx"
#include "MPIProcessorGroup.hxx"
#include "CommInterface.hxx"

#include <mpi.h>
#include <iostream>
#include <stdlib.h>

class MPI2ParaMEDMEMTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE( MPI2ParaMEDMEMTest );
  CPPUNIT_TEST( testBasicMPI2_1 );
  CPPUNIT_TEST_SUITE_END();
public:
  void testBasicMPI2_1();
};

using namespace ParaMEDMEM;

void MPI2ParaMEDMEMTest::testBasicMPI2_1()
{
  int lsize, lrank, gsize, grank;
  MPI_Comm gcom;
  std::string service = "SERVICE";
  std::ostringstream meshfilename, meshname;
  ParaMEDMEM::ParaMESH *paramesh=0;
  ParaMEDMEM::MEDCouplingUMesh* mesh;
  ParaMEDMEM::ParaFIELD *parafield=0;
  ParaMEDMEM::CommInterface* interface;
  ParaMEDMEM::MPIProcessorGroup* source, *target;
  
  MPI_Comm_size( MPI_COMM_WORLD, &lsize );
  MPI_Comm_rank( MPI_COMM_WORLD, &lrank );
  if(lsize!=3)
    {
      CPPUNIT_ASSERT(false);
      return;
    }

  /* Connection to remote programm */
  MPI2Connector *mpio = new MPI2Connector;
  gcom = mpio->remoteMPI2Connect(service);
  
  MPI_Comm_size( gcom, &gsize );
  MPI_Comm_rank( gcom, &grank );
  if(gsize!=5)
    {
      CPPUNIT_ASSERT(false);
      return;
    }

  interface = new ParaMEDMEM::CommInterface;
  source = new ParaMEDMEM::MPIProcessorGroup(*interface,0,gsize-lsize-1,gcom);
  target = new ParaMEDMEM::MPIProcessorGroup(*interface,gsize-lsize,gsize-1,gcom);

  const double targetCoordsAll[3][16]={{0.7,1.45,0.7,1.65,0.9,1.65,0.9,1.45,  1.1,1.4,1.1,1.6,1.3,1.6,1.3,1.4},
                                       {0.7,-0.6,0.7,0.7,0.9,0.7,0.9,-0.6,  1.1,-0.7,1.1,0.6,1.3,0.6,1.3,-0.7},
                                       {0.7,-1.55,0.7,-1.35,0.9,-1.35,0.9,-1.55,  1.1,-1.65,1.1,-1.45,1.3,-1.45,1.3,-1.65}};
  int conn4All[8]={0,1,2,3,4,5,6,7};
  double targetResults[3][2]={{34.,34.},{38.333333333333336,42.666666666666664},{47.,47.}};

  std::ostringstream stream; stream << "targetmesh2D proc " << grank-(gsize-lsize);
  mesh=MEDCouplingUMesh::New(stream.str().c_str(),2);
  mesh->allocateCells(2);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn4All);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn4All+4);
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(8,2);
  const double *targetCoords=targetCoordsAll[grank-(gsize-lsize)];
  std::copy(targetCoords,targetCoords+16,myCoords->getPointer());
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  paramesh=new ParaMESH (mesh,*target,"target mesh");
  ParaMEDMEM::ComponentTopology comptopo;
  parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);

  ParaMEDMEM::InterpKernelDEC dec(*source,*target);
  parafield->getField()->setNature(ConservativeVolumic);

  dec.setMethod("P0");
  dec.attachLocalField(parafield);
  dec.synchronize();
  dec.setForcedRenormalization(false);
  dec.recvData();
  const double *res=parafield->getField()->getArray()->getConstPointer();
  const double *expected=targetResults[grank-(gsize-lsize)];
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[0],res[0],1e-13);
  CPPUNIT_ASSERT_DOUBLES_EQUAL(expected[1],res[1],1e-13);
  /* Deconnection of remote programm */
  mpio->remoteMPI2Disconnect(service);
  /* clean-up */
  delete mpio;
  delete parafield;
  mesh->decrRef();
  delete paramesh;
  delete source;
  delete target;
  delete interface;
}

CPPUNIT_TEST_SUITE_REGISTRATION( MPI2ParaMEDMEMTest );

#include "MPIMainTest.hxx"
