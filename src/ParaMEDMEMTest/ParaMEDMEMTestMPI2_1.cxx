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

using namespace MEDCoupling;

void MPI2ParaMEDMEMTest::testBasicMPI2_1()
{
  int lsize, lrank, gsize, grank;
  MPI_Comm gcom;
  std::string service = "SERVICE";
  std::ostringstream meshfilename, meshname;
  MEDCoupling::ParaMESH *paramesh=0;
  MEDCoupling::MEDCouplingUMesh *mesh;
  MEDCoupling::ParaFIELD *parafield=0;
  MEDCoupling::CommInterface *interface;
  MEDCoupling::MPIProcessorGroup *source, *target;

  MPI_Comm_size( MPI_COMM_WORLD, &lsize );
  MPI_Comm_rank( MPI_COMM_WORLD, &lrank );
  if(lsize!=2)
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
  interface = new MEDCoupling::CommInterface;
  source = new MEDCoupling::MPIProcessorGroup(*interface,0,lsize-1,gcom);
  target = new MEDCoupling::MPIProcessorGroup(*interface,lsize,gsize-1,gcom);

  const double sourceCoordsAll[2][8]={{0.4,0.5,0.4,1.5,1.6,1.5,1.6,0.5},
                                      {0.3,-0.5,1.6,-0.5,1.6,-1.5,0.3,-1.5}};
  
  int conn4All[8]={0,1,2,3,4,5,6,7};
  
  std::ostringstream stream; stream << "sourcemesh2D proc " << grank;
  mesh=MEDCouplingUMesh::New(stream.str().c_str(),2);
  mesh->allocateCells(2);
  mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn4All);
  mesh->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(4,2);
  const double *sourceCoords=sourceCoordsAll[grank];
  std::copy(sourceCoords,sourceCoords+8,myCoords->getPointer());
  mesh->setCoords(myCoords);
  myCoords->decrRef();
  paramesh=new ParaMESH(mesh,*source,"source mesh");
  MEDCoupling::ComponentTopology comptopo;
  parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
  double *value=parafield->getField()->getArray()->getPointer();
  value[0]=34+13*((double)grank);

  MEDCoupling::InterpKernelDEC dec(*source,*target);
  parafield->getField()->setNature(IntensiveMaximum);


  dec.setMethod("P0");
  dec.attachLocalField(parafield);
  dec.synchronize();
  dec.setForcedRenormalization(false);
  dec.sendData();
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
