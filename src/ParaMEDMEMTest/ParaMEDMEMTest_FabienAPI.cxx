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
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "InterpKernelDEC.hxx"
#include "MEDCouplingUMesh.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ComponentTopology.hxx"

#include <set>

using namespace MEDCoupling;

void ParaMEDMEMTest::testFabienAPI1()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=3)
    return ;
  int procs_source_c[1]={0};
  std::set<int> procs_source(procs_source_c,procs_source_c+1);
  int procs_target_c[1]={1};
  std::set<int> procs_target(procs_target_c,procs_target_c+1);
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
  //
  MEDCoupling::InterpKernelDEC *dec=new MEDCoupling::InterpKernelDEC(procs_source,procs_target);
  if(dec->isInSourceSide())
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
      paramesh=new ParaMESH(mesh,*dec->getSourceGrp(),"source mesh");
      parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafield->getField()->setNature(IntensiveMaximum);
      double *vals=parafield->getField()->getArray()->getPointer();
      vals[0]=7.;
    }
  if(dec->isInTargetSide())
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
      paramesh=new ParaMESH(mesh,*dec->getTargetGrp(),"target mesh");
      parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafield->getField()->setNature(IntensiveMaximum);
    }
  dec->attachLocalField(parafield);
  dec->synchronize();
  dec->sendRecvData();
  if(dec->isInTargetSide())
    {
      const double *valsToTest=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsToTest[0],7.,1e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsToTest[1],7.,1e-14);
    }
  //
  delete parafield;
  delete paramesh;
  if(mesh)
    mesh->decrRef();
  delete dec;
  MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * Idem testFabienAPI1 except that procs are shuffled. Test of the good management of group translation in newly created communicator.
 */
void ParaMEDMEMTest::testFabienAPI2()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=3)
    return ;
  int procs_source_c[1]={2};//difference with testFabienAPI1
  std::set<int> procs_source(procs_source_c,procs_source_c+1);
  int procs_target_c[1]={1};
  std::set<int> procs_target(procs_target_c,procs_target_c+1);
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
  //
  MEDCoupling::InterpKernelDEC *dec=new MEDCoupling::InterpKernelDEC(procs_source,procs_target);
  if(dec->isInSourceSide())
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
      paramesh=new ParaMESH(mesh,*dec->getSourceGrp(),"source mesh");
      parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafield->getField()->setNature(IntensiveMaximum);
      double *vals=parafield->getField()->getArray()->getPointer();
      vals[0]=7.;
    }
  if(dec->isInTargetSide())
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
      paramesh=new ParaMESH(mesh,*dec->getTargetGrp(),"target mesh");
      parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafield->getField()->setNature(IntensiveMaximum);
    }
  dec->attachLocalField(parafield);
  dec->synchronize();
  dec->sendRecvData();
  if(dec->isInTargetSide())
    {
      const double *valsToTest=parafield->getField()->getArray()->getConstPointer();
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsToTest[0],7.,1e-14);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(valsToTest[1],7.,1e-14);
    }
  //
  delete parafield;
  delete paramesh;
  if(mesh)
    mesh->decrRef();
  delete dec;
  MPI_Barrier(MPI_COMM_WORLD);
}
