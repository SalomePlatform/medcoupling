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

using namespace ParaMEDMEM;

void ParaMEDMEMTest::testFabienAPI1()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=3)
    return ;
  int procs_source[1]={0};
  int procs_target[1]={1};
  //
  ParaMEDMEM::MEDCouplingUMesh *mesh=0;
  ParaMEDMEM::ParaMESH *paramesh=0;
  ParaMEDMEM::ParaFIELD *parafield=0;
  //
  ParaMEDMEM::CommInterface interface;
  //
  MPI_Barrier(MPI_COMM_WORLD);
  double targetCoords[8]={ 0.,0., 1., 0., 0., 1., 1., 1. };
  CommInterface comm;
  //
  ParaMEDMEM::InterpKernelDEC *dec=new ParaMEDMEM::InterpKernelDEC(procs_source,procs_source+1,procs_target,procs_target+1);
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
      ParaMEDMEM::ComponentTopology comptopo;
      paramesh=new ParaMESH(mesh,*dec->getSourceGrp(),"source mesh");
      parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafield->getField()->setNature(ConservativeVolumic);
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
      ParaMEDMEM::ComponentTopology comptopo;
      paramesh=new ParaMESH(mesh,*dec->getTargetGrp(),"target mesh");
      parafield=new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
      parafield->getField()->setNature(ConservativeVolumic);
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
