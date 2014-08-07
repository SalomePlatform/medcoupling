// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
#include "OverlapDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ComponentTopology.hxx"

#include "MEDCouplingUMesh.hxx"

#include <set>

void ParaMEDMEMTest::testOverlapDEC1()
{
  std::string srcM("P0");
  std::string targetM("P0");
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (size != 3) return ;
   
  int nproc = 3;
  std::set<int> procs;
  
  for (int i=0; i<nproc; i++)
    procs.insert(i);
  
  ParaMEDMEM::CommInterface interface;

  ParaMEDMEM::OverlapDEC dec(procs);

  ParaMEDMEM::MEDCouplingUMesh* meshS=0;
  ParaMEDMEM::MEDCouplingUMesh* meshT=0;
  ParaMEDMEM::ParaMESH* parameshS=0;
  ParaMEDMEM::ParaMESH* parameshT=0;
  ParaMEDMEM::ParaFIELD* parafieldS=0;
  ParaMEDMEM::ParaFIELD* parafieldT=0;
  
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    {
      const double coordsS[10]={0.,0.,0.5,0.,1.,0.,0.,0.5,0.5,0.5};
      const double coordsT[6]={0.,0.,1.,0.,1.,1.};
      meshS=ParaMEDMEM::MEDCouplingUMesh::New();
      meshS->setMeshDimension(2);
      ParaMEDMEM::DataArrayDouble *myCoords=ParaMEDMEM::DataArrayDouble::New();
      myCoords->alloc(5,2);
      std::copy(coordsS,coordsS+10,myCoords->getPointer());
      meshS->setCoords(myCoords);
      myCoords->decrRef();
      int connS[7]={0,3,4,1, 1,4,2};
      meshS->allocateCells(2);
      meshS->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,connS);
      meshS->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connS+4);
      meshS->finishInsertingCells();
      ParaMEDMEM::ComponentTopology comptopo;
      parameshS=new ParaMEDMEM::ParaMESH(meshS,*dec.getGrp(),"source mesh");
      parafieldS=new ParaMEDMEM::ParaFIELD(ParaMEDMEM::ON_CELLS,ParaMEDMEM::NO_TIME,parameshS,comptopo);
      parafieldS->getField()->setNature(ParaMEDMEM::ConservativeVolumic);//IntegralGlobConstraint
      double *valsS=parafieldS->getField()->getArray()->getPointer();
      valsS[0]=7.; valsS[1]=8.;
      //
      meshT=ParaMEDMEM::MEDCouplingUMesh::New();
      meshT->setMeshDimension(2);
      myCoords=ParaMEDMEM::DataArrayDouble::New();
      myCoords->alloc(3,2);
      std::copy(coordsT,coordsT+6,myCoords->getPointer());
      meshT->setCoords(myCoords);
      myCoords->decrRef();
      int connT[3]={0,2,1};
      meshT->allocateCells(1);
      meshT->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT);
      meshT->finishInsertingCells();
      parameshT=new ParaMEDMEM::ParaMESH(meshT,*dec.getGrp(),"target mesh");
      parafieldT=new ParaMEDMEM::ParaFIELD(ParaMEDMEM::ON_CELLS,ParaMEDMEM::NO_TIME,parameshT,comptopo);
      parafieldT->getField()->setNature(ParaMEDMEM::ConservativeVolumic);//IntegralGlobConstraint
      double *valsT=parafieldT->getField()->getArray()->getPointer();
      valsT[0]=7.;
    }
  //
  if(rank==1)
    {
      const double coordsS[10]={1.,0.,0.5,0.5,1.,0.5,0.5,1.,1.,1.};
      const double coordsT[6]={0.,0.,0.5,0.5,0.,1.};
      meshS=ParaMEDMEM::MEDCouplingUMesh::New();
      meshS->setMeshDimension(2);
      ParaMEDMEM::DataArrayDouble *myCoords=ParaMEDMEM::DataArrayDouble::New();
      myCoords->alloc(5,2);
      std::copy(coordsS,coordsS+10,myCoords->getPointer());
      meshS->setCoords(myCoords);
      myCoords->decrRef();
      int connS[7]={0,1,2, 1,3,4,2};
      meshS->allocateCells(2);
      meshS->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connS);
      meshS->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,connS+3);
      meshS->finishInsertingCells();
      ParaMEDMEM::ComponentTopology comptopo;
      parameshS=new ParaMEDMEM::ParaMESH(meshS,*dec.getGrp(),"source mesh");
      parafieldS=new ParaMEDMEM::ParaFIELD(ParaMEDMEM::ON_CELLS,ParaMEDMEM::NO_TIME,parameshS,comptopo);
      parafieldS->getField()->setNature(ParaMEDMEM::ConservativeVolumic);//IntegralGlobConstraint
      double *valsS=parafieldS->getField()->getArray()->getPointer();
      valsS[0]=9.; valsS[1]=11.;
      //
      meshT=ParaMEDMEM::MEDCouplingUMesh::New();
      meshT->setMeshDimension(2);
      myCoords=ParaMEDMEM::DataArrayDouble::New();
      myCoords->alloc(3,2);
      std::copy(coordsT,coordsT+6,myCoords->getPointer());
      meshT->setCoords(myCoords);
      myCoords->decrRef();
      int connT[3]={0,2,1};
      meshT->allocateCells(1);
      meshT->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT);
      meshT->finishInsertingCells();
      parameshT=new ParaMEDMEM::ParaMESH(meshT,*dec.getGrp(),"target mesh");
      parafieldT=new ParaMEDMEM::ParaFIELD(ParaMEDMEM::ON_CELLS,ParaMEDMEM::NO_TIME,parameshT,comptopo);
      parafieldT->getField()->setNature(ParaMEDMEM::ConservativeVolumic);//IntegralGlobConstraint
      double *valsT=parafieldT->getField()->getArray()->getPointer();
      valsT[0]=8.;
    }
  //
  if(rank==2)
    {
      const double coordsS[8]={0.,0.5, 0.5,0.5, 0.,1., 0.5,1.};
      const double coordsT[6]={0.5,0.5,0.,1.,1.,1.};
      meshS=ParaMEDMEM::MEDCouplingUMesh::New();
      meshS->setMeshDimension(2);
      ParaMEDMEM::DataArrayDouble *myCoords=ParaMEDMEM::DataArrayDouble::New();
      myCoords->alloc(4,2);
      std::copy(coordsS,coordsS+8,myCoords->getPointer());
      meshS->setCoords(myCoords);
      myCoords->decrRef();
      int connS[4]={0,2,3,1};
      meshS->allocateCells(1);
      meshS->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,connS);
      meshS->finishInsertingCells();
      ParaMEDMEM::ComponentTopology comptopo;
      parameshS=new ParaMEDMEM::ParaMESH(meshS,*dec.getGrp(),"source mesh");
      parafieldS=new ParaMEDMEM::ParaFIELD(ParaMEDMEM::ON_CELLS,ParaMEDMEM::NO_TIME,parameshS,comptopo);
      parafieldS->getField()->setNature(ParaMEDMEM::ConservativeVolumic);//IntegralGlobConstraint
      double *valsS=parafieldS->getField()->getArray()->getPointer();
      valsS[0]=10.;
      //
      meshT=ParaMEDMEM::MEDCouplingUMesh::New();
      meshT->setMeshDimension(2);
      myCoords=ParaMEDMEM::DataArrayDouble::New();
      myCoords->alloc(3,2);
      std::copy(coordsT,coordsT+6,myCoords->getPointer());
      meshT->setCoords(myCoords);
      myCoords->decrRef();
      int connT[3]={0,1,2};
      meshT->allocateCells(1);
      meshT->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT);
      meshT->finishInsertingCells();
      parameshT=new ParaMEDMEM::ParaMESH(meshT,*dec.getGrp(),"target mesh");
      parafieldT=new ParaMEDMEM::ParaFIELD(ParaMEDMEM::ON_CELLS,ParaMEDMEM::NO_TIME,parameshT,comptopo);
      parafieldT->getField()->setNature(ParaMEDMEM::ConservativeVolumic);//IntegralGlobConstraint
      double *valsT=parafieldT->getField()->getArray()->getPointer();
      valsT[0]=9.;
    }
  dec.attachSourceLocalField(parafieldS);
  dec.attachTargetLocalField(parafieldT);
  dec.synchronize();
  dec.sendRecvData(true);
  //
  if(rank==0)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(8.75,parafieldT->getField()->getArray()->getIJ(0,0),1e-12);
    }
  if(rank==1)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(8.5,parafieldT->getField()->getArray()->getIJ(0,0),1e-12);
    }
  if(rank==2)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(10.5,parafieldT->getField()->getArray()->getIJ(0,0),1e-12);
    }
  delete parafieldS;
  delete parafieldT;
  delete parameshS;
  delete parameshT;
  meshS->decrRef();
  meshT->decrRef();

  MPI_Barrier(MPI_COMM_WORLD);
}

