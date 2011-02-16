void ParaMEDMEMTest::testOverlapDEC()
{
  std::string srcM(srcMeth);
  std::string targetM(targetMeth);
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  if (size != 3) return ;
   
  int nproc = 3;
  set<int> procs;
  
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
      meshS=MEDCouplingUMesh::New();
      meshS->setMeshDimension(2);
      DataArrayDouble *myCoords=DataArrayDouble::New();
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
      parameshS=new ParaMESH(meshS,*dec.getGrp(),"source mesh");
      parafieldS=new ParaFIELD(ON_CELLS,NO_TIME,parameshS,comptopo);
      parafieldS->getField()->setNature(ConservativeVolumic);
      double *vals=parafieldS->getField()->getArray()->getPointer();
      vals[0]=7.; vals[1]=8.;
      //
      meshT=MEDCouplingUMesh::New();
      meshT->setMeshDimension(2);
      myCoords=DataArrayDouble::New();
      myCoords->alloc(3,2);
      std::copy(coordsT,coordsT+6,myCoords->getPointer());
      meshT->setCoords(myCoords);
      myCoords->decrRef();
      int connT[3]={0,2,1};
      meshT->allocateCells(1);
      meshT->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT);
      meshT->finishInsertingCells();
      parameshT=new ParaMESH(meshT,*dec.getGrp(),"target mesh");
      parafieldT=new ParaFIELD(ON_CELLS,NO_TIME,parameshT,comptopo);
      parafieldT->getField()->setNature(ConservativeVolumic);
    }
  //
  if(rank==1)
    {
      const double coordsS[10]={1.,0.,0.5,0.5,1.,0.5,0.5,1.,1.,1.};
      const double coordsT[6]={0.,0.,0.5,0.5,0.,1.};
      meshS=MEDCouplingUMesh::New();
      meshS->setMeshDimension(2);
      DataArrayDouble *myCoords=DataArrayDouble::New();
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
      parameshS=new ParaMESH(meshS,*dec.getGrp(),"source mesh");
      parafieldS=new ParaFIELD(ON_CELLS,NO_TIME,parameshS,comptopo);
      parafieldS->getField()->setNature(ConservativeVolumic);
      double *vals=parafieldS->getField()->getArray()->getPointer();
      vals[0]=9.; vals[1]=11.;
      //
      meshT=MEDCouplingUMesh::New();
      meshT->setMeshDimension(2);
      myCoords=DataArrayDouble::New();
      myCoords->alloc(3,2);
      std::copy(coordsT,coordsT+6,myCoords->getPointer());
      meshT->setCoords(myCoords);
      myCoords->decrRef();
      int connT[3]={0,2,1};
      meshT->allocateCells(1);
      meshT->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT);
      meshT->finishInsertingCells();
      parameshT=new ParaMESH(meshT,*dec.getGrp(),"target mesh");
      parafieldT=new ParaFIELD(ON_CELLS,NO_TIME,parameshT,comptopo);
      parafieldT->getField()->setNature(ConservativeVolumic);
    }
  //
  if(rank==2)
    {
      const double coordsS[8]={0.,0.5, 0.5,0.5, 0.,1., 0.5,1.};
      const double coordsT[6]={0.5,0.5,0.,1.,1.,1.};
      meshS=MEDCouplingUMesh::New();
      meshS->setMeshDimension(2);
      DataArrayDouble *myCoords=DataArrayDouble::New();
      myCoords->alloc(4,2);
      std::copy(coordsS,coordsS+8,myCoords->getPointer());
      meshS->setCoords(myCoords);
      myCoords->decrRef();
      int connS[4]={0,2,3,1};
      meshS->allocateCells(1);
      meshS->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,connS);
      meshS->finishInsertingCells();
      ParaMEDMEM::ComponentTopology comptopo;
      parameshS=new ParaMESH(meshS,*dec.getGrp(),"source mesh");
      parafieldS=new ParaFIELD(ON_CELLS,NO_TIME,parameshS,comptopo);
      parafieldS->getField()->setNature(ConservativeVolumic);
      double *vals=parafieldS->getField()->getArray()->getPointer();
      vals[0]=10.;
      //
      meshT=MEDCouplingUMesh::New();
      meshT->setMeshDimension(2);
      myCoords=DataArrayDouble::New();
      myCoords->alloc(3,2);
      std::copy(coordsT,coordsT+6,myCoords->getPointer());
      meshT->setCoords(myCoords);
      myCoords->decrRef();
      int connT[3]={0,1,2};
      meshT->allocateCells(1);
      meshT->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT);
      meshT->finishInsertingCells();
      parameshT=new ParaMESH(meshT,*dec.getGrp(),"target mesh");
      parafieldT=new ParaFIELD(ON_CELLS,NO_TIME,parameshT,comptopo);
      parafieldT->getField()->setNature(ConservativeVolumic);
    }
  dec.attachSourceLocalField(parafieldS);
  dec.attachTargetLocalField(parafieldT);
  dec.synchronize();
  dec.sendRecvData(true);
  //
  if(rank==0)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(9.75,parafieldT->getField()->getArray()->getIJ(0,0));
    }
  if(rank==1)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(8.5,parafieldT->getField()->getArray()->getIJ(0,0));
    }
  if(rank==2)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(10.5,parafieldT->getField()->getArray()->getIJ(0,0));
    }

  delete parafieldS;
  delete parafieldT;
  delete parameshS;
  delete parameshT;
  meshS->decrRef();
  meshT->decrRef();

  MPI_Barrier(MPI_COMM_WORLD);
}

