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
#include "OverlapDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ComponentTopology.hxx"

#include "MEDCouplingUMesh.hxx"

#include <set>

using namespace std;

#include "MCAuto.hxx"
#include "MEDLoader.hxx"
#include "MEDLoaderBase.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingRemapper.hxx"

using namespace MEDCoupling;

typedef  MCAuto<MEDCouplingUMesh> MUMesh;
typedef  MCAuto<MEDCouplingFieldDouble> MFDouble;
typedef  MCAuto<DataArrayDouble> DADouble;

//void ParaMEDMEMTest::testOverlapDEC_LMEC_seq()
//{
//  //  T_SC_Trio_src.med  -- "SupportOf_"
//  //  T_SC_Trio_dst.med  -- "SupportOf_T_SC_Trio"
//  //  h_TH_Trio_src.med  -- "SupportOf_"
//  //  h_TH_Trio_dst.med  -- "SupportOf_h_TH_Trio"
//  string rep("/export/home/adrien/support/antoine_LMEC/");
//  string src_mesh_nam(rep + string("T_SC_Trio_src.med"));
//  string tgt_mesh_nam(rep + string("T_SC_Trio_dst.med"));
////  string src_mesh_nam(rep + string("h_TH_Trio_src.med"));
////  string tgt_mesh_nam(rep + string("h_TH_Trio_dst.med"));
//  MUMesh src_mesh=ReadUMeshFromFile(src_mesh_nam,"SupportOf_",0);
//  MUMesh tgt_mesh=ReadUMeshFromFile(tgt_mesh_nam,"SupportOf_T_SC_Trio",0);
////  MUMesh tgt_mesh=ReadUMeshFromFile(tgt_mesh_nam,"SupportOf_h_TH_Trio",0);
//
//  MFDouble srcField = MEDCouplingFieldDouble::New(ON_CELLS, ONE_TIME);
//  srcField->setMesh(src_mesh);
//  DataArrayDouble * dad = DataArrayDouble::New(); dad->alloc(src_mesh->getNumberOfCells(),1);
//  dad->fillWithValue(1.0);
//  srcField->setArray(dad);
//  srcField->setNature(IntensiveMaximum);
//
//  MEDCouplingRemapper remap;
//  remap.setOrientation(2); // always consider surface intersections as absolute areas.
//  remap.prepare(src_mesh, tgt_mesh, "P0P0");
//  MFDouble tgtField = remap.transferField(srcField, 1.0e+300);
//  tgtField->setName("result");
//  string out_nam(rep + string("adrien.med"));
//  WriteField(out_nam,tgtField, true);
//  cout << "wrote: " << out_nam << "\n";
//  double integ1 = 0.0, integ2 = 0.0;
//  srcField->integral(true, &integ1);
//  tgtField->integral(true, &integ2);
////  tgtField->reprQuickOverview(cout);
//  CPPUNIT_ASSERT_DOUBLES_EQUAL(integ1,integ2,1e-8);
//
//  dad->decrRef();
//}
//
//void ParaMEDMEMTest::testOverlapDEC_LMEC_para()
//{
//  using namespace MEDCoupling;
//
//  int size;
//  int rank;
//  MPI_Comm_size(MPI_COMM_WORLD,&size);
//  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
//
//  if (size != 1) return ;
//
//  int nproc = 1;
//  std::set<int> procs;
//
//  for (int i=0; i<nproc; i++)
//    procs.insert(i);
//
//  CommInterface interface;
//  OverlapDEC dec(procs);
//
//  ParaMESH* parameshS=0;
//  ParaMESH* parameshT=0;
//  ParaFIELD* parafieldS=0;
//  ParaFIELD* parafieldT=0;
//  MFDouble srcField;
//
//  // **** FILE LOADING
//  //  T_SC_Trio_src.med  -- "SupportOf_"
//  //  T_SC_Trio_dst.med  -- "SupportOf_T_SC_Trio"
//  //  h_TH_Trio_src.med  -- "SupportOf_"
//  //  h_TH_Trio_dst.med  -- "SupportOf_h_TH_Trio"
//  string rep("/export/home/adrien/support/antoine_LMEC/");
//  string src_mesh_nam(rep + string("T_SC_Trio_src.med"));
//  string tgt_mesh_nam(rep + string("T_SC_Trio_dst.med"));
//
//
//  MPI_Barrier(MPI_COMM_WORLD);
//  if(rank==0)
//    {
//    //  string src_mesh_nam(rep + string("h_TH_Trio_src.med"));
//    //  string tgt_mesh_nam(rep + string("h_TH_Trio_dst.med"));
//      MUMesh src_mesh=ReadUMeshFromFile(src_mesh_nam,"SupportOf_",0);
//      MUMesh tgt_mesh=ReadUMeshFromFile(tgt_mesh_nam,"SupportOf_T_SC_Trio",0);
//    //  MUMesh tgt_mesh=ReadUMeshFromFile(tgt_mesh_nam,"SupportOf_h_TH_Trio",0);
//
//      // **** SOURCE
//      srcField = MEDCouplingFieldDouble::New(ON_CELLS, ONE_TIME);
//      srcField->setMesh(src_mesh);
//      DataArrayDouble * dad = DataArrayDouble::New(); dad->alloc(src_mesh->getNumberOfCells(),1);
//      dad->fillWithValue(1.0);
//      srcField->setArray(dad);
//      srcField->setNature(IntensiveMaximum);
//
//      ComponentTopology comptopo;
//      parameshS = new ParaMESH(src_mesh,*dec.getGroup(),"source mesh");
//      parafieldS = new ParaFIELD(ON_CELLS,ONE_TIME,parameshS,comptopo);
//      parafieldS->getField()->setNature(IntensiveMaximum);//ExtensiveConservation
//      parafieldS->getField()->setArray(dad);
//
//      // **** TARGET
//      parameshT=new ParaMESH(tgt_mesh,*dec.getGroup(),"target mesh");
//      parafieldT=new ParaFIELD(ON_CELLS,ONE_TIME,parameshT,comptopo);
//      parafieldT->getField()->setNature(IntensiveMaximum);//ExtensiveConservation
//      parafieldT->getField()->getArray()->fillWithValue(1.0e300);
////      valsT[0]=7.;
//    }
//  dec.setOrientation(2);
//  dec.attachSourceLocalField(parafieldS);
//  dec.attachTargetLocalField(parafieldT);
//  dec.synchronize();
//  dec.sendRecvData(true);
//  //
//  if(rank==0)
//    {
//      double integ1 = 0.0, integ2 = 0.0;
//      MEDCouplingFieldDouble * tgtField;
//
//      srcField->integral(true, &integ1);
//      tgtField = parafieldT->getField();
////      tgtField->reprQuickOverview(cout);
//      tgtField->integral(true, &integ2);
//      tgtField->setName("result");
//      string out_nam(rep + string("adrien_para.med"));
//      WriteField(out_nam,tgtField, true);
//      cout << "wrote: " << out_nam << "\n";
//      CPPUNIT_ASSERT_DOUBLES_EQUAL(integ1,integ2,1e-8);
//    }
//  delete parafieldS;
//  delete parafieldT;
//  delete parameshS;
//  delete parameshT;
//
//  MPI_Barrier(MPI_COMM_WORLD);
//}
//
void prepareData1(int rank, NatureOfField nature,
                  MEDCouplingFieldDouble *& fieldS, MEDCouplingFieldDouble *& fieldT)
{
  if(rank==0)
    {
      const double coordsS[10]={0.,0.,0.5,0.,1.,0.,0.,0.5,0.5,0.5};
      const double coordsT[6]={0.,0.,1.,0.,1.,1.};
      MUMesh meshS=MEDCouplingUMesh::New();
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
      fieldS = MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
      DADouble arr = DataArrayDouble::New(); arr->alloc(meshS->getNumberOfCells(), 1);
      fieldS->setMesh(meshS); fieldS->setArray(arr);
      fieldS->setNature(nature);
      double *valsS=fieldS->getArray()->getPointer();
      valsS[0]=7.; valsS[1]=8.;
      //
      MUMesh meshT=MEDCouplingUMesh::New();
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
      fieldT = MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
      DADouble arr2 = DataArrayDouble::New(); arr2->alloc(meshT->getNumberOfCells(), 1);
      fieldT->setMesh(meshT);  fieldT->setArray(arr2);
      fieldT->setNature(nature);
      double *valsT=fieldT->getArray()->getPointer();
      valsT[0]=7.;
    }
  //
  if(rank==1)
    {
      const double coordsS[10]={1.,0.,0.5,0.5,1.,0.5,0.5,1.,1.,1.};
      const double coordsT[6]={0.,0.,0.5,0.5,0.,1.};
      MUMesh meshS=MEDCouplingUMesh::New();
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
      fieldS = MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
      DADouble arr = DataArrayDouble::New(); arr->alloc(meshS->getNumberOfCells(), 1);
      fieldS->setMesh(meshS); fieldS->setArray(arr);
      fieldS->setNature(nature);
      double *valsS=fieldS->getArray()->getPointer();
      valsS[0]=9.; valsS[1]=11.;
      //
      MUMesh meshT=MEDCouplingUMesh::New();
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
      fieldT = MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
      DADouble arr2 = DataArrayDouble::New(); arr2->alloc(meshT->getNumberOfCells(), 1);
      fieldT->setMesh(meshT);  fieldT->setArray(arr2);
      fieldT->setNature(nature);
      double *valsT=fieldT->getArray()->getPointer();
      valsT[0]=8.;
    }
  //
  if(rank==2)
    {
      const double coordsS[8]={0.,0.5, 0.5,0.5, 0.,1., 0.5,1.};
      const double coordsT[6]={0.5,0.5,0.,1.,1.,1.};
      MUMesh meshS=MEDCouplingUMesh::New();
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
      fieldS = MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
      DADouble arr = DataArrayDouble::New(); arr->alloc(meshS->getNumberOfCells(), 1);
      fieldS->setMesh(meshS); fieldS->setArray(arr);
      fieldS->setNature(nature);
      double *valsS=fieldS->getArray()->getPointer();
      valsS[0]=10.;
      //
      MUMesh meshT=MEDCouplingUMesh::New();
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
      fieldT = MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
      DADouble arr2 = DataArrayDouble::New(); arr2->alloc(meshT->getNumberOfCells(), 1);
      fieldT->setMesh(meshT); fieldT->setArray(arr2);
      fieldT->setNature(nature);
      double *valsT=fieldT->getArray()->getPointer();
      valsT[0]=9.;
    }
}

void prepareData2_buildOneSquare(MEDCouplingUMesh* & meshS_0, MEDCouplingUMesh* & meshT_0)
{
  const double coords[10] = {0.0,0.0,  0.0,1.0,  1.0,1.0,  1.0,0.0, 0.5,0.5};
  meshS_0 = MEDCouplingUMesh::New("source", 2);
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(5,2);
  std::copy(coords,coords+10,myCoords->getPointer());
  meshS_0->setCoords(myCoords);  myCoords->decrRef();
  int connS[4]={0,1,2,3};
  meshS_0->allocateCells(2);
  meshS_0->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,connS);
  //
  meshT_0 = MEDCouplingUMesh::New("target", 2);
  myCoords=DataArrayDouble::New();
  myCoords->alloc(5,2);
  std::copy(coords,coords+10,myCoords->getPointer());
  meshT_0->setCoords(myCoords);
  myCoords->decrRef();
  int connT[12]={0,1,4,  1,2,4,  2,3,4,  3,0,4};
  meshT_0->allocateCells(4);
  meshT_0->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT);
  meshT_0->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT+3);
  meshT_0->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT+6);
  meshT_0->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,connT+9);
}

/**
 * Prepare five (detached) QUAD4 disposed like this:
 *   (0)  (1)  (2)
 *   (3)  (4)
 *
 * On the target side the global mesh is identical except that each QUAD4 is split in 4 TRI3 (along the diagonals).
 * This is a case for two procs:
 *    - proc #0 has source squares 0,1,2 and target squares 0,3 (well, sets of TRI3s actually)
 *    - proc #1 has source squares 3,4 and target squares 1,2,4
 */
void prepareData2(int rank, ProcessorGroup * grp, NatureOfField nature,
                  MEDCouplingUMesh *& meshS, MEDCouplingUMesh *& meshT,
                  ParaMESH*& parameshS, ParaMESH*& parameshT,
                  ParaFIELD*& parafieldS, ParaFIELD*& parafieldT,
                  bool stripPartOfSource=false,
                  int fieldCompoNum=1)
{
  MEDCouplingUMesh *meshS_0 = 0, *meshT_0 = 0;
  prepareData2_buildOneSquare(meshS_0, meshT_0);

  if(rank==0)
    {
      const double tr1[] = {1.5, 0.0};
      MEDCouplingUMesh *meshS_1 = static_cast<MEDCouplingUMesh*>(meshS_0->deepCopy());
      meshS_1->translate(tr1);
      const double tr2[] = {3.0, 0.0};
      MEDCouplingUMesh *meshS_2 = static_cast<MEDCouplingUMesh*>(meshS_0->deepCopy());
      meshS_2->translate(tr2);

      std::vector<const MEDCouplingUMesh*> vec;
      vec.push_back(meshS_0);vec.push_back(meshS_1);
      if (!stripPartOfSource)
        vec.push_back(meshS_2);
      meshS = MEDCouplingUMesh::MergeUMeshes(vec);
      meshS_1->decrRef(); meshS_2->decrRef();

      ComponentTopology comptopo(fieldCompoNum);
      parameshS=new ParaMESH(meshS, *grp,"source mesh");
      parafieldS=new ParaFIELD(ON_CELLS,ONE_TIME,parameshS,comptopo);
      parafieldS->getField()->setNature(nature);
      double *valsS=parafieldS->getField()->getArray()->getPointer();
      for(int i=0; i < fieldCompoNum; i++)
        {
          valsS[i] = 1. * (10^i);
          valsS[fieldCompoNum+i] = 2. * (10^i);
          if (!stripPartOfSource)
            {
              valsS[2*fieldCompoNum+i] = 3. * (10^i);
            }
        }

      //
      const double tr3[] = {0.0, -1.5};
      MEDCouplingUMesh *meshT_3 = static_cast<MEDCouplingUMesh*>(meshT_0->deepCopy());
      meshT_3->translate(tr3);
      vec.clear();
      vec.push_back(meshT_0);vec.push_back(meshT_3);
      meshT = MEDCouplingUMesh::MergeUMeshes(vec);
      meshT_3->decrRef();

      parameshT=new ParaMESH(meshT,*grp,"target mesh");
      parafieldT=new ParaFIELD(ON_CELLS,ONE_TIME,parameshT,comptopo);
      parafieldT->getField()->setNature(nature);
    }
  //
  if(rank==1)
    {
      const double tr3[] = {0.0, -1.5};
      MEDCouplingUMesh *meshS_3 = static_cast<MEDCouplingUMesh*>(meshS_0->deepCopy());
      meshS_3->translate(tr3);
      const double tr4[] = {1.5, -1.5};
      MEDCouplingUMesh *meshS_4 = static_cast<MEDCouplingUMesh*>(meshS_0->deepCopy());
      meshS_4->translate(tr4);

      std::vector<const MEDCouplingUMesh*> vec;
      vec.push_back(meshS_3);vec.push_back(meshS_4);
      meshS = MEDCouplingUMesh::MergeUMeshes(vec);
      meshS_3->decrRef(); meshS_4->decrRef();

      ComponentTopology comptopo(fieldCompoNum);
      parameshS=new ParaMESH(meshS, *grp,"source mesh");
      parafieldS=new ParaFIELD(ON_CELLS,ONE_TIME,parameshS,comptopo);
      parafieldS->getField()->setNature(nature);
      double *valsS=parafieldS->getField()->getArray()->getPointer();
      for(int i=0; i < fieldCompoNum; i++)
        {
          valsS[i] = 4. * (10^i);
          valsS[fieldCompoNum+i] = 5. * (10^i);
        }

      //
      const double tr5[] = {1.5, 0.0};
      MEDCouplingUMesh *meshT_1 = static_cast<MEDCouplingUMesh*>(meshT_0->deepCopy());
      meshT_1->translate(tr5);
      const double tr6[] = {3.0, 0.0};
      MEDCouplingUMesh *meshT_2 = static_cast<MEDCouplingUMesh*>(meshT_0->deepCopy());
      meshT_2->translate(tr6);
      const double tr7[] = {1.5, -1.5};
      MEDCouplingUMesh *meshT_4 = static_cast<MEDCouplingUMesh*>(meshT_0->deepCopy());
      meshT_4->translate(tr7);

      vec.clear();
      vec.push_back(meshT_1);vec.push_back(meshT_2);vec.push_back(meshT_4);
      meshT = MEDCouplingUMesh::MergeUMeshes(vec);
      meshT_1->decrRef(); meshT_2->decrRef(); meshT_4->decrRef();

      parameshT=new ParaMESH(meshT,*grp,"target mesh");
      parafieldT=new ParaFIELD(ON_CELLS,ONE_TIME,parameshT,comptopo);
      parafieldT->getField()->setNature(nature);
    }
  meshS_0->decrRef();
  meshT_0->decrRef();
}

/*! Test case from the official doc of the OverlapDEC.
 *  WARNING: bounding boxes might be tweaked here to make the case more interesting (i.e. to avoid an all to all exchange
 *  between all procs).
 */
void testOverlapDEC_generic(int workSharingAlgo, double bbAdj)
{
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //  char hostname[256];
  //  printf("(%d) PID %d on localhost ready for attach\n", rank, getpid());
  //  fflush(stdout);

//    if (rank == 0)
//      {
//        int i=1, j=0;
//        while (i!=0)
//          j=2;
//      }

  if (size != 3) return ;
  int nproc = 3;
  std::set<int> procs;
  for (int i=0; i<nproc; i++)
    procs.insert(i);

  CommInterface interface;
  OverlapDEC dec(procs);
  MEDCouplingFieldDouble * mcfieldS=0, *mcfieldT=0;

  prepareData1(rank, IntensiveMaximum, mcfieldS, mcfieldT);

  // See comment in the caller:
  dec.setBoundingBoxAdjustmentAbs(bbAdj);
  dec.setWorkSharingAlgo(workSharingAlgo);  // just to ease debugging

  dec.attachSourceLocalField(mcfieldS);
  dec.attachTargetLocalField(mcfieldT);
  dec.synchronize();
//  dec.debugPrintWorkSharing(std::cout);
  dec.sendRecvData(true);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if(rank==0)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(8.75,mcfieldT->getArray()->getIJ(0,0),1e-12);
    }
  if(rank==1)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(8.5,mcfieldT->getArray()->getIJ(0,0),1e-12);
    }
  if(rank==2)
    {
      CPPUNIT_ASSERT_DOUBLES_EQUAL(10.5,mcfieldT->getArray()->getIJ(0,0),1e-12);
    }

  mcfieldS->decrRef();
  mcfieldT->decrRef();

  MPI_Barrier(MPI_COMM_WORLD);
}

void ParaMEDMEMTest::testOverlapDEC1()
{
  /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   *  HACK ON BOUNDING BOX TO MAKE THIS CASE SIMPLE AND USABLE IN DEBUG
   * Bounding boxes are slightly smaller than should be, thus localizing the work to be done
   * and avoiding every proc talking to everyone else.
   * Obviously this is NOT a good idea to do this in production code :-)
   * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   */
  testOverlapDEC_generic(0,-1.0e-12);
}

void ParaMEDMEMTest::testOverlapDEC1_bis()
{
   // Same BB hack as above
  testOverlapDEC_generic(1,-1.0e-12);
}

void ParaMEDMEMTest::testOverlapDEC1_ter()
{
   // Same BB hack as above
  testOverlapDEC_generic(2, -1.0e-12);
}


/*!
 * Same as testOverlapDEC1() but with regular bounding boxes. If you're looking for a nice debug case,
 * testOverlapDEC1() is identical in terms of geometry and field values, and more appropriate.
 */
void ParaMEDMEMTest::testOverlapDEC2()
{
  testOverlapDEC_generic(0,1.0e-12);
}

void ParaMEDMEMTest::testOverlapDEC2_bis()
{
  testOverlapDEC_generic(1,1.0e-12);
}

void ParaMEDMEMTest::testOverlapDEC2_ter()
{
  testOverlapDEC_generic(2,1.0e-12);
}


/*! Test focused on the mapping of cell IDs.
 * (i.e. when only part of the source/target mesh is transmitted)
 */
void ParaMEDMEMTest::testOverlapDEC3()
{
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int nproc = 2;
  if (size != nproc) return ;
  std::set<int> procs;
  for (int i=0; i<nproc; i++)
    procs.insert(i);

  CommInterface interface;
  OverlapDEC dec(procs);
  ProcessorGroup * grp = dec.getGroup();
  MEDCouplingUMesh* meshS=0, *meshT=0;
  ParaMESH* parameshS=0, *parameshT=0;
  ParaFIELD* parafieldS=0, *parafieldT=0;

  prepareData2(rank, grp, IntensiveMaximum, meshS, meshT, parameshS, parameshT, parafieldS, parafieldT);

  dec.attachSourceLocalField(parafieldS);
  dec.attachTargetLocalField(parafieldT);
  dec.synchronize();
  dec.sendRecvData(true);
  //
  MEDCouplingFieldDouble * resField = parafieldT->getField();
  if(rank==0)
    {
      CPPUNIT_ASSERT_EQUAL(8, (int)resField->getNumberOfTuples());
      for(int i=0;i<4;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0,resField->getArray()->getIJ(i,0),1e-12);
      for(int i=4;i<8;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0,resField->getArray()->getIJ(i,0),1e-12);
    }
  if(rank==1)
    {
      CPPUNIT_ASSERT_EQUAL(12, (int)resField->getNumberOfTuples());
      for(int i=0;i<4;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0,resField->getArray()->getIJ(i,0),1e-12);
      for(int i=4;i<8;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(3.0,resField->getArray()->getIJ(i,0),1e-12);
      for(int i=8;i<12;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0,resField->getArray()->getIJ(i,0),1e-12);
    }
  delete parafieldS;
  delete parafieldT;
  delete parameshS;
  delete parameshT;
  meshS->decrRef();
  meshT->decrRef();

  MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * Tests:
 *  - default value
 *  - multi-component fields
 */
void ParaMEDMEMTest::testOverlapDEC4()
{
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  int nproc = 2;
  if (size != nproc) return ;
  std::set<int> procs;
  for (int i=0; i<nproc; i++)
    procs.insert(i);

  CommInterface interface;
  OverlapDEC dec(procs);
  ProcessorGroup * grp = dec.getGroup();
  MEDCouplingUMesh* meshS=0, *meshT=0;
  ParaMESH* parameshS=0, *parameshT=0;
  ParaFIELD* parafieldS=0, *parafieldT=0;

  // As before, except than one of the source cell is removed, and that the field now has 2 components
  prepareData2(rank, grp, IntensiveMaximum, meshS, meshT, parameshS, parameshT, parafieldS, parafieldT,
               true, 2);
//  if (rank == 1)
//    {
//      int i=1, j=0;
//      while (i!=0)
//        j=2;
//    }

  dec.attachSourceLocalField(parafieldS);
  dec.attachTargetLocalField(parafieldT);
  double defVal = -300.0;
  dec.setDefaultValue(defVal);
  dec.synchronize();
  dec.sendRecvData(true);
  //
  MEDCouplingFieldDouble * resField = parafieldT->getField();
  if(rank==0)
    {
      CPPUNIT_ASSERT_EQUAL(8, (int)resField->getNumberOfTuples());
      for(int i=0;i<4;i++)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(1.0,resField->getArray()->getIJ(i*2,0),1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(10.0,resField->getArray()->getIJ(i*2+1,0),1e-12);
        }
      for(int i=4;i<8;i++)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(4.0,resField->getArray()->getIJ(i*2,0),1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(40.0,resField->getArray()->getIJ(i*2+1,0),1e-12);
        }
    }
  if(rank==1)
    {
      CPPUNIT_ASSERT_EQUAL(12, (int)resField->getNumberOfTuples());
      for(int i=0;i<4;i++)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(2.0,resField->getArray()->getIJ(i*2,0),1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(20.0,resField->getArray()->getIJ(i*2+1,0),1e-12);
        }
      // Default value should be here:
      for(int i=4;i<8;i++)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(defVal,resField->getArray()->getIJ(i*2,0),1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(defVal,resField->getArray()->getIJ(i*2+1,0),1e-12);
        }
      for(int i=8;i<12;i++)
        {
          CPPUNIT_ASSERT_DOUBLES_EQUAL(5.0,resField->getArray()->getIJ(i*2,0),1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(50.0,resField->getArray()->getIJ(i*2+1,0),1e-12);
        }
    }
  delete parafieldS;
  delete parafieldT;
  delete parameshS;
  delete parameshT;
  meshS->decrRef();
  meshT->decrRef();

  MPI_Barrier(MPI_COMM_WORLD);
}

