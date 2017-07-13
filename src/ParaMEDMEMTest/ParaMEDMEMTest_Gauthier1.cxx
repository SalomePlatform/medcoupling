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
#include "DEC.hxx"
#include "InterpKernelDEC.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ComponentTopology.hxx"
#include "BlockTopology.hxx"

#include <set>
#include <time.h>
#include <iostream>
#include <assert.h>
#include <string>
#include <math.h>

using namespace std;
using namespace MEDCoupling;
using namespace ICoCo;

void afficheGauthier1(const ParaFIELD& field, const double *vals, int lgth)
{
  const DataArrayDouble *valsOfField(field.getField()->getArray());
  CPPUNIT_ASSERT_EQUAL(lgth,(int)valsOfField->getNumberOfTuples());
  for (int ele=0;ele<valsOfField->getNumberOfTuples();ele++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vals[ele],valsOfField->getIJ(ele,0),1e-12);
}

MEDCouplingUMesh *init_quadGauthier1(int is_master)
{
  MCAuto<MEDCouplingUMesh> m(MEDCouplingUMesh::New("champ_quad",2));
  MCAuto<DataArrayDouble> coo(DataArrayDouble::New());
  if(is_master)
    {
      const double dataCoo[24]={0,0,0,1,0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,1,1,1,1,1};
      coo->alloc(8,3);
      std::copy(dataCoo,dataCoo+24,coo->getPointer());
      const int conn[8]={0,1,3,2,4,5,7,6};
      m->allocateCells(2);
      m->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
      m->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+4);
    }
  else
    {
      coo->alloc(0,3);
      m->allocateCells(0);
    }
  m->setCoords(coo);
  return m.retn();
}

MEDCouplingUMesh *init_triangleGauthier1(int is_master)
{
  MCAuto<MEDCouplingUMesh> m(MEDCouplingUMesh::New("champ_triangle",2));
  MCAuto<DataArrayDouble> coo(DataArrayDouble::New());
  if(is_master)
    {
      const double dataCoo[24]={0,0,0,1,0,0,0,0,1,1,0,1,0,1,0,1,1,0,0,1,1,1,1,1};
      coo->alloc(8,3);
      std::copy(dataCoo,dataCoo+24,coo->getPointer());
      const int conn[12]={0,1,2,1,2,3,4,5,7,4,6,7};
      m->allocateCells(2);
      for(int i=0;i<4;i++)
        m->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+3*i);
    }
  else
    {
      coo->alloc(0,3);
      m->allocateCells(0);
    }
  m->setCoords(coo);
  return m.retn();
}


void ParaMEDMEMTest::testGauthier1()
{
  int num_cas=0;
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  
  int is_master=0;

  CommInterface comm;
  set<int> emetteur_ids;
  set<int> recepteur_ids;
  emetteur_ids.insert(0);
  if(size!=4)
    return;
  recepteur_ids.insert(1);
  if (size >2) 
    recepteur_ids.insert(2);
  if (size >2) 
    emetteur_ids.insert(3);
  if ((rank==0)||(rank==1)) 
    is_master=1;
  
  MPIProcessorGroup recepteur_group(comm,recepteur_ids);
  MPIProcessorGroup emetteur_group(comm,emetteur_ids);

  string cas;
  if (recepteur_group.containsMyRank())
    {
      cas="recepteur";
      //freopen("recpeteur.out","w",stdout);
      //freopen("recepteur.err","w",stderr);
    }
  else
    {
      cas="emetteur";
      // freopen("emetteur.out","w",stdout);
      //freopen("emetteur.err","w",stderr);
    }
  double expected[8][4]={
    {1.,1.,1.,1.},
    {40., 40., 1., 1.},
    {1.,1.,1e200,1e200},
    {40.,1.,1e200,1e200},
    {1.,1.,1.,1.},
    {40.,1.,1.,1.},
    {1.,1.,1e200,1e200},
    {20.5,1.,1e200,1e200}
  };
  int expectedLgth[8]={4,4,2,2,4,4,2,2};
  
  for (int send=0;send<2;send++)
    for (int rec=0;rec<2;rec++)
      {
        InterpKernelDEC dec_emetteur(emetteur_group, recepteur_group);
        MEDCoupling::ParaFIELD *champ_emetteur(0),*champ_recepteur(0);
        MEDCoupling::ParaMESH *paramesh(0);
        MCAuto<MEDCouplingUMesh> mesh;
        dec_emetteur.setOrientation(2);
        if (send==0)
          {
            mesh=init_quadGauthier1(is_master);
          }
        else
          {
            mesh=init_triangleGauthier1(is_master);
          }
        paramesh=new MEDCoupling::ParaMESH(mesh,recepteur_group.containsMyRank()?recepteur_group:emetteur_group,"emetteur mesh");
        MEDCoupling::ComponentTopology comptopo;
        champ_emetteur=new MEDCoupling::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
        champ_emetteur->getField()->setNature(IntensiveMaximum);
        champ_emetteur->setOwnSupport(true);
        if (rec==0)
          {
            mesh=init_triangleGauthier1(is_master);
          }
        else
          {
            mesh=init_quadGauthier1(is_master);
          }
        paramesh=new MEDCoupling::ParaMESH(mesh,recepteur_group.containsMyRank()?recepteur_group:emetteur_group,"recepteur mesh");
        champ_recepteur=new MEDCoupling::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
        champ_recepteur->getField()->setNature(IntensiveMaximum);
        champ_recepteur->setOwnSupport(true);
        if (cas=="emetteur") 
          {
            champ_emetteur->getField()->getArray()->fillWithValue(1.);
          }
  
  
        MPI_Barrier(MPI_COMM_WORLD);

        //clock_t clock0= clock ();
        int compti=0;

        bool init=true; // first time step ??
        bool stop=false;
        //boucle sur les pas de quads
        while (!stop) {
  
          compti++;
          //clock_t clocki= clock ();
          //cout << compti << " CLOCK " << (clocki-clock0)*1.e-6 << endl; 
          for (int non_unif=0;non_unif<2;non_unif++)
            {
              if (cas=="emetteur") 
                {
                  if (non_unif)
                    if(rank!=3)
                      champ_emetteur->getField()->getArray()->setIJ(0,0,40);
                }
              //bool ok=false; // Is the time interval successfully solved ?
    
              // Loop on the time interval tries
              if(1) {
      

                if (cas=="emetteur")
                  dec_emetteur.attachLocalField(champ_emetteur);
                else
                  dec_emetteur.attachLocalField(champ_recepteur);


                if(init) dec_emetteur.synchronize();
                init=false;

                if (cas=="emetteur") {
                  //    affiche(champ_emetteur);
                  dec_emetteur.sendData();
                }
                else if (cas=="recepteur")
                  {
                    dec_emetteur.recvData();
                    if (is_master)
                      afficheGauthier1(*champ_recepteur,expected[num_cas],expectedLgth[num_cas]);
                  }
                else
                  throw 0;
                MPI_Barrier(MPI_COMM_WORLD);
              }
              stop=true;
              num_cas++;
            }
        }
        delete champ_emetteur;
        delete champ_recepteur;
      }
}

void ParaMEDMEMTest::testGauthier2()
{
  std::cout << "testGauthier2\n";
  double valuesExpected1[2]={0.,0.};
  double valuesExpected2[2]={0.95,0.970625};
  
  double valuesExpected30[]={0., 0., 0.05, 0., 0., 0.15, 0., 0., 0.25, 0., 0., 0.35, 0., 0., 0.45, 0., 0., 0.55, 0., 0., 0.65, 0., 0., 0.75, 0., 0., 0.85, 0., 0., 0.95};
  double valuesExpected31[]={0.,  0.,  0.029375,  0.,  0.,  0.029375,  0.,  0.,  0.1,  0.,  0.,  0.1,  0.,  0.,  0.2,  0.,  0.,  0.2,  0.,  0.,  0.3,  0.,  0.,  0.3,  0.,  0.,  0.4,  0.,  0.,  0.4,  0.,  0.,  0.5,  0.,  0.,  0.5,  0.,  0.,  0.6,  0.,  0.,  0.6,  0.,  0.,  0.7,  0.,  0.,  0.7,  0.,  0.,  0.8,  0.,  0.,  0.8,  0.,  0.,  0.9,  0.,  0.,  0.9,  0.,  0.,  0.970625,  0.,  0.,  0.970625 };

  double *valuesExpected3[2]={valuesExpected30,valuesExpected31};

  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if (size <2)
    return ;
  CommInterface comm;
  set<int> Genepi_ids;
  set<int> entree_chaude_ids;
  Genepi_ids.insert(0);
  for (int i=1;i<size;i++)
    entree_chaude_ids.insert(i);
  for (int type=0;type<2;type++)
    {
      MPIProcessorGroup entree_chaude_group(comm,entree_chaude_ids);
      MPIProcessorGroup Genepi_group(comm,Genepi_ids);

      MEDCoupling::ParaFIELD *vitesse(0);
      InterpKernelDEC dec_vit_in_chaude(entree_chaude_group, Genepi_group);

      if ( entree_chaude_group.containsMyRank())
        {
          MCAuto<MEDCouplingUMesh> mesh(MEDCouplingUMesh::New("mesh",2));
          MCAuto<DataArrayDouble> arr(DataArrayDouble::New()); arr->alloc(63,3);
          const double cooData[189]={0.,0.,0.,0.5,0.,0.,0.5,0.05,0.,0.,0.1,0.,0.5,0.1,0.,0.5,0.15,0.,0.,0.2,0.,0.5,0.2,0.,0.5,0.25,0.,0.,0.3,0.,0.5,0.3,0.,0.5,0.35,0.,0.,0.4,0.,0.5,0.4,0.,0.5,0.45,0.,0.,0.5,0.,0.5,0.5,0.,0.5,0.55,0.,0.,0.6,0.,0.5,0.6,0.,0.5,0.65,0.,0.,0.7,0.,0.5,0.7,0.,0.5,0.75,0.,0.,0.8,0.,0.5,0.8,0.,0.5,0.85,0.,0.,0.9,0.,0.5,0.9,0.,0.5,0.95,0.,1.,0.,0.,1.,0.1,0.,1.,0.2,0.,1.,0.3,0.,1.,0.4,0.,1.,0.5,0.,1.,0.6,0.,1.,0.7,0.,1.,0.8,0.,1.,0.9,0.,1.,0.05,0.,1.,0.15,0.,1.,0.25,0.,1.,0.35,0.,1.,0.45,0.,1.,0.55,0.,1.,0.65,0.,1.,0.75,0.,1.,0.85,0.,1.,0.95,0.,1.,1.,0.,0.,1.,0.,0.5,1.,0.,0.,0.05,0.,0.,0.15,0.,0.,0.25,0.,0.,0.35,0.,0.,0.45,0.,0.,0.55,0.,0.,0.65,0.,0.,0.75,0.,0.,0.85,0.,0.,0.95,0.};
          std::copy(cooData,cooData+189,arr->getPointer());
          mesh->setCoords(arr);
          mesh->allocateCells(80);
          const int conn[240]={0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,2,1,31,5,4,32,8,7,33,11,10,34,14,13,35,17,16,36,20,19,37,23,22,38,26,25,39,29,28,30,40,2,31,41,5,32,42,8,33,43,11,34,44,14,35,45,17,36,46,20,37,47,23,38,48,26,39,49,29,31,2,40,32,5,41,33,8,42,34,11,43,35,14,44,36,17,45,37,20,46,38,23,47,39,26,48,50,29,49,3,2,4,6,5,7,9,8,10,12,11,13,15,14,16,18,17,19,21,20,22,24,23,25,27,26,28,51,29,52,31,4,2,32,7,5,33,10,8,34,13,11,35,16,14,36,19,17,37,22,20,38,25,23,39,28,26,50,52,29,0,2,53,3,5,54,6,8,55,9,11,56,12,14,57,15,17,58,18,20,59,21,23,60,24,26,61,27,29,62,3,53,2,6,54,5,9,55,8,12,56,11,15,57,14,18,58,17,21,59,20,24,60,23,27,61,26,51,62,29};
          for(int i=0;i<80;i++)
            mesh->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+3*i);
          MCAuto<MEDCouplingFieldDouble> f(MEDCouplingFieldDouble::New(ON_NODES,ONE_TIME));
          const double valsOfField[189]={0.,0.,0.,0.,0.,0.,0.,0.,0.05,0.,0.,0.1,0.,0.,0.1,0.,0.,0.15,0.,0.,0.2,0.,0.,0.2,0.,0.,0.25,0.,0.,0.3,0.,0.,0.3,0.,0.,0.35,0.,0.,0.4,0.,0.,0.4,0.,0.,0.45,0.,0.,0.5,0.,0.,0.5,0.,0.,0.55,0.,0.,0.6,0.,0.,0.6,0.,0.,0.65,0.,0.,0.7,0.,0.,0.7,0.,0.,0.75,0.,0.,0.8,0.,0.,0.8,0.,0.,0.85,0.,0.,0.9,0.,0.,0.9,0.,0.,0.95,0.,0.,0.,0.,0.,0.1,0.,0.,0.2,0.,0.,0.3,0.,0.,0.4,0.,0.,0.5,0.,0.,0.6,0.,0.,0.7,0.,0.,0.8,0.,0.,0.9,0.,0.,0.05,0.,0.,0.15,0.,0.,0.25,0.,0.,0.35,0.,0.,0.45,0.,0.,0.55,0.,0.,0.65,0.,0.,0.75,0.,0.,0.85,0.,0.,0.95,0.,0.,1.,0.,0.,1.,0.,0.,1.,0.,0.,0.05,0.,0.,0.15,0.,0.,0.25,0.,0.,0.35,0.,0.,0.45,0.,0.,0.55,0.,0.,0.65,0.,0.,0.75,0.,0.,0.85,0.,0.,0.95};
          f->setMesh(mesh); f->setName("VITESSE_P1_OUT");
          arr=DataArrayDouble::New(); arr->alloc(63,3);
          std::copy(valsOfField,valsOfField+189,arr->getPointer());
          f->setArray(arr); f->setNature(IntensiveMaximum);
          MEDCoupling::ParaMESH *paramesh(new MEDCoupling::ParaMESH(mesh,entree_chaude_group,"emetteur mesh"));
          vitesse=new MEDCoupling::ParaFIELD(f,paramesh,entree_chaude_group);
          vitesse->setOwnSupport(true);
          dec_vit_in_chaude.setMethod("P1");
        }
      else
        {
          MCAuto<MEDCouplingUMesh> mesh(MEDCouplingUMesh::New("mesh",2));
          MCAuto<DataArrayDouble> arr(DataArrayDouble::New()); arr->alloc(22,3);
          const double cooData[66]={0,0,0,1,0,0,0,0.1,0,1,0.1,0,0,0.2,0,1,0.2,0,0,0.3,0,1,0.3,0,0,0.4,0,1,0.4,0,0,0.5,0,1,0.5,0,0,0.6,0,1,0.6,0,0,0.7,0,1,0.7,0,0,0.8,0,1,0.8,0,0,0.9,0,1,0.9,0,0,1,0,1,1,0};
          std::copy(cooData,cooData+66,arr->getPointer());
          mesh->setCoords(arr);
          mesh->allocateCells(10);
          const int conn[40]={0,1,3,2,2,3,5,4,4,5,7,6,6,7,9,8,8,9,11,10,10,11,13,12,12,13,15,14,14,15,17,16,16,17,19,18,18,19,21,20};
          for(int i=0;i<10;i++)
            mesh->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+4*i);
          MCAuto<MEDCouplingFieldDouble> f(MEDCouplingFieldDouble::New(type==0?ON_CELLS:ON_NODES,ONE_TIME));
          f->setMesh(mesh); f->setName("vitesse_in_chaude");
          arr=DataArrayDouble::New(); arr->alloc(f->getNumberOfTuplesExpected()*3); arr->fillWithZero(); arr->rearrange(3);
          f->setArray(arr); f->setNature(IntensiveMaximum);
          MEDCoupling::ParaMESH *paramesh(new MEDCoupling::ParaMESH(mesh,Genepi_group,"recepteur mesh"));
          vitesse=new MEDCoupling::ParaFIELD(f,paramesh,Genepi_group);
          vitesse->setOwnSupport(true);
          dec_vit_in_chaude.setMethod(f->getDiscretization()->getRepr());
        }

      dec_vit_in_chaude.attachLocalField(vitesse);
      
      dec_vit_in_chaude.synchronize();
  
  
      // Envois - receptions
      if (entree_chaude_group.containsMyRank())
        {
          dec_vit_in_chaude.sendData();
        }
      else
        {
          dec_vit_in_chaude.recvData(); 
        }
      if ( !entree_chaude_group.containsMyRank() )
        {
          double pmin=1e38, pmax=-1e38;
          const double *p(vitesse->getField()->getArray()->begin());
          for(std::size_t i=0;i<vitesse->getField()->getArray()->getNbOfElems();i++,p++)
            {
              if (*p<pmin) pmin=*p;
              if (*p>pmax) pmax=*p;
            }
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected1[type],pmin,1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[type],pmax,1e-12);
      
          int nbCompo(vitesse->getField()->getNumberOfComponents());
          p=vitesse->getField()->getArray()->begin();
          for(int i=0;i<vitesse->getField()->getNumberOfTuples();i++)
            for(int c=0;c<nbCompo;c++,p++)
              CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected3[type][i*nbCompo+c],*p,1e-12);
        }
      delete vitesse;
    }
}

void ParaMEDMEMTest::testGauthier3_1()
{
  testGauthier3_GEN(true,4);
}

void ParaMEDMEMTest::testGauthier3_2()
{
  testGauthier3_GEN(false,4);
}

void ParaMEDMEMTest::testGauthier3_3()
{
  testGauthier3_GEN(true,5);
}

void ParaMEDMEMTest::testGauthier3_4()
{
  int size;
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  if(size!=5)
    return;

  // Should throw since the two groups (source/target) do not form a partition of
  // all the procs.
  CPPUNIT_ASSERT_THROW(testGauthier3_GEN(false,5), INTERP_KERNEL::Exception);
}


/*!
 * Non regression test testing copy constructor of InterpKernelDEC. 
 */
void ParaMEDMEMTest::testGauthier3_GEN(bool withIDs, int nprocs)
{
  int num_cas=0;
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  int is_master=0;

  CommInterface comm;
  set<int> emetteur_ids;
  set<int> recepteur_ids;
  emetteur_ids.insert(0);
  if(size!=nprocs)
    return;
  recepteur_ids.insert(1);

  recepteur_ids.insert(size-2);

  emetteur_ids.insert(size-1);
  if ((rank==0)||(rank==1))
    is_master=1;

  MPIProcessorGroup recepteur_group(comm,recepteur_ids);
  MPIProcessorGroup emetteur_group(comm,emetteur_ids);

  string cas;
  if (recepteur_group.containsMyRank())
    {
      cas="recepteur";
      //freopen("recpeteur.out","w",stdout);
      //freopen("recepteur.err","w",stderr);
    }
  else
    {
      if (emetteur_group.containsMyRank())
        cas="emetteur";
      else
        cas="vide";
      // freopen("emetteur.out","w",stdout);
      //freopen("emetteur.err","w",stderr);
    }

  double expected[8][4]={
    {1.,1.,1.,1.},
    {40., 40., 1., 1.},
    {1.,1.,1e200,1e200},
    {40.,1.,1e200,1e200},
    {1.,1.,1.,1.},
    {40.,1.,1.,1.},
    {1.,1.,1e200,1e200},
    {20.5,1.,1e200,1e200}
  };
  int expectedLgth[8]={4,4,2,2,4,4,2,2};

  for (int send=0;send<2;send++)
    for (int rec=0;rec<2;rec++)
      {
        std::vector<InterpKernelDEC> decu(1);
        if (withIDs)
          decu[0] = InterpKernelDEC(emetteur_ids,recepteur_ids);
        else
          decu[0] = InterpKernelDEC(emetteur_group,recepteur_group);
        InterpKernelDEC& dec_emetteur=decu[0];
        MEDCoupling::ParaFIELD *champ_emetteur(0),*champ_recepteur(0);
        MEDCoupling::ParaMESH *paramesh(0);
        MCAuto<MEDCouplingUMesh> mesh;
        dec_emetteur.setOrientation(2);
        if (send==0)
          {
            mesh=init_quadGauthier1(is_master);
          }
        else
          {
            mesh=init_triangleGauthier1(is_master);
          }
        if (cas!="vide")
          {
            paramesh=new MEDCoupling::ParaMESH(mesh,recepteur_group.containsMyRank()?recepteur_group:emetteur_group,"emetteur mesh");
            MEDCoupling::ComponentTopology comptopo;
            champ_emetteur=new MEDCoupling::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
            champ_emetteur->getField()->setNature(IntensiveMaximum);
            champ_emetteur->setOwnSupport(true);
            if (rec==0)
              {
                mesh=init_triangleGauthier1(is_master);
              }
            else
              {
                mesh=init_quadGauthier1(is_master);
              }
            paramesh=new MEDCoupling::ParaMESH(mesh,recepteur_group.containsMyRank()?recepteur_group:emetteur_group,"recepteur mesh");
            champ_recepteur=new MEDCoupling::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
            champ_recepteur->getField()->setNature(IntensiveMaximum);
            champ_recepteur->setOwnSupport(true);
            if (cas=="emetteur")
              {
                champ_emetteur->getField()->getArray()->fillWithValue(1.);
              }
          }

        MPI_Barrier(MPI_COMM_WORLD);

        //clock_t clock0= clock ();
        int compti=0;

        bool init=true; // first time step ??
        bool stop=false;
        //boucle sur les pas de quads
        while (!stop) {

            compti++;
            //clock_t clocki= clock ();
            //cout << compti << " CLOCK " << (clocki-clock0)*1.e-6 << endl;
            for (int non_unif=0;non_unif<2;non_unif++)
              {
                if (cas=="emetteur")
                  {
                    if (non_unif)
                      if(rank!=3)
                        champ_emetteur->getField()->getArray()->setIJ(0,0,40);
                  }
                //bool ok=false; // Is the time interval successfully solved ?

                // Loop on the time interval tries
                if(1) {


                    if (cas=="emetteur")
                      dec_emetteur.attachLocalField(champ_emetteur);
                    else
                      dec_emetteur.attachLocalField(champ_recepteur);


                    if(init) dec_emetteur.synchronize();
                    init=false;

                    if (cas=="emetteur") {
                        //    affiche(champ_emetteur);
                        dec_emetteur.sendData();
                    }
                    else if (cas=="recepteur")
                      {
                        dec_emetteur.recvData();
                        if (is_master)
                          afficheGauthier1(*champ_recepteur,expected[num_cas],expectedLgth[num_cas]);
                      }
                    // else
                    //   throw 0;
                    MPI_Barrier(MPI_COMM_WORLD);
                }
                stop=true;
                num_cas++;
              }
        }
        delete champ_emetteur;
        delete champ_recepteur;
      }
}

/*!
 * This test is the parallel version of MEDCouplingBasicsTest.test3D1DOnP1P0_1 test.
 */
void ParaMEDMEMTest::testGauthier4()
{
  //
  const double sourceCoords[19*3]={0.5,0.5,0.1,0.5,0.5,1.2,0.5,0.5,1.6,0.5,0.5,1.8,0.5,0.5,2.43,0.5,0.5,2.55,0.5,0.5,4.1,0.5,0.5,4.4,0.5,0.5,4.9,0.5,0.5,5.1,0.5,0.5,7.6,0.5,0.5,7.7,0.5,0.5,8.2,0.5,0.5,8.4,0.5,0.5,8.6,0.5,0.5,8.8,0.5,0.5,9.2,0.5,0.5,9.6,0.5,0.5,11.5};
  const int sourceConn[18*2]={0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10,11,11,12,12,13,13,14,14,15,15,16,16,17,17,18};
  const double sourceVals[19]={0.49,2.8899999999999997,7.29,13.69,22.09,32.49,44.89,59.29,75.69,94.09, 114.49,136.89,161.29,187.69,216.09,246.49,278.89,313.29,349.69};
  const double targetCoords0[20*3]={0.,0.,0.,1.,0.,0.,0.,1.,0.,1.,1.,0.,0.,0.,1.,1.,0.,1.,0.,1.,1.,1.,1.,1.,0.,0.,2.,1.,0.,2.,0.,1.,2.,1.,1.,2.,0.,0.,3.,1.,0.,3.,0.,1.,3.,1.,1.,3.,0.,0.,4.,1.,0.,4.,0.,1.,4.,1.,1.,4.};
  const int targetConn0[8*4]={1,0,2,3,5,4,6,7,5,4,6,7,9,8,10,11,9,8,10,11,13,12,14,15,13,12,14,15,17,16,18,19};
  const double targetCoords1[28*3]={0.,0.,4.,1.,0.,4.,0.,1.,4.,1.,1.,4.,0.,0.,5.,1.,0.,5.,0.,1.,5.,1.,1.,5.,0.,0.,6.,1.,0.,6.,0.,1.,6.,1.,1.,6.,0.,0.,7.,1.,0.,7.,0.,1.,7.,1.,1.,7.,0.,0.,8.,1.,0.,8.,0.,1.,8.,1.,1.,8.,0.,0.,9.,1.,0.,9.,0.,1.,9.,1.,1.,9.,0.,0.,10.,1.,0.,10.,0.,1.,10.,1.,1.,10.};
  const int targetConn1[8*6]={1,0,2,3,5,4,6,7,5,4,6,7,9,8,10,11,9,8,10,11,13,12,14,15,13,12,14,15,17,16,18,19,17,16,18,19,21,20,22,23,21,20,22,23,25,24,26,27};
  //
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  //
  if(size!=3)
    return ;
  int nproc_source = 1;
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  //
  MEDCoupling::MEDCouplingUMesh *mesh=0;
  MEDCoupling::ParaMESH *paramesh=0;
  MEDCoupling::ParaFIELD* parafield=0;
  //
  MEDCoupling::CommInterface interface;
  //
  ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  //
  MPI_Barrier(MPI_COMM_WORLD);
  if(source_group->containsMyRank())
    {
      std::ostringstream stream; stream << "sourcemesh2D proc " << rank;
      mesh=MEDCouplingUMesh::New(stream.str().c_str(),1);
      mesh->allocateCells();
      for(int i=0;i<18;i++)
        mesh->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,sourceConn+2*i);
      mesh->finishInsertingCells();
      DataArrayDouble *myCoords=DataArrayDouble::New();
      myCoords->alloc(19,3);
      std::copy(sourceCoords,sourceCoords+19*3,myCoords->getPointer());
      mesh->setCoords(myCoords);
      myCoords->decrRef();
      paramesh=new ParaMESH(mesh,*source_group,"source mesh");
      MEDCoupling::ComponentTopology comptopo;
      parafield = new ParaFIELD(ON_NODES,NO_TIME,paramesh,comptopo);
      double *value=parafield->getField()->getArray()->getPointer();
      std::copy(sourceVals,sourceVals+19,value);
    }
  else
    {
      if(rank==1)
        {
          std::ostringstream stream; stream << "targetmesh2D proc " << rank-nproc_source;
          mesh=MEDCouplingUMesh::New(stream.str().c_str(),3);
          mesh->allocateCells();
          for(int i=0;i<4;i++)
            mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,targetConn0+8*i);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(20,3);
          std::copy(targetCoords0,targetCoords0+20*3,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          paramesh=new ParaMESH (mesh,*target_group,"target mesh");
          MEDCoupling::ComponentTopology comptopo;
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
        }
      else if(rank==2)
        {
          std::ostringstream stream; stream << "targetmesh2D proc " << rank-nproc_source;
          mesh=MEDCouplingUMesh::New(stream.str().c_str(),3);
          mesh->allocateCells();
          for(int i=0;i<6;i++)
            mesh->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8,targetConn1+8*i);
          mesh->finishInsertingCells();
          DataArrayDouble *myCoords=DataArrayDouble::New();
          myCoords->alloc(28,3);
          std::copy(targetCoords1,targetCoords1+28*3,myCoords->getPointer());
          mesh->setCoords(myCoords);
          myCoords->decrRef();
          paramesh=new ParaMESH (mesh,*target_group,"target mesh");
          MEDCoupling::ComponentTopology comptopo;
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
        }
    }
  //test 1 - primaire -> secondaire
  MEDCoupling::InterpKernelDEC dec(*source_group,*target_group);
  dec.setIntersectionType(INTERP_KERNEL::PointLocator);
  parafield->getField()->setNature(IntensiveMaximum);//very important
  if (source_group->containsMyRank())
    { 
      dec.setMethod("P1");
      dec.attachLocalField(parafield);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.sendData();
    }
  else
    {
      dec.setMethod("P0");
      dec.attachLocalField(parafield);
      dec.synchronize();
      dec.setForcedRenormalization(false);
      dec.recvData();
      const double *res(parafield->getField()->getArray()->getConstPointer());
      if(rank==1)
        {
          const double expected0[4]={0.49,7.956666666666667,27.29,0.};
          for(int i=0;i<4;i++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expected0[i],res[i],1e-13);
        }
      else
        {
          const double expected1[6]={59.95666666666667,94.09,0.,125.69,202.89,296.09};
          for(int i=0;i<6;i++)
            CPPUNIT_ASSERT_DOUBLES_EQUAL(expected1[i],res[i],1e-13);
        }
    }
  MPI_Barrier(MPI_COMM_WORLD);
  if (source_group->containsMyRank())
    {
      dec.recvData();
      const double expected2[19]={0.49,7.956666666666667,7.956666666666667,7.956666666666667,27.29,27.29,59.95666666666667,59.95666666666667,59.95666666666667,94.09,125.69,125.69,202.89,202.89,202.89,202.89,296.09,296.09,0.};
      const double *res(parafield->getField()->getArray()->getConstPointer());
      for(int i=0;i<19;i++)
        CPPUNIT_ASSERT_DOUBLES_EQUAL(expected2[i],res[i],1e-13);
    }
  else
    {
      dec.sendData();
    }
  delete parafield;
  mesh->decrRef();
  delete paramesh;
  delete self_group;
  delete target_group;
  delete source_group;
  //
  MPI_Barrier(MPI_COMM_WORLD);
}
