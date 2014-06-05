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
#include "DEC.hxx"
#include "InterpKernelDEC.hxx"
#include "ICoCoTrioField.hxx"
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
using namespace ParaMEDMEM;
using namespace ICoCo;

void afficheGauthier1(const ParaFIELD& field, const double *vals, int lgth)
{
  const DataArrayDouble *valsOfField(field.getField()->getArray());
  CPPUNIT_ASSERT_EQUAL(lgth,valsOfField->getNumberOfTuples());
  for (int ele=0;ele<valsOfField->getNumberOfTuples();ele++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vals[ele],valsOfField->getIJ(ele,0),1e-12);
}

MEDCouplingUMesh *init_quadGauthier1(int is_master)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m(MEDCouplingUMesh::New("champ_quad",2));
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coo(DataArrayDouble::New());
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
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> m(MEDCouplingUMesh::New("champ_triangle",2));
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coo(DataArrayDouble::New());
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
        ParaMEDMEM::ParaFIELD *champ_emetteur(0),*champ_recepteur(0);
        ParaMEDMEM::ParaMESH *paramesh(0);
        MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mesh;
        dec_emetteur.setOrientation(2);
        if (send==0)
          {
            mesh=init_quadGauthier1(is_master);
          }
        else
          {
            mesh=init_triangleGauthier1(is_master);
          }
        paramesh=new ParaMEDMEM::ParaMESH(mesh,recepteur_group.containsMyRank()?recepteur_group:emetteur_group,"emetteur mesh");
        ParaMEDMEM::ComponentTopology comptopo;
        champ_emetteur=new ParaMEDMEM::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
        champ_emetteur->getField()->setNature(ConservativeVolumic);
        champ_emetteur->setOwnSupport(true);
        if (rec==0)
          {
            mesh=init_triangleGauthier1(is_master);
          }
        else
          {
            mesh=init_quadGauthier1(is_master);
          }
        paramesh=new ParaMEDMEM::ParaMESH(mesh,recepteur_group.containsMyRank()?recepteur_group:emetteur_group,"recepteur mesh");
        champ_recepteur=new ParaMEDMEM::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
        champ_recepteur->getField()->setNature(ConservativeVolumic);
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
  const char save_vit_in_2[]="VITESSE_P1_OUT\n1\n2\n3\n63\n3\n80\n0\n 0 1 2\n 3 4 5\n 6 7 8\n 9 10 11\n 12 13 14\n 15 16 17\n 18 19 20\n 21 22 23\n 24 25 26\n 27 28 29\n 30 2 1\n 31 5 4\n 32 8 7\n 33 11 10\n 34 14 13\n 35 17 16\n 36 20 19\n 37 23 22\n 38 26 25\n 39 29 28\n 30 40 2\n 31 41 5\n 32 42 8\n 33 43 11\n 34 44 14\n 35 45 17\n 36 46 20\n 37 47 23\n 38 48 26\n 39 49 29\n 31 2 40\n 32 5 41\n 33 8 42\n 34 11 43\n 35 14 44\n 36 17 45\n 37 20 46\n 38 23 47\n 39 26 48\n 50 29 49\n 3 2 4\n 6 5 7\n 9 8 10\n 12 11 13\n 15 14 16\n 18 17 19\n 21 20 22\n 24 23 25\n 27 26 28\n 51 29 52\n 31 4 2\n 32 7 5\n 33 10 8\n 34 13 11\n 35 16 14\n 36 19 17\n 37 22 20\n 38 25 23\n 39 28 26\n 50 52 29\n 0 2 53\n 3 5 54\n 6 8 55\n 9 11 56\n 12 14 57\n 15 17 58\n 18 20 59\n 21 23 60\n 24 26 61\n 27 29 62\n 3 53 2\n 6 54 5\n 9 55 8\n 12 56 11\n 15 57 14\n 18 58 17\n 21 59 20\n 24 60 23\n 27 61 26\n 51 62 29\n 0 0 0\n 0.5 0 0\n 0.5 0.05 0\n 0 0.1 0\n 0.5 0.1 0\n 0.5 0.15 0\n 0 0.2 0\n 0.5 0.2 0\n 0.5 0.25 0\n 0 0.3 0\n 0.5 0.3 0\n 0.5 0.35 0\n 0 0.4 0\n 0.5 0.4 0\n 0.5 0.45 0\n 0 0.5 0\n 0.5 0.5 0\n 0.5 0.55 0\n 0 0.6 0\n 0.5 0.6 0\n 0.5 0.65 0\n 0 0.7 0\n 0.5 0.7 0\n 0.5 0.75 0\n 0 0.8 0\n 0.5 0.8 0\n 0.5 0.85 0\n 0 0.9 0\n 0.5 0.9 0\n 0.5 0.95 0\n 1 0 0\n 1 0.1 0\n 1 0.2 0\n 1 0.3 0\n 1 0.4 0\n 1 0.5 0\n 1 0.6 0\n 1 0.7 0\n 1 0.8 0\n 1 0.9 0\n 1 0.05 0\n 1 0.15 0\n 1 0.25 0\n 1 0.35 0\n 1 0.45 0\n 1 0.55 0\n 1 0.65 0\n 1 0.75 0\n 1 0.85 0\n 1 0.95 0\n 1 1 0\n 0 1 0\n 0.5 1 0\n 0 0.05 0\n 0 0.15 0\n 0 0.25 0\n 0 0.35 0\n 0 0.45 0\n 0 0.55 0\n 0 0.65 0\n 0 0.75 0\n 0 0.85 0\n 0 0.95 0\n2.9268\n3.1707\n3\n1\n 0 0 0\n 0 0 0\n 0 0 0.05\n 0 0 0.1\n 0 0 0.1\n 0 0 0.15\n 0 0 0.2\n 0 0 0.2\n 0 0 0.25\n 0 0 0.3\n 0 0 0.3\n 0 0 0.35\n 0 0 0.4\n 0 0 0.4\n 0 0 0.45\n 0 0 0.5\n 0 0 0.5\n 0 0 0.55\n 0 0 0.6\n 0 0 0.6\n 0 0 0.65\n 0 0 0.7\n 0 0 0.7\n 0 0 0.75\n 0 0 0.8\n 0 0 0.8\n 0 0 0.85\n 0 0 0.9\n 0 0 0.9\n 0 0 0.95\n 0 0 0\n 0 0 0.1\n 0 0 0.2\n 0 0 0.3\n 0 0 0.4\n 0 0 0.5\n 0 0 0.6\n 0 0 0.7\n 0 0 0.8\n 0 0 0.9\n 0 0 0.05\n 0 0 0.15\n 0 0 0.25\n 0 0 0.35\n 0 0 0.45\n 0 0 0.55\n 0 0 0.65\n 0 0 0.75\n 0 0 0.85\n 0 0 0.95\n 0 0 1\n 0 0 1\n 0 0 1\n 0 0 0.05\n 0 0 0.15\n 0 0 0.25\n 0 0 0.35\n 0 0 0.45\n 0 0 0.55\n 0 0 0.65\n 0 0 0.75\n 0 0 0.85\n 0 0 0.95\n1\n";

  const char save_vit_out_0_2[]="vitesse_in_chaude\n0\n2\n3\n22\n4\n10\n-1081737852\n 0 1 3 2\n 2 3 5 4\n 4 5 7 6\n 6 7 9 8\n 8 9 11 10\n 10 11 13 12\n 12 13 15 14\n 14 15 17 16\n 16 17 19 18\n 18 19 21 20\n 0 0 0\n 1 0 0\n 0 0.1 0\n 1 0.1 0\n 0 0.2 0\n 1 0.2 0\n 0 0.3 0\n 1 0.3 0\n 0 0.4 0\n 1 0.4 0\n 0 0.5 0\n 1 0.5 0\n 0 0.6 0\n 1 0.6 0\n 0 0.7 0\n 1 0.7 0\n 0 0.8 0\n 1 0.8 0\n 0 0.9 0\n 1 0.9 0\n 0 1 0\n 1 1 0\n2.9268\n3.1707\n3\n1\n 0 0 0.05\n 0 0 0.15\n 0 0 0.25\n 0 0 0.35\n 0 0 0.45\n 0 0 0.55\n 0 0 0.65\n 0 0 0.75\n 0 0 0.85\n 0 0 0.95\n0\n";
  const char save_vit_out_1_2[]="vitesse_in_chaude\n1\n2\n3\n22\n4\n10\n-1081737852\n 0 1 3 2\n 2 3 5 4\n 4 5 7 6\n 6 7 9 8\n 8 9 11 10\n 10 11 13 12\n 12 13 15 14\n 14 15 17 16\n 16 17 19 18\n 18 19 21 20\n 0 0 0\n 1 0 0\n 0 0.1 0\n 1 0.1 0\n 0 0.2 0\n 1 0.2 0\n 0 0.3 0\n 1 0.3 0\n 0 0.4 0\n 1 0.4 0\n 0 0.5 0\n 1 0.5 0\n 0 0.6 0\n 1 0.6 0\n 0 0.7 0\n 1 0.7 0\n 0 0.8 0\n 1 0.8 0\n 0 0.9 0\n 1 0.9 0\n 0 1 0\n 1 1 0\n2.9268\n3.1707\n3\n1\n 0 0 0.029375\n 0 0 0.029375\n 0 0 0.1\n 0 0 0.1\n 0 0 0.2\n 0 0 0.2\n 0 0 0.3\n 0 0 0.3\n 0 0 0.4\n 0 0 0.4\n 0 0 0.5\n 0 0 0.5\n 0 0 0.6\n 0 0 0.6\n 0 0 0.7\n 0 0 0.7\n 0 0 0.8\n 0 0 0.8\n 0 0 0.9\n 0 0 0.9\n 0 0 0.970625\n 0 0 0.970625\n0\n";

  const char *save_vit_outs[2]={save_vit_out_1_2,save_vit_out_0_2};

  const char save_vit_out_1_0[]="vitesse_in_chaude\n1\n2\n3\n22\n4\n10\n-1081737852\n 0 1 3 2\n 2 3 5 4\n 4 5 7 6\n 6 7 9 8\n 8 9 11 10\n 10 11 13 12\n 12 13 15 14\n 14 15 17 16\n 16 17 19 18\n 18 19 21 20\n 0 0 0\n 1 0 0\n 0 0.1 0\n 1 0.1 0\n 0 0.2 0\n 1 0.2 0\n 0 0.3 0\n 1 0.3 0\n 0 0.4 0\n 1 0.4 0\n 0 0.5 0\n 1 0.5 0\n 0 0.6 0\n 1 0.6 0\n 0 0.7 0\n 1 0.7 0\n 0 0.8 0\n 1 0.8 0\n 0 0.9 0\n 1 0.9 0\n 0 1 0\n 1 1 0\n2.9268\n3.1707\n3\n1\n 0 0 0.029375\n 0 0 0.029375\n 0 0 0.1\n 0 0 0.1\n 0 0 0.2\n 0 0 0.2\n 0 0 0.3\n 0 0 0.3\n 0 0 0.4\n 0 0 0.4\n 0 0 0.5\n 0 0 0.5\n 0 0 0.6\n 0 0 0.6\n 0 0 0.7\n 0 0 0.7\n 0 0 0.8\n 0 0 0.8\n 0 0 0.9\n 0 0 0.9\n 0 0 0.970625\n 0 0 0.970625\n0\n";
  
  const char save_vit_in[]="VITESSE_P1_OUT\n1\n2\n3\n63\n3\n80\n0\n 0 1 2\n 3 4 5\n 6 7 8\n 9 10 11\n 12 13 14\n 15 16 17\n 18 19 20\n 21 22 23\n 24 25 26\n 27 28 29\n 30 2 1\n 31 5 4\n 32 8 7\n 33 11 10\n 34 14 13\n 35 17 16\n 36 20 19\n 37 23 22\n 38 26 25\n 39 29 28\n 30 40 2\n 31 41 5\n 32 42 8\n 33 43 11\n 34 44 14\n 35 45 17\n 36 46 20\n 37 47 23\n 38 48 26\n 39 49 29\n 31 2 40\n 32 5 41\n 33 8 42\n 34 11 43\n 35 14 44\n 36 17 45\n 37 20 46\n 38 23 47\n 39 26 48\n 50 29 49\n 3 2 4\n 6 5 7\n 9 8 10\n 12 11 13\n 15 14 16\n 18 17 19\n 21 20 22\n 24 23 25\n 27 26 28\n 51 29 52\n 31 4 2\n 32 7 5\n 33 10 8\n 34 13 11\n 35 16 14\n 36 19 17\n 37 22 20\n 38 25 23\n 39 28 26\n 50 52 29\n 0 2 53\n 3 5 54\n 6 8 55\n 9 11 56\n 12 14 57\n 15 17 58\n 18 20 59\n 21 23 60\n 24 26 61\n 27 29 62\n 3 53 2\n 6 54 5\n 9 55 8\n 12 56 11\n 15 57 14\n 18 58 17\n 21 59 20\n 24 60 23\n 27 61 26\n 51 62 29\n 0 0 0\n 0.5 0 0\n 0.5 0.05 0\n 0 0.1 0\n 0.5 0.1 0\n 0.5 0.15 0\n 0 0.2 0\n 0.5 0.2 0\n 0.5 0.25 0\n 0 0.3 0\n 0.5 0.3 0\n 0.5 0.35 0\n 0 0.4 0\n 0.5 0.4 0\n 0.5 0.45 0\n 0 0.5 0\n 0.5 0.5 0\n 0.5 0.55 0\n 0 0.6 0\n 0.5 0.6 0\n 0.5 0.65 0\n 0 0.7 0\n 0.5 0.7 0\n 0.5 0.75 0\n 0 0.8 0\n 0.5 0.8 0\n 0.5 0.85 0\n 0 0.9 0\n 0.5 0.9 0\n 0.5 0.95 0\n 1 0 0\n 1 0.1 0\n 1 0.2 0\n 1 0.3 0\n 1 0.4 0\n 1 0.5 0\n 1 0.6 0\n 1 0.7 0\n 1 0.8 0\n 1 0.9 0\n 1 0.05 0\n 1 0.15 0\n 1 0.25 0\n 1 0.35 0\n 1 0.45 0\n 1 0.55 0\n 1 0.65 0\n 1 0.75 0\n 1 0.85 0\n 1 0.95 0\n 1 1 0\n 0 1 0\n 0.5 1 0\n 0 0.05 0\n 0 0.15 0\n 0 0.25 0\n 0 0.35 0\n 0 0.45 0\n 0 0.55 0\n 0 0.65 0\n 0 0.75 0\n 0 0.85 0\n 0 0.95 0\n2.9268\n3.1707\n3\n1\n 0 0 0\n 0 0 0\n 0 0 0.05\n 0 0 0.1\n 0 0 0.1\n 0 0 0.15\n 0 0 0.2\n 0 0 0.2\n 0 0 0.25\n 0 0 0.3\n 0 0 0.3\n 0 0 0.35\n 0 0 0.4\n 0 0 0.4\n 0 0 0.45\n 0 0 0.5\n 0 0 0.5\n 0 0 0.55\n 0 0 0.6\n 0 0 0.6\n 0 0 0.65\n 0 0 0.7\n 0 0 0.7\n 0 0 0.75\n 0 0 0.8\n 0 0 0.8\n 0 0 0.85\n 0 0 0.9\n 0 0 0.9\n 0 0 0.95\n 0 0 0\n 0 0 0.1\n 0 0 0.2\n 0 0 0.3\n 0 0 0.4\n 0 0 0.5\n 0 0 0.6\n 0 0 0.7\n 0 0 0.8\n 0 0 0.9\n 0 0 0.05\n 0 0 0.15\n 0 0 0.25\n 0 0 0.35\n 0 0 0.45\n 0 0 0.55\n 0 0 0.65\n 0 0 0.75\n 0 0 0.85\n 0 0 0.95\n 0 0 1\n 0 0 1\n 0 0 1\n 0 0 0.05\n 0 0 0.15\n 0 0 0.25\n 0 0 0.35\n 0 0 0.45\n 0 0 0.55\n 0 0 0.65\n 0 0 0.75\n 0 0 0.85\n 0 0 0.95\n1\n";

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

      TrioField vitesse;
      InterpKernelDEC dec_vit_in_chaude(entree_chaude_group, Genepi_group);

      if ( entree_chaude_group.containsMyRank())
        {
          istringstream save_vit(save_vit_in);
          vitesse.restore(save_vit);
        }
      else
        {
          istringstream save_vit(save_vit_out_1_0);
          vitesse.restore(save_vit);
          vitesse._has_field_ownership=false;
      
          if (vitesse._field)
            {
              delete [] vitesse._field;
              // cette ligne est super importante sinon c'est tout faux !!!!!!!
              vitesse._field=0;
            }
          // pour tester P1->P0
          vitesse._type=type;  
        }
  
      if (vitesse._type==1)
        dec_vit_in_chaude.setMethod("P1");
  
  

      dec_vit_in_chaude.attachLocalField((ICoCo::Field*) &vitesse);
      
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
      if (entree_chaude_group.containsMyRank() )
        {
          if (1)
            {
              ostringstream save_vit(save_vit_in_2);
              vitesse.save(save_vit);
            }
        }
      else
        {
      
          double pmin=1e38, pmax=-1e38;
      
          for(int i=0;i<vitesse.nb_values()*vitesse._nb_field_components;i++)
            {
              double p=*(vitesse._field+i);
              if (p<pmin) pmin=p;
              if (p>pmax) pmax=p;
            }
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected1[type],pmin,1e-12);
          CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected2[type],pmax,1e-12);
      
          ostringstream save_vit(save_vit_outs[type]);
          vitesse.save(save_vit);

          for(int i=0;i<vitesse.nb_values();i++)
            {
              for(int c=0;c<vitesse._nb_field_components;c++)
                {
                  double p=vitesse._field[i*vitesse._nb_field_components+c];
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(valuesExpected3[type][i*vitesse._nb_field_components+c],p,1e-12);
                }
            }
      
        }
    }
}

/*!
 * Non regression test testing copy constructor of InterpKernelDEC. 
 */
void ParaMEDMEMTest::testGauthier3()
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
        std::vector<InterpKernelDEC> decu(1);
        decu[0]=InterpKernelDEC(emetteur_group,recepteur_group);
        InterpKernelDEC& dec_emetteur=decu[0];
        ParaMEDMEM::ParaFIELD *champ_emetteur(0),*champ_recepteur(0);
        ParaMEDMEM::ParaMESH *paramesh(0);
        MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> mesh;
        dec_emetteur.setOrientation(2);
        if (send==0)
          {
            mesh=init_quadGauthier1(is_master);
          }
        else
          {
            mesh=init_triangleGauthier1(is_master);
          }
        paramesh=new ParaMEDMEM::ParaMESH(mesh,recepteur_group.containsMyRank()?recepteur_group:emetteur_group,"emetteur mesh");
        ParaMEDMEM::ComponentTopology comptopo;
        champ_emetteur=new ParaMEDMEM::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
        champ_emetteur->getField()->setNature(ConservativeVolumic);
        champ_emetteur->setOwnSupport(true);
        if (rec==0)
          {
            mesh=init_triangleGauthier1(is_master);
          }
        else
          {
            mesh=init_quadGauthier1(is_master);
          }
        paramesh=new ParaMEDMEM::ParaMESH(mesh,recepteur_group.containsMyRank()?recepteur_group:emetteur_group,"recepteur mesh");
        champ_recepteur=new ParaMEDMEM::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
        champ_recepteur->getField()->setNature(ConservativeVolumic);
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
  ParaMEDMEM::MEDCouplingUMesh *mesh=0;
  ParaMEDMEM::ParaMESH *paramesh=0;
  ParaMEDMEM::ParaFIELD* parafield=0;
  //
  ParaMEDMEM::CommInterface interface;
  //
  ProcessorGroup* self_group = new ParaMEDMEM::MPIProcessorGroup(interface,self_procs);
  ProcessorGroup* target_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_target);
  ProcessorGroup* source_group = new ParaMEDMEM::MPIProcessorGroup(interface,procs_source);
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
      ParaMEDMEM::ComponentTopology comptopo;
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
          ParaMEDMEM::ComponentTopology comptopo;
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
          ParaMEDMEM::ComponentTopology comptopo;
          parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);
        }
    }
  //test 1 - primaire -> secondaire
  ParaMEDMEM::InterpKernelDEC dec(*source_group,*target_group);
  dec.setIntersectionType(INTERP_KERNEL::PointLocator);
  parafield->getField()->setNature(ConservativeVolumic);//very important
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
