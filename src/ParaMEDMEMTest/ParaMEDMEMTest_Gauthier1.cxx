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

#include "ParaMEDMEMTest.hxx"
#include <cppunit/TestAssert.h>

#include <string>
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "DEC.hxx"
#include "InterpKernelDEC.hxx"
#include <set>
#include <time.h>
#include "ICoCoTrioField.hxx"
#include <iostream>
#include <assert.h>
#include <math.h>

using namespace std;
using namespace ParaMEDMEM;
using namespace ICoCo;

void afficheGauthier1( const TrioField&   field, const double *vals, int lgth)
{
  CPPUNIT_ASSERT_EQUAL(lgth,field._nb_elems);
  for (int ele=0;ele<field._nb_elems;ele++)
    CPPUNIT_ASSERT_DOUBLES_EQUAL(vals[ele],field._field[ele],1e-12);
}

void remplit_coordGauthier1(double* coords)
{
  double angle,epaisseur;
  angle=0*45*(asin(1.)/90);
  epaisseur=1e-0;
  coords[0*3+0]=0.;
  coords[0*3+1]=0.;
  coords[0*3+2]=0.;
  
  coords[1*3+0]=cos(angle);
  coords[1*3+1]=0.;
  coords[1*3+2]=sin(angle);
  
    
  coords[2*3+0]=-sin(angle);
  coords[2*3+1]=0.;
  coords[2*3+2]=cos(angle);
  
  for (int d=0;d<3;d++)
    coords[3*3+d]=coords[1*3+d]+ coords[2*3+d];
  
  for (int i=4;i<8;i++)
    {
      for (int d=0;d<3;d++)
        coords[i*3+d]=coords[(i-4)*3+d];
      coords[i*3+1]+=epaisseur;
    }

}

void init_quadGauthier1(TrioField& champ_quad,int is_master)
{
  
  champ_quad.setName("champ_quad");
  champ_quad._space_dim=3;
  champ_quad._mesh_dim=2;
  champ_quad._nodes_per_elem=4;
  champ_quad._itnumber=0;
  champ_quad._time1=0;
  champ_quad._time2=1;
  champ_quad._nb_field_components=1;

  if (is_master)
    {
      champ_quad._nbnodes=8;
      champ_quad._nb_elems=2;
      
      champ_quad._coords=new double[champ_quad._nbnodes*champ_quad._space_dim];
      //memcpy(afield._coords,sommets.addr(),champ_quad._nbnodes*champ_quad._space_dim*sizeof(double));
      
      remplit_coordGauthier1(champ_quad._coords);
  
  
      champ_quad._connectivity=new int[champ_quad._nb_elems*champ_quad._nodes_per_elem];
      champ_quad._connectivity[0*champ_quad._nodes_per_elem+0]=0;
      champ_quad._connectivity[0*champ_quad._nodes_per_elem+1]=1;
      champ_quad._connectivity[0*champ_quad._nodes_per_elem+2]=3;
      champ_quad._connectivity[0*champ_quad._nodes_per_elem+3]=2;
      champ_quad._connectivity[1*champ_quad._nodes_per_elem+0]=4;
      champ_quad._connectivity[1*champ_quad._nodes_per_elem+1]=5;
      champ_quad._connectivity[1*champ_quad._nodes_per_elem+2]=7;
      champ_quad._connectivity[1*champ_quad._nodes_per_elem+3]=6;
      
    }
  else
    {
      champ_quad._nbnodes=0;
      champ_quad._nb_elems=0;
      champ_quad._coords=new double[champ_quad._nbnodes*champ_quad._space_dim];
    
    }
  champ_quad._has_field_ownership=false;
  champ_quad._field=0;
  //champ_quad._field=new double[champ_quad._nb_elems];
  //  assert(champ_quad._nb_field_components==1);
}
void init_triangleGauthier1(TrioField& champ_triangle,int is_master)
{
   
  champ_triangle.setName("champ_triangle");
  champ_triangle._space_dim=3;
  champ_triangle._mesh_dim=2;
  champ_triangle._nodes_per_elem=3;
  champ_triangle._itnumber=0;
  champ_triangle._time1=0;
  champ_triangle._time2=1;
  champ_triangle._nb_field_components=1;

  if (is_master)
    {
      champ_triangle._nb_elems=4;
      champ_triangle._nbnodes=8;
    
      champ_triangle._coords=new double[champ_triangle._nbnodes*champ_triangle._space_dim];
      //memcpy(afield._coords,sommets.addr(),champ_triangle._nbnodes*champ_triangle._space_dim*sizeof(double));
      remplit_coordGauthier1(champ_triangle._coords);
      
      champ_triangle._connectivity=new int[champ_triangle._nb_elems*champ_triangle._nodes_per_elem];
      champ_triangle._connectivity[0*champ_triangle._nodes_per_elem+0]=0;
      champ_triangle._connectivity[0*champ_triangle._nodes_per_elem+1]=1;
      champ_triangle._connectivity[0*champ_triangle._nodes_per_elem+2]=2;
      champ_triangle._connectivity[1*champ_triangle._nodes_per_elem+0]=1;
      champ_triangle._connectivity[1*champ_triangle._nodes_per_elem+1]=2;
      champ_triangle._connectivity[1*champ_triangle._nodes_per_elem+2]=3;
      
      champ_triangle._connectivity[2*champ_triangle._nodes_per_elem+0]=4;
      champ_triangle._connectivity[2*champ_triangle._nodes_per_elem+1]=5;
      champ_triangle._connectivity[2*champ_triangle._nodes_per_elem+2]=7;
      champ_triangle._connectivity[3*champ_triangle._nodes_per_elem+0]=4;
      champ_triangle._connectivity[3*champ_triangle._nodes_per_elem+1]=6;
      champ_triangle._connectivity[3*champ_triangle._nodes_per_elem+2]=7;
    }
  else
    {
      champ_triangle._nb_elems=0;
      champ_triangle._nbnodes=0;
      champ_triangle._coords=new double[champ_triangle._nbnodes*champ_triangle._space_dim];
    
    }
  champ_triangle._has_field_ownership=false;
  // champ_triangle._field=new double[champ_triangle._nb_elems];
  champ_triangle._field=0;
  
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
        dec_emetteur.setOrientation(2);
        TrioField champ_emetteur, champ_recepteur;
   
        if (send==0)
          init_quadGauthier1(champ_emetteur,is_master);
        else
          init_triangleGauthier1(champ_emetteur,is_master);
        if (rec==0)
          init_triangleGauthier1(champ_recepteur,is_master);
        else
          init_quadGauthier1(champ_recepteur,is_master);
  
        if (cas=="emetteur") 
          {
            champ_emetteur._field=new double[champ_emetteur._nb_elems];
            for (int ele=0;ele<champ_emetteur._nb_elems;ele++)
              champ_emetteur._field[ele]=1;
      
            champ_emetteur._has_field_ownership=true;
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
              // if (champ_recepteur._field)
              //   delete [] champ_recepteur._field;
              champ_recepteur._field=0;
              // champ_recepteur._has_field_ownership=false;
  

  
              if (cas=="emetteur") 
                {
                  if (non_unif)
                    if(rank!=3)
                      champ_emetteur._field[0]=40;
                }
              //bool ok=false; // Is the time interval successfully solved ?
    
              // Loop on the time interval tries
              if(1) {
      

                if (cas=="emetteur")
                  dec_emetteur.attachLocalField((ICoCo::Field*) &champ_emetteur);
                else
                  dec_emetteur.attachLocalField((ICoCo::Field*) &champ_recepteur);


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
                      afficheGauthier1(champ_recepteur,expected[num_cas],expectedLgth[num_cas]);
                  }
                else
                  throw 0;
                MPI_Barrier(MPI_COMM_WORLD);
              }
              stop=true;
              num_cas++;
            }
          // destruction des champs, des DEC, et des tableaux associ�s
        }
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
        decu[0]=InterpKernelDEC(emetteur_group, recepteur_group);
        InterpKernelDEC& dec_emetteur=decu[0];

        //InterpKernelDEC dec_emetteur(emetteur_group, recepteur_group);
        dec_emetteur.setOrientation(2);
        TrioField champ_emetteur, champ_recepteur;
   
        if (send==0)
          init_quadGauthier1(champ_emetteur,is_master);
        else
          init_triangleGauthier1(champ_emetteur,is_master);
        if (rec==0)
          init_triangleGauthier1(champ_recepteur,is_master);
        else
          init_quadGauthier1(champ_recepteur,is_master);
  
        if (cas=="emetteur") 
          {
            champ_emetteur._field=new double[champ_emetteur._nb_elems];
            for (int ele=0;ele<champ_emetteur._nb_elems;ele++)
              champ_emetteur._field[ele]=1;
      
            champ_emetteur._has_field_ownership=true;
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
              // if (champ_recepteur._field)
              //   delete [] champ_recepteur._field;
              champ_recepteur._field=0;
              // champ_recepteur._has_field_ownership=false;
  

  
              if (cas=="emetteur") 
                {
                  if (non_unif)
                    if(rank!=3)
                      champ_emetteur._field[0]=40;
                }
              //bool ok=false; // Is the time interval successfully solved ?
    
              // Loop on the time interval tries
              if(1) {
      

                if (cas=="emetteur")
                  dec_emetteur.attachLocalField((ICoCo::Field*) &champ_emetteur);
                else
                  dec_emetteur.attachLocalField((ICoCo::Field*) &champ_recepteur);


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
                      afficheGauthier1(champ_recepteur,expected[num_cas],expectedLgth[num_cas]);
                  }
                else
                  throw 0;
                MPI_Barrier(MPI_COMM_WORLD);
              }
              stop=true;
              num_cas++;
            }
          // destruction des champs, des DEC, et des tableaux associ�s
        }
      }
}
