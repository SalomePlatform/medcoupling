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

using namespace std;
using namespace ParaMEDMEM;
using namespace ICoCo;

typedef enum {sync_and,sync_or} synctype;
void synchronize_bool(bool& stop, synctype s)
{
  int my_stop;
  int my_stop_temp = stop?1:0;
  if (s==sync_and)
    MPI_Allreduce(&my_stop_temp,&my_stop,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD);
  else if (s==sync_or)
    MPI_Allreduce(&my_stop_temp,&my_stop,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD);
  stop =(my_stop==1);
}

void synchronize_dt(double& dt)
{
  double dttemp=dt;
  MPI_Allreduce(&dttemp,&dt,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
}


void affiche( const TrioField&   field)
{
  cout <<field.getName()<<endl;
  for (int ele=0;ele<field._nb_elems;ele++)
    cout <<ele <<": "<<field._field[ele]<<endl;;
  
}

void remplit_coord(double* coords)
{
  coords[0*3+0]=0.;
  coords[0*3+1]=0.;
  coords[0*3+2]=0.;
  
  coords[1*3+0]=1.;
  coords[1*3+1]=0.;
  coords[1*3+2]=0.;
  
    
  coords[2*3+0]=0.;
  coords[2*3+1]=0.;
  coords[2*3+2]=1.;
  
  coords[3*3+0]=1.;
  coords[3*3+1]=0.;
  coords[3*3+2]=1.;
  
  for (int i=4;i<8;i++)
    {
      for (int d=0;d<3;d++)
        coords[i*3+d]=coords[(i-4)*3+d];
      coords[i*3+1]+=1e-5;
    }

}

void init_quad(TrioField& champ_quad)
{
  
  champ_quad.setName("champ_quad");
  champ_quad._space_dim=3;
  champ_quad._mesh_dim=2;
  champ_quad._nbnodes=8;
  champ_quad._nodes_per_elem=4;
  champ_quad._nb_elems=2;
  champ_quad._itnumber=0;
  champ_quad._time1=0;
  champ_quad._time2=1;
  champ_quad._nb_field_components=1;
  
  champ_quad._coords=new double[champ_quad._nbnodes*champ_quad._space_dim];
  //memcpy(afield._coords,sommets.addr(),champ_quad._nbnodes*champ_quad._space_dim*sizeof(double));
  
  remplit_coord(champ_quad._coords);
  
  
  champ_quad._connectivity=new int[champ_quad._nb_elems*champ_quad._nodes_per_elem];
  champ_quad._connectivity[0*champ_quad._nodes_per_elem+0]=0;
  champ_quad._connectivity[0*champ_quad._nodes_per_elem+1]=1;
  champ_quad._connectivity[0*champ_quad._nodes_per_elem+2]=3;
  champ_quad._connectivity[0*champ_quad._nodes_per_elem+3]=2;
  champ_quad._connectivity[1*champ_quad._nodes_per_elem+0]=4;
  champ_quad._connectivity[1*champ_quad._nodes_per_elem+1]=5;
  champ_quad._connectivity[1*champ_quad._nodes_per_elem+2]=7;
  champ_quad._connectivity[1*champ_quad._nodes_per_elem+3]=6;
  
  
  champ_quad._has_field_ownership=false;
  champ_quad._field=0;
  //champ_quad._field=new double[champ_quad._nb_elems];
  //  assert(champ_quad._nb_field_components==1);
}
void init_triangle(TrioField& champ_triangle)
{
   
  champ_triangle.setName("champ_triangle");
  champ_triangle._space_dim=3;
  champ_triangle._mesh_dim=2;
  champ_triangle._nbnodes=8;
  champ_triangle._nodes_per_elem=3;
  champ_triangle._nb_elems=4;
  champ_triangle._itnumber=0;
  champ_triangle._time1=0;
  champ_triangle._time2=1;
  champ_triangle._nb_field_components=1;
    
  champ_triangle._coords=new double[champ_triangle._nbnodes*champ_triangle._space_dim];
  //memcpy(afield._coords,sommets.addr(),champ_triangle._nbnodes*champ_triangle._space_dim*sizeof(double));
  remplit_coord(champ_triangle._coords);

  champ_triangle._connectivity=new int[champ_triangle._nb_elems*champ_triangle._nodes_per_elem];
  champ_triangle._connectivity[0*champ_triangle._nodes_per_elem+0]=0;
  champ_triangle._connectivity[0*champ_triangle._nodes_per_elem+1]=1;
  champ_triangle._connectivity[0*champ_triangle._nodes_per_elem+2]=2;
  champ_triangle._connectivity[1*champ_triangle._nodes_per_elem+0]=1;
  champ_triangle._connectivity[1*champ_triangle._nodes_per_elem+1]=3;
  champ_triangle._connectivity[1*champ_triangle._nodes_per_elem+2]=2;
  
  champ_triangle._connectivity[2*champ_triangle._nodes_per_elem+0]=4;
  champ_triangle._connectivity[2*champ_triangle._nodes_per_elem+1]=5;
  champ_triangle._connectivity[2*champ_triangle._nodes_per_elem+2]=7;
  champ_triangle._connectivity[3*champ_triangle._nodes_per_elem+0]=4;
  champ_triangle._connectivity[3*champ_triangle._nodes_per_elem+1]=7;
  champ_triangle._connectivity[3*champ_triangle._nodes_per_elem+2]=6;
  
  champ_triangle._has_field_ownership=false;
  // champ_triangle._field=new double[champ_triangle._nb_elems];
  champ_triangle._field=0;
}

void ParaMEDMEMTest::testICocoTrio1()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //the test is meant to run on five processors
  if (size !=2) return ;
  
  CommInterface comm;
  set<int> emetteur_ids;
  set<int> recepteur_ids;
  emetteur_ids.insert(0);
  recepteur_ids.insert(1);

  MPIProcessorGroup recepteur_group(comm,recepteur_ids);
  MPIProcessorGroup emetteur_group(comm,emetteur_ids);


  string cas;
  if (recepteur_group.containsMyRank())
    {
      cas="recepteur";
      
    }
  else
    cas="emetteur";

  InterpKernelDEC dec_emetteur(emetteur_group, recepteur_group);

  TrioField champ_emetteur, champ_recepteur;
   
  init_triangle(champ_emetteur);
  //init_triangle(champ_emetteur);
  init_quad(champ_recepteur);
  //init_emetteur(champ_recepteur);
  
  if (cas=="emetteur") 
    {
      champ_emetteur._field=new double[champ_emetteur._nb_elems];
      for (int ele=0;ele<champ_emetteur._nb_elems;ele++)
        champ_emetteur._field[ele]=1;
      
      champ_emetteur._has_field_ownership=true;
    }
  
  
  MPI_Barrier(MPI_COMM_WORLD);

  clock_t clock0= clock ();
  int compti=0;

  bool init=true; // first time step ??
  bool stop=false;
  //boucle sur les pas de quads
  while (!stop) {
  
    compti++;
    clock_t clocki= clock ();
    cout << compti << " CLOCK " << (clocki-clock0)*1.e-6 << endl; 
    for (int non_unif=0;non_unif<2;non_unif++)
      {
        // if (champ_recepteur._field)
        //   delete [] champ_recepteur._field;
        champ_recepteur._field=0;
        // champ_recepteur._has_field_ownership=false;
  

  
        if (cas=="emetteur") 
          if (non_unif)
            champ_emetteur._field[0]=40;
        //bool ok=false; // Is the time interval successfully solved ?
    
        // Loop on the time interval tries
        if(1)
          {
            if (cas=="emetteur")
              dec_emetteur.attachLocalField((ICoCo::Field*) &champ_emetteur);
            else
              dec_emetteur.attachLocalField((ICoCo::Field*) &champ_recepteur);
            
            dec_emetteur.setNature(ConservativeVolumic);
            
            if(init)
              dec_emetteur.synchronize();
            init=false;
            
            if (cas=="emetteur")
              {
                dec_emetteur.sendData();
                affiche(champ_emetteur);
              }
            else if (cas=="recepteur")
              {
                dec_emetteur.recvData();
                affiche(champ_recepteur);
              }
            else
              throw 0;
          }
        stop=true;
      }
  }
}
