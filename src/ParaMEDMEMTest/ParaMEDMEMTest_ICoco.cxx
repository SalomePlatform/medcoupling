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
#include "ComponentTopology.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "InterpKernelDEC.hxx"

#include "MEDCouplingUMesh.hxx"

#include <set>
#include <string>
#include <time.h>
#include <iostream>
#include <assert.h>

using namespace std;
using namespace MEDCoupling;
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


void affiche(const ParaFIELD& field)
{
  cout <<field.getField()->getName()<<endl;
  const double *vals(field.getField()->getArray()->begin());
  for(int ele=0;ele<field.getField()->getNumberOfTuples();ele++)
    cout << ele <<": "<< vals[ele] << endl;
}

MEDCouplingUMesh *init_quad()
{
  MCAuto<MEDCouplingUMesh> m(MEDCouplingUMesh::New("champ_quad",2));
  MCAuto<DataArrayDouble> coo(DataArrayDouble::New());
  const double dataCoo[24]={0.,0.,0.,1.,0.,0.,0.,0.,1.,1.,0.,1.,0.,1e-05,0.,1.,1e-05,0.,0.,1e-05,1.,1.,1e-05,1.};
  coo->alloc(8,3);
  std::copy(dataCoo,dataCoo+24,coo->getPointer());
  const int conn[8]={0,1,3,2,4,5,7,6};
  m->allocateCells(2);
  m->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn);
  m->insertNextCell(INTERP_KERNEL::NORM_QUAD4,4,conn+4);
  m->setCoords(coo);
  return m.retn();
}

MEDCouplingUMesh *init_triangle()
{
  MCAuto<MEDCouplingUMesh> m(MEDCouplingUMesh::New("champ_triangle",2));
  MCAuto<DataArrayDouble> coo(DataArrayDouble::New());
  const double dataCoo[24]={0.,0.,0.,1.,0.,0.,0.,0.,1.,1.,0.,1.,0.,1e-05,0.,1.,1e-05,0.,0.,1e-05,1.,1.,1e-05,1.};
  coo->alloc(8,3);
  std::copy(dataCoo,dataCoo+24,coo->getPointer());
  const int conn[12]={0,1,2,1,2,3,4,5,7,4,6,7};
  m->allocateCells(4);
  for(int i=0;i<4;i++)
    m->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,conn+3*i);
  m->setCoords(coo);
  return m.retn();
}

void ParaMEDMEMTest::testICoco1()
{
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  //the test is meant to run on 2 processors
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
    cas="recepteur";
  else
    cas="emetteur";

  InterpKernelDEC dec_emetteur(emetteur_group,recepteur_group);
  dec_emetteur.setOrientation(2);
  MEDCoupling::ParaFIELD *champ_emetteur(0),*champ_recepteur(0);
  MEDCoupling::ParaMESH *paramesh(0);
  if (cas=="emetteur") 
    {
      MCAuto<MEDCoupling::MEDCouplingUMesh> mesh_emetteur(init_triangle());
      paramesh=new MEDCoupling::ParaMESH(mesh_emetteur,emetteur_group,"emetteur mesh");
      MEDCoupling::ComponentTopology comptopo;
      champ_emetteur=new MEDCoupling::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
      champ_emetteur->getField()->setNature(IntensiveMaximum);
      champ_emetteur->setOwnSupport(true);
      champ_emetteur->getField()->getArray()->fillWithValue(1.);
    }
  else
    {
      MCAuto<MEDCoupling::MEDCouplingUMesh> mesh_recepteur(init_quad());
      paramesh=new MEDCoupling::ParaMESH(mesh_recepteur,recepteur_group,"recepteur mesh");
      MEDCoupling::ComponentTopology comptopo;
      champ_recepteur=new MEDCoupling::ParaFIELD(ON_CELLS,ONE_TIME,paramesh,comptopo);
      champ_recepteur->getField()->setNature(IntensiveMaximum);
      champ_recepteur->setOwnSupport(true);
    }
  
  
  MPI_Barrier(MPI_COMM_WORLD);

  clock_t clock0(clock());
  int compti=0;

  bool init(true),stop(false);
  //boucle sur les pas de quads
  while(!stop)
    {
      compti++;
      clock_t clocki= clock ();
      cout << compti << " CLOCK " << (clocki-clock0)*1.e-6 << endl; 
      for (int non_unif=0;non_unif<2;non_unif++)
        {
          if (cas=="emetteur") 
            if (non_unif)
              champ_emetteur->getField()->getArray()->setIJ(0,0,40.);
          //bool ok=false; // Is the time interval successfully solved ?
    
          // Loop on the time interval tries
          if (cas=="emetteur")
            dec_emetteur.attachLocalField(champ_emetteur);
          else
            dec_emetteur.attachLocalField(champ_recepteur);
            
          if(init)
            dec_emetteur.synchronize();
          init=false;
            
          if (cas=="emetteur")
            {
              dec_emetteur.sendData();
              affiche(*champ_emetteur);
            }
          else if (cas=="recepteur")
            {
              dec_emetteur.recvData();
              affiche(*champ_recepteur);
            }
          else
            throw 0;
        }
      stop=true;
    }
  delete champ_recepteur;
  delete champ_emetteur;
}
