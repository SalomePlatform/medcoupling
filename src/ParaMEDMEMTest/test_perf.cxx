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

#include <time.h>
#include <sys/times.h>
#include <sys/time.h>
#include "ParaMEDMEMTest.hxx"
#include <cppunit/TestAssert.h>

#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "Topology.hxx"
#include "DEC.hxx"
#include "MxN_Mapping.hxx"
#include "InterpKernelDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ComponentTopology.hxx"
#include "ICoCoMEDField.hxx"
#include "MEDLoader.hxx"
 
#include <string>
#include <cstring>

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES

#ifndef CLK_TCK 
#include <unistd.h>
#define CLK_TCK sysconf(_SC_CLK_TCK);
#endif 

using namespace std;
using namespace MEDCoupling;
 
void testInterpKernelDEC_2D(const string& filename1, const string& meshname1,
                            const string& filename2, const string& meshname2,
                            int nproc_source, double epsilon, bool tri, bool all);
void get_time( float *telps, float *tuser, float *tsys, float *tcpu );

int main(int argc, char *argv[])
{
  string filename1, filename2;
  string meshname1, meshname2;
  int nproc_source=1, rank;
  double epsilon=1.e-6;
  int count=0;
  bool tri=false;
  bool all=false;

  MPI_Init(&argc,&argv);

  for(int i=1;i<argc;i++){
    if( strcmp(argv[i],"-f1") == 0 ){
      filename1 = argv[++i];
      count++;
    }
    else if( strcmp(argv[i],"-f2") == 0 ){
      filename2 = argv[++i];
      count++;
    }
    else if( strcmp(argv[i],"-m1") == 0 ){
      meshname1 = argv[++i];
      count++;
    }
    else if( strcmp(argv[i],"-m2") == 0 ){
      meshname2 = argv[++i];
      count++;
    }
    else if( strcmp(argv[i],"-ns") == 0 ){
      nproc_source = atoi(argv[++i]);
    }
    else if( strcmp(argv[i],"-eps") == 0 ){
      epsilon = atof(argv[++i]);
    }
    else if( strcmp(argv[i],"-tri") == 0 ){
      tri = true;
    }
    else if( strcmp(argv[i],"-all") == 0 ){
      all = true;
    }
  }

  if( count != 4 ){
    cout << "usage test_perf -f1 filename1 -m1 meshname1 -f2 filename2 -m2 meshname2 (-ns nproc_source -eps epsilon -tri -all)" << endl;
    exit(0);
  }

  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  testInterpKernelDEC_2D(filename1,meshname1,filename2,meshname2,nproc_source,epsilon,tri,all);

  MPI_Finalize();
}

void testInterpKernelDEC_2D(const string& filename_xml1, const string& meshname1,
                            const string& filename_xml2, const string& meshname2,
                            int nproc_source, double epsilon, bool tri, bool all)
{
  float tcpu, tcpu_u, tcpu_s, telps;
  int size;
  int rank;
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
 
  set<int> self_procs;
  set<int> procs_source;
  set<int> procs_target;
  
  for (int i=0; i<nproc_source; i++)
    procs_source.insert(i);
  for (int i=nproc_source; i<size; i++)
    procs_target.insert(i);
  self_procs.insert(rank);
  
  MEDCoupling::CommInterface interface;
    
  MEDCoupling::ProcessorGroup* self_group = new MEDCoupling::MPIProcessorGroup(interface,self_procs);
  MEDCoupling::ProcessorGroup* target_group = new MEDCoupling::MPIProcessorGroup(interface,procs_target);
  MEDCoupling::ProcessorGroup* source_group = new MEDCoupling::MPIProcessorGroup(interface,procs_source);
  
  //loading the geometry for the source group

  MEDCoupling::InterpKernelDEC dec (*source_group,*target_group);
  if(tri)
    dec.setIntersectionType(INTERP_KERNEL::Triangulation);
  else
    dec.setIntersectionType(INTERP_KERNEL::Convex);

  MEDCoupling::MEDCouplingUMesh* mesh;
  MEDCoupling::ParaMESH* paramesh;
  MEDCoupling::ParaFIELD* parafield;
  ICoCo::MEDField* icocofield ;
  
  // To remove tmp files from disk
  ParaMEDMEMTest_TmpFilesRemover aRemover;
  
  MPI_Barrier(MPI_COMM_WORLD);
  if (source_group->containsMyRank()){
    string master = filename_xml1;
      
    ostringstream strstream;
    if( nproc_source == 1 )
      strstream <<master<<".med";
    else
      strstream <<master<<rank+1<<".med";

    ostringstream meshname ;
    if( nproc_source == 1 )
      meshname<< meshname1;
    else
      meshname<< meshname1<<"_"<< rank+1;
      
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    mesh=ReadUMeshFromFile(strstream.str().c_str(),meshname.str().c_str(),0);
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    if( rank == 0 )
      cout << "IO : Telapse = " << telps << " TuserCPU = " << tcpu_u << " TsysCPU = " << tcpu_s << " TCPU = " << tcpu << endl;
    mesh->incrRef();
    
    paramesh=new ParaMESH (mesh,*source_group,"source mesh");
    
    MEDCoupling::ComponentTopology comptopo;
    parafield = new ParaFIELD(ON_CELLS, NO_TIME, paramesh, comptopo);

    int nb_local=mesh->getNumberOfCells();
    double *value=parafield->getField()->getArray()->getPointer();
    for(int ielem=0; ielem<nb_local;ielem++)
      value[ielem]=1.0;
    
    icocofield=new ICoCo::MEDField(parafield->getField());
     
    dec.attachLocalField(icocofield);
  }
  
  //loading the geometry for the target group
  if (target_group->containsMyRank()){
    string master= filename_xml2;
    ostringstream strstream;
    if( (size-nproc_source) == 1 )
      strstream << master<<".med";
    else
      strstream << master<<(rank-nproc_source+1)<<".med";
    ostringstream meshname ;
    if( (size-nproc_source) == 1 )
      meshname<< meshname2;
    else
      meshname<< meshname2<<"_"<<rank-nproc_source+1;
      
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    mesh = ReadUMeshFromFile(strstream.str().c_str(),meshname.str().c_str(),0);
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    mesh->incrRef();

    paramesh=new ParaMESH (mesh,*target_group,"target mesh");
    MEDCoupling::ComponentTopology comptopo;
    parafield = new ParaFIELD(ON_CELLS,NO_TIME,paramesh, comptopo);

    int nb_local=mesh->getNumberOfCells();
    double *value=parafield->getField()->getArray()->getPointer();
    for(int ielem=0; ielem<nb_local;ielem++)
      value[ielem]=0.0;
    icocofield=new ICoCo::MEDField(parafield->getField());
      
    dec.attachLocalField(icocofield);
  }
    
  
  //attaching a DEC to the source group 
  double field_before_int;
  double field_after_int;
  
  if (source_group->containsMyRank()){ 
    field_before_int = parafield->getVolumeIntegral(0,true);
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    dec.synchronize();
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    if( rank == 0 )
      cout << "SYNCHRONIZE : Telapse = " << telps << " TuserCPU = " << tcpu_u << " TsysCPU = " << tcpu_s << " TCPU = " << tcpu << endl;
    cout<<"DEC usage"<<endl;
    dec.setForcedRenormalization(false);
    if(all)
      dec.setAllToAllMethod(PointToPoint);

    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    dec.sendData();
    
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    if( rank == 0 )
      cout << "SEND DATA : Telapse = " << telps << " TuserCPU = " << tcpu_u << " TsysCPU = " << tcpu_s << " TCPU = " << tcpu << endl;
    dec.recvData();
     
    field_after_int = parafield->getVolumeIntegral(0,true);
//    CPPUNIT_ASSERT_DOUBLES_EQUAL(field_before_int, field_after_int, epsilon);
      
  }
  
  //attaching a DEC to the target group
  if (target_group->containsMyRank()){
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    dec.synchronize();
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    dec.setForcedRenormalization(false);
    if(all)
      dec.setAllToAllMethod(PointToPoint);

    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    dec.recvData();
    get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
    dec.sendData();
  }
  
  get_time( &telps, &tcpu_u, &tcpu_s, &tcpu );
  if( rank == 0 )
    cout << "RECV DATA : Telapse = " << telps << " TuserCPU = " << tcpu_u << " TsysCPU = " << tcpu_s << " TCPU = " << tcpu << endl;

  delete source_group;
  delete target_group;
  delete self_group;
  delete paramesh;
  delete parafield;
  mesh->decrRef() ;
  delete icocofield;

  MPI_Barrier(MPI_COMM_WORLD);
  cout << "end of InterpKernelDEC_2D test"<<endl;
}

void get_time( float *telps, float *tuser, float *tsys, float *tcpu )
{

  /* Variables declaration */
  static time_t zsec = 0;
  static long zusec = 0;
  time_t nsec;
  long nusec;
  static clock_t zclock = 0;
  clock_t nclock;
  static clock_t zuser = 0;
  static clock_t zsys = 0;
  clock_t nuser, nsys;

  struct timeval tp;
  struct timezone tzp;
  struct tms local;

  MPI_Barrier(MPI_COMM_WORLD);

  /* Elapsed time reading */

  gettimeofday(&tp,&tzp);
  nsec = tp.tv_sec;
  nusec = tp.tv_usec;
  *telps = (float)(nsec-zsec) + (float)(nusec-zusec)/(float)CLOCKS_PER_SEC;
  
  zsec = nsec;
  zusec = nusec;

  /* User and system CPU time reading */

  times(&local);
  nuser = local.tms_utime;
  nsys = local.tms_stime;
  *tuser = (float)(nuser-zuser) / (float)CLK_TCK;
  *tsys = (float)(nsys-zsys) / (float)CLK_TCK;

  zuser = nuser;
  zsys = nsys;

  /* CPU time reading */

  nclock = clock();
  *tcpu = (float)(nclock-zclock) / (float)CLOCKS_PER_SEC;
  zclock = nclock;

}


