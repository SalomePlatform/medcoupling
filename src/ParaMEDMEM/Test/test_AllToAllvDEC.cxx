//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#include <mpi.h>
#include <string>
#include <vector>
#include <map>
#include <iostream>

#include "MPIAccessDECTest.hxx"
#include <cppunit/TestAssert.h>

//#include "CommInterface.hxx"
//#include "ProcessorGroup.hxx"
//#include "MPIProcessorGroup.hxx"
#include "MPI_AccessDEC.hxx"

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES

using namespace std;
using namespace ParaMEDMEM;

void MPIAccessDECTest::test_AllToAllvDECSynchronousPointToPoint() {
  test_AllToAllvDEC( false ) ;
}
void MPIAccessDECTest::test_AllToAllvDECAsynchronousPointToPoint() {
  test_AllToAllvDEC( true ) ;
}

static void chksts( int sts , int myrank , ParaMEDMEM::MPI_Access mpi_access ) {
  char msgerr[MPI_MAX_ERROR_STRING] ;
  int lenerr ;
  if ( sts != MPI_SUCCESS ) {
    mpi_access.Error_String(sts, msgerr, &lenerr) ;
    cout << "test_AllToAllvDEC" << myrank << " lenerr " << lenerr << " "
         << msgerr << endl ;
    ostringstream strstream ;
    strstream << "==========================================================="
              << "test_AllToAllvDEC" << myrank << " KO"
              << "==========================================================="
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  return ;
}

void MPIAccessDECTest::test_AllToAllvDEC( bool Asynchronous ) {

  cout << "test_AllToAllvDEC" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 || size > 11 ) {
    ostringstream strstream ;
    strstream << "usage :" << endl
              << "mpirun -np <nbprocs> test_AllToAllvDEC" << endl
              << " (nbprocs >=2)" << endl
              << "test must be runned with more than 1 proc and less than 12 procs"
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

//  int Asynchronous = atoi(argv[1]);

  cout << "test_AllToAllvDEC" << myrank << endl ;

  ParaMEDMEM::CommInterface interface ;
  std::set<int> sourceprocs;
  std::set<int> targetprocs;
  int i ;
  for ( i = 0 ; i < size/2 ; i++ ) {
     sourceprocs.insert(i);
  }
  for ( i = size/2 ; i < size ; i++ ) {
     targetprocs.insert(i);
  }

  ParaMEDMEM::MPIProcessorGroup* sourcegroup = new ParaMEDMEM::MPIProcessorGroup(interface,sourceprocs) ;
  ParaMEDMEM::MPIProcessorGroup* targetgroup = new ParaMEDMEM::MPIProcessorGroup(interface,targetprocs) ;

  MPI_AccessDEC * MPIAccessDEC = new MPI_AccessDEC( *sourcegroup , *targetgroup ,
                                                    Asynchronous ) ;
  
  MPI_Access * mpi_access = MPIAccessDEC->MPIAccess() ;

#define maxreq 100
#define datamsglength 10

//  int sts ;
  int *sendcounts = new int[size] ;
  int *sdispls = new int[size] ;
  int *recvcounts = new int[size] ;
  int *rdispls = new int[size] ;
  for ( i = 0 ; i < size ; i++ ) {
     sendcounts[i] = datamsglength-i;
     sdispls[i] = i*datamsglength ;
     recvcounts[i] = datamsglength-myrank;
     rdispls[i] = i*datamsglength ;
  }
  int * recvbuf = new int[datamsglength*size] ;

  int ireq ;
  for ( ireq = 0 ; ireq < maxreq ; ireq++ ) {
    int * sendbuf = new int[datamsglength*size] ;
//    int * sendbuf = (int *) malloc( sizeof(int)*datamsglength*size) ;
    int j ;
    for ( j = 0 ; j < datamsglength*size ; j++ ) {
       sendbuf[j] = myrank*1000000 + ireq*1000 + j ;
       recvbuf[j] = -1 ;
    }

    MPIAccessDEC->AllToAllv( sendbuf, sendcounts , sdispls , MPI_INT ,
	                     recvbuf, recvcounts , rdispls , MPI_INT ) ;

//    cout << "test_AllToAllvDEC" << myrank << " recvbuf before CheckSent" ;
//    for ( i = 0 ; i < datamsglength*size ; i++ ) {
//       cout << " " << recvbuf[i] ;
//    }
//    cout << endl ;

//    cout << "test_AllToAllvDEC" << myrank << " sendbuf " << sendbuf << endl ;
//    MPIAccessDEC->CheckSent() ;

    int nRecvReq = mpi_access->RecvRequestIdsSize() ;
//    cout << "test_AllToAllvDEC" << myrank << " WaitAllRecv " << nRecvReq << " Requests" << endl ;
    int *ArrayOfRecvRequests = new int[nRecvReq] ;
    int nReq = mpi_access->RecvRequestIds( nRecvReq, ArrayOfRecvRequests ) ;
    mpi_access->WaitAll( nReq , ArrayOfRecvRequests ) ;
    mpi_access->DeleteRequests( nReq , ArrayOfRecvRequests ) ;
    delete [] ArrayOfRecvRequests ;

//    cout << "test_AllToAllvDEC" << myrank << " recvbuf" ;
//    for ( i = 0 ; i < datamsglength*size ; i++ ) {
//       cout << " " << recvbuf[i] ;
//    }
//    cout << endl ;
  }

//  cout << "test_AllToAllvDEC" << myrank << " final CheckSent" << endl ;
//  MPIAccessDEC->CheckSent() ;

  int nSendReq = mpi_access->SendRequestIdsSize() ;
  cout << "test_AllToAllvDEC" << myrank << " final SendRequestIds " << nSendReq << " SendRequests"
       << endl ;
  if ( nSendReq ) {
    int *ArrayOfSendRequests = new int[nSendReq] ;
    int nReq = mpi_access->SendRequestIds( nSendReq, ArrayOfSendRequests ) ;
    mpi_access->WaitAll( nReq , ArrayOfSendRequests ) ;
    delete [] ArrayOfSendRequests ;
  }

  int nRecvReq = mpi_access->RecvRequestIdsSize() ;
  if ( nRecvReq ) {
    ostringstream strstream ;
    strstream << "test_AllToAllvDEC" << myrank << " final RecvRequestIds " << nRecvReq
              << " RecvRequests # 0 Error" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    cout << "test_AllToAllvDEC" << myrank << " final RecvRequestIds " << nRecvReq
         << " RecvRequests = 0 OK" << endl ;
  }

  mpi_access->Barrier() ;

  delete sourcegroup ;
  delete targetgroup ;
  delete MPIAccessDEC ;
  delete [] sendcounts ;
  delete [] sdispls ;
  delete [] recvcounts ;
  delete [] rdispls ;
  delete [] recvbuf ;

//  MPI_Finalize();

  cout << "test_AllToAllvDEC" << myrank << " OK" << endl ;

  return ;
}




