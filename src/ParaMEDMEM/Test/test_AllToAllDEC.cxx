#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <mpi.h>

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

void MPIAccessDECTest::test_AllToAllDECSynchronousPointToPoint() {
  test_AllToAllDEC( false ) ;
}
void MPIAccessDECTest::test_AllToAllDECAsynchronousPointToPoint() {
  test_AllToAllDEC( true ) ;
}

static void chksts( int sts , int myrank , ParaMEDMEM::MPI_Access mpi_access ) {
  char msgerr[MPI_MAX_ERROR_STRING] ;
  int lenerr ;
  if ( sts != MPI_SUCCESS ) {
    mpi_access.Error_String(sts, msgerr, &lenerr) ;
    cout << "test" << myrank << " lenerr " << lenerr << " "
         << msgerr << endl ;
    ostringstream strstream ;
    strstream << "===========================================================" << endl
              << "test_AllToAllDEC" << myrank << " KO" << endl
              << "==========================================================="
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  return ;
}

void MPIAccessDECTest::test_AllToAllDEC( bool Asynchronous ) {

  cout << "test_AllToAllDEC" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 || size > 11 ) {
    ostringstream strstream ;
    strstream << "usage :" << endl
              << "mpirun -np <nbprocs> test_AllToAllDEC" << endl
              << " (nbprocs >=2)" << endl
              << "test must be runned with more than 1 proc and less than 12 procs"
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

//  int Asynchronous = atoi(argv[1]);

  cout << "test_AllToAllDEC" << myrank << endl ;

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

//  MPI_AccessDEC * MPIAccessDEC = new MPI_AccessDEC( *sourcegroup , *targetgroup ,
//                                                    NULL , Asynchronous ) ;
  MPI_AccessDEC * MPIAccessDEC = new MPI_AccessDEC( *sourcegroup , *targetgroup ,
                                                    Asynchronous ) ;
  
  MPI_Access * mpi_access = MPIAccessDEC->MPIAccess() ;

#define maxreq 100
#define datamsglength 10

//  int sts ;
  int sendcount = datamsglength ;
  int recvcount = datamsglength ;
  int * recvbuf = new int[datamsglength*size] ;

  int ireq ;
  for ( ireq = 0 ; ireq < maxreq ; ireq++ ) {
    int * sendbuf = new int[datamsglength*size] ;
//    cout << "test_AllToAllDEC" << myrank << " ireq " << ireq << " RecvRequestIdsSize "
//         << mpi_access->RecvRequestIdsSize() << endl ;
//    int * sendbuf = (int *) malloc( sizeof(int)*datamsglength*size) ;
    int j ;
    for ( j = 0 ; j < datamsglength*size ; j++ ) {
       sendbuf[j] = myrank*1000000 + ireq*1000 + j ;
       recvbuf[j] = -1 ;
    }

    MPIAccessDEC->AllToAll( sendbuf, sendcount , MPI_INT ,
	                    recvbuf, recvcount , MPI_INT ) ;

//    cout << "test_AllToAllDEC" << myrank << " recvbuf before CheckSent" ;
//    for ( i = 0 ; i < datamsglength*size ; i++ ) {
//       cout << " " << recvbuf[i] ;
//    }
//    cout << endl ;

//    cout << "test_AllToAllDEC" << myrank << " sendbuf " << sendbuf << endl ;
//    MPIAccessDEC->CheckSent() ;

    int nRecvReq = mpi_access->RecvRequestIdsSize() ;
//    cout << "test_AllToAllDEC" << myrank << " WaitAllRecv " << nRecvReq << " Requests" << endl ;
    int *ArrayOfRecvRequests = new int[nRecvReq] ;
    int nReq = mpi_access->RecvRequestIds( nRecvReq, ArrayOfRecvRequests ) ;
    mpi_access->WaitAll( nReq , ArrayOfRecvRequests ) ;
    mpi_access->DeleteRequests( nReq , ArrayOfRecvRequests ) ;
    delete [] ArrayOfRecvRequests ;
//    cout << "test_AllToAllDEC" << myrank << " RecvRequestIdsSize " << mpi_access->RecvRequestIdsSize()
//         << " after WaitAll" << endl ;

//    cout << "test_AllToAllDEC" << myrank << " recvbuf" ;
//    for ( i = 0 ; i < datamsglength*size ; i++ ) {
//       cout << " " << recvbuf[i] ;
//    }
//    cout << endl ;
  }

//  cout << "test_AllToAllDEC" << myrank << " final CheckSent" << endl ;
//  MPIAccessDEC->CheckSent() ;

  int nSendReq = mpi_access->SendRequestIdsSize() ;
  cout << "test_AllToAllDEC" << myrank << " final SendRequestIds " << nSendReq << " SendRequests"
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
    strstream << "test_AllToAllDEC" << myrank << " final RecvRequestIds " << nRecvReq
              << " RecvRequests # 0 Error" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    cout << "test_AllToAllDEC" << myrank << " final RecvRequestIds " << nRecvReq
         << " RecvRequests = 0 OK" << endl ;
  }

  mpi_access->Barrier() ;

  delete sourcegroup ;
  delete targetgroup ;
  delete MPIAccessDEC ;
  delete [] recvbuf ;

//  MPI_Finalize();

  cout << "test_AllToAllDEC" << myrank << " OK" << endl ;

  return ;
}




