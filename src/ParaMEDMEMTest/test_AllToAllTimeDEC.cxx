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
#include "MPIAccessDEC.hxx"
#include "LinearTimeInterpolator.hxx"

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES

using namespace std;
using namespace MEDCoupling;

void MPIAccessDECTest::test_AllToAllTimeDECSynchronousPointToPoint() {
  test_AllToAllTimeDEC( false ) ;
}
void MPIAccessDECTest::test_AllToAllTimeDECAsynchronousPointToPoint() {
  test_AllToAllTimeDEC( true ) ;
}

static void chksts( int sts , int myrank , MEDCoupling::MPIAccess * mpi_access ) {
  char msgerr[MPI_MAX_ERROR_STRING] ;
  int lenerr ;
  if ( sts != MPI_SUCCESS ) {
    mpi_access->errorString(sts, msgerr, &lenerr) ;
    debugStream << "test_AllToAllTimeDEC" << myrank << " lenerr " << lenerr << " "
         << msgerr << endl ;
    ostringstream strstream ;
    strstream << "==========================================================="
              << "test_AllToAllTimeDEC" << myrank << " KO"
              << "==========================================================="
              << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  return ;
}

void MPIAccessDECTest::test_AllToAllTimeDEC( bool Asynchronous ) {

  debugStream << "test_AllToAllTimeDEC" << endl ;

  //  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 || size > 11 ) {
    ostringstream strstream ;
    strstream << "usage :" << endl
              << "mpirun -np <nbprocs> test_AllToAllTimeDEC" << endl
              << " (nbprocs >=2)" << endl
              << "test must be runned with more than 1 proc and less than 12 procs"
              << endl ;
    cerr << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  //  int Asynchronous = atoi(argv[1]);

  debugStream << "test_AllToAllTimeDEC" << myrank << " Asynchronous " << Asynchronous << endl ;

  MEDCoupling::CommInterface interface ;
  std::set<int> sourceprocs;
  std::set<int> targetprocs;
  int i ;
  for ( i = 0 ; i < size/2 ; i++ ) {
    sourceprocs.insert(i);
  }
  for ( i = size/2 ; i < size ; i++ ) {
    targetprocs.insert(i);
  }

  MEDCoupling::MPIProcessorGroup* sourcegroup = new MEDCoupling::MPIProcessorGroup(interface,sourceprocs) ;
  MEDCoupling::MPIProcessorGroup* targetgroup = new MEDCoupling::MPIProcessorGroup(interface,targetprocs) ;

  //  LinearTimeInterpolator * aLinearInterpDEC = new LinearTimeInterpolator( 0.5 ) ;
  MPIAccessDEC * MyMPIAccessDEC = new MPIAccessDEC( *sourcegroup , *targetgroup ,
                                                    Asynchronous ) ;
  //                                                    Asynchronous , LinearInterp , 0.5 ) ;
  MyMPIAccessDEC->setTimeInterpolator( LinearTimeInterp ) ;
  MPIAccess * mpi_access = MyMPIAccessDEC->getMPIAccess() ;

  debugStream << "test_AllToAllTimeDEC" << myrank << " Barrier :" << endl ;
  mpi_access->barrier() ;
  debugStream << "test_AllToAllTimeDEC" << myrank << " Barrier done" << endl ;
  
#define maxproc 11
#define maxreq 10000
#define datamsglength 10

  int sts ;
  int sendcount = datamsglength ;
  int recvcount = datamsglength ;

  double time = 0 ;
  //  double deltatime[maxproc] = {1.,2.1,3.2,4.3,5.4,6.5,7.6,8.7,9.8,10.9,11.} ;
  double deltatime[maxproc] = {1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.} ;
  double maxtime = maxreq ;
  double nextdeltatime = deltatime[myrank] ;
  //  MyMPIAccessDEC->InitTime( time , deltatime[myrank] , maxtime ) ;
  //  for ( time = 0 ; time <= maxtime ; time+=deltatime[myrank] ) {
  for ( time = 0 ; time <= maxtime && nextdeltatime != 0 ; time+=nextdeltatime ) {
    if ( time != 0 ) {
      nextdeltatime = deltatime[myrank] ;
      if ( time+nextdeltatime > maxtime ) {
        nextdeltatime = 0 ;
      }
      //       MyMPIAccessDEC->NextTime( nextdeltatime ) ;
    }
    MyMPIAccessDEC->setTime( time , nextdeltatime ) ;
    debugStream << "test_AllToAllTimeDEC" << myrank << "=====TIME " << time << "=====DELTATIME "
         << nextdeltatime << "=====MAXTIME " << maxtime << " ======" << endl ; 
    int * sendbuf = new int[datamsglength*size] ;
    //     int * sendbuf = (int *) malloc(sizeof(int)*datamsglength*size) ;
    int * recvbuf = new int[datamsglength*size] ;
    int j ;
    for ( j = 0 ; j < datamsglength*size ; j++ ) {
      sendbuf[j] = myrank*1000000 + (j/datamsglength)*1000 + j ;
      recvbuf[j] = -1 ;
    }

    int sts = MyMPIAccessDEC->allToAllTime( sendbuf, sendcount , MPI_INT ,
                                            recvbuf, recvcount , MPI_INT ) ;
    chksts( sts , myrank , mpi_access ) ;

    //     debugStream << "test_AllToAllTimeDEC" << myrank << " recvbuf before CheckSent" ;
    //     for ( i = 0 ; i < datamsglength*size ; i++ ) {
    //        debugStream << " " << recvbuf[i] ;
    //     }
    //     debugStream << endl ;

    //     debugStream << "test_AllToAllTimeDEC" << myrank << " sendbuf " << sendbuf << endl ;
    //     MyMPIAccessDEC->CheckSent() ;

    int nRecvReq = mpi_access->recvRequestIdsSize() ;
    if ( nRecvReq != 0 ) {
      ostringstream strstream ;
      strstream << "=============================================================" << endl
                << "test_AllToAllTimeDEC" << myrank << " WaitAllRecv " << nRecvReq << " Requests # 0 ERROR"
                << endl << "============================================================="
                << endl ;
      int *ArrayOfRecvRequests = new int[nRecvReq] ;
      int nReq = mpi_access->recvRequestIds( nRecvReq, ArrayOfRecvRequests ) ;
      mpi_access->waitAll( nReq , ArrayOfRecvRequests ) ;
      delete [] ArrayOfRecvRequests ;
      debugStream << strstream.str() << endl ;
      CPPUNIT_FAIL( strstream.str() ) ;
    }

    //     debugStream << "test_AllToAllTimeDEC" << myrank << " recvbuf" << endl ;
    bool badrecvbuf = false ;
    for ( i = 0 ; i < datamsglength*size ; i++ ) {
      if ( recvbuf[i] != (i/datamsglength)*1000000 + myrank*1000 +
           myrank*datamsglength+(i%datamsglength) ) {
        badrecvbuf = true ;
        debugStream << "test_AllToAllTimeDEC" << myrank << " recvbuf[" << i << "] "
             << recvbuf[i] << " # " << (i/datamsglength)*1000000 + myrank*1000 +
          myrank*datamsglength+(i%datamsglength) << endl ;
      }
      else if ( badrecvbuf ) {
        debugStream << "test_AllToAllTimeDEC" << myrank << " recvbuf[" << i << "] "
             << recvbuf[i] << " == " << (i/datamsglength)*1000000 + myrank*1000 +
          myrank*datamsglength+(i%datamsglength) << endl ;
      }
    }
    if ( badrecvbuf ) {
      ostringstream strstream ;
      strstream << "==============================================================" << endl
                << "test_AllToAllTimeDEC" << myrank << " badrecvbuf"
                << endl << "============================================================="
                << endl ;
      debugStream << strstream.str() << endl ;
      CPPUNIT_FAIL( strstream.str() ) ;
    }
    delete [] recvbuf ;
  }

  debugStream << "test_AllToAllTimeDEC" << myrank << " final CheckSent" << endl ;
  sts = MyMPIAccessDEC->checkSent() ;
  if ( sts != MPI_SUCCESS ) {
    ostringstream strstream ;
    strstream << "================================================================" << endl
              << "test_AllToAllTimeDEC" << myrank << " final CheckSent ERROR"
              << endl << "================================================================"
              << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  int nSendReq = mpi_access->sendRequestIdsSize() ;
  debugStream << "test_AllToAllTimeDEC" << myrank << " final SendRequestIds " << nSendReq << " SendRequests"
       << endl ;
  if ( nSendReq ) {
    int *ArrayOfSendRequests = new int[nSendReq] ;
    int nReq = mpi_access->sendRequestIds( nSendReq, ArrayOfSendRequests ) ;
    mpi_access->waitAll( nReq , ArrayOfSendRequests ) ;
    delete [] ArrayOfSendRequests ;
  }

  int nRecvReq = mpi_access->recvRequestIdsSize() ;
  if ( nRecvReq ) {
    ostringstream strstream ;
    strstream << "===============================================================" << endl
              << "test_AllToAllTimeDEC" << myrank << " RecvRequestIds " << nRecvReq
              << " RecvRequests # 0 Error"
              << endl << "==============================================================="
              << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    debugStream << "test_AllToAllTimeDEC" << myrank << " RecvRequestIds " << nRecvReq
         << " RecvRequests = 0 OK" << endl ;
  }

  debugStream << "test_AllToAllTimeDEC" << myrank << " Barrier :" << endl ;
  mpi_access->barrier() ;
  debugStream << "test_AllToAllTimeDEC" << myrank << " Barrier done" << endl ;

  delete sourcegroup ;
  delete targetgroup ;
  //  delete aLinearInterpDEC ;
  delete MyMPIAccessDEC ;

  //  MPI_Finalize();

  debugStream << "test_AllToAllTimeDEC" << myrank << " OK" << endl ;

  return ;
}




