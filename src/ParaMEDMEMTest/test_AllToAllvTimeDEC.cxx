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
#include <ctime>

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

void MPIAccessDECTest::test_AllToAllvTimeDECSynchronousNative() {
  test_AllToAllvTimeDEC( false , true ) ;
}
void MPIAccessDECTest::test_AllToAllvTimeDECSynchronousPointToPoint() {
  test_AllToAllvTimeDEC( false , false ) ;
}
void MPIAccessDECTest::test_AllToAllvTimeDECAsynchronousPointToPoint() {
  test_AllToAllvTimeDEC( true , false ) ;
}

static void chksts( int sts , int myrank , MEDCoupling::MPIAccess * mpi_access ) {
  char msgerr[MPI_MAX_ERROR_STRING] ;
  int lenerr ;
  if ( sts != MPI_SUCCESS ) {
    mpi_access->errorString(sts, msgerr, &lenerr) ;
    debugStream << "test_AllToAllvTimeDEC" << myrank << " lenerr " << lenerr << " "
         << msgerr << endl ;
    ostringstream strstream ;
    strstream << "==========================================================="
              << "test_AllToAllvTimeDEC" << myrank << " KO"
              << "==========================================================="
              << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  return ;
}

void MPIAccessDECTest::test_AllToAllvTimeDEC( bool Asynchronous , bool UseMPINative ) {

  debugStream << "test_AllToAllvTimeDEC" << endl ;

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

  //  int Asynchronous = atoi(argv[1]) ;
  int UseMPI_Alltoallv = UseMPINative ;
  //  if ( argc == 3 ) {
  //    UseMPI_Alltoallv = atoi(argv[2]) ;
  //  }

  debugStream << "test_AllToAllvTimeDEC" << myrank << " Asynchronous " << Asynchronous
       << " UseMPI_Alltoallv " << UseMPI_Alltoallv << endl ;

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

  //  TimeInterpolator * aLinearInterpDEC = new LinearTimeInterpolator( 0.5 ) ;
  MPIAccessDEC * MyMPIAccessDEC = new MPIAccessDEC( *sourcegroup , *targetgroup ,
                                                    Asynchronous ) ;
  //                                                    Asynchronous , LinearInterp , 0.5 ) ;
  MyMPIAccessDEC->setTimeInterpolator( LinearTimeInterp , 0.5 ) ;
  MPIAccess * mpi_access = MyMPIAccessDEC->getMPIAccess() ;

  debugStream << "test_AllToAllvTimeDEC" << myrank << " Barrier :" << endl ;
  mpi_access->barrier() ;
  debugStream << "test_AllToAllvTimeDEC" << myrank << " Barrier done" << endl ;

#define maxproc 11
#define maxreq 10000
#define datamsglength 10

  int sts ;
  int *sendcounts = new int[size] ;
  int *sdispls = new int[size] ;
  int *recvcounts = new int[size] ;
  int *rdispls = new int[size] ;
  int *sendtimecounts = new int[size] ;
  int *stimedispls = new int[size] ;
  int *recvtimecounts = new int[size] ;
  int *rtimedispls = new int[size] ;
  for ( i = 0 ; i < size ; i++ ) {
    sendcounts[i] = datamsglength-i ;
    sdispls[i] = i*datamsglength ;
    recvcounts[i] = datamsglength-myrank ;
    rdispls[i] = i*datamsglength ;
    sendtimecounts[i] = 1 ;
    stimedispls[i] = 0 ;
    recvtimecounts[i] = 1 ;
    rtimedispls[i] = i ;
    //rtimedispls[i] = i*mpi_access->TimeExtent() ;
  }

  double timeLoc = 0 ;
  double deltatime[maxproc] = {1.,2.1,3.2,4.3,5.4,6.5,7.6,8.7,9.8,10.9,11.} ;
  double maxtime ;
  double nextdeltatime = deltatime[myrank] ;
  if ( UseMPI_Alltoallv ) {
    maxtime = maxreq*nextdeltatime - 0.1 ;
  }
  else {
    maxtime = maxreq ;
    //    MyMPIAccessDEC->InitTime( time , nextdeltatime , maxtime ) ;
  }
  time_t begintime = time(NULL) ;
  //  for ( time = 0 ; time <= maxtime ; time+=deltatime[myrank] ) {
  for ( timeLoc = 0 ; timeLoc <= maxtime && nextdeltatime != 0 ; timeLoc+=nextdeltatime ) {
    nextdeltatime = deltatime[myrank] ;
    if ( timeLoc != 0 ) {
      nextdeltatime = deltatime[myrank] ;
      if ( timeLoc+nextdeltatime > maxtime ) {
        nextdeltatime = 0 ;
      }
      //       MyMPIAccessDEC->NextTime( nextdeltatime ) ;
    }
    MyMPIAccessDEC->setTime( timeLoc , nextdeltatime ) ;
    debugStream << "test_AllToAllvTimeDEC" << myrank << "=====TIME " << time << "=====DELTATIME "
         << nextdeltatime << "=====MAXTIME " << maxtime << " ======" << endl ; 
    int * sendbuf = new int[datamsglength*size] ;
    //     int * sendbuf = (int *) malloc(sizeof(int)*datamsglength*size) ;
    int * recvbuf = new int[datamsglength*size] ;
    int j ;
    for ( j = 0 ; j < datamsglength*size ; j++ ) {
      sendbuf[j] = myrank*1000000 + (j/datamsglength)*1000 + j ;
      recvbuf[j] = -1 ;
    }

    if ( UseMPI_Alltoallv ) {
      const MPI_Comm* comm = MyMPIAccessDEC->getComm();
      TimeMessage * aSendTimeMessage = new TimeMessage ;
      aSendTimeMessage->time = timeLoc ;
      //       aSendTimeMessage->deltatime = deltatime[myrank] ;
      aSendTimeMessage->deltatime = nextdeltatime ;
      //       aSendTimeMessage->maxtime = maxtime ;
      aSendTimeMessage->tag = (int ) (timeLoc/deltatime[myrank]) ;
      TimeMessage * aRecvTimeMessage = new TimeMessage[size] ;
      interface.allToAllV(aSendTimeMessage, sendtimecounts , stimedispls ,
                          mpi_access->timeType() ,
                          aRecvTimeMessage, recvtimecounts , rtimedispls ,
                          mpi_access->timeType() , *comm ) ;
      //       for ( j = 0 ; j < size ; j++ ) {
      //          debugStream << "test_AllToAllvTimeDEC" << myrank << " TimeMessage received " << j << " "
      //               << aRecvTimeMessage[j] << endl ;
      //       }
      delete aSendTimeMessage ;
      delete [] aRecvTimeMessage ;
      interface.allToAllV(sendbuf, sendcounts , sdispls , MPI_INT ,
                          recvbuf, recvcounts , rdispls , MPI_INT , *comm ) ;
      //       free(sendbuf) ;
      delete [] sendbuf ;
    }
    else {
      int sts = MyMPIAccessDEC->allToAllvTime( sendbuf, sendcounts , sdispls , MPI_INT ,
                                               recvbuf, recvcounts , rdispls , MPI_INT ) ;
      chksts( sts , myrank , mpi_access ) ;
    }

    //     debugStream << "test_AllToAllvTimeDEC" << myrank << " recvbuf before CheckSent" ;
    //     for ( i = 0 ; i < datamsglength*size ; i++ ) {
    //        debugStream << " " << recvbuf[i] ;
    //     }
    //     debugStream << endl ;

    //     debugStream << "test_AllToAllvTimeDEC" << myrank << " sendbuf " << sendbuf << endl ;
    //     MyMPIAccessDEC->CheckSent() ;

    int nRecvReq = mpi_access->recvRequestIdsSize() ;
    if ( nRecvReq != 0 ) {
      ostringstream strstream ;
      strstream << "=============================================================" << endl
                << "test_AllToAllvTimeDEC" << myrank << " WaitAllRecv " << nRecvReq << " Requests # 0 ERROR"
                << endl << "============================================================="
                << endl ;
      int *ArrayOfRecvRequests = new int[nRecvReq] ;
      int nReq = mpi_access->recvRequestIds( nRecvReq, ArrayOfRecvRequests ) ;
      mpi_access->waitAll( nReq , ArrayOfRecvRequests ) ;
      delete [] ArrayOfRecvRequests ;
      debugStream << strstream.str() << endl ;
      CPPUNIT_FAIL( strstream.str() ) ;
    }

    //     debugStream << "test_AllToAllvTimeDEC" << myrank << " check of recvbuf" << endl ;
    bool badrecvbuf = false ;
    for ( i = 0 ; i < size ; i++ ) {
      int j ;
      for ( j = 0 ; j < datamsglength ; j++ ) {
        int index = i*datamsglength+j ;
        if ( j < recvcounts[i] ) {
          if ( recvbuf[index] != (index/datamsglength)*1000000 + myrank*1000 +
               myrank*datamsglength+(index%datamsglength) ) {
            badrecvbuf = true ;
            debugStream << "test_AllToAllvTimeDEC" << myrank << " recvbuf[" << index << "] "
                 << recvbuf[index] << " # " << (index/datamsglength)*1000000 +
              myrank*1000 +
              myrank*datamsglength+(index%datamsglength) << endl ;
          }
          else if ( badrecvbuf ) {
            debugStream << "test_AllToAllvTimeDEC" << myrank << " recvbuf[" << index << "] "
                 << recvbuf[index] << " == " << (index/datamsglength)*1000000 +
              myrank*1000 +
              myrank*datamsglength+(index%datamsglength) << endl ;
          }
        }
        else if ( recvbuf[index] != -1 ) {
          badrecvbuf = true ;
          debugStream << "test_AllToAllvTimeDEC" << myrank << " recvbuf[" << index << "] "
               << recvbuf[index] << " # -1" << endl ;
        }
      }
    }
    if ( badrecvbuf ) {
      ostringstream strstream ;
      strstream << "==============================================================" << endl
                << "test_AllToAllvTimeDEC" << myrank << " badrecvbuf"
                << endl << "============================================================="
                << endl ;
      debugStream << strstream.str() << endl ;
      CPPUNIT_FAIL( strstream.str() ) ;
    }
    delete [] recvbuf ;
  }

  debugStream << "test_AllToAllvTimeDEC" << myrank << " Barrier :" << endl ;
  mpi_access->barrier() ;
  debugStream << "test_AllToAllvTimeDEC" << myrank << " Barrier done" << endl ;

  debugStream << "test_AllToAllvTimeDEC" << myrank << " CheckFinalSent" << endl ;
  sts = MyMPIAccessDEC->checkFinalSent() ;
  if ( sts != MPI_SUCCESS ) {
    ostringstream strstream ;
    strstream << "================================================================" << endl
              << "test_AllToAllvTimeDEC" << myrank << " final CheckSent ERROR"
              << endl << "================================================================"
              << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  debugStream << "test_AllToAllvTimeDEC" << myrank << " CheckFinalRecv" << endl ;
  sts = MyMPIAccessDEC->checkFinalRecv() ;
  if ( sts != MPI_SUCCESS ) {
    ostringstream strstream ;
    strstream << "================================================================" << endl
              << "test_AllToAllvTimeDEC" << myrank << " CheckFinalRecv ERROR"
              << endl << "================================================================"
              << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  int nRecvReq = mpi_access->recvRequestIdsSize() ;
  if ( nRecvReq ) {
    ostringstream strstream ;
    strstream << "===============================================================" << endl
              << "test_AllToAllvTimeDEC" << myrank << " RecvRequestIds " << nRecvReq
              << " RecvRequests # 0 Error"
              << endl << "==============================================================="
              << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    debugStream << "test_AllToAllvTimeDEC" << myrank << " RecvRequestIds " << nRecvReq
         << " RecvRequests = 0 OK" << endl ;
  }

  time_t endtime = time(NULL) ;
  debugStream << "test_AllToAllvTimeDEC" << myrank << " begintime " << begintime << " endtime " << endtime
       << " elapse " << endtime-begintime << " " << maxtime/deltatime[myrank]
       << " calls to AllToAll" << endl ;

  debugStream << "test_AllToAllvTimeDEC" << myrank << " Barrier :" << endl ;
  mpi_access->barrier() ;
  debugStream << "test_AllToAllvTimeDEC" << myrank << " Barrier done" << endl ;

  delete sourcegroup ;
  delete targetgroup ;
  delete MyMPIAccessDEC ;
  //  delete aLinearInterpDEC ;

  delete [] sendcounts ;
  delete [] sdispls ;
  delete [] recvcounts ;
  delete [] rdispls ;
  delete [] sendtimecounts ;
  delete [] stimedispls ;
  delete [] recvtimecounts ;
  delete [] rtimedispls ;

  //  MPI_Finalize();

  endtime = time(NULL) ;

  debugStream << "test_AllToAllvTimeDEC" << myrank << " OK begintime " << begintime << " endtime " << endtime
       << " elapse " << endtime-begintime << " " << maxtime/deltatime[myrank]
       << " calls to AllToAll" << endl ;

  return ;
}




