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

#include <math.h>
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
using namespace ParaMEDMEM;

void MPIAccessDECTest::test_AllToAllvTimeDoubleDECSynchronousPointToPoint() {
  test_AllToAllvTimeDoubleDEC( false ) ;
}
void MPIAccessDECTest::test_AllToAllvTimeDoubleDECAsynchronousPointToPoint() {
  test_AllToAllvTimeDoubleDEC( true ) ;
}

static void chksts( int sts , int myrank , ParaMEDMEM::MPIAccess * mpi_access ) {
  char msgerr[MPI_MAX_ERROR_STRING] ;
  int lenerr ;
  if ( sts != MPI_SUCCESS ) {
    mpi_access->errorString(sts, msgerr, &lenerr) ;
    cout << "test" << myrank << " lenerr " << lenerr << " "
         << msgerr << endl ;
    ostringstream strstream ;
    strstream << "==========================================================="
              << "test" << myrank << " KO"
              << "==========================================================="
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  return ;
}

void MPIAccessDECTest::test_AllToAllvTimeDoubleDEC( bool Asynchronous ) {

  cout << "test_AllToAllvTimeDoubleDEC" << endl ;

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
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

//  int Asynchronous = atoi(argv[1]) ;

  cout << "test_AllToAllvTimeDoubleDEC" << myrank << " Asynchronous " << Asynchronous << endl ;

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

//  TimeInterpolator * aLinearInterpDEC = new LinearTimeInterpolator( 0 ) ;
  MPIAccessDEC * MyMPIAccessDEC = new MPIAccessDEC( *sourcegroup , *targetgroup ,
                                                    Asynchronous ) ;
//                                                    Asynchronous , LinearInterp , 0.5 ) ;
  MyMPIAccessDEC->setTimeInterpolator( LinearTimeInterp ) ;
  MPIAccess * mpi_access = MyMPIAccessDEC->getMPIAccess() ;

  cout << "test_AllToAllvTimeDoubleDEC" << myrank << " Barrier :" << endl ;
  mpi_access->barrier() ;

#define maxproc 11
#define maxreq 100
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
  }

  double timeLoc[maxproc] ;
  double deltatime[maxproc] = {1.,2.1,3.2,4.3,5.4,6.5,7.6,8.7,9.8,10.9,11.} ;
  double maxtime[maxproc] ;
  double nextdeltatime[maxproc] ;
  for ( i = 0 ; i < size ; i++ ) {
     timeLoc[i] = 0 ;
     maxtime[i] = maxreq ;
     nextdeltatime[i] = deltatime[i] ;
  }
  time_t begintime = time(NULL) ;
  for ( timeLoc[myrank] = 0 ; timeLoc[myrank] <= maxtime[myrank] && nextdeltatime[myrank] != 0 ;
        timeLoc[myrank]+=nextdeltatime[myrank] ) {
//local and target times
     int target ;
     for ( target = 0 ; target < size ; target++ ) {
        nextdeltatime[target] = deltatime[target] ;
        if ( timeLoc[target] != 0 ) {
          if ( timeLoc[target]+nextdeltatime[target] > maxtime[target] ) {
            nextdeltatime[target] = 0 ;
          }
        }
        if ( target != myrank ) {
          while ( timeLoc[myrank] >= timeLoc[target] ) {
               timeLoc[target] += deltatime[target] ;
          }
        }
     }
     MyMPIAccessDEC->setTime( timeLoc[myrank] , nextdeltatime[myrank] ) ;
     cout << "test" << myrank << "=====TIME " << timeLoc[myrank] << "=====DELTATIME "
          << nextdeltatime[myrank] << "=====MAXTIME " << maxtime[myrank] << " ======"
          << endl ; 
     double * sendbuf = new double[datamsglength*size] ;
//     double * sendbuf = (double *) malloc(sizeof(double)*datamsglength*size) ;
     double * recvbuf = new double[datamsglength*size] ;
     int j ;
     //cout << "test_AllToAllvTimeDoubleDEC" << myrank << " sendbuf" ;
     for ( target = 0 ; target < size ; target++ ) {
        for ( j = 0 ; j < datamsglength ; j++ ) {
           //sendbuf[j] = myrank*10000 + (j/datamsglength)*100 + j ;
           sendbuf[target*datamsglength+j] = myrank*1000000 + target*10000 +
                                             (timeLoc[myrank]/deltatime[myrank])*100 + j ;
           //cout << " " << (int ) sendbuf[target*datamsglength+j] ;
           recvbuf[target*datamsglength+j] = -1 ;
        }
        //cout << endl ;
     }

     int sts = MyMPIAccessDEC->allToAllvTime( sendbuf, sendcounts , sdispls , MPI_DOUBLE ,
                                            recvbuf, recvcounts , rdispls , MPI_DOUBLE ) ;
     chksts( sts , myrank , mpi_access ) ;

//     cout << "test_AllToAllvTimeDoubleDEC" << myrank << " recvbuf before CheckSent" ;
//     for ( i = 0 ; i < datamsglength*size ; i++ ) {
//        cout << " " << recvbuf[i] ;
//     }
//     cout << endl ;

     int nRecvReq = mpi_access->recvRequestIdsSize() ;
     if ( nRecvReq != 0 ) {
       ostringstream strstream ;
       strstream << "=============================================================" << endl
                 << "test_AllToAllvTimeDoubleDEC" << myrank << " WaitAllRecv "
                 << nRecvReq << " Requests # 0 ERROR"
                 << endl << "============================================================"
                 << endl ;
       int *ArrayOfRecvRequests = new int[nRecvReq] ;
       int nReq = mpi_access->recvRequestIds( nRecvReq, ArrayOfRecvRequests ) ;
       mpi_access->waitAll( nReq , ArrayOfRecvRequests ) ;
       delete [] ArrayOfRecvRequests ;
       cout << strstream.str() << endl ;
       CPPUNIT_FAIL( strstream.str() ) ;
     }

//     cout << "test_AllToAllvTimeDoubleDEC" << myrank << " check of recvbuf" << endl ;
     bool badrecvbuf = false ;
     for ( target = 0 ; target < size ; target++ ) {
        int j ;
        for ( j = 0 ; j < datamsglength ; j++ ) {
           int index = target*datamsglength+j ;
           if ( j < recvcounts[target] ) {
             if ( fabs(recvbuf[index] - (target*1000000 + myrank*10000 +
                  (timeLoc[target]/deltatime[target])*100 + j)) > 101) {
               badrecvbuf = true ;
               cout << "test_AllToAllvTimeDoubleDEC" << myrank << " target " << target << " timeLoc[target] "
                    << timeLoc[target] << " recvbuf[" << index << "] " << (int ) recvbuf[index]
                    << " # " << (int ) (target*1000000 +
                       myrank*10000 + (timeLoc[target]/deltatime[target])*100 + j)
                    << endl ;
             }
             else if ( badrecvbuf ) {
               cout << "test_AllToAllvTimeDoubleDEC" << myrank << " recvbuf[" << index << "] "
                    << recvbuf[index] << " ~= " << (int ) (target*1000000 +
                       myrank*10000 + (timeLoc[target]/deltatime[target])*100 + j) << endl ;
             }
           }
           else if ( recvbuf[index] != -1 ) {
             badrecvbuf = true ;
             cout << "test_AllToAllvTimeDoubleDEC" << myrank << " recvbuf[" << index << "] "
                  << recvbuf[index] << " # -1" << endl ;
           }
        }
     }
     if ( badrecvbuf ) {
       ostringstream strstream ;
       strstream << "==================================================================" << endl
                 << "test_AllToAllvTimeDoubleDEC" << myrank << " badrecvbuf"
                 << endl << "=================================================================="
                 << endl ;
       cout << strstream.str() << endl ;
       CPPUNIT_FAIL( strstream.str() ) ;
     }
     delete [] recvbuf ;
  }

  cout << "test_AllToAllvTimeDoubleDEC" << myrank << " Barrier :" << endl ;
  mpi_access->barrier() ;

  cout << "test_AllToAllvTimeDoubleDEC" << myrank << " CheckFinalSent" << endl ;
  sts = MyMPIAccessDEC->checkFinalSent() ;
  if ( sts != MPI_SUCCESS ) {
    ostringstream strstream ;
    strstream << "=================================================================" << endl
              << "test_AllToAllvTimeDoubleDEC" << myrank << " CheckFinalSent ERROR"
              << endl << "================================================================="
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  cout << "test_AllToAllvTimeDoubleDEC" << myrank << " CheckFinalRecv" << endl ;
  sts = MyMPIAccessDEC->checkFinalRecv() ;
  if ( sts != MPI_SUCCESS ) {
    ostringstream strstream ;
    strstream << "=================================================================" << endl
              << "test_AllToAllvTimeDoubleDEC" << myrank << " CheckFinalRecv ERROR"
              << endl << "================================================================"
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  int nRecvReq = mpi_access->recvRequestIdsSize() ;
  if ( nRecvReq ) {
    ostringstream strstream ;
    strstream << "===============================================================" << endl
              << "test_AllToAllvTimeDoubleDEC" << myrank << " RecvRequestIds " << nRecvReq
              << " RecvRequests # 0 Error"
              << endl << "==============================================================="
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    cout << "test_AllToAllvTimeDoubleDEC" << myrank << " RecvRequestIds " << nRecvReq
         << " RecvRequests = 0 OK" << endl ;
  }

  time_t endtime = time(NULL) ;
  cout << "test_AllToAllvTimeDoubleDEC" << myrank << " begintime " << begintime << " endtime " << endtime
       << " elapse " << endtime-begintime << " " << maxtime[myrank]/deltatime[myrank]
       << " calls to AllToAll" << endl ;

  cout << "test" << myrank << " Barrier :" << endl ;
  mpi_access->barrier() ;

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

  cout << "test_AllToAllvTimeDoubleDEC" << myrank << " OK begintime " << begintime << " endtime " << endtime
       << " elapse " << endtime-begintime << " " << maxtime[myrank]/deltatime[myrank]
       << " calls to AllToAll" << endl ;

  return ;
}




