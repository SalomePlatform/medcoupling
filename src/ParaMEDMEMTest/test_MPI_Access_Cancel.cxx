// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <mpi.h>

#ifndef WIN32
#include <unistd.h>
#endif

#include "MPIAccessTest.hxx"
#include <cppunit/TestAssert.h>

//#include "CommInterface.hxx"
//#include "ProcessorGroup.hxx"
//#include "MPIProcessorGroup.hxx"
#include "MPIAccess.hxx"

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES

using namespace std;
using namespace ParaMEDMEM;

void MPIAccessTest::test_MPI_Access_Cancel() {

  cout << "test_MPI_Access_Cancel" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 ) {
    ostringstream strstream ;
    strstream << "test_MPI_Access_Cancel must be runned with 2 procs" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  cout << "test_MPI_Access_Cancel" << myrank << endl ;

  ParaMEDMEM::CommInterface interface ;

  ParaMEDMEM::MPIProcessorGroup* group = new ParaMEDMEM::MPIProcessorGroup(interface) ;

  ParaMEDMEM::MPIAccess mpi_access( group ) ;

  if ( myrank >= 2 ) {
    mpi_access.barrier() ;
    delete group ;
    return ;
  }

  int target = 1 - myrank ;
  int intsendbuf[5] ;
  double doublesendbuf[10] ;
  int RequestId[10] ;
  int sts ;
  int i , j ;
  for ( j = 0 ; j < 3 ; j++ ) {
     for ( i = 0 ; i < 10 ; i++ ) {
        cout << "test" << myrank << " ============================ i " << i
             << "============================" << endl ;
        if ( myrank == 0 ) {
          if ( i < 5 ) {
            intsendbuf[i] = i ;
            sts = mpi_access.ISend(&intsendbuf[i],1,MPI_INT,target, RequestId[i]) ;
            cout << "test" << myrank << " Send MPI_INT RequestId " << RequestId[i]
                 << endl ;
          }
          else {
            doublesendbuf[i] = i ;
            sts = mpi_access.ISend(&doublesendbuf[i],1,MPI_DOUBLE,target,
                                   RequestId[i]) ;
            cout << "test" << myrank << " Send MPI_DOUBLE RequestId " << RequestId[i]
                 << endl ;
          }
        }
        else {
          int flag = false ;
          while ( !flag ) {
               int source, tag, outcount ;
               MPI_Datatype datatype ;
               sts = mpi_access.IProbe(target, source, tag, datatype, outcount,
                                       flag ) ;
               if ( flag ) {
                 cout << "test" << myrank << " " << i << " IProbe target " << target
                      << " source " << source << " tag " << tag
                      << " outcount " << outcount << " flag " << flag << endl ;
               }
               else {
                 cout << "test" << myrank << " flag " << flag << endl ;
                 sleep( 1 ) ;
               }
               if ( flag ) {
                 int recvbuf ;
                 sts = mpi_access.IRecv(&recvbuf,outcount,MPI_INT,source,
                                        RequestId[i] ) ;
                 if ( datatype == MPI_INT ) {
                   int source, tag, error, outcount ;
                   mpi_access.wait( RequestId[i] ) ;
                   mpi_access.status( RequestId[i], source, tag, error, outcount,
                                      true ) ;
                   if ( (outcount != 1) | (recvbuf != i) ) {
                     ostringstream strstream ;
                     strstream << "======================================================"
                               << endl << "test" << myrank << " outcount " << outcount
                               << " recvbuf " << recvbuf << " KO" << endl
                               << "======================================================"
                               << endl ;
                     cout << strstream.str() << endl ;
                     CPPUNIT_FAIL( strstream.str() ) ;
                   }
                   cout << "========================================================"
                        << endl << "test" << myrank << " outcount " << outcount
                        << " recvbuf " << recvbuf << " OK" << endl
                        << "========================================================"
                        << endl ;
                 }
               }
          }
        }
        char msgerr[MPI_MAX_ERROR_STRING] ;
        int lenerr ;
        mpi_access.errorString(sts, msgerr, &lenerr) ;
        cout << "test" << myrank << " lenerr " << lenerr << " "
             << msgerr << endl ;
        if ( sts != MPI_SUCCESS ) {
          ostringstream strstream ;
          strstream << "==========================================================="
                    << endl << "test" << myrank << " KO"
                    << "==========================================================="
                    << endl ;
          cout << strstream.str() << endl ;
          CPPUNIT_FAIL( strstream.str() ) ;
        }
        mpi_access.check() ;
     }

     if ( myrank != 0 ) {
       int iprobe ;
       for ( iprobe = 5 ; iprobe < 10 ; iprobe++ ) {
          cout << "test" << myrank << " ============================ iprobe "
               << iprobe << "============================" << endl ;
          int source, tag, outcount ;
          MPI_Datatype datatype ;
          int probeflag = false ;
          while ( !probeflag ) {
               sts = mpi_access.IProbe( target, source, tag, datatype, outcount,
                                        probeflag ) ;
               char msgerr[MPI_MAX_ERROR_STRING] ;
               int lenerr ;
               mpi_access.errorString(sts, msgerr, &lenerr) ;
               cout << "test" << myrank << " IProbe iprobe " << iprobe
                    << " target " << target << " probeflag " << probeflag
                    << " tag " << tag << " outcount " << outcount << " datatype "
                    << datatype << " lenerr " << lenerr << " " << msgerr << endl ;
               if ( sts != MPI_SUCCESS ) {
                 ostringstream strstream ;
                 strstream << "=========================================================="
                           << endl << "test" << myrank << " IProbe KO iprobe " << iprobe
                           << endl
                           << "=========================================================="
                           << endl ;
                 cout << strstream.str() << endl ;
                 CPPUNIT_FAIL( strstream.str() ) ;
               }
               if ( !probeflag ) {
                 //cout << "========================================================"
                 //     << endl << "test" << myrank << " IProbe KO(OK) iprobe " << iprobe
                 //     << " probeflag " << probeflag << endl
                 //     << "========================================================"
                 //     << endl ;
               }
               else {
                 cout << "test" << myrank << " " << iprobe << " IProbe target "
                      << target << " source " << source << " tag " << tag
                      << " outcount " << outcount << " probeflag " << probeflag
                      << endl ;
                 if ( datatype != MPI_DOUBLE ) {
                   ostringstream strstream ;
                   strstream << "========================================================"
                             << endl << "test" << myrank << " MPI_DOUBLE KO" << endl
                             << "========================================================"
                             << endl ;
                   cout << strstream.str() << endl ;
                   CPPUNIT_FAIL( strstream.str() ) ;
                 }
                 else {
                   int flag ;
                   sts = mpi_access.cancel( source, tag, datatype, outcount, flag ) ;
                   if ( sts != MPI_SUCCESS || !flag ) {
                     mpi_access.errorString(sts, msgerr, &lenerr) ;
                     cout << "======================================================"
                          << endl << "test" << myrank << " lenerr " << lenerr << " "
                          << msgerr << endl << "test" << myrank
                          << " Cancel PendingIrecv KO flag " << flag << " iprobe "
                          << iprobe << " Irecv completed" << endl
                          << "======================================================"
                          << endl ;
                     //return 1 ;
                   }
                   else {
                     cout << "======================================================"
                          << endl << "test" << myrank
                          << " Cancel PendingIrecv OK RequestId " << " flag "
                          << flag << " iprobe " << iprobe << endl
                          << "======================================================"
                          << endl ;
                   }
                 }
                 int Reqtarget, Reqtag, Reqerror, Reqoutcount ;
                 mpi_access.status( RequestId[iprobe], Reqtarget, Reqtag, Reqerror,
                                    Reqoutcount, true ) ;
                 cout << "test" << myrank << " Status Reqtarget "<< Reqtarget
                      << " Reqtag " << Reqtag << " Reqoutcount " << Reqoutcount
                      << endl ;
                 int Reqflag ;
                 sts = mpi_access.cancel( RequestId[iprobe] , Reqflag ) ;
                 cout << "test" << myrank << " " << iprobe
                      << " Cancel Irecv done Reqtarget " << Reqtarget
                      << " Reqtag " << Reqtag << " Reqoutcount " << Reqoutcount
                      << " Reqflag " << Reqflag << endl ;
                 if ( sts != MPI_SUCCESS || !Reqflag ) {
                   mpi_access.errorString(sts, msgerr, &lenerr) ;
                   ostringstream strstream ;
                   strstream << "========================================================"
                             << endl << "test" << myrank << " lenerr " << lenerr << " "
                             << msgerr << endl << "test" << myrank
                             << " Cancel Irecv KO Reqflag " << Reqflag << " iprobe "
                             << iprobe << endl
                             << "========================================================"
                             << endl ;
                   cout << strstream.str() << endl ;
                   CPPUNIT_FAIL( strstream.str() ) ;
                 }
                 else {
                   cout << "========================================================"
                        << endl << "test" << myrank
                        << " Cancel Irecv OK RequestId " << RequestId[iprobe]
                        << " Reqflag " << Reqflag << " iprobe " << iprobe << endl
                        << "========================================================"
                        << endl ;
                   probeflag = Reqflag ;
                 }
               }
          }
       }
     }
     mpi_access.waitAll(10,RequestId) ;
     mpi_access.deleteRequests(10,RequestId) ;
  }

  int source, tag, outcount, flag ;
  MPI_Datatype datatype ;
  sts = mpi_access.IProbe(target, source, tag, datatype, outcount, flag ) ;
  char msgerr[MPI_MAX_ERROR_STRING] ;
  int lenerr ;
  mpi_access.errorString(sts, msgerr, &lenerr) ;
  cout << "test" << myrank << " lenerr " << lenerr << " "
       << msgerr << endl ;
  if ( sts != MPI_SUCCESS || flag ) {
    ostringstream strstream ;
    strstream << "==========================================================="
              << endl << "test" << myrank << " IProbe KO flag " << flag
              << " remaining unread/cancelled message :" << endl
              << " source " << source << " tag " << tag << endl
              << "==========================================================="
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  mpi_access.testAll(10,RequestId,flag) ;
  mpi_access.waitAll(10,RequestId) ;
  mpi_access.deleteRequests(10,RequestId) ;
  mpi_access.testAll(10,RequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  mpi_access.check() ;

  mpi_access.barrier() ;

  delete group ;

//  MPI_Finalize();

  cout << "test" << myrank << " OK" << endl ;

  return ;
}




