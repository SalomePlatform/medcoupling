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

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <mpi.h>

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

void MPIAccessTest::test_MPI_Access_Time() {

  cout << "test_MPI_Access_Time" << endl ;

  //  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 ) {
    ostringstream strstream ;
    strstream << "test_MPI_Access_Time must be runned with 2 procs" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  cout << "test_MPI_Access_Time" << myrank << endl ;

  ParaMEDMEM::CommInterface interface ;

  ParaMEDMEM::MPIProcessorGroup* group = new ParaMEDMEM::MPIProcessorGroup(interface) ;

  ParaMEDMEM::MPIAccess mpi_access( group ) ;

#define maxreq 10

  if ( myrank >= 2 ) {
    cout << "test_MPI_Access_Time_0 rank" << myrank << " --> mpi_access->Barrier" << endl ;
    mpi_access.barrier() ;
    cout << "test_MPI_Access_Time_0 rank" << myrank << " <-- mpi_access->Barrier" << endl ;
    delete group ;
    cout << "test_MPI_Access_Time" << myrank << " OK" << endl ;
    return ;
  }

  int target = 1 - myrank ;
  int SendTimeRequestId[maxreq] ;
  int RecvTimeRequestId[maxreq] ;
  int SendRequestId[maxreq] ;
  int RecvRequestId[maxreq] ;
  int sts ;
  int sendbuf[maxreq] ;
  int recvbuf[maxreq] ;
  int i = 0 ;
  ParaMEDMEM::TimeMessage aSendTimeMsg[maxreq] ;
  ParaMEDMEM::TimeMessage aRecvTimeMsg[maxreq] ;
  double t ;
  double dt = 1. ;
  double maxt = 10. ;
  for ( t = 0 ; t < maxt ; t = t+dt ) {
    if ( myrank == 0 ) {
      aSendTimeMsg[i].time = t ;
      aSendTimeMsg[i].deltatime = dt ;
      //aSendTimeMsg[i].maxtime = maxt ;
      //sts = mpi_access.ISend( &aSendTimeMsg , mpi_access.timeExtent() ,
      sts = mpi_access.ISend( &aSendTimeMsg[i] , 1 ,
                              mpi_access.timeType() , target ,
                              SendTimeRequestId[i]) ;
      cout << "test" << myrank << " ISend RequestId " << SendTimeRequestId[i]
           << " tag " << mpi_access.sendMPITag(target) << endl ;
      sendbuf[i] = i ;
      sts = mpi_access.ISend(&sendbuf[i],1,MPI_INT,target, SendRequestId[i]) ;
      cout << "test" << myrank << " ISend RequestId " << SendRequestId[i]
           << " tag " << mpi_access.sendMPITag(target) << endl ;
    }
    else {
      //sts = mpi_access.IRecv( &aRecvTimeMsg , mpi_access.timeExtent() ,
      sts = mpi_access.IRecv( &aRecvTimeMsg[i] , 1 ,
                              mpi_access.timeType() , target ,
                              RecvTimeRequestId[i]) ;
      cout << "test" << myrank << " IRecv RequestId " << RecvTimeRequestId[i]
           << " tag " << mpi_access.recvMPITag(target) << endl ;
      sts = mpi_access.IRecv(&recvbuf[i],1,MPI_INT,target, RecvRequestId[i]) ;
      cout << "test" << myrank << " IRecv RequestId " << RecvRequestId[i]
           << " tag " << mpi_access.recvMPITag(target) << endl ;
    }
    int j ;
    for (j = 0 ; j <= i ; j++) {
      int flag ;
      if ( myrank == 0 ) {
        mpi_access.test( SendTimeRequestId[j], flag ) ;
      }
      else {
        mpi_access.test( RecvTimeRequestId[j], flag ) ;
      }
      if ( flag ) {
        int target,source, tag, error, outcount ;
        if ( myrank == 0 ) {
          mpi_access.status( SendTimeRequestId[j], target, tag, error, outcount,
                             true ) ;
          cout << "test" << myrank << " Test(Send TimeRequestId " << SendTimeRequestId[j]
               << ") : target " << target << " tag " << tag << " error " << error
               << " flag " << flag << aSendTimeMsg[j] << endl ;
        }
        else {
          mpi_access.status( RecvTimeRequestId[j], source, tag, error, outcount,
                             true ) ;
          cout << "test" << myrank << " Test(Recv TimeRequestId "
               << RecvTimeRequestId[j] << ") : source " << source << " tag " << tag
               << " error " << error << " outcount " << outcount
               << " flag " << flag << aRecvTimeMsg[j] << endl ;
          if ( (outcount != 1) | (aRecvTimeMsg[j].time != j) ) {
            ostringstream strstream ;
            strstream << "==========================================================="
                      << endl << "test" << myrank << " outcount " << outcount << " KO"
                      << " RecvTimeRequestId " << RecvTimeRequestId[j] << endl
                      << "==========================================================="
                      << endl ;
            cout << strstream.str() << endl ;
            CPPUNIT_FAIL( strstream.str() ) ;
          }
          else {
            cout << "==========================================================="
                 << endl << "test" << myrank << " outcount " << outcount
                 << " RecvTimeRequestId " << RecvTimeRequestId[j] << " OK" << endl
                 << "==========================================================="
                 << endl ;
          }
        }
      }
      if ( myrank == 0 ) {
        mpi_access.test( SendRequestId[j], flag ) ;
      }
      else {
        mpi_access.test( RecvRequestId[j], flag ) ;
      }
      if ( flag ) {
        int target,source, tag, error, outcount ;
        if ( myrank == 0 ) {
          mpi_access.status( SendRequestId[j], target, tag, error, outcount,
                             true ) ;
          cout << "test" << myrank << " Test(Send RequestId " << SendRequestId[j]
               << ") : target " << target << " tag " << tag << " error " << error
               << " flag " << flag << endl ;
        }
        else {
          mpi_access.status( RecvRequestId[j], source, tag, error, outcount,
                             true ) ;
          cout << "test" << myrank << " Test(Recv RequestId "
               << RecvRequestId[j] << ") : source " << source << " tag " << tag
               << " error " << error << " outcount " << outcount
               << " flag " << flag << endl ;
          if ( (outcount != 1) | (recvbuf[j] != j) ) {
            ostringstream strstream ;
            strstream << "==========================================================="
                      << endl << "test" << myrank << " outcount "
                      << outcount << " recvbuf " << recvbuf[j] << " KO" << endl
                      << "==========================================================="
                      << endl ;
            cout << strstream.str() << endl ;
            CPPUNIT_FAIL( strstream.str() ) ;
          }
          else {
            cout << "==========================================================="
                 << endl << "test" << myrank << " outcount " << outcount
                 << " RequestId " << RecvRequestId[j] << " OK" << endl
                 << "==========================================================="
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
                << "test" << myrank << " KO"
                << "==========================================================="
                << endl ;
      cout << strstream.str() << endl ;
      CPPUNIT_FAIL( strstream.str() ) ;
    }
    i = i + 1 ;
  }

  mpi_access.check() ;
  if ( myrank == 0 ) {
    mpi_access.waitAll(maxreq, SendTimeRequestId) ;
    mpi_access.deleteRequests(maxreq, SendTimeRequestId) ;
    mpi_access.waitAll(maxreq, SendRequestId) ;
    mpi_access.deleteRequests(maxreq, SendRequestId) ;
  }
  else {
    mpi_access.waitAll(maxreq, RecvTimeRequestId) ;
    mpi_access.deleteRequests(maxreq, RecvTimeRequestId) ;
    mpi_access.waitAll(maxreq, RecvRequestId) ;
    mpi_access.deleteRequests(maxreq, RecvRequestId) ;
  }
  mpi_access.check() ;

  if ( myrank == 0 ) {
    int sendrequests[2*maxreq] ;
    int sendreqsize = mpi_access.sendRequestIds( target , 2*maxreq , sendrequests ) ;
    if ( sendreqsize != 0 ) {
      ostringstream strstream ;
      strstream << "=========================================================" << endl
                << "test" << myrank << " sendreqsize " << sendreqsize << " KO" << endl
                << "=========================================================" << endl ;
      cout << strstream.str() << endl ;
      CPPUNIT_FAIL( strstream.str() ) ;
    }
    else {
      cout << "=========================================================" << endl
           << "test" << myrank << " sendreqsize " << sendreqsize << " OK" << endl
           << "=========================================================" << endl ;
    }
  }
  else {
    int recvrequests[2*maxreq] ;
    int recvreqsize = mpi_access.sendRequestIds( target , 2*maxreq , recvrequests ) ;
    if ( recvreqsize != 0 ) {
      ostringstream strstream ;
      strstream << "=========================================================" << endl
                << "test" << myrank << " recvreqsize " << recvreqsize << " KO" << endl
                << "=========================================================" << endl ;
      cout << strstream.str() << endl ;
      CPPUNIT_FAIL( strstream.str() ) ;
    }
    else {
      cout << "=========================================================" << endl
           << "test" << myrank << " recvreqsize " << recvreqsize << " OK" << endl
           << "=========================================================" << endl ;
    }
  }

  cout << "test_MPI_Access_Time_0 rank" << myrank << " --> mpi_access->Barrier" << endl ;
  mpi_access.barrier() ;
  cout << "test_MPI_Access_Time_0 rank" << myrank << " <-- mpi_access->Barrier" << endl ;

  delete group ;

  //  MPI_Finalize();

  cout << "test_MPI_Access_Time" << myrank << " OK" << endl ;

  return ;
}




