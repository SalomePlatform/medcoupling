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

void MPIAccessTest::test_MPI_Access_ISend_IRecv_BottleNeck() {

  cout << "test_MPI_Access_ISend_IRecv_BottleNeck" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 ) {
    ostringstream strstream ;
    strstream << "test_MPI_Access_ISend_IRecv_BottleNeck must be runned with 2 procs"
              << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  cout << "test_MPI_Access_ISend_IRecv_BottleNeck" << myrank << endl ;

  ParaMEDMEM::CommInterface interface ;

  ParaMEDMEM::MPIProcessorGroup* group = new ParaMEDMEM::MPIProcessorGroup(interface) ;

  ParaMEDMEM::MPIAccess mpi_access( group ) ;

#define maxreq 10000

  if ( myrank >= 2 ) {
    mpi_access.barrier() ;
    delete group ;
    return ;
  }

  int target = 1 - myrank ;
  int SendRequestId[maxreq] ;
  int RecvRequestId[maxreq] ;
  int sts ;
  int sendbuf[maxreq] ;
  int recvbuf[maxreq] ;
  int i ;
  for ( i = 0 ; i < maxreq ; i++ ) {
     if ( myrank == 0 ) {
       sendbuf[i] = i ;
       sts = mpi_access.ISend(sendbuf,i,MPI_INT,target, SendRequestId[i]) ;
       cout << "test" << myrank << " ISend RequestId " << SendRequestId[i]
            << " tag " << mpi_access.sendMPITag(target) << endl ;
     }
     else {
       //sleep( 1 ) ;
       sts = mpi_access.IRecv(recvbuf,i,MPI_INT,target, RecvRequestId[i]) ;
       cout << "test" << myrank << " IRecv RequestId " << RecvRequestId[i]
            << " tag " << mpi_access.recvMPITag(target) << endl ;
       int recvreqsize = mpi_access.recvRequestIdsSize() ;
       int * recvrequests = new int[ recvreqsize ] ;
       recvreqsize = mpi_access.recvRequestIds( target , recvreqsize , recvrequests ) ;
       int j ;
       for (j = 0 ; j < recvreqsize ; j++) {
          int flag ;
          mpi_access.test( recvrequests[j], flag ) ;
          if ( flag ) {
            int source, tag, error, outcount ;
            mpi_access.status( recvrequests[j], source, tag, error, outcount,
                               true ) ;
            cout << "test" << myrank << " Test(Recv RequestId "
                 << recvrequests[j] << ") : source " << source << " tag " << tag
                 << " error " << error << " outcount " << outcount
                 << " flag " << flag << " : DeleteRequest" << endl ;
            mpi_access.deleteRequest( recvrequests[j] ) ;
          }
          else {
//            cout << "test" << myrank << " Test(Recv RequestId "
//                 << recvrequests[j] << ") flag " << flag << endl ;
          }
       }
       delete [] recvrequests ;
     }
     if ( sts != MPI_SUCCESS ) {
       char msgerr[MPI_MAX_ERROR_STRING] ;
       int lenerr ;
       mpi_access.errorString(sts, msgerr, &lenerr) ;
       cout << "test" << myrank << " lenerr " << lenerr << " "
            << msgerr << endl ;
     }

     if ( sts != MPI_SUCCESS ) {
       ostringstream strstream ;
       strstream << "==========================================================="
                 << "test" << myrank << " KO"
                 << "==========================================================="
                 << endl ;
       cout << strstream.str() << endl ;
       CPPUNIT_FAIL( strstream.str() ) ;
     }
  }

  mpi_access.check() ;
  if ( myrank == 0 ) {
    int size = mpi_access.sendRequestIdsSize() ;
    cout << "test" << myrank << " before WaitAll sendreqsize " << size << endl ;
    mpi_access.waitAll(maxreq, SendRequestId) ;
    size = mpi_access.sendRequestIdsSize() ;
    cout << "test" << myrank << " after WaitAll sendreqsize " << size << endl ;
    int * ArrayOfSendRequests = new int[ size ] ;
    int nSendRequest = mpi_access.sendRequestIds( size , ArrayOfSendRequests ) ;
    int i ;
    for ( i = 0 ; i < nSendRequest ; i++ ) {
       mpi_access.deleteRequest( ArrayOfSendRequests[i] ) ;
    }
    delete [] ArrayOfSendRequests ;
  }
  else {
    int size = mpi_access.recvRequestIdsSize() ;
    cout << "test" << myrank << " before WaitAll recvreqsize " << size << endl ;
    mpi_access.waitAll(maxreq, RecvRequestId) ;
    size = mpi_access.recvRequestIdsSize() ;
    cout << "test" << myrank << " after WaitAll recvreqsize " << size << endl ;
    int * ArrayOfRecvRequests = new int[ size ] ;
    int nRecvRequest = mpi_access.recvRequestIds( size , ArrayOfRecvRequests ) ;
    int i ;
    for ( i = 0 ; i < nRecvRequest ; i++ ) {
       mpi_access.deleteRequest( ArrayOfRecvRequests[i] ) ;
    }
    delete [] ArrayOfRecvRequests ;
  }
  mpi_access.check() ;

  if ( myrank == 0 ) {
    int sendrequests[maxreq] ;
    int sendreqsize = mpi_access.sendRequestIds( target , maxreq , sendrequests ) ;
    int i ;
    if ( sendreqsize != 0 ) {
      ostringstream strstream ;
      strstream << "=========================================================" << endl
                << "test" << myrank << " sendreqsize " << sendreqsize << " KO" << endl
                << "=========================================================" << endl ;
      cout << strstream.str() << endl ;
      for ( i = 0 ; i < sendreqsize ; i++ ) {
         cout << "test" << myrank << " sendrequests[ " << i << " ] = "
              << sendrequests[i] << endl ;
      }
      CPPUNIT_FAIL( strstream.str() ) ;
    }
    else {
      cout << "=========================================================" << endl
           << "test" << myrank << " sendreqsize " << sendreqsize << " OK" << endl
           << "=========================================================" << endl ;
    }
  }
  else {
    int recvrequests[maxreq] ;
    int recvreqsize = mpi_access.recvRequestIds( target , maxreq , recvrequests ) ;
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

  mpi_access.barrier() ;

  delete group ;

//  MPI_Finalize();

  cout << "test" << myrank << " OK" << endl ;

  return ;
}




