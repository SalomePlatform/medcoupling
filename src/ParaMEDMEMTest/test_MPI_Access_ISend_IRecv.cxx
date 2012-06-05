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

void MPIAccessTest::test_MPI_Access_ISend_IRecv() {

  cout << "test_MPI_Access_ISend_IRecv" << endl ;

  //  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 ) {
    cout << "test_MPI_Access_ISend_IRecv must be runned with 2 procs" << endl ;
    CPPUNIT_FAIL("test_MPI_Access_ISend_IRecv must be runned with 2 procs") ;
  }

  cout << "test_MPI_Access_ISend_IRecv" << myrank << endl ;

  ParaMEDMEM::CommInterface interface ;

  ParaMEDMEM::MPIProcessorGroup* group = new ParaMEDMEM::MPIProcessorGroup(interface) ;

  ParaMEDMEM::MPIAccess mpi_access( group ) ;

#define maxreq 100

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
      sts = mpi_access.ISend(&sendbuf[i],1,MPI_INT,target, SendRequestId[i]) ;
      cout << "test" << myrank << " ISend RequestId " << SendRequestId[i]
           << " tag " << mpi_access.sendMPITag(target) << endl ;
    }
    else {
      sts = mpi_access.IRecv(&recvbuf[i],1,MPI_INT,target, RecvRequestId[i]) ;
      cout << "test" << myrank << " IRecv RequestId " << RecvRequestId[i]
           << " tag " << mpi_access.recvMPITag(target) << endl ;
    }
    int j ;
    for (j = 0 ; j <= i ; j++) {
      int flag ;
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
          //else {
          //  cout << "==========================================================="
          //       << endl << "test" << myrank << " outcount " << outcount
          //       << " RequestId " << RecvRequestId[j] << " recvbuf "
          //       << recvbuf[j] << " OK" << endl
          //       << "==========================================================="
          //       << endl ;
          //}
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
  }

  mpi_access.check() ;
  if ( myrank == 0 ) {
    mpi_access.waitAll(maxreq, SendRequestId) ;
    mpi_access.deleteRequests(maxreq, SendRequestId) ;
  }
  else {
    mpi_access.waitAll(maxreq, RecvRequestId) ;
    mpi_access.deleteRequests(maxreq, RecvRequestId) ;
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
    int recvreqsize = mpi_access.sendRequestIds( target , maxreq , recvrequests ) ;
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




