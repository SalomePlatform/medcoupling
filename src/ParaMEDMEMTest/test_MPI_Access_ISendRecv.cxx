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
using namespace MEDCoupling;

void MPIAccessTest::test_MPI_Access_ISendRecv() {

  debugStream << "test_MPI_Access_ISendRecv" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 ) {
      cerr << "test_MPI_Access_ISendRecv must be runned with 2 procs" << endl ;
    //CPPUNIT_FAIL("test_MPI_Access_ISendRecv must be runned with 2 procs") ;
    return;
  }

  debugStream << "test_MPI_Access_ISendRecv" << myrank << endl ;

  MEDCoupling::CommInterface interface ;

  MEDCoupling::MPIProcessorGroup* group = new MEDCoupling::MPIProcessorGroup(interface) ;

  MEDCoupling::MPIAccess mpi_access( group ) ;

  if ( myrank >= 2 ) {
    mpi_access.barrier() ;
    delete group ;
    return ;
  }

  int target = 1 - myrank ;
  int SendRequestId[10] ;
  int RecvRequestId[10] ;
  int sendbuf[10] ;
  int recvbuf[10] ;
  int sts ;
  int i ;
  for ( i = 0 ; i < 10 ; i++ ) {
     sendbuf[i] = i ;
     sts = mpi_access.ISendRecv(&sendbuf[i],1,MPI_INT,target, SendRequestId[i],
                                &recvbuf[i],1,MPI_INT,target, RecvRequestId[i]) ;
     debugStream << "test" << myrank << " Send sendRequestId " << SendRequestId[i]
          << " tag " << mpi_access.sendMPITag(target)
          << " recvRequestId " << RecvRequestId[i]
          << " tag " << mpi_access.recvMPITag(target) << endl ;
     char msgerr[MPI_MAX_ERROR_STRING] ;
     int lenerr ;
     mpi_access.errorString(sts, msgerr, &lenerr) ;
     debugStream << "test" << myrank << " lenerr " << lenerr
          << " " << msgerr << endl ;

     if ( sts != MPI_SUCCESS ) {
       ostringstream strstream ;
       strstream << "==========================================================="
                 << "test" << myrank << " KO"
                 << "==========================================================="
                 << endl ;
       debugStream << strstream.str() << endl ;
       CPPUNIT_FAIL( strstream.str() ) ;
     }
     int j ;
     for (j = 0 ; j <= i ; j++) {
        int flag ;
        if ( j < i ) {
          debugStream << "test" << myrank << " " << j << " -> Test-Send("<< SendRequestId[j]
               << ")" << endl ;
          mpi_access.test( SendRequestId[j], flag ) ;
          if ( flag ) {
            int target, tag, error, outcount ;
              mpi_access.status( SendRequestId[j], target, tag, error, outcount,
                                 true ) ;
              debugStream << "test" << myrank << " Send RequestId " << SendRequestId[j]
                   << " target " << target << " tag " << tag << " error " << error
                   << endl ;
            mpi_access.deleteRequest( SendRequestId[j] ) ;
          }
        }
        debugStream << "test" << myrank << " " << j << " -> Test-Recv("<< SendRequestId[j]
             << ")" << endl ;
        mpi_access.test( RecvRequestId[j], flag ) ;
        if ( flag ) {
          int source, tag, error, outcount ;
          mpi_access.status( RecvRequestId[j], source, tag, error, outcount,
                             true ) ;
          debugStream << "test" << myrank << " Recv RequestId" << j << " "
               << RecvRequestId[j] << " source " << source << " tag " << tag
               << " error " << error << " outcount " << outcount << endl ;
          if ( (outcount != 1) | (recvbuf[j] != j) ) {
             ostringstream strstream ;
             strstream << "==========================================================="
                       << "test" << myrank << " outcount "
                       << outcount << " recvbuf[ " << j << " ] " << recvbuf[j] << " KO"
                       << "==========================================================="
                       << endl ;
            debugStream << strstream.str() << endl ;
            CPPUNIT_FAIL( strstream.str() ) ;
          }
        }
     }
     mpi_access.errorString(sts, msgerr, &lenerr) ;
     debugStream << "test" << myrank << " lenerr " << lenerr << " "
          << msgerr << endl ;
     if(MPI_ACCESS_VERBOSE) mpi_access.check() ;
  }

  int flag ;
  mpi_access.testAll(10,SendRequestId,flag) ;
  mpi_access.waitAll(10,SendRequestId) ;
  mpi_access.deleteRequests(10,SendRequestId) ;
  mpi_access.testAll(10,SendRequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  mpi_access.testAll(10,RecvRequestId,flag) ;
  mpi_access.waitAll(10,RecvRequestId) ;
  mpi_access.deleteRequests(10,RecvRequestId) ;
  mpi_access.testAll(10,RecvRequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  if(MPI_ACCESS_VERBOSE) mpi_access.check() ;

  int sendrequests[10] ;
  int sendreqsize = mpi_access.sendRequestIds( target , 10 , sendrequests ) ;
  if ( sendreqsize != 0 ) {
    ostringstream strstream ;
    strstream << "=========================================================" << endl
              << "test" << myrank << " sendreqsize " << sendreqsize << " KO" << endl
              << "=========================================================" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    debugStream << "=========================================================" << endl
         << "test" << myrank << " sendreqsize " << sendreqsize << " OK" << endl
         << "=========================================================" << endl ;
  }
  int recvrequests[10] ;
  int recvreqsize = mpi_access.sendRequestIds( target , 10 , recvrequests ) ;
  if ( recvreqsize != 0 ) {
    ostringstream strstream ;
    strstream << "=========================================================" << endl
              << "test" << myrank << " recvreqsize " << recvreqsize << " KO" << endl
              << "=========================================================" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    debugStream << "=========================================================" << endl
         << "test" << myrank << " recvreqsize " << recvreqsize << " OK" << endl
         << "=========================================================" << endl ;
  }

  mpi_access.barrier() ;

  delete group ;

//  MPI_Finalize();

  debugStream << "test" << myrank << " OK" << endl ;

  return ;
}




