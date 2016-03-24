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

void MPIAccessTest::test_MPI_Access_SendRecv() {

  debugStream << "MPIAccessTest::test_MPI_Access_SendRecv" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 ) {
      cerr << "MPIAccessTest::test_MPI_Access_SendRecv must be runned with 2 procs" << endl ;
    //CPPUNIT_FAIL("test_MPI_Access_SendRecv must be runned with 2 procs") ;
    return;
  }

  debugStream << "MPIAccessTest::test_MPI_Access_SendRecv" << myrank << endl ;

  MEDCoupling::CommInterface interface ;

  MEDCoupling::MPIProcessorGroup* group = new MEDCoupling::MPIProcessorGroup(interface) ;

  MEDCoupling::MPIAccess mpi_access( group ) ;

  if ( myrank >= 2 ) {
    mpi_access.barrier() ;
    delete group ;
    return ;
  }

  int target = 1 - myrank ;
  int sendRequestId[10] ;
  int recvRequestId[10] ;
  int sts ;
  int i ;
  for ( i = 0 ; i < 10 ; i++ ) {
     int recvbuf ;
     int outcount ;
     if ( i & 1 ) {
       outcount = -1 ;
       sts = mpi_access.sendRecv(&i,1,MPI_INT,target, sendRequestId[i],
                                 &recvbuf,1,MPI_INT,target, recvRequestId[i],
                                 &outcount) ;
     }
     else {
       sts = mpi_access.sendRecv(&i,1,MPI_INT,target, sendRequestId[i],
                                 &recvbuf,1,MPI_INT,target, recvRequestId[i]) ;
//       outcount = mpi_access.MPIOutCount( recvRequestId[i] ) ;
       outcount = 1 ;
     }
     debugStream << "test" << myrank << " Send sendRequestId " << sendRequestId[i]
          << " tag " << mpi_access.sendMPITag(target)
          << " recvRequestId " << recvRequestId[i]
          << " tag " << mpi_access.recvMPITag(target)
          << " outcount " << outcount << " MPIOutCount "
          << mpi_access.MPIOutCount( recvRequestId[i] ) << endl ;
     if ( (outcount != 1) | (recvbuf != i) ) {
       ostringstream strstream ;
       strstream << "==========================================================="
                 << "test" << myrank << " outcount " << outcount
                 << " recvbuf " << recvbuf << " KO"
                 << "==========================================================="
                 << endl ;
       debugStream << strstream.str() << endl ;
       CPPUNIT_FAIL( strstream.str() ) ;
     }
     char msgerr[MPI_MAX_ERROR_STRING] ;
     int lenerr ;
     mpi_access.errorString(sts, msgerr, &lenerr) ;
     debugStream << "test" << myrank << " lenerr " << lenerr << " "
          << msgerr << endl ;

     if ( sts != MPI_SUCCESS ) {
       ostringstream strstream ;
       strstream << "==========================================================="
                 << "test" << myrank << " KO"
                 << "==========================================================="
                 << endl ;
       debugStream << strstream.str() << endl ;
       CPPUNIT_FAIL( strstream.str() ) ;
     }
     if(MPI_ACCESS_VERBOSE) mpi_access.check() ;
  }

  int flag ;
  mpi_access.testAll(10,sendRequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  mpi_access.waitAll(10,sendRequestId) ;
  mpi_access.testAll(10,recvRequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  mpi_access.waitAll(10,recvRequestId) ;
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

  mpi_access.barrier() ;

  delete group ;

//  MPI_Finalize();

  debugStream << "test" << myrank << " OK" << endl ;

  return ;
}




