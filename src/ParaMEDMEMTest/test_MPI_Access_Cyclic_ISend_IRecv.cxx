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

void MPIAccessTest::test_MPI_Access_Cyclic_ISend_IRecv() {

  debugStream << "test_MPI_Access_Cyclic_ISend_IRecv" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 3 ) {
      cerr << "test_MPI_Access_Cyclic_ISend_IRecv must be runned with 3 procs" << endl ;
    //CPPUNIT_FAIL("test_MPI_Access_Cyclic_ISend_IRecv must be runned with 3 procs") ;
    return;
  }

  debugStream << "test_MPI_Access_Cyclic_ISend_IRecv" << myrank << endl ;

  MEDCoupling::CommInterface interface ;

  MEDCoupling::MPIProcessorGroup* group = new MEDCoupling::MPIProcessorGroup(interface) ;

  MEDCoupling::MPIAccess mpi_access( group ) ;

#define maxsend 100

  if ( myrank >= 3 ) {
    mpi_access.barrier() ;
    delete group ;
    return ;
  }

  int alltarget[3] = {1 , 2 , 0 } ;
  int allsource[3] = {2 , 0 , 1 } ;
  int SendRequestId[maxsend] ;
  int RecvRequestId[maxsend] ;
  int sendbuf[maxsend] ;
  int recvbuf[maxsend] ;
  int sts ;
  int i = 0 ;
  if ( myrank == 0 ) {
    sendbuf[i] = i ;
    sts = mpi_access.ISend(&sendbuf[i],1,MPI_INT,alltarget[myrank],
                           SendRequestId[i]) ;
    debugStream << "test" << myrank << " Send RequestId " << SendRequestId[i]
         << " tag " << mpi_access.sendMPITag(alltarget[myrank]) << endl ;
  }
  for ( i = 0 ; i < maxsend ; i++ ) {
     recvbuf[i] = -1 ;
     sts = mpi_access.IRecv(&recvbuf[i],1,MPI_INT,allsource[myrank],
                            RecvRequestId[i]) ;
     debugStream << "test" << myrank << " Recv RequestId " << RecvRequestId[i]
          << " tag " << mpi_access.recvMPITag(allsource[myrank]) << endl ;
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
            strstream << "====================================================="
                      << endl << "test" << myrank << " outcount "
                      << outcount << " recvbuf[ " << j << " ] " << recvbuf[j] << " KO"
                      << endl << "====================================================="
                      << endl ;
            debugStream << strstream.str() << endl ;
            CPPUNIT_FAIL( strstream.str() ) ;
          }
        }
     }
     if ( myrank == 0 ) {
       if ( i != maxsend-1 ) {
         sendbuf[i+1] = i + 1 ;
         sts = mpi_access.ISend(&sendbuf[i+1],1,MPI_INT,alltarget[myrank],
                                SendRequestId[i+1]) ;
         debugStream << "test" << myrank << " Send RequestId " << SendRequestId[i+1]
              << " tag " << mpi_access.sendMPITag(alltarget[myrank]) << endl ;
       }
     }
     else {
       sendbuf[i] = i ;
       sts = mpi_access.ISend(&sendbuf[i],1,MPI_INT,alltarget[myrank],
                              SendRequestId[i]) ;
       debugStream << "test" << myrank << " Send RequestId " << SendRequestId[i]
            << " tag " << mpi_access.sendMPITag(alltarget[myrank]) << endl ;
     }
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
     if(MPI_ACCESS_VERBOSE) mpi_access.check() ;
  }

  int flag ;
  mpi_access.testAll(maxsend,SendRequestId,flag) ;
  mpi_access.testAll(maxsend,RecvRequestId,flag) ;
  mpi_access.waitAll(maxsend,SendRequestId) ;
  mpi_access.deleteRequests(maxsend,SendRequestId) ;
  mpi_access.waitAll(maxsend,RecvRequestId) ;
  mpi_access.deleteRequests(maxsend,RecvRequestId) ;
  if(MPI_ACCESS_VERBOSE) mpi_access.check() ;
  mpi_access.testAll(maxsend,SendRequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "=========================================================" << endl
              << "test" << myrank << " TestAllSendflag " << flag << " KO" << endl
              << "=========================================================" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    debugStream << "=========================================================" << endl
         << "test" << myrank << " TestAllSendflag " << flag << " OK" << endl
         << "=========================================================" << endl ;
  }
  mpi_access.testAll(maxsend,RecvRequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "=========================================================" << endl
              << "test" << myrank << " TestAllRecvflag " << flag << " KO" << endl
              << "=========================================================" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    debugStream << "=========================================================" << endl
         << "test" << myrank << " TestAllRecvflag " << flag << " OK" << endl
         << "=========================================================" << endl ;
  }

  int sendrequests[maxsend] ;
  int sendreqsize = mpi_access.sendRequestIds( alltarget[myrank] , maxsend ,
                                               sendrequests ) ;
  if ( sendreqsize != 0 ) {
    ostringstream strstream ;
    strstream << "=========================================================" << endl
              << "test" << myrank << " sendreqsize " << sendreqsize << " KO" << endl
              << "=========================================================" << endl ;
    debugStream << strstream.str() << endl ;
    int source, tag, error, outcount ;
    mpi_access.status(sendrequests[0], source, tag, error, outcount, true) ;
    debugStream << "test" << myrank << " RequestId " << sendrequests[0]
         << " source " << source << " tag " << tag << " error " << error
         << " outcount " << outcount << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  else {
    debugStream << "=========================================================" << endl
         << "test" << myrank << " sendreqsize " << sendreqsize << " OK" << endl
         << "=========================================================" << endl ;
  }
  int recvrequests[maxsend] ;
  int recvreqsize = mpi_access.sendRequestIds( allsource[myrank] , maxsend ,
                                               recvrequests ) ;
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




