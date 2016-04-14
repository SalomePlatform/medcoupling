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

void MPIAccessTest::test_MPI_Access_Cyclic_Send_Recv() {

  debugStream << "test_MPI_Access_Cyclic_Send_Recv" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 3 ) {
      cerr << "test_MPI_Access_Send_Recv must be runned with 3 procs" << endl ;
    //CPPUNIT_FAIL("test_MPI_Access_Send_Recv must be runned with 3 procs") ;
    return;
  }

  debugStream << "test_MPI_Access_Cyclic_Send_Recv" << myrank << endl ;

  MEDCoupling::CommInterface interface ;

  MEDCoupling::MPIProcessorGroup* group = new MEDCoupling::MPIProcessorGroup(interface) ;

  MEDCoupling::MPIAccess mpi_access( group ) ;

  if ( myrank >= 3 ) {
    mpi_access.barrier() ;
    delete group ;
    return ;
  }

  int alltarget[3] = {1 , 2 , 0 } ;
  int allsource[3] = {2 , 0 , 1 } ;
  int RequestId[10] ;
  int sts ;
  int i = 0 ;
  if ( myrank == 0 ) {
    sts = mpi_access.send(&i,1,MPI_INT,alltarget[myrank], RequestId[i]) ;
    debugStream << "test" << myrank << " Send RequestId " << RequestId[i]
         << " tag " << mpi_access.sendMPITag(alltarget[myrank]) << endl ;
  }
  for ( i = 0 ; i < 10 ; i++ ) {
     int recvbuf ;
     int outcount ;
     if ( i & 1 ) {
       outcount = 0 ;
       sts = mpi_access.recv(&recvbuf,1,MPI_INT,allsource[myrank], RequestId[i],
                             &outcount) ;
     }
     else {
       sts = mpi_access.recv(&recvbuf,1,MPI_INT,allsource[myrank], RequestId[i]) ;
       outcount = 1 ;
     }
     //int source, tag, error, outcount ;
     //mpi_access.Status( RequestId[i], source, tag, error, outcount, true) ;
     debugStream << "test" << myrank << " Recv RequestId " << RequestId[i]
          << " tag " << mpi_access.recvMPITag(allsource[myrank])
          << " outcount " << outcount << endl ;
     if ( (outcount != 1) | (recvbuf != i) ) {
       ostringstream strstream ;
       strstream << "==========================================================="
                 << "test" << myrank << " outcount "
                 << outcount << " recvbuf " << recvbuf << " KO"
                 << "==========================================================="
                 << endl ;
       debugStream << strstream.str() << endl ;
       CPPUNIT_FAIL( strstream.str() ) ;
     }
     if ( myrank == 0 ) {
       if ( i != 9 ) {
         int ii = i + 1 ;
         sts = mpi_access.send(&ii,1,MPI_INT,alltarget[myrank], RequestId[i]) ;
         debugStream << "test" << myrank << " Send RequestId " << RequestId[i]
              << " tag " << mpi_access.sendMPITag(alltarget[myrank]) << endl ;
       }
     }
     else {
       sts = mpi_access.send(&i,1,MPI_INT,alltarget[myrank], RequestId[i]) ;
       debugStream << "test" << myrank << " Send RequestId " << RequestId[i]
            << " tag " << mpi_access.sendMPITag(alltarget[myrank]) << endl ;
     }
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
     if(MPI_ACCESS_VERBOSE) mpi_access.check() ;
  }

  int flag ;
  mpi_access.testAll(10,RequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  mpi_access.waitAll(10,RequestId) ;
  if(MPI_ACCESS_VERBOSE) mpi_access.check() ;

  int sendrequests[10] ;
  int sendreqsize = mpi_access.sendRequestIds( alltarget[myrank] , 10 ,
                                               sendrequests ) ;
  if ( sendreqsize != 0 ) {
    ostringstream strstream ;
    strstream << "=========================================================" << endl
              << "test" << myrank << " sendreqsize " << sendreqsize << " KO" << endl
              << "=========================================================" << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  int recvrequests[10] ;
  int recvreqsize = mpi_access.sendRequestIds( allsource[myrank] , 10 ,
                                               recvrequests ) ;
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




