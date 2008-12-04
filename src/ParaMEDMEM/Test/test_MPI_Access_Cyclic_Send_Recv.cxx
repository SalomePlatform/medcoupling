//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
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
#include "MPI_Access.hxx"

// use this define to enable lines, execution of which leads to Segmentation Fault
#define ENABLE_FAULTS

// use this define to enable CPPUNIT asserts and fails, showing bugs
#define ENABLE_FORCED_FAILURES

using namespace std;
using namespace ParaMEDMEM;

void MPIAccessTest::test_MPI_Access_Cyclic_Send_Recv() {

  cout << "test_MPI_Access_Cyclic_Send_Recv" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 3 ) {
    cout << "test_MPI_Access_Send_Recv must be runned with 3 procs" << endl ;
    CPPUNIT_FAIL("test_MPI_Access_Send_Recv must be runned with 3 procs") ;
  }

  cout << "test_MPI_Access_Cyclic_Send_Recv" << myrank << endl ;

  ParaMEDMEM::CommInterface interface ;

  ParaMEDMEM::MPIProcessorGroup* group = new ParaMEDMEM::MPIProcessorGroup(interface) ;

  ParaMEDMEM::MPI_Access mpi_access( group ) ;

  if ( myrank >= 3 ) {
    mpi_access.Barrier() ;
    delete group ;
    return ;
  }

  int alltarget[3] = {1 , 2 , 0 } ;
  int allsource[3] = {2 , 0 , 1 } ;
  int RequestId[10] ;
  int sts ;
  int i = 0 ;
  if ( myrank == 0 ) {
    sts = mpi_access.Send(&i,1,MPI_INT,alltarget[myrank], RequestId[i]) ;
    cout << "test" << myrank << " Send RequestId " << RequestId[i]
         << " tag " << mpi_access.SendMPITag(alltarget[myrank]) << endl ;
  }
  for ( i = 0 ; i < 10 ; i++ ) {
     int recvbuf ;
     int outcount ;
     if ( i & 1 ) {
       outcount = 0 ;
       sts = mpi_access.Recv(&recvbuf,1,MPI_INT,allsource[myrank], RequestId[i],
                             &outcount) ;
     }
     else {
       sts = mpi_access.Recv(&recvbuf,1,MPI_INT,allsource[myrank], RequestId[i]) ;
       outcount = 1 ;
     }
     //int source, tag, error, outcount ;
     //mpi_access.Status( RequestId[i], source, tag, error, outcount, true) ;
     cout << "test" << myrank << " Recv RequestId " << RequestId[i]
          << " tag " << mpi_access.RecvMPITag(allsource[myrank])
          << " outcount " << outcount << endl ;
     if ( (outcount != 1) | (recvbuf != i) ) {
       ostringstream strstream ;
       strstream << "==========================================================="
                 << "test" << myrank << " outcount "
                 << outcount << " recvbuf " << recvbuf << " KO"
                 << "==========================================================="
                 << endl ;
       cout << strstream.str() << endl ;
       CPPUNIT_FAIL( strstream.str() ) ;
     }
     if ( myrank == 0 ) {
       if ( i != 9 ) {
         int ii = i + 1 ;
         sts = mpi_access.Send(&ii,1,MPI_INT,alltarget[myrank], RequestId[i]) ;
         cout << "test" << myrank << " Send RequestId " << RequestId[i]
              << " tag " << mpi_access.SendMPITag(alltarget[myrank]) << endl ;
       }
     }
     else {
       sts = mpi_access.Send(&i,1,MPI_INT,alltarget[myrank], RequestId[i]) ;
       cout << "test" << myrank << " Send RequestId " << RequestId[i]
            << " tag " << mpi_access.SendMPITag(alltarget[myrank]) << endl ;
     }
     char msgerr[MPI_MAX_ERROR_STRING] ;
     int lenerr ;
     mpi_access.Error_String(sts, msgerr, &lenerr) ;
     cout << "test" << myrank << " lenerr " << lenerr
          << " " << msgerr << endl ;

     if ( sts != MPI_SUCCESS ) {
       ostringstream strstream ;
       strstream << "==========================================================="
                 << "test" << myrank << " KO"
                 << "==========================================================="
                 << endl ;
       cout << strstream.str() << endl ;
       CPPUNIT_FAIL( strstream.str() ) ;
     }
     mpi_access.Check() ;
  }

  int flag ;
  mpi_access.TestAll(10,RequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  mpi_access.WaitAll(10,RequestId) ;
  mpi_access.Check() ;

  int sendrequests[10] ;
  int sendreqsize = mpi_access.SendRequestIds( alltarget[myrank] , 10 ,
                                               sendrequests ) ;
  if ( sendreqsize != 0 ) {
    ostringstream strstream ;
    strstream << "=========================================================" << endl
              << "test" << myrank << " sendreqsize " << sendreqsize << " KO" << endl
              << "=========================================================" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  int recvrequests[10] ;
  int recvreqsize = mpi_access.SendRequestIds( allsource[myrank] , 10 ,
                                               recvrequests ) ;
  if ( recvreqsize != 0 ) {
    ostringstream strstream ;
    strstream << "=========================================================" << endl
              << "test" << myrank << " recvreqsize " << recvreqsize << " KO" << endl
              << "=========================================================" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  mpi_access.Barrier() ;

  delete group ;

//  MPI_Finalize();

  cout << "test" << myrank << " OK" << endl ;

  return ;
}




