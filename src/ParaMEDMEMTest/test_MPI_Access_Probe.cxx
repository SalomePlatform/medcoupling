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

void MPIAccessTest::test_MPI_Access_Probe() {

  cout << "test_MPI_Access_Probe" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 ) {
    cout << "test_MPI_Access_Probe must be runned with 2 procs" << endl ;
    CPPUNIT_FAIL("test_MPI_Access_Probe must be runned with 2 procs") ;
  }

  cout << "test_MPI_Access_Probe" << myrank << endl ;

  ParaMEDMEM::CommInterface interface ;

  ParaMEDMEM::MPIProcessorGroup* group = new ParaMEDMEM::MPIProcessorGroup(interface) ;

  ParaMEDMEM::MPIAccess mpi_access( group ) ;

  if ( myrank >= 2 ) {
    mpi_access.barrier() ;
    delete group ;
    return ;
  }

  int target = 1 - myrank ;
  int RequestId[10] ;
  int sts ;
  int i ;
  for ( i = 0 ; i < 10 ; i++ ) {
     if ( myrank == 0 ) {
       sts = mpi_access.send(&i,1,MPI_INT,target, RequestId[i]) ;
       cout << "test" << myrank << " Send RequestId " << RequestId[i]
            << endl ;
     }
     else {
       int source, tag, outcount ;
       MPI_Datatype datatype ;
       sts = mpi_access.probe(target, source, tag, datatype, outcount ) ;
       cout << "test" << myrank << " Probe target " << target << " source " << source
            << " tag " << tag << " outcount " << outcount << endl ;
       int recvbuf ;
       sts = mpi_access.recv(&recvbuf,outcount,datatype,source, RequestId[i],
                             &outcount) ;
       if ( (outcount != 1) | (recvbuf != i) ) {
         ostringstream strstream ;
         strstream << "==========================================================="
                   << "test" << myrank << " outcount " << outcount
                   << " recvbuf " << recvbuf << " KO"
                   << "==========================================================="
                   << endl ;
         cout << strstream.str() << endl ;
         CPPUNIT_FAIL( strstream.str() ) ;
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
     mpi_access.check() ;
  }
  int flag ;
  mpi_access.testAll(10,RequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  mpi_access.waitAll(10,RequestId) ;
  mpi_access.check() ;

  mpi_access.barrier() ;

  delete group ;

//  MPI_Finalize();

  cout << "test" << myrank << " OK" << endl ;

  return ;
}




