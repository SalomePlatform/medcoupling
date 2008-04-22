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

void MPIAccessTest::test_MPI_Access_ISendRecv() {

  cout << "test_MPI_Access_ISendRecv" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 ) {
    cout << "test_MPI_Access_ISendRecv must be runned with 2 procs" << endl ;
    CPPUNIT_FAIL("test_MPI_Access_ISendRecv must be runned with 2 procs") ;
  }

  cout << "test_MPI_Access_ISendRecv" << myrank << endl ;

  ParaMEDMEM::CommInterface interface ;

  ParaMEDMEM::MPIProcessorGroup* group = new ParaMEDMEM::MPIProcessorGroup(interface) ;

  ParaMEDMEM::MPI_Access mpi_access( group ) ;

  if ( myrank >= 2 ) {
    mpi_access.Barrier() ;
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
     cout << "test" << myrank << " Send sendRequestId " << SendRequestId[i]
          << " tag " << mpi_access.SendMPITag(target)
          << " recvRequestId " << RecvRequestId[i]
          << " tag " << mpi_access.RecvMPITag(target) << endl ;
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
     int j ;
     for (j = 0 ; j <= i ; j++) {
        int flag ;
        if ( j < i ) {
          cout << "test" << myrank << " " << j << " -> Test-Send("<< SendRequestId[j]
               << ")" << endl ;
          mpi_access.Test( SendRequestId[j], flag ) ;
          if ( flag ) {
            int target, tag, error, outcount ;
              mpi_access.Status( SendRequestId[j], target, tag, error, outcount,
                                 true ) ;
              cout << "test" << myrank << " Send RequestId " << SendRequestId[j]
                   << " target " << target << " tag " << tag << " error " << error
                   << endl ;
            mpi_access.DeleteRequest( SendRequestId[j] ) ;
          }
        }
        cout << "test" << myrank << " " << j << " -> Test-Recv("<< SendRequestId[j]
             << ")" << endl ;
        mpi_access.Test( RecvRequestId[j], flag ) ;
        if ( flag ) {
          int source, tag, error, outcount ;
          mpi_access.Status( RecvRequestId[j], source, tag, error, outcount,
                             true ) ;
          cout << "test" << myrank << " Recv RequestId" << j << " "
               << RecvRequestId[j] << " source " << source << " tag " << tag
               << " error " << error << " outcount " << outcount << endl ;
          if ( (outcount != 1) | (recvbuf[j] != j) ) {
             ostringstream strstream ;
             strstream << "==========================================================="
                       << "test" << myrank << " outcount "
                       << outcount << " recvbuf[ " << j << " ] " << recvbuf[j] << " KO"
                       << "==========================================================="
                       << endl ;
            cout << strstream.str() << endl ;
            CPPUNIT_FAIL( strstream.str() ) ;
          }
        }
     }
     mpi_access.Error_String(sts, msgerr, &lenerr) ;
     cout << "test" << myrank << " lenerr " << lenerr << " "
          << msgerr << endl ;
     mpi_access.Check() ;
  }

  int flag ;
  mpi_access.TestAll(10,SendRequestId,flag) ;
  mpi_access.WaitAll(10,SendRequestId) ;
  mpi_access.DeleteRequests(10,SendRequestId) ;
  mpi_access.TestAll(10,SendRequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }

  mpi_access.TestAll(10,RecvRequestId,flag) ;
  mpi_access.WaitAll(10,RecvRequestId) ;
  mpi_access.DeleteRequests(10,RecvRequestId) ;
  mpi_access.TestAll(10,RecvRequestId,flag) ;
  if ( !flag ) {
    ostringstream strstream ;
    strstream << "test" << myrank << " flag " << flag << " KO" << endl ;
    cout << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
  mpi_access.Check() ;

  int sendrequests[10] ;
  int sendreqsize = mpi_access.SendRequestIds( target , 10 , sendrequests ) ;
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
  int recvrequests[10] ;
  int recvreqsize = mpi_access.SendRequestIds( target , 10 , recvrequests ) ;
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

  mpi_access.Barrier() ;

  delete group ;

//  MPI_Finalize();

  cout << "test" << myrank << " OK" << endl ;

  return ;
}




