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

void chksts( int sts , int myrank , MEDCoupling::MPIAccess * mpi_access ) {
  char msgerr[MPI_MAX_ERROR_STRING] ;
  int lenerr ;
  if ( sts != MPI_SUCCESS ) {
    mpi_access->errorString(sts, msgerr, &lenerr) ;
    debugStream << "test" << myrank << " lenerr " << lenerr << " "
         << msgerr << endl ;
    ostringstream strstream ;
    strstream << "==========================================================="
              << "test" << myrank << " KO"
              << "==========================================================="
              << endl ;
    debugStream << strstream.str() << endl ;
    CPPUNIT_FAIL( strstream.str() ) ;
  }
return ;
}

void MPIAccessTest::test_MPI_Access_Time_0() {

  debugStream << "test_MPI_Access_Time_0" << endl ;

//  MPI_Init(&argc, &argv) ; 

  int size ;
  int myrank ;
  MPI_Comm_size(MPI_COMM_WORLD,&size) ;
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank) ;

  if ( size < 2 ) {
    ostringstream strstream ;
    strstream << "usage :" << endl
              << "mpirun -np <nbprocs> test_MPI_Access_Time_0" <<endl
              << " nbprocs =2" << endl
              << "test must be runned with 2 procs" << endl ;
    cerr << strstream.str() << endl ;
    //CPPUNIT_FAIL( strstream.str() ) ;
    return;
  }

#define maxreq 100

  double t ;
  double dt[2] = {2., 1.} ;
  double maxt = maxreq/dt[myrank] ;

  debugStream << "test_MPI_Access_Time_0 rank" << myrank << endl ;

  MEDCoupling::CommInterface interface ;

  MEDCoupling::MPIProcessorGroup* group = new MEDCoupling::MPIProcessorGroup(interface) ;

  MEDCoupling::MPIAccess * mpi_access = new MEDCoupling::MPIAccess( group ) ;

  if ( myrank >= 2 ) {
    debugStream << "test_MPI_Access_Time_0 rank" << myrank << " --> mpi_access->barrier" << endl ;
    mpi_access->barrier() ;
    debugStream << "test_MPI_Access_Time_0 rank" << myrank << " <-- mpi_access->barrier" << endl ;
    debugStream << "test_MPI_Access_Time_0 rank" << myrank << " --> mpi_access->barrier" << endl ;
    mpi_access->barrier() ;
    debugStream << "test_MPI_Access_Time_0 rank" << myrank << " <-- mpi_access->barrier" << endl ;
    delete group ;
    delete mpi_access ;
    debugStream << "test_MPI_Access_Time" << myrank << " OK" << endl ;
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
  MEDCoupling::TimeMessage aSendTimeMsg[maxreq] ;
  int lasttime = -1 ;
  MEDCoupling::TimeMessage RecvTimeMessages[maxreq+1] ;
  MEDCoupling::TimeMessage *aRecvTimeMsg = &RecvTimeMessages[1] ;
//  mpi_access->Trace() ;
  int istep = 0 ;
  for ( t = 0 ; t < maxt ; t = t+dt[myrank] ) {
     debugStream << "test" << myrank << " ==========================TIME " << t
          << " ==========================" << endl ;
     if ( myrank == 0 ) {
       aSendTimeMsg[istep].time = t ;
       aSendTimeMsg[istep].deltatime = dt[myrank] ;
       //aSendTimeMsg[istep].maxtime = maxt ;
       if ( t+dt[myrank] >= maxt ) {
         aSendTimeMsg[istep].deltatime = 0 ;
       }
       sts = mpi_access->ISend( &aSendTimeMsg[istep] , 1 ,
                               mpi_access->timeType() , target ,
                               SendTimeRequestId[istep]) ;
       debugStream << "test" << myrank << " ISend TimeRequestId " << SendTimeRequestId[istep]
            << " tag " << mpi_access->MPITag(SendTimeRequestId[istep]) << endl ;
       chksts( sts , myrank , mpi_access ) ;
       sendbuf[istep] = istep ;
       sts = mpi_access->ISend(&sendbuf[istep],1,MPI_INT,target, SendRequestId[istep]) ;
       debugStream << "test" << myrank << " ISend Data RequestId " << SendRequestId[istep]
            << " tag " << mpi_access->MPITag(SendRequestId[istep]) << endl ;
       chksts( sts , myrank , mpi_access ) ;
//CheckSent
//=========
       int sendrequests[2*maxreq] ;
       int sendreqsize = mpi_access->sendRequestIds( target , 2*maxreq ,
                                                    sendrequests ) ;
       int j , flag ;
       for ( j = 0 ; j < sendreqsize ; j++ ) {
          sts = mpi_access->test( sendrequests[j] , flag ) ;
          chksts( sts , myrank , mpi_access ) ;
          if ( flag ) {
            mpi_access->deleteRequest( sendrequests[j] ) ;
            debugStream << "test" << myrank << " " << j << ". " << sendrequests[j]
                 << " sendrequest deleted" << endl ;
          }
       }
     }
     else {
//InitRecv
//========
       if ( t == 0 ) {
         aRecvTimeMsg[lasttime].time = 0 ;
         sts = mpi_access->IRecv( &aRecvTimeMsg[lasttime+1] , 1 ,
                                 mpi_access->timeType() ,
                                 target , RecvTimeRequestId[lasttime+1]) ;
         debugStream << "test" << myrank << " t == 0 IRecv TimeRequestId "
              << RecvTimeRequestId[lasttime+1]
              << " MPITag " << mpi_access->MPITag( RecvTimeRequestId[lasttime+1] )
              << " MPICompleted "
              << mpi_access->MPICompleted( RecvTimeRequestId[lasttime+1] ) << endl ;
         chksts( sts , myrank , mpi_access ) ;
       }
       else {
         debugStream << "test" << myrank << " t # 0 lasttime " << lasttime << endl ;
//InitialOutTime
//==============
         bool outtime = false ;
         if ( lasttime != -1 ) {
           if ( t <= aRecvTimeMsg[lasttime-1].time ) {
             ostringstream strstream ;
             strstream << "==========================================================="
                       << endl << "test" << myrank << " t " << t << " <= "
                       << "aRecvTimeMsg[ " << lasttime << "-1 ].time "
                       << aRecvTimeMsg[lasttime-1].time << " KO" << endl
                       << "==========================================================="
                       << endl ;
             debugStream << strstream.str() << endl ;
             CPPUNIT_FAIL( strstream.str() ) ;
           }
           else {
             debugStream << "==========================================================="
                  << endl << "test" << myrank << " t " << t << " > "
                  << "aRecvTimeMsg[ " << lasttime << "-1 ].time "
                  << aRecvTimeMsg[lasttime-1].time << " OK" << endl
                  << "==========================================================="
                  << endl ;
           }
           //outtime = ((aRecvTimeMsg[lasttime].time +
           //            aRecvTimeMsg[lasttime].deltatime) >=
           //           aRecvTimeMsg[lasttime].maxtime) ;
           outtime = aRecvTimeMsg[lasttime].deltatime == 0 ;
         }
// CheckRecv - CheckTime
// On a lasttime tel que :
// aRecvTimeMsg[ lasttime-1 ].time < T(i-1) <= aRecvTimeMsg[ lasttime ].time
// On cherche lasttime tel que :
// aRecvTimeMsg[ lasttime-1 ].time < T(i) <= aRecvTimeMsg[ lasttime ].time
         if ( t <= aRecvTimeMsg[lasttime].time ) {
           outtime = false ;
         }
         debugStream << "test" << myrank << " while outtime( " << outtime << " && t " << t
              << " > aRecvTimeMsg[ " << lasttime << " ] "
              << aRecvTimeMsg[lasttime].time << " )" << endl ;
         while ( !outtime && (t > aRecvTimeMsg[lasttime].time) ) {
              lasttime += 1 ;
//TimeMessage
//===========
              sts = mpi_access->wait( RecvTimeRequestId[lasttime] ) ;
              chksts( sts , myrank , mpi_access ) ;
              debugStream << "test" << myrank << " Wait done RecvTimeRequestId "
                   << RecvTimeRequestId[lasttime] << " lasttime " << lasttime
                   << " tag " << mpi_access->MPITag(RecvTimeRequestId[lasttime])
                   << aRecvTimeMsg[lasttime] << endl ;
              if ( lasttime == 0 ) {
                aRecvTimeMsg[lasttime-1] = aRecvTimeMsg[lasttime] ;
              }
              mpi_access->deleteRequest( RecvTimeRequestId[lasttime] ) ;

              double deltatime = aRecvTimeMsg[lasttime].deltatime ;
              //double maxtime = aRecvTimeMsg[lasttime].maxtime ;
              double nexttime = aRecvTimeMsg[lasttime].time + deltatime ;
              debugStream << "test" << myrank << " t " << t << " lasttime " << lasttime
                   << " deltatime " << deltatime
                   << " nexttime " << nexttime << endl ;
              //if ( nexttime < maxtime && t > nexttime ) {
              if ( deltatime != 0 && t > nexttime ) {
//CheckRecv :
//=========   
                //while ( nexttime < maxtime && t > nexttime ) {
                while ( deltatime != 0 && t > nexttime ) {
                     int source, MPITag, outcount ;
                     MPI_Datatype datatype ;
                     sts = mpi_access->probe( target , source, MPITag, datatype,
                                             outcount ) ;
                     chksts( sts , myrank , mpi_access ) ;
// Cancel DataMessages jusqu'a un TimeMessage
                     int cancelflag ;
                     while ( !mpi_access->isTimeMessage( MPITag ) ) {
                          sts = mpi_access->cancel( source, MPITag, datatype, outcount ,
                          //sts = mpi_access->cancel( source, datatype, outcount ,
                                                   //RecvRequestId[lasttime] ,
                                                   cancelflag ) ;
                          debugStream << "test" << myrank << " Recv TO CANCEL RequestId "
                               << RecvRequestId[lasttime]
                               << " tag " << mpi_access->recvMPITag( target )
                               << " cancelflag " << cancelflag << endl ;
                          chksts( sts , myrank , mpi_access ) ;
                          sts = mpi_access->probe( target , source, MPITag, datatype,
                                                  outcount ) ;
                          chksts( sts , myrank , mpi_access ) ;
                     }
//On peut avancer en temps
                     nexttime += deltatime ;
                     //if ( nexttime < maxtime && t > nexttime ) {
                     if ( deltatime != 0 && t > nexttime ) {
// Cancel du TimeMessage
                       sts = mpi_access->cancel( source, MPITag, datatype, outcount ,
                       //sts = mpi_access->cancel( source, datatype, outcount ,
                                                //RecvRequestId[lasttime] ,
                                                cancelflag ) ;
                       debugStream << "test" << myrank << " Time TO CANCEL RequestId "
                            << RecvRequestId[lasttime]
                            << " tag " << mpi_access->recvMPITag( target )
                            << " cancelflag " << cancelflag << endl ;
                       chksts( sts , myrank , mpi_access ) ;
                     }
                }
              }
              else {
//DoRecv
//======
                debugStream << "test" << myrank << " Recv target " << target
                     << " lasttime " << lasttime
                     << " lasttime-1 " << aRecvTimeMsg[lasttime-1]
                     << " lasttime " << aRecvTimeMsg[lasttime]
                     << endl ;
                sts = mpi_access->recv(&recvbuf[lasttime],1,MPI_INT,target,
                                       RecvRequestId[lasttime]) ;
                debugStream << "test" << myrank << " Recv RequestId "
                     << RecvRequestId[lasttime]
                     << " tag " << mpi_access->recvMPITag( target )
                     << endl ;
                chksts( sts , myrank , mpi_access ) ;
              }
              //outtime = ((aRecvTimeMsg[lasttime].time +
              //            aRecvTimeMsg[lasttime].deltatime) >=
              //           aRecvTimeMsg[lasttime].maxtime) ;
              outtime = aRecvTimeMsg[lasttime].deltatime == 0 ;
              if ( !outtime ) {
// Une lecture asynchrone d'un message temps a l'avance
                sts = mpi_access->IRecv( &aRecvTimeMsg[lasttime+1] , 1 ,
                                        mpi_access->timeType() , target ,
                                        RecvTimeRequestId[lasttime+1]) ;
                debugStream << "test" << myrank << " IRecv TimeRequestId "
                     << RecvTimeRequestId[lasttime+1] << " MPITag "
                     << mpi_access->MPITag( RecvTimeRequestId[lasttime+1] )
                     << " MPICompleted "
                     << mpi_access->MPICompleted( RecvTimeRequestId[lasttime+1] )
                     << endl ;
                chksts( sts , myrank , mpi_access ) ;
              }
              else if ( t <= aRecvTimeMsg[lasttime].time ) {
                outtime = false ;
              }
         }
         
         //printf("DEBUG t %.15f Msg[lasttime-1] %.15f Msg[lasttime] %.15f \n",t,
         //       aRecvTimeMsg[lasttime-1].time,aRecvTimeMsg[lasttime].time) ;
         if ( ((t <= aRecvTimeMsg[lasttime-1].time) ||
               (t > aRecvTimeMsg[lasttime].time)) && !outtime ) {
           ostringstream strstream ;
           strstream << "==========================================================="
                     << endl << "test" << myrank << " t " << t << " <= "
                     << "aRecvTimeMsg[ " << lasttime << "-1 ].time "
                     << aRecvTimeMsg[lasttime-1].time << " ou t " << t << " > "
                     << "aRecvTimeMsg[ " << lasttime << " ].time "
                     << aRecvTimeMsg[lasttime].time << endl
                     << " ou bien outtime " << outtime << " KO RequestTimeIds "
                     << RecvTimeRequestId[lasttime-1] << " " << RecvTimeRequestId[lasttime]
                     << " RequestIds "
                     << RecvRequestId[lasttime-1] << " " << RecvRequestId[lasttime] << endl
                     << "==========================================================="
                     << endl ;
           debugStream << strstream.str() << endl ;
           CPPUNIT_FAIL( strstream.str() ) ;
         }
         else {
           debugStream << "==========================================================="
                << endl << "test" << myrank 
                << " aRecvTimeMsg[ " << lasttime << "-1 ].time "
                << aRecvTimeMsg[lasttime-1].time << " < t " << t << " <= "
                << "aRecvTimeMsg[ " << lasttime << " ].time "
                << aRecvTimeMsg[lasttime].time << endl
                << " ou bien outtime " << outtime << " OK RequestTimeIds "
                << RecvTimeRequestId[lasttime-1] << " " << RecvTimeRequestId[lasttime]
                << " RequestIds "
                << RecvRequestId[lasttime-1] << " " << RecvRequestId[lasttime] << endl
                << "==========================================================="
                << endl ;
         }
       }
     }
     chksts( sts , myrank , mpi_access ) ;
     istep = istep + 1 ;
  }

  debugStream << "test" << myrank << " Barrier :" << endl ;
  mpi_access->barrier() ;

  if (MPI_ACCESS_VERBOSE) mpi_access->check() ;

  if ( myrank == 0 ) {
//CheckFinalSent
//==============
    debugStream << "test" << myrank << " CheckFinalSent :" << endl ;
    int sendrequests[2*maxreq] ;
    int sendreqsize = mpi_access->sendRequestIds( target , 2*maxreq , sendrequests ) ;
    int j ;
    for ( j = 0 ; j < sendreqsize ; j++ ) {
       sts = mpi_access->wait( sendrequests[j] ) ;
       chksts( sts , myrank , mpi_access ) ;
       mpi_access->deleteRequest( sendrequests[j] ) ;
       debugStream << "test" << myrank << " " << j << ". " << sendrequests[j] << " deleted"
            << endl ;
    }
  }
  else {
    debugStream << "test" << myrank << " CheckFinalRecv :" << endl ;
    int recvrequests[2*maxreq] ;
    int recvreqsize = mpi_access->recvRequestIds( target , 2*maxreq , recvrequests ) ;
    int cancelflag ;
    int j ;
    for ( j = 0 ; j < recvreqsize ; j++ ) {
       sts = mpi_access->cancel( recvrequests[j] , cancelflag ) ;
       chksts( sts , myrank , mpi_access ) ;
       mpi_access->deleteRequest( recvrequests[j] ) ;
       debugStream << "test" << myrank << " " << j << ". " << recvrequests[j] << " deleted"
            << " cancelflag " << cancelflag << endl ;
    }
    int source, MPITag, outcount , flag ;
    MPI_Datatype datatype ;
    sts = mpi_access->IProbe( target , source, MPITag, datatype,
                             outcount , flag ) ;
    chksts( sts , myrank , mpi_access ) ;
    while ( flag ) {
         sts = mpi_access->cancel( source, MPITag, datatype, outcount ,
         //sts = mpi_access->cancel( source, datatype, outcount ,
                                  //RecvRequestId[lasttime] ,
                                  cancelflag ) ;
         debugStream << "test" << myrank << " TO CANCEL RequestId "
              << RecvRequestId[lasttime]
              << " tag " << mpi_access->recvMPITag( target )
              << " cancelflag " << cancelflag << endl ;
         chksts( sts , myrank , mpi_access ) ;
         sts = mpi_access->IProbe( target , source, MPITag, datatype,
                                  outcount , flag ) ;
         chksts( sts , myrank , mpi_access ) ;
    }
  }
  if(MPI_ACCESS_VERBOSE) mpi_access->check() ;

  if ( myrank == 0 ) {
    int sendrequests[2*maxreq] ;
    int sendreqsize = mpi_access->sendRequestIds( target , 2*maxreq , sendrequests ) ;
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
  }
  else {
    int recvrequests[2*maxreq] ;
    int recvreqsize = mpi_access->recvRequestIds( target , 2*maxreq , recvrequests ) ;
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
  }

  int i ;
  for ( i = 0 ; i <= lasttime ; i++ ) {
     debugStream << "test" << myrank << " " << i << ". RecvTimeMsg "
          << aRecvTimeMsg[i].time << " recvbuf " << recvbuf[i] << endl ;
  }

  debugStream << "test_MPI_Access_Time_0 rank" << myrank << " --> mpi_access->barrier" << endl ;
  mpi_access->barrier() ;
  debugStream << "test_MPI_Access_Time_0 rank" << myrank << " <-- mpi_access->barrier" << endl ;

  delete group ;
  delete mpi_access ;

//  MPI_Finalize();

  debugStream << "test" << myrank << " OK" << endl ;

  return ;
}




