
#include "MPI_AccessDEC.hxx"

using namespace std;

namespace ParaMEDMEM {    

/*!
This constructor creates an MPI_AccessDEC which has \a local_group as a working side 
and  \a distant_group as an idle side. 
The constructor must be called synchronously on all processors of both processor groups.

\param local_group working side ProcessorGroup
\param distant_group lazy side ProcessorGroup
\param Asynchronous Communication mode (default asynchronous)
\param nStepBefore Number of Time step needed for the interpolation before current time
\param nStepAfter Number of Time step needed for the interpolation after current time

*/
//MPI_AccessDEC::MPI_AccessDEC( ProcessorGroup& local_group,
//                              ProcessorGroup& distant_group,
//                              TimeInterpolator * aTimeInterpolator,
//                              bool Asynchronous ) {
MPI_AccessDEC::MPI_AccessDEC( const ProcessorGroup& local_group,
                              const ProcessorGroup& distant_group,
                              bool Asynchronous ) {
  _local_group = &local_group ;
  _distant_group = &distant_group ;
  ProcessorGroup * union_group = _local_group->fuse( *_distant_group ) ;  
  int i ;
  std::set<int> procs;
  for ( i = 0 ; i < union_group->size() ; i++ ) {
     procs.insert(i) ;
  }
  MPIProcessorGroup *mpilg = (MPIProcessorGroup *)_local_group;
  _MPI_union_group = new ParaMEDMEM::MPIProcessorGroup( union_group->getCommInterface() , procs,
							mpilg->getWorldComm());
  delete union_group ;
  _MyRank = _MPI_union_group->myRank() ;
  _GroupSize = _MPI_union_group->size() ;
  _MPIAccess = new MPI_Access( _MPI_union_group ) ;
  _Asynchronous = Asynchronous ;
  _TimeMessages = new vector< vector< TimeMessage > > ;
  _TimeMessages->resize( _GroupSize ) ;
  _OutOfTime = new vector< bool > ;
  _OutOfTime->resize( _GroupSize ) ;
  _DataMessagesRecvCount = new vector< int > ;
  _DataMessagesRecvCount->resize( _GroupSize ) ;
  for ( i = 0 ; i < _GroupSize ; i++ ) {
     (*_OutOfTime)[i] = false ;
     (*_DataMessagesRecvCount)[i] = 0 ;
  }
  _DataMessagesType = new vector< MPI_Datatype > ;
  _DataMessagesType->resize( _GroupSize ) ;
  _DataMessages = new vector< vector< void * > > ;
  _DataMessages->resize( _GroupSize ) ;
  _TimeInterpolator = NULL ;
  _MapOfSendBuffers = new map< int , SendBuffStruct * > ;
  cout << "MPI_AccessDEC" << _MyRank << " Asynchronous " << _Asynchronous << endl ;
}

MPI_AccessDEC::~MPI_AccessDEC() {
  CheckFinalSent() ;
  CheckFinalRecv() ;
  delete _MPI_union_group ;
  delete _MPIAccess ;
  if ( _TimeInterpolator )
    delete _TimeInterpolator ;
  if ( _TimeMessages )
    delete _TimeMessages ;
  if ( _OutOfTime )
    delete _OutOfTime ;
  if ( _DataMessagesRecvCount )
    delete _DataMessagesRecvCount ;
  if ( _DataMessagesType )
    delete _DataMessagesType ;
  if ( _DataMessages )
    delete _DataMessages ;
  if ( _MapOfSendBuffers )
    delete _MapOfSendBuffers ;
} 

void MPI_AccessDEC::SetTimeInterpolator( TimeInterpolationMethod aTimeInterp ,
                                         double InterpPrecision, int nStepBefore,
                                         int nStepAfter ) {
  cout << "MPI_AccessDEC::SetTimeInterpolator" << _MyRank << " Asynchronous "
       << _Asynchronous << " TimeInterpolationMethod " << aTimeInterp
       << " InterpPrecision " << InterpPrecision << " nStepBefore " << nStepBefore
       << " nStepAfter " << nStepAfter << endl ;
  if ( _TimeInterpolator )
    delete _TimeInterpolator ;
	switch ( aTimeInterp ) {
    case WithoutTimeInterp :
      _TimeInterpolator = NULL ;
      _nStepBefore = 0 ;
      _nStepAfter = 0 ;
      break ;
    case LinearTimeInterp :
      _TimeInterpolator = new LinearTimeInterpolator( InterpPrecision , nStepBefore ,
                                                      nStepAfter ) ;
      _nStepBefore = nStepBefore ;
      _nStepAfter = nStepAfter ;
      int i ;
      for ( i = 0 ; i < _GroupSize ; i++ ) {
         (*_TimeMessages)[ i ].resize( _nStepBefore + _nStepAfter ) ;
         (*_DataMessages)[ i ].resize( _nStepBefore + _nStepAfter ) ;
         int j ;
         for ( j = 0 ; j < _nStepBefore + _nStepAfter ; j++ ) {
            (*_TimeMessages)[ i ][ j ].time = -1 ;
            (*_TimeMessages)[ i ][ j ].deltatime = -1 ;
            (*_DataMessages)[ i ][ j ] = NULL ;
         }
      }
      break ;
  }
}

/*!
Send sendcount datas from sendbuf[offset] with type sendtype to target of IntraCommunicator
(Internal Protected method)

Returns the request identifier SendRequestId

*/
int MPI_AccessDEC::Send( void* sendbuf, int sendcount , int offset ,
                         MPI_Datatype sendtype , int target , int &SendRequestId ) {
  int sts ;
  if ( _Asynchronous ) {
    if ( sendtype == MPI_INT ) {
      sts = _MPIAccess->ISend( &((int *) sendbuf)[offset] , sendcount , sendtype ,
                               target , SendRequestId ) ;
    }
    else {
      sts = _MPIAccess->ISend( &((double *) sendbuf)[offset] , sendcount , sendtype ,
                               target , SendRequestId ) ;
    }
  }
  else {
    if ( sendtype == MPI_INT ) {
      sts = _MPIAccess->Send( &((int *) sendbuf)[offset] , sendcount , sendtype ,
                              target , SendRequestId ) ;
    }
    else {
      sts = _MPIAccess->Send( &((double *) sendbuf)[offset] , sendcount , sendtype ,
                              target , SendRequestId ) ;
    }
  }
  return sts ;
}

/*!
Receive recvcount datas to recvbuf[offset] with type recvtype from target of IntraCommunicator
(Internal Protected method)

Returns the request identifier RecvRequestId

*/
int MPI_AccessDEC::Recv( void* recvbuf, int recvcount , int offset ,
                         MPI_Datatype recvtype , int target , int &RecvRequestId ) {
  int sts ;
  if ( _Asynchronous ) {
    if ( recvtype == MPI_INT ) {
      sts = _MPIAccess->IRecv( &((int *) recvbuf)[offset] , recvcount , recvtype ,
                               target , RecvRequestId ) ;
    }
    else {
      sts = _MPIAccess->IRecv( &((double *) recvbuf)[offset] , recvcount , recvtype ,
                               target , RecvRequestId ) ;
    }
  }
  else {
    if ( recvtype == MPI_INT ) {
      sts = _MPIAccess->Recv( &((int *) recvbuf)[offset] , recvcount , recvtype ,
                              target , RecvRequestId ) ;
    }
    else {
      sts = _MPIAccess->Recv( &((double *) recvbuf)[offset] , recvcount , recvtype ,
                              target , RecvRequestId ) ;
    }
  }
  return sts ;
}

/*!
Send sendcount datas from sendbuf[offset] with type sendtype to target of IntraCommunicator
Receive recvcount datas to recvbuf[offset] with type recvtype from target of IntraCommunicator
(Internal Protected method)

Returns the request identifier SendRequestId
Returns the request identifier RecvRequestId

*/
int MPI_AccessDEC::SendRecv( void* sendbuf, int sendcount , int sendoffset ,
                             MPI_Datatype sendtype ,
                             void* recvbuf, int recvcount , int recvoffset ,
                             MPI_Datatype recvtype , int target ,
                             int &SendRequestId , int &RecvRequestId ) {
  int sts ;
  if ( _Asynchronous ) {
    if ( sendtype == MPI_INT ) {
      if ( recvtype == MPI_INT ) {
        sts = _MPIAccess->ISendRecv( &((int *) sendbuf)[sendoffset] , sendcount ,
                                     sendtype , target , SendRequestId ,
                                     &((int *) recvbuf)[recvoffset] , recvcount ,
                                     recvtype , target , RecvRequestId ) ;
      }
      else {
        sts = _MPIAccess->ISendRecv( &((int *) sendbuf)[sendoffset] , sendcount ,
                                     sendtype , target , SendRequestId ,
                                     &((double *) recvbuf)[recvoffset] ,
                                     recvcount , recvtype , target , RecvRequestId ) ;
      }
    }
    else {
      if ( recvtype == MPI_INT ) {
        sts = _MPIAccess->ISendRecv( &((double *) sendbuf)[sendoffset] , sendcount ,
                                     sendtype , target , SendRequestId ,
                                     &((int *) recvbuf)[recvoffset] ,
                                     recvcount , recvtype , target , RecvRequestId ) ;
      }
      else {
        sts = _MPIAccess->ISendRecv( &((double *) sendbuf)[sendoffset] , sendcount ,
                                     sendtype , target , SendRequestId ,
                                     &((double *) recvbuf)[recvoffset] ,
                                     recvcount , recvtype , target , RecvRequestId ) ;
      }
    }
  }
  else {
    if ( sendtype == MPI_INT ) {
      if ( recvtype == MPI_INT ) {
        sts = _MPIAccess->SendRecv( &((int *) sendbuf)[sendoffset] , sendcount ,
                                    sendtype , target , SendRequestId ,
                                    &((int *) recvbuf)[recvoffset] , recvcount ,
                                    recvtype , target , RecvRequestId ) ;
      }
      else {
        sts = _MPIAccess->SendRecv( &((int *) sendbuf)[sendoffset] , sendcount ,
                                    sendtype , target , SendRequestId ,
                                    &((double *) recvbuf)[recvoffset] ,
                                    recvcount , recvtype , target , RecvRequestId ) ;
      }
    }
    else {
      if ( recvtype == MPI_INT ) {
        sts = _MPIAccess->SendRecv( &((double *) sendbuf)[sendoffset] , sendcount ,
                                    sendtype , target , SendRequestId ,
                                    &((int *) recvbuf)[recvoffset] ,
                                    recvcount , recvtype , target , RecvRequestId ) ;
      }
      else {
        cout << "SendRecv" << _MyRank << " target " << target << " sendbuf "
             << &((double *) sendbuf)[sendoffset] << " sendcount " << sendcount
             << " recvbuf " << &((double *) recvbuf)[recvoffset] << " recvcount "
             << recvcount << endl ;
        sts = _MPIAccess->SendRecv( &((double *) sendbuf)[sendoffset] , sendcount ,
                                    sendtype , target , SendRequestId ,
                                    &((double *) recvbuf)[recvoffset] ,
                                    recvcount , recvtype , target , RecvRequestId ) ;
      }
    }
  }
  return sts ;
}

/*!
Send sendcount datas from sendbuf[offset] with type sendtype to all targets of IntraCommunicator
Receive recvcount datas to recvbuf[offset] with type recvtype from all targets of IntraCommunicator

*/
int MPI_AccessDEC::AllToAll( void* sendbuf, int sendcount, MPI_Datatype sendtype ,
	                           void* recvbuf, int recvcount, MPI_Datatype recvtype ) {
  if ( _TimeInterpolator ) {
    return AllToAllTime( sendbuf, sendcount, sendtype , recvbuf, recvcount, recvtype ) ;
  }
  int sts ;
  int target ;
  int sendoffset = 0 ;
  int recvoffset = 0 ;
  int SendRequestId ;
  int RecvRequestId ;

//Free of SendBuffers 
  if ( _Asynchronous ) {
    CheckSent() ;
  }

//DoSend + DoRecv : SendRecv
  SendBuffStruct * aSendDataStruct = NULL ;
//  cout << "AllToAll" << _MyRank << " sendbuf " << sendbuf << " recvbuf " << recvbuf << endl ;
  if ( _Asynchronous && sendbuf ) {
    aSendDataStruct = new SendBuffStruct ;
    aSendDataStruct->SendBuffer = sendbuf ;
    aSendDataStruct->Counter = 0 ;
    aSendDataStruct->DataType = sendtype ;
  }
  for ( target = 0 ; target < _GroupSize ; target++ ) {
     sts = SendRecv( sendbuf , sendcount , sendoffset , sendtype ,
                     recvbuf , recvcount , recvoffset , recvtype ,
                     target , SendRequestId , RecvRequestId ) ;
     if ( _Asynchronous && sendbuf && sendcount ) {
       aSendDataStruct->Counter += 1 ;
       (*_MapOfSendBuffers)[ SendRequestId ] = aSendDataStruct ;
     }
     sendoffset += sendcount ;
     recvoffset += recvcount ;
  }
  if ( !_Asynchronous && sendbuf ) {
//    cout << "AllToAll" << _MyRank << " free of sendbuf " << sendbuf << endl ;
//    free( sendbuf ) ;
    if ( sendtype == MPI_INT ) {
      delete [] (int *) sendbuf ;
    }
    else {
      delete [] (double *) sendbuf ;
    }
  }
  return sts ;
}

/*!
Send sendcounts[target] datas from sendbuf[sdispls[target]] with type sendtype to all targets of IntraCommunicator
Receive recvcounts[target] datas to recvbuf[rdispls[target]] with type recvtype from all targets of IntraCommunicator

*/
int MPI_AccessDEC::AllToAllv( void* sendbuf, int* sendcounts, int* sdispls,
                              MPI_Datatype sendtype ,
	                            void* recvbuf, int* recvcounts, int* rdispls,
                              MPI_Datatype recvtype ) {
  if ( _TimeInterpolator ) {
    return AllToAllvTime( sendbuf, sendcounts, sdispls, sendtype ,
                          recvbuf, recvcounts, rdispls, recvtype ) ;
  }
  int sts ;
  int target ;
  int SendRequestId ;
  int RecvRequestId ;

//Free of SendBuffers 
  if ( _Asynchronous ) {
    CheckSent() ;
  }

//DoSend + DoRecv : SendRecv
  SendBuffStruct * aSendDataStruct = NULL ;
//  cout << "AllToAllv" << _MyRank << " sendbuf " << sendbuf << " recvbuf " << recvbuf
//       << endl ;
  if ( _Asynchronous && sendbuf ) {
    aSendDataStruct = new SendBuffStruct ;
    aSendDataStruct->SendBuffer = sendbuf ;
    aSendDataStruct->Counter = 0 ;
    aSendDataStruct->DataType = sendtype ;
  }
  for ( target = 0 ; target < _GroupSize ; target++ ) {
     if ( sendcounts[target] || recvcounts[target] ) {
       sts = SendRecv( sendbuf , sendcounts[target] , sdispls[target] , sendtype ,
                       recvbuf , recvcounts[target] , rdispls[target] , recvtype ,
                       target , SendRequestId , RecvRequestId ) ;
       if ( _Asynchronous && sendbuf && sendcounts[target]) {
         aSendDataStruct->Counter += 1 ;
         (*_MapOfSendBuffers)[ SendRequestId ] = aSendDataStruct ;
       }
     }
  }
  if ( !_Asynchronous && sendbuf ) {
//    cout << "AllToAllv" << _MyRank << " free of sendbuf " << sendbuf << endl ;
//    free( sendbuf ) ;
    if ( sendtype == MPI_INT ) {
      delete [] (int *) sendbuf ;
    }
    else {
      delete [] (double *) sendbuf ;
    }
  }
  return sts ;
}

/*
MPI_AccessDEC and the management of SendBuffers :
=================================================

. In the collective communications collectives we send only parts of
  the same buffer to each "target". So in asynchronous mode it is
  necessary that all parts are free before to delete/free the
  buffer.

. We assume that buffers are allocated with a new double[]. so a
  delete [] is done.

. The structure SendBuffStruct permit to keep the adress of the buffer
  and to manage a reference counter of that buffer. It contains
  also MPI_Datatype for the delete [] (double *) ... when the counter
  is null.

. The map _MapOfSendBuffers etablish the correspondance between each
  RequestId given by a MPI_Access->ISend(...) and a SendBuffStruct
  for each "target" of a part of the buffer.

. All that concerns only asynchronous Send. In synchronous mode,
  we delete senbuf just after the Send.
*/

/*
MPI_AccessDEC and the management of RecvBuffers :
=================================================

If there is no interpolation, no special action is done.

With interpolation for each target :
------------------------------------
. We have _TimeMessages[target] which is a vector of TimesMessages.
  We have 2 TimesMessages in our case with a linear interpolation.
  They contain the previous time(t0)/deltatime and the last
  time(t1)/deltatime.

. We have _DataMessages[target] which is a vector of DatasMessages.
  We have 2 DatasMessages in our case with a linear interpolation.
  They contain the previous datas at time(t0)/deltatime and at last
  time(t1)/deltatime.

. At time _t(t*) of current processus we do the interpolation of
  the values of the 2 DatasMessages which are returned in the part of
  recvbuf corresponding to the target with t0 < t* <= t1.

. Because of the difference of "deltatimes" between processes, we
  may have t0 < t1 < t* and there is an extrapolation.

. The vectors _OutOfTime, _DataMessagesRecvCount and _DataMessagesType
  contain for each target true if t* > last t1, recvcount and
  MPI_Datatype for the finalize of messages at the end.
*/

/*!
Send a TimeMessage to all targets of IntraCommunicator
Receive the TimeMessages from targets of IntraCommunicator if necessary.

Send sendcount datas from sendbuf[offset] with type sendtype to all targets of IntraCommunicator
Returns recvcount datas to recvbuf[offset] with type recvtype after an interpolation
with datas received from all targets of IntraCommunicator.

*/
int MPI_AccessDEC::AllToAllTime( void* sendbuf, int sendcount , MPI_Datatype sendtype ,
	                               void* recvbuf, int recvcount , MPI_Datatype recvtype ) {
  int sts ;
  int target ;
  int sendoffset = 0 ;
  int SendTimeRequestId ;
  int SendDataRequestId ;

  if ( _TimeInterpolator == NULL ) {
    return MPI_ERR_OTHER ;
  }

//Free of SendBuffers 
  if ( _Asynchronous ) {
    CheckSent() ;
  }

//DoSend : Time + SendBuff
  SendBuffStruct * aSendTimeStruct = NULL ;
  SendBuffStruct * aSendDataStruct = NULL ;
//  cout << "AllToAllTime" << _MyRank << " sendbuf " << sendbuf << " recvbuf " << recvbuf
//       << endl ;
  if ( sendbuf && sendcount ) {
    TimeMessage * aSendTimeMessage = new TimeMessage ;
    if ( _Asynchronous ) {
      aSendTimeStruct = new SendBuffStruct ;
      aSendTimeStruct->SendBuffer = aSendTimeMessage ;
      aSendTimeStruct->Counter = 0 ;
      aSendTimeStruct->DataType = _MPIAccess->TimeType() ;
      aSendDataStruct = new SendBuffStruct ;
      aSendDataStruct->SendBuffer = sendbuf ;
      aSendDataStruct->Counter = 0 ;
      aSendDataStruct->DataType = sendtype ;
    }
    aSendTimeMessage->time = _t ;
    aSendTimeMessage->deltatime = _dt ;
    for ( target = 0 ; target < _GroupSize ; target++ ) {
       sts = Send( aSendTimeMessage , 1 , 0 , _MPIAccess->TimeType() , target ,
                   SendTimeRequestId ) ;
//       cout << "Send" << _MyRank << " --> target " << target << " _t " << _t
//            << " SendTimeRequestId " << SendTimeRequestId << endl ;
//            << " aSendTimeStruct/aSendTimeMessage " << aSendTimeStruct << "/"
//            << aSendTimeMessage << endl ;
       sts = Send( sendbuf , sendcount , sendoffset , sendtype , target , SendDataRequestId ) ;
//       cout << "Send" << _MyRank << " --> target " << target << " SendDataRequestId "
//            << SendDataRequestId << " aSendDataStruct/sendbuf " << aSendDataStruct << "/"
//            << sendbuf << endl ;
       if ( _Asynchronous ) {
         aSendTimeStruct->Counter += 1 ;
         (*_MapOfSendBuffers)[ SendTimeRequestId ] = aSendTimeStruct ;
         aSendDataStruct->Counter += 1 ;
         (*_MapOfSendBuffers)[ SendDataRequestId ] = aSendDataStruct ;
       }
       sendoffset += sendcount ;
    }
    if ( !_Asynchronous ) {
//      cout << "SynchronousAllToAllTime" << _MyRank << " free of SendTimeMessage & sendbuf "
//           << sendbuf << endl ;
      delete aSendTimeMessage ;
//      free( sendbuf ) ;
      if ( sendtype == MPI_INT ) {
        delete [] (int *) sendbuf ;
      }
      else {
        delete [] (double *) sendbuf ;
      }
    }
  }

//CheckTime + DoRecv + DoInterp
  if ( recvbuf && recvcount ) {
    for ( target = 0 ; target < _GroupSize ; target++ ) {
       int recvsize = recvcount*_MPIAccess->Extent( recvtype ) ;
//       bool OutOfTime ;
//       CheckTime( recvcount , recvsize , recvtype , target , OutOfTime) ;
       CheckTime( recvcount , recvtype , target , false ) ;
       //cout << "AllToAllTime" << _MyRank << " memcpy to recvbuf " << recvbuf
       //     << "+target" << target << "*" << recvsize << " bytes" << endl ;
//===========================================================================
//TODO : it is assumed actually that we have only 1 timestep before nad after
//===========================================================================
       if ( _TimeInterpolator && (*_TimeMessages)[target][0].time != -1 ) {
//       if ( _TimeInterpolator && !OutOfTime ) {
           if ( (*_OutOfTime)[target] ) {
             cout << " =====================================================" << endl
                  << "Recv" << _MyRank << " <-- target " << target << " t0 "
                  << (*_TimeMessages)[target][0].time << " < t1 "
                  << (*_TimeMessages)[target][1].time << " < t* " << _t << endl
                  << " =====================================================" << endl ;
           }
         if ( recvtype == MPI_INT ) {
           _TimeInterpolator->DoInterp( (*_TimeMessages)[target][0].time,
                                        (*_TimeMessages)[target][1].time, _t, recvcount ,
                                        _nStepBefore, _nStepAfter,
                                        (int **) &(*_DataMessages)[target][0],
                                        (int **) &(*_DataMessages)[target][1],
                                        &((int *)recvbuf)[target*recvcount] ) ;
         }
         else {
           _TimeInterpolator->DoInterp( (*_TimeMessages)[target][0].time,
                                        (*_TimeMessages)[target][1].time, _t, recvcount ,
                                        _nStepBefore, _nStepAfter,
                                        (double **) &(*_DataMessages)[target][0],
                                        (double **) &(*_DataMessages)[target][1],
                                        &((double *)recvbuf)[target*recvcount] ) ;
         }
       }
       else {
         char * buffdest = (char *) recvbuf ;
         char * buffsrc = (char *) (*_DataMessages)[target][1] ;
         memcpy( &buffdest[target*recvsize] , buffsrc , recvsize ) ;
       }
    }
  }

  return sts ;
}

int MPI_AccessDEC::AllToAllvTime( void* sendbuf, int* sendcounts, int* sdispls,
                                  MPI_Datatype sendtype ,
	                                void* recvbuf, int* recvcounts, int* rdispls,
                                  MPI_Datatype recvtype ) {
  int sts ;
  int target ;
  int SendTimeRequestId ;
  int SendDataRequestId ;

  if ( _TimeInterpolator == NULL ) {
    return MPI_ERR_OTHER ;
  }

//Free of SendBuffers 
  if ( _Asynchronous ) {
    CheckSent() ;
  }

/*
. DoSend :
  + We create a TimeMessage (look at that structure in MPI_Access).
  + If we are in asynchronous mode, we create two structures SendBuffStruct
    aSendTimeStruct and aSendDataStruct that we fill.
  + We fill the structure aSendTimeMessage with time/deltatime of
    the current process. "deltatime" must be nul if it is the last step of
    Time.
  + After that for each "target", we Send the TimeMessage and the part
    of sendbuf corresponding to that target.
  + If we are in asynchronous mode, we increment the counter and we add
     aSendTimeStruct and aSendDataStruct to _MapOfSendBuffers with the
    identifiers SendTimeRequestId and SendDataRequestId returned by
    MPI_Access->Send(...).
  + And if we are in synchronous mode we delete the SendMessages.
*/
//DoSend : Time + SendBuff
  SendBuffStruct * aSendTimeStruct = NULL ;
  SendBuffStruct * aSendDataStruct = NULL ;
//  cout << "AllToAllvTime" << _MyRank << " sendbuf " << sendbuf << " recvbuf " << recvbuf
//       << endl ;
  if ( sendbuf ) {
    TimeMessage * aSendTimeMessage = new TimeMessage ;
    if ( _Asynchronous ) {
      aSendTimeStruct = new SendBuffStruct ;
      aSendTimeStruct->SendBuffer = aSendTimeMessage ;
      aSendTimeStruct->Counter = 0 ;
      aSendTimeStruct->DataType = _MPIAccess->TimeType() ;
      aSendDataStruct = new SendBuffStruct ;
      aSendDataStruct->SendBuffer = sendbuf ;
      aSendDataStruct->Counter = 0 ;
      aSendDataStruct->DataType = sendtype ;
    }
    aSendTimeMessage->time = _t ;
    aSendTimeMessage->deltatime = _dt ;
    for ( target = 0 ; target < _GroupSize ; target++ ) {
       if ( sendcounts[target] ) {
         sts = Send( aSendTimeMessage , 1 , 0 , _MPIAccess->TimeType() , target ,
                     SendTimeRequestId ) ;
//         cout << "AllToAllvTime" << _MyRank << " TimeMessage target " << target
//              << " SendTimeRequestId " << SendTimeRequestId << " MPITag "
//              << _MPIAccess->SendMPITag(target) << endl ;
//         cout << "Send" << _MyRank << " --> target " << target << " _t " << _t
//              << " SendTimeRequestId " << SendTimeRequestId << endl ;
//              << " aSendTimeStruct/aSendTimeMessage " << aSendTimeStruct << "/"
//              << aSendTimeMessage << endl ;
         //int i ;
         //cout << "Send" << _MyRank << " --> target " << target
         //     << "sendcount" << sendcounts[target] << endl ;
         //for ( i = 0 ; i < sendcounts[target] ; i++ ) {
         //   cout << " " << ((double *) sendbuf)[sdispls[target]+i] ;
         //}
         //cout << endl ;
         sts = Send( sendbuf , sendcounts[target] , sdispls[target] , sendtype , target ,
                     SendDataRequestId ) ;
//         cout << "AllToAllvTime" << _MyRank << " DataMessage target " << target
//              << " SendDataRequestId " << SendDataRequestId << " MPITag "
//              << _MPIAccess->SendMPITag(target) << endl ;
         //cout << "Send" << _MyRank << " --> target " << target << " SendDataRequestId "
         //     << SendDataRequestId << " aSendDataStruct/sendbuf " << aSendDataStruct << "/"
         //     << sendbuf << " sendcount " << sendcounts[target] << " sdispls "
         //     << sdispls[target] << endl ;
         //int i ;
         //cout << "Send" << _MyRank << " --> target " << target << " count "
         //     << sendcounts[target] << endl ;
         //for ( i = 0 ; i < sendcounts[target] ; i++ ) {
         //   cout << " " << ((int *)sendbuf)[sdispls[target]+i] ;
         //}
         //cout << endl ;
         if ( _Asynchronous ) {
           aSendTimeStruct->Counter += 1 ;
           (*_MapOfSendBuffers)[ SendTimeRequestId ] = aSendTimeStruct ;
           //cout << "TimeSent" << _MyRank << " Request " << SendTimeRequestId
           //     << " _MapOfSendBuffers->SendBuffer "
           //     << (*_MapOfSendBuffers)[ SendTimeRequestId ]->SendBuffer
           //     << " Counter " << (*_MapOfSendBuffers)[ SendTimeRequestId ]->Counter
           //     << endl ;
           aSendDataStruct->Counter += 1 ;
           (*_MapOfSendBuffers)[ SendDataRequestId ] = aSendDataStruct ;
           //cout << "DataSent" << _MyRank << " Request " << SendDataRequestId
           //     << " _MapOfSendBuffers->SendBuffer "
           //     << (*_MapOfSendBuffers)[ SendDataRequestId ]->SendBuffer
           //     << " Counter " << (*_MapOfSendBuffers)[ SendDataRequestId ]->Counter
           //     << endl ;
         }
       }
    }
    if ( !_Asynchronous ) {
      //cout << "SynchronousAllToAllv" << _MyRank << " free of SendTimeMessage & sendbuf "
      //     << sendbuf << endl ;
      delete aSendTimeMessage ;
//      free( sendbuf ) ;
      if ( sendtype == MPI_INT ) {
        delete [] (int *) sendbuf ;
      }
      else {
        delete [] (double *) sendbuf ;
      }
    }
  }

/*
. CheckTime + DoRecv + DoInterp
  + For each target we call CheckTime
  + If there is a TimeInterpolator and if the TimeMessage of the target
    is not the first, we call the interpolator which return its
    results in the part of the recv buffer corresponding to the "target".
  + If not, there is a copy of received datas for that first step of time
    in the part of the recv buffer corresponding to the "target".
*/
//CheckTime + DoRecv + DoInterp
  if ( recvbuf ) {
    for ( target = 0 ; target < _GroupSize ; target++ ) {
       if ( recvcounts[target] ) {
         int recvsize = recvcounts[target]*_MPIAccess->Extent( recvtype ) ;
//         bool OutOfTime ;
//         CheckTime( recvcounts[target] , recvsize , recvtype , target , OutOfTime ) ;
         CheckTime( recvcounts[target] , recvtype , target , false ) ;
         //cout << "AllToAllvTime" << _MyRank << " memcpy to recvbuf " << recvbuf
         //     << "+target" << target << "*" << recvsize << " bytes" << endl ;
//===========================================================================
//TODO : it is assumed actually that we have only 1 timestep before nad after
//===========================================================================
         if ( _TimeInterpolator && (*_TimeMessages)[target][0].time != -1 ) {
//         if ( _TimeInterpolator && !OutOfTime ) {
             if ( (*_OutOfTime)[target] ) {
               cout << " =====================================================" << endl
                    << "Recv" << _MyRank << " <-- target " << target << " t0 "
                    << (*_TimeMessages)[target][0].time << " < t1 "
                    << (*_TimeMessages)[target][1].time << " < t* " << _t << endl
                    << " =====================================================" << endl ;
             }
           if ( recvtype == MPI_INT ) {
             _TimeInterpolator->DoInterp( (*_TimeMessages)[target][0].time,
                                          (*_TimeMessages)[target][1].time, _t,
                                          recvcounts[target] , _nStepBefore, _nStepAfter,
                                          (int **) &(*_DataMessages)[target][0],
                                          (int **) &(*_DataMessages)[target][1],
                                          &((int *)recvbuf)[rdispls[target]] ) ;
           }
           else {
             _TimeInterpolator->DoInterp( (*_TimeMessages)[target][0].time,
                                          (*_TimeMessages)[target][1].time, _t,
                                          recvcounts[target] , _nStepBefore, _nStepAfter,
                                          (double **) &(*_DataMessages)[target][0],
                                          (double **) &(*_DataMessages)[target][1],
                                          &((double *)recvbuf)[rdispls[target]] ) ;
           }
         }
         else {
           char * buffdest = (char *) recvbuf ;
           char * buffsrc = (char *) (*_DataMessages)[target][1] ;
           memcpy( &buffdest[rdispls[target]*_MPIAccess->Extent( recvtype )] , buffsrc ,
                   recvsize ) ;
         }
         //int i ;
         //cout << "Recv" << _MyRank << " --> target " << target << " recvbuf " << recvbuf
         //     << " recvbuf " << recvbuf << " &recvbuf[rdispls[target]] "
         //     << (&((int *) recvbuf)[rdispls[target]]) << endl ;
         //for ( i = 0 ; i < recvcounts[target] ; i++ ) {
         //   cout << " " << ((int *) (*_DataMessages)[target][1])[i] ;
         //}
         //cout << endl ;
         //cout << "Recv" << _MyRank << " <-- target " << target << " count "
         //     << recvcounts[target] << endl ;
         //for ( i = 0 ; i < recvcounts[target] ; i++ ) {
         //   cout << " " << (&((int *) recvbuf)[rdispls[target]])[i] ;
         //}
         //cout << endl ;
       }
    }
  }

  return sts ;
}

/*
. CheckTime(recvcount , recvtype , target , UntilEnd)
  + At the beginning, we read the first TimeMessage in
    &(*_TimeMessages)[target][1] and the first DataMessage
    in the allocated buffer (*_DataMessages)[target][1].
  + deltatime of TimesMessages must be nul if it is the last one.
  + While : _t(t*) is the current time of the processus.
    "while _t(t*) is greater than the time of the "target"
     (*_TimeMessages)[target][1].time and
     (*_TimeMessages)[target][1].deltatime is not nul",
    So at the end of the while we have :
     _t(t*) <= (*_TimeMessages)[target][1].time with
     _t(t*) > (*_TimeMessages)[target][0].time
    or we have the last TimeMessage of the "target".
  + If it is the finalization of the recv of TimeMessages and
    DataMessages (UntilEnd value is true), we execute the while
    until (*_TimeMessages)[target][1].deltatime is nul.
  + In the while :
    We copy the last TimeMessage in the previoud TimeMessage and
      we read a new TimeMessage
    We delete the previous DataMessage.
    We copy the last DataMessage pointer in the previous one.
    We allocate a new last DataMessage buffer
      (*_DataMessages)[target][1] and we read the corresponding
      datas in that buffe.
  + If the current time of the current process is greater than the
    last time (*_TimeMessages)[target][1].time du target, we give
    a true value to (*_OutOfTime)[target].
    (*_TimeMessages)[target][1].deltatime is nul.
*/
int MPI_AccessDEC::CheckTime( int recvcount , MPI_Datatype recvtype , int target ,
                              bool UntilEnd ) {
//                              int target , bool &OutOfTime ) {
  //cout << "CheckTime" << _MyRank << " time " << _t << " deltatime " << _dt << endl ;

  int sts = MPI_SUCCESS ;
  int RecvTimeRequestId ;
  int RecvDataRequestId ;
//Pour l'instant on cherche _TimeMessages[target][0] < _t <= _TimeMessages[target][1]
//===========================================================================
//TODO : it is assumed actually that we have only 1 timestep before and after
//       instead of _nStepBefore and _nStepAfter ...
//===========================================================================
  (*_DataMessagesRecvCount)[target] = recvcount ;
  (*_DataMessagesType)[target] = recvtype ;
//  (*_OutOfTime)[target] = false ;
//Actually we need 1 timestep before and after
//  if ( _steptime <= 1 ) {
  if ( (*_TimeMessages)[target][1].time == -1 ) {
    (*_TimeMessages)[target][0] = (*_TimeMessages)[target][1] ;
    sts = Recv( &(*_TimeMessages)[target][1] , 1 , _MPIAccess->TimeType() ,
                target , RecvTimeRequestId ) ;
    (*_DataMessages)[target][0] = (*_DataMessages)[target][1] ;
//    (*_DataMessages)[target][1] = malloc( recvsize ) ;
    if ( recvtype == MPI_INT ) {
      (*_DataMessages)[target][1] = new int[recvcount] ;
    }
    else {
      (*_DataMessages)[target][1] = new double[recvcount] ;
    }
    sts = Recv( (*_DataMessages)[target][1] , recvcount , recvtype , target ,
                RecvDataRequestId ) ;
  }
  else {
    while ( ( _t > (*_TimeMessages)[target][1].time || UntilEnd ) &&
            (*_TimeMessages)[target][1].deltatime != 0 ) {
         //cout << "CheckTime" << _MyRank << " TimeMessage target " << target << " _t "
         //     << _t << " > " << (*_TimeMessages)[target][1].time << " et "
         //     << (*_TimeMessages)[target][1].deltatime
         //     << " != 0 ==> Recv" << endl ;
         (*_TimeMessages)[target][0] = (*_TimeMessages)[target][1] ;
         sts = Recv( &(*_TimeMessages)[target][1] , 1 , _MPIAccess->TimeType() ,
                     target , RecvTimeRequestId ) ;
         if ( UntilEnd ) {
           cout << "CheckTime" << _MyRank << " TimeMessage target " << target
                << " RecvTimeRequestId " << RecvTimeRequestId << " MPITag "
                << _MPIAccess->RecvMPITag(target) << endl ;
         }
         if ( recvtype == MPI_INT ) {
           delete [] (int *) (*_DataMessages)[target][0] ;
         }
         else {
           delete [] (double *) (*_DataMessages)[target][0] ;
         }
         (*_DataMessages)[target][0] = (*_DataMessages)[target][1] ;
         if ( recvtype == MPI_INT ) {
           (*_DataMessages)[target][1] = new int[recvcount] ;
         }
         else {
           (*_DataMessages)[target][1] = new double[recvcount] ;
         }
         sts = Recv( (*_DataMessages)[target][1] , recvcount , recvtype , target ,
                     RecvDataRequestId ) ;
         if ( UntilEnd ) {
           cout << "CheckTime" << _MyRank << " DataMessage target " << target
                << " RecvDataRequestId " << RecvDataRequestId << " MPITag "
                << _MPIAccess->RecvMPITag(target) << endl ;
         }
    }

    if ( _t > (*_TimeMessages)[target][0].time &&
         _t <= (*_TimeMessages)[target][1].time ) {
//Ok
    }
    else {
      (*_OutOfTime)[target] = true ;
    }
  }

  //cout << "CheckTime" << _MyRank << " TimeMessage target " << target << " Time "
  //     << (*_TimeMessages)[target][1].time << " DataMessage" ;
  //int i ;
  //void * buff = (*_DataMessages)[target][1] ;
  //for ( i = 0 ; i < recvcount ; i++ ) {
  //   cout << " " << ((int *) buff)[i] ;
  //}
  //cout << endl ;
  return sts ;
}

/*
. CheckSent() :
  + call  SendRequestIds of MPI_Access in order to get all
    RequestIds of SendMessages of all "targets".
  + For each RequestId, CheckSent call "Test" of MPI_Access in order
    to know if the buffer is "free" (flag = true). If it is the
    FinalCheckSent (WithWait = true), we call Wait instead of Test.
  + If the buffer is "free", the counter of the structure SendBuffStruct
    (from _MapOfSendBuffers) is decremented.
  + If that counter is nul we delete the TimeMessage or the
    SendBuffer according to the DataType.
  + And we delete the structure SendBuffStruct before the suppression
    (erase) of that item of _MapOfSendBuffers
*/
int MPI_AccessDEC::CheckSent(bool WithWait) {
  int sts = MPI_SUCCESS ;
  int flag = WithWait ;
  int size = _MPIAccess->SendRequestIdsSize() ;
  int * ArrayOfSendRequests = new int[ size ] ;
  int nSendRequest = _MPIAccess->SendRequestIds( size , ArrayOfSendRequests ) ;
  bool SendTrace = false ;
//  if ( nSendRequest > 2*_GroupSize || (_t > 1800 && _MyRank== 0) ) {
//    SendTrace = true ;
//    _MPIAccess->Trace() ;
//    if ( nSendRequest > 3*_GroupSize ) {
//      WithWait = true ;
//    }
//  }
//  cout << "CheckSent" << _MyRank << " nSendRequest " << nSendRequest << " SendTrace "
//       << SendTrace << " WithWait " << WithWait << " :" << endl ;
  int i ;
  for ( i = 0 ; i < nSendRequest ; i++ ) {
     if ( WithWait ) {
       cout << "CheckSent" << _MyRank << " " << i << "./" << nSendRequest
            << " SendRequestId " << ArrayOfSendRequests[i] << " MPITarget "
            << _MPIAccess->MPITarget(ArrayOfSendRequests[i]) << " MPITag "
            << _MPIAccess->MPITag(ArrayOfSendRequests[i]) << " Wait :" << endl ;
       sts = _MPIAccess->Wait( ArrayOfSendRequests[i] ) ;
     }
     else {
       sts = _MPIAccess->Test( ArrayOfSendRequests[i] , flag ) ;
     }
     if ( flag ) {
       _MPIAccess->DeleteRequest( ArrayOfSendRequests[i] ) ;
       if ( SendTrace ) {
         cout << "CheckSent" << _MyRank << " " << i << "./" << nSendRequest
              << " SendRequestId " << ArrayOfSendRequests[i]
              << " flag " << flag
              << " Counter " << (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->Counter
              << " DataType " << (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->DataType
              << endl ;
       }
       (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->Counter -= 1 ;
       if ( SendTrace ) {
         if ( (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->DataType == 
              _MPIAccess->TimeType() ) {
           cout << "CheckTimeSent" << _MyRank << " Request " ;
         }
         else {
           cout << "CheckDataSent" << _MyRank << " Request " ;
         }
         cout << ArrayOfSendRequests[i]
              << " _MapOfSendBuffers->SendBuffer "
              << (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->SendBuffer
              << " Counter " << (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->Counter
              << endl ;
       }
       if ( (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->Counter  == 0 ) {
         if ( SendTrace ) {
           cout << "CheckSent" << _MyRank << " SendRequestId " << ArrayOfSendRequests[i]
                << " Counter " << (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->Counter
                << " flag " << flag << " SendBuffer "
                << (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->SendBuffer
                << " deleted. Erase in _MapOfSendBuffers :" << endl ;
         }
         if ( (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->DataType ==
              _MPIAccess->TimeType() ) {
           delete (TimeMessage * ) (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->SendBuffer ;
         }
         else {
           if ( (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->DataType == MPI_INT ) {
             delete [] (int *) (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->SendBuffer ;
           }
           else {
             delete [] (double *) (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->SendBuffer ;
           }
         }
         delete (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ] ;
       }
       if ( SendTrace ) {
         cout << "CheckSent" << _MyRank << " Erase in _MapOfSendBuffers SendRequestId "
              << ArrayOfSendRequests[i] << endl ;
       }
       (*_MapOfSendBuffers).erase( ArrayOfSendRequests[i] ) ;
     }
     else if ( SendTrace ) {
       cout << "CheckSent" << _MyRank << " " << i << "./" << nSendRequest
            << " SendRequestId " << ArrayOfSendRequests[i]
            << " flag " << flag
            << " Counter " << (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->Counter
            << " DataType " << (*_MapOfSendBuffers)[ ArrayOfSendRequests[i] ]->DataType
            << endl ;
     }
  }
  if ( SendTrace ) {
    _MPIAccess->Check() ;
  }
  delete [] ArrayOfSendRequests ;
  return sts ;
}

int MPI_AccessDEC::CheckFinalRecv() {
  int sts = MPI_SUCCESS ;
  if ( _TimeInterpolator ) {
    int target ;
    for ( target = 0 ; target < _GroupSize ; target++ ) {
       if ( (*_DataMessages)[target][0] != NULL ) {
         sts = CheckTime( (*_DataMessagesRecvCount)[target] , (*_DataMessagesType)[target] ,
                          target , true ) ;
         if ( (*_DataMessagesType)[target] == MPI_INT ) {
           delete [] (int *) (*_DataMessages)[target][0] ;
         }
         else {
           delete [] (double *) (*_DataMessages)[target][0] ;
         }
         (*_DataMessages)[target][0] = NULL ;
         if ( (*_DataMessages)[target][1] != NULL ) {
           if ( (*_DataMessagesType)[target] == MPI_INT ) {
             delete [] (int *) (*_DataMessages)[target][1] ;
           }
           else {
             delete [] (double *) (*_DataMessages)[target][1] ;
           }
           (*_DataMessages)[target][1] = NULL ;
         }
       }
     }
  }
//  return _MPIAccess->CancelAll() ;
  return sts ;
}

ostream & operator<< (ostream & f ,const TimeInterpolationMethod & interpolationmethod ) {
  switch (interpolationmethod) {
  case WithoutTimeInterp :
    f << " WithoutTimeInterpolation ";
    break;
  case LinearTimeInterp :
    f << " LinearTimeInterpolation ";
    break;
  default :
    f << " UnknownTimeInterpolation ";
    break;
  }

  return f;
}


}

