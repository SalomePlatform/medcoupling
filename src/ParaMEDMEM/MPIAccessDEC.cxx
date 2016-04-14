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

#include "MPIAccessDEC.hxx"

#include <cstring>

using namespace std;

namespace MEDCoupling
{    

  /*!
    This constructor creates an MPIAccessDEC which has \a source_group as a working side 
    and  \a target_group as an idle side. 
    The constructor must be called synchronously on all processors of both processor groups.

    \param source_group working side ProcessorGroup
    \param target_group lazy side ProcessorGroup
    \param Asynchronous Communication mode (default asynchronous)
    \param nStepBefore Number of Time step needed for the interpolation before current time
    \param nStepAfter Number of Time step needed for the interpolation after current time

  */

  MPIAccessDEC::MPIAccessDEC( const ProcessorGroup& source_group,
                              const ProcessorGroup& target_group,
                              bool Asynchronous )
  {

    ProcessorGroup * union_group = source_group.fuse(target_group) ;  
    int i ;
    std::set<int> procs;
    for ( i = 0 ; i < union_group->size() ; i++ )
      {
        procs.insert(i) ;
      }
    MPIProcessorGroup *mpilg = static_cast<MPIProcessorGroup *>(const_cast<ProcessorGroup *>(&source_group));
    _MPI_union_group = new MEDCoupling::MPIProcessorGroup( union_group->getCommInterface(),procs,mpilg->getWorldComm());
    delete union_group ;
    _my_rank = _MPI_union_group->myRank() ;
    _group_size = _MPI_union_group->size() ;
    _MPI_access = new MPIAccess( _MPI_union_group ) ;
    _asynchronous = Asynchronous ;
    _time_messages = new vector< vector< TimeMessage > > ;
    _time_messages->resize( _group_size ) ;
    _out_of_time = new vector< bool > ;
    _out_of_time->resize( _group_size ) ;
    _data_messages_recv_count = new vector< int > ;
    _data_messages_recv_count->resize( _group_size ) ;
    for ( i = 0 ; i < _group_size ; i++ )
      {
        (*_out_of_time)[i] = false ;
        (*_data_messages_recv_count)[i] = 0 ;
      }
    _data_messages_type = new vector< MPI_Datatype > ;
    _data_messages_type->resize( _group_size ) ;
    _data_messages = new vector< vector< void * > > ;
    _data_messages->resize( _group_size ) ;
    _time_interpolator = NULL ;
    _map_of_send_buffers = new map< int , SendBuffStruct * > ;
  }

  MPIAccessDEC::~MPIAccessDEC()
  {
    checkFinalSent() ;
    checkFinalRecv() ;
    delete _MPI_union_group ;
    delete _MPI_access ;
    if ( _time_interpolator )
      delete _time_interpolator ;
    if ( _time_messages )
      delete _time_messages ;
    if ( _out_of_time )
      delete _out_of_time ;
    if ( _data_messages_recv_count )
      delete _data_messages_recv_count ;
    if ( _data_messages_type )
      delete _data_messages_type ;
    if ( _data_messages )
      delete _data_messages ;
    if ( _map_of_send_buffers )
      delete _map_of_send_buffers ;
  } 

  void MPIAccessDEC::setTimeInterpolator( TimeInterpolationMethod aTimeInterp ,
                                          double InterpPrecision, int nStepBefore,
                                          int nStepAfter )
  {
    if ( _time_interpolator )
      delete _time_interpolator ;
    switch ( aTimeInterp )
      {
      case WithoutTimeInterp :
        _time_interpolator = NULL ;
        _n_step_before = 0 ;
        _n_step_after = 0 ;
        break ;
      case LinearTimeInterp :
        _time_interpolator = new LinearTimeInterpolator( InterpPrecision , nStepBefore ,
                                                         nStepAfter ) ;
        _n_step_before = nStepBefore ;
        _n_step_after = nStepAfter ;
        int i ;
        for ( i = 0 ; i < _group_size ; i++ )
          {
            (*_time_messages)[ i ].resize( _n_step_before + _n_step_after ) ;
            (*_data_messages)[ i ].resize( _n_step_before + _n_step_after ) ;
            int j ;
            for ( j = 0 ; j < _n_step_before + _n_step_after ; j++ )
              {
                (*_time_messages)[ i ][ j ].time = -1 ;
                (*_time_messages)[ i ][ j ].deltatime = -1 ;
                (*_data_messages)[ i ][ j ] = NULL ;
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
  int MPIAccessDEC::send( void* sendbuf, int sendcount , int offset ,
                          MPI_Datatype sendtype , int target , int &SendRequestId )
  {
    int sts ;
    if ( _asynchronous )
      {
        if ( sendtype == MPI_INT )
          {
            sts = _MPI_access->ISend( &((int *) sendbuf)[offset] , sendcount , sendtype ,
                                      target , SendRequestId ) ;
          }
        else
          {
            sts = _MPI_access->ISend( &((double *) sendbuf)[offset] , sendcount , sendtype ,
                                      target , SendRequestId ) ;
          }
      }
    else
      {
        if ( sendtype == MPI_INT )
          {
            sts = _MPI_access->send( &((int *) sendbuf)[offset] , sendcount , sendtype ,
                                     target , SendRequestId ) ;
          }
        else
          {
            sts = _MPI_access->send( &((double *) sendbuf)[offset] , sendcount , sendtype ,
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
  int MPIAccessDEC::recv( void* recvbuf, int recvcount , int offset ,
                          MPI_Datatype recvtype , int target , int &RecvRequestId )
  {
    int sts ;
    if ( _asynchronous )
      {
        if ( recvtype == MPI_INT )
          {
            sts = _MPI_access->IRecv( &((int *) recvbuf)[offset] , recvcount , recvtype ,
                                      target , RecvRequestId ) ;
          }
        else
          {
            sts = _MPI_access->IRecv( &((double *) recvbuf)[offset] , recvcount , recvtype ,
                                      target , RecvRequestId ) ;
          }
      }
    else
      {
        if ( recvtype == MPI_INT )
          {
            sts = _MPI_access->recv( &((int *) recvbuf)[offset] , recvcount , recvtype ,
                                     target , RecvRequestId ) ;
          }
        else
          {
            sts = _MPI_access->recv( &((double *) recvbuf)[offset] , recvcount , recvtype ,
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
  int MPIAccessDEC::sendRecv( void* sendbuf, int sendcount , int sendoffset ,
                              MPI_Datatype sendtype ,
                              void* recvbuf, int recvcount , int recvoffset ,
                              MPI_Datatype recvtype , int target ,
                              int &SendRequestId , int &RecvRequestId )
  {
    int sts ;
    if ( _asynchronous )
      {
        if ( sendtype == MPI_INT )
          {
            if ( recvtype == MPI_INT )
              {
                sts = _MPI_access->ISendRecv( &((int *) sendbuf)[sendoffset] , sendcount ,
                                              sendtype , target , SendRequestId ,
                                              &((int *) recvbuf)[recvoffset] , recvcount ,
                                              recvtype , target , RecvRequestId ) ;
              }
            else
              {
                sts = _MPI_access->ISendRecv( &((int *) sendbuf)[sendoffset] , sendcount ,
                                              sendtype , target , SendRequestId ,
                                              &((double *) recvbuf)[recvoffset] ,
                                              recvcount , recvtype , target , RecvRequestId ) ;
              }
          }
        else
          {
            if ( recvtype == MPI_INT )
              {
                sts = _MPI_access->ISendRecv( &((double *) sendbuf)[sendoffset] , sendcount ,
                                              sendtype , target , SendRequestId ,
                                              &((int *) recvbuf)[recvoffset] ,
                                              recvcount , recvtype , target , RecvRequestId ) ;
              }
            else
              {
                sts = _MPI_access->ISendRecv( &((double *) sendbuf)[sendoffset] , sendcount ,
                                              sendtype , target , SendRequestId ,
                                              &((double *) recvbuf)[recvoffset] ,
                                              recvcount , recvtype , target , RecvRequestId ) ;
              }
          }
      }
    else
      {
        if ( sendtype == MPI_INT )
          {
            if ( recvtype == MPI_INT )
              {
                sts = _MPI_access->sendRecv( &((int *) sendbuf)[sendoffset] , sendcount ,
                                             sendtype , target , SendRequestId ,
                                             &((int *) recvbuf)[recvoffset] , recvcount ,
                                             recvtype , target , RecvRequestId ) ;
              }
            else
              {
                sts = _MPI_access->sendRecv( &((int *) sendbuf)[sendoffset] , sendcount ,
                                             sendtype , target , SendRequestId ,
                                             &((double *) recvbuf)[recvoffset] ,
                                             recvcount , recvtype , target , RecvRequestId ) ;
              }
          }
        else
          {
            if ( recvtype == MPI_INT )
              {
                sts = _MPI_access->sendRecv( &((double *) sendbuf)[sendoffset] , sendcount ,
                                             sendtype , target , SendRequestId ,
                                             &((int *) recvbuf)[recvoffset] ,
                                             recvcount , recvtype , target , RecvRequestId ) ;
              }
            else
              {
                sts = _MPI_access->sendRecv( &((double *) sendbuf)[sendoffset] , sendcount ,
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
  int MPIAccessDEC::allToAll( void* sendbuf, int sendcount, MPI_Datatype sendtype ,
                              void* recvbuf, int recvcount, MPI_Datatype recvtype )
  {
    if ( _time_interpolator )
      {
        return allToAllTime( sendbuf, sendcount, sendtype , recvbuf, recvcount, recvtype ) ;
      }
    int sts ;
    int target ;
    int sendoffset = 0 ;
    int recvoffset = 0 ;
    int SendRequestId ;
    int RecvRequestId ;

    //Free of SendBuffers 
    if ( _asynchronous )
      checkSent() ;

    //DoSend + DoRecv : SendRecv
    SendBuffStruct * aSendDataStruct = NULL ;
    if ( _asynchronous && sendbuf )
      {
        aSendDataStruct = new SendBuffStruct ;
        aSendDataStruct->SendBuffer = sendbuf ;
        aSendDataStruct->Counter = 0 ;
        aSendDataStruct->DataType = sendtype ;
      }
    for ( target = 0 ; target < _group_size ; target++ )
      {
        sts = sendRecv( sendbuf , sendcount , sendoffset , sendtype ,
                        recvbuf , recvcount , recvoffset , recvtype ,
                        target , SendRequestId , RecvRequestId ) ;
        if ( _asynchronous && sendbuf && sendcount )
          {
            aSendDataStruct->Counter += 1 ;
            (*_map_of_send_buffers)[ SendRequestId ] = aSendDataStruct ;
          }
        sendoffset += sendcount ;
        recvoffset += recvcount ;
      }
    if ( !_asynchronous && sendbuf )
      {
        if ( sendtype == MPI_INT )
          {
            delete [] (int *) sendbuf ;
          }
        else
          {
            delete [] (double *) sendbuf ;
          }
      }
    return sts ;
  }

  /*!
    Send sendcounts[target] datas from sendbuf[sdispls[target]] with type sendtype to all targets of IntraCommunicator
    Receive recvcounts[target] datas to recvbuf[rdispls[target]] with type recvtype from all targets of IntraCommunicator

  */
  int MPIAccessDEC::allToAllv( void* sendbuf, int* sendcounts, int* sdispls,
                               MPI_Datatype sendtype ,
                               void* recvbuf, int* recvcounts, int* rdispls,
                               MPI_Datatype recvtype )
  {
    if ( _time_interpolator )
      {
        return allToAllvTime( sendbuf, sendcounts, sdispls, sendtype ,
                              recvbuf, recvcounts, rdispls, recvtype ) ;
      }
    int sts ;
    int target ;
    int SendRequestId ;
    int RecvRequestId ;

    //Free of SendBuffers 
    if ( _asynchronous )
      {
        checkSent() ;
      }

    //DoSend + DoRecv : SendRecv
    SendBuffStruct * aSendDataStruct = NULL ;
    if ( _asynchronous && sendbuf )
      {
        aSendDataStruct = new SendBuffStruct ;
        aSendDataStruct->SendBuffer = sendbuf ;
        aSendDataStruct->Counter = 0 ;
        aSendDataStruct->DataType = sendtype ;
      }
    for ( target = 0 ; target < _group_size ; target++ )
      {
        if ( sendcounts[target] || recvcounts[target] )
          {
            sts = sendRecv( sendbuf , sendcounts[target] , sdispls[target] , sendtype ,
                            recvbuf , recvcounts[target] , rdispls[target] , recvtype ,
                            target , SendRequestId , RecvRequestId ) ;
            if ( _asynchronous && sendbuf && sendcounts[target])
              {
                aSendDataStruct->Counter += 1 ;
                (*_map_of_send_buffers)[ SendRequestId ] = aSendDataStruct ;
              }
          }
      }
    if ( !_asynchronous && sendbuf )
      {
        if ( sendtype == MPI_INT )
          {
            delete [] (int *) sendbuf ;
          }
        else
          {
            delete [] (double *) sendbuf ;
          }
      }
    return sts ;
  }

  /*
    MPIAccessDEC and the management of SendBuffers :
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
    MPIAccessDEC and the management of RecvBuffers :
    =================================================

    If there is no interpolation, no special action is done.

    With interpolation for each target :
    ------------------------------------
    . We have _time_messages[target] which is a vector of TimesMessages.
    We have 2 TimesMessages in our case with a linear interpolation.
    They contain the previous time(t0)/deltatime and the last
    time(t1)/deltatime.

    . We have _data_messages[target] which is a vector of DatasMessages.
    We have 2 DatasMessages in our case with a linear interpolation.
    They contain the previous datas at time(t0)/deltatime and at last
    time(t1)/deltatime.

    . At time _t(t*) of current processus we do the interpolation of
    the values of the 2 DatasMessages which are returned in the part of
    recvbuf corresponding to the target with t0 < t* <= t1.

    . Because of the difference of "deltatimes" between processes, we
    may have t0 < t1 < t* and there is an extrapolation.

    . The vectors _out_of_time, _DataMessagesRecvCount and _DataMessagesType
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
  int MPIAccessDEC::allToAllTime( void* sendbuf, int sendcount , MPI_Datatype sendtype ,
                                  void* recvbuf, int recvcount , MPI_Datatype recvtype )
  {
    int sts ;
    int target ;
    int sendoffset = 0 ;
    int SendTimeRequestId ;
    int SendDataRequestId ;

    if ( _time_interpolator == NULL )
      {
        return MPI_ERR_OTHER ;
      }

    //Free of SendBuffers 
    if ( _asynchronous )
      {
        checkSent() ;
      }

    //DoSend : Time + SendBuff
    SendBuffStruct * aSendTimeStruct = NULL ;
    SendBuffStruct * aSendDataStruct = NULL ;
    if ( sendbuf && sendcount )
      {
        TimeMessage * aSendTimeMessage = new TimeMessage ;
        if ( _asynchronous )
          {
            aSendTimeStruct = new SendBuffStruct ;
            aSendTimeStruct->SendBuffer = aSendTimeMessage ;
            aSendTimeStruct->Counter = 0 ;
            aSendTimeStruct->DataType = _MPI_access->timeType() ;
            aSendDataStruct = new SendBuffStruct ;
            aSendDataStruct->SendBuffer = sendbuf ;
            aSendDataStruct->Counter = 0 ;
            aSendDataStruct->DataType = sendtype ;
          }
        aSendTimeMessage->time = _t ;
        aSendTimeMessage->deltatime = _dt ;
        for ( target = 0 ; target < _group_size ; target++ )
          {
            sts = send( aSendTimeMessage , 1 , 0 , _MPI_access->timeType() , target ,
                        SendTimeRequestId ) ;
            sts = send( sendbuf , sendcount , sendoffset , sendtype , target , SendDataRequestId ) ;
            if ( _asynchronous )
              {
                aSendTimeStruct->Counter += 1 ;
                (*_map_of_send_buffers)[ SendTimeRequestId ] = aSendTimeStruct ;
                aSendDataStruct->Counter += 1 ;
                (*_map_of_send_buffers)[ SendDataRequestId ] = aSendDataStruct ;
              }
            sendoffset += sendcount ;
          }
        if ( !_asynchronous )
          {
            delete aSendTimeMessage ;
            if ( sendtype == MPI_INT )
              {
                delete [] (int *) sendbuf ;
              }
            else
              {
                delete [] (double *) sendbuf ;
              }
          }
      }

    //CheckTime + DoRecv + DoInterp
    if ( recvbuf && recvcount )
      {
        for ( target = 0 ; target < _group_size ; target++ )
          {
            int recvsize = recvcount*_MPI_access->extent( recvtype ) ;
            checkTime( recvcount , recvtype , target , false ) ;
            //===========================================================================
            //TODO : it is assumed actually that we have only 1 timestep before and after
            //===========================================================================
            if ( _time_interpolator && (*_time_messages)[target][0].time != -1 )
              {
                if ( (*_out_of_time)[target] )
                  {
                    cout << " =====================================================" << endl
                         << "Recv" << _my_rank << " <-- target " << target << " t0 "
                         << (*_time_messages)[target][0].time << " < t1 "
                         << (*_time_messages)[target][1].time << " < t* " << _t << endl
                         << " =====================================================" << endl ;
                  }
                if ( recvtype == MPI_INT )
                  {
                    _time_interpolator->doInterp( (*_time_messages)[target][0].time,
                                                  (*_time_messages)[target][1].time, _t, recvcount ,
                                                  _n_step_before, _n_step_after,
                                                  (int **) &(*_data_messages)[target][0],
                                                  (int **) &(*_data_messages)[target][1],
                                                  &((int *)recvbuf)[target*recvcount] ) ;
                  }
                else
                  {
                    _time_interpolator->doInterp( (*_time_messages)[target][0].time,
                                                  (*_time_messages)[target][1].time, _t, recvcount ,
                                                  _n_step_before, _n_step_after,
                                                  (double **) &(*_data_messages)[target][0],
                                                  (double **) &(*_data_messages)[target][1],
                                                  &((double *)recvbuf)[target*recvcount] ) ;
                  }
              }
            else
              {
                char * buffdest = (char *) recvbuf ;
                char * buffsrc = (char *) (*_data_messages)[target][1] ;
                memcpy( &buffdest[target*recvsize] , buffsrc , recvsize ) ;
              }
          }
      }

    return sts ;
  }

  int MPIAccessDEC::allToAllvTime( void* sendbuf, int* sendcounts, int* sdispls,
                                   MPI_Datatype sendtype ,
                                   void* recvbuf, int* recvcounts, int* rdispls,
                                   MPI_Datatype recvtype )
  {
    int sts ;
    int target ;
    int SendTimeRequestId ;
    int SendDataRequestId ;

    if ( _time_interpolator == NULL )
      {
        return MPI_ERR_OTHER ;
      }

    //Free of SendBuffers 
    if ( _asynchronous )
      {
        checkSent() ;
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
    if ( sendbuf )
      {
        TimeMessage * aSendTimeMessage = new TimeMessage ;
        if ( _asynchronous )
          {
            aSendTimeStruct = new SendBuffStruct ;
            aSendTimeStruct->SendBuffer = aSendTimeMessage ;
            aSendTimeStruct->Counter = 0 ;
            aSendTimeStruct->DataType = _MPI_access->timeType() ;
            aSendDataStruct = new SendBuffStruct ;
            aSendDataStruct->SendBuffer = sendbuf ;
            aSendDataStruct->Counter = 0 ;
            aSendDataStruct->DataType = sendtype ;
          }
        aSendTimeMessage->time = _t ;
        aSendTimeMessage->deltatime = _dt ;
        for ( target = 0 ; target < _group_size ; target++ )
          {
            if ( sendcounts[target] )
              {
                sts = send( aSendTimeMessage , 1 , 0 , _MPI_access->timeType() , target ,
                            SendTimeRequestId ) ;
                sts = send( sendbuf , sendcounts[target] , sdispls[target] , sendtype , target ,
                            SendDataRequestId ) ;
                if ( _asynchronous )
                  {
                    aSendTimeStruct->Counter += 1 ;
                    (*_map_of_send_buffers)[ SendTimeRequestId ] = aSendTimeStruct ;
                    aSendDataStruct->Counter += 1 ;
                    (*_map_of_send_buffers)[ SendDataRequestId ] = aSendDataStruct ;
                  }
              }
          }
        if ( !_asynchronous )
          {
            delete aSendTimeMessage ;
            if ( sendtype == MPI_INT )
              {
                delete [] (int *) sendbuf ;
              }
            else
              {
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
    if ( recvbuf )
      {
        for ( target = 0 ; target < _group_size ; target++ )
          {
            if ( recvcounts[target] )
              {
                int recvsize = recvcounts[target]*_MPI_access->extent( recvtype ) ;
                checkTime( recvcounts[target] , recvtype , target , false ) ;
                //===========================================================================
                //TODO : it is assumed actually that we have only 1 timestep before nad after
                //===========================================================================
                if ( _time_interpolator && (*_time_messages)[target][0].time != -1 )
                  {
                    if ( (*_out_of_time)[target] )
                      {
                        cout << " =====================================================" << endl
                             << "Recv" << _my_rank << " <-- target " << target << " t0 "
                             << (*_time_messages)[target][0].time << " < t1 "
                             << (*_time_messages)[target][1].time << " < t* " << _t << endl
                             << " =====================================================" << endl ;
                      }
                    if ( recvtype == MPI_INT )
                      {
                        _time_interpolator->doInterp( (*_time_messages)[target][0].time,
                                                      (*_time_messages)[target][1].time, _t,
                                                      recvcounts[target] , _n_step_before, _n_step_after,
                                                      (int **) &(*_data_messages)[target][0],
                                                      (int **) &(*_data_messages)[target][1],
                                                      &((int *)recvbuf)[rdispls[target]] ) ;
                      }
                    else
                      {
                        _time_interpolator->doInterp( (*_time_messages)[target][0].time,
                                                      (*_time_messages)[target][1].time, _t,
                                                      recvcounts[target] , _n_step_before, _n_step_after,
                                                      (double **) &(*_data_messages)[target][0],
                                                      (double **) &(*_data_messages)[target][1],
                                                      &((double *)recvbuf)[rdispls[target]] ) ;
                      }
                  }
                else
                  {
                    char * buffdest = (char *) recvbuf ;
                    char * buffsrc = (char *) (*_data_messages)[target][1] ;
                    memcpy( &buffdest[rdispls[target]*_MPI_access->extent( recvtype )] , buffsrc ,
                            recvsize ) ;
                  }
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
  int MPIAccessDEC::checkTime( int recvcount , MPI_Datatype recvtype , int target ,
                               bool UntilEnd )
  {
    int sts = MPI_SUCCESS ;
    int RecvTimeRequestId ;
    int RecvDataRequestId ;
    //Pour l'instant on cherche _time_messages[target][0] < _t <= _time_messages[target][1]
    //===========================================================================
    //TODO : it is assumed actually that we have only 1 timestep before and after
    //       instead of _n_step_before and _n_step_after ...
    //===========================================================================
    (*_data_messages_recv_count)[target] = recvcount ;
    (*_data_messages_type)[target] = recvtype ;
    if ( (*_time_messages)[target][1].time == -1 )
      {
        (*_time_messages)[target][0] = (*_time_messages)[target][1] ;
        sts = recv( &(*_time_messages)[target][1] , 1 , _MPI_access->timeType() ,
                    target , RecvTimeRequestId ) ;
        (*_data_messages)[target][0] = (*_data_messages)[target][1] ;
        if ( recvtype == MPI_INT )
          {
            (*_data_messages)[target][1] = new int[recvcount] ;
          }
        else
          {
            (*_data_messages)[target][1] = new double[recvcount] ;
          }
        sts = recv( (*_data_messages)[target][1] , recvcount , recvtype , target ,
                    RecvDataRequestId ) ;
      }
    else
      {
        while ( ( _t > (*_time_messages)[target][1].time || UntilEnd ) &&
                (*_time_messages)[target][1].deltatime != 0 )
          {
            (*_time_messages)[target][0] = (*_time_messages)[target][1] ;
            sts = recv( &(*_time_messages)[target][1] , 1 , _MPI_access->timeType() ,
                        target , RecvTimeRequestId ) ;
            if ( UntilEnd )
              {
                cout << "CheckTime" << _my_rank << " TimeMessage target " << target
                     << " RecvTimeRequestId " << RecvTimeRequestId << " MPITag "
                     << _MPI_access->recvMPITag(target) << endl ;
              }
            if ( recvtype == MPI_INT )
              {
                delete [] (int *) (*_data_messages)[target][0] ;
              }
            else
              {
                delete [] (double *) (*_data_messages)[target][0] ;
              }
            (*_data_messages)[target][0] = (*_data_messages)[target][1] ;
            if ( recvtype == MPI_INT )
              {
                (*_data_messages)[target][1] = new int[recvcount] ;
              }
            else
              {
                (*_data_messages)[target][1] = new double[recvcount] ;
              }
            sts = recv( (*_data_messages)[target][1] , recvcount , recvtype , target ,
                        RecvDataRequestId ) ;
            if ( UntilEnd )
              {
                cout << "CheckTime" << _my_rank << " DataMessage target " << target
                     << " RecvDataRequestId " << RecvDataRequestId << " MPITag "
                     << _MPI_access->recvMPITag(target) << endl ;
              }
          }

        if ( _t > (*_time_messages)[target][0].time &&
             _t <= (*_time_messages)[target][1].time )
          {
          }
        else
          {
            (*_out_of_time)[target] = true ;
          }
      }
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
  int MPIAccessDEC::checkSent(bool WithWait)
  {
    int sts = MPI_SUCCESS ;
    int flag = WithWait ;
    int size = _MPI_access->sendRequestIdsSize() ;
    int * ArrayOfSendRequests = new int[ size ] ;
    int nSendRequest = _MPI_access->sendRequestIds( size , ArrayOfSendRequests ) ;
    bool SendTrace = false ;
    int i ;
    for ( i = 0 ; i < nSendRequest ; i++ )
      {
        if ( WithWait )
          {
            if (SendTrace)
              {
                cout << "CheckSent" << _my_rank << " " << i << "./" << nSendRequest
                    << " SendRequestId " << ArrayOfSendRequests[i] << " MPITarget "
                    << _MPI_access->MPITarget(ArrayOfSendRequests[i]) << " MPITag "
                    << _MPI_access->MPITag(ArrayOfSendRequests[i]) << " Wait :" << endl ;
              }
            sts = _MPI_access->wait( ArrayOfSendRequests[i] ) ;
          }
        else
          {
            sts = _MPI_access->test( ArrayOfSendRequests[i] , flag ) ;
          }
        if ( flag )
          {
            _MPI_access->deleteRequest( ArrayOfSendRequests[i] ) ;
            if ( SendTrace )
              {
                cout << "CheckSent" << _my_rank << " " << i << "./" << nSendRequest
                     << " SendRequestId " << ArrayOfSendRequests[i]
                     << " flag " << flag
                     << " Counter " << (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->Counter
                     << " DataType " << (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->DataType
                     << endl ;
              }
            (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->Counter -= 1 ;
            if ( SendTrace )
              {
                if ( (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->DataType == 
                     _MPI_access->timeType() )
                  {
                    cout << "CheckTimeSent" << _my_rank << " Request " ;
                  }
                else
                  {
                    cout << "CheckDataSent" << _my_rank << " Request " ;
                  }
                cout << ArrayOfSendRequests[i]
                     << " _map_of_send_buffers->SendBuffer "
                     << (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->SendBuffer
                     << " Counter " << (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->Counter
                     << endl ;
              }
            if ( (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->Counter  == 0 )
              {
                if ( SendTrace )
                  {
                    cout << "CheckSent" << _my_rank << " SendRequestId " << ArrayOfSendRequests[i]
                         << " Counter " << (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->Counter
                         << " flag " << flag << " SendBuffer "
                         << (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->SendBuffer
                         << " deleted. Erase in _map_of_send_buffers :" << endl ;
                  }
                if ( (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->DataType ==
                     _MPI_access->timeType() )
                  {
                    delete (TimeMessage * ) (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->SendBuffer ;
                  }
                else
                  {
                    if ( (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->DataType == MPI_INT )
                      {
                        delete [] (int *) (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->SendBuffer ;
                      }
                    else
                      {
                        delete [] (double *) (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->SendBuffer ;
                      }
                  }
                delete (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ] ;
              }
            if ( SendTrace )
              {
                cout << "CheckSent" << _my_rank << " Erase in _map_of_send_buffers SendRequestId "
                     << ArrayOfSendRequests[i] << endl ;
              }
            (*_map_of_send_buffers).erase( ArrayOfSendRequests[i] ) ;
          }
        else if ( SendTrace )
          {
            cout << "CheckSent" << _my_rank << " " << i << "./" << nSendRequest
                 << " SendRequestId " << ArrayOfSendRequests[i]
                 << " flag " << flag
                 << " Counter " << (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->Counter
                 << " DataType " << (*_map_of_send_buffers)[ ArrayOfSendRequests[i] ]->DataType
                 << endl ;
          }
      }
    if ( SendTrace )
      {
        _MPI_access->check() ;
      }
    delete [] ArrayOfSendRequests ;
    return sts ;
  }

  int MPIAccessDEC::checkFinalRecv()
  {
    int sts = MPI_SUCCESS ;
    if ( _time_interpolator )
      {
        int target ;
        for ( target = 0 ; target < _group_size ; target++ )
          {
            if ( (*_data_messages)[target][0] != NULL )
              {
                sts = checkTime( (*_data_messages_recv_count)[target] , (*_data_messages_type)[target] ,
                                 target , true ) ;
                if ( (*_data_messages_type)[target] == MPI_INT )
                  {
                    delete [] (int *) (*_data_messages)[target][0] ;
                  }
                else
                  {
                    delete [] (double *) (*_data_messages)[target][0] ;
                  }
                (*_data_messages)[target][0] = NULL ;
                if ( (*_data_messages)[target][1] != NULL )
                  {
                    if ( (*_data_messages_type)[target] == MPI_INT )
                      {
                        delete [] (int *) (*_data_messages)[target][1] ;
                      }
                    else
                      {
                        delete [] (double *) (*_data_messages)[target][1] ;
                      }
                    (*_data_messages)[target][1] = NULL ;
                  }
              }
          }
      }
    return sts ;
  }

  ostream & operator<< (ostream & f ,const TimeInterpolationMethod & interpolationmethod )
  {
    switch (interpolationmethod)
      {
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
