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

#include "MPIAccess.hxx"
#include "InterpolationUtils.hxx"

#include <iostream>

using namespace std;

namespace MEDCoupling
{
  /**!
    \anchor MPIAccess-det
    \class MPIAccess

    The class \a MPIAccess is the gateway to the MPI library.
    It is a helper class that gathers the calls to the MPI
    library that are made in the ParaMEDMEM library. This gathering
    allows easier gathering of information about the communication
    in the library. With MPIAccess, tags are managed automatically
    and asynchronous operations are easier.

    It is typically called after the MPI_Init() call in a program. It is afterwards passed as a parameter to the constructors of ParaMEDMEM objects so that they access the MPI library via the MPIAccess.

    As an example, the following code initializes a processor group made of the zero processor.

    \verbatim
    #include "MPIAccess.hxx"
    #include "ProcessorGroup.hxx"

    int main(int argc, char** argv)
    {
    //initialization
    MPI_Init(&argc, &argv);
    MEDCoupling::CommInterface comm_interface;

    //setting up a processor group with proc 0
    set<int> procs;
    procs.insert(0);
    MEDCoupling::ProcessorGroup group(procs, comm_interface);

    MEDCoupling::MPIAccess mpi_access(group);

    //cleanup
    MPI_Finalize();
    }
    \endverbatim
  */


  /*! Creates a MPIAccess that is based on the processors included in \a ProcessorGroup.
    This class may be called for easier use of MPI API.

    \param ProcessorGroup MPIProcessorGroup object giving access to group management
    \param BaseTag and MaxTag define the range of tags to be used.
    Tags are managed by MPIAccess. They are cyclically incremented.
    When there is a Send or a Receive operation there is a new RequestId tag returned
    to the caller. That RequestId may be used to manage the operation Wait, Check of
    status etc... The MPITag internally managed by MPIAccess is used as "tag" argument
    in MPI call.
  */

  MPIAccess::MPIAccess(MPIProcessorGroup * ProcessorGroup, int BaseTag, int MaxTag) :
    _comm_interface( ProcessorGroup->getCommInterface() ) ,
    _intra_communicator( ProcessorGroup->getComm() )
  {
    void *v ;
    int mpitagub ;
    int flag ;
    //MPI_Comm_get_attr does not run with _IntraCommunicator ???
    //MPI_Comm_get_attr(*_IntraCommunicator,MPID_TAG_UB,&mpitagub,&flag) ;
    MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_TAG_UB,&v,&flag) ;
    mpitagub=*(reinterpret_cast<int*>(v));
    if ( BaseTag != 0 )
      BaseTag = (BaseTag/MODULO_TAG)*MODULO_TAG ;
    if ( MaxTag == 0 )
      MaxTag = (mpitagub/MODULO_TAG-1)*MODULO_TAG ;
    MPI_Comm_rank( *_intra_communicator, &_my_rank ) ;
    if ( !flag | (BaseTag < 0) | (BaseTag >= MaxTag) | (MaxTag > mpitagub) )
      throw INTERP_KERNEL::Exception("wrong call to MPIAccess constructor");

    _processor_group = ProcessorGroup ;
    _processor_group_size = _processor_group->size() ;
    _trace = false ;

    _base_request = -1 ;
    _max_request = std::numeric_limits<int>::max() ;
    _request = _base_request ;
    
    _base_MPI_tag = BaseTag ;
    _max_MPI_tag = MaxTag ;
    
    _send_request = new int[ _processor_group_size ] ;
    _recv_request = new int[ _processor_group_size ] ;

    _send_requests.resize( _processor_group_size ) ;
    _recv_requests.resize( _processor_group_size ) ;

    _send_MPI_tag = new int[ _processor_group_size ] ;
    _recv_MPI_Tag = new int[ _processor_group_size ] ;
    int i ;
    for (i = 0 ; i < _processor_group_size ; i++ )
      {
        _send_request[ i ] = _max_request ;
        _recv_request[ i ] = _max_request ;
        _send_requests[ i ].resize(0) ;
        _recv_requests[ i ].resize(0) ;
        _send_MPI_tag[ i ] = _max_MPI_tag ;
        _recv_MPI_Tag[ i ] = _max_MPI_tag ;
      }
    MPI_Datatype array_of_types[3] ;
    array_of_types[0] = MPI_DOUBLE ;
    array_of_types[1] = MPI_DOUBLE ;
    array_of_types[2] = MPI_INT ;
    int array_of_blocklengths[3] ;
    array_of_blocklengths[0] = 1 ;
    array_of_blocklengths[1] = 1 ;
    array_of_blocklengths[2] = 1 ;
    MPI_Aint array_of_displacements[3] ;
    array_of_displacements[0] = 0 ;
    array_of_displacements[1] = sizeof(double) ;
    array_of_displacements[2] = 2*sizeof(double) ;
    MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements,
                    array_of_types, &_MPI_TIME) ;
    MPI_Type_commit(&_MPI_TIME) ;
  }

  MPIAccess::~MPIAccess()
  {
    delete [] _send_request ;
    delete [] _recv_request ;
    delete [] _send_MPI_tag ;
    delete [] _recv_MPI_Tag ;
    MPI_Type_free(&_MPI_TIME) ;
  }

  /*
    MPIAccess and "RequestIds" :
    ============================

    . WARNING : In the specification document, the distinction
    between "MPITags" and "RequestIds" is not clear. "MPITags"
    are arguments of calls to MPI. "RequestIds" does not concern
    calls to MPI. "RequestIds" are named "tag"as arguments in/out
    in the MPIAccess API in the specification documentation.
    But in the implementation we have the right name RequestId (or
    RecvRequestId/SendRequestId).

    . When we have a MPI write/read request via MPIAccess, we get
    an identifier "RequestId".
    That identifier matches a  structure RequestStruct of
    MPIAccess. The access to that structure is done with the map
    "_MapOfRequestStruct".
    That structure RequestStruct give the possibility to manage
    the structures MPI_Request and MPI_Status * of MPI. It give
    also the possibility to get informations about that request :
    target, send/recv, tag, [a]synchronous, type, outcount.

    . That identifier is used to control an asynchronous request
    via MPIAccess : Wait, Test, Probe, etc...

    . In practise "RequestId" is simply an integer fo the interval
    [0 , 2**32-1]. There is only one such a cyclic for
    [I]Sends and [I]Recvs.

    . That "RequestIds" and their associated structures give an easy
    way to manage asynchronous communications.
    For example we have mpi_access->Wait( int RequestId ) instead of
    MPI_Wait(MPI_Request *request, MPI_Status *status).

    . The API of MPIAccess may give the "SendRequestIds" of a "target",
    the "RecvRequestIds" from a "source" or the "SendRequestIds" of
    all "targets" or the "RecvRequestIds" of all "sources".
    That avoid to manage them in Presentation-ParaMEDMEM.
  */

  int MPIAccess::newRequest( MPI_Datatype datatype, int tag , int destsourcerank ,
                             bool fromsourcerank , bool asynchronous )
  {
    RequestStruct *mpiaccessstruct = new RequestStruct;
    mpiaccessstruct->MPITag = tag ;
    mpiaccessstruct->MPIDatatype = datatype ;
    mpiaccessstruct->MPITarget = destsourcerank ;
    mpiaccessstruct->MPIIsRecv = fromsourcerank ;
    MPI_Status *aStatus = new MPI_Status ;
    mpiaccessstruct->MPIStatus = aStatus ;
    mpiaccessstruct->MPIAsynchronous = asynchronous ;
    mpiaccessstruct->MPICompleted = !asynchronous ;
    mpiaccessstruct->MPIOutCount = -1 ;
    if ( !asynchronous )
      {
        mpiaccessstruct->MPIRequest = MPI_REQUEST_NULL ;
        mpiaccessstruct->MPIStatus->MPI_SOURCE = destsourcerank ;
        mpiaccessstruct->MPIStatus->MPI_TAG = tag ;
        mpiaccessstruct->MPIStatus->MPI_ERROR = MPI_SUCCESS ;
      }
    if ( _request == _max_request )
      _request = _base_request ;
    _request += 1 ;
    _map_of_request_struct[_request] = mpiaccessstruct ;
    if ( fromsourcerank )
      _recv_request[ destsourcerank ] = _request;
    else
      _send_request[ destsourcerank ] = _request;
    if ( _trace )
      cout << "NewRequest" << _my_rank << "( " << _request << " ) "
           << mpiaccessstruct << endl ;
    return _request ;
  }

  /*
    MPIAccess and "tags" (or "MPITags") :
    =====================================

    . The constructor give the possibility to choose an interval of
    tags to use : [BaseTag , MaxTag].
    The default is [ 0 , MPI_TAG_UB], MPI_TAG_UB being the maximum
    value in an implementation of MPI (minimum 32767 = 2**15-1).
    On awa with the implementation lam MPI_TAG_UB value is
    7353944. The norm MPI specify that value is the same in all
    processes started by mpirun.
    In the case of the use of the same IntraCommunicator in a process
    for several distinct data flows (or for several IntraCommunicators
    with common processes), that permits to avoid ambiguity
    and may help debug.

    . In MPIAccess the tags have two parts (#define MODULO_TAG 10) :
    + The last decimal digit decimal correspond to MPI_DataType ( 1 for
    TimeMessages, 2 for MPI_INT and 3 for MPI_DOUBLE)
    + The value of other digits correspond to a circular number for each
    message.
    + A TimeMessage and the associated DataMessage have the same number
    (but the types are different and the tags also).

    . For a Send of a message from a process "source" to a process
    "target", we have _send_MPI_tag[target] in the process
    source (it contains the last "tag" used for the Send of a
    message to the process target).
    And in the process "target" which receive that message, we have
    _recv_MPI_Tag[source] (it contains the last "tag" used for the Recv
    of messages from the process source).
    Naturally in the MPI norm the values of that tags must be the same.
  */
  int MPIAccess::newSendTag( MPI_Datatype datatype, int destrank , int method ,
                             bool asynchronous, int &RequestId )
  {
    int tag ;
    tag = incrTag( _send_MPI_tag[destrank] ) ;
    tag = valTag( tag, method ) ;
    _send_MPI_tag[ destrank ] = tag ;
    RequestId = newRequest( datatype, tag, destrank , false , asynchronous ) ;
    _send_request[ destrank ] = RequestId ;
    _send_requests[ destrank ].push_back( RequestId ) ;
    return tag ;
  }

  int MPIAccess::newRecvTag( MPI_Datatype datatype, int sourcerank , int method ,
                             bool asynchronous, int &RequestId )
  {
    int tag ;
    tag = incrTag( _recv_MPI_Tag[sourcerank] ) ;
    tag = valTag( tag, method ) ;
    _recv_MPI_Tag[ sourcerank ] = tag ;
    RequestId = newRequest( datatype, tag , sourcerank , true , asynchronous ) ;
    _recv_request[ sourcerank ] = RequestId ;
    _recv_requests[ sourcerank ].push_back( RequestId ) ;
    return tag ;
  }

  // Returns the number of all SendRequestIds that may be used to allocate
  // ArrayOfSendRequests for the call to SendRequestIds
  int MPIAccess::sendRequestIdsSize()
  {
    int size = 0;
    for (int i = 0 ; i < _processor_group_size ; i++ )
      size += _send_requests[ i ].size() ;
    return size ;
  }

  // Returns in ArrayOfSendRequests with the dimension "size" all the
  // SendRequestIds
  int MPIAccess::sendRequestIds(int size, int *ArrayOfSendRequests)
  {
    int destrank ;
    int i = 0 ;
    for ( destrank = 0 ; destrank < _processor_group_size ; destrank++ )
      {
        list< int >::const_iterator iter ;
        for (iter = _send_requests[ destrank ].begin() ; iter != _send_requests[destrank].end() ; iter++ )
          ArrayOfSendRequests[i++] = *iter ;
      }
    return i ;
  }

  // Returns the number of all RecvRequestIds that may be used to allocate
  // ArrayOfRecvRequests for the call to RecvRequestIds
  int MPIAccess::recvRequestIdsSize()
  {
    int size = 0 ;
    for (int i = 0 ; i < _processor_group_size ; i++ )
      size += _recv_requests[ i ].size() ;
    return size ;
  }

  // Returns in ArrayOfRecvRequests with the dimension "size" all the
  // RecvRequestIds
  int MPIAccess::recvRequestIds(int size, int *ArrayOfRecvRequests)
  {
    int sourcerank ;
    int i = 0 ;
    for ( sourcerank = 0 ; sourcerank < _processor_group_size ; sourcerank++ )
      {
        list< int >::const_iterator iter ;
        for (iter = _recv_requests[ sourcerank ].begin() ; iter != _recv_requests[sourcerank].end() ; iter++ )
          ArrayOfRecvRequests[i++] = *iter ;
      }
    return i ;
  }

  // Returns in ArrayOfSendRequests with the dimension "size" all the
  // SendRequestIds to a destination rank
  int MPIAccess::sendRequestIds(int destrank, int size, int *ArrayOfSendRequests)
  {
    if (size < (int)_send_requests[destrank].size() )
      throw INTERP_KERNEL::Exception("wrong call to MPIAccess::SendRequestIds");
    int i = 0 ;
    list< int >::const_iterator iter ;
    for (iter = _send_requests[ destrank ].begin() ; iter != _send_requests[destrank].end() ; iter++ )
      ArrayOfSendRequests[i++] = *iter ;
    return _send_requests[destrank].size() ;
  }

  // Returns in ArrayOfRecvRequests with the dimension "size" all the
  // RecvRequestIds from a sourcerank
  int MPIAccess::recvRequestIds(int sourcerank, int size, int *ArrayOfRecvRequests)
  {
    if (size < (int)_recv_requests[sourcerank].size() )
      throw INTERP_KERNEL::Exception("wrong call to MPIAccess::RecvRequestIds");
    int i = 0 ;
    list< int >::const_iterator iter ;
    _recv_requests[ sourcerank ] ;
    for (iter = _recv_requests[ sourcerank ].begin() ; iter != _recv_requests[sourcerank].end() ; iter++ )
      ArrayOfRecvRequests[i++] = *iter ;
    return _recv_requests[sourcerank].size() ;
  }

  // Send in synchronous mode count values of type datatype from buffer to target
  // (returns RequestId identifier even if the corresponding structure is deleted :
  // it is only in order to have the same signature as the asynchronous mode)
  int MPIAccess::send(void* buffer, int count, MPI_Datatype datatype, int target, int &RequestId)
  {
    int sts = MPI_SUCCESS ;
    RequestId = -1 ;
    if ( count )
      {
        _MessageIdent aMethodIdent = methodId( datatype ) ;
        int MPItag = newSendTag( datatype, target , aMethodIdent , false , RequestId ) ;
        if ( aMethodIdent == _message_time )
          {
            TimeMessage *aTimeMsg = (TimeMessage *) buffer ;
            aTimeMsg->tag = MPItag ;
          }
        deleteRequest( RequestId ) ;
        sts = _comm_interface.send(buffer, count, datatype, target, MPItag,
                                  *_intra_communicator ) ;
        if ( _trace )
          cout << "MPIAccess::Send" << _my_rank << " SendRequestId "
               << RequestId << " count " << count << " target " << target
               << " MPItag " << MPItag << endl ;
      }
    return sts ;
  }

  // Receive (read) in synchronous mode count values of type datatype in buffer from source
  // (returns RequestId identifier even if the corresponding structure is deleted :
  // it is only in order to have the same signature as the asynchronous mode)
  // The output argument OutCount is optionnal : *OutCount <= count
  int MPIAccess::recv(void* buffer, int count, MPI_Datatype datatype, int source, int &RequestId, int *OutCount)
  {
    int sts = MPI_SUCCESS ;
    RequestId = -1 ;
    if ( OutCount != NULL )
      *OutCount = -1 ;
    if ( count )
      {
        _MessageIdent aMethodIdent = methodId( datatype ) ;
        int MPItag = newRecvTag( datatype, source , aMethodIdent , false , RequestId ) ;
        sts =  _comm_interface.recv(buffer, count, datatype, source, MPItag,
                                   *_intra_communicator , MPIStatus( RequestId ) ) ;
        int outcount = 0 ;
        if ( sts == MPI_SUCCESS )
          {
            MPI_Datatype datatype = MPIDatatype( RequestId ) ;
            _comm_interface.getCount(MPIStatus( RequestId ), datatype, &outcount ) ;
            setMPIOutCount( RequestId , outcount ) ;
            setMPICompleted( RequestId , true ) ;
            deleteStatus( RequestId ) ;
          }
        if ( OutCount != NULL )
          *OutCount = outcount ;
        if ( _trace )
          cout << "MPIAccess::Recv" << _my_rank << " RecvRequestId "
               << RequestId << " count " << count << " source " << source
               << " MPItag " << MPItag << endl ;
        deleteRequest( RequestId ) ;
      }
    return sts ;
  }

  // Send in asynchronous mode count values of type datatype from buffer to target
  // Returns RequestId identifier.
  int MPIAccess::ISend(void* buffer, int count, MPI_Datatype datatype, int target, int &RequestId)
  {
    int sts = MPI_SUCCESS ;
    RequestId = -1 ;
    if ( count )
      {
        _MessageIdent aMethodIdent = methodId( datatype ) ;
        int MPItag = newSendTag( datatype, target , aMethodIdent , true , RequestId ) ;
        if ( aMethodIdent == _message_time )
          {
            TimeMessage *aTimeMsg = (TimeMessage *) buffer ;
            aTimeMsg->tag = MPItag ;
          }
        MPI_Request *aSendRequest = MPIRequest( RequestId ) ;
        if ( _trace )
          {
            cout << "MPIAccess::ISend" << _my_rank << " ISendRequestId "
                 << RequestId << " count " << count << " target " << target
                 << " MPItag " << MPItag << endl ;
            if ( MPItag == 1 )
              cout << "MPIAccess::ISend" << _my_rank << " time "
                   << ((TimeMessage *)buffer)->time << " " << ((TimeMessage *)buffer)->deltatime
                   << endl ;
          }
        sts = _comm_interface.Isend(buffer, count, datatype, target, MPItag,
                                   *_intra_communicator , aSendRequest) ;
      }
    return sts ;
  }

  // Receive (read) in asynchronous mode count values of type datatype in buffer from source
  // returns RequestId identifier.
  int MPIAccess::IRecv(void* buffer, int count, MPI_Datatype datatype, int source, int &RequestId)
  {
    int sts = MPI_SUCCESS ;
    RequestId = -1 ;
    if ( count )
      {
        _MessageIdent aMethodIdent = methodId( datatype ) ;
        int MPItag = newRecvTag( datatype, source , aMethodIdent , true , RequestId ) ;
        MPI_Request *aRecvRequest = MPIRequest( RequestId ) ;
        if ( _trace )
          {
            cout << "MPIAccess::IRecv" << _my_rank << " IRecvRequestId "
                 << RequestId << " count " << count << " source " << source
                 << " MPItag " << MPItag << endl ;
            if ( MPItag == 1 )
              cout << "MPIAccess::ISend" << _my_rank << " time "
                   << ((TimeMessage *)buffer)->time << " " << ((TimeMessage *)buffer)->deltatime
                   << endl ;
          }
        sts = _comm_interface.Irecv(buffer, count, datatype, source, MPItag,
                                   *_intra_communicator , aRecvRequest) ;
      }
    return sts ;
  }

  // Perform a Send and a Recv in synchronous mode
  int MPIAccess::sendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                          int dest, int &SendRequestId,
                          void* recvbuf, int recvcount, MPI_Datatype recvtype,
                          int source, int &RecvRequestId, int *OutCount)
  {
    int sts = MPI_SUCCESS ;
    SendRequestId = -1 ;
    RecvRequestId = -1 ;
    if ( recvcount )
      sts = IRecv(recvbuf, recvcount, recvtype, source, RecvRequestId) ;
    int outcount = -1 ;
    if ( _trace )
      cout << "MPIAccess::SendRecv" << _my_rank << " IRecv RecvRequestId "
           << RecvRequestId << endl ;
    if ( sts == MPI_SUCCESS )
      {
        if ( sendcount )
          sts = send(sendbuf, sendcount, sendtype, dest, SendRequestId) ;
        if ( _trace )
          cout << "MPIAccess::SendRecv" << _my_rank << " Send SendRequestId "
               << SendRequestId << endl ;
        if ( sts == MPI_SUCCESS && recvcount )
          {
            sts = wait( RecvRequestId ) ;
            outcount = MPIOutCount( RecvRequestId ) ;
            if ( _trace )
              cout << "MPIAccess::SendRecv" << _my_rank << " IRecv RecvRequestId "
                   << RecvRequestId << " outcount " << outcount << endl ;
          }
      }
    if ( OutCount != NULL )
      {
        *OutCount = outcount ;
        if ( _trace )
          cout << "MPIAccess::SendRecv" << _my_rank << " *OutCount = " << *OutCount
               << endl ;
      }
    deleteRequest( RecvRequestId ) ;
    return sts ;
  }

  // Perform a Send and a Recv in asynchronous mode
  int MPIAccess::ISendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                           int dest, int &SendRequestId,
                           void* recvbuf, int recvcount, MPI_Datatype recvtype,
                           int source, int &RecvRequestId)
  {
    int sts = MPI_SUCCESS ;
    SendRequestId = -1 ;
    RecvRequestId = -1 ;
    if ( recvcount )
      sts = IRecv(recvbuf, recvcount, recvtype, source, RecvRequestId) ;
    if ( sts == MPI_SUCCESS )
      if ( sendcount )
        sts = ISend(sendbuf, sendcount, sendtype, dest, SendRequestId) ;
    return sts ;
  }

  // Perform a wait of a Send or Recv asynchronous Request
  // Do nothing for a synchronous Request
  // Manage MPI_Request * and MPI_Status * structure
  int MPIAccess::wait( int RequestId )
  {
    int status = MPI_SUCCESS ;
    if ( !MPICompleted( RequestId ) )
      {
        if ( *MPIRequest( RequestId ) != MPI_REQUEST_NULL )
          {
            if ( _trace )
              cout << "MPIAccess::Wait" << _my_rank << " -> wait( " << RequestId
                   << " ) MPIRequest " << MPIRequest( RequestId ) << " MPIStatus "
                   << MPIStatus( RequestId ) << " MPITag " << MPITag( RequestId )
                   << " MPIIsRecv " << MPIIsRecv( RequestId ) << endl ;
            status = _comm_interface.wait(MPIRequest( RequestId ), MPIStatus( RequestId )) ;
          }
        else
          {
            if ( _trace )
              cout << "MPIAccess::Wait" << _my_rank << " MPIRequest == MPI_REQUEST_NULL"
                   << endl ;
          }
        setMPICompleted( RequestId , true ) ;
        if ( MPIIsRecv( RequestId ) && MPIStatus( RequestId ) )
          {
            MPI_Datatype datatype = MPIDatatype( RequestId ) ;
            int outcount ;
            status = _comm_interface.getCount(MPIStatus( RequestId ), datatype,
                                             &outcount ) ;
            if ( status == MPI_SUCCESS )
              {
                setMPIOutCount( RequestId , outcount ) ;
                deleteStatus( RequestId ) ;
                if ( _trace )
                  cout << "MPIAccess::Wait" << _my_rank << " RequestId " << RequestId
                       << "MPIIsRecv " << MPIIsRecv( RequestId ) << " outcount " << outcount
                       << endl ;
              }
            else
              {
                if ( _trace )
                  cout << "MPIAccess::Wait" << _my_rank << " MPIIsRecv "
                       << MPIIsRecv( RequestId ) << " outcount " << outcount << endl ;
              }
          }
        else
          {
            if ( _trace )
              cout << "MPIAccess::Wait" << _my_rank << " MPIIsRecv " << MPIIsRecv( RequestId )
                   << " MPIOutCount " << MPIOutCount( RequestId ) << endl ;
          }
      }
    if ( _trace )
      cout << "MPIAccess::Wait" << _my_rank << " RequestId " << RequestId
           << " Request " << MPIRequest( RequestId )
           << " Status " << MPIStatus( RequestId ) << " MPICompleted "
           << MPICompleted( RequestId ) << " MPIOutCount " << MPIOutCount( RequestId )
           << endl ;
    return status ;
  }

  // Perform a "test" of a Send or Recv asynchronous Request
  // If the request is done, returns true in the flag argument
  // If the request is not finished, returns false in the flag argument
  // Do nothing for a synchronous Request
  // Manage MPI_request * and MPI_status * structure
  int MPIAccess::test(int RequestId, int &flag)
  {
    int status = MPI_SUCCESS ;
    flag = MPICompleted( RequestId ) ;
    if ( _trace )
      cout << "MPIAccess::Test" << _my_rank << " flag " << flag ;
    if ( MPIIsRecv( RequestId ) )
      {
        if ( _trace )
          cout << " Recv" ;
      }
    else
      {
        if ( _trace )
          cout << " Send" ;
      }
    if( _trace )
      cout << "Request" << RequestId << " " << MPIRequest( RequestId )
           << " Status " << MPIStatus( RequestId ) << endl ;
    if ( !flag )
      {
        if ( *MPIRequest( RequestId ) != MPI_REQUEST_NULL )
          {
            if ( _trace )
              cout << "MPIAccess::Test" << _my_rank << " -> test( " << RequestId
                   << " ) MPIRequest " << MPIRequest( RequestId )
                   << " MPIStatus " << MPIStatus( RequestId )
                   << " MPITag " << MPITag( RequestId )
                   << " MPIIsRecv " << MPIIsRecv( RequestId ) << endl ;
            status = _comm_interface.test(MPIRequest( RequestId ), &flag,
                                         MPIStatus( RequestId )) ;
          }
        else
          {
            if ( _trace )
              cout << "MPIAccess::Test" << _my_rank << " MPIRequest == MPI_REQUEST_NULL"
                   << endl ;
          }
        if ( flag )
          {
            setMPICompleted( RequestId , true ) ;
            if ( MPIIsRecv( RequestId ) && MPIStatus( RequestId ) )
              {
                int outcount ;
                MPI_Datatype datatype = MPIDatatype( RequestId ) ;
                status = _comm_interface.getCount( MPIStatus( RequestId ), datatype,
                                                  &outcount ) ;
                if ( status == MPI_SUCCESS )
                  {
                    setMPIOutCount( RequestId , outcount ) ;
                    deleteStatus( RequestId ) ;
                    if ( _trace )
                      cout << "MPIAccess::Test" << _my_rank << " MPIIsRecv "
                           << MPIIsRecv( RequestId ) << " outcount " << outcount << endl ;
                  }
                else
                  {
                    if ( _trace )
                      cout << "MPIAccess::Test" << _my_rank << " MPIIsRecv "
                           << MPIIsRecv( RequestId ) << " outcount " << outcount << endl ;
                  }
              }
            else
              {
                if ( _trace )
                  cout << "MPIAccess::Test" << _my_rank << " MPIIsRecv "
                       << MPIIsRecv( RequestId ) << " MPIOutCount "
                       << MPIOutCount( RequestId ) << endl ;
              }
          }
      }
    if ( _trace )
      cout << "MPIAccess::Test" << _my_rank << " RequestId " << RequestId
           << " flag " << flag << " MPICompleted " << MPICompleted( RequestId )
           << " MPIOutCount " << MPIOutCount( RequestId ) << endl ;
    return status ;
  }

  int MPIAccess::waitAny(int count, int *array_of_RequestIds, int &RequestId)
  {
    int status = MPI_ERR_OTHER ;
    RequestId = -1 ;
    cout << "MPIAccess::WaitAny not yet implemented" << endl ;
    return status ;
  }

  int MPIAccess::testAny(int count, int *array_of_RequestIds, int &RequestId, int &flag)
  {
    int status = MPI_ERR_OTHER ;
    RequestId = -1 ;
    flag = 0 ;
    cout << "MPIAccess::TestAny not yet implemented" << endl ;
    return status ;
  }
  
  // Perform a wait of each Send or Recv asynchronous Request of the array 
  // array_of_RequestIds of size "count".
  // That array may be filled with a call to SendRequestIdsSize or RecvRequestIdsSize
  // Do nothing for a synchronous Request
  // Manage MPI_Request * and MPI_Status * structure
  int MPIAccess::waitAll(int count, int *array_of_RequestIds)
  {
    if ( _trace )
      cout << "WaitAll" << _my_rank << " : count " << count << endl ;
    int status ;
    int retstatus = MPI_SUCCESS ;
    int i ;
    for ( i = 0 ; i < count ; i++ )
      {
        if ( _trace )
          cout << "WaitAll" << _my_rank << " " << i << " -> Wait( "
               << array_of_RequestIds[i] << " )" << endl ;
        status = wait( array_of_RequestIds[i] ) ;
        if ( status != MPI_SUCCESS )
          retstatus = status ;
      }
    if ( _trace )
      cout << "EndWaitAll" << _my_rank << endl ;
    return retstatus ;
  }

  // Perform a "test" of each Send or Recv asynchronous Request of the array 
  // array_of_RequestIds of size "count".
  // That array may be filled with a call to SendRequestIdsSize or RecvRequestIdsSize
  // If all requests are done, returns true in the flag argument
  // If all requests are not finished, returns false in the flag argument
  // Do nothing for a synchronous Request
  // Manage MPI_Request * and MPI_Status * structure
  int MPIAccess::testAll(int count, int *array_of_RequestIds, int &flag)
  {
    if ( _trace )
      cout << "TestAll" << _my_rank << " : count " << count << endl ;
    int status ;
    int retstatus = MPI_SUCCESS ;
    bool retflag = true ;
    int i ;
    for ( i = 0 ; i < count ; i++ )
      {
        status = test( array_of_RequestIds[i] , flag ) ;
        retflag = retflag && (flag != 0) ;
        if ( status != MPI_SUCCESS )
          retstatus = status ;
      }
    flag = retflag ;
    if ( _trace )
      cout << "EndTestAll" << _my_rank << endl ;
    return retstatus ;
  }

  int MPIAccess::waitSome(int count, int *array_of_RequestIds, int outcount,
                          int *outarray_of_RequestIds)
  {
    int status = MPI_ERR_OTHER ;
    cout << "MPIAccess::WaitSome not yet implemented" << endl ;
    return status ;
  }

  int MPIAccess::testSome(int count, int *array_of_RequestIds, int outcounts,
                          int *outarray_of_RequestIds)
  {
    int status = MPI_ERR_OTHER ;
    cout << "MPIAccess::TestSome not yet implemented" << endl ;
    return status ;
  }
  
  // Probe checks if a message is available for read from FromSource rank.
  // Returns the corresponding source, MPITag, datatype and outcount
  // Probe is a blocking call which wait until a message is available
  int MPIAccess::probe(int FromSource, int &source, int &MPITag,
                       MPI_Datatype &myDatatype, int &outcount)
  {
    MPI_Status aMPIStatus ;
    int sts =  _comm_interface.probe( FromSource, MPI_ANY_TAG,
                                     *_intra_communicator , &aMPIStatus ) ;
    if ( sts == MPI_SUCCESS )
      {
        source = aMPIStatus.MPI_SOURCE ;
        MPITag = aMPIStatus.MPI_TAG ;
        int MethodId = (MPITag % MODULO_TAG) ;
        myDatatype = datatype( (MEDCoupling::_MessageIdent) MethodId ) ;
        _comm_interface.getCount(&aMPIStatus, myDatatype, &outcount ) ;
        if ( _trace )
          cout << "MPIAccess::Probe" << _my_rank << " FromSource " << FromSource
               << " source " << source << " MPITag " << MPITag << " MethodId "
               << MethodId << " datatype " << myDatatype << " outcount " << outcount
               << endl ;
      }
    else
      {
        source = -1 ;
        MPITag = -1 ;
        myDatatype = 0 ;
        outcount = -1 ;
      }
    return sts ;
  }

  // IProbe checks if a message is available for read from FromSource rank.
  // If there is a message available, returns the corresponding source,
  // MPITag, datatype and outcount with flag = true
  // If not, returns flag = false
  int MPIAccess::IProbe(int FromSource, int &source, int &MPITag,
                        MPI_Datatype &myDataType, int &outcount, int &flag)
  {
    MPI_Status aMPIStatus ;
    int sts =  _comm_interface.Iprobe( FromSource, MPI_ANY_TAG,
                                      *_intra_communicator , &flag,
                                      &aMPIStatus ) ;
    if ( sts == MPI_SUCCESS && flag )
      {
        source = aMPIStatus.MPI_SOURCE ;
        MPITag = aMPIStatus.MPI_TAG ;
        int MethodId = (MPITag % MODULO_TAG) ;
        myDataType = datatype( (MEDCoupling::_MessageIdent) MethodId ) ;
        _comm_interface.getCount(&aMPIStatus, myDataType, &outcount ) ;
        if ( _trace )
          cout << "MPIAccess::IProbe" << _my_rank << " FromSource " << FromSource
               << " source " << source << " MPITag " << MPITag << " MethodId "
               << MethodId << " datatype " << myDataType << " outcount " << outcount
               << " flag " << flag << endl ;
      }
    else
      {
        source = -1 ;
        MPITag = -1 ;
        myDataType = 0 ;
        outcount = -1 ;
      }
    return sts ;
  }

  // Cancel concerns a "posted" asynchronous IRecv
  // Returns flag = true if the receiving request was successfully canceled
  // Returns flag = false if the receiving request was finished but not canceled
  // Use cancel, wait and test_cancelled of the MPI API
  int MPIAccess::cancel( int RecvRequestId, int &flag )
  {
    flag = 0 ;
    int sts = _comm_interface.cancel( MPIRequest( RecvRequestId ) ) ;
    if ( sts == MPI_SUCCESS )
      {
        sts = _comm_interface.wait( MPIRequest( RecvRequestId ) ,
                                   MPIStatus( RecvRequestId ) ) ;
        if ( sts == MPI_SUCCESS )
          sts = _comm_interface.testCancelled( MPIStatus( RecvRequestId ) , &flag ) ;
      }
    return sts ;
  }

  // Cancel concerns a "pending" receiving message (without IRecv "posted")
  // Returns flag = true if the message was successfully canceled
  // Returns flag = false if the receiving request was finished but not canceled
  // Use Irecv, cancel, wait and test_cancelled of the MPI API
  int MPIAccess::cancel( int source, int theMPITag, MPI_Datatype datatype, int outcount, int &flag )
  {
    int sts ;
    MPI_Aint extent, lbound ;
    flag = 0 ;
    sts =  MPI_Type_get_extent( datatype , &lbound, &extent ) ;
    if ( sts == MPI_SUCCESS )
      {
        void * recvbuf = malloc( extent*outcount ) ;
        MPI_Request aRecvRequest ;
        if ( _trace )
          cout << "MPIAccess::Cancel" << _my_rank << " Irecv extent " << extent
               << " datatype " << datatype << " source " << source << " theMPITag "
               << theMPITag << endl ;
        sts = _comm_interface.Irecv( recvbuf, outcount, datatype, source, theMPITag,
                                    *_intra_communicator , &aRecvRequest ) ;
        if ( sts == MPI_SUCCESS )
          {
            sts = _comm_interface.cancel( &aRecvRequest ) ;
            if ( _trace )
              cout << "MPIAccess::Cancel" << _my_rank << " theMPITag " << theMPITag
                   << " cancel done" << endl ;
            if ( sts == MPI_SUCCESS )
              {
                MPI_Status aStatus ;
                if ( _trace )
                  cout << "MPIAccess::Cancel" << _my_rank << " wait" << endl ;
                sts = _comm_interface.wait( &aRecvRequest , &aStatus ) ;
                if ( sts == MPI_SUCCESS )
                  {
                    if ( _trace )
                      cout << "MPIAccess::Cancel" << _my_rank << " test_cancelled" << endl ;
                    sts = _comm_interface.testCancelled( &aStatus , &flag ) ;
                  }
              }
          }
        if ( _trace && datatype == timeType() )
          cout << "MPIAccess::Cancel" << _my_rank << " time "
               << ((TimeMessage *) recvbuf)->time << " "
               << ((TimeMessage *) recvbuf)->deltatime << endl ;
        free( recvbuf ) ;
      }
    if ( _trace )
      cout << "MPIAccess::Cancel" << _my_rank << " flag " << flag << endl ;
    return sts ;
  }


  // CancelAll concerns all "pending" receiving message (without IRecv "posted")
  // CancelAll use IProbe and Cancel (see obove)
  int MPIAccess::cancelAll()
  {
    int sts = MPI_SUCCESS ;
    int target ;
    int source ;
    int MPITag ;
    MPI_Datatype datatype ;
    int outcount ;
    int flag ;
    for ( target = 0 ; target < _processor_group_size ; target++ )
      {
        sts = IProbe(target, source, MPITag, datatype, outcount, flag) ;
        if ( sts == MPI_SUCCESS && flag )
          {
            sts = cancel(source, MPITag, datatype, outcount, flag) ;
            if ( _trace )
              cout << "MPIAccess::CancelAll" << _my_rank << " source " << source
                   << " MPITag " << MPITag << " datatype " << datatype
                   << " outcount " << outcount << " Cancel flag " << flag << endl ;
            if ( sts != MPI_SUCCESS )
              break ;
          }
        else if ( sts != MPI_SUCCESS )
          break ;
      }
    return sts ;
  }

  // Same as barrier of MPI API
  int MPIAccess::barrier()
  {
    int status = _comm_interface.barrier( *_intra_communicator ) ;
    return status ;
  }

  // Same as Error_string of MPI API
  int MPIAccess::errorString(int errorcode, char *string, int *resultlen) const
  {
    return _comm_interface.errorString( errorcode, string, resultlen) ;
  }
  
  // Returns source, tag, error and outcount corresponding to receiving RequestId
  // By default the corresponding structure of RequestId is deleted
  int MPIAccess::status(int RequestId, int &source, int &tag, int &error,
                        int &outcount, bool keepRequestStruct)
  {
    MPI_Status *myStatus = MPIStatus( RequestId ) ;
    if ( _trace )
      cout << "MPIAccess::status" << _my_rank << " RequestId " << RequestId
           << " status " << myStatus << endl ;
    if ( myStatus != NULL && MPIAsynchronous( RequestId ) &&
         MPICompleted( RequestId ) )
      {
        if ( MPIIsRecv( RequestId ) )
          {
            source = myStatus->MPI_SOURCE ;
            tag = myStatus->MPI_TAG ;
            error = myStatus->MPI_ERROR ;
            MPI_Datatype datatype = MPIDatatype( RequestId ) ;
            _comm_interface.getCount(myStatus, datatype, &outcount ) ;
            if ( _trace )
              cout << "MPIAccess::status" << _my_rank << " RequestId " << RequestId
                   << " status " << myStatus << " outcount " << outcount << endl ;
            setMPIOutCount( RequestId , outcount ) ;
          }
        else
          {
            source = MPITarget( RequestId ) ;
            tag = MPITag( RequestId ) ;
            error = 0 ;
            outcount = MPIOutCount( RequestId ) ;
          }
        if ( !keepRequestStruct )
          deleteRequest( RequestId ) ;
        return MPI_SUCCESS ;
      }
    else
      {
        source = MPITarget( RequestId ) ;
        tag = MPITag( RequestId ) ;
        error = 0 ;
        outcount = MPIOutCount( RequestId ) ;
      }
    return MPI_SUCCESS ;
  }
  
  int MPIAccess::requestFree( MPI_Request *request )
  {
    return _comm_interface.requestFree( request ) ;
  }
  
  // Print all informations of all known requests for debugging purpose
  void MPIAccess::check() const
  {
    int i = 0 ;
    map< int , RequestStruct * >::const_iterator MapOfRequestStructiterator ;
    cout << "MPIAccess::Check" << _my_rank << "_map_of_request_struct_size "
         << _map_of_request_struct.size() << endl ;
    for ( MapOfRequestStructiterator = _map_of_request_struct.begin() ;
          MapOfRequestStructiterator != _map_of_request_struct.end() ;
          MapOfRequestStructiterator++ )
      {
        if ( MapOfRequestStructiterator->second != NULL )
          {
            cout << "    Check" << _my_rank << " " << i << ". Request"
                 << MapOfRequestStructiterator->first << "-->" ;
            if ( (MapOfRequestStructiterator->second)->MPIAsynchronous )
              cout << "I" ;
            if ( (MapOfRequestStructiterator->second)->MPIIsRecv )
              cout << "Recv from " ;
            else
              cout << "Send to " ;
            cout << (MapOfRequestStructiterator->second)->MPITarget
                 << " MPITag " << (MapOfRequestStructiterator->second)->MPITag
                 << " DataType " << (MapOfRequestStructiterator->second)->MPIDatatype
                 << " Request " << (MapOfRequestStructiterator->second)->MPIRequest
                 << " Status " << (MapOfRequestStructiterator->second)->MPIStatus
                 << " Completed " << (MapOfRequestStructiterator->second)->MPICompleted
                 << endl ;
          }
        i++ ;
      }
  }

  // Returns the MPI size of a TimeMessage
  MPI_Aint MPIAccess::timeExtent() const
  {
    MPI_Aint aextent, lbound ;
    MPI_Type_get_extent( _MPI_TIME , &lbound, &aextent ) ;
    return aextent ;
  }

  // Returns the MPI size of a MPI_INT
  MPI_Aint MPIAccess::intExtent() const
  {
    MPI_Aint aextent, lbound ;
    MPI_Type_get_extent( MPI_INT , &lbound, &aextent ) ;
    return aextent ;
  }

  // Returns the MPI size of a MPI_DOUBLE
  MPI_Aint MPIAccess::doubleExtent() const
  {
    MPI_Aint aextent, lbound ;
    MPI_Type_get_extent( MPI_DOUBLE , &lbound, &aextent ) ;
    return aextent ;
  }

  // Outputs fields of a TimeMessage structure
  ostream & operator<< (ostream & f ,const TimeMessage & aTimeMsg )
  {
    f << " time " << aTimeMsg.time << " deltatime " << aTimeMsg.deltatime
      << " tag " << aTimeMsg.tag ;
    return f;
  }

  // Outputs the DataType coded in a Tag
  ostream & operator<< (ostream & f ,const _MessageIdent & methodtype )
  {
    switch (methodtype)
      {
      case _message_time :
        f << " MethodTime ";
        break;
      case _message_int :
        f << " MPI_INT ";
        break;
      case _message_double :
        f << " MPI_DOUBLE ";
        break;
      default :
        f << " UnknownMethodType ";
        break;
      }
    return f;
  }
}
