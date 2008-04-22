#include <iostream>

#include "MPI_Access.hxx"
#include "MEDMEM_Exception.hxx"

using namespace std;

namespace ParaMEDMEM
{
/*! \defgroup mpi_access MPI_Access
	Class \a MPI_Access is the gateway to the MPI library.
	It is a helper class that gathers the calls to the MPI
	library that are made in the ParaMEDMEM library. This gathering
	allows easier gathering of information about the communication
	in the library. With MPI_Access, tags are managed automatically
	and asynchronous operations are easier.

It is typically called after the MPI_Init() call in a program. It is afterwards passed as a parameter to the constructors of ParaMEDMEM objects so that they access the MPI library via the MPI_Access.

As an example, the following code initializes a processor group made of the zero processor.

\verbatim
#include "MPI_Access.hxx"
#include "ProcessorGroup.hxx"

int main(int argc, char** argv)
{
  //initialization
  MPI_Init(&argc, &argv);
  ParaMEDMEM::CommInterface comm_interface;

  //setting up a processor group with proc 0
	set<int> procs;
	procs.insert(0);
  ParaMEDMEM::ProcessorGroup group(procs, comm_interface);

  ParaMEDMEM::MPI_Access mpi_access(group);

	//cleanup
	MPI_Finalize();
}
\endverbatim
*/


	/*! Creates a MPI_Access that is based on the processors included in \a ProcessorGroup.
This class may be called for easier use of MPI API.

\param ProcessorGroup MPIProcessorGroup object giving access to group management
\param BaseTag and MaxTag define the range of tags to be used.
Tags are managed by MPI_Access. They are cyclically incremented.
When there is a Send or a Receive operation there is a new RequestId tag returned
to the caller. That RequestId may be used to manage the operation Wait, Check of
status etc... The MPITag internally managed by MPI_Access is used as "tag" argument
in MPI call.
	*/

MPI_Access::MPI_Access(MPIProcessorGroup * ProcessorGroup, int BaseTag, int MaxTag) :
    _CommInterface( ProcessorGroup->getCommInterface() ) ,
    _IntraCommunicator( ProcessorGroup->getComm() ) {
  int mpitagub ;
  int flag ;
//MPI_Attr_get does not run with _IntraCommunicator ???
  //MPI_Attr_get(*_IntraCommunicator,MPI_TAG_UB,&mpitagub,&flag) ;
  MPI_Attr_get(MPI_COMM_WORLD,MPI_TAG_UB,&mpitagub,&flag) ;
  if ( BaseTag != 0 ) {
    BaseTag = (BaseTag/ModuloTag)*ModuloTag ;
  }
  if ( MaxTag == 0 ) {
    MaxTag = (mpitagub/ModuloTag-1)*ModuloTag ;
  }
  MPI_Comm_rank( *_IntraCommunicator, &_MyRank ) ;
  cout << "MPI_Access::MPI_Access" << _MyRank << " this " << this << " BaseTag " << BaseTag
       << " MaxTag " << MaxTag << " mpitagub " << mpitagub << " (minimum 32767) "
       << " flag " << flag << endl ;
  if ( !flag | (BaseTag < 0) | (BaseTag >= MaxTag) | (MaxTag > mpitagub) )
    throw MEDMEM::MEDEXCEPTION("wrong call to MPI_Access constructor");

  _ProcessorGroup = ProcessorGroup ;
  _ProcessorGroupSize = _ProcessorGroup->size() ;
  _Trace = false ;
  
  cout << "MPI_Access::MPI_Access" << _MyRank << " _ProcessorGroupSize "
       << _ProcessorGroupSize << endl ;

  _BaseRequest = -1 ;
  _MaxRequest = 4294967295 ;
  _Request = _BaseRequest ;

  _BaseMPITag = BaseTag ;
  _MaxMPITag = MaxTag ;

  _SendRequest = new int[ _ProcessorGroupSize ] ;
  _RecvRequest = new int[ _ProcessorGroupSize ] ;

  _SendRequests.resize( _ProcessorGroupSize ) ;
  _RecvRequests.resize( _ProcessorGroupSize ) ;

  _SendMPITag = new int[ _ProcessorGroupSize ] ;
  _RecvMPITag = new int[ _ProcessorGroupSize ] ;
  int i ;
  for (i = 0 ; i < _ProcessorGroupSize ; i++ ) {
     _SendRequest[ i ] = _MaxRequest ;
     _RecvRequest[ i ] = _MaxRequest ;
     _SendRequests[ i ].resize(0) ;
     _RecvRequests[ i ].resize(0) ;
     _SendMPITag[ i ] = _MaxMPITag ;
     _RecvMPITag[ i ] = _MaxMPITag ;
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
  MPI_Type_struct(3, array_of_blocklengths, array_of_displacements,
                  array_of_types, &_MPI_TIME) ;
  MPI_Type_commit(&_MPI_TIME) ;
}

MPI_Access::~MPI_Access() {
  cout << "MPI_Access::~MPI_Access" << _MyRank << " this " << this << endl ;
  delete [] _SendRequest ;
  delete [] _RecvRequest ;
  delete [] _SendMPITag ;
  delete [] _RecvMPITag ;
  cout << "End of MPI_Access::~MPI_Access" << _MyRank << " this " << this << endl ;
}

/*
MPI_Access and "RequestIds" :
============================

. WARNING : In the specification document, the distinction
  between "MPITags" and "RequestIds" is not clear. "MPITags"
  are arguments of calls to MPI. "RequestIds" does not concern
  calls to MPI. "RequestIds" are named "tag"as arguments in/out
  in the MPI_Access API in the specification documentation.
  But in the implementation we have the right name RequestId (or
  RecvRequestId/SendRequestId).

. When we have a MPI write/read request via MPI_Access, we get
    an identifier "RequestId".
  That identifier matches a  structure RequestStruct of
    MPI_Access. The access to that structure is done with the map
    "_MapOfRequestStruct".
  That structure RequestStruct give the possibility to manage
    the structures MPI_Request and MPI_Status * of MPI. It give
    also the possibility to get informations about that request :
    target, send/recv, tag, [a]synchronous, type, outcount.

. That identifier is used to control an asynchronous request
    via MPI_Access : Wait, Test, Probe, etc...

. In practise "RequestId" is simply an integer fo the interval
    [0 , 2**32-1]. There is only one such a cyclic for
    [I]Sends and [I]Recvs.

. That "RequestIds" and their associated structures give an easy
    way to manage asynchronous communications.
    For example we have mpi_access->Wait( int RequestId ) instead of
    MPI_Wait(MPI_Request *request, MPI_Status *status).

. The API of MPI_Access may give the "SendRequestIds" of a "target",
    the "RecvRequestIds" from a "source" or the "SendRequestIds" of
    all "targets" or the "RecvRequestIds" of all "sources".
  That avoid to manage them in Presentation-ParaMEDMEM.
*/

int MPI_Access::NewRequest( MPI_Datatype datatype, int tag , int destsourcerank ,
                            bool fromsourcerank , bool asynchronous ) {
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
  if ( !asynchronous ) {
    mpiaccessstruct->MPIRequest = MPI_REQUEST_NULL ;
    mpiaccessstruct->MPIStatus->MPI_SOURCE = destsourcerank ;
    mpiaccessstruct->MPIStatus->MPI_TAG = tag ;
    mpiaccessstruct->MPIStatus->MPI_ERROR = MPI_SUCCESS ;
  }
  //cout << "NewRequest old _Request" << _MyRank << " " << _Request
  // << " _MaxRequest " << _MaxRequest << endl ;
  if ( _Request == _MaxRequest ) {
    _Request = _BaseRequest ;
  }
  _Request += 1 ;
  _MapOfRequestStruct[_Request] = mpiaccessstruct ;
  if ( fromsourcerank ) {
    _RecvRequest[ destsourcerank ] = _Request ;
  }
  else {
    _SendRequest[ destsourcerank ] = _Request ;
  }
  if ( _Trace )
    cout << "NewRequest" << _MyRank << "( " << _Request << " ) "
         << mpiaccessstruct << endl ;
  return _Request ;
}

/*
MPI_Access and "tags" (or "MPITags") :
=====================================

. The constructor give the possibility to choose an interval of
    tags to use : [BaseTag , MaxTag].
  The default is [ 0 , MPI_TAG_UB], MPI_TAG_UB being the maximum
    value in an implementation of MPI (minimum 32767 = 2**15-1).
  On awa with the implementation lam MPI_TAG_UB value is
    7353944. The norma MPI specify that value is the same in all
    processes started by mpirun.
  In the case of the use of the same IntraCommunicator in a process
    for several distinct data flows (or for several IntraCommunicators
    with common processes), that permits to avoid ambibuity
    and may help debug.

. In MPI_Access the tags have two parts (#define ModuloTag 10) :
  + The last decimal digit decimal correspond to MPI_DataType ( 1 for
    TimeMessages, 2 for MPI_INT and 3 for MPI_DOUBLE)
  + The value of other digits correspond to a circular numero for each
    message.
  + A TimeMessage and the associated DataMessage have the same numero
    (but the types are different and the tags also).

. For a Send of a message from a process "source" to a process
    "target", we have _SendMPITag[target] in the process
    source (it contains the last "tag" used for the Send of a pour l'envoi de
    message to the process target).
  And in the process "target" which receive that message, we have
    _RecvMPITag[source] (it contains the last "tag" used for the Recv
    of messages from the process source).
  Naturally in the MPI norma the values of that tags must be the same.
*/
int MPI_Access::NewSendTag( MPI_Datatype datatype, int destrank , int method ,
                            bool asynchronous, int &RequestId ) {
  int tag ;
//  if ( method == _MessageTime ) {
//    tag = method ;
//  }
//  else {
    tag = IncrTag( _SendMPITag[destrank] ) ;
    tag = ValTag( tag, method ) ;
    _SendMPITag[ destrank ] = tag ;
//  }
  RequestId = NewRequest( datatype, tag, destrank , false , asynchronous ) ;
  _SendRequest[ destrank ] = RequestId ;
  _SendRequests[ destrank ].push_back( RequestId ) ;
  return tag ;
}

int MPI_Access::NewRecvTag( MPI_Datatype datatype, int sourcerank , int method ,
                            bool asynchronous, int &RequestId ) {
  int tag ;
//  if ( method == _MessageTime ) {
//    tag = method ;
//  }
//  else {
    tag = IncrTag( _RecvMPITag[sourcerank] ) ;
    tag = ValTag( tag, method ) ;
    _RecvMPITag[ sourcerank ] = tag ;
//  }
  RequestId = NewRequest( datatype, tag , sourcerank , true , asynchronous ) ;
  _RecvRequest[ sourcerank ] = RequestId ;
  _RecvRequests[ sourcerank ].push_back( RequestId ) ;
  return tag ;
}

// Returns the number of all SendRequestIds that may be used to allocate
// ArrayOfSendRequests for the call to SendRequestIds
int MPI_Access::SendRequestIdsSize() {
  int size = 0 ;
  int i ;
  for (i = 0 ; i < _ProcessorGroupSize ; i++ ) {
     size += _SendRequests[ i ].size() ;
  }
  return size ;
}

// Returns in ArrayOfSendRequests with the dimension "size" all the
// SendRequestIds
int MPI_Access::SendRequestIds(int size, int *ArrayOfSendRequests) {
  int destrank ;
  int i = 0 ;
  for ( destrank = 0 ; destrank < _ProcessorGroupSize ; destrank++ ) {
     list< int >::const_iterator iter ;
     for (iter = _SendRequests[ destrank ].begin() ; iter != _SendRequests[destrank].end() ; iter++ ) {
        ArrayOfSendRequests[i++] = *iter ;
     }
  }
  return i ;
}

// Returns the number of all RecvRequestIds that may be used to allocate
// ArrayOfRecvRequests for the call to RecvRequestIds
int MPI_Access::RecvRequestIdsSize() {
  int size = 0 ;
  int i ;
  for (i = 0 ; i < _ProcessorGroupSize ; i++ ) {
     size += _RecvRequests[ i ].size() ;
  }
  return size ;
}

// Returns in ArrayOfRecvRequests with the dimension "size" all the
// RecvRequestIds
int MPI_Access::RecvRequestIds(int size, int *ArrayOfRecvRequests) {
  int sourcerank ;
  int i = 0 ;
  for ( sourcerank = 0 ; sourcerank < _ProcessorGroupSize ; sourcerank++ ) {
     list< int >::const_iterator iter ;
     for (iter = _RecvRequests[ sourcerank ].begin() ; iter != _RecvRequests[sourcerank].end() ; iter++ ) {
        ArrayOfRecvRequests[i++] = *iter ;
     }
  }
  return i ;
}

// Returns in ArrayOfSendRequests with the dimension "size" all the
// SendRequestIds to a destination rank
int MPI_Access::SendRequestIds(int destrank, int size, int *ArrayOfSendRequests) {
  if (size < _SendRequests[destrank].size() ) throw MEDMEM::MEDEXCEPTION("wrong call to MPI_Access::SendRequestIds");
  int i = 0 ;
  list< int >::const_iterator iter ;
  for (iter = _SendRequests[ destrank ].begin() ; iter != _SendRequests[destrank].end() ; iter++ ) {
     ArrayOfSendRequests[i++] = *iter ;
  }
  return _SendRequests[destrank].size() ;
}

// Returns in ArrayOfRecvRequests with the dimension "size" all the
// RecvRequestIds from a sourcerank
int MPI_Access::RecvRequestIds(int sourcerank, int size, int *ArrayOfRecvRequests) {
  if (size < _RecvRequests[sourcerank].size() ) throw MEDMEM::MEDEXCEPTION("wrong call to MPI_Access::RecvRequestIds");
  int i = 0 ;
  list< int >::const_iterator iter ;
  _RecvRequests[ sourcerank ] ;
  for (iter = _RecvRequests[ sourcerank ].begin() ; iter != _RecvRequests[sourcerank].end() ; iter++ ) {
     ArrayOfRecvRequests[i++] = *iter ;
  }
  return _RecvRequests[sourcerank].size() ;
}

// Send in synchronous mode count values of type datatype from buffer to target
// (returns RequestId identifier even if the corresponding structure is deleted :
// it is only in order to have the same signature as the asynchronous mode)
int MPI_Access::Send(void* buffer, int count, MPI_Datatype datatype, int target,
                     int &RequestId) {
  int sts = MPI_SUCCESS ;
  RequestId = -1 ;
  if ( count ) {
    _MessageIdent aMethodIdent = MethodId( datatype ) ;
    int MPItag = NewSendTag( datatype, target , aMethodIdent , false , RequestId ) ;
    if ( aMethodIdent == _MessageTime ) {
      TimeMessage *aTimeMsg = (TimeMessage *) buffer ;
      aTimeMsg->tag = MPItag ;
    }
    DeleteRequest( RequestId ) ;
    sts = _CommInterface.send(buffer, count, datatype, target, MPItag,
                              *_IntraCommunicator ) ;
    if ( _Trace )
      cout << "MPI_Access::Send" << _MyRank << " SendRequestId "
           << RequestId << " count " << count << " target " << target
           << " MPItag " << MPItag << endl ;
  }
  return sts ;
}

// Receive (read) in synchronous mode count values of type datatype in buffer from source
// (returns RequestId identifier even if the corresponding structure is deleted :
// it is only in order to have the same signature as the asynchronous mode)
// The output argument OutCount is optionnal : *OutCount <= count
int MPI_Access::Recv(void* buffer, int count, MPI_Datatype datatype, int source,
                     int &RequestId, int *OutCount) {
  int sts = MPI_SUCCESS ;
  RequestId = -1 ;
  if ( OutCount != NULL ) {
    *OutCount = -1 ;
  }
  if ( count ) {
    _MessageIdent aMethodIdent = MethodId( datatype ) ;
    int MPItag = NewRecvTag( datatype, source , aMethodIdent , false , RequestId ) ;
    sts =  _CommInterface.recv(buffer, count, datatype, source, MPItag,
                               *_IntraCommunicator , MPIStatus( RequestId ) ) ;
    int outcount = 0 ;
    if ( sts == MPI_SUCCESS ) {
      MPI_Datatype datatype = MPIDatatype( RequestId ) ;
      _CommInterface.get_count(MPIStatus( RequestId ), datatype, &outcount ) ;
      SetMPIOutCount( RequestId , outcount ) ;
      SetMPICompleted( RequestId , true ) ;
      DeleteStatus( RequestId ) ;
    }
    if ( OutCount != NULL ) {
      *OutCount = outcount ;
    }
    if ( _Trace )
      cout << "MPI_Access::Recv" << _MyRank << " RecvRequestId "
           << RequestId << " count " << count << " source " << source
           << " MPItag " << MPItag << endl ;
    DeleteRequest( RequestId ) ;
  }
  return sts ;
}

// Send in asynchronous mode count values of type datatype from buffer to target
// Returns RequestId identifier.
int MPI_Access::ISend(void* buffer, int count, MPI_Datatype datatype, int target,
                      int &RequestId) {
  int sts = MPI_SUCCESS ;
  RequestId = -1 ;
  if ( count ) {
    _MessageIdent aMethodIdent = MethodId( datatype ) ;
    int MPItag = NewSendTag( datatype, target , aMethodIdent , true , RequestId ) ;
    if ( aMethodIdent == _MessageTime ) {
      TimeMessage *aTimeMsg = (TimeMessage *) buffer ;
      aTimeMsg->tag = MPItag ;
    }
    MPI_Request *aSendRequest = MPIRequest( RequestId ) ;
    if ( _Trace ) {
      cout << "MPI_Access::ISend" << _MyRank << " ISendRequestId "
           << RequestId << " count " << count << " target " << target
           << " MPItag " << MPItag << endl ;
      if ( MPItag == 1 )
        cout << "MPI_Access::ISend" << _MyRank << " time "
             << ((TimeMessage *)buffer)->time << " " << ((TimeMessage *)buffer)->deltatime
             << endl ;
    }
    sts = _CommInterface.Isend(buffer, count, datatype, target, MPItag,
                                *_IntraCommunicator , aSendRequest) ;
  }
  return sts ;
}

// Receive (read) in asynchronous mode count values of type datatype in buffer from source
// returns RequestId identifier.
int MPI_Access::IRecv(void* buffer, int count, MPI_Datatype datatype, int source,
                      int &RequestId) {
  int sts = MPI_SUCCESS ;
  RequestId = -1 ;
  if ( count ) {
    _MessageIdent aMethodIdent = MethodId( datatype ) ;
    int MPItag = NewRecvTag( datatype, source , aMethodIdent , true , RequestId ) ;
    MPI_Request *aRecvRequest = MPIRequest( RequestId ) ;
    if ( _Trace ) {
      cout << "MPI_Access::IRecv" << _MyRank << " IRecvRequestId "
           << RequestId << " count " << count << " source " << source
           << " MPItag " << MPItag << endl ;
      if ( MPItag == 1 )
        cout << "MPI_Access::ISend" << _MyRank << " time "
             << ((TimeMessage *)buffer)->time << " " << ((TimeMessage *)buffer)->deltatime
             << endl ;
    }
    sts = _CommInterface.Irecv(buffer, count, datatype, source, MPItag,
                                *_IntraCommunicator , aRecvRequest) ;
  }
  return sts ;
}

// Perform a Send and a Recv in synchronous mode
int MPI_Access::SendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                         int dest, int &SendRequestId,
                         void* recvbuf, int recvcount, MPI_Datatype recvtype,
                         int source, int &RecvRequestId, int *OutCount) {
  int sts = MPI_SUCCESS ;
  SendRequestId = -1 ;
  RecvRequestId = -1 ;
  if ( recvcount ) {
    sts = IRecv(recvbuf, recvcount, recvtype, source, RecvRequestId) ;
  }
  int outcount = -1 ;
  if ( _Trace )
    cout << "MPI_Access::SendRecv" << _MyRank << " IRecv RecvRequestId "
         << RecvRequestId << endl ;
  if ( sts == MPI_SUCCESS ) {
    if ( sendcount ) {
      sts = Send(sendbuf, sendcount, sendtype, dest, SendRequestId) ;
    }
    if ( _Trace )
      cout << "MPI_Access::SendRecv" << _MyRank << " Send SendRequestId "
           << SendRequestId << endl ;
    if ( sts == MPI_SUCCESS && recvcount ) {
      sts = Wait( RecvRequestId ) ;
      outcount = MPIOutCount( RecvRequestId ) ;
      if ( _Trace )
        cout << "MPI_Access::SendRecv" << _MyRank << " IRecv RecvRequestId "
             << RecvRequestId << " outcount " << outcount << endl ;
    }
  }
  if ( OutCount != NULL ) {
    *OutCount = outcount ;
    if ( _Trace )
      cout << "MPI_Access::SendRecv" << _MyRank << " *OutCount = " << *OutCount
           << endl ;
  }
  DeleteRequest( RecvRequestId ) ;
  return sts ;
}

// Perform a Send and a Recv in asynchronous mode
int MPI_Access::ISendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                          int dest, int &SendRequestId,
                          void* recvbuf, int recvcount, MPI_Datatype recvtype,
                          int source, int &RecvRequestId) {
  int sts = MPI_SUCCESS ;
  SendRequestId = -1 ;
  RecvRequestId = -1 ;
  if ( recvcount ) {
    sts = IRecv(recvbuf, recvcount, recvtype, source, RecvRequestId) ;
  }
  if ( sts == MPI_SUCCESS ) {
    if ( sendcount ) {
      sts = ISend(sendbuf, sendcount, sendtype, dest, SendRequestId) ;
    }
  }
  return sts ;
}

// Perform a wait of a Send or Recv asynchronous Request
// Do nothing for a synchronous Request
// Manage MPI_Request * and MPI_Status * structure
int MPI_Access::Wait( int RequestId ) {
  int status = MPI_SUCCESS ;
  if ( !MPICompleted( RequestId ) ) {
    if ( *MPIRequest( RequestId ) != MPI_REQUEST_NULL ) {
      if ( _Trace )
        cout << "MPI_Access::Wait" << _MyRank << " -> wait( " << RequestId
             << " ) MPIRequest " << MPIRequest( RequestId ) << " MPIStatus "
             << MPIStatus( RequestId ) << " MPITag " << MPITag( RequestId )
             << " MPIIsRecv " << MPIIsRecv( RequestId ) << endl ;
      status = _CommInterface.wait(MPIRequest( RequestId ), MPIStatus( RequestId )) ;
    }
    else {
      if ( _Trace )
        cout << "MPI_Access::Wait" << _MyRank << " MPIRequest == MPI_REQUEST_NULL"
             << endl ;
    }
    SetMPICompleted( RequestId , true ) ;
    if ( MPIIsRecv( RequestId ) && MPIStatus( RequestId ) ) {
      MPI_Datatype datatype = MPIDatatype( RequestId ) ;
      int outcount ;
      status = _CommInterface.get_count(MPIStatus( RequestId ), datatype,
                                        &outcount ) ;
      if ( status == MPI_SUCCESS ) {
        SetMPIOutCount( RequestId , outcount ) ;
        DeleteStatus( RequestId ) ;
        if ( _Trace )
          cout << "MPI_Access::Wait" << _MyRank << " RequestId " << RequestId
               << "MPIIsRecv " << MPIIsRecv( RequestId ) << " outcount " << outcount
               << endl ;
      }
      else {
        if ( _Trace )
          cout << "MPI_Access::Wait" << _MyRank << " MPIIsRecv "
               << MPIIsRecv( RequestId ) << " outcount " << outcount << endl ;
      }
    }
    else {
      if ( _Trace )
        cout << "MPI_Access::Wait" << _MyRank << " MPIIsRecv " << MPIIsRecv( RequestId )
             << " MPIOutCount " << MPIOutCount( RequestId ) << endl ;
    }
  }
  if ( _Trace )
    cout << "MPI_Access::Wait" << _MyRank << " RequestId " << RequestId
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
// Manage MPI_Request * and MPI_Status * structure
int MPI_Access::Test(int RequestId, int &flag) {
  int status = MPI_SUCCESS ;
  flag = MPICompleted( RequestId ) ;
  if ( _Trace )
    cout << "MPI_Access::Test" << _MyRank << " flag " << flag ;
  if ( MPIIsRecv( RequestId ) ) {
    if ( _Trace )
      cout << " Recv" ;
  }
  else {
    if ( _Trace )
      cout << " Send" ;
  }
  if ( _Trace )
    cout << "Request" << RequestId << " " << MPIRequest( RequestId )
          << " Status " << MPIStatus( RequestId ) << endl ;
  if ( !flag ) {
    if ( *MPIRequest( RequestId ) != MPI_REQUEST_NULL ) {
      if ( _Trace )
        cout << "MPI_Access::Test" << _MyRank << " -> test( " << RequestId
             << " ) MPIRequest " << MPIRequest( RequestId )
             << " MPIStatus " << MPIStatus( RequestId )
             << " MPITag " << MPITag( RequestId )
             << " MPIIsRecv " << MPIIsRecv( RequestId ) << endl ;
      status = _CommInterface.test(MPIRequest( RequestId ), &flag,
                                   MPIStatus( RequestId )) ;
    }
    else {
      if ( _Trace )
        cout << "MPI_Access::Test" << _MyRank << " MPIRequest == MPI_REQUEST_NULL"
             << endl ;
    }
    if ( flag ) {
      SetMPICompleted( RequestId , true ) ;
      if ( MPIIsRecv( RequestId ) && MPIStatus( RequestId ) ) {
        int outcount ;
        MPI_Datatype datatype = MPIDatatype( RequestId ) ;
        status = _CommInterface.get_count( MPIStatus( RequestId ), datatype,
                                           &outcount ) ;
        if ( status == MPI_SUCCESS ) {
          SetMPIOutCount( RequestId , outcount ) ;
          DeleteStatus( RequestId ) ;
          if ( _Trace )
            cout << "MPI_Access::Test" << _MyRank << " MPIIsRecv "
                 << MPIIsRecv( RequestId ) << " outcount " << outcount << endl ;
        }
        else {
          if ( _Trace )
            cout << "MPI_Access::Test" << _MyRank << " MPIIsRecv "
                 << MPIIsRecv( RequestId ) << " outcount " << outcount << endl ;
        }
      }
      else {
        if ( _Trace )
          cout << "MPI_Access::Test" << _MyRank << " MPIIsRecv "
               << MPIIsRecv( RequestId ) << " MPIOutCount "
               << MPIOutCount( RequestId ) << endl ;
      }
    }
  }
  if ( _Trace )
    cout << "MPI_Access::Test" << _MyRank << " RequestId " << RequestId
         << " flag " << flag << " MPICompleted " << MPICompleted( RequestId )
         << " MPIOutCount " << MPIOutCount( RequestId ) << endl ;
  return status ;
}

int MPI_Access::WaitAny(int count, int *array_of_RequestIds, int &RequestId) {
  int status = MPI_ERR_OTHER ;
  RequestId = -1 ;
  cout << "MPI_Access::WaitAny not yet implemented" << endl ;
  return status ;
}

int MPI_Access::TestAny(int count, int *array_of_RequestIds, int &RequestId,
                        int &flag) {
  int status = MPI_ERR_OTHER ;
  RequestId = -1 ;
  flag = 0 ;
  cout << "MPI_Access::TestAny not yet implemented" << endl ;
  return status ;
}

// Perform a wait of each Send or Recv asynchronous Request of the array 
// array_of_RequestIds of size "count".
// That array may be filled with a call to SendRequestIdsSize or RecvRequestIdsSize
// Do nothing for a synchronous Request
// Manage MPI_Request * and MPI_Status * structure
int MPI_Access::WaitAll(int count, int *array_of_RequestIds) {
  if ( _Trace )
    cout << "WaitAll" << _MyRank << " : count " << count << endl ;
  int status ;
  int retstatus = MPI_SUCCESS ;
  int i ;
  for ( i = 0 ; i < count ; i++ ) {
     if ( _Trace )
       cout << "WaitAll" << _MyRank << " " << i << " -> Wait( "
            << array_of_RequestIds[i] << " )" << endl ;
     status = Wait( array_of_RequestIds[i] ) ;
     if ( status != MPI_SUCCESS ) {
       retstatus = status ;
     }
  }
  if ( _Trace )
    cout << "EndWaitAll" << _MyRank << endl ;
  return retstatus ;
}

// Perform a "test" of each Send or Recv asynchronous Request of the array 
// array_of_RequestIds of size "count".
// That array may be filled with a call to SendRequestIdsSize or RecvRequestIdsSize
// If all requests are done, returns true in the flag argument
// If all requests are not finished, returns false in the flag argument
// Do nothing for a synchronous Request
// Manage MPI_Request * and MPI_Status * structure
int MPI_Access::TestAll(int count, int *array_of_RequestIds, int &flag) {
  if ( _Trace )
    cout << "TestAll" << _MyRank << " : count " << count << endl ;
  int status ;
  int retstatus = MPI_SUCCESS ;
  bool retflag = true ;
  int i ;
  for ( i = 0 ; i < count ; i++ ) {
     status = Test( array_of_RequestIds[i] , flag ) ;
     retflag = retflag && (flag != 0) ;
     if ( status != MPI_SUCCESS ) {
       retstatus = status ;
     }
  }
  flag = retflag ;
  if ( _Trace )
    cout << "EndTestAll" << _MyRank << endl ;
  return retstatus ;
}

int MPI_Access::WaitSome(int count, int *array_of_RequestIds, int outcount,
                         int *outarray_of_RequestIds) {
  int status = MPI_ERR_OTHER ;
  cout << "MPI_Access::WaitSome not yet implemented" << endl ;
  return status ;
}

int MPI_Access::TestSome(int count, int *array_of_RequestIds, int outcounts,
                         int *outarray_of_RequestIds) {
  int status = MPI_ERR_OTHER ;
  cout << "MPI_Access::TestSome not yet implemented" << endl ;
  return status ;
}

// Probe checks if a message is available for read from FromSource rank.
// Returns the corresponding source, MPITag, datatype and outcount
// Probe is a blocking call which wait until a message is available
int MPI_Access::Probe(int FromSource, int &source, int &MPITag,
                      MPI_Datatype &datatype, int &outcount) {
  MPI_Status aMPIStatus ;
  int sts =  _CommInterface.probe( FromSource, MPI_ANY_TAG,
                                   *_IntraCommunicator , &aMPIStatus ) ;
  if ( sts == MPI_SUCCESS ) {
    source = aMPIStatus.MPI_SOURCE ;
    MPITag = aMPIStatus.MPI_TAG ;
    int MethodId = (MPITag % ModuloTag) ;
    datatype = Datatype( (ParaMEDMEM::_MessageIdent) MethodId ) ;
    _CommInterface.get_count(&aMPIStatus, datatype, &outcount ) ;
    if ( _Trace )
      cout << "MPI_Access::Probe" << _MyRank << " FromSource " << FromSource
           << " source " << source << " MPITag " << MPITag << " MethodId "
           << MethodId << " datatype " << datatype << " outcount " << outcount
           << endl ;
  }
  else {
    source = -1 ;
    MPITag = -1 ;
    datatype = 0 ;
    outcount = -1 ;
  }
  return sts ;
}

// IProbe checks if a message is available for read from FromSource rank.
// If there is a message available, returns the corresponding source,
// MPITag, datatype and outcount with flag = true
// If not, returns flag = false
int MPI_Access::IProbe(int FromSource, int &source, int &MPITag,
                       MPI_Datatype &datatype, int &outcount, int &flag) {
  MPI_Status aMPIStatus ;
  int sts =  _CommInterface.Iprobe( FromSource, MPI_ANY_TAG,
                                    *_IntraCommunicator , &flag,
                                    &aMPIStatus ) ;
  if ( sts == MPI_SUCCESS && flag ) {
    source = aMPIStatus.MPI_SOURCE ;
    MPITag = aMPIStatus.MPI_TAG ;
    int MethodId = (MPITag % ModuloTag) ;
    datatype = Datatype( (ParaMEDMEM::_MessageIdent) MethodId ) ;
    _CommInterface.get_count(&aMPIStatus, datatype, &outcount ) ;
    if ( _Trace )
      cout << "MPI_Access::IProbe" << _MyRank << " FromSource " << FromSource
           << " source " << source << " MPITag " << MPITag << " MethodId "
           << MethodId << " datatype " << datatype << " outcount " << outcount
           << " flag " << flag << endl ;
  }
  else {
    source = -1 ;
    MPITag = -1 ;
    datatype = 0 ;
    outcount = -1 ;
  }
  return sts ;
}

// Cancel concerns a "posted" asynchronous IRecv
// Returns flag = true if the receiving request was successfully canceled
// Returns flag = false if the receiving request was finished but not canceled
// Use cancel, wait and test_cancelled of the MPI API
int MPI_Access::Cancel( int RecvRequestId, int &flag ) {
  flag = 0 ;
  int sts = _CommInterface.cancel( MPIRequest( RecvRequestId ) ) ;
  if ( sts == MPI_SUCCESS ) {
    sts = _CommInterface.wait( MPIRequest( RecvRequestId ) ,
                               MPIStatus( RecvRequestId ) ) ;
    if ( sts == MPI_SUCCESS ) {
      sts = _CommInterface.test_cancelled( MPIStatus( RecvRequestId ) , &flag ) ;
    }
  }
  return sts ;
}

// Cancel concerns a "pending" receiving message (without IRecv "posted")
// Returns flag = true if the message was successfully canceled
// Returns flag = false if the receiving request was finished but not canceled
// Use Irecv, cancel, wait and test_cancelled of the MPI API
int MPI_Access::Cancel( int source, int theMPITag, MPI_Datatype datatype,
                        int outcount, int &flag ) {
  int sts ;
  MPI_Aint extent ;
  flag = 0 ;
  sts =  MPI_Type_extent( datatype , &extent ) ;
  if ( sts == MPI_SUCCESS ) {
    void * recvbuf = malloc( extent*outcount ) ;
    MPI_Request aRecvRequest ;
    if ( _Trace )
      cout << "MPI_Access::Cancel" << _MyRank << " Irecv extent " << extent
           << " datatype " << datatype << " source " << source << " theMPITag "
           << theMPITag << endl ;
    sts = _CommInterface.Irecv( recvbuf, outcount, datatype, source, theMPITag,
                                *_IntraCommunicator , &aRecvRequest ) ;
    if ( sts == MPI_SUCCESS ) {
      sts = _CommInterface.cancel( &aRecvRequest ) ;
      if ( _Trace )
        cout << "MPI_Access::Cancel" << _MyRank << " theMPITag " << theMPITag
             << " cancel done" << endl ;
      if ( sts == MPI_SUCCESS ) {
        MPI_Status aStatus ;
        if ( _Trace )
          cout << "MPI_Access::Cancel" << _MyRank << " wait" << endl ;
        sts = _CommInterface.wait( &aRecvRequest , &aStatus ) ;
        if ( sts == MPI_SUCCESS ) {
          if ( _Trace )
            cout << "MPI_Access::Cancel" << _MyRank << " test_cancelled" << endl ;
          sts = _CommInterface.test_cancelled( &aStatus , &flag ) ;
        }
      }
    }
    if ( _Trace && datatype == TimeType() )
      cout << "MPI_Access::Cancel" << _MyRank << " time "
           << ((TimeMessage *) recvbuf)->time << " "
           << ((TimeMessage *) recvbuf)->deltatime << endl ;
    free( recvbuf ) ;
  }
  if ( _Trace )
    cout << "MPI_Access::Cancel" << _MyRank << " flag " << flag << endl ;
  return sts ;
}


// CancelAll concerns all "pending" receiving message (without IRecv "posted")
// CancelAll use IProbe and Cancel (see obove)
int MPI_Access::CancelAll() {
  int sts = MPI_SUCCESS ;
  int target ;
  int source ;
  int MPITag ;
  MPI_Datatype datatype ;
  int outcount ;
  int flag ;
  for ( target = 0 ; target < _ProcessorGroupSize ; target++ ) {
     sts = IProbe(target, source, MPITag, datatype, outcount, flag) ;
     if ( sts == MPI_SUCCESS && flag ) {
       sts = Cancel(source, MPITag, datatype, outcount, flag) ;
       if ( _Trace )
         cout << "MPI_Access::CancelAll" << _MyRank << " source " << source
              << " MPITag " << MPITag << " datatype " << datatype
              << " outcount " << outcount << " Cancel flag " << flag << endl ;
       if ( sts != MPI_SUCCESS ) {
         break ;
       }
     }
     else if ( sts != MPI_SUCCESS ) {
       break ;
     }
  }
  return sts ;
}

// Same as barrier of MPI API
int MPI_Access::Barrier() {
  int status = _CommInterface.barrier( *_IntraCommunicator ) ;
  return status ;
}

// Same as Error_String of MPI API
int MPI_Access::Error_String(int errorcode, char *string, int *resultlen) const {
  return _CommInterface.error_string( errorcode, string, resultlen) ;
}

// Returns source, tag, error and outcount corresponding to receiving RequestId
// By default the corresponding structure of RequestId is deleted
int MPI_Access::Status(int RequestId, int &source, int &tag, int &error,
                       int &outcount, bool keepRequestStruct) {
  MPI_Status *status = MPIStatus( RequestId ) ;
  if ( _Trace )
    cout << "MPI_Access::Status" << _MyRank << " RequestId " << RequestId
         << " status " << status << endl ;
  if ( status != NULL && MPIAsynchronous( RequestId ) &&
       MPICompleted( RequestId ) ) {
    if ( MPIIsRecv( RequestId ) ) {
      source = status->MPI_SOURCE ;
      tag = status->MPI_TAG ;
      error = status->MPI_ERROR ;
      MPI_Datatype datatype = MPIDatatype( RequestId ) ;
      _CommInterface.get_count(status, datatype, &outcount ) ;
      if ( _Trace )
        cout << "MPI_Access::Status" << _MyRank << " RequestId " << RequestId
             << " status " << status << " outcount " << outcount << endl ;
      SetMPIOutCount( RequestId , outcount ) ;
    }
    else {
      source = MPITarget( RequestId ) ;
      tag = MPITag( RequestId ) ;
      error = 0 ;
      outcount = MPIOutCount( RequestId ) ;
    }
    if ( !keepRequestStruct ) {
      DeleteRequest( RequestId ) ;
    }
    return MPI_SUCCESS ;
  }
  else {
    source = MPITarget( RequestId ) ;
    tag = MPITag( RequestId ) ;
    error = 0 ;
    outcount = MPIOutCount( RequestId ) ;
  }
  return MPI_SUCCESS ;
}

int MPI_Access::Request_Free( MPI_Request *request ) {
  return _CommInterface.request_free( request ) ;
}

// Print all informations of all known requests for debugging purpose
void MPI_Access::Check() const {
  int i = 0 ;
  map< int , RequestStruct * >::const_iterator MapOfRequestStructiterator ;
  cout << "MPI_Access::Check" << _MyRank << "_MapOfRequestStructSize "
       << _MapOfRequestStruct.size() << endl ;
  for ( MapOfRequestStructiterator = _MapOfRequestStruct.begin() ;
        MapOfRequestStructiterator != _MapOfRequestStruct.end() ;
        MapOfRequestStructiterator++ ) {
     if ( MapOfRequestStructiterator->second == NULL ) {
//       cout << " MapOfRequestStructiterator->second == NULL" << endl ;
     }
     else {
       cout << "    Check" << _MyRank << " " << i << ". Request"
            << MapOfRequestStructiterator->first << "-->" ;
       if ( (MapOfRequestStructiterator->second)->MPIAsynchronous ) {
         cout << "I" ;
       }
       if ( (MapOfRequestStructiterator->second)->MPIIsRecv ) {
         cout << "Recv from " ;
       }
       else {
         cout << "Send to " ;
       }
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
  cout << "EndCheck" << _MyRank << endl ;
}

// Outputs fields of a TimeMessage structure
ostream & operator<< (ostream & f ,const TimeMessage & aTimeMsg ) {
  f << " time " << aTimeMsg.time << " deltatime " << aTimeMsg.deltatime
    << " tag " << aTimeMsg.tag ;

  return f;
}

// Outputs the DataType coded in a Tag
ostream & operator<< (ostream & f ,const _MessageIdent & methodtype ) {
  switch (methodtype) {
  case _MessageTime :
    f << " MethodTime ";
    break;
  case _MessageInt :
    f << " MPI_INT ";
    break;
  case _MessageDouble :
    f << " MPI_DOUBLE ";
    break;
  default :
    f << " UnknownMethodType ";
    break;
  }

  return f;
}

}

