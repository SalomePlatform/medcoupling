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

#ifndef __MPIACCESS_HXX__
#define __MPIACCESS_HXX__

#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"

#include <map>
#include <list>
#include <vector>
#include <iostream>

namespace MEDCoupling
{
  typedef struct
  {
    double time ;
    double deltatime ;
    int tag ;
  } TimeMessage;
  
  static MPI_Request mpirequestnull = MPI_REQUEST_NULL ;
  enum _MessageIdent { _message_unknown, _message_time, _message_int, _message_double } ;

  class MPIAccess
  {
  private:
    struct RequestStruct
    {
      int MPITarget ;
      bool MPIIsRecv ;
      int MPITag ;
      bool MPIAsynchronous ;
      bool MPICompleted ;
      MPI_Datatype MPIDatatype ;
      MPI_Request MPIRequest ;
      MPI_Status *MPIStatus ;
      int MPIOutCount ;
    };
  public:
    MPIAccess(MPIProcessorGroup * ProcessorGroup, int BaseTag=0, int MaxTag=0) ;
    virtual ~MPIAccess() ;

    void trace( bool trace = true ) ;

    void deleteRequest( int RequestId ) ;
    void deleteRequests(int size , int *ArrayOfSendRequests ) ;

    int sendMPITag(int destrank) ;
    int recvMPITag(int sourcerank) ;

    int sendRequestIdsSize() ;
    int sendRequestIds(int size, int *ArrayOfSendRequests) ;
    int recvRequestIdsSize() ;
    int recvRequestIds(int size, int *ArrayOfRecvRequests) ;

    int sendRequestIdsSize(int destrank) ;
    int sendRequestIds(int destrank, int size, int *ArrayOfSendRequests) ;
    int recvRequestIdsSize(int sourcerank) ;
    int recvRequestIds(int sourcerank, int size, int *ArrayOfRecvRequests) ;

    int send(void* buffer, int count, MPI_Datatype datatype, int target,
             int &RequestId) ;
    int ISend(void* buffer, int count, MPI_Datatype datatype, int target,
              int &RequestId) ;
    int recv(void* buffer, int count, MPI_Datatype datatype, int source,
             int &RequestId, int *OutCount=NULL) ;
    int IRecv(void* buffer, int count, MPI_Datatype datatype, int source,
              int &RequestId) ;
    int sendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
                 int &SendRequestId, void* recvbuf, int recvcount,
                 MPI_Datatype recvtype, int source,
                 int &RecvRequestId, int *OutCount=NULL) ;
    int ISendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
                  int &SendRequestId, void* recvbuf, int recvcount,
                  MPI_Datatype recvtype, int source, int &RecvRequestId) ;

    int wait(int RequestId) ;
    int test(int RequestId, int &flag) ;
    int waitAny(int count, int *array_of_RequestIds, int &RequestId) ;
    int testAny(int count, int *array_of_RequestIds, int &RequestId, int &flag)  ;
    int waitAll(int count, int *array_of_RequestIds) ;
    int testAll(int count, int *array_of_RequestIds, int &flag)  ;
    int waitSome(int count, int *array_of_RequestIds, int outcount,
                 int *outarray_of_RequestIds) ;
    int testSome(int count, int *array_of_RequestIds, int outcounts,
                 int *outarray_of_RequestIds) ;
    int probe(int FromSource, int &source, int &MPITag, MPI_Datatype &datatype,
              int &outcount) ;
    int IProbe(int FromSource, int &source, int &MPITag, MPI_Datatype &datatype,
               int &outcount, int &flag) ;
    int cancel( int RecvRequestId, int &flag ) ; 
    int cancel( int source, int MPITag, MPI_Datatype datatype, int outcount,
                int &flag ) ;
    int cancelAll() ;
    int barrier() ;
    int errorString(int errorcode, char *string, int *resultlen) const ;
    int status(int RequestId, int &source, int &tag, int &error, int &outcount,
               bool keepRequestStruct=false) ;
    int requestFree( MPI_Request *request ) ;

    void check() const ;

    MPI_Datatype timeType() const ;
    bool isTimeMessage( int MPITag ) const ;
    MPI_Aint timeExtent() const ;
    MPI_Aint intExtent() const ;
    MPI_Aint doubleExtent() const ;
    MPI_Aint extent( MPI_Datatype datatype ) const ;

    int MPITag( int RequestId ) ;
    int MPITarget( int RequestId ) ;
    bool MPIIsRecv( int RequestId ) ;
    bool MPIAsynchronous( int RequestId ) ;
    bool MPICompleted( int RequestId ) ;
    MPI_Datatype MPIDatatype( int RequestId ) ;
    int MPIOutCount( int RequestId ) ;

  private:
    int newRequest( MPI_Datatype datatype, int tag , int destsourcerank ,
                    bool fromsourcerank , bool asynchronous ) ;
    int newSendTag( MPI_Datatype datatype, int destrank , int method ,
                    bool asynchronous, int &RequestId ) ;
    int newRecvTag( MPI_Datatype datatype, int sourcerank , int method ,
                    bool asynchronous, int &RequestId ) ;
    int incrTag( int prevtag ) ;
    int valTag( int tag, int method ) ;

    void deleteSendRecvRequest( int RequestId ) ;

    void deleteStatus( int RequestId ) ;

    MPI_Request *MPIRequest( int RequestId ) ;
    MPI_Status *MPIStatus( int RequestId ) ;
    void setMPICompleted( int RequestId , bool completed ) ;
    void setMPIOutCount( int RequestId , int outcount ) ;
    void clearMPIStatus( int RequestId ) ;

    _MessageIdent methodId( MPI_Datatype datatype ) const ;
    MPI_Datatype datatype( _MessageIdent aMethodIdent ) const ;
  private:
    const CommInterface &_comm_interface ;
    const MPI_Comm* _intra_communicator ;
    MPIProcessorGroup * _processor_group ;
    int _processor_group_size ;
    int _my_rank ;
    bool _trace ;
    int _base_request ;
    int _max_request ;
    int _request ;
    int * _send_request ;
    int * _recv_request ;
    std::vector< std::list< int > > _send_requests ;
    std::vector< std::list< int > > _recv_requests ;
    int _base_MPI_tag ;
    int _max_MPI_tag ;
    int * _send_MPI_tag ;
    int * _recv_MPI_Tag ;
    MPI_Datatype _MPI_TIME ;
    static const int MODULO_TAG=10;
    std::map< int , RequestStruct * > _map_of_request_struct ;

  };

  inline void MPIAccess::trace( bool atrace )
  {
    _trace = atrace ;
  }

  // Delete the structure Request corresponding to RequestId identifier after
  // the deletion of the structures MPI_Request * and MPI_Status *
  // remove it from _MapOfRequestStruct (erase)
  inline void MPIAccess::deleteRequest( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      {
        if ( _trace )
          std::cout << "MPIAccess::DeleteRequest" << _my_rank << "( " << RequestId << " ) "
                    << aRequestStruct << " MPIRequest " << aRequestStruct->MPIRequest
                    << " MPIIsRecv " << aRequestStruct->MPIIsRecv << std::endl ;
        if ( _map_of_request_struct[RequestId]->MPIRequest != MPI_REQUEST_NULL )
          requestFree( &_map_of_request_struct[RequestId]->MPIRequest ) ;
        deleteSendRecvRequest( RequestId ) ;
        deleteStatus( RequestId ) ;
        _map_of_request_struct.erase( RequestId ) ;
        delete aRequestStruct ;
      }
    else
      {
        if ( _trace )
          std::cout << "MPIAccess::DeleteRequest" << _my_rank << "( " << RequestId
                    << " ) Request not found" << std::endl ;
      }
  }

  // Delete all requests of the array ArrayOfSendRequests
  inline void MPIAccess::deleteRequests(int size , int *ArrayOfSendRequests )
  {
    for (int i = 0 ; i < size ; i++ )
      deleteRequest( ArrayOfSendRequests[i] ) ;
  }

  // Returns the last MPITag of the destination rank destrank
  inline int MPIAccess::sendMPITag(int destrank)
  {
    return _send_MPI_tag[destrank] ;
  }

  // Returns the last MPITag of the source rank sourcerank
  inline int MPIAccess::recvMPITag(int sourcerank)
  {
    return _recv_MPI_Tag[sourcerank] ;
  }

  // Returns the number of all SendRequestIds matching a destination rank. It may be
  // used to allocate ArrayOfSendRequests for the call to SendRequestIds
  inline int MPIAccess::sendRequestIdsSize(int destrank)
  {
    return _send_requests[destrank].size() ;
  }

  // Returns the number of all RecvRequestIds matching a source rank. It may be
  // used to allocate ArrayOfRecvRequests for the call to RecvRequestIds
  inline int MPIAccess::recvRequestIdsSize(int sourcerank)
  {
    return _recv_requests[sourcerank].size() ;
  }

  // Returns the MPI_Datatype (registered in MPI in the constructor with
  // MPI_Type_struct and MPI_Type_commit) for TimeMessages
  inline MPI_Datatype MPIAccess::timeType() const
  {
    return _MPI_TIME ;
  }
  
  // Returns true if the tag MPITag corresponds to a TimeMessage
  inline bool MPIAccess::isTimeMessage( int aMPITag ) const
  {
    return ((aMPITag%MODULO_TAG) == _message_time) ;
  }

  // Returns the MPI size of the MPI_Datatype datatype
  inline MPI_Aint MPIAccess::extent( MPI_Datatype adatatype ) const
  {
    if ( adatatype == _MPI_TIME )
      return timeExtent() ;
    if ( adatatype == MPI_INT )
      return intExtent() ;
    if ( adatatype == MPI_DOUBLE )
      return doubleExtent() ;
    return 0 ;
  }
  
  // Returns the MPITag of the request corresponding to RequestId identifier
  inline int MPIAccess::MPITag( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      return aRequestStruct->MPITag ;
    return -1 ;
  }
  
  // Returns the MPITarget of the request corresponding to RequestId identifier
  inline int MPIAccess::MPITarget( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      return aRequestStruct->MPITarget ;
    return -1 ;
  }

  // Returns true if the request corresponding to RequestId identifier was [I]Recv
  inline bool MPIAccess::MPIIsRecv( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      return aRequestStruct->MPIIsRecv ;
    return false ;
  }

  // Returns true if the request corresponding to RequestId identifier was asynchronous
  inline bool MPIAccess::MPIAsynchronous( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      return aRequestStruct->MPIAsynchronous ;
    return false ;
  }
  
  // Returns true if the request corresponding to RequestId identifier was completed
  inline bool MPIAccess::MPICompleted( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      return aRequestStruct->MPICompleted;
    return true ;
  }

  // Returns the MPI_datatype  of the request corresponding to RequestId identifier
  inline MPI_Datatype MPIAccess::MPIDatatype( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      return aRequestStruct->MPIDatatype;
    return MPI_DATATYPE_NULL;
  }

  // Returns the size of the receiving message of the request corresponding to
  // RequestId identifier
  inline int MPIAccess::MPIOutCount( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      return aRequestStruct->MPIOutCount;
    return 0 ;
  }

  // Increments the previous tag value (cyclically)
  // Look at MPIAccess::NewSendTag/NewRecvTag in MPIAccess.cxx
  inline int MPIAccess::incrTag( int prevtag )
  {
    int tag;
    if ( (prevtag % MODULO_TAG) == _message_time )
      tag = ((prevtag/MODULO_TAG)*MODULO_TAG);
    else
      tag = ((prevtag/MODULO_TAG + 1)*MODULO_TAG);
    if ( tag > _max_MPI_tag )
      tag = _base_MPI_tag ;
    return tag ;
  }

  // Returns the MPITag with the method-type field
  // Look at MPIAccess::NewSendTag/NewRecvTag in MPIAccess.cxx
  inline int MPIAccess::valTag( int tag, int method )
  {
    return ((tag/MODULO_TAG)*MODULO_TAG) + method;
  }
  
  // Remove a Request identifier from the list _RecvRequests/_SendRequests for
  // the corresponding target.
  inline void MPIAccess::deleteSendRecvRequest( int RequestId )
  {
    if ( _trace )
      std::cout << "MPIAccess::DeleteSendRecvRequest" << _my_rank
                << "( " << RequestId << " ) " << std::endl ;
    if ( MPIIsRecv( RequestId ) )
      _recv_requests[ MPITarget( RequestId ) ].remove( RequestId );
    else
      _send_requests[ MPITarget( RequestId ) ].remove( RequestId );
  }

  // Delete the MPI structure MPI_status * of a ReaquestId
  inline void MPIAccess::deleteStatus( int RequestId )
  {
    if ( _map_of_request_struct[RequestId]->MPIStatus != NULL )
      {
        delete _map_of_request_struct[RequestId]->MPIStatus ;
        clearMPIStatus( RequestId ) ;
      }
  }

  // Returns the MPI structure MPI_request * of a RequestId
  inline MPI_Request * MPIAccess::MPIRequest( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      return &aRequestStruct->MPIRequest;
    return &mpirequestnull ;
  }
  
  // Returns the MPI structure MPI_status * of a RequestId
  inline MPI_Status * MPIAccess::MPIStatus( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ];
    if ( aRequestStruct )
      return aRequestStruct->MPIStatus;
    return NULL ;
  }

  // Set the MPICompleted field of the structure Request corresponding to RequestId
  // identifier with the value completed
  inline void MPIAccess::setMPICompleted( int RequestId , bool completed )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      aRequestStruct->MPICompleted = completed;
  }

  // Set the MPIOutCount field of the structure Request corresponding to RequestId
  // identifier with the value outcount
  inline void MPIAccess::setMPIOutCount( int RequestId , int outcount )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      aRequestStruct->MPIOutCount = outcount;
  }

  // Nullify the MPIStatusfield of the structure Request corresponding to RequestId
  // identifier
  inline void MPIAccess::clearMPIStatus( int RequestId )
  {
    struct RequestStruct *aRequestStruct = _map_of_request_struct[ RequestId ] ;
    if ( aRequestStruct )
      aRequestStruct->MPIStatus = NULL ;
  }

  // Returns the _MessageIdent enum value corresponding to the MPI_Datatype datatype
  // Look at MPIAccess::NewSendTag/NewRecvTag in MPIAccess.cxx
  inline _MessageIdent MPIAccess::methodId( MPI_Datatype adatatype ) const
  {
    _MessageIdent aMethodIdent ;
    if ( adatatype == _MPI_TIME )
      aMethodIdent = _message_time;
    else if ( adatatype == MPI_INT )
      aMethodIdent = _message_int ;
    else if ( adatatype == MPI_DOUBLE )
      aMethodIdent = _message_double ;
    else
      aMethodIdent = _message_unknown ;
    return aMethodIdent ;
  }
  
  // Returns the MPI_Datatype corresponding to the _MessageIdent enum aMethodIdent
  inline MPI_Datatype MPIAccess::datatype( _MessageIdent aMethodIdent ) const
  {
    MPI_Datatype aDataType ;
    switch( aMethodIdent )
      {
      case _message_time :
        aDataType = _MPI_TIME ;
        break ;
      case _message_int :
        aDataType = MPI_INT ;
        break ;
      case _message_double :
        aDataType = MPI_DOUBLE ;
        break ;
      default :
        aDataType = (MPI_Datatype) -1 ;
        break ;
      }
    return aDataType ;
  }

  std::ostream & operator<< (std::ostream &,const _MessageIdent &);

  std::ostream & operator<< (std::ostream &,const TimeMessage &);

}

#endif
