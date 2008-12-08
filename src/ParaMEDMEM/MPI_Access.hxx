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
#ifndef MPI_ACCESS_HXX_
#define MPI_ACCESS_HXX_

#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"

#include <map>
#include <list>
#include <vector>
#include <iostream>

namespace ParaMEDMEM
{
typedef struct { double time ;
                 double deltatime ;
                 int tag ; } TimeMessage ;

static MPI_Request mpirequestnull = MPI_REQUEST_NULL ;
enum _MessageIdent { _MessageUnknown, _MessageTime, _MessageInt, _MessageDouble } ;

class MPI_Access {
  public:
#define ModuloTag 10
    MPI_Access(MPIProcessorGroup * ProcessorGroup, int BaseTag=0,
               int MaxTag=0) ;
    virtual ~MPI_Access() ;

    void Trace( bool trace = true ) ;

    void DeleteRequest( int RequestId ) ;
    void DeleteRequests(int size , int *ArrayOfSendRequests ) ;

    int SendMPITag(int destrank) ;
    int RecvMPITag(int sourcerank) ;

    int SendRequestIdsSize() ;
    int SendRequestIds(int size, int *ArrayOfSendRequests) ;
    int RecvRequestIdsSize() ;
    int RecvRequestIds(int size, int *ArrayOfRecvRequests) ;

    int SendRequestIdsSize(int destrank) ;
    int SendRequestIds(int destrank, int size, int *ArrayOfSendRequests) ;
    int RecvRequestIdsSize(int sourcerank) ;
    int RecvRequestIds(int sourcerank, int size, int *ArrayOfRecvRequests) ;

    int Send(void* buffer, int count, MPI_Datatype datatype, int target,
             int &RequestId) ;
    int ISend(void* buffer, int count, MPI_Datatype datatype, int target,
              int &RequestId) ;
    int Recv(void* buffer, int count, MPI_Datatype datatype, int source,
             int &RequestId, int *OutCount=NULL) ;
    int IRecv(void* buffer, int count, MPI_Datatype datatype, int source,
              int &RequestId) ;
    int SendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
                 int &SendRequestId, void* recvbuf, int recvcount,
                 MPI_Datatype recvtype, int source,
                 int &RecvRequestId, int *OutCount=NULL) ;
    int ISendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, int dest,
                  int &SendRequestId, void* recvbuf, int recvcount,
                  MPI_Datatype recvtype, int source, int &RecvRequestId) ;

    int Wait(int RequestId) ;
    int Test(int RequestId, int &flag) ;
    int WaitAny(int count, int *array_of_RequestIds, int &RequestId) ;
    int TestAny(int count, int *array_of_RequestIds, int &RequestId, int &flag)  ;
    int WaitAll(int count, int *array_of_RequestIds) ;
    int TestAll(int count, int *array_of_RequestIds, int &flag)  ;
    int WaitSome(int count, int *array_of_RequestIds, int outcount,
                 int *outarray_of_RequestIds) ;
    int TestSome(int count, int *array_of_RequestIds, int outcounts,
                 int *outarray_of_RequestIds) ;
    int Probe(int FromSource, int &source, int &MPITag, MPI_Datatype &datatype,
              int &outcount) ;
    int IProbe(int FromSource, int &source, int &MPITag, MPI_Datatype &datatype,
               int &outcount, int &flag) ;
    int Cancel( int RecvRequestId, int &flag ) ; 
    int Cancel( int source, int MPITag, MPI_Datatype datatype, int outcount,
                int &flag ) ;
    int CancelAll() ;
    int Barrier() ;
    int Error_String(int errorcode, char *string, int *resultlen) const ;
    int Status(int RequestId, int &source, int &tag, int &error, int &outcount,
               bool keepRequestStruct=false) ;
    int Request_Free( MPI_Request *request ) ;

    void Check() const ;

    MPI_Datatype TimeType() const ;
    bool IsTimeMessage( int MPITag ) const ;
    MPI_Aint TimeExtent() const ;
    MPI_Aint IntExtent() const ;
    MPI_Aint DoubleExtent() const ;
    MPI_Aint Extent( MPI_Datatype datatype ) const ;

    int MPITag( int RequestId ) ;
    int MPITarget( int RequestId ) ;
    bool MPIIsRecv( int RequestId ) ;
    bool MPIAsynchronous( int RequestId ) ;
    bool MPICompleted( int RequestId ) ;
    MPI_Datatype MPIDatatype( int RequestId ) ;
    int MPIOutCount( int RequestId ) ;

  private:
    int NewRequest( MPI_Datatype datatype, int tag , int destsourcerank ,
                    bool fromsourcerank , bool asynchronous ) ;
    int NewSendTag( MPI_Datatype datatype, int destrank , int method ,
                    bool asynchronous, int &RequestId ) ;
    int NewRecvTag( MPI_Datatype datatype, int sourcerank , int method ,
                    bool asynchronous, int &RequestId ) ;
    int IncrTag( int prevtag ) ;
    int ValTag( int tag, int method ) ;

    void DeleteSendRecvRequest( int RequestId ) ;

    void DeleteStatus( int RequestId ) ;

    MPI_Request *MPIRequest( int RequestId ) ;
    MPI_Status *MPIStatus( int RequestId ) ;
    void SetMPICompleted( int RequestId , bool completed ) ;
    void SetMPIOutCount( int RequestId , int outcount ) ;
    void ClearMPIStatus( int RequestId ) ;

    _MessageIdent MethodId( MPI_Datatype datatype ) const ;
    MPI_Datatype Datatype( _MessageIdent aMethodIdent ) const ;

    const CommInterface &_CommInterface ;
    const MPI_Comm* _IntraCommunicator ;
    MPIProcessorGroup * _ProcessorGroup ;
    int _ProcessorGroupSize ;
    int _MyRank ;
    bool _Trace ;

    int _BaseRequest ;
    int _MaxRequest ;
    int _Request ;
    int * _SendRequest ;
    int * _RecvRequest ;
    vector< list< int > > _SendRequests ;
    vector< list< int > > _RecvRequests ;

    int _BaseMPITag ;
    int _MaxMPITag ;
    int * _SendMPITag ;
    int * _RecvMPITag ;

    MPI_Datatype _MPI_TIME ;

    struct RequestStruct { int MPITarget ;
                           bool MPIIsRecv ;
                           int MPITag ;
                           bool MPIAsynchronous ;
                           bool MPICompleted ;
                           MPI_Datatype MPIDatatype ;
                           MPI_Request MPIRequest ;
                           MPI_Status *MPIStatus ;
                           int MPIOutCount ; } ;
    map< int , RequestStruct * > _MapOfRequestStruct ;

};

  inline void MPI_Access::Trace( bool trace ) {
         _Trace = trace ;
  }

// Delete the structure Request corresponding to RequestId identifier after
// the deletion of the structures MPI_Request * and MPI_Status *
// remove it from _MapOfRequestStruct (erase)
  inline void MPI_Access::DeleteRequest( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           if ( _Trace )
             std::cout << "MPI_Access::DeleteRequest" << _MyRank << "( " << RequestId << " ) "
                  << aRequestStruct << " MPIRequest " << aRequestStruct->MPIRequest
                       << " MPIIsRecv " << aRequestStruct->MPIIsRecv << std::endl ;
           if ( _MapOfRequestStruct[RequestId]->MPIRequest != MPI_REQUEST_NULL ) {
             Request_Free( &_MapOfRequestStruct[RequestId]->MPIRequest ) ;
           }
           DeleteSendRecvRequest( RequestId ) ;
           DeleteStatus( RequestId ) ;
           _MapOfRequestStruct.erase( RequestId ) ;
           delete aRequestStruct ;
         }
         else {
           if ( _Trace )
             std::cout << "MPI_Access::DeleteRequest" << _MyRank << "( " << RequestId
                       << " ) Request not found" << std::endl ;
         }
  }
// Delete all requests of the array ArrayOfSendRequests
  inline void MPI_Access::DeleteRequests(int size , int *ArrayOfSendRequests ) {
         int i ;
         for ( i = 0 ; i < size ; i++ ) {
            DeleteRequest( ArrayOfSendRequests[i] ) ;
         }
  }

// Returns the last MPITag of the destination rank destrank
  inline int MPI_Access::SendMPITag(int destrank) {
         return _SendMPITag[destrank] ;
  }
// Returns the last MPITag of the source rank sourcerank
  inline int MPI_Access::RecvMPITag(int sourcerank) {
         return _RecvMPITag[sourcerank] ;
  }

// Returns the number of all SendRequestIds matching a destination rank. It may be
// used to allocate ArrayOfSendRequests for the call to SendRequestIds
  inline int MPI_Access::SendRequestIdsSize(int destrank) {
         return _SendRequests[destrank].size() ;
  }
// Returns the number of all RecvRequestIds matching a source rank. It may be
// used to allocate ArrayOfRecvRequests for the call to RecvRequestIds
  inline int MPI_Access::RecvRequestIdsSize(int sourcerank) {
         return _RecvRequests[sourcerank].size() ;
  }
// Returns the MPI_Datatype (registered in MPI in the constructor with
// MPI_Type_struct and MPI_Type_commit) for TimeMessages
  inline MPI_Datatype MPI_Access::TimeType() const {
         return _MPI_TIME ;
  }
// Returns true if the tag MPITag corresponds to a TimeMessage
  inline bool MPI_Access::IsTimeMessage( int MPITag ) const {
         return ((MPITag%ModuloTag) == _MessageTime) ; } ;
// Returns the MPI size of a TimeMessage
  inline MPI_Aint MPI_Access::TimeExtent() const {
         MPI_Aint extent ;
         MPI_Type_extent( _MPI_TIME , &extent ) ;
         return extent ;
  }
// Returns the MPI size of a MPI_INT
  inline MPI_Aint MPI_Access::IntExtent() const {
         MPI_Aint extent ;
         MPI_Type_extent( MPI_INT , &extent ) ;
         return extent ;
  }
// Returns the MPI size of a MPI_DOUBLE
  inline MPI_Aint MPI_Access::DoubleExtent() const {
         MPI_Aint extent ;
         MPI_Type_extent( MPI_DOUBLE , &extent ) ;
         return extent ;
  }
// Returns the MPI size of the MPI_Datatype datatype
  inline MPI_Aint MPI_Access::Extent( MPI_Datatype datatype ) const {
         if ( datatype == _MPI_TIME )
           return TimeExtent() ;
         if ( datatype == MPI_INT )
           return IntExtent() ;
         if ( datatype == MPI_DOUBLE )
           return DoubleExtent() ;
         return 0 ;
  }

// Returns the MPITag of the request corresponding to RequestId identifier
  inline int MPI_Access::MPITag( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           return aRequestStruct->MPITag ;
         }
         return -1 ;
  }
// Returns the MPITarget of the request corresponding to RequestId identifier
  inline int MPI_Access::MPITarget( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           return aRequestStruct->MPITarget ;
         }
         return -1 ;
  }
// Returns true if the request corresponding to RequestId identifier was [I]Recv
  inline bool MPI_Access::MPIIsRecv( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           return aRequestStruct->MPIIsRecv ;
         }
         return false ;
  }
// Returns true if the request corresponding to RequestId identifier was asynchronous
  inline bool MPI_Access::MPIAsynchronous( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           return aRequestStruct->MPIAsynchronous ;
         }
         return false ;
  }
// Returns true if the request corresponding to RequestId identifier was completed
  inline bool MPI_Access::MPICompleted( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           return aRequestStruct->MPICompleted ;
         }
         return true ;
  }
// Returns the MPI_Datatype  of the request corresponding to RequestId identifier
  inline MPI_Datatype MPI_Access::MPIDatatype( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           return aRequestStruct->MPIDatatype ; ;
         }
         return (MPI_Datatype ) NULL ;
  }
// Returns the size of the receiving message of the request corresponding to
// RequestId identifier
  inline int MPI_Access::MPIOutCount( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           return aRequestStruct->MPIOutCount ;
         }
         return 0 ;
  }

// Increments the previous tag value (cyclically)
// Look at MPI_Access::NewSendTag/NewRecvTag in MPI_Access.cxx
  inline int MPI_Access::IncrTag( int prevtag ) {
         int tag ;
         if ( (prevtag % ModuloTag) == _MessageTime ) {
           tag = ((prevtag/ModuloTag)*ModuloTag) ;
         }
         else {
           tag = ((prevtag/ModuloTag + 1)*ModuloTag) ;
         }
         if ( tag > _MaxMPITag )
           tag = _BaseMPITag ;
         return tag ;
  }
// Returns the MPITag with the method-type field
// Look at MPI_Access::NewSendTag/NewRecvTag in MPI_Access.cxx
  inline int MPI_Access::ValTag( int tag, int method ) {
         return ((tag/ModuloTag)*ModuloTag) + method ;
  }

// Remove a Request identifier from the list _RecvRequests/_SendRequests for
// the corresponding target.
  inline void MPI_Access::DeleteSendRecvRequest( int RequestId ) {
         if ( _Trace )
           std::cout << "MPI_Access::DeleteSendRecvRequest" << _MyRank
                     << "( " << RequestId << " ) " << std::endl ;
         if ( MPIIsRecv( RequestId ) ) {
           _RecvRequests[ MPITarget( RequestId ) ].remove( RequestId ) ;
         }
         else {
           _SendRequests[ MPITarget( RequestId ) ].remove( RequestId ) ;
         }
  }

// Delete the MPI structure MPI_Status * of a ReaquestId
  inline void MPI_Access::DeleteStatus( int RequestId ) {
         if ( _MapOfRequestStruct[RequestId]->MPIStatus != NULL ) {
           delete _MapOfRequestStruct[RequestId]->MPIStatus ;
           ClearMPIStatus( RequestId ) ;
         }
  }

// Returns the MPI structure MPI_Request * of a RequestId
  inline MPI_Request * MPI_Access::MPIRequest( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         //cout << "MPIRequest" << _MyRank << "(" << RequestId
         //     << " ) aRequestStruct->MPIRequest " << aRequestStruct->MPIRequest
         //     << endl ;
         if ( aRequestStruct ) {
           return &aRequestStruct->MPIRequest ;
         }
         return &mpirequestnull ;
  }
// Returns the MPI structure MPI_Status * of a RequestId
  inline MPI_Status * MPI_Access::MPIStatus( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         //cout << "MPIStatus" << _MyRank << "(" << RequestId
         //     << " ) aRequestStruct->MPIStatus " << aRequestStruct->MPIStatus
         //     << endl ;
         if ( aRequestStruct ) {
           return aRequestStruct->MPIStatus ; ;
         }
         return NULL ;
  }
// Set the MPICompleted field of the structure Request corresponding to RequestId
// identifier with the value completed
  inline void MPI_Access::SetMPICompleted( int RequestId , bool completed ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           aRequestStruct->MPICompleted = completed ;
         }
  }
// Set the MPIOutCount field of the structure Request corresponding to RequestId
// identifier with the value outcount
  inline void MPI_Access::SetMPIOutCount( int RequestId , int outcount ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           aRequestStruct->MPIOutCount = outcount ;
         }
  }
// Nullify the MPIStatusfield of the structure Request corresponding to RequestId
// identifier
  inline void MPI_Access::ClearMPIStatus( int RequestId ) {
         struct RequestStruct *aRequestStruct = _MapOfRequestStruct[ RequestId ] ;
         if ( aRequestStruct ) {
           aRequestStruct->MPIStatus = NULL ;
         }
  }

// Returns the _MessageIdent enum value corresponding to the MPI_Datatype datatype
// Look at MPI_Access::NewSendTag/NewRecvTag in MPI_Access.cxx
  inline _MessageIdent MPI_Access::MethodId( MPI_Datatype datatype ) const {
         _MessageIdent aMethodIdent ;
         if ( datatype == _MPI_TIME ) {
           aMethodIdent = _MessageTime ;
         }
         else if ( datatype == MPI_INT ) {
           aMethodIdent = _MessageInt ;
         }
         else if ( datatype == MPI_DOUBLE ) {
           aMethodIdent = _MessageDouble ;
         }
         else {
           aMethodIdent = _MessageUnknown ;
         }
         return aMethodIdent ;
  }
// Returns the MPI_Datatype corresponding to the _MessageIdent enum aMethodIdent
  inline MPI_Datatype MPI_Access::Datatype( _MessageIdent aMethodIdent ) const {
         MPI_Datatype aDataType ;
         switch( aMethodIdent ) {
           case _MessageTime :
             aDataType = _MPI_TIME ;
             break ;
           case _MessageInt :
             aDataType = MPI_INT ;
             break ;
           case _MessageDouble :
             aDataType = MPI_DOUBLE ;
             break ;
           default :
             aDataType = (MPI_Datatype) -1 ;
             break ;
         }
         return aDataType ;
  }

ostream & operator<< (ostream &,const _MessageIdent &);

ostream & operator<< (ostream &,const TimeMessage &);

}

#endif /*MPI_ACCESS_HXX_*/
