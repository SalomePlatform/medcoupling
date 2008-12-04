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
#ifndef MPI_ACESSDEC_HXX_
#define MPI_ACESSDEC_HXX_

#include <map>
#include <iostream>

#include "MPI_Access.hxx"
#include "DEC.hxx"
#include "LinearTimeInterpolator.hxx"

namespace ParaMEDMEM {

	//  typedef enum{WithoutTimeInterp,LinearTimeInterp} TimeInterpolationMethod;

  class MPI_AccessDEC {

    public:  
      MPI_AccessDEC( const ProcessorGroup& local_group, const ProcessorGroup& distant_group,
                     bool Asynchronous = true ) ;
      //MPI_AccessDEC( ProcessorGroup& local_group, ProcessorGroup& distant_group,
      //               TimeInterpolator * aTimeInterpolator ,
      //               bool Asynchronous = true ) ;
      virtual ~MPI_AccessDEC();
      MPI_Access * MPIAccess() ;
      const MPI_Comm* GetComm() ;

      void Asynchronous( bool Asynchronous = true ) ;
      void SetTimeInterpolator( TimeInterpolationMethod anInterp , double InterpPrecision=0 ,
                                int nStepBefore=1, int nStepAfter=1 ) ;

      void SetTime( double t ) ;
      void SetTime( double t , double dt ) ;
      bool OutOfTime( int target ) ;

      int Send( void* sendbuf, int sendcount , MPI_Datatype sendtype , int target ) ;
      int Recv( void* recvbuf, int recvcount , MPI_Datatype recvtype , int target ) ;
      int Recv( void* recvbuf, int recvcount , MPI_Datatype recvtype , int target ,
                int &RecvRequestId , bool Asynchronous=false ) ;
      int SendRecv( void* sendbuf, int sendcount , MPI_Datatype sendtype ,
                    void* recvbuf, int recvcount , MPI_Datatype recvtype , int target ) ;

      int AllToAll( void* sendbuf, int sendcount, MPI_Datatype sendtype ,
	            void* recvbuf, int recvcount, MPI_Datatype recvtype ) ;
      int AllToAllv( void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype sendtype ,
	             void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype recvtype ) ;

      int AllToAllTime( void* sendbuf, int sendcount , MPI_Datatype sendtype ,
	                void* recvbuf, int recvcount , MPI_Datatype recvtype ) ;
      int AllToAllvTime( void* sendbuf, int* sendcounts, int* sdispls,
                         MPI_Datatype sendtype ,
	                 void* recvbuf, int* recvcounts, int* rdispls,
                         MPI_Datatype recvtype ) ;
//      int CheckTime( int recvcount , int recvsize , MPI_Datatype recvtype , int target ,
//                     bool &OutOfTime ) ;
      int CheckTime( int recvcount , MPI_Datatype recvtype , int target , bool UntilEnd ) ;
      int CheckSent(bool WithWait=false) ;
      int CheckFinalSent() {
          return CheckSent( true ) ; } ;
      int CheckFinalRecv() ;

    protected:
//      int Send( void* sendbuf, int sendcount , int sendoffset , MPI_Datatype sendtype ,
//                int target ) ;
//      int Recv( void* recvbuf, int recvcount , int recvoffset , MPI_Datatype recvtype ,
//                int target ) ;

      int Send( void* sendbuf, int sendcount , int sendoffset , MPI_Datatype sendtype ,
                int target, int &SendRequestId ) ;
      int Recv( void* recvbuf, int recvcount , int recvoffset , MPI_Datatype recvtype ,
                int target, int &RecvRequestId ) ;
      int SendRecv( void* sendbuf, int sendcount , int sendoffset ,
                    MPI_Datatype sendtype , 
                    void* recvbuf, int recvcount , int recvoffset ,
                    MPI_Datatype recvtype , int target ,
                    int &SendRequestId ,int &RecvRequestId ) ;

    private :
      bool                   _Asynchronous ;
      const ProcessorGroup * _local_group ;
      const ProcessorGroup * _distant_group ;
      MPIProcessorGroup    * _MPI_union_group ;

      TimeInterpolator     * _TimeInterpolator ;
      int                    _nStepBefore ;
      int                    _nStepAfter ;

      int                    _MyRank ;
      int                    _GroupSize ;
      MPI_Access           * _MPIAccess ;

// Current time and deltatime of current process
      double                 _t ;
      double                 _dt ;

// TimeMessages from each target _TimeMessages[target][Step] : TimeMessage
      vector< vector< TimeMessage > > *_TimeMessages ;
// Corresponding DataMessages from each target _DataMessages[target][~TimeStep]
      vector< bool >             * _OutOfTime ;
      vector< int >              * _DataMessagesRecvCount ;
      vector< MPI_Datatype >     * _DataMessagesType ;
      vector< vector< void * > > *_DataMessages ;

      typedef struct { void * SendBuffer ;
                       int Counter ;
                       MPI_Datatype DataType ; } SendBuffStruct ;
// int RequestId -> SendBuffStruct :
      map< int ,  SendBuffStruct * > *_MapOfSendBuffers ;

  };

  inline MPI_Access * MPI_AccessDEC::MPIAccess() {
         return _MPIAccess ; } ;
  inline const MPI_Comm* MPI_AccessDEC::GetComm() {
         return _MPI_union_group->getComm() ; } ;

  inline void MPI_AccessDEC::Asynchronous( bool Asynchronous ) {
         _Asynchronous = Asynchronous ; } ;

  inline void MPI_AccessDEC::SetTime( double t ) {
         _t = t ; _dt = -1 ; } ;
  inline void MPI_AccessDEC::SetTime( double t , double dt ) {
         _t = t ; _dt = dt ; } ;
  inline bool MPI_AccessDEC::OutOfTime( int target ) {
         return (*_OutOfTime)[target] ; } ;

  inline int MPI_AccessDEC::Send( void* sendbuf, int sendcount , MPI_Datatype sendtype ,
                                  int target ) {
         int SendRequestId ;
         int sts ;
         if ( _Asynchronous ) {
           sts = _MPIAccess->ISend( sendbuf , sendcount , sendtype , target ,
                                    SendRequestId ) ;
         }
         else {
           sts = _MPIAccess->Send( sendbuf , sendcount , sendtype , target ,
                                   SendRequestId ) ;
           if ( sts == MPI_SUCCESS ) {
             free( sendbuf ) ;
           }
         }
         return sts ; } ;
  inline int MPI_AccessDEC::Recv( void* recvbuf, int recvcount , MPI_Datatype recvtype ,
                                  int target )  {
         int RecvRequestId ;
         int sts ;
         if ( _Asynchronous ) {
           sts = _MPIAccess->IRecv( recvbuf , recvcount , recvtype , target ,
                                    RecvRequestId ) ;
         }
         else {
           sts = _MPIAccess->Recv( recvbuf , recvcount , recvtype , target ,
                                   RecvRequestId ) ;
         }
         return sts ; } ;
  inline int MPI_AccessDEC::Recv( void* recvbuf, int recvcount , MPI_Datatype recvtype ,
                                  int target ,  int &RecvRequestId , bool Asynchronous )  {
         int sts ;
         if ( Asynchronous ) {
           sts = _MPIAccess->IRecv( recvbuf , recvcount , recvtype , target ,
                                    RecvRequestId ) ;
         }
         else {
           sts = _MPIAccess->Recv( recvbuf , recvcount , recvtype , target ,
                                   RecvRequestId ) ;
         }
         return sts ; } ;
  inline int MPI_AccessDEC::SendRecv( void* sendbuf, int sendcount , MPI_Datatype sendtype ,
                                      void* recvbuf, int recvcount , MPI_Datatype recvtype ,
                                      int target ) {
         int SendRequestId ;
         int RecvRequestId ;
         int sts ;
         if ( _Asynchronous ) {
           sts = _MPIAccess->ISendRecv( sendbuf , sendcount , sendtype , target ,
                                        SendRequestId ,
                                        recvbuf , recvcount , recvtype , target ,
                                        RecvRequestId ) ;
         }
         else {
           sts = _MPIAccess->SendRecv( sendbuf , sendcount , sendtype , target ,
                                       SendRequestId ,
                                       recvbuf , recvcount , recvtype , target ,
                                       RecvRequestId ) ;
         }
         return sts ; } ;

//  inline int MPI_AccessDEC::Send( void* sendbuf, int sendcount , int sendoffset ,
//                                  MPI_Datatype sendtype , int target ) {
//         int SendRequestId ;
//         return Send( sendbuf , sendcount , sendoffset , sendtype , target ,
//                      SendRequestId ) ; } ;
//  inline int MPI_AccessDEC::Recv( void* recvbuf, int recvcount , int recvoffset ,
//                                  MPI_Datatype recvtype , int target )  {
//         int RecvRequestId ;
//         return Recv( recvbuf , recvcount , recvoffset , recvtype , target ,
//                      RecvRequestId ) ; } ;


ostream & operator<< (ostream &,const TimeInterpolationMethod &);

}

#endif
