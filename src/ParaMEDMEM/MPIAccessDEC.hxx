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

#ifndef __MPIACCESSDEC_HXX__
#define __MPIACCESSDEC_HXX__

#include "MPIAccess.hxx"
#include "DEC.hxx"
#include "LinearTimeInterpolator.hxx"

#include <map>
#include <iostream>

namespace MEDCoupling
{
  /*
   * Internal class, not part of the public API.
   *
   * Another gateway to the MPI library?
   */
  class MPIAccessDEC
  {
  public:  
    MPIAccessDEC( const ProcessorGroup& local_group, const ProcessorGroup& distant_group,
                  bool Asynchronous = true );
    virtual ~MPIAccessDEC();
    MPIAccess * getMPIAccess() { return _MPI_access; }
    const MPI_Comm* getComm() { return _MPI_union_group->getComm(); }
    void asynchronous( bool Asynchronous = true ) { _asynchronous = Asynchronous; }
    void setTimeInterpolator( TimeInterpolationMethod anInterp , double InterpPrecision=0 ,
                              int n_step_before=1, int nStepAfter=1 );

    void setTime( double t ) { _t = t; _dt = -1; }
    void setTime( double t , double dt ) { _t = t; _dt = dt; }
    bool outOfTime( int target ) { return (*_out_of_time)[target]; }

    int send( void* sendbuf, int sendcount , MPI_Datatype sendtype , int target );
    int recv( void* recvbuf, int recvcount , MPI_Datatype recvtype , int target );
    int recv( void* recvbuf, int recvcount , MPI_Datatype recvtype , int target ,
              int &RecvRequestId , bool Asynchronous=false );
    int sendRecv( void* sendbuf, int sendcount , MPI_Datatype sendtype ,
                  void* recvbuf, int recvcount , MPI_Datatype recvtype , int target );

    int allToAll( void* sendbuf, int sendcount, MPI_Datatype sendtype ,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype );
    int allToAllv( void* sendbuf, int* sendcounts, int* sdispls, MPI_Datatype sendtype ,
                   void* recvbuf, int* recvcounts, int* rdispls, MPI_Datatype recvtype );

    int allToAllTime( void* sendbuf, int sendcount , MPI_Datatype sendtype ,
                      void* recvbuf, int recvcount , MPI_Datatype recvtype );
    int allToAllvTime( void* sendbuf, int* sendcounts, int* sdispls,
                       MPI_Datatype sendtype ,
                       void* recvbuf, int* recvcounts, int* rdispls,
                       MPI_Datatype recvtype );
    int checkTime( int recvcount , MPI_Datatype recvtype , int target , bool UntilEnd );
    int checkSent(bool WithWait=false);
    int checkFinalSent() { return checkSent( true ); }
    int checkFinalRecv();
  protected:
    int send( void* sendbuf, int sendcount , int sendoffset , MPI_Datatype sendtype ,
              int target, int &SendRequestId );
    int recv( void* recvbuf, int recvcount , int recvoffset , MPI_Datatype recvtype ,
              int target, int &RecvRequestId );
    int sendRecv( void* sendbuf, int sendcount , int sendoffset ,
                  MPI_Datatype sendtype , 
                  void* recvbuf, int recvcount , int recvoffset ,
                  MPI_Datatype recvtype , int target ,
                  int &SendRequestId ,int &RecvRequestId );
  private :
    bool _asynchronous;
    MPIProcessorGroup* _MPI_union_group;

    TimeInterpolator* _time_interpolator;
    int _n_step_before;
    int _n_step_after;

    int _my_rank;
    int _group_size;
    MPIAccess* _MPI_access;

    // Current time and deltatime of current process
    double _t;
    double _dt;

    // TimeMessages from each target _TimeMessages[target][Step] : TimeMessage
    std::vector< std::vector< TimeMessage > > *_time_messages;
    // Corresponding DataMessages from each target _DataMessages[target][~TimeStep]
    std::vector< bool >* _out_of_time;
    std::vector< int >* _data_messages_recv_count;
    std::vector< MPI_Datatype >* _data_messages_type;
    std::vector< std::vector< void * > >* _data_messages;

    typedef struct
    {
      void * SendBuffer;
      int Counter;
      MPI_Datatype DataType; }
      SendBuffStruct;
    std::map< int ,  SendBuffStruct * > *_map_of_send_buffers;
  };

  inline int MPIAccessDEC::send( void* sendbuf, int sendcount , MPI_Datatype sendtype , int target )
  {
    int SendRequestId;
    int sts;
    if ( _asynchronous )
      {
        sts = _MPI_access->ISend( sendbuf , sendcount , sendtype , target ,
                                  SendRequestId );
      }
    else
      {
        sts = _MPI_access->send( sendbuf , sendcount , sendtype , target ,
                                 SendRequestId );
        if ( sts == MPI_SUCCESS )
          free( sendbuf );
      }
    return sts;
  }

  inline int MPIAccessDEC::recv( void* recvbuf, int recvcount , MPI_Datatype recvtype , int target )
  {
    int RecvRequestId;
    int sts;
    if ( _asynchronous )
      sts = _MPI_access->IRecv( recvbuf , recvcount , recvtype , target , RecvRequestId );
    else
      sts = _MPI_access->recv( recvbuf , recvcount , recvtype , target ,  RecvRequestId );
    return sts;
  }

  inline int MPIAccessDEC::recv( void* recvbuf, int recvcount , MPI_Datatype recvtype ,
                                 int target ,  int &RecvRequestId , bool Asynchronous )
  {
    int sts;
    if ( Asynchronous )
      sts = _MPI_access->IRecv( recvbuf , recvcount , recvtype , target ,
                                RecvRequestId );
    else
      sts = _MPI_access->recv( recvbuf , recvcount , recvtype , target ,
                               RecvRequestId );
    return sts;
  }
  
  inline int MPIAccessDEC::sendRecv( void* sendbuf, int sendcount , MPI_Datatype sendtype ,
                                     void* recvbuf, int recvcount , MPI_Datatype recvtype ,
                                     int target )
  {
    int SendRequestId;
    int RecvRequestId;
    int sts;
    if ( _asynchronous )
      sts = _MPI_access->ISendRecv( sendbuf , sendcount , sendtype , target ,
                                    SendRequestId ,
                                    recvbuf , recvcount , recvtype , target ,
                                    RecvRequestId );
    else
      sts = _MPI_access->sendRecv( sendbuf , sendcount , sendtype , target ,
                                   SendRequestId ,
                                   recvbuf , recvcount , recvtype , target ,
                                   RecvRequestId );
    return sts;
  }

  std::ostream & operator<< (std::ostream &,const TimeInterpolationMethod &);
}

#endif
