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

#ifndef __COMMINTERFACE_HXX__
#define __COMMINTERFACE_HXX__

#include <mpi.h>
namespace MEDCoupling
{

  class CommInterface
  {
  public:
    CommInterface(){}
    virtual ~CommInterface(){}
    int worldSize() const {
      int size;
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      return size;}
    int commSize(MPI_Comm comm, int* size) const { return MPI_Comm_size(comm,size); }
    int commRank(MPI_Comm comm, int* rank) const { return MPI_Comm_rank(comm,rank); }
    int commGroup(MPI_Comm comm, MPI_Group* group) const { return MPI_Comm_group(comm, group); }
    int groupIncl(MPI_Group group, int size, int* ranks, MPI_Group* group_output) const { return MPI_Group_incl(group, size, ranks, group_output); }
    int commCreate(MPI_Comm comm, MPI_Group group, MPI_Comm* comm_output) const { return MPI_Comm_create(comm,group,comm_output); }
    int groupFree(MPI_Group* group) const { return MPI_Group_free(group); }
    int commFree(MPI_Comm* comm) const { return MPI_Comm_free(comm); }

    int send(void* buffer, int count, MPI_Datatype datatype, int target, int tag, MPI_Comm comm) const { return MPI_Send(buffer,count, datatype, target, tag, comm); }
    int recv(void* buffer, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status) const { return MPI_Recv(buffer,count, datatype, source, tag, comm, status); }
    int sendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
                 int dest, int sendtag, void* recvbuf, int recvcount, 
                 MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
                 MPI_Status* status) { return MPI_Sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount, recvtype, source, recvtag, comm,status); }

    int Isend(void* buffer, int count, MPI_Datatype datatype, int target,
              int tag, MPI_Comm comm, MPI_Request *request) const { return MPI_Isend(buffer,count, datatype, target, tag, comm, request); }
    int Irecv(void* buffer, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request* request) const { return MPI_Irecv(buffer,count, datatype, source, tag, comm, request); }

    int wait(MPI_Request *request, MPI_Status *status) const { return MPI_Wait(request, status); }
    int test(MPI_Request *request, int *flag, MPI_Status *status) const { return MPI_Test(request, flag, status); }
    int requestFree(MPI_Request *request) const { return MPI_Request_free(request); }
    int waitany(int count, MPI_Request *array_of_requests, int *index, MPI_Status *status) const { return MPI_Waitany(count, array_of_requests, index, status); }
    int testany(int count, MPI_Request *array_of_requests, int *index, int *flag, MPI_Status *status) const { return MPI_Testany(count, array_of_requests, index, flag, status); }
    int waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_status) const { return MPI_Waitall(count, array_of_requests, array_of_status); }
    int testall(int count, MPI_Request *array_of_requests, int *flag, MPI_Status *array_of_status) const { return MPI_Testall(count, array_of_requests, flag, array_of_status); }
    int waitsome(int incount, MPI_Request *array_of_requests,int *outcount, int *array_of_indices, MPI_Status *array_of_status) const { return MPI_Waitsome(incount, array_of_requests, outcount, array_of_indices, array_of_status); }
    int testsome(int incount, MPI_Request *array_of_requests, int *outcount,
                 int *array_of_indices, MPI_Status *array_of_status) const { return MPI_Testsome(incount, array_of_requests, outcount, array_of_indices, array_of_status); }
    int probe(int source, int tag, MPI_Comm comm, MPI_Status *status) const { return MPI_Probe(source, tag, comm, status) ; }
    int Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status) const { return MPI_Iprobe(source, tag, comm, flag, status) ; }
    int cancel(MPI_Request *request) const { return MPI_Cancel(request); }
    int testCancelled(MPI_Status *status, int *flag) const { return MPI_Test_cancelled(status, flag); }
    int barrier(MPI_Comm comm) const { return MPI_Barrier(comm); }
    int errorString(int errorcode, char *string, int *resultlen) const { return MPI_Error_string(errorcode, string, resultlen); }
    int getCount(MPI_Status *status, MPI_Datatype datatype, int *count) const { return MPI_Get_count(status, datatype, count); }

    int broadcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) const { return MPI_Bcast(buffer, count,  datatype, root, comm); }
    int allGather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm) const { return MPI_Allgather(sendbuf,sendcount, sendtype, recvbuf, recvcount, recvtype, comm); }
    int allToAll(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype,
                 MPI_Comm comm) const { return MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm); }
    int allToAllV(void* sendbuf, int* sendcounts, int* senddispls,
                  MPI_Datatype sendtype, void* recvbuf, int* recvcounts,
                  int* recvdispls, MPI_Datatype recvtype, 
                  MPI_Comm comm) const { return MPI_Alltoallv(sendbuf, sendcounts, senddispls, sendtype, recvbuf, recvcounts, recvdispls, recvtype, comm); }

    int reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm) const { return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm); }
    int allReduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) const { return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm); }
  };
}

#endif /*COMMINTERFACE_HXX_*/
