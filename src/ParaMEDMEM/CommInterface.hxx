// Copyright (C) 2007-2020  CEA/DEN, EDF R&D
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

#pragma once

#include "ParaIdType.hxx"
#include "MEDCouplingMemArray.hxx"

#include <mpi.h>

#include <memory>
#include <numeric>

namespace MEDCoupling
{
  template<class T>
  struct ParaTraits
  {
    using EltType = T;
  };
  
  template<>
  struct ParaTraits<double>
  {
    static MPI_Datatype MPIDataType;
  };

  template<>
  struct ParaTraits<Int32>
  {
    static MPI_Datatype MPIDataType;
  };

  template<>
  struct ParaTraits<Int64>
  {
    static MPI_Datatype MPIDataType;
  };

  class CommInterface
  {
  public:
    CommInterface() { }
    virtual ~CommInterface() { }
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
    int gather(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, int recvcount, MPI_Datatype recvtype, int root, MPI_Comm comm) const { return MPI_Gather(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,root,comm); }
    int gatherV(const void* sendbuf, int sendcount, MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root, MPI_Comm comm) const { return MPI_Gatherv(sendbuf,sendcount,sendtype,recvbuf,recvcounts,displs,recvtype,root,comm); }
    int allGather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm) const { return MPI_Allgather(sendbuf,sendcount, sendtype, recvbuf, recvcount, recvtype, comm); }
    int allGatherV(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
                   const int displs[], MPI_Datatype recvtype, MPI_Comm comm) const { return MPI_Allgatherv(sendbuf,sendcount,sendtype,recvbuf,recvcounts,displs,recvtype,comm); }
    int allToAll(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype,
                 MPI_Comm comm) const { return MPI_Alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, comm); }
    int allToAllV(const void* sendbuf, int* sendcounts, int* senddispls,
                  MPI_Datatype sendtype, void* recvbuf, int* recvcounts,
                  int* recvdispls, MPI_Datatype recvtype, 
                  MPI_Comm comm) const { return MPI_Alltoallv(sendbuf, sendcounts, senddispls, sendtype, recvbuf, recvcounts, recvdispls, recvtype, comm); }

    int reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype,
               MPI_Op op, int root, MPI_Comm comm) const { return MPI_Reduce(sendbuf, recvbuf, count, datatype, op, root, comm); }
    int allReduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) const { return MPI_Allreduce(sendbuf, recvbuf, count, datatype, op, comm); }
  public:
    void gatherArrays(MPI_Comm comm, int root, const DataArrayIdType *array, std::vector< MCAuto<DataArrayIdType> >& arraysOut) const;
    void allGatherArrays(MPI_Comm comm, const DataArrayIdType *array, std::vector< MCAuto<DataArrayIdType> >& arraysOut) const;
    int allGatherArrays(MPI_Comm comm, const DataArrayIdType *array, std::unique_ptr<mcIdType[]>& result, std::unique_ptr<mcIdType[]>& resultIndex) const;
    void allToAllArrays(MPI_Comm comm, const std::vector< MCAuto<DataArrayIdType> >& arrays, std::vector< MCAuto<DataArrayIdType> >& arraysOut) const;
    void allToAllArrays(MPI_Comm comm, const std::vector< MCAuto<DataArrayDouble> >& arrays, std::vector< MCAuto<DataArrayDouble> >& arraysOut) const;
    void allToAllArrays(MPI_Comm comm, const std::vector< MCAuto<DataArrayDouble> >& arrays, MCAuto<DataArrayDouble>& arraysOut) const;
    
    template<class T>
    int gatherArraysT(MPI_Comm comm, int root, const typename Traits<T>::ArrayType *array, std::unique_ptr<T[]>& result, std::unique_ptr<mcIdType[]>& resultIndex, int& rank) const
    {
      using DataArrayT = typename Traits<T>::ArrayType;
      int size;
      this->commSize(comm,&size);
      rank = -1;
      this->commRank(comm,&rank);
      std::unique_ptr<mcIdType[]> nbOfElems;
      if(rank==root)
        nbOfElems.reset(new mcIdType[size]);
      mcIdType nbOfCellsRequested(array->getNumberOfTuples());
      this->gather(&nbOfCellsRequested,1,MPI_ID_TYPE,nbOfElems.get(),1,MPI_ID_TYPE,root,comm);
      std::unique_ptr<int[]> nbOfElemsInt,offsetsIn;
      if(rank==root)
      {
        mcIdType nbOfCellIdsSum(std::accumulate(nbOfElems.get(),nbOfElems.get()+size,0));
        result.reset(new T[nbOfCellIdsSum]);
        nbOfElemsInt = CommInterface::ToIntArray<mcIdType>(nbOfElems,size);
        offsetsIn = CommInterface::ComputeOffset(nbOfElemsInt,size);
      }
      this->gatherV(array->begin(),nbOfCellsRequested,ParaTraits<T>::MPIDataType,result.get(),nbOfElemsInt.get(),offsetsIn.get(),ParaTraits<T>::MPIDataType,root,comm);
      if(rank==root)
      {
        resultIndex = ComputeOffsetFull<mcIdType>(nbOfElems,size);
      }
      return size;
    }

    template<class T>
    void gatherArraysT2(MPI_Comm comm, int root, const typename Traits<T>::ArrayType *array, std::vector< MCAuto<typename Traits<T>::ArrayType> >& arraysOut) const
    {
      using DataArrayT = typename Traits<T>::ArrayType;
      std::unique_ptr<T[]> result;
      std::unique_ptr<mcIdType[]> resultIndex;
      int rank(-1);
      int size(this->gatherArraysT<T>(comm,root,array,result,resultIndex,rank));
      arraysOut.resize(size);
      for(int i = 0 ; i < size ; ++i)
      {
        arraysOut[i] = DataArrayT::New();
        if(rank == root)
        {
          mcIdType nbOfEltPack(resultIndex[i+1]-resultIndex[i]);
          arraysOut[i]->alloc(nbOfEltPack,1);
          std::copy(result.get()+resultIndex[i],result.get()+resultIndex[i+1],arraysOut[i]->getPointer());
        }
      }
    }

    template<class T>
    int allGatherArraysT(MPI_Comm comm, const typename Traits<T>::ArrayType *array, std::unique_ptr<T[]>& result, std::unique_ptr<mcIdType[]>& resultIndex) const
    {
      using DataArrayT = typename Traits<T>::ArrayType;
      int size;
      this->commSize(comm,&size);
      std::unique_ptr<mcIdType[]> nbOfElems(new mcIdType[size]);
      mcIdType nbOfCellsRequested(array->getNumberOfTuples());
      this->allGather(&nbOfCellsRequested,1,MPI_ID_TYPE,nbOfElems.get(),1,MPI_ID_TYPE,comm);
      mcIdType nbOfCellIdsSum(std::accumulate(nbOfElems.get(),nbOfElems.get()+size,0));
      result.reset(new T[nbOfCellIdsSum]);
      std::unique_ptr<int[]> nbOfElemsInt( CommInterface::ToIntArray<mcIdType>(nbOfElems,size) );
      std::unique_ptr<int[]> offsetsIn( CommInterface::ComputeOffset(nbOfElemsInt,size) );
      this->allGatherV(array->begin(),nbOfCellsRequested,ParaTraits<T>::MPIDataType,result.get(),nbOfElemsInt.get(),offsetsIn.get(),ParaTraits<T>::MPIDataType,comm);
      resultIndex = ComputeOffsetFull<mcIdType>(nbOfElems,size);
      return size;
    }

    template<class T>
    void allGatherArraysT2(MPI_Comm comm, const typename Traits<T>::ArrayType *array, std::vector< MCAuto<typename Traits<T>::ArrayType> >& arraysOut) const
    {
      using DataArrayT = typename Traits<T>::ArrayType;
      std::unique_ptr<T[]> result;
      std::unique_ptr<mcIdType[]> resultIndex;
      int size(this->allGatherArraysT<T>(comm,array,result,resultIndex));
      arraysOut.resize(size);
      for(int i = 0 ; i < size ; ++i)
      {
        arraysOut[i] = DataArrayT::New();
        mcIdType nbOfEltPack(resultIndex[i+1]-resultIndex[i]);
        arraysOut[i]->alloc(nbOfEltPack,1);
        std::copy(result.get()+resultIndex[i],result.get()+resultIndex[i+1],arraysOut[i]->getPointer());
      }
    }

    template<class T>
    int allToAllArraysT2(MPI_Comm comm, const std::vector< MCAuto<typename Traits<T>::ArrayType> >& arrays, MCAuto<typename Traits<T>::ArrayType>& arrayOut, std::unique_ptr<mcIdType[]>& nbOfElems2, mcIdType& nbOfComponents) const
    {
      using DataArrayT = typename Traits<T>::ArrayType;
      int size;
      this->commSize(comm,&size);
      if( arrays.size() != ToSizeT(size) )
        throw INTERP_KERNEL::Exception("AllToAllArrays : internal error ! Invalid size of input array.");
        
      std::vector< const DataArrayT *> arraysBis(FromVecAutoToVecOfConst<DataArrayT>(arrays));
      std::unique_ptr<mcIdType[]> nbOfElems3(new mcIdType[size]);
      nbOfElems2.reset(new mcIdType[size]);
      nbOfComponents = std::numeric_limits<mcIdType>::max();
      for(int curRk = 0 ; curRk < size ; ++curRk)
      {
        mcIdType curNbOfCompo( ToIdType( arrays[curRk]->getNumberOfComponents() ) );
        if(nbOfComponents != std::numeric_limits<mcIdType>::max())
        {
          if( nbOfComponents != curNbOfCompo )
            throw INTERP_KERNEL::Exception("AllToAllArrays : internal error ! Nb of components is not homogeneous !");
        }
        else
        {
          nbOfComponents = curNbOfCompo;
        }
        nbOfElems3[curRk] = arrays[curRk]->getNbOfElems();
      }
      this->allToAll(nbOfElems3.get(),1,MPI_ID_TYPE,nbOfElems2.get(),1,MPI_ID_TYPE,comm);
      mcIdType nbOfCellIdsSum(std::accumulate(nbOfElems2.get(),nbOfElems2.get()+size,0));
      arrayOut = DataArrayT::New();
      arrayOut->alloc(nbOfCellIdsSum/nbOfComponents,nbOfComponents);
      std::unique_ptr<int[]> nbOfElemsInt( CommInterface::ToIntArray<mcIdType>(nbOfElems3,size) ),nbOfElemsOutInt( CommInterface::ToIntArray<mcIdType>(nbOfElems2,size) );
      std::unique_ptr<int[]> offsetsIn( CommInterface::ComputeOffset(nbOfElemsInt,size) ), offsetsOut( CommInterface::ComputeOffset(nbOfElemsOutInt,size) );
      {
        MCAuto<DataArrayT> arraysAcc(DataArrayT::Aggregate(arraysBis));
        this->allToAllV(arraysAcc->begin(),nbOfElemsInt.get(),offsetsIn.get(),ParaTraits<T>::MPIDataType,
                        arrayOut->getPointer(),nbOfElemsOutInt.get(),offsetsOut.get(),ParaTraits<T>::MPIDataType,comm);
      }
      return size;
    }

    template<class T>
    void allToAllArraysT(MPI_Comm comm, const std::vector< MCAuto<typename Traits<T>::ArrayType> >& arrays, std::vector< MCAuto<typename Traits<T>::ArrayType> >& arraysOut) const
    {
      using DataArrayT = typename Traits<T>::ArrayType;
      MCAuto<DataArrayT> cellIdsFromProcs;
      std::unique_ptr<mcIdType[]> nbOfElems2;
      mcIdType nbOfComponents(0);
      int size(this->allToAllArraysT2<T>(comm,arrays,cellIdsFromProcs,nbOfElems2,nbOfComponents));
      std::unique_ptr<mcIdType[]> offsetsOutIdType( CommInterface::ComputeOffset(nbOfElems2,size) );
      // build output arraysOut by spliting cellIdsFromProcs into parts
      arraysOut.resize(size);
      for(int curRk = 0 ; curRk < size ; ++curRk)
      {
        arraysOut[curRk] = DataArrayT::NewFromArray(cellIdsFromProcs->begin()+offsetsOutIdType[curRk],cellIdsFromProcs->begin()+offsetsOutIdType[curRk]+nbOfElems2[curRk]);
        arraysOut[curRk]->rearrange(nbOfComponents);
      }
    }
  public:

    /*!
    * \a counts is expected to be an array of array length. This method returns an array of split array.
    */
    static std::unique_ptr<mcIdType[]> SplitArrayOfLength(const std::unique_ptr<mcIdType[]>& counts, std::size_t countsSz, int rk, int size)
    {
      std::unique_ptr<mcIdType[]> ret(new mcIdType[countsSz]);
      for(std::size_t i=0;i<countsSz;++i)
      {
        mcIdType a,b;
        DataArray::GetSlice(0,counts[i],1,rk,size,a,b);
        ret[i] = b-a;
      }
      return ret;
    }

    /*!
    * Helper of alltoallv and allgatherv
    */
    template<class T>
    static std::unique_ptr<int []> ToIntArray(const std::unique_ptr<T []>& arr, std::size_t size)
    {
      std::unique_ptr<int []> ret(new int[size]);
      std::copy(arr.get(),arr.get()+size,ret.get());
      return ret;
    }
    
    /*!
    * Helper of alltoallv and allgatherv
    */
    template<class T>
    static std::unique_ptr<T []> ComputeOffset(const std::unique_ptr<T []>& counts, std::size_t sizeOfCounts)
    {
      std::unique_ptr<T []> ret(new T[sizeOfCounts]);
      ret[0] = static_cast<T>(0);
      for(std::size_t i = 1 ; i < sizeOfCounts ; ++i)
      {
        ret[i] = ret[i-1] + counts[i-1];
      }
      return ret;
    }

    /*!
    * Helper of alltoallv and allgatherv
    */
    template<class T>
    static std::unique_ptr<T []> ComputeOffsetFull(const std::unique_ptr<T []>& counts, std::size_t sizeOfCounts)
    {
      std::unique_ptr<T []> ret(new T[sizeOfCounts+1]);
      ret[0] = static_cast<T>(0);
      for(std::size_t i = 1 ; i < sizeOfCounts+1 ; ++i)
      {
        ret[i] = ret[i-1] + counts[i-1];
      }
      return ret;
    }
  };
}
