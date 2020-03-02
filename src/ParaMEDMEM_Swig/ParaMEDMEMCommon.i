// Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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

%include std_set.i

%template() std::set<int>;

%{
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "Topology.hxx"
#include "MPIProcessorGroup.hxx"
#include "DEC.hxx"
#include "InterpKernelDEC.hxx"
#include "NonCoincidentDEC.hxx"
#include "StructuredCoincidentDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ICoCoMEDField.hxx"
#include "ComponentTopology.hxx"
#include "ParaUMesh.hxx"

using namespace INTERP_KERNEL;
using namespace MEDCoupling;
using namespace ICoCo;
%}

%include "InterpolationOptions.hxx"
%include "ProcessorGroup.hxx"
%include "DECOptions.hxx"
%include "ParaMESH.hxx"
%include "ParaFIELD.hxx"
%include "MPIProcessorGroup.hxx"
%include "ComponentTopology.hxx"
%include "DEC.hxx"
%include "DisjointDEC.hxx"
%include "InterpKernelDEC.hxx"
%include "StructuredCoincidentDEC.hxx"

%include "ICoCoField.hxx"
%rename(ICoCoMEDField) ICoCo::MEDField;
%include "ICoCoMEDField.hxx"

%newobject MEDCoupling::ParaUMesh::getCellIdsLyingOnNodes;

%nodefaultctor;

namespace MEDCoupling
{
  class CommInterface
  {
  public:
    CommInterface();
    virtual ~CommInterface();
    int worldSize() const;
    int commSize(MPI_Comm comm, int* size) const;
    int commRank(MPI_Comm comm, int* rank) const;
    int commGroup(MPI_Comm comm, MPI_Group* group) const;
    int groupIncl(MPI_Group group, int size, int* ranks, MPI_Group* group_output) const;
    int commCreate(MPI_Comm comm, MPI_Group group, MPI_Comm* comm_output) const;
    int groupFree(MPI_Group* group) const;
    int commFree(MPI_Comm* comm) const;

    int send(void* buffer, int count, MPI_Datatype datatype, int target, int tag, MPI_Comm comm) const;
    int recv(void* buffer, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status* status) const;
    int sendRecv(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
                 int dest, int sendtag, void* recvbuf, int recvcount, 
                 MPI_Datatype recvtype, int source, int recvtag, MPI_Comm comm,
                 MPI_Status* status);

    int Isend(void* buffer, int count, MPI_Datatype datatype, int target,
              int tag, MPI_Comm comm, MPI_Request *request) const;
    int Irecv(void* buffer, int count, MPI_Datatype datatype, int source,
              int tag, MPI_Comm comm, MPI_Request* request) const;

    int wait(MPI_Request *request, MPI_Status *status) const;
    int test(MPI_Request *request, int *flag, MPI_Status *status) const;
    int requestFree(MPI_Request *request) const;
    int waitany(int count, MPI_Request *array_of_requests, int *index, MPI_Status *status) const;
    int testany(int count, MPI_Request *array_of_requests, int *index, int *flag, MPI_Status *status) const;
    int waitall(int count, MPI_Request *array_of_requests, MPI_Status *array_of_status) const { return MPI_Waitall(count, array_of_requests, array_of_status); }
    int testall(int count, MPI_Request *array_of_requests, int *flag, MPI_Status *array_of_status) const;
    int waitsome(int incount, MPI_Request *array_of_requests,int *outcount, int *array_of_indices, MPI_Status *array_of_status) const;
    int testsome(int incount, MPI_Request *array_of_requests, int *outcount,
                 int *array_of_indices, MPI_Status *array_of_status) const;
    int probe(int source, int tag, MPI_Comm comm, MPI_Status *status) const;
    int Iprobe(int source, int tag, MPI_Comm comm, int *flag, MPI_Status *status) const;
    int cancel(MPI_Request *request) const;
    int testCancelled(MPI_Status *status, int *flag) const;
    int barrier(MPI_Comm comm) const;
    int errorString(int errorcode, char *string, int *resultlen) const;
    int getCount(MPI_Status *status, MPI_Datatype datatype, int *count) const;

    int broadcast(void* buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm) const;
    int allGather(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                  void* recvbuf, int recvcount, MPI_Datatype recvtype,
                  MPI_Comm comm) const;
    int allToAll(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                 void* recvbuf, int recvcount, MPI_Datatype recvtype,
                 MPI_Comm comm) const;
    int allToAllV(void* sendbuf, int* sendcounts, int* senddispls,
                  MPI_Datatype sendtype, void* recvbuf, int* recvcounts,
                  int* recvdispls, MPI_Datatype recvtype, 
                  MPI_Comm comm) const;

    int reduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) const;
    int allReduce(void* sendbuf, void* recvbuf, int count, MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) const;
  };

  class ParaUMesh
  {
  public:
    ParaUMesh(MEDCouplingUMesh *mesh, DataArrayIdType *globalCellIds, DataArrayIdType *globalNodeIds);
    %extend
    {
      DataArrayIdType *getCellIdsLyingOnNodes(const DataArrayIdType *globalNodeIds, bool fullyIn) const
      { 
        MCAuto<DataArrayIdType> ret(self->getCellIdsLyingOnNodes(globalNodeIds,fullyIn));
        return ret.retn();
      }
    }
  };
}

/* This object can be used only if MED_ENABLE_FVM is defined*/
#ifdef MED_ENABLE_FVM
class NonCoincidentDEC : public DEC
{
public:
  NonCoincidentDEC(ProcessorGroup& source, ProcessorGroup& target);
};
#endif

%extend MEDCoupling::ParaMESH
{
  PyObject *getGlobalNumberingCell2() const
  {
    const mcIdType *tmp=self->getGlobalNumberingCell();
    mcIdType size=self->getCellMesh()->getNumberOfCells();
    PyObject *ret=PyList_New(size);
    for(mcIdType i=0;i<size;i++)
      PyList_SetItem(ret,i,PyInt_FromLong(tmp[i])); 
    return ret;
  }

  PyObject *getGlobalNumberingFace2() const
  {
    const mcIdType *tmp=self->getGlobalNumberingFace();
    mcIdType size=self->getFaceMesh()->getNumberOfCells();
    PyObject *ret=PyList_New(size);
    for(mcIdType i=0;i<size;i++)
      PyList_SetItem(ret,i,PyInt_FromLong(tmp[i])); 
    return ret;
  }

  PyObject *getGlobalNumberingNode2() const
  {
    const mcIdType *tmp=self->getGlobalNumberingNode();
    mcIdType size=self->getCellMesh()->getNumberOfNodes();
    PyObject *ret=PyList_New(size);
    for(mcIdType i=0;i<size;i++)
      PyList_SetItem(ret,i,PyInt_FromLong(tmp[i])); 
    return ret;
  }
}
