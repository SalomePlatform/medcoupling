// Copyright (C) 2017-2020  CEA/DEN, EDF R&D
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
#include "ParaSkyLineArray.hxx"

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

%newobject MEDCoupling::ParaUMesh::New;
%newobject MEDCoupling::ParaUMesh::getMesh;
%newobject MEDCoupling::ParaUMesh::getGlobalCellIds;
%newobject MEDCoupling::ParaUMesh::getGlobalNodeIds;
%newobject MEDCoupling::ParaUMesh::getCellIdsLyingOnNodes;
%newobject MEDCoupling::ParaUMesh::redistributeCells;
%newobject MEDCoupling::ParaUMesh::redistributeCellField;
%newobject MEDCoupling::ParaUMesh::redistributeNodeField;
%newobject MEDCoupling::ParaSkyLineArray::New;
%newobject MEDCoupling::ParaSkyLineArray::equiRedistribute;
%newobject MEDCoupling::ParaSkyLineArray::getSkyLineArray;
%newobject MEDCoupling::ParaSkyLineArray::getGlobalIdsArray;

%feature("unref") ParaSkyLineArray "$this->decrRef();"
%feature("unref") ParaUMesh "$this->decrRef();"

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
    %extend
    {
      PyObject *allGatherArrays(const DataArrayIdType *array) const
      {
        std::vector< MCAuto<DataArrayIdType> > ret;
        self->allGatherArrays(MPI_COMM_WORLD,array,ret);
        return convertFromVectorAutoObjToPyObj<DataArrayIdType>(ret,SWIGTITraits<mcIdType>::TI);
      }

      PyObject *allToAllArrays(PyObject *arrays) const
      {
        std::vector< DataArrayIdType * > arraysIn;
        std::vector< MCAuto<DataArrayIdType> > arrayOut;
        convertFromPyObjVectorOfObj<MEDCoupling::DataArrayIdType*>(arrays,SWIGTITraits<mcIdType>::TI,"DataArrayIdType",arraysIn);
        std::vector< MCAuto<DataArrayIdType> > arraysIn2(FromVecToVecAuto<DataArrayIdType>(arraysIn));
        self->allToAllArrays(MPI_COMM_WORLD,arraysIn2,arrayOut);
        return convertFromVectorAutoObjToPyObj<DataArrayIdType>(arrayOut,SWIGTITraits<mcIdType>::TI);
      }
    }
  };

  class ParaUMesh : public RefCountObject
  {
  public:
    static ParaUMesh *New(MEDCouplingUMesh *mesh, DataArrayIdType *globalCellIds, DataArrayIdType *globalNodeIds);
    ParaUMesh *redistributeCells(const DataArrayIdType *globalCellIds) const;
    DataArrayIdType *redistributeCellField(const DataArrayIdType *globalCellIds, const DataArrayIdType *fieldValueToRed) const;
    DataArrayDouble *redistributeCellField(const DataArrayIdType *globalCellIds, const DataArrayDouble *fieldValueToRed) const;
    DataArrayIdType *redistributeNodeField(const DataArrayIdType *globalCellIds, const DataArrayIdType *fieldValueToRed) const;
    DataArrayDouble *redistributeNodeField(const DataArrayIdType *globalCellIds, const DataArrayDouble *fieldValueToRed) const;
    %extend
    {
      ParaUMesh(MEDCouplingUMesh *mesh, DataArrayIdType *globalCellIds, DataArrayIdType *globalNodeIds)
      {
        return ParaUMesh::New(mesh,globalCellIds,globalNodeIds);
      }

      MEDCouplingUMesh *getMesh()
      {
        MEDCouplingUMesh *ret(self->getMesh());
        if(ret) ret->incrRef();
        return ret;
      }

      DataArrayIdType *getGlobalCellIds()
      {
        DataArrayIdType *ret(self->getGlobalCellIds());
        if(ret) ret->incrRef();
        return ret;
      }

      DataArrayIdType *getGlobalNodeIds()
      {
        DataArrayIdType *ret(self->getGlobalNodeIds());
        if(ret) ret->incrRef();
        return ret;
      }

      DataArrayIdType *getCellIdsLyingOnNodes(const DataArrayIdType *globalNodeIds, bool fullyIn) const
      { 
        MCAuto<DataArrayIdType> ret(self->getCellIdsLyingOnNodes(globalNodeIds,fullyIn));
        return ret.retn();
      }
    }
  };

  class ParaSkyLineArray : public RefCountObject
  {
  public:
    static ParaSkyLineArray *New(MEDCouplingSkyLineArray *ska, DataArrayIdType *globalIds);
    %extend
    {
      ParaSkyLineArray(MEDCouplingSkyLineArray *ska, DataArrayIdType *globalIds)
      {
        return ParaSkyLineArray::New(ska,globalIds);
      }

      ParaSkyLineArray *equiRedistribute(mcIdType nbOfEntities) const
      {
        MCAuto<ParaSkyLineArray> ret(self->equiRedistribute(nbOfEntities));
        return ret.retn();
      }

      MEDCouplingSkyLineArray *getSkyLineArray() const
      {
        MEDCouplingSkyLineArray *ret(self->getSkyLineArray());
        if(ret)
          ret->incrRef();
        return ret;
      }
      
      DataArrayIdType *getGlobalIdsArray() const
      {
        DataArrayIdType *ret(self->getGlobalIdsArray());
        if(ret)
          ret->incrRef();
        return ret;
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
