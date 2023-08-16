// Copyright (C) 2017-2023  CEA, EDF
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
#include "OverlapDEC.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "ICoCoMEDDoubleField.hxx"
#include "ICoCoMEDIntField.hxx"
#include "ComponentTopology.hxx"
#include "ParaUMesh.hxx"
#include "ParaSkyLineArray.hxx"
#include "ParaDataArray.hxx"

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
%include "StructuredCoincidentDEC.hxx"

%newobject MEDCoupling::ParaUMesh::New;
%newobject MEDCoupling::ParaUMesh::getMesh;
%newobject MEDCoupling::ParaUMesh::getGlobalCellIds;
%newobject MEDCoupling::ParaUMesh::getGlobalNodeIds;
%newobject MEDCoupling::ParaUMesh::getCellIdsLyingOnNodes;
%newobject MEDCoupling::ParaUMesh::redistributeCells;
%newobject MEDCoupling::ParaUMesh::redistributeCellField;
%newobject MEDCoupling::ParaUMesh::redistributeNodeField;
%newobject MEDCoupling::ParaDataArrayInt32::New;
%newobject MEDCoupling::ParaDataArrayInt32::buildComplement;
%newobject MEDCoupling::ParaDataArrayInt64::New;
%newobject MEDCoupling::ParaDataArrayInt64::buildComplement;
%newobject MEDCoupling::ParaSkyLineArray::New;
%newobject MEDCoupling::ParaSkyLineArray::equiRedistribute;
%newobject MEDCoupling::ParaSkyLineArray::getSkyLineArray;
%newobject MEDCoupling::ParaSkyLineArray::getGlobalIdsArray;

%newobject MEDCoupling::InterpKernelDEC::_NewWithPG_internal;
%newobject MEDCoupling::InterpKernelDEC::_NewWithComm_internal;
%newobject MEDCoupling::InterpKernelDEC::retrieveNonFetchedIds;
%newobject MEDCoupling::OverlapDEC::_NewWithComm_internal;

%feature("unref") ParaSkyLineArray "$this->decrRef();"
%feature("unref") ParaUMesh "$this->decrRef();"
%feature("unref") ParaDataArrayInt32 "$this->decrRef();"
%feature("unref") ParaDataArrayInt64 "$this->decrRef();"

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

  class ParaDataArray : public RefCountObject
  {
  };

  class ParaDataArrayInt32 : public ParaDataArray
  {
  public:
    static ParaDataArrayInt32 *New(DataArrayInt32 *seqDa);
    DataArrayIdType *buildComplement(int nbOfElems) const;
    %extend
    {
      ParaDataArrayInt32(DataArrayInt32 *seqDa)
      {
        return ParaDataArrayInt32::New(seqDa);
      }
    }
  };

  class ParaDataArrayInt64 : public ParaDataArray
  {
  public:
    static ParaDataArrayInt64 *New(DataArrayInt64 *seqDa);
    DataArrayIdType *buildComplement(long nbOfElems) const;
    %extend
    {
      ParaDataArrayInt64(DataArrayInt64 *seqDa)
      {
        return ParaDataArrayInt64::New(seqDa);
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

  /* This object can be used only if MED_ENABLE_FVM is defined*/
  #ifdef MED_ENABLE_FVM
  class NonCoincidentDEC : public DEC
  {
  public:
    NonCoincidentDEC(ProcessorGroup& source, ProcessorGroup& target);
  };
  #endif

  class InterpKernelDEC : public DisjointDEC, public INTERP_KERNEL::InterpolationOptions
  {
    public:
      InterpKernelDEC();
      InterpKernelDEC(ProcessorGroup& source_group, ProcessorGroup& target_group);
      InterpKernelDEC(const std::set<int>& src_ids, const std::set<int>& trg_ids); // hide last optional parameter!
      virtual ~InterpKernelDEC();
      void release();

      void synchronize();
      void synchronizeWithDefaultValue(double val);
      void recvData();
      void recvData(double time);
      void sendData();
      void sendData(double time , double deltatime);
      void prepareSourceDE();
      void prepareTargetDE();

      %extend {
        // Provides a direct ctor for which the communicator can be passed with "MPI._addressof(the_com)":
        InterpKernelDEC(const std::set<int>& src_ids, const std::set<int>& trg_ids, long long comm_ptr)
        {
            return new InterpKernelDEC(src_ids, trg_ids, *((MPI_Comm*)comm_ptr));
        }

        // This one should really not be called directly by the user since it still has an interface with a pointer to MPI_Comm 
        // which Swig doesn't handle nicely.
        // It is just here to provide a constructor taking a **pointer** to a comm - See pythoncode below.
        static InterpKernelDEC* _NewWithPG_internal(ProcessorGroup& source_group, ProcessorGroup& target_group)
        {
          return new InterpKernelDEC(source_group,target_group);
        }

        static InterpKernelDEC* _NewWithComm_internal(const std::set<int>& src_ids, const std::set<int>& trg_ids, long long another_comm)
        {
          return new InterpKernelDEC(src_ids,trg_ids, *(MPI_Comm*)another_comm); // I know, ugly cast ...
        }

        DataArrayIdType *retrieveNonFetchedIds() const
        {
          MCAuto<DataArrayIdType> ret = self->retrieveNonFetchedIds();
          return ret.retn();
        }
      }
  };

  class OverlapDEC : public DEC, public INTERP_KERNEL::InterpolationOptions
  {
      public:
        OverlapDEC(const std::set<int>& procIds);  // hide optional param comm
        virtual ~OverlapDEC();
        void release();

        void sendRecvData(bool way=true);
        void sendData();
        void recvData();
        void synchronize();
        void attachSourceLocalField(ParaFIELD *field, bool ownPt=false);
        void attachTargetLocalField(ParaFIELD *field, bool ownPt=false);
        void attachSourceLocalField(MEDCouplingFieldDouble *field);
        void attachTargetLocalField(MEDCouplingFieldDouble *field);
        void attachSourceLocalField(ICoCo::MEDDoubleField *field);
        void attachTargetLocalField(ICoCo::MEDDoubleField *field);
        ProcessorGroup *getGroup();
        bool isInGroup() const;

        void setDefaultValue(double val);
        void setWorkSharingAlgo(int method);

        void debugPrintWorkSharing(std::ostream & ostr) const;

        %extend {
          OverlapDEC(const std::set<int>& ids, long long comm_ptr)
          {
             return new OverlapDEC(ids, *((MPI_Comm*)comm_ptr));
          }

          // This one should really not be called directly by the user since it still has an interface with a pointer to MPI_Comm 
          // which Swig doesn't handle nicely.
          // It is just here to provide a constructor taking a **pointer** to a comm - See pythoncode below.
          static OverlapDEC* _NewWithComm_internal(const std::set<int>& ids, long long another_comm)
          {
            return new OverlapDEC(ids, *(MPI_Comm*)another_comm); // I know, ugly cast ...
          }
        }
   };

} // end namespace MEDCoupling

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

%pythoncode %{
if MEDCouplingUse64BitIDs():
  ParaDataArrayInt = ParaDataArrayInt64
else:
  ParaDataArrayInt = ParaDataArrayInt32
%}

%pythoncode %{

# And here we use mpi4py ability to provide its internal (C++) pointer to the communicator:
# NB: doing a proper typemap from MPI_Comm from Python to C++ requires the inclusion of mpi4py headers and .i file ... an extra dependency ...
def _IKDEC_WithComm_internal(src_procs, tgt_procs, mpicomm=None):
    from mpi4py import MPI
    # Check iterable:
    try:
        s, t = [el for el in src_procs], [el for el in tgt_procs]
    except:
        s, t = None, None
    msg =  "InterpKernelDEC: invalid type in ctor arguments! Possible signatures are:\n"
    msg += "   - InterpKernelDEC(ProcessorGroup, ProcessorGroup)\n"
    msg += "   - InterpKernelDEC(<iterable>, <iterable>)\n"
    msg += "   - InterpKernelDEC(<iterable>, <iterable>, MPI_Comm*) : WARNING here the address of the communicator should be passed with MPI._addressof(the_com)\n"
    msg += "   - InterpKernelDEC.New(ProcessorGroup, ProcessorGroup)\n"
    msg += "   - InterpKernelDEC.New(<iterable>, <iterable>)\n"
    msg += "   - InterpKernelDEC.New(<iterable>, <iterable>, MPI_Comm)\n"
    if mpicomm is None:
        if isinstance(src_procs, ProcessorGroup) and isinstance(tgt_procs, ProcessorGroup):
            return InterpKernelDEC._NewWithPG_internal(src_procs, tgt_procs)
        elif not s is None:  # iterable
            return InterpKernelDEC._NewWithComm_internal(s, t, MPI._addressof(MPI.COMM_WORLD))
        else:
            raise InterpKernelException(msg)
    else:
        if s is None: raise InterpKernelException(msg)  # must be iterable
        return InterpKernelDEC._NewWithComm_internal(s, t, MPI._addressof(mpicomm))

def _ODEC_WithComm_internal(procs, mpicomm=None):
    from mpi4py import MPI
    # Check iterable:
    try:
        g = [el for el in procs]
    except:
        msg =  "OverlapDEC: invalid type in ctor arguments! Possible signatures are:\n"
        msg += "   - OverlapDEC.New(<iterable>)\n"
        msg += "   - OverlapDEC.New(<iterable>, MPI_Comm)\n"
        msg += "   - OverlapDEC(<iterable>)\n"
        msg += "   - OverlapDEC(<iterable>, MPI_Comm*) : WARNING here the address of the communicator should be passed with MPI._addressof(the_com)\n"
        raise InterpKernelException(msg)
    if mpicomm is None:
        return OverlapDEC(g)
    else:
        return OverlapDEC._NewWithComm_internal(g, MPI._addressof(mpicomm))

InterpKernelDEC.New = _IKDEC_WithComm_internal
OverlapDEC.New = _ODEC_WithComm_internal

%}
