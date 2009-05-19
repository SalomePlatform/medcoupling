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
#include "GlobalizerMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "CommInterface.hxx"

using namespace std;

namespace ParaMEDMEM
{
  GlobalizerMesh::GlobalizerMesh(const MPI_Comm *comm, MEDCouplingFieldDouble *localField):_comm(comm),_local_field(localField)
  {
    if(_local_field)
      _local_field->incrRef();
  }
  
  GlobalizerMesh::~GlobalizerMesh()
  {
    if(_local_field)
      _local_field->decrRef();
  }

  NatureOfField GlobalizerMesh::getLocalNature() const
  {
    return _local_field->getNature();
  }

  GlobalizerMeshWorkingSide::GlobalizerMeshWorkingSide(const MPI_Comm *comm, MEDCouplingFieldDouble *localField,
                                                       const std::string& distantMeth, const std::vector<int>& lazyProcs):GlobalizerMesh(comm,localField),_distant_method(distantMeth),_lazy_procs(lazyProcs)
  {
  }

  GlobalizerMeshWorkingSide::~GlobalizerMeshWorkingSide()
  {
  }

  const std::vector<int>& GlobalizerMeshWorkingSide::getProcIdsInInteraction() const
  {
    return _lazy_procs;
  }

  /*!
   * connected with GlobalizerMeshLazySide::recvFromWorkingSide
   */
  void GlobalizerMeshWorkingSide::sendSumToLazySide(const std::vector< std::vector<int> >& distantLocEltIds, const std::vector< std::vector<double> >& partialSumRelToDistantIds)
  {
    int procId=0;
    CommInterface comm;
    for(vector<int>::const_iterator iter=_lazy_procs.begin();iter!=_lazy_procs.end();iter++,procId++)
      {
        const vector<int>& eltIds=distantLocEltIds[procId];
        const vector<double>& valued=partialSumRelToDistantIds[procId];
        int lgth=eltIds.size();
        comm.send(&lgth,1,MPI_INT,*iter,1114,*_comm);
        comm.send((void *)&eltIds[0],lgth,MPI_INT,*iter,1115,*_comm);
        comm.send((void *)&valued[0],lgth,MPI_DOUBLE,*iter,1116,*_comm);
      }
  }

  /*!
   * connected with GlobalizerMeshLazySide::sendToWorkingSide
   */
  void GlobalizerMeshWorkingSide::recvSumFromLazySide(std::vector< std::vector<double> >& globalSumRelToDistantIds)
  {
    int procId=0;
    CommInterface comm;
    MPI_Status status;
    for(vector<int>::const_iterator iter=_lazy_procs.begin();iter!=_lazy_procs.end();iter++,procId++)
      {
        std::vector<double>& vec=globalSumRelToDistantIds[procId];
        comm.recv(&vec[0],vec.size(),MPI_DOUBLE,*iter,1117,*_comm,&status);
      }
  }

  GlobalizerMeshLazySide::GlobalizerMeshLazySide(const MPI_Comm *comm, MEDCouplingFieldDouble *localField, const std::vector<int>& computeProcs):GlobalizerMesh(comm,localField),_compute_procs(computeProcs)
  {
  }

  GlobalizerMeshLazySide::~GlobalizerMeshLazySide()
  {
  }

  /*!
   * connected with GlobalizerMeshWorkingSide::sendSumToLazySide
   */
  void GlobalizerMeshLazySide::recvFromWorkingSide()
  {
    _values_added.resize(_local_field->getNumberOfTuples());
    int procId=0;
    CommInterface comm;
    _ids_per_working_proc.resize(_compute_procs.size());
    MPI_Status status;
    for(vector<int>::const_iterator iter=_compute_procs.begin();iter!=_compute_procs.end();iter++,procId++)
      {
        int lgth;
        comm.recv(&lgth,1,MPI_INT,*iter,1114,*_comm,&status);
        vector<int>& ids=_ids_per_working_proc[procId];
        ids.resize(lgth);
        vector<double> values(lgth);
        comm.recv(&ids[0],lgth,MPI_INT,*iter,1115,*_comm,&status);
        comm.recv(&values[0],lgth,MPI_DOUBLE,*iter,1116,*_comm,&status);
        for(int i=0;i<lgth;i++)
          _values_added[ids[i]]+=values[i];
      }
  }

  /*!
   * connected with GlobalizerMeshWorkingSide::recvSumFromLazySide
   */
  void GlobalizerMeshLazySide::sendToWorkingSide()
  {
    int procId=0;
    CommInterface comm;
    for(vector<int>::const_iterator iter=_compute_procs.begin();iter!=_compute_procs.end();iter++,procId++)
      {
        vector<int>& ids=_ids_per_working_proc[procId];
        vector<double> valsToSend(ids.size());
        vector<double>::iterator iter3=valsToSend.begin();
        for(vector<int>::const_iterator iter2=ids.begin();iter2!=ids.end();iter2++,iter3++)
          *iter3=_values_added[*iter2];
        comm.send(&valsToSend[0],ids.size(),MPI_DOUBLE,*iter,1117,*_comm);
        ids.clear();
      }
    _ids_per_working_proc.clear();
  }
}

