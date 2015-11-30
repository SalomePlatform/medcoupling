// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#include "OverlapElementLocator.hxx"

#include "CommInterface.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ParaFIELD.hxx"
#include "ParaMESH.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "OverlapInterpolationMatrix.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "DirectedBoundingBox.hxx"
#include "InterpKernelAutoPtr.hxx"

#include <limits>

using namespace std;

namespace ParaMEDMEM 
{ 
  const int OverlapElementLocator::START_TAG_MESH_XCH = 1140;

  OverlapElementLocator::OverlapElementLocator(const ParaFIELD *sourceField, const ParaFIELD *targetField,
                                               const ProcessorGroup& group, double epsAbs)
    : _local_source_field(sourceField),
      _local_target_field(targetField),
      _local_source_mesh(0),
      _local_target_mesh(0),
      _domain_bounding_boxes(0),
      _group(group),
      _epsAbs(epsAbs)
  { 
    if(_local_source_field)
      _local_source_mesh=_local_source_field->getSupport()->getCellMesh();
    if(_local_target_field)
      _local_target_mesh=_local_target_field->getSupport()->getCellMesh();
    _comm=getCommunicator();
    computeBoundingBoxesAndTodoList();
  }

  OverlapElementLocator::~OverlapElementLocator()
  {
    delete [] _domain_bounding_boxes;
  }

  const MPI_Comm *OverlapElementLocator::getCommunicator() const
  {
    const MPIProcessorGroup* group=static_cast<const MPIProcessorGroup*>(&_group);
    return group->getComm();
  }

  void OverlapElementLocator::computeBoundingBoxesAndTodoList()
  {
    CommInterface comm_interface=_group.getCommInterface();
    const MPIProcessorGroup* group=static_cast<const MPIProcessorGroup*> (&_group);
    _local_space_dim=0;
    if(_local_source_mesh)
      _local_space_dim=_local_source_mesh->getSpaceDimension();
    else
      _local_space_dim=_local_target_mesh->getSpaceDimension();
    //
    const MPI_Comm* comm = group->getComm();
    int bbSize=2*2*_local_space_dim;//2 (for source/target) 2 (min/max)
    _domain_bounding_boxes=new double[bbSize*_group.size()];
    INTERP_KERNEL::AutoPtr<double> minmax=new double[bbSize];
    //Format minmax : Xmin_src,Xmax_src,Ymin_src,Ymax_src,Zmin_src,Zmax_src,Xmin_trg,Xmax_trg,Ymin_trg,Ymax_trg,Zmin_trg,Zmax_trg
    if(_local_source_mesh)
      _local_source_mesh->getBoundingBox(minmax);
    else
      {
        for(int i=0;i<_local_space_dim;i++)
          {
            minmax[i*2]=std::numeric_limits<double>::max();
            minmax[i*2+1]=-std::numeric_limits<double>::max();
          }
      }
    if(_local_target_mesh)
      _local_target_mesh->getBoundingBox(minmax+2*_local_space_dim);
    else
      {
        for(int i=0;i<_local_space_dim;i++)
          {
            minmax[i*2+2*_local_space_dim]=std::numeric_limits<double>::max();
            minmax[i*2+1+2*_local_space_dim]=-std::numeric_limits<double>::max();
          }
      }
    comm_interface.allGather(minmax, bbSize, MPI_DOUBLE,
                             _domain_bounding_boxes,bbSize, MPI_DOUBLE, 
                             *comm);
  
    // Computation of all pairs needing an interpolation pairs are duplicated now !
    
    _proc_pairs.clear();//first is source second is target
    _proc_pairs.resize(_group.size());
    for(int i=0;i<_group.size();i++)
      for(int j=0;j<_group.size();j++)
        {
          if(intersectsBoundingBox(i,j))
            {
            _proc_pairs[i].push_back(j);
            }
        }
    // OK now let's assigning as balanced as possible, job to each proc of group
    std::vector< std::vector< ProcCouple > > pairsToBeDonePerProc(_group.size());
    int i=0;
    for(std::vector< std::vector< int > >::const_iterator it1=_proc_pairs.begin();it1!=_proc_pairs.end();it1++,i++)
      for(std::vector< int >::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
        {
          if(pairsToBeDonePerProc[i].size()<=pairsToBeDonePerProc[*it2].size())//it includes the fact that i==*it2
            pairsToBeDonePerProc[i].push_back(ProcCouple(i,*it2));
          else
            pairsToBeDonePerProc[*it2].push_back(ProcCouple(i,*it2));
        }
    //Keeping todo list of current proc. _to_do_list contains a set of pair where at least _group.myRank() appears once.
    //This proc will be in charge to perform interpolation of any of element of '_to_do_list'
    //If _group.myRank()==myPair.first, current proc should fetch target mesh of myPair.second (if different from _group.myRank()).
    //If _group.myRank()==myPair.second, current proc should fetch source mesh of myPair.second.
    
    int myProcId=_group.myRank();
    _to_do_list=pairsToBeDonePerProc[myProcId];

#ifdef DEC_DEBUG
    std::stringstream scout;
    scout << "(" << myProcId << ") my TODO list is: ";
    for (std::vector< ProcCouple >::const_iterator itdbg=_to_do_list.begin(); itdbg!=_to_do_list.end(); itdbg++)
      scout << "(" << (*itdbg).first << "," << (*itdbg).second << ")";
    std::cout << scout.str() << "\n";
#endif

    // Feeding now '_procs_to_send*'. A same id can appears twice. The second parameter in pair means what
    // to send true=source, false=target
    _procs_to_send_mesh.clear();
    _procs_to_send_field.clear();
    for(int i=_group.size()-1;i>=0;i--)
      {
        const std::vector< ProcCouple >& anRemoteProcToDoList=pairsToBeDonePerProc[i];
        for(std::vector< ProcCouple >::const_iterator it=anRemoteProcToDoList.begin();it!=anRemoteProcToDoList.end();it++)
          {
            if((*it).first==myProcId)
              {
                if(i!=myProcId)
                  _procs_to_send_mesh.push_back(Proc_SrcOrTgt(i,true));
                _procs_to_send_field.push_back((*it).second);
              }
            if((*it).second==myProcId)
              if(i!=myProcId)
                _procs_to_send_mesh.push_back(Proc_SrcOrTgt(i,false));
          }
      }
  }

  /*!
   * The aim of this method is to perform the communication to get data corresponding to '_to_do_list' attribute.
   * The principle is the following : if proc n1 and n2 need to perform a cross sending with n1<n2, then n1 will send first and receive then.
   */
  void OverlapElementLocator::exchangeMeshes(OverlapInterpolationMatrix& matrix)
  {
    int myProcId=_group.myRank();
    //starting to receive every procs whose id is lower than myProcId.
    std::vector< ProcCouple > toDoListForFetchRemaining;
    for(std::vector< ProcCouple >::const_iterator it=_to_do_list.begin();it!=_to_do_list.end();it++)
      {
        int first = (*it).first, second = (*it).second;
        if(first!=second)
          {
            if(first==myProcId)
              {
                if(second<myProcId)
                  receiveRemoteMeshFrom(second,false);
                else
                  toDoListForFetchRemaining.push_back(ProcCouple(first,second));
              }
            else
              {//(*it).second==myProcId
                if(first<myProcId)
                  receiveRemoteMeshFrom(first,true);
                else
                  toDoListForFetchRemaining.push_back(ProcCouple(first,second));
              }
          }
      }
    //sending source or target mesh to remote procs
    for(std::vector< Proc_SrcOrTgt >::const_iterator it2=_procs_to_send_mesh.begin();it2!=_procs_to_send_mesh.end();it2++)
      sendLocalMeshTo((*it2).first,(*it2).second,matrix);
    //fetching remaining meshes
    for(std::vector< ProcCouple >::const_iterator it=toDoListForFetchRemaining.begin();it!=toDoListForFetchRemaining.end();it++)
      {
        if((*it).first!=(*it).second)
          {
            if((*it).first==myProcId)
              receiveRemoteMeshFrom((*it).second,false);
            else//(*it).second==myProcId
              receiveRemoteMeshFrom((*it).first,true);
          }
      }
  }
  
  std::string OverlapElementLocator::getSourceMethod() const
  {
    return _local_source_field->getField()->getDiscretization()->getStringRepr();
  }

  std::string OverlapElementLocator::getTargetMethod() const
  {
    return _local_target_field->getField()->getDiscretization()->getStringRepr();
  }

  const MEDCouplingPointSet *OverlapElementLocator::getSourceMesh(int procId) const
  {
    int myProcId=_group.myRank();
    if(myProcId==procId)
      return _local_source_mesh;
    Proc_SrcOrTgt p(procId,true);
    std::map<Proc_SrcOrTgt, AutoMCPointSet >::const_iterator it=_remote_meshes.find(p);
    return (*it).second;
  }

  const DataArrayInt *OverlapElementLocator::getSourceIds(int procId) const
  {
    int myProcId=_group.myRank();
    if(myProcId==procId)
      return 0;
    Proc_SrcOrTgt p(procId,true);
    std::map<Proc_SrcOrTgt, AutoDAInt >::const_iterator it=_remote_elems.find(p);
    return (*it).second;
  }

  const MEDCouplingPointSet *OverlapElementLocator::getTargetMesh(int procId) const
  {
    int myProcId=_group.myRank();
    if(myProcId==procId)
      return _local_target_mesh;
    Proc_SrcOrTgt p(procId,false);
    std::map<Proc_SrcOrTgt, AutoMCPointSet >::const_iterator it=_remote_meshes.find(p);
    return (*it).second;
  }

  const DataArrayInt *OverlapElementLocator::getTargetIds(int procId) const
  {
    int myProcId=_group.myRank();
    if(myProcId==procId)
      return 0;
    Proc_SrcOrTgt p(procId,false);
    std::map<Proc_SrcOrTgt, AutoDAInt >::const_iterator it=_remote_elems.find(p);
    return (*it).second;
  }

  bool OverlapElementLocator::isInMyTodoList(int i, int j) const
  {
    ProcCouple cpl = std::make_pair(i, j);
    return std::find(_to_do_list.begin(), _to_do_list.end(), cpl)!=_to_do_list.end();
  }

  bool OverlapElementLocator::intersectsBoundingBox(int isource, int itarget) const
  {
    const double *source_bb=_domain_bounding_boxes+isource*2*2*_local_space_dim;
    const double *target_bb=_domain_bounding_boxes+itarget*2*2*_local_space_dim+2*_local_space_dim;

    for (int idim=0; idim < _local_space_dim; idim++)
      {
        bool intersects = (target_bb[idim*2]<source_bb[idim*2+1]+_epsAbs)
          && (source_bb[idim*2]<target_bb[idim*2+1]+_epsAbs);
        if (!intersects)
          return false; 
      }
    return true;
  }

  /*!
   * This methods sends (part of) local source if 'sourceOrTarget'==True to proc 'procId'.
   * This methods sends (part of) local target if 'sourceOrTarget'==False to proc 'procId'.
   *
   * This method prepares the matrix too, for matrix assembling and future matrix-vector computation.
   */
  void OverlapElementLocator::sendLocalMeshTo(int procId, bool sourceOrTarget, OverlapInterpolationMatrix& matrix) const
  {
   //int myProcId=_group.myRank();
   const double *distant_bb=0;
   MEDCouplingPointSet *local_mesh=0;
   const ParaFIELD *field=0;
   if(sourceOrTarget)//source for local mesh but target for distant mesh
     {
       distant_bb=_domain_bounding_boxes+procId*2*2*_local_space_dim+2*_local_space_dim;
       local_mesh=_local_source_mesh;
       field=_local_source_field;
     }
   else//target for local but source for distant
     {
       distant_bb=_domain_bounding_boxes+procId*2*2*_local_space_dim;
       local_mesh=_local_target_mesh;
       field=_local_target_field;
     }
   AutoDAInt elems=local_mesh->getCellsInBoundingBox(distant_bb,getBoundingBoxAdjustment());
   DataArrayInt *old2new_map;
   MEDCouplingPointSet *send_mesh=static_cast<MEDCouplingPointSet *>(field->getField()->buildSubMeshData(elems->begin(),elems->end(),old2new_map));
   if(sourceOrTarget)
     matrix.keepTracksOfSourceIds(procId,old2new_map);
   else
     matrix.keepTracksOfTargetIds(procId,old2new_map);
   sendMesh(procId,send_mesh,old2new_map);
   send_mesh->decrRef();
   old2new_map->decrRef();
  }

  /*!
   * This method recieves source remote mesh on proc 'procId' if sourceOrTarget==True
   * This method recieves target remote mesh on proc 'procId' if sourceOrTarget==False
   */
  void OverlapElementLocator::receiveRemoteMeshFrom(int procId, bool sourceOrTarget)
  {
    DataArrayInt *old2new_map=0;
    MEDCouplingPointSet *m=0;
    receiveMesh(procId,m,old2new_map);
    Proc_SrcOrTgt p(procId,sourceOrTarget);
    _remote_meshes[p]=m;
    _remote_elems[p]=old2new_map;
  }

  void OverlapElementLocator::sendMesh(int procId, const MEDCouplingPointSet *mesh, const DataArrayInt *idsToSend) const
  {
    CommInterface comInterface=_group.getCommInterface();

    // First stage : exchanging sizes
    vector<double> tinyInfoLocalD;//tinyInfoLocalD not used for the moment
    vector<int> tinyInfoLocal;
    vector<string> tinyInfoLocalS;
    mesh->getTinySerializationInformation(tinyInfoLocalD,tinyInfoLocal,tinyInfoLocalS);
    const MPI_Comm *comm=getCommunicator();
    //
    int lgth[2];
    lgth[0]=tinyInfoLocal.size();
    lgth[1]=idsToSend->getNbOfElems();
    comInterface.send(&lgth,2,MPI_INT,procId,START_TAG_MESH_XCH,*_comm);
    comInterface.send(&tinyInfoLocal[0],tinyInfoLocal.size(),MPI_INT,procId,START_TAG_MESH_XCH+1,*comm);
    //
    DataArrayInt *v1Local=0;
    DataArrayDouble *v2Local=0;
    mesh->serialize(v1Local,v2Local);
    comInterface.send(v1Local->getPointer(),v1Local->getNbOfElems(),MPI_INT,procId,START_TAG_MESH_XCH+2,*comm);
    comInterface.send(v2Local->getPointer(),v2Local->getNbOfElems(),MPI_DOUBLE,procId,START_TAG_MESH_XCH+3,*comm);
    //finished for mesh, ids now
    comInterface.send(const_cast<int *>(idsToSend->getConstPointer()),lgth[1],MPI_INT,procId,START_TAG_MESH_XCH+4,*comm);
    //
    v1Local->decrRef();
    v2Local->decrRef();
  }

  void OverlapElementLocator::receiveMesh(int procId, MEDCouplingPointSet* &mesh, DataArrayInt *&ids) const
  {
    int lgth[2];
    MPI_Status status;
    const MPI_Comm *comm=getCommunicator();
    CommInterface comInterface=_group.getCommInterface();
    comInterface.recv(lgth,2,MPI_INT,procId,START_TAG_MESH_XCH,*_comm,&status);
    std::vector<int> tinyInfoDistant(lgth[0]);
    ids=DataArrayInt::New();
    ids->alloc(lgth[1],1);
    comInterface.recv(&tinyInfoDistant[0],lgth[0],MPI_INT,procId,START_TAG_MESH_XCH+1,*comm,&status);
    mesh=MEDCouplingPointSet::BuildInstanceFromMeshType((MEDCouplingMeshType)tinyInfoDistant[0]);
    std::vector<std::string> unusedTinyDistantSts;
    vector<double> tinyInfoDistantD(1);//tinyInfoDistantD not used for the moment
    DataArrayInt *v1Distant=DataArrayInt::New();
    DataArrayDouble *v2Distant=DataArrayDouble::New();
    mesh->resizeForUnserialization(tinyInfoDistant,v1Distant,v2Distant,unusedTinyDistantSts);
    comInterface.recv(v1Distant->getPointer(),v1Distant->getNbOfElems(),MPI_INT,procId,START_TAG_MESH_XCH+2,*comm,&status);
    comInterface.recv(v2Distant->getPointer(),v2Distant->getNbOfElems(),MPI_DOUBLE,procId,START_TAG_MESH_XCH+3,*comm,&status);
    mesh->unserialization(tinyInfoDistantD,tinyInfoDistant,v1Distant,v2Distant,unusedTinyDistantSts);
    //finished for mesh, ids now
    comInterface.recv(ids->getPointer(),lgth[1],MPI_INT,procId,1144,*comm,&status);
    //
    v1Distant->decrRef();
    v2Distant->decrRef();
  }
}
