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

#include <mpi.h>
#include "CommInterface.hxx"
#include "ElementLocator.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ParaFIELD.hxx"
#include "ParaMESH.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MCAuto.hxx"
#include "DirectedBoundingBox.hxx"

#include <map>
#include <set>
#include <limits>

using namespace std;

//#define USE_DIRECTED_BB

namespace MEDCoupling 
{ 
  ElementLocator::ElementLocator(const ParaFIELD& sourceField,
                                 const ProcessorGroup& distant_group,
                                 const ProcessorGroup& local_group)
    : _local_para_field(sourceField),
      _local_cell_mesh(sourceField.getSupport()->getCellMesh()),
      _local_face_mesh(sourceField.getSupport()->getFaceMesh()),
      _distant_group(distant_group),
      _local_group(local_group)
  { 
    _union_group = _local_group.fuse(distant_group);
    _computeBoundingBoxes();
    _comm=getCommunicator();
  }

  ElementLocator::~ElementLocator()
  {
    delete _union_group;
    delete [] _domain_bounding_boxes;
  }

  const MPI_Comm *ElementLocator::getCommunicator() const
  {
    MPIProcessorGroup* group=static_cast<MPIProcessorGroup*> (_union_group);
    return group->getComm();
  }

  NatureOfField ElementLocator::getLocalNature() const
  {
    return _local_para_field.getField()->getNature();
  }


  /*! Procedure for exchanging a mesh between a distant proc and a local processor
   \param idistantrank  proc id on distant group
   \param distant_mesh on return , points to a local reconstruction of
          the distant mesh
   \param distant_ids on return, contains a vector defining a correspondence
          between the distant ids and the ids of the local reconstruction
  */
  void ElementLocator::exchangeMesh(int idistantrank,
                                    MEDCouplingPointSet*& distant_mesh,
                                    int*& distant_ids)
  {
    int rank = _union_group->translateRank(&_distant_group,idistantrank);

    if (find(_distant_proc_ids.begin(), _distant_proc_ids.end(),rank)==_distant_proc_ids.end())
      return;
   
    MCAuto<DataArrayInt> elems;
#ifdef USE_DIRECTED_BB
    INTERP_KERNEL::DirectedBoundingBox dbb;
    double* distant_bb = _domain_bounding_boxes+rank*dbb.dataSize(_local_cell_mesh_space_dim);
    dbb.setData(distant_bb);
    elems=_local_cell_mesh->getCellsInBoundingBox(dbb,getBoundingBoxAdjustment());
#else
    double* distant_bb = _domain_bounding_boxes+rank*2*_local_cell_mesh_space_dim;
    elems=_local_cell_mesh->getCellsInBoundingBox(distant_bb,getBoundingBoxAdjustment());
#endif
    
    DataArrayInt *distant_ids_send;
    MEDCouplingPointSet *send_mesh = (MEDCouplingPointSet *)_local_para_field.getField()->buildSubMeshData(elems->begin(),elems->end(),distant_ids_send);
    _exchangeMesh(send_mesh, distant_mesh, idistantrank, distant_ids_send, distant_ids);
    distant_ids_send->decrRef();
    
    if(send_mesh)
      send_mesh->decrRef();
  }

  void ElementLocator::exchangeMethod(const std::string& sourceMeth, int idistantrank, std::string& targetMeth)
  {
    CommInterface comm_interface=_union_group->getCommInterface();
    MPIProcessorGroup* group=static_cast<MPIProcessorGroup*> (_union_group);
    const MPI_Comm* comm=(group->getComm());
    MPI_Status status;
    // it must be converted to union numbering before communication
    int idistRankInUnion = group->translateRank(&_distant_group,idistantrank);
    char *recv_buffer=new char[4];
    std::vector<char> send_buffer(4);
    std::copy(sourceMeth.begin(),sourceMeth.end(),send_buffer.begin());
    comm_interface.sendRecv(&send_buffer[0], 4, MPI_CHAR,idistRankInUnion, 1112,
                            recv_buffer, 4, MPI_CHAR,idistRankInUnion, 1112,
                            *comm, &status);
    targetMeth=recv_buffer;
    delete [] recv_buffer;
  }



  /*!
   Compute bounding boxes
  */
  void ElementLocator::_computeBoundingBoxes()
  {
    CommInterface comm_interface =_union_group->getCommInterface();
    MPIProcessorGroup* group=static_cast<MPIProcessorGroup*> (_union_group);
    const MPI_Comm* comm = group->getComm();
    _local_cell_mesh_space_dim = -1;
    if(_local_cell_mesh->getMeshDimension() != -1)
      _local_cell_mesh_space_dim=_local_cell_mesh->getSpaceDimension();
    int *spaceDimForAll=new int[_union_group->size()];
    comm_interface.allGather(&_local_cell_mesh_space_dim, 1, MPI_INT,
                             spaceDimForAll,1, MPI_INT, 
                             *comm);
    _local_cell_mesh_space_dim=*std::max_element(spaceDimForAll,spaceDimForAll+_union_group->size());
    _is_m1d_corr=((*std::min_element(spaceDimForAll,spaceDimForAll+_union_group->size()))==-1);
    for(int i=0;i<_union_group->size();i++)
      if(spaceDimForAll[i]!=_local_cell_mesh_space_dim && spaceDimForAll[i]!=-1)
        throw INTERP_KERNEL::Exception("Spacedim not matches !");
    delete [] spaceDimForAll;
#ifdef USE_DIRECTED_BB
    INTERP_KERNEL::DirectedBoundingBox dbb;
    int bbSize = dbb.dataSize(_local_cell_mesh_space_dim);
    _domain_bounding_boxes = new double[bbSize*_union_group->size()];
    if(_local_cell_mesh->getMeshDimension() != -1)
      dbb = INTERP_KERNEL::DirectedBoundingBox(_local_cell_mesh->getCoords()->getPointer(),
                                               _local_cell_mesh->getNumberOfNodes(),
                                               _local_cell_mesh_space_dim);
    std::vector<double> dbbData = dbb.getData();
    if ( dbbData.size() < bbSize ) dbbData.resize(bbSize,0);
    double * minmax= &dbbData[0];
#else
    int bbSize = 2*_local_cell_mesh_space_dim;
    _domain_bounding_boxes = new double[bbSize*_union_group->size()];
    double * minmax=new double [bbSize];
    if(_local_cell_mesh->getMeshDimension() != -1)
      _local_cell_mesh->getBoundingBox(minmax);
    else
      for(int i=0;i<_local_cell_mesh_space_dim;i++)
        {
          minmax[i*2]=-std::numeric_limits<double>::max();
          minmax[i*2+1]=std::numeric_limits<double>::max();
        }
#endif

    comm_interface.allGather(minmax, bbSize, MPI_DOUBLE,
                             _domain_bounding_boxes,bbSize, MPI_DOUBLE, 
                             *comm);
  
    for (int i=0; i< _distant_group.size(); i++)
      {
        int rank=_union_group->translateRank(&_distant_group,i);

        if (_intersectsBoundingBox(rank))
          {
            _distant_proc_ids.push_back(rank);
          }
      }
#ifdef USE_DIRECTED_BB
#else
    delete [] minmax;
#endif
  }



  /*!
   * Intersect local bounding box with a given distant bounding box on "irank"
   */
  bool ElementLocator::_intersectsBoundingBox(int irank)
  {
#ifdef USE_DIRECTED_BB
    INTERP_KERNEL::DirectedBoundingBox local_dbb, distant_dbb;
    local_dbb.setData( _domain_bounding_boxes+_union_group->myRank()*local_dbb.dataSize( _local_cell_mesh_space_dim ));
    distant_dbb.setData( _domain_bounding_boxes+irank*distant_dbb.dataSize( _local_cell_mesh_space_dim ));
    return !local_dbb.isDisjointWith( distant_dbb );
#else
    double*  local_bb = _domain_bounding_boxes+_union_group->myRank()*2*_local_cell_mesh_space_dim;
    double*  distant_bb =  _domain_bounding_boxes+irank*2*_local_cell_mesh_space_dim;

    const double eps = 1e-12;
    for (int idim=0; idim < _local_cell_mesh_space_dim; idim++)
      {
        bool intersects = (distant_bb[idim*2]<local_bb[idim*2+1]+eps)
          && (local_bb[idim*2]<distant_bb[idim*2+1]+eps);
        if (!intersects) return false; 
      }
    return true;
#endif
  } 


  /*!
   *  Exchange mesh data
   */
  void ElementLocator::_exchangeMesh( MEDCouplingPointSet* local_mesh,
                                      MEDCouplingPointSet*& distant_mesh,
                                      int iproc_distant,
                                      const DataArrayInt* distant_ids_send,
                                      int*& distant_ids_recv)
  {
    CommInterface comm_interface=_union_group->getCommInterface();
  
    // First stage : exchanging sizes
    // ------------------------------
    vector<double> tinyInfoLocalD,tinyInfoDistantD(1);//not used for the moment
    vector<int> tinyInfoLocal,tinyInfoDistant;
    vector<string> tinyInfoLocalS;
    //Getting tiny info of local mesh to allow the distant proc to initialize and allocate
    //the transmitted mesh.
    local_mesh->getTinySerializationInformation(tinyInfoLocalD,tinyInfoLocal,tinyInfoLocalS);
    tinyInfoLocal.push_back(distant_ids_send->getNumberOfTuples());
    tinyInfoDistant.resize(tinyInfoLocal.size());
    std::fill(tinyInfoDistant.begin(),tinyInfoDistant.end(),0);
    MPIProcessorGroup* group=static_cast<MPIProcessorGroup*> (_union_group);
    const MPI_Comm* comm=group->getComm();
    MPI_Status status; 
    
    // iproc_distant is the number of proc in distant group
    // it must be converted to union numbering before communication
    int iprocdistant_in_union = group->translateRank(&_distant_group,
                                                     iproc_distant);
    
    comm_interface.sendRecv(&tinyInfoLocal[0], tinyInfoLocal.size(), MPI_INT, iprocdistant_in_union, 1112,
                            &tinyInfoDistant[0], tinyInfoDistant.size(), MPI_INT,iprocdistant_in_union,1112,
                            *comm, &status);
    DataArrayInt *v1Local=0;
    DataArrayDouble *v2Local=0;
    DataArrayInt *v1Distant=DataArrayInt::New();
    DataArrayDouble *v2Distant=DataArrayDouble::New();
    //serialization of local mesh to send data to distant proc.
    local_mesh->serialize(v1Local,v2Local);
    //Building the right instance of copy of distant mesh.
    MEDCouplingPointSet *distant_mesh_tmp=MEDCouplingPointSet::BuildInstanceFromMeshType((MEDCouplingMeshType)tinyInfoDistant[0]);
    std::vector<std::string> unusedTinyDistantSts;
    distant_mesh_tmp->resizeForUnserialization(tinyInfoDistant,v1Distant,v2Distant,unusedTinyDistantSts);
    int nbLocalElems=0;
    int nbDistElem=0;
    int *ptLocal=0;
    int *ptDist=0;
    if(v1Local)
      {
        nbLocalElems=v1Local->getNbOfElems();
        ptLocal=v1Local->getPointer();
      }
    if(v1Distant)
      {
        nbDistElem=v1Distant->getNbOfElems();
        ptDist=v1Distant->getPointer();
      }
    comm_interface.sendRecv(ptLocal, nbLocalElems, MPI_INT,
                            iprocdistant_in_union, 1111,
                            ptDist, nbDistElem, MPI_INT,
                            iprocdistant_in_union,1111,
                            *comm, &status);
    nbLocalElems=0;
    double *ptLocal2=0;
    double *ptDist2=0;
    if(v2Local)
      {
        nbLocalElems=v2Local->getNbOfElems();
        ptLocal2=v2Local->getPointer();
      }
    nbDistElem=0;
    if(v2Distant)
      {
        nbDistElem=v2Distant->getNbOfElems();
        ptDist2=v2Distant->getPointer();
      }
    comm_interface.sendRecv(ptLocal2, nbLocalElems, MPI_DOUBLE,
                            iprocdistant_in_union, 1112,
                            ptDist2, nbDistElem, MPI_DOUBLE,
                            iprocdistant_in_union, 1112, 
                            *comm, &status);
    //
    distant_mesh=distant_mesh_tmp;
    //finish unserialization
    distant_mesh->unserialization(tinyInfoDistantD,tinyInfoDistant,v1Distant,v2Distant,unusedTinyDistantSts);
    //
    distant_ids_recv=new int[tinyInfoDistant.back()];
    comm_interface.sendRecv(const_cast<void *>(reinterpret_cast<const void *>(distant_ids_send->getConstPointer())),tinyInfoLocal.back(), MPI_INT,
                            iprocdistant_in_union, 1113,
                            distant_ids_recv,tinyInfoDistant.back(), MPI_INT,
                            iprocdistant_in_union,1113,
                            *comm, &status);
    if(v1Local)
      v1Local->decrRef();
    if(v2Local)
      v2Local->decrRef();
    if(v1Distant)
      v1Distant->decrRef();
    if(v2Distant)
      v2Distant->decrRef();
  }
  
  /*!
   * connected with ElementLocator::sendPolicyToWorkingSideL
   */
  void ElementLocator::recvPolicyFromLazySideW(std::vector<int>& policy)
  {
    policy.resize(_distant_proc_ids.size());
    int procId=0;
    CommInterface comm;
    MPI_Status status;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        int toRecv;
        comm.recv((void *)&toRecv,1,MPI_INT,*iter,1120,*_comm,&status);
        policy[procId]=toRecv;
      }
  }

  /*!
   * connected with ElementLocator::recvFromWorkingSideL
   */
  void ElementLocator::sendSumToLazySideW(const std::vector< std::vector<int> >& distantLocEltIds, const std::vector< std::vector<double> >& partialSumRelToDistantIds)
  {
    int procId=0;
    CommInterface comm;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        const vector<int>& eltIds=distantLocEltIds[procId];
        const vector<double>& valued=partialSumRelToDistantIds[procId];
        int lgth=eltIds.size();
        comm.send(&lgth,1,MPI_INT,*iter,1114,*_comm);
        comm.send(const_cast<void *>(reinterpret_cast<const void *>(&eltIds[0])),lgth,MPI_INT,*iter,1115,*_comm);
        comm.send(const_cast<void *>(reinterpret_cast<const void *>(&valued[0])),lgth,MPI_DOUBLE,*iter,1116,*_comm);
      }
  }

  /*!
   * connected with ElementLocator::sendToWorkingSideL
   */
  void ElementLocator::recvSumFromLazySideW(std::vector< std::vector<double> >& globalSumRelToDistantIds)
  {
    int procId=0;
    CommInterface comm;
    MPI_Status status;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        std::vector<double>& vec=globalSumRelToDistantIds[procId];
        comm.recv(&vec[0],vec.size(),MPI_DOUBLE,*iter,1117,*_comm,&status);
      }
  }

  /*!
   * connected with ElementLocator::recvLocalIdsFromWorkingSideL
   */
  void ElementLocator::sendLocalIdsToLazyProcsW(const std::vector< std::vector<int> >& distantLocEltIds)
  {
    int procId=0;
    CommInterface comm;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        const vector<int>& eltIds=distantLocEltIds[procId];
        int lgth=eltIds.size();
        comm.send(&lgth,1,MPI_INT,*iter,1121,*_comm);
        comm.send(const_cast<void *>(reinterpret_cast<const void *>(&eltIds[0])),lgth,MPI_INT,*iter,1122,*_comm);
      }
  }

  /*!
   * connected with ElementLocator::sendGlobalIdsToWorkingSideL
   */
  void ElementLocator::recvGlobalIdsFromLazyProcsW(const std::vector< std::vector<int> >& distantLocEltIds, std::vector< std::vector<int> >& globalIds)
  {
    int procId=0;
    CommInterface comm;
    MPI_Status status;
    globalIds.resize(_distant_proc_ids.size());
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        const std::vector<int>& vec=distantLocEltIds[procId];
        std::vector<int>& global=globalIds[procId];
        global.resize(vec.size());
        comm.recv(&global[0],vec.size(),MPI_INT,*iter,1123,*_comm,&status);
      }
  }
  
  /*!
   * connected with ElementLocator::sendCandidatesGlobalIdsToWorkingSideL
   */
  void ElementLocator::recvCandidatesGlobalIdsFromLazyProcsW(std::vector< std::vector<int> >& globalIds)
  {
    int procId=0;
    CommInterface comm;
    MPI_Status status;
    globalIds.resize(_distant_proc_ids.size());
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        std::vector<int>& global=globalIds[procId];
        int lgth;
        comm.recv(&lgth,1,MPI_INT,*iter,1132,*_comm,&status);
        global.resize(lgth);
        comm.recv(&global[0],lgth,MPI_INT,*iter,1133,*_comm,&status);
      }
  }
  
  /*!
   * connected with ElementLocator::recvSumFromWorkingSideL
   */
  void ElementLocator::sendPartialSumToLazyProcsW(const std::vector<int>& distantGlobIds, const std::vector<double>& sum)
  {
    int procId=0;
    CommInterface comm;
    int lgth=distantGlobIds.size();
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        comm.send(&lgth,1,MPI_INT,*iter,1124,*_comm);
        comm.send(const_cast<void *>(reinterpret_cast<const void *>(&distantGlobIds[0])),lgth,MPI_INT,*iter,1125,*_comm);
        comm.send(const_cast<void *>(reinterpret_cast<const void *>(&sum[0])),lgth,MPI_DOUBLE,*iter,1126,*_comm);
      }
  }

  /*!
   * connected with ElementLocator::recvCandidatesForAddElementsL
   */
  void ElementLocator::sendCandidatesForAddElementsW(const std::vector<int>& distantGlobIds)
  {
    int procId=0;
    CommInterface comm;
    int lgth=distantGlobIds.size();
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        comm.send(const_cast<void *>(reinterpret_cast<const void *>(&lgth)),1,MPI_INT,*iter,1128,*_comm);
        comm.send(const_cast<void *>(reinterpret_cast<const void *>(&distantGlobIds[0])),lgth,MPI_INT,*iter,1129,*_comm);
      }
  }
  
  /*!
   * connected with ElementLocator::sendAddElementsToWorkingSideL
   */
  void ElementLocator::recvAddElementsFromLazyProcsW(std::vector<std::vector<int> >& elementsToAdd)
  {
    int procId=0;
    CommInterface comm;
    MPI_Status status;
    int lgth=_distant_proc_ids.size();
    elementsToAdd.resize(lgth);
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        int locLgth;
        std::vector<int>& eltToFeed=elementsToAdd[procId];
        comm.recv(&locLgth,1,MPI_INT,*iter,1130,*_comm,&status);
        eltToFeed.resize(locLgth);
        comm.recv(&eltToFeed[0],locLgth,MPI_INT,*iter,1131,*_comm,&status);
      }
  }

  /*!
   * connected with ElementLocator::recvPolicyFromLazySideW
   */
  int ElementLocator::sendPolicyToWorkingSideL()
  {
    CommInterface comm;
    int toSend;
    DataArrayInt *isCumulative=_local_para_field.returnCumulativeGlobalNumbering();
    if(isCumulative)
      {
        toSend=CUMULATIVE_POLICY;
        isCumulative->decrRef();
      }
    else
      toSend=NO_POST_TREATMENT_POLICY;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++)
      comm.send(&toSend,1,MPI_INT,*iter,1120,*_comm);
    return toSend;
  }

  /*!
   * connected with ElementLocator::sendSumToLazySideW
   */
  void ElementLocator::recvFromWorkingSideL()
  {
    _values_added.resize(_local_para_field.getField()->getNumberOfTuples());
    int procId=0;
    CommInterface comm;
    _ids_per_working_proc.resize(_distant_proc_ids.size());
    MPI_Status status;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
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
   * connected with ElementLocator::recvSumFromLazySideW
   */
  void ElementLocator::sendToWorkingSideL()
  {
    int procId=0;
    CommInterface comm;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        vector<int>& ids=_ids_per_working_proc[procId];
        vector<double> valsToSend(ids.size());
        vector<double>::iterator iter3=valsToSend.begin();
        for(vector<int>::const_iterator iter2=ids.begin();iter2!=ids.end();iter2++,iter3++)
          *iter3=_values_added[*iter2];
        comm.send(&valsToSend[0],ids.size(),MPI_DOUBLE,*iter,1117,*_comm);
        //ids.clear();
      }
    //_ids_per_working_proc.clear();
  }

  /*!
   * connected with ElementLocator::sendLocalIdsToLazyProcsW
   */
  void ElementLocator::recvLocalIdsFromWorkingSideL()
  {
    int procId=0;
    CommInterface comm;
    _ids_per_working_proc.resize(_distant_proc_ids.size());
    MPI_Status status;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        int lgth;
        vector<int>& ids=_ids_per_working_proc[procId];
        comm.recv(&lgth,1,MPI_INT,*iter,1121,*_comm,&status);
        ids.resize(lgth);
        comm.recv(&ids[0],lgth,MPI_INT,*iter,1122,*_comm,&status);
      }
  }

  /*!
   * connected with ElementLocator::recvGlobalIdsFromLazyProcsW
   */
  void ElementLocator::sendGlobalIdsToWorkingSideL()
  {
    int procId=0;
    CommInterface comm;
    DataArrayInt *globalIds=_local_para_field.returnGlobalNumbering();
    const int *globalIdsC=globalIds->getConstPointer();
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        const vector<int>& ids=_ids_per_working_proc[procId];
        vector<int> valsToSend(ids.size());
        vector<int>::iterator iter1=valsToSend.begin();
        for(vector<int>::const_iterator iter2=ids.begin();iter2!=ids.end();iter2++,iter1++)
          *iter1=globalIdsC[*iter2];
        comm.send(&valsToSend[0],ids.size(),MPI_INT,*iter,1123,*_comm);
      }
    if(globalIds)
      globalIds->decrRef();
  }

  /*!
   * connected with ElementLocator::sendPartialSumToLazyProcsW
   */
  void ElementLocator::recvSumFromWorkingSideL()
  {
    int procId=0;
    int wProcSize=_distant_proc_ids.size();
    CommInterface comm;
    _ids_per_working_proc.resize(wProcSize);
    _values_per_working_proc.resize(wProcSize);
    MPI_Status status;
    std::map<int,double> sums;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        int lgth;
        comm.recv(&lgth,1,MPI_INT,*iter,1124,*_comm,&status);
        vector<int>& ids=_ids_per_working_proc[procId];
        vector<double>& vals=_values_per_working_proc[procId];
        ids.resize(lgth);
        vals.resize(lgth);
        comm.recv(&ids[0],lgth,MPI_INT,*iter,1125,*_comm,&status);
        comm.recv(&vals[0],lgth,MPI_DOUBLE,*iter,1126,*_comm,&status);
        vector<int>::const_iterator iter1=ids.begin();
        vector<double>::const_iterator iter2=vals.begin();
        for(;iter1!=ids.end();iter1++,iter2++)
          sums[*iter1]+=*iter2;
      }
    //assign sum to prepare sending to working side
    for(procId=0;procId<wProcSize;procId++)
      {
        vector<int>& ids=_ids_per_working_proc[procId];
        vector<double>& vals=_values_per_working_proc[procId];
        vector<int>::const_iterator iter1=ids.begin();
        vector<double>::iterator iter2=vals.begin();
        for(;iter1!=ids.end();iter1++,iter2++)
          *iter2=sums[*iter1];
        ids.clear();
      }
  }

  /*!
   * Foreach working procs Wi compute and push it in _ids_per_working_proc3,
   * if it exist, local id of nodes that are in interaction with an another lazy proc than this
   * and that exists in this \b but with no interaction with this.
   * The computation is performed here. sendAddElementsToWorkingSideL is only in charge to send
   * precomputed _ids_per_working_proc3 attribute.
   * connected with ElementLocator::sendCandidatesForAddElementsW
   */
  void ElementLocator::recvCandidatesForAddElementsL()
  {
    int procId=0;
    int wProcSize=_distant_proc_ids.size();
    CommInterface comm;
    _ids_per_working_proc3.resize(wProcSize);
    MPI_Status status;
    std::map<int,double> sums;
    DataArrayInt *globalIds=_local_para_field.returnGlobalNumbering();
    const int *globalIdsC=globalIds->getConstPointer();
    int nbElts=globalIds->getNumberOfTuples();
    std::set<int> globalIdsS(globalIdsC,globalIdsC+nbElts);
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        const std::vector<int>& ids0=_ids_per_working_proc[procId];
        int lgth0=ids0.size();
        std::set<int> elts0;
        for(int i=0;i<lgth0;i++)
          elts0.insert(globalIdsC[ids0[i]]);
        int lgth;
        comm.recv(&lgth,1,MPI_INT,*iter,1128,*_comm,&status);
        vector<int> ids(lgth);
        comm.recv(&ids[0],lgth,MPI_INT,*iter,1129,*_comm,&status);
        set<int> ids1(ids.begin(),ids.end());
        ids.clear();
        set<int> tmp5,tmp6;
        set_intersection(globalIdsS.begin(),globalIdsS.end(),ids1.begin(),ids1.end(),inserter(tmp5,tmp5.begin()));
        set_difference(tmp5.begin(),tmp5.end(),elts0.begin(),elts0.end(),inserter(tmp6,tmp6.begin()));
        std::vector<int>& ids2=_ids_per_working_proc3[procId];
        ids2.resize(tmp6.size());
        std::copy(tmp6.begin(),tmp6.end(),ids2.begin());
        //global->local
        for(std::vector<int>::iterator iter2=ids2.begin();iter2!=ids2.end();iter2++)
          *iter2=std::find(globalIdsC,globalIdsC+nbElts,*iter2)-globalIdsC;
      }
    if(globalIds)
      globalIds->decrRef();
  }

  /*!
   * connected with ElementLocator::recvAddElementsFromLazyProcsW
   */
  void ElementLocator::sendAddElementsToWorkingSideL()
  {
    int procId=0;
    CommInterface comm;
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        const std::vector<int>& vals=_ids_per_working_proc3[procId];
        int size=vals.size();
        comm.send(const_cast<void *>(reinterpret_cast<const void *>(&size)),1,MPI_INT,*iter,1130,*_comm);
        comm.send(const_cast<void *>(reinterpret_cast<const void *>(&vals[0])),size,MPI_INT,*iter,1131,*_comm);
      }
  }

  /*!
   * This method sends to working side Wi only nodes in interaction with Wi \b and located on boundary, to reduce number.
   * connected with ElementLocator::recvCandidatesGlobalIdsFromLazyProcsW
   */
  void ElementLocator::sendCandidatesGlobalIdsToWorkingSideL()
  { 
    int procId=0;
    CommInterface comm;
    DataArrayInt *globalIds=_local_para_field.returnGlobalNumbering();
    const int *globalIdsC=globalIds->getConstPointer();
    MCAuto<DataArrayInt> candidates=_local_para_field.getSupport()->getCellMesh()->findBoundaryNodes();
    for(int *iter1=candidates->getPointer();iter1!=candidates->getPointer()+candidates->getNumberOfTuples();iter1++)
      (*iter1)=globalIdsC[*iter1];
    std::set<int> candidatesS(candidates->begin(),candidates->end());
    for(vector<int>::const_iterator iter=_distant_proc_ids.begin();iter!=_distant_proc_ids.end();iter++,procId++)
      {
        const vector<int>& ids=_ids_per_working_proc[procId];
        vector<int> valsToSend(ids.size());
        vector<int>::iterator iter1=valsToSend.begin();
        for(vector<int>::const_iterator iter2=ids.begin();iter2!=ids.end();iter2++,iter1++)
          *iter1=globalIdsC[*iter2];
        std::set<int> tmp2(valsToSend.begin(),valsToSend.end());
        std::vector<int> tmp3;
        set_intersection(candidatesS.begin(),candidatesS.end(),tmp2.begin(),tmp2.end(),std::back_insert_iterator< std::vector<int> >(tmp3));
        int lgth=tmp3.size();
        comm.send(&lgth,1,MPI_INT,*iter,1132,*_comm);
        comm.send(&tmp3[0],lgth,MPI_INT,*iter,1133,*_comm);
      }
    if(globalIds)
      globalIds->decrRef();
  }
}
