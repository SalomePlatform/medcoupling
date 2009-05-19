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

#include <map>
#include <set>
#include <limits>

using namespace std;

namespace ParaMEDMEM 
{ 
  ElementLocator::ElementLocator(const ParaFIELD& sourceField,
                                 const ProcessorGroup& distant_group)
    : _local_para_field(sourceField),
      _local_cell_mesh(sourceField.getSupport()->getCellMesh()),
      _local_face_mesh(sourceField.getSupport()->getFaceMesh()),
      _local_group(*sourceField.getSupport()->getBlockTopology()->getProcGroup()),
      _distant_group(distant_group)
  { 
    _union_group = _local_group.fuse(distant_group);
    _computeBoundingBoxes();
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

  // ==========================================================================
  // Procedure for exchanging mesh between a distant proc and a local processor
  // param idistantrank  proc id on distant group
  // param distant_mesh on return , points to a local reconstruction of
  //  the distant mesh
  // param distant_ids on return, contains a vector defining a correspondence
  // between the distant ids and the ids of the local reconstruction 
  // ==========================================================================
  void ElementLocator::exchangeMesh(int idistantrank,
                                    MEDCouplingPointSet*& distant_mesh,
                                    int*& distant_ids)
  {
    int dim  = _local_cell_mesh->getSpaceDimension();
    int rank = _union_group->translateRank(&_distant_group,idistantrank);

    if (find(_distant_proc_ids.begin(), _distant_proc_ids.end(),rank)==_distant_proc_ids.end())
      {
        return;
      }
   
    vector<int> elems;
    double* distant_bb =  _domain_bounding_boxes+rank*2*dim;
    _local_cell_mesh->giveElemsInBoundingBox(distant_bb,getBoundingBoxAdjustment(),elems);
    
    DataArrayInt *distant_ids_send;
    MEDCouplingPointSet *send_mesh = (MEDCouplingPointSet *)_local_para_field.getField()->buildSubMeshData(&elems[0],&elems[elems.size()],distant_ids_send);
    _exchangeMesh(send_mesh, distant_mesh, idistantrank, distant_ids_send, distant_ids);
    distant_ids_send->decrRef();
    
    if(send_mesh)
      send_mesh->decrRef();
#if 0
    int* distant_ids_send=0;
    //send_mesh contains null pointer if elems is empty
    MEDCouplingPointSet* send_mesh = _local_cell_mesh->buildPartOfMySelf(&elems[0],&elems[elems.size()],false);
    // Constituting an array containing the ids of the elements that are 
    // going to be sent to the distant subdomain.
    // This array  enables the correct redistribution of the data when the
    // interpolated field is transmitted to the target array

    if (elems.size()>0)
      {
        distant_ids_send = new int[elems.size()];
        std::copy(elems.begin(),elems.end(),distant_ids_send);
      }
    _exchangeMesh(send_mesh, distant_mesh, idistantrank, distant_ids_send, distant_ids);
    delete[] distant_ids_send;
    if(send_mesh)
      send_mesh->decrRef();
#endif
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


  // ======================
  // Compute bounding boxes
  // ======================

  void ElementLocator::_computeBoundingBoxes()
  {
    CommInterface comm_interface =_union_group->getCommInterface();
    int dim = _local_cell_mesh->getSpaceDimension();
    _domain_bounding_boxes = new double[2*dim*_union_group->size()];
    double * minmax=new double [2*dim];
    _local_cell_mesh->getBoundingBox(minmax);

    MPIProcessorGroup* group=static_cast<MPIProcessorGroup*> (_union_group);
    const MPI_Comm* comm = group->getComm();
    comm_interface.allGather(minmax, 2*dim, MPI_DOUBLE,
                             _domain_bounding_boxes,2*dim, MPI_DOUBLE, 
                             *comm);
  
    for (int i=0; i< _distant_group.size(); i++)
      {
        int rank= _union_group->translateRank(&_distant_group,i);

        if (_intersectsBoundingBox(rank))
          {
            _distant_proc_ids.push_back(rank);
          }
      }
    delete[] minmax;
  }


  // =============================================
  // Intersect Bounding Box (with a given "irank")
  // =============================================
  bool ElementLocator::_intersectsBoundingBox(int irank)
  {
    int dim=_local_cell_mesh->getSpaceDimension();
    double*  local_bb = _domain_bounding_boxes+_union_group->myRank()*2*dim;
    double*  distant_bb =  _domain_bounding_boxes+irank*2*dim;

    for (int idim=0; idim < _local_cell_mesh->getSpaceDimension(); idim++)
      {
        const double eps =  1e-12;
        bool intersects = (distant_bb[idim*2]<local_bb[idim*2+1]+eps)
          && (local_bb[idim*2]<distant_bb[idim*2+1]+eps);
        if (!intersects) return false; 
      }
    return true;
  } 

  // ======================
  // Exchanging meshes data
  // ======================
  void ElementLocator::_exchangeMesh( MEDCouplingPointSet* local_mesh,
                                      MEDCouplingPointSet*& distant_mesh,
                                      int iproc_distant,
                                      const DataArrayInt* distant_ids_send,
                                      int*& distant_ids_recv)
  {
    CommInterface comm_interface=_union_group->getCommInterface();
  
    // First stage : exchanging sizes
    // ------------------------------
    vector<int> tinyInfoLocal,tinyInfoDistant;
    vector<string> tinyInfoLocalS;
    //Getting tiny info of local mesh to allow the distant proc to initialize and allocate
    //the transmitted mesh.
    local_mesh->getTinySerializationInformation(tinyInfoLocal,tinyInfoLocalS);
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
    MEDCouplingPointSet *distant_mesh_tmp=MEDCouplingPointSet::buildInstanceFromMeshType((MEDCouplingMeshType)tinyInfoDistant[0]);
    std::vector<std::string> unusedTinyDistantSts;
    distant_mesh_tmp->resizeForUnserialization(tinyInfoDistant,v1Distant,v2Distant,unusedTinyDistantSts);
    comm_interface.sendRecv(v1Local->getPointer(), v1Local->getNbOfElems(), MPI_INT,
                            iprocdistant_in_union, 1111,
                            v1Distant->getPointer(), v1Distant->getNbOfElems(), MPI_INT,
                            iprocdistant_in_union,1111,
                            *comm, &status);
    comm_interface.sendRecv(v2Local->getPointer(), v2Local->getNbOfElems(), MPI_DOUBLE,
                            iprocdistant_in_union, 1112,
                            v2Distant->getPointer(), v2Distant->getNbOfElems(), MPI_DOUBLE,
                            iprocdistant_in_union, 1112, 
                            *comm, &status);
    if(v1Distant->getNbOfElems()>0)
      {
        distant_mesh=distant_mesh_tmp;
        //finish unserialization
        distant_mesh->unserialization(tinyInfoDistant,v1Distant,v2Distant,unusedTinyDistantSts);
      }
    else
      distant_mesh_tmp->decrRef();
    distant_ids_recv=new int[tinyInfoDistant.back()];
    comm_interface.sendRecv((void *)distant_ids_send->getConstPointer(),tinyInfoLocal.back(), MPI_INT,
                            iprocdistant_in_union, 1113,
                            distant_ids_recv,tinyInfoDistant.back(), MPI_INT,
                            iprocdistant_in_union,1113,
                            *comm, &status);
    v1Local->decrRef();
    v2Local->decrRef();
    v1Distant->decrRef();
    v2Distant->decrRef();
  }
}
