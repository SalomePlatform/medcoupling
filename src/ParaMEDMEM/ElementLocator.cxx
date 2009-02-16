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
#include "ParaMESH.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"

#include <map>
#include <set>
#include <limits>

using namespace std;

namespace ParaMEDMEM 
{ 
  ElementLocator::ElementLocator(const ParaMESH& sourceMesh,
                                 const ProcessorGroup& distant_group)
    : _local_para_mesh(sourceMesh),
      _local_cell_mesh(sourceMesh.getCellMesh()),
      _local_face_mesh(sourceMesh.getFaceMesh()),
      _local_group(*sourceMesh.getBlockTopology()->getProcGroup()),
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

  // ==========================================================================
  // Procedure for exchanging mesh between a distant proc and a local processor
  // param idistantrank  proc id on distant group
  // param distant_mesh on return , points to a local reconstruction of
  //  the distant mesh
  // param distant_ids on return, contains a vector defining a correspondence
  // between the distant ids and the ids of the local reconstruction 
  // ==========================================================================
  void ElementLocator::exchangeMesh(int idistantrank,
                                    MEDCouplingUMesh*& distant_mesh,
                                    int*& distant_ids)
  {
    int dim  = _local_cell_mesh->getSpaceDimension();
    int rank = _union_group->translateRank(&_distant_group,idistantrank);

    if (find(_distant_proc_ids.begin(), _distant_proc_ids.end(),rank)==_distant_proc_ids.end())
      {
        return;
      }
   
    set <int> elems;
    double* distant_bb =  _domain_bounding_boxes+rank*2*dim;
    double* elem_bb=new double[2*dim];

    //defining pointers to med
    const int* conn      = _local_cell_mesh->getNodalConnectivity()->getPointer() ;
    const int* conn_index= _local_cell_mesh->getNodalConnectivityIndex()->getPointer();
    const double* coords = _local_cell_mesh->getCoords()->getPointer() ;
   
    for ( int ielem=0; ielem<_local_cell_mesh->getNumberOfCells() ; ielem++)
      {
        for (int i=0; i<dim; i++)
          {
            elem_bb[i*2]=std::numeric_limits<double>::max();
            elem_bb[i*2+1]=-std::numeric_limits<double>::max();
          }

        for (int inode=conn_index[ielem]+1; inode<conn_index[ielem+1]; inode++)//+1 due to offset of cell type.
          {
            int node= conn[inode];
     
            for (int idim=0; idim<dim; idim++)
              {
                if ( coords[node*dim+idim] < elem_bb[idim*2] )
                  {
                    elem_bb[idim*2] = coords[node*dim+idim] ;
                  }
                if ( coords[node*dim+idim] > elem_bb[idim*2+1] )
                  {
                    elem_bb[idim*2+1] = coords[node*dim+idim] ;
                  }
              }
          }
        if (_intersectsBoundingBox(elem_bb, distant_bb, dim))
          {
            elems.insert(ielem);
          }
      }
    //send_mesh contains null pointer if elems is empty
    MEDCouplingUMesh* send_mesh = _meshFromElems(elems);
    
    // Constituting an array containing the ids of the elements that are 
    // going to be sent to the distant subdomain.
    // This array  enables the correct redistribution of the data when the
    // interpolated field is transmitted to the target array

    int* distant_ids_send=0;
    if (elems.size()>0)
      {
        distant_ids_send = new int[elems.size()];
        int index=0;
        for (std::set<int>::const_iterator iter = elems.begin(); iter!= elems.end(); iter++)
          {
            distant_ids_send[index]=*iter;
            index++;
          }
      }
    _exchangeMesh(send_mesh, distant_mesh, idistantrank, distant_ids_send, distant_ids);
    delete[] distant_ids_send;
    delete[] elem_bb;
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


  // ======================
  // Compute bounding boxes
  // ======================

  void ElementLocator::_computeBoundingBoxes()
  {
    CommInterface comm_interface =_union_group->getCommInterface();
    int dim = _local_cell_mesh->getSpaceDimension();
    _domain_bounding_boxes = new double[2*dim*_union_group->size()];
    const double* coords = _local_cell_mesh->getCoords()->getPointer() ;
 
    int nbnodes =  _local_cell_mesh->getNumberOfNodes();
    double * minmax=new double [2*dim];
    for (int idim=0; idim<dim; idim++)
      {
        minmax[idim*2]=std::numeric_limits<double>::max();
        minmax[idim*2+1]=-std::numeric_limits<double>::max();
      } 

    for (int i=0; i<nbnodes; i++)
      {
        for (int idim=0; idim<dim;idim++)
          {
            if ( minmax[idim*2] > coords[i*dim+idim] )
              {
                minmax[idim*2] = coords[i*dim+idim] ;
              }
            if ( minmax[idim*2+1] < coords[i*dim+idim] )
              {
                minmax[idim*2+1] = coords[i*dim+idim] ;
              }
          }
      }

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

  // =============================================
  // Intersect Bounding Box given 2 Bounding Boxes
  // =============================================
  bool ElementLocator::_intersectsBoundingBox(double* bb1, double* bb2, int dim)
  {
    double bbtemp[2*dim];
    double deltamax=0.0;
    double adjustment_eps=getBoundingBoxAdjustment();

    for (int i=0; i< dim; i++)
      {
        double delta = bb1[2*i+1]-bb1[2*i];
        if ( delta > deltamax )
          {
            deltamax = delta ;
          }
        //    deltamax = (delta>deltamax)?delta:deltamax;
      }
    for (int i=0; i<dim; i++)
      {
        bbtemp[i*2]=bb1[i*2]-deltamax*adjustment_eps;
        bbtemp[i*2+1]=bb1[i*2+1]+deltamax*adjustment_eps;
      }
  
    for (int idim=0; idim < dim; idim++)
      {
        bool intersects = (bbtemp[idim*2]<bb2[idim*2+1])
          && (bb2[idim*2]<bbtemp[idim*2+1]) ;
        if (!intersects) return false; 
      }
    return true;
  }


  // ======================
  // Exchanging meshes data
  // ======================
  void ElementLocator::_exchangeMesh( MEDCouplingUMesh* local_mesh,
                                      MEDCouplingUMesh*& distant_mesh,
                                      int iproc_distant,
                                      const int* distant_ids_send,
                                      int*& distant_ids_recv)
  {
    CommInterface comm_interface=_union_group->getCommInterface();
  
    // First stage : exchanging sizes
    // ------------------------------

    int* send_buffer = new int[5];
    int* recv_buffer = new int[5];
 
    //treatment for non-empty mesh
    int nbconn=0;
    int nbelems=0;

    if (local_mesh !=0)
      {
        nbelems = local_mesh->getNumberOfCells();
        nbconn = local_mesh->getMeshLength();
        send_buffer[0] = local_mesh->getSpaceDimension();
        send_buffer[1] = local_mesh->getMeshDimension();
        send_buffer[2] = local_mesh->getNumberOfNodes();
        send_buffer[3] = nbelems;
        send_buffer[4] = nbconn;
      }
    else
      {
        for (int i=0; i<5; i++)
          {
            send_buffer[i]=0;
          }
      }

    MPIProcessorGroup* group=static_cast<MPIProcessorGroup*> (_union_group);
    const MPI_Comm* comm=(group->getComm());
    MPI_Status status; 

    // iproc_distant is the number of proc in distant group
    // it must be converted to union numbering before communication
    int iprocdistant_in_union = group->translateRank(&_distant_group,
                                                     iproc_distant);

    comm_interface.sendRecv(send_buffer, 5, MPI_INT, iprocdistant_in_union, 1112,
                            recv_buffer, 5, MPI_INT,iprocdistant_in_union,1112,
                            *comm, &status);

    int distant_space_dim = recv_buffer[0];
    int distant_mesh_dim  = recv_buffer[1];
    int distant_nnodes    = recv_buffer[2];
    int distant_nb_elems  = recv_buffer[3];
    int distant_nb_conn   = recv_buffer[4];
  
    delete[] send_buffer;
    delete[] recv_buffer;
  
    // Second stage : exchanging connectivity buffers
    // ----------------------------------------------

    int nb_integers = nbconn + 2*nbelems + 1;
    send_buffer     = new int[nb_integers];
    const int* conn = 0;
    const int* global_numbering=0;
    int* ptr_buffer = send_buffer;   

    if (local_mesh != 0)
      {
        conn = local_mesh->getNodalConnectivity()->getPointer();
      
        global_numbering = local_mesh->getNodalConnectivityIndex()->getPointer() ;
      
        //copying the data in the integer buffer
      
        memcpy(ptr_buffer, global_numbering,  (nbelems+1)*sizeof(int));
        ptr_buffer += nbelems+1;
        memcpy(ptr_buffer,conn, nbconn*sizeof(int));
        ptr_buffer += nbconn;
        memcpy(ptr_buffer, distant_ids_send,  nbelems*sizeof(int));
      }

    // Preparing the recv buffers
    int nb_recv_integers = distant_nb_conn + 2*distant_nb_elems + 1 ;
    recv_buffer=new int[nb_recv_integers];
  
    // Exchanging  integer buffer
    comm_interface.sendRecv(send_buffer, nb_integers, MPI_INT,
                            iprocdistant_in_union, 1111,
                            recv_buffer, nb_recv_integers, MPI_INT,
                            iprocdistant_in_union,1111,
                            *comm, &status);
                 
    if ( nb_integers>0 )
      {
        delete[] send_buffer;
      }

    // Third stage : exchanging coordinates  
    // ------------------------------------

    int nb_recv_floats = distant_space_dim*distant_nnodes;
    int nb_send_floats = 0;
    double* coords=0;
 
    if ( local_mesh!=0 )
      {
        nb_send_floats = local_mesh->getSpaceDimension()
          * local_mesh->getNumberOfNodes();
        coords = local_mesh->getCoords()->getPointer();
      }
  
    DataArrayDouble* myCoords=DataArrayDouble::New();
    myCoords->alloc(distant_nnodes,distant_space_dim);

    comm_interface.sendRecv(coords, nb_send_floats, MPI_DOUBLE,
                            iprocdistant_in_union, 1112,
                            myCoords->getPointer(), nb_recv_floats, MPI_DOUBLE,
                            iprocdistant_in_union, 1112, 
                            *group->getComm(), &status);
  

    // Reconstructing an image of the distant mesh locally
  
    if ( nb_recv_integers>0 && distant_space_dim !=0 ) 
      {
        MEDCouplingUMesh* meshing = MEDCouplingUMesh::New() ;

        // Coordinates
        meshing->setCoords(myCoords) ;
        myCoords->decrRef();
        // Connectivity

        int *work=recv_buffer;
        DataArrayInt* myConnecIndex=DataArrayInt::New();
        myConnecIndex->alloc(distant_nb_elems+1,1);
        memcpy(myConnecIndex->getPointer(), work, (distant_nb_elems+1)*sizeof(int));
        work += distant_nb_elems + 1 ;
    
        DataArrayInt* myConnec=DataArrayInt::New();
        myConnec->alloc(distant_nb_conn,1);
        memcpy(myConnec->getPointer(), work, (distant_nb_conn)*sizeof(int));
        work+=distant_nb_conn;
        meshing->setConnectivity(myConnec, myConnecIndex) ;
        myConnec->decrRef();
        myConnecIndex->decrRef();

        // correspondence between the distant ids and the ids of
        // the local reconstruction

        distant_ids_recv=new int [distant_nb_elems];
        for (int i=0; i<distant_nb_elems; i++)
          {
            distant_ids_recv[i]=*work++;
          }

        // Mesh dimension
        meshing->setMeshDimension(distant_mesh_dim);

        distant_mesh=meshing;  
        delete[] recv_buffer;
      }

  }


  // ==============
  // _meshFromElems
  // ==============

  MEDCouplingUMesh* ElementLocator::_meshFromElems(set<int>& elems)
  {
    //returns null pointer if there are no elems in the mesh
    if ( elems.size()==0 ) return 0;

    // Defining pointers
    const int* conn_mesh =
      const_cast<int*> (_local_cell_mesh->getNodalConnectivity()->getPointer());

    const int* conn_index =
      const_cast<int*> (_local_cell_mesh->getNodalConnectivityIndex()->getPointer());

    const double* coords =
      const_cast<double*> ( _local_cell_mesh->getCoords()->getPointer());

    set<int> nodes;
    int nbconn=0;
    for (set<int>::const_iterator iter=elems.begin(); iter!=elems.end(); iter++)
      {
        // Conn_index : C-like Addresses
        for (int inode=conn_index[*iter]+1; inode<conn_index[*iter+1]; inode++)
          {
            nodes.insert(conn_mesh[inode]);
            nbconn++ ;
          }
      }

    map<int,int> big2small;
    int i=0;
    for (set<int>::const_iterator iter=nodes.begin(); iter!=nodes.end(); iter++)
      {
        big2small[*iter]=i;
        i++;
      }

    // Memory allocate
    DataArrayInt *conn = DataArrayInt::New() ;
    conn->alloc(nbconn+elems.size(),1) ;
    int *connPtr=conn->getPointer();

    DataArrayInt * connIndex = DataArrayInt::New() ;
    connIndex->alloc(elems.size()+1,1) ;
    int* connIndexPtr=connIndex->getPointer();

    DataArrayDouble *new_coords = DataArrayDouble::New() ;
    new_coords->alloc(nodes.size(), _local_cell_mesh->getSpaceDimension()) ;
    double *new_coords_ptr = new_coords->getPointer();

    // New connectivity table
    int index=0;
    int mainIndex=0;
    for (set<int>::const_iterator iter=elems.begin(); iter!=elems.end(); iter++,mainIndex++)
      {
        connIndexPtr[mainIndex]=index;
        connPtr[index++]=conn_mesh[conn_index[*iter]];
        for (int inode = conn_index[*iter]+1; inode < conn_index[*iter+1]; inode++)
          {
            connPtr[index]=big2small[conn_mesh[inode]] ; // C-like number
            index++;
          } 
      }
    connIndexPtr[mainIndex]=index;
    // Coordinates
    index=0;
    for (set<int>::const_iterator iter=nodes.begin(); iter!=nodes.end(); iter++)
      {
        int dim = _local_cell_mesh->getSpaceDimension();
        for (int i=0; i<dim;i++)
          {
            new_coords_ptr[index]=coords[(*iter)*dim+i];
            index++;
          }
      }
  
    // Initialize
    MEDCouplingUMesh* meshing = MEDCouplingUMesh::New() ;
    meshing->setCoords(new_coords) ;
    new_coords->decrRef();
    meshing->setConnectivity(conn, connIndex) ;
    conn->decrRef();
    connIndex->decrRef();
    meshing->setMeshDimension(_local_cell_mesh->getMeshDimension());

    return meshing;
  }
}
