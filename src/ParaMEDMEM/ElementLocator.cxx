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
#include "ParaSUPPORT.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Meshing.hxx"

#include <set>
#include <limits>
using namespace std;


namespace ParaMEDMEM 
{ 
 
ElementLocator::ElementLocator(const ParaMESH& mesh, const ProcessorGroup& distant_group) 
:_local_mesh(mesh.getMesh()),
_local_group(*mesh.getBlockTopology()->getProcGroup()),
 _distant_group(distant_group)
{ 
  _union_group = _local_group.fuse(distant_group);
  _computeBoundingBoxes();
}

ElementLocator::ElementLocator(const ParaSUPPORT& support, const ProcessorGroup& distant_group)
:_local_group(*support.getTopology()->getProcGroup()),
_distant_group(distant_group),
_union_group(_local_group.fuse(distant_group))
{
  throw ("Element Locator SUPPORT constructor not implemented yet");
}

ElementLocator::~ElementLocator()
{
  delete _union_group;
  delete [] _domain_bounding_boxes;
}


/*! Procedure for exchanging mesh between a distant proc and a local processor
\param idistantrank  proc id on distant group
\param distant_mesh on return , points to a local reconstruction of the distant mesh
\param distant_ids on return, contains a vector defining a correspondence between the distant ids and the ids of the local reconstruction 
*/

void ElementLocator::exchangeMesh(int idistantrank, MEDMEM::MESH*& distant_mesh, int*& distant_ids)
{
  int dim=_local_mesh->getSpaceDimension();
   int rank= _union_group->translateRank(&_distant_group,idistantrank);
   if (find(_distant_proc_ids.begin(), _distant_proc_ids.end(),rank)==_distant_proc_ids.end())
     return;
   
   set <int> elems;
   double* distant_bb =  _domain_bounding_boxes+rank*2*dim;
   double* elem_bb=new double[2*dim];

   //defining pointers to med
   const int* conn=_local_mesh->getConnectivity(MED_EN::MED_FULL_INTERLACE,
                                               MED_EN::MED_NODAL,
                                               MED_EN::MED_CELL,
                                               MED_EN::MED_ALL_ELEMENTS);
   const int* conn_index= _local_mesh->getConnectivityIndex(
                                               MED_EN::MED_NODAL,
                                               MED_EN::MED_CELL);
   const double* coords = _local_mesh->getCoordinates(MED_EN::MED_FULL_INTERLACE);
   
   for (int ielem=0; ielem<_local_mesh->getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);ielem++)
     {
       for (int i=0; i<dim; i++)
	 {
	   elem_bb[i*2]=std::numeric_limits<double>::max();
	   elem_bb[i*2+1]=-std::numeric_limits<double>::max();
	 }
       for (int inode=conn_index[ielem]; inode<conn_index[ielem+1]; inode++)
	 {
	   int node= conn [inode-1]-1; // - 1 because of MED numbering starting at 1
	   
	   for (int idim=0; idim<dim; idim++)
	     {
	       elem_bb[idim*2]=coords[node*dim+idim]<elem_bb[idim*2]? coords[node*dim+idim]: elem_bb[idim*2];
	       elem_bb[idim*2+1]=coords[node*dim+idim]>elem_bb[idim*2+1]? coords[node*dim+idim]: elem_bb[idim*2+1];
	     }
	 }
       if (_intersectsBoundingBox(elem_bb, distant_bb, dim))
	 {
	   elems.insert(ielem+1);
	 }
     }
   
	 //send_mesh contains null pointer if elems is empty
   MEDMEM::MESH* send_mesh= _meshFromElems(elems);
   
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
  delete send_mesh;
}

void ElementLocator::_computeBoundingBoxes()
{
  CommInterface comm_interface=_union_group->getCommInterface();
  int dim = _local_mesh->getSpaceDimension();
  _domain_bounding_boxes = new double[2*dim*_union_group->size()];
  const double* coords = _local_mesh->getCoordinates(MED_EN::MED_FULL_INTERLACE);
 
  int nbnodes =  _local_mesh->getNumberOfNodes();
  double * minmax=new double [2*dim];
  for (int idim=0; idim<dim; idim++)
  {
    minmax[idim*2]=std::numeric_limits<double>::max();
    minmax[idim*2+1]=-std::numeric_limits<double>::max();
  } 
  for (int i=0; i<nbnodes; i++)
    for (int idim=0; idim<dim;idim++)
    {
      minmax[idim*2]=(minmax[idim*2]<coords[i*dim+idim]?minmax[idim*2]:coords[i*dim+idim]);
      minmax[idim*2+1]=(minmax[idim*2+1]>coords[i*dim+idim]?minmax[idim*2+1]:coords[i*dim+idim]);
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
      _distant_proc_ids.push_back(rank);
  }
  delete[] minmax;
}

bool ElementLocator::_intersectsBoundingBox(int irank)
{
  int dim=_local_mesh->getSpaceDimension();
  double*  local_bb = _domain_bounding_boxes+_union_group->myRank()*2*dim;
  double*  distant_bb =  _domain_bounding_boxes+irank*2*dim;
  for (int idim=0; idim < _local_mesh->getSpaceDimension(); idim++)
  {
    const double eps =  1e-12;
    bool intersects = distant_bb[idim*2]<local_bb[idim*2+1]+eps && local_bb[idim*2]<distant_bb[idim*2+1]+eps;
    if (!intersects) return false; 
  }
  return true;
} 

bool ElementLocator::_intersectsBoundingBox(double* bb1, double* bb2, int dim)
{
	double bbtemp[2*dim];
	double deltamax=0.0;
	double adjustment_eps=getBoundingBoxAdjustment();
	for (int i=0; i< dim; i++)
		{
			double delta = bb1[2*i+1]-bb1[2*i];
			deltamax = (delta>deltamax)?delta:deltamax;
		}
	for (int i=0; i<dim; i++)
		{
			bbtemp[i*2]=bb1[i*2]-deltamax*adjustment_eps;
			bbtemp[i*2+1]=bb1[i*2+1]+deltamax*adjustment_eps;
		}
	
  for (int idim=0; idim < dim; idim++)
  {
    bool intersects = bbtemp[idim*2]<bb2[idim*2+1] && bb2[idim*2]<bbtemp[idim*2+1];
    if (!intersects) return false; 
  }
  return true;
}

void ElementLocator::_exchangeMesh(MEDMEM::MESH* local_mesh, MEDMEM::MESH*& distant_mesh, int iproc_distant, const int* distant_ids_send, int*& distant_ids_recv)
{
 
  CommInterface comm_interface=_union_group->getCommInterface();
  
  // First stage : exchanging sizes
  int* send_buffer = new int[6];
  int* recv_buffer = new int[6];
 
	//treatment for non-empty mesh
	int nbtypes=0;
	int nbconn=0;
	int nbelems=0;

  if (local_mesh !=0)
  {
		nbtypes= local_mesh->getNumberOfTypes(MED_EN::MED_CELL);
		nbconn = local_mesh->getConnectivityLength(MED_EN::MED_FULL_INTERLACE, MED_EN::MED_NODAL,MED_EN::MED_CELL, MED_EN::MED_ALL_ELEMENTS);
		nbelems = local_mesh->getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);
 
    send_buffer[0] = local_mesh->getSpaceDimension();
    send_buffer[1] = local_mesh->getMeshDimension();
    send_buffer[2] = local_mesh->getNumberOfNodes();
    send_buffer[3] = nbelems;
    send_buffer[4] = nbtypes;
    send_buffer[5] = nbconn;
  }
  else
  {
     for (int i=0; i<6; i++)
      send_buffer[i]=0;
  }
 MPIProcessorGroup* group=static_cast<MPIProcessorGroup*> (_union_group);
 const MPI_Comm* comm=(group->getComm());
 MPI_Status status; 
 // iproc_distant is the number of proc in distant group
 // it must be converted to union numbering before communication
 int iprocdistant_in_union = group->translateRank(&_distant_group, iproc_distant);
  comm_interface.sendRecv(send_buffer, 6, MPI_INT, iprocdistant_in_union, 1112, recv_buffer, 6, MPI_INT,iprocdistant_in_union,1112, *comm, &status);
   int distant_space_dim = recv_buffer[0];
  int distant_mesh_dim = recv_buffer[1];
  int distant_nnodes = recv_buffer[2];
  int distant_nb_elems = recv_buffer[3];
  int distant_nb_types = recv_buffer[4];
  int distant_nb_conn = recv_buffer[5];
  
  delete[] send_buffer;
  delete[] recv_buffer;
  
  //Second stage : exchanging connectivity buffers
  int nb_integers = nbtypes*2+nbconn+1+nbelems;
  send_buffer = new int[nb_integers];
	const MED_EN::medGeometryElement* types=0;
  const int* conn = 0;
	const int* global_numbering=0;
	int* ptr_buffer=send_buffer;   
	if (local_mesh!=0)
		{
			conn=local_mesh->getConnectivity(MED_EN::MED_FULL_INTERLACE, MED_EN::MED_NODAL, 
																			 MED_EN::MED_CELL, MED_EN::MED_ALL_ELEMENTS);
			
			global_numbering = local_mesh->getGlobalNumberingIndex(MED_EN::MED_CELL);
			types = local_mesh->getTypes(MED_EN::MED_CELL);
			
			//copying the data in the integer buffer
			
			memcpy(ptr_buffer, types, nbtypes*sizeof(int));
			ptr_buffer+=nbtypes;
			memcpy(ptr_buffer, global_numbering,  (nbtypes+1)*sizeof(int));
			ptr_buffer+=nbtypes+1;
			memcpy(ptr_buffer,conn, nbconn*sizeof(int));
			ptr_buffer+=nbconn;
			memcpy(ptr_buffer, distant_ids_send,  nbelems*sizeof(int));
		}
	//preparing the recv buffers
	int nb_recv_integers = 2*distant_nb_types+1+distant_nb_conn+distant_nb_elems;
	recv_buffer=new int[nb_recv_integers];
  
	//exchanging  integer buffer
  comm_interface.sendRecv(send_buffer, nb_integers, MPI_INT, iprocdistant_in_union, 1111,
                           recv_buffer, nb_recv_integers, MPI_INT, iprocdistant_in_union,1111,
                           *comm, &status);
                 
  if (nb_integers>0) delete[] send_buffer;

  //Third stage : exchanging coordinates  
  int nb_recv_floats = distant_space_dim*distant_nnodes;
	int nb_send_floats=0;
	double* coords=0;
	double* recv_coords=0;
 
	if (local_mesh!=0)
		{
			nb_send_floats =  local_mesh->getSpaceDimension()* local_mesh->getNumberOfNodes();
			coords = const_cast<double*> (local_mesh->getCoordinates(MED_EN::MED_FULL_INTERLACE));
		}
	
	if (nb_recv_floats>0)
		recv_coords = new double[nb_recv_floats];

	comm_interface.sendRecv(coords, nb_send_floats, MPI_DOUBLE, iprocdistant_in_union, 1112,
                           recv_coords, nb_recv_floats, MPI_DOUBLE, iprocdistant_in_union, 1112, 
                           *group->getComm(), &status);
  
  //Reconstructing an image of the distant mesh locally
  
  if (nb_recv_integers>0 && distant_space_dim !=0) 
  {
    MEDMEM::MESHING* meshing = new MEDMEM::MESHING ();
    int* recv_buffer_ptr = recv_buffer;
    meshing->setCoordinates(distant_space_dim, distant_nnodes, recv_coords, "CARTESIAN", MED_EN::MED_FULL_INTERLACE);
		
    meshing->setNumberOfTypes(distant_nb_types,MED_EN::MED_CELL);

    // converting the types from int to medGeometryElement
    MED_EN::medGeometryElement* types_vector = new MED_EN::medGeometryElement[distant_nb_types];
    for (int i=0; i<distant_nb_types; i++)
      types_vector[i]=(MED_EN::medGeometryElement)recv_buffer_ptr[i];
    meshing->setTypes(types_vector, MED_EN::MED_CELL);
    delete[] types_vector;

    recv_buffer_ptr+=distant_nb_types;
    int* nbtypes = new int[distant_nb_types];
    for (int i=0; i<distant_nb_types; i++)
      nbtypes[i]=recv_buffer_ptr[i+1]-recv_buffer_ptr[i];
    recv_buffer_ptr+=distant_nb_types+1;
    meshing->setNumberOfElements( nbtypes, MED_EN::MED_CELL);
                      
    for (int i=0; i<distant_nb_types; i++)
    {
      meshing->setConnectivity(recv_buffer_ptr, MED_EN::MED_CELL,recv_buffer[i]);
      recv_buffer_ptr+=nbtypes[i]*(recv_buffer[i]%100);
    }
    distant_ids_recv=new int [distant_nb_elems];
    for (int i=0; i<distant_nb_elems; i++)
    {
      distant_ids_recv[i]=*recv_buffer_ptr++;
    }
    meshing->setMeshDimension(distant_mesh_dim);

    distant_mesh=meshing;  
    delete[] recv_buffer;
    delete[] nbtypes;
  }
  if (nb_recv_floats >0)
     delete[] recv_coords; // the coordinates are present if the recv_buffer is not empty
  
}

MEDMEM::MESH* ElementLocator::_meshFromElems(set<int>& elems)
{
	//returns null pointer if there are no elems in the mesh
	if (elems.size()==0) return 0;

  //defining pointers to med
   const int* conn_mesh=_local_mesh->getConnectivity(MED_EN::MED_FULL_INTERLACE,
                                               MED_EN::MED_NODAL,
                                               MED_EN::MED_CELL,
                                               MED_EN::MED_ALL_ELEMENTS);
   const int* conn_index= _local_mesh->getConnectivityIndex(
                                               MED_EN::MED_NODAL,
                                               MED_EN::MED_CELL);
   const double* coords = _local_mesh->getCoordinates(MED_EN::MED_FULL_INTERLACE);
  set<int> nodes;
  int nbconn=0;
  map<MED_EN::medGeometryElement,int> nbelems_per_type;
  for (set<int>::const_iterator iter=elems.begin(); iter!=elems.end(); iter++)
  {
    for (int inode = conn_index[*iter-1]-1; inode < conn_index[*iter]-1; inode++)
      nodes.insert(conn_mesh[inode]);
    nbconn+=conn_index[*iter]-conn_index[*iter-1];
    MED_EN::medGeometryElement type = _local_mesh->getElementType(MED_EN::MED_CELL,*iter);
    nbelems_per_type[type]++;
  }
  int* small2big=new int[nodes.size()];
  map<int,int> big2small;
  int i=0;
  for (set<int>::const_iterator iter=nodes.begin(); iter!=nodes.end(); iter++)
  {
    small2big[i]=*iter;
    big2small[*iter]=i;
    i++;
  }
  int* conn= new int[nbconn];
  double* new_coords = new double[nodes.size()*_local_mesh->getSpaceDimension()]; 
  int index=0;
  for (set<int>::const_iterator iter=elems.begin(); iter!=elems.end(); iter++)
  {
    for (int inode = conn_index[*iter-1]-1; inode < conn_index[*iter]-1; inode++)
    {
      conn[index]=big2small[conn_mesh[inode]]+1;
      index++;
    } 
  }
  index=0;
  for (set<int>::const_iterator iter=nodes.begin(); iter!=nodes.end(); iter++)
  {
    int dim = _local_mesh->getSpaceDimension();
    for (int i=0; i<dim;i++)
    {
      new_coords[index]=coords[(*iter-1)*dim+i];
      index++;
    }
  }
  
  int* nbtypes=new int[nbelems_per_type.size()];
  MED_EN::medGeometryElement* new_types = new MED_EN::medGeometryElement[nbelems_per_type.size()];
  index=0;
  for (map<MED_EN::medGeometryElement,int>::const_iterator iter= nbelems_per_type.begin(); iter!=nbelems_per_type.end(); iter++)
    {
      nbtypes[index]=iter->second;
      new_types[index]=iter->first;
      index++;
    }
  MEDMEM::MESHING* meshing = new MEDMEM::MESHING();
  meshing->setCoordinates(_local_mesh->getSpaceDimension(), nodes.size(), new_coords, string("CARTESIAN"), MED_EN::MED_FULL_INTERLACE);
  meshing->setNumberOfTypes(nbelems_per_type.size(),MED_EN::MED_CELL);
  meshing->setTypes(new_types,MED_EN::MED_CELL);
  meshing->setNumberOfElements(nbtypes,MED_EN::MED_CELL);

	int dimmax=0;
 int* conn_ptr= conn;
  for (int i=0; i<nbelems_per_type.size(); i++)
  {
    meshing->setConnectivity(conn_ptr, MED_EN::MED_CELL,new_types[i]);
    conn_ptr+=nbtypes[i]*(new_types[i]%100);
		if (new_types[i]/100>dimmax) dimmax=new_types[i]/100;
  }
	meshing->setMeshDimension(dimmax);
  delete [] small2big;
  delete [] nbtypes;
  delete [] conn;
  delete [] new_coords;
  delete [] new_types;
  return meshing;
} 

}
