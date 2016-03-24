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
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"
#include "ParaFIELD.hxx"
#include "MPIProcessorGroup.hxx"
#include "DEC.hxx"
#include "NonCoincidentDEC.hxx"

extern "C" {
#include <fvm_parall.h>
#include <fvm_nodal.h>
#include <fvm_nodal_append.h>
#include <fvm_locator.h>
}

namespace MEDCoupling
{

  /*!
    \anchor NonCoincidentDEC-det
    \class NonCoincidentDEC

    \c NonCoincidentDEC enables non-conservative remapping of fields
    between two parallel codes. 
    The computation is possible for 3D meshes and 2D meshes.
    It is not available for 3D surfaces.

    The computation enables fast parallel localization, and is based on a point in element search, followed
    by a field evaluation at the point location. Thus, it is typically
    faster than the \ref InterpKernelDEC-det "InterpKernelDEC" which uses a
    \ref InterpKerRemapGlobal "conservative remapping" (i.e. the same algorithms of volume
    intersection as in the \ref remapper "sequential remapper")
    It is particularly true for the initialisation phase (synchronize() method)
    which has a significant computation cost in \ref InterpKernelDEC-det.

    In the present version, only fields lying on elements ("P0") are considered.
    The value is estimated by locating the barycenter of the target
    side cell in a source cell and sending the value of this source cell 
    as the value of the target cell.

    \image html NonCoincident_small.png "Example showing the transfer from a field based on a quadrangular mesh to a triangular mesh. The triangle barycenters are computed and located in the quadrangles. In a P0-P0 interpolation, the value on the quadrangle is then applied to the triangles whose barycenter lies within."

    \image latex NonCoincident_small.eps "Example showing the transfer from a field based on a quadrangular mesh to a triangular mesh. The triangle barycenters are computed and located in the quadrangles. In a P0-P0 interpolation, the value on the quadrangle is then applied to the triangles whose barycenter lies within."

    A typical use of NonCoincidentDEC encompasses two distinct phases :
    - A setup phase during which the intersection volumes are computed and the communication structures are setup. This corresponds to calling the NonCoincidentDEC::synchronize() method.
    - A use phase during which the remappings are actually performed. This corresponds to the calls to sendData() and recvData() which actually trigger the data exchange. The data exchange are synchronous in the current version of the library so that recvData() and sendData() calls must be synchronized on code A and code B processor groups. 

    The following code excerpt illutrates a typical use of the NonCoincidentDEC class.

    \code
    ...
    NonCoincidentDEC dec(groupA, groupB);
    dec.attachLocalField(field);
    dec.synchronize();
    if (groupA.containsMyRank())
    dec.recvData();
    else if (groupB.containsMyRank())
    dec.sendData();
    ...
    \endcode

    Computing the field on the receiving side can be expressed in terms 
    of a matrix-vector product : \f$ \phi_t=W.\phi_s\f$, with \f$ \phi_t
    \f$ the field on the target side and \f$ \phi_s \f$ the field on 
    the source side.
    In the P0-P0 case, this matrix is a plain rectangular matrix with one 
    non-zero element per row (with value 1). For instance, in the above figure, the matrix is :
    \f[

    \begin{tabular}{|cccc|}
    1 & 0 & 0 & 0\\
    0 & 0 & 1 & 0\\
    1 & 0 & 0 & 0\\
    0 & 0 & 1 & 0\\
    \end{tabular}
    \f]
  */

  fvm_nodal_t*  medmemMeshToFVMMesh(const MEDMEM::MESH* mesh)
  {
    // create an FVM structure from the paramesh structure
    std::string meshName(mesh->getName());//this line avoid that mesh->getName() object killed before fvm_nodal_create read the const char *.
    fvm_nodal_t * fvm_nodal = fvm_nodal_create(meshName.c_str(),mesh->getMeshDimension());
      
    //loop on cell types
    int nbtypes = mesh->getNumberOfTypes(MED_EN::MED_CELL);
    const MED_EN::medGeometryElement* types = mesh->getTypes(MED_EN::MED_CELL);
    for (int itype=0; itype<nbtypes; itype++)
      {
        fvm_element_t fvm_type;
        switch (types[itype]) 
          {
          case MED_EN::MED_TRIA3 :
            fvm_type=FVM_FACE_TRIA;
            break;
          case MED_EN::MED_QUAD4 :
            fvm_type=FVM_FACE_QUAD;
            break;
          case MED_EN::MED_TETRA4 :
            fvm_type=FVM_CELL_TETRA;
            break;
          case MED_EN::MED_HEXA8 :
            fvm_type=FVM_CELL_HEXA;
            break;
          default:
            throw MEDEXCEPTION(" MED type  conversion to fvm is not handled yet.");
            break;

          }

        fvm_lnum_t nbelems = mesh->getNumberOfElements(MED_EN::MED_CELL, types[itype]);
        fvm_lnum_t* conn = new fvm_lnum_t[nbelems*(types[itype]%100)];
        const int* mesh_conn =mesh->getConnectivity(MED_EN::MED_FULL_INTERLACE,MED_EN::MED_NODAL, MED_EN::MED_CELL, types[itype]);
        for (int i=0; i<nbelems*(types[itype]%100); i++)
          conn[i]=mesh_conn[i]; 
        //swapping trias
        if (types[itype]==MED_EN::MED_TRIA3)
          {
            for (int i=0; i<nbelems;i++)
              {
                int tmp=conn[3*i];
                conn[3*i]=mesh_conn[3*i+1];
                conn[3*i+1]=tmp;
              }
          }
        //swapping tetras
        if (types[itype]==MED_EN::MED_TETRA4)
          {
            for (int i=0; i<nbelems;i++)
              {
                int tmp=conn[4*i];
                conn[4*i]=mesh_conn[4*i+1];
                conn[4*i+1]=tmp;
              }
          }
        fvm_nodal_append_by_transfer(fvm_nodal, nbelems, fvm_type,0,0,0,conn,0);
         
        int nbnodes= mesh->getNumberOfNodes();
        int spacedim=mesh->getSpaceDimension();
        fvm_coord_t* coords = new fvm_coord_t[nbnodes*spacedim];
        const double* mesh_coords=mesh->getCoordinates(MED_EN::MED_FULL_INTERLACE);
        for (int i=0; i<nbnodes*spacedim; i++)
          coords[i]=mesh_coords[i];                  
        fvm_nodal_transfer_vertices(fvm_nodal,coords);
      }
    return fvm_nodal;
  }
  
  fvm_nodal_t*  medmemSupportToFVMMesh(const MEDMEM::SUPPORT* support)
  {

    // create an FVM structure from the paramesh structure
    std::string supportName(support->getName());//this line avoid that support->getName() object killed before fvm_nodal_create read the const char *.
    fvm_nodal_t * fvm_nodal = fvm_nodal_create(supportName.c_str(),1);
      
    const MEDMEM::MESH* mesh= support->getMesh();
      
    //loop on cell types
    MED_EN::medEntityMesh entity = support->getEntity();
      
    int nbtypes = support->getNumberOfTypes();
    const MED_EN::medGeometryElement* types = support->getTypes();
    int ioffset=0;
    const int* type_offset = support->getNumberIndex();
      
    //browsing through all types
    for (int itype=0; itype<nbtypes; itype++)
      {
        fvm_element_t fvm_type;
        switch (types[itype]) 
          {
          case MED_EN::MED_TRIA3 :
            fvm_type=FVM_FACE_TRIA;
            break;
          case MED_EN::MED_QUAD4 :
            fvm_type=FVM_FACE_QUAD;
            break;
          case MED_EN::MED_TETRA4 :
            fvm_type=FVM_CELL_TETRA;
            break;
          case MED_EN::MED_HEXA8 :
            fvm_type=FVM_CELL_HEXA;
            break;
          default:
            throw MEDEXCEPTION(" MED type  conversion to fvm is not handled yet.");
            break;

          }
        fvm_lnum_t nbelems = support->getNumberOfElements(types[itype]);
         
        //for a partial support, defining the element numbers that are taken into
        //account in the support
        fvm_lnum_t* elem_numbers=0;
        if (!support->isOnAllElements())
          {
            elem_numbers = const_cast<fvm_lnum_t*> (support->getNumber(types[itype]));
           
            //creating work arrays to store list of elems for partial suports
            if (itype>0)
              {
                fvm_lnum_t* temp = new int[nbelems];
                for (int i=0; i< nbelems; i++)
                  temp[i] = elem_numbers [i]-ioffset;
                ioffset+=type_offset[itype];
                elem_numbers = temp;
              }
          }
        //retrieving original mesh connectivity
        fvm_lnum_t* conn = const_cast<fvm_lnum_t*> (mesh->getConnectivity(MED_EN::MED_FULL_INTERLACE,MED_EN::MED_NODAL,entity, types[itype]));
       
        // adding the elements to the FVM structure 
        fvm_nodal_append_by_transfer(fvm_nodal, nbelems, fvm_type,0,0,0,conn,elem_numbers);
         
        //cleaning work arrays (for partial supports)
        if (!support->isOnAllElements() && itype>0)
          delete[] elem_numbers;
      
      }
    return fvm_nodal;
  }
  
  NonCoincidentDEC::NonCoincidentDEC()
  {  
  }

  /*! Constructor of a non coincident \ref para-dec "DEC" with
   * a source group on which lies a field lying on a mesh and a 
   * target group on which lies a mesh.
   * 
   * \param source_group ProcessorGroup on the source side
   * \param target_group ProcessorGroup on the target side 
   */
  
  NonCoincidentDEC::NonCoincidentDEC(ProcessorGroup& source_group,
                                     ProcessorGroup& target_group)
    :DEC(source_group, target_group)
  {}
                                   
  NonCoincidentDEC::~NonCoincidentDEC()
  {
  }

  /*! Synchronization process. Calling this method 
   * synchronizes the topologies so that the target side
   * gets the information which enable it to fetch the field value 
   * from the source side.
   * A typical call is : 
   \verbatim
   NonCoincidentDEC dec(source_group,target_group);
   dec.attachLocalField(field);
   dec.synchronize();
   \endverbatim
  */
  void NonCoincidentDEC::synchronize()
  {
  
    //initializing FVM parallel environment
    const MPI_Comm* comm=dynamic_cast<const MPIProcessorGroup*> (_union_group)->getComm();
    fvm_parall_set_mpi_comm(*const_cast<MPI_Comm*> (comm));
  
  
    //setting up the communication DEC on both sides
  
    if (_source_group->containsMyRank())
      {
        MEDMEM::MESH* mesh = _local_field->getField()->getSupport()->getMesh();
        fvm_nodal_t* source_nodal = MEDCoupling::medmemMeshToFVMMesh(mesh);
      
        int target_size = _target_group->size()  ;
        int start_rank=  _source_group->size();
        const MPI_Comm* comm = (dynamic_cast<const MPIProcessorGroup*> (_union_group))->getComm();
      
        _locator =  fvm_locator_create(1e-6,
                                       *comm,
                                       target_size,
                                       start_rank);
      
        fvm_locator_set_nodal(_locator,
                              source_nodal,
                              mesh->getSpaceDimension(),
                              0,
                              NULL,
                              0);

      
        _nb_distant_points = fvm_locator_get_n_dist_points(_locator);
        _distant_coords = fvm_locator_get_dist_coords(_locator);
        _distant_locations = fvm_locator_get_dist_locations(_locator);
           
      }
    if (_target_group->containsMyRank())
      {
        MEDMEM::MESH* mesh = _local_field->getField()->getSupport()->getMesh();
      
        fvm_nodal_t* target_nodal = MEDCoupling::medmemMeshToFVMMesh(mesh);
        int source_size = _source_group->size();
        int start_rank=  0 ;
        const MPI_Comm* comm = (dynamic_cast<const MPIProcessorGroup*> (_union_group))->getComm();
      
        _locator = fvm_locator_create(1e-6,
                                      *comm,
                                      source_size,
                                      start_rank);
        int nbcells = mesh->getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);
        const MEDMEM::SUPPORT* support=_local_field->getField()->getSupport();
        MEDMEM::FIELD<double>* barycenter_coords = mesh->getBarycenter(support);
        const double* coords = barycenter_coords->getValue();
        fvm_locator_set_nodal(_locator,
                              target_nodal,
                              mesh->getSpaceDimension(),
                              nbcells,
                              NULL,
                              coords);  
        delete barycenter_coords;
      }
  }


  /*! This method is called on the target group in order to 
   * trigger the retrieveal of field data. It must 
   * be called synchronously with a sendData() call on 
   * the source group.
   */
  void NonCoincidentDEC::recvData()
  {
    int nbelems = _local_field->getField()->getSupport()->getMesh()->getNumberOfElements(MED_EN::MED_CELL, MED_EN::MED_ALL_ELEMENTS);
    int nbcomp =  _local_field->getField()->getNumberOfComponents();
    double* values = new double [nbelems*nbcomp];
    fvm_locator_exchange_point_var(_locator,
                                   0,
                                   values,
                                   0,
                                   sizeof(double),
                                   nbcomp,
                                   0);
    _local_field->getField()->setValue(values);
    if (_forced_renormalization_flag)
      renormalizeTargetField();
    delete[]values;
  }

  /*! This method is called on the source group in order to 
   * send field data. It must be called synchronously with 
   * a recvData() call on 
   * the target group.
   */
  void NonCoincidentDEC::sendData()
  {
    const double* values=_local_field->getField()->getValue();
    int nbcomp = _local_field->getField()->getNumberOfComponents();
    double* distant_values = new double [_nb_distant_points*nbcomp];

    //cheap interpolation :  the value of the cell is transfered to the point
    for (int i=0; i<_nb_distant_points; i++)
      for (int j=0; j <nbcomp; j++)
        distant_values[i*nbcomp+j]=values[(_distant_locations[i]-1)*nbcomp+j];
  
    fvm_locator_exchange_point_var(_locator,
                                   distant_values,
                                   0,
                                   0,
                                   sizeof(double),
                                   nbcomp,
                                   0);

    delete [] distant_values;
    if (_forced_renormalization_flag)
      renormalizeTargetField();

  }
}
