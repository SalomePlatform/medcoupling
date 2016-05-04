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
#include "ParaMESH.hxx"
#include "DEC.hxx"
#include "InterpolationMatrix.hxx"
#include "InterpKernelDEC.hxx"
#include "ElementLocator.hxx"

namespace MEDCoupling
{  

  /*!
    \anchor InterpKernelDEC-det
    \class InterpKernelDEC

    \section InterpKernelDEC-over Overview

    The InterpKernelDEC enables the \ref InterpKerRemapGlobal "remapping" (or interpolation) of fields between
    two parallel codes.

    The projection
    methodology is based on the algorithms of %INTERP_KERNEL, that is to say, they work in a similar fashion than
    what the \ref remapper "sequential remapper" does. The following \ref discretization "projection methods"
    are supported: P0->P0 (the most common case), P1->P0, P0->P1.

    The computation is possible for 3D meshes, 2D meshes, and 3D-surface
    meshes. Dimensions must be identical for code A and code B (for instance, though it could be
    desirable, it is not yet possible to couple 3D surfaces with 2D surfaces).

    The name "InterpKernelDEC" comes from the fact that this class uses exactly the same algorithms
    as the sequential remapper. Both this class and the sequential
    \ref MEDCoupling::MEDCouplingRemapper "MEDCouplingRemapper" are built on top of the %INTERP_KERNEL
    algorithms (notably the computation of the intersection volumes).

    Among the important properties inherited from the parent abstract class \ref DisjointDEC-det "DisjointDEC",
    the two \ref MPIProcessorGroup-det "processor groups" (source and target) must have a void intersection.

    \image html NonCoincident_small.png "Transfer of a field supported by a quadrangular mesh to a triangular mesh".

    \image latex NonCoincident_small.eps "Transfer of a field supported by a quadrangular mesh to a triangular mesh"

    In the figure above we see the transfer of a field based on a quadrangular mesh to a new field supported by
    a triangular mesh. In a P0-P0 interpolation, to obtain the value on a triangle, the values on the
    quadrangles are weighted by their intersection area and summed.

    A typical use of InterpKernelDEC encompasses two distinct phases :
    - A setup phase during which the intersection volumes are computed and the communication structures are
    setup. This corresponds to calling the InterpKernelDEC::synchronize() method.
    - A running phase during which the projections are actually performed. This corresponds to the calls to
    sendData() and recvData() which actually trigger the data exchange. The data exchange are synchronous
    in the current version of the library so that recvData() and sendData() calls must be synchronized
    on code A and code B processor groups.

    The following code excerpt illustrates a typical use of the InterpKernelDEC class.

    \code
    ...
    InterpKernelDEC dec(groupA, groupB);
    dec.attachLocalField(field);
    dec.synchronize();
    if (groupA.containsMyRank())
    dec.recvData();
    else if (groupB.containsMyRank())
    dec.sendData();
    ...
    \endcode
    A \ref InterpKerRemapGlobal "remapping" of the field from the source mesh to the target mesh is performed by
    the function synchronise(), which computes the interpolation matrix.

    Computing the field on the receiving side can be expressed in terms of a matrix-vector product :
    \f$ \phi_t=W.\phi_s\f$, with \f$ \phi_t \f$ the field on the target side and \f$ \phi_s \f$ the field
    on the source side.
    When remapping a 3D surface to another 3D surface, a projection phase is necessary to match elements
    from both sides. Care must be taken when defining this projection to obtain a
    \ref InterpKerRemapGlobal "conservative remapping".

    In the P0-P0 case, this matrix is a plain rectangular matrix with coefficients equal to the
    intersection areas between triangle and quadrangles. For instance, in the above figure, the matrix
    is :

    \f[
    \begin{tabular}{|cccc|}
    0.72 & 0 & 0.2 & 0 \\
    0.46 & 0 & 0.51 & 0.03\\
    0.42 & 0.53 & 0 & 0.05\\
    0 & 0 & 0.92 & 0.05 \\
    \end{tabular}
    \f]

    \section InterpKernelDEC-options Options
    On top of the usual \ref MEDCoupling::DECOptions "DEC options", the options supported by %InterpKernelDEC objects are
    related to the underlying \ref InterpKerIntersectors "intersector class".
    All the options available in the intersector objects are
    available for the %InterpKernelDEC object. The various options available for  intersectors can
    be reviewed in \ref InterpKerIntersectors.
 
    For instance :
    \verbatim
    InterpKernelDEC dec(source_group, target_group);
    dec.attachLocalField(field);
    dec.setDoRotate(false);
    dec.setPrecision(1e-12);
    dec.synchronize();
    \endverbatim

    \warning{  Options must be set before calling the synchronize method. }
  */
  
  InterpKernelDEC::InterpKernelDEC():
    DisjointDEC(),
    _nb_distant_points(0), _distant_coords(0),
    _distant_locations(0), _interpolation_matrix(0)
  {  
  }

  /*!
    This constructor creates an InterpKernelDEC which has \a source_group as a working side 
    and  \a target_group as an idle side. All the processors will actually participate, but intersection computations will be performed on the working side during the \a synchronize() phase.
    The constructor must be called synchronously on all processors of both processor groups.
    The source group and target group MUST form a partition of all the procs within the communicator passed as 'world_comm'
    when building the group.

    \param source_group working side ProcessorGroup
    \param target_group lazy side ProcessorGroup

  */
  InterpKernelDEC::InterpKernelDEC(ProcessorGroup& source_group, ProcessorGroup& target_group):
    DisjointDEC(source_group, target_group),
    _nb_distant_points(0), _distant_coords(0),
    _distant_locations(0), _interpolation_matrix(0)
  {

  }


  /*!
   * Creates an InterpKernelDEC from a set of source procs IDs and target group IDs.
   * The difference with the ctor using groups is that the set of procs might not cover entirely MPI_COMM_WORLD
   * (a sub-communicator holding the union of source and target procs is recreated internally).
   */
  InterpKernelDEC::InterpKernelDEC(const std::set<int>& src_ids, const std::set<int>& trg_ids,
                                   const MPI_Comm& world_comm):
    DisjointDEC(src_ids,trg_ids,world_comm),
    _nb_distant_points(0), _distant_coords(0),
    _distant_locations(0), _interpolation_matrix(0)
  {
  }

  InterpKernelDEC::~InterpKernelDEC()
  {
    if (_interpolation_matrix !=0)
      delete _interpolation_matrix;
  } 

  /*! 
    \brief Synchronization process for exchanging topologies.

    This method prepares all the structures necessary for sending data from a processor group to the other. It uses the mesh
    underlying the fields that have been set with attachLocalField method.
    It works in four steps :
    -# Bounding boxes are computed for each sub-domain,
    -# The lazy side mesh parts that are likely to intersect the working side local processor are sent to the working side,
    -# The working side calls the interpolation kernel to compute the intersection between local and imported mesh.
    -# The lazy side is updated so that it knows the structure of the data that will be sent by
    the working side during a \a sendData() call.

  */
  void InterpKernelDEC::synchronize()
  {
    if(!isInUnion())
      return ;
    delete _interpolation_matrix;
    _interpolation_matrix = new InterpolationMatrix (_local_field, *_source_group,*_target_group,*this,*this); 

    //setting up the communication DEC on both sides  
    if (_source_group->containsMyRank())
      {
        //locate the distant meshes
        ElementLocator locator(*_local_field, *_target_group, *_source_group);
        //transfering option from InterpKernelDEC to ElementLocator   
        locator.copyOptions(*this);
        MEDCouplingPointSet* distant_mesh=0; 
        int* distant_ids=0;
        std::string distantMeth;
        for (int i=0; i<_target_group->size(); i++)
          {
            //        int idistant_proc = (i+_source_group->myRank())%_target_group->size();
            int idistant_proc=i;

            //gathers pieces of the target meshes that can intersect the local mesh
            locator.exchangeMesh(idistant_proc,distant_mesh,distant_ids);
            if (distant_mesh !=0)
              {
                locator.exchangeMethod(_method,idistant_proc,distantMeth);
                //adds the contribution of the distant mesh on the local one
                int idistant_proc_in_union=_union_group->translateRank(_target_group,idistant_proc);
                //std::cout <<"add contribution from proc "<<idistant_proc_in_union<<" to proc "<<_union_group->myRank()<<std::endl;
                _interpolation_matrix->addContribution(*distant_mesh,idistant_proc_in_union,distant_ids,_method,distantMeth);
                distant_mesh->decrRef();
                delete [] distant_ids;
                distant_mesh=0;
                distant_ids=0;
              }
          }
       _interpolation_matrix->finishContributionW(locator);
      }

    if (_target_group->containsMyRank())
      {
        ElementLocator locator(*_local_field, *_source_group, *_target_group);
        //transfering option from InterpKernelDEC to ElementLocator
        locator.copyOptions(*this);
        MEDCouplingPointSet* distant_mesh=0;
        int* distant_ids=0;
        for (int i=0; i<_source_group->size(); i++)
          {
            //        int idistant_proc = (i+_target_group->myRank())%_source_group->size();
            int  idistant_proc=i;
            //gathers pieces of the target meshes that can intersect the local mesh
            locator.exchangeMesh(idistant_proc,distant_mesh,distant_ids);
            //std::cout << " Data sent from "<<_union_group->myRank()<<" to source proc "<< idistant_proc<<std::endl;
            if (distant_mesh!=0)
              {
                std::string distantMeth;
                locator.exchangeMethod(_method,idistant_proc,distantMeth);
                distant_mesh->decrRef();
                delete [] distant_ids;
                distant_mesh=0;
                distant_ids=0;
              }
          }
        _interpolation_matrix->finishContributionL(locator);
      }
    _interpolation_matrix->prepare();
  }


  /*!
    Receives the data whether the processor is on the working side or on the lazy side. It must match a \a sendData() call on the other side.
  */
  void InterpKernelDEC::recvData()
  {
    if (_source_group->containsMyRank())
      _interpolation_matrix->transposeMultiply(*_local_field->getField());
    else if (_target_group->containsMyRank())
      {
        _interpolation_matrix->multiply(*_local_field->getField());
        if (getForcedRenormalization())
          renormalizeTargetField(getMeasureAbsStatus());
      }
  }


  /*!
    Receives the data at time \a time in asynchronous mode. The value of the field
    will be time-interpolated from the field values received.
    \param time time at which the value is desired
  */
  void InterpKernelDEC::recvData( double time )
  {
    _interpolation_matrix->getAccessDEC()->setTime(time);
    recvData() ;
  }

  /*!
    Sends the data whether the processor is on the working side or on the lazy side.
    It must match a recvData() call on the other side.
  */
  void InterpKernelDEC::sendData()
  {
    if (_source_group->containsMyRank())
      {
    
        _interpolation_matrix->multiply(*_local_field->getField());
        if (getForcedRenormalization())
          renormalizeTargetField(getMeasureAbsStatus());
    
      }
    else if (_target_group->containsMyRank())
      _interpolation_matrix->transposeMultiply(*_local_field->getField());
  }

  /*!
    Sends the data available at time \a time in asynchronous mode. 
    \param time time at which the value is available
    \param deltatime time interval between the value presently sent and the next one. 
  */
  void InterpKernelDEC::sendData( double time , double deltatime )
  {
    _interpolation_matrix->getAccessDEC()->setTime(time,deltatime);
    sendData() ;
  }
  
}
