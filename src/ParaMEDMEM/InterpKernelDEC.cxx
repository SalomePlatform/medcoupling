// Copyright (C) 2007-2023  CEA/DEN, EDF R&D
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
  InterpKernelDEC::InterpKernelDEC():
    DisjointDEC(),
    _interpolation_matrix(0)
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
    _interpolation_matrix(0)
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
    _interpolation_matrix(0)
  {
  }

  InterpKernelDEC::~InterpKernelDEC()
  {
    release();
  }

  void InterpKernelDEC::release()
  {
    if (_interpolation_matrix != nullptr)
      delete _interpolation_matrix;
    _interpolation_matrix = nullptr;
    DisjointDEC::cleanInstance();
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
        //transferring option from InterpKernelDEC to ElementLocator   
        locator.copyOptions(*this);
        MEDCouplingPointSet* distant_mesh=0; 
        mcIdType* distant_ids=0;
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
        //transferring option from InterpKernelDEC to ElementLocator
        locator.copyOptions(*this);
        MEDCouplingPointSet* distant_mesh=0;
        mcIdType* distant_ids=0;
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
