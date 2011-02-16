//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#include "OverlapDEC.hxx"
#include "CommInterface.hxx"
#include "ParaFIELD.hxx"
#include "MPIProcessorGroup.hxx"
#include "OverlapElementLocator.hxx"
#include "OverlapInterpolationMatrix.hxx"

namespace ParaMEDMEM
{
  OverlapDEC::OverlapDEC(const std::set<int>& procIds, const MPI_Comm& world_comm):_own_group(true),_interpolation_matrix(0),
                                                                                   _source_field(0),_own_source_field(false),
                                                                                   _target_field(0),_own_target_field(false)
  {
    ParaMEDMEM::CommInterface comm;
    int *ranks_world=new int[procIds.size()]; // ranks of sources and targets in world_comm
    std::copy(procIds.begin(),procIds.end(),ranks_world);
    MPI_Group group,world_group;
    comm.commGroup(world_comm,&world_group);
    comm.groupIncl(world_group,procIds.size(),ranks_world,&group);
    delete [] ranks_world;
    MPI_Comm theComm;
    comm.commCreate(world_comm,group,&theComm);
    if(theComm==MPI_COMM_NULL)
      {
        _group=0;
        return ;
      }
    std::set<int> idsUnion;
    for(std::size_t i=0;i<procIds.size();i++)
      idsUnion.insert(i);
    _group=new MPIProcessorGroup(comm,idsUnion,theComm);
  }

  OverlapDEC::~OverlapDEC()
  {
    if(_own_group)
      delete _group;
    if(_own_source_field)
      delete _source_field;
    if(_own_target_field)
      delete _target_field;
  }

  void OverlapDEC::sendRecvData(bool way)
  {
  }
  
  void OverlapDEC::synchronize()
  {
    if(!isInGroup())
      return ;
    delete _interpolation_matrix;
    _interpolation_matrix=new OverlapInterpolationMatrix(_source_field,_target_field,*_group,*this,*this);
    OverlapElementLocator locator(_source_field,_target_field,*_group);
  }

  void OverlapDEC::attachSourceLocalField(ParaFIELD *field, bool ownPt)
  {
    if(!isInGroup())
      return ;
    if(_own_source_field)
      delete _source_field;
    _source_field=field;
    _own_source_field=ownPt;
  }

  void OverlapDEC::attachTargetLocalField(ParaFIELD *field, bool ownPt)
  {
    if(!isInGroup())
      return ;
    if(_own_target_field)
      delete _target_field;
    _target_field=field;
    _own_target_field=ownPt;
  }

  bool OverlapDEC::isInGroup() const
  {
    if(!_group)
      return false;
    return _group->containsMyRank();
  }
}
