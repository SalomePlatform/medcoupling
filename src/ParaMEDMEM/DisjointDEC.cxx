// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#include "DisjointDEC.hxx"
#include "CommInterface.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"
#include "ParaFIELD.hxx"
#include "ParaMESH.hxx"
#include "ICoCoField.hxx"
#include "ICoCoMEDField.hxx"
#include "ICoCoTrioField.hxx"
#include "MPIProcessorGroup.hxx"

#include <cmath>
#include <iostream>

/*! \defgroup dec DEC
 *
 * \section decintroduction Introduction
 *
 * Interface class for creation of a link between two 
 * processor groups for exhanging mesh or field data.
 * The \c DEC is defined by attaching a field on the receiving or on the 
 * sending side. 
 * On top of attaching a \c ParaMEDMEM::FIELD, it is possible to
 * attach a ICoCo::Field. This class is an abstract class that enables 
 * coupling of codes that respect the ICoCo interface \ref icoco. It has two implementations:
 * one for codes that express their fields as \ref medoupling fields (ICoCo::MEDField) and one
 * for codes that express their fields as Trio/U fields.
 * 
 * \section dec_options DEC Options
 * Options supported by \c DEC objects are
 *
 * <TABLE BORDER=1 >
 * <TR><TD>Option</TD><TD>Description</TD><TD>Default value</TD></TR>
 * <TR><TD>ForcedRenormalization</TD><TD>After receiving data, the target field is renormalized so that L2-norms of the source and target fields match.</TD><TD> false </TD></TR>
 *</TABLE>


 The following code excerpt shows how to set options for an object that inherits from \c DEC :

 \code
 InterpKernelDEC dec(source_group,target_group);
 dec.setOptions("ForcedRenormalization",true);
 dec.attachLocalField(field);
 dec.synchronize();
 if (source_group.containsMyRank())
 dec.sendData();
 else
 dec.recvData();
 \endcode
*/

namespace ParaMEDMEM
{


  /*! \addtogroup dec
    @{ 
  */
  DisjointDEC::DisjointDEC(ProcessorGroup& source_group, ProcessorGroup& target_group):_local_field(0), 
                                                                                       _source_group(&source_group),
                                                                                       _target_group(&target_group),
                                                                                       _owns_field(false),
                                                                                       _owns_groups(false),
                                                                                       _icoco_field(0)
  {
    _union_group = source_group.fuse(target_group);  
  }
  
  DisjointDEC::DisjointDEC(const DisjointDEC& s):DEC(s),_local_field(0),_union_group(0),_source_group(0),_target_group(0),_owns_field(false),_owns_groups(false),_icoco_field(0)
  {
    copyInstance(s);
  }

  DisjointDEC & DisjointDEC::operator=(const DisjointDEC& s)
  {
    cleanInstance();
    copyInstance(s);
    return *this;
     
  }

  void DisjointDEC::copyInstance(const DisjointDEC& other)
  {
    DEC::copyFrom(other);
    if(other._target_group)
      {
        _target_group=other._target_group->deepCpy();
        _owns_groups=true;
      }
    if(other._source_group)
      {
        _source_group=other._source_group->deepCpy();
        _owns_groups=true;
      }
    if (_source_group && _target_group)
      _union_group = _source_group->fuse(*_target_group);
  }

  DisjointDEC::DisjointDEC(const std::set<int>& source_ids, const std::set<int>& target_ids, const MPI_Comm& world_comm):_local_field(0), 
                                                                                                                         _owns_field(false),
                                                                                                                         _owns_groups(true),
                                                                                                                         _icoco_field(0)
  {
    ParaMEDMEM::CommInterface comm;
    // Create the list of procs including source and target
    std::set<int> union_ids; // source and target ids in world_comm
    union_ids.insert(source_ids.begin(),source_ids.end());
    union_ids.insert(target_ids.begin(),target_ids.end());
    if(union_ids.size()!=(source_ids.size()+target_ids.size()))
      throw INTERP_KERNEL::Exception("DisjointDEC constructor : source_ids and target_ids overlap partially or fully. This type of DEC does not support it ! OverlapDEC class could be the solution !");
    int* union_ranks_world=new int[union_ids.size()]; // ranks of sources and targets in world_comm
    std::copy(union_ids.begin(), union_ids.end(), union_ranks_world);

    // Create a communicator on these procs
    MPI_Group union_group,world_group;
    comm.commGroup(world_comm,&world_group);
    comm.groupIncl(world_group,union_ids.size(),union_ranks_world,&union_group);
    MPI_Comm union_comm;
    comm.commCreate(world_comm,union_group,&union_comm);
    delete[] union_ranks_world;

    if (union_comm==MPI_COMM_NULL)
      { // This process is not in union
        _source_group=0;
        _target_group=0;
        _union_group=0;
        return;
      }

    // Translate source_ids and target_ids from world_comm to union_comm
    int* source_ranks_world=new int[source_ids.size()]; // ranks of sources in world_comm
    std::copy(source_ids.begin(), source_ids.end(),source_ranks_world);
    int* source_ranks_union=new int[source_ids.size()]; // ranks of sources in union_comm
    int* target_ranks_world=new int[target_ids.size()]; // ranks of targets in world_comm
    std::copy(target_ids.begin(), target_ids.end(),target_ranks_world);
    int* target_ranks_union=new int[target_ids.size()]; // ranks of targets in union_comm
    MPI_Group_translate_ranks(world_group,source_ids.size(),source_ranks_world,union_group,source_ranks_union);
    MPI_Group_translate_ranks(world_group,target_ids.size(),target_ranks_world,union_group,target_ranks_union);
    std::set<int> source_ids_union;
    for (int i=0;i<(int)source_ids.size();i++)
      source_ids_union.insert(source_ranks_union[i]);
    std::set<int> target_ids_union;
    for (int i=0;i<(int)target_ids.size();i++)
      target_ids_union.insert(target_ranks_union[i]);
    delete [] source_ranks_world;
    delete [] source_ranks_union;
    delete [] target_ranks_world;
    delete [] target_ranks_union;

    // Create the MPIProcessorGroups
    _source_group = new MPIProcessorGroup(comm,source_ids_union,union_comm);
    _target_group = new MPIProcessorGroup(comm,target_ids_union,union_comm);
    _union_group = _source_group->fuse(*_target_group);

  }

  DisjointDEC::~DisjointDEC()
  {
    cleanInstance();
  }

  void DisjointDEC::cleanInstance()
  {
    if(_owns_field)
      {
        delete _local_field;
      }
    _local_field=0;
    _owns_field=false;
    if(_owns_groups)
      {
        delete _source_group;
        delete _target_group;
      }
    _owns_groups=false;
    _source_group=0;
    _target_group=0;
    delete _icoco_field;
    _icoco_field=0;
    delete _union_group;
    _union_group=0;
  }

  void DisjointDEC::setNature(NatureOfField nature)
  {
    if(_local_field)
      _local_field->getField()->setNature(nature);
  }

  /*! Attaches a local field to a DEC.
    If the processor is on the receiving end of the DEC, the field
    will be updated by a recvData() call.
    Reversely, if the processor is on the sending end, the field will be read, possibly transformed, and sent appropriately to the other side.
  */
  void DisjointDEC::attachLocalField(const ParaFIELD* field, bool ownPt) 
  {
    if(!isInUnion())
      return ;
    if(_owns_field)
      delete _local_field;
    _local_field=field;
    _owns_field=ownPt;
    _comm_interface=&(field->getTopology()->getProcGroup()->getCommInterface());
    compareFieldAndMethod();
  }

  /*! Attaches a local field to a DEC. The method will test whether the processor
    is on the source or the target side and will associate the mesh underlying the 
    field to the local side.

    If the processor is on the receiving end of the DEC, the field
    will be updated by a recvData() call.
    Reversely, if the processor is on the sending end, the field will be read, possibly transformed,
    and sent appropriately to the other side.
  */

  void DisjointDEC::attachLocalField(MEDCouplingFieldDouble* field) 
  {
    if(!isInUnion())
      return ;
    ProcessorGroup* local_group;
    if (_source_group->containsMyRank())
      local_group=_source_group;
    else if (_target_group->containsMyRank())
      local_group=_target_group;
    else
      throw INTERP_KERNEL::Exception("Invalid procgroup for field attachment to DEC");
    ParaMESH *paramesh=new ParaMESH((MEDCouplingPointSet *)field->getMesh(),*local_group,field->getMesh()->getName());
    ParaFIELD *tmp=new ParaFIELD(field, paramesh, *local_group);
    tmp->setOwnSupport(true);
    attachLocalField(tmp,true);
    //_comm_interface=&(local_group->getCommInterface());
  }

  /*! 
    Attaches a local field to a DEC.
    If the processor is on the receiving end of the DEC, the field
    will be updated by a recvData() call.
    Reversely, if the processor is on the sending end, the field will be read, possibly transformed, and sent appropriately to the other side.
    The field type is a generic ICoCo Field, so that the DEC can couple a number of different fields :
    - a ICoCo::MEDField, that is created from a MEDCoupling structure
    - a ICOCo::TrioField, that is created from tables extracted from a TRIO-U structure.
    
  */
  void DisjointDEC::attachLocalField(const ICoCo::Field* field)
  {
    if(!isInUnion())
      return ;
    const ICoCo::MEDField* medfield=dynamic_cast<const ICoCo::MEDField*> (field);
    if(medfield !=0)
      {
        attachLocalField(medfield->getField());
        return;
      }
    const ICoCo::TrioField* triofield=dynamic_cast<const ICoCo::TrioField*> (field);
    if (triofield !=0)
      {
        ProcessorGroup* localgroup;
        if (_source_group->containsMyRank())
          localgroup=_source_group;
        else
          localgroup=_target_group;
        delete _icoco_field;
        
        _icoco_field=new ICoCo::MEDField(*const_cast<ICoCo::TrioField* >(triofield));
        attachLocalField(_icoco_field);
        return;
      }
    throw INTERP_KERNEL::Exception("incompatible field type");
  }
  
  /*!
    Computes the field norm over its support 
    on the source side and renormalizes the field on the target side
    so that the norms match.

    \f[
    I_{source}=\sum_{i=1}^{n_{source}}V_{i}.|\Phi^{source}_{i}|^2,
    \f]

    \f[
    I_{target}=\sum_{i=1}^{n_{target}}V_{i}.|\Phi^{target}_{i}|^2,
    \f]
    
    \f[
    \Phi^{target}:=\Phi^{target}.\sqrt{I_{source}/I_{target}}.
    \f]

  */
  void DisjointDEC::renormalizeTargetField(bool isWAbs)
  {
    if (_source_group->containsMyRank())
      for (int icomp=0; icomp<_local_field->getField()->getArray()->getNumberOfComponents(); icomp++)
        {
          double total_norm = _local_field->getVolumeIntegral(icomp+1,isWAbs);
          double source_norm = total_norm;
          _comm_interface->broadcast(&source_norm, 1, MPI_DOUBLE, 0,* dynamic_cast<MPIProcessorGroup*>(_union_group)->getComm());

        }
    if (_target_group->containsMyRank())
      {
        for (int icomp=0; icomp<_local_field->getField()->getArray()->getNumberOfComponents(); icomp++)
          {
            double total_norm = _local_field->getVolumeIntegral(icomp+1,isWAbs);
            double source_norm=total_norm;
            _comm_interface->broadcast(&source_norm, 1, MPI_DOUBLE, 0,* dynamic_cast<MPIProcessorGroup*>(_union_group)->getComm());

            if (fabs(total_norm)>1e-100)
              _local_field->getField()->applyLin(source_norm/total_norm,0.0,icomp+1);
          }
      }
  }
  /*! @} */

  bool DisjointDEC::isInSourceSide() const
  {
    if(!_source_group)
      return false;
    return _source_group->containsMyRank();
  }

  bool DisjointDEC::isInTargetSide() const
  {
    if(!_target_group)
      return false;
    return _target_group->containsMyRank();
  }

  bool DisjointDEC::isInUnion() const
  {
    if(!_union_group)
      return false;
    return _union_group->containsMyRank();
  }

  void DisjointDEC::compareFieldAndMethod() const throw(INTERP_KERNEL::Exception)
  {
    if (_local_field)
      {
        TypeOfField entity = _local_field->getField()->getTypeOfField();
        if ( getMethod() == "P0" )
          {
            if ( entity != ON_CELLS )
              throw INTERP_KERNEL::Exception("Field support and interpolation method mismatch."
                                             " For P0 interpolation, field must be on MED_CELL's");
          }
        else if ( getMethod() == "P1" )
          {
            if ( entity != ON_NODES )
              throw INTERP_KERNEL::Exception("Field support and interpolation method mismatch."
                                             " For P1 interpolation, field must be on MED_NODE's");
          }
        else if ( getMethod() == "P1d" )
          {
            if ( entity != ON_CELLS )
              throw INTERP_KERNEL::Exception("Field support and interpolation method mismatch."
                                             " For P1d interpolation, field must be on MED_CELL's");
            if ( _target_group->containsMyRank() )
              throw INTERP_KERNEL::Exception("Projection to P1d field not supported");
          }
        else
          {
            throw INTERP_KERNEL::Exception("Unknown interpolation method. Possible methods: P0, P1, P1d");
          }
      }
  }

  /*!
    If way==true, source procs call sendData() and target procs call recvData().
    if way==false, it's the other way round.
  */
  void DisjointDEC::sendRecvData(bool way)
  {
    if(!isInUnion())
      return;
    if(isInSourceSide())
      {
        if(way)
          sendData();
        else
          recvData();
      }
    else if(isInTargetSide())
      {
        if(way)
          recvData();
        else
          sendData();
      }
  }
}
