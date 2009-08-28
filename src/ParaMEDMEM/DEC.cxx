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
#include "CommInterface.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"
#include "ParaFIELD.hxx"
#include "DEC.hxx"
#include "ICoCoField.hxx"
#include "ICoCoMEDField.hxx"
#include "ICoCoTrioField.hxx"
#include "MPIProcessorGroup.hxx"

#include <cmath>

/*! \defgroup dec DEC
 *
 * \section decintroduction Introduction
 *
 * Interface class for creation of a link between two 
 * processor groups for exhanging mesh or field data.
 * The DEC is defined by attaching a field on the receiving or on the 
 * sending side. 
 * On top of attaching a ParaMEDMEM::FIELD, it is possible to
 * attach a ICoCo::Field. This class is an abstract class that enables 
 * coupling of codes that respect the ICoCo interface \ref icoco. It has two implementations:
 * one for codes that express their fields as MEDCoupling fields (ICoCo::MEDField) and one
 * for codes that express their fields as Trio/U fields.
 * 
 * \section dec_options DEC Options
 * Options supported by DEC objects are
 *
 * <TABLE BORDER=1 >
 * <TR><TD>Option</TD><TD>Description</TD><TD>Default value</TD></TR>
 * <TR><TD>ForcedRenormalization</TD><TD>After receiving data, the target field is renormalized so that L2-norms of the source and target fields match.</TD><TD> false </TD></TR>
 *</TABLE>


 The following code excerpt shows how to set options for an object that inherits from DEC :

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
  DEC::DEC(ProcessorGroup& source_group, ProcessorGroup& target_group):_local_field(0), 
                                                                       _source_group(&source_group),
                                                                       _target_group(&target_group),
                                                                       _owns_field(false),
                                                                       _icoco_field(0)
  {
    _union_group = source_group.fuse(target_group);  
  }

  DEC::~DEC()
  {
    //    delete _union_group;
    if(_owns_field)
      delete _local_field;
    delete _icoco_field;
    delete _union_group;
  }

  void DEC::setNature(NatureOfField nature)
  {
    if(_local_field)
      _local_field->getField()->setNature(nature);
  }

  /*! Attaches a local field to a DEC.
    If the processor is on the receiving end of the DEC, the field
    will be updated by a recvData() call.
    Reversely, if the processor is on the sending end, the field will be read, possibly transformed, and sent appropriately to the other side.
  */
  void DEC::attachLocalField(const ParaFIELD* field) 
  {
    _local_field=field;
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

  void DEC::attachLocalField(MEDCouplingFieldDouble* field) 
  {
    ProcessorGroup* local_group;
    if (_source_group->containsMyRank())
      local_group=_source_group;
    else if (_target_group->containsMyRank())
      local_group=_target_group;
    else
      throw INTERP_KERNEL::Exception("Invalid procgroup for field attachment to DEC");
        
    _local_field= new ParaFIELD(field, *local_group);
    _owns_field=true;
    _comm_interface=&(local_group->getCommInterface());
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
  
  void DEC::attachLocalField(const ICoCo::Field* field)
  {
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
        
        _icoco_field=new ICoCo::MEDField(*const_cast<ICoCo::TrioField* >(triofield), *localgroup);
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
  void DEC::renormalizeTargetField(bool isWAbs)
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

  void DEC::compareFieldAndMethod() const throw(INTERP_KERNEL::Exception)
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
}
