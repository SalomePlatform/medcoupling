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
#ifndef __DEC_HXX__
#define __DEC_HXX__

#include "MEDCouplingFieldDouble.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "DECOptions.hxx"

namespace ICoCo
{
  class Field;
}

namespace ParaMEDMEM
{
  class ProcessorGroup;
  class ParaFIELD;
  class CommInterface;
  class DEC : public DECOptions
  {
  public:
    DEC():_local_field(0) { }
    DEC(ProcessorGroup& source_group, ProcessorGroup& target_group);
    void setNature(NatureOfField nature);
    void attachLocalField( MEDCouplingFieldDouble* field);
    void attachLocalField(const ParaFIELD* field, bool ownPt=false);
    void attachLocalField(const ICoCo::Field* field);
    
    virtual void prepareSourceDE()=0;
    virtual void prepareTargetDE()=0;
    virtual void recvData()=0;
    virtual void sendData()=0;
    virtual void synchronize()=0;
    virtual ~DEC();
    virtual void computeProcGroup() { }
    void renormalizeTargetField(bool isWAbs);
  protected:
    void compareFieldAndMethod() const throw(INTERP_KERNEL::Exception);
  protected:
    const ParaFIELD* _local_field;
    //! Processor group representing the union of target and source processors
    ProcessorGroup* _union_group;
    ProcessorGroup* _source_group;
    ProcessorGroup* _target_group;
    
    const CommInterface* _comm_interface;
    bool _owns_field;
  private:
    ICoCo::Field* _icoco_field;
  };
}

#endif
