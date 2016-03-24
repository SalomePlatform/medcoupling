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

#ifndef __PARAFIELD_HXX__
#define __PARAFIELD_HXX__

#include "MEDCouplingRefCountObject.hxx"
#include "ComponentTopology.hxx"

namespace MEDCoupling
{
  class DataArrayInt;
  class ParaMESH;
  class ProcessorGroup;
  class MEDCouplingFieldDouble;
  class ComponentTopology;
  class Topology;

  class ParaFIELD
  {
  public:
    ParaFIELD(TypeOfField type, TypeOfTimeDiscretization td, ParaMESH* mesh, const ComponentTopology& component_topology); 
    ParaFIELD(MEDCouplingFieldDouble* field, ParaMESH *sup, const ProcessorGroup& group);
    virtual ~ParaFIELD();

    void synchronizeTarget( MEDCoupling::ParaFIELD* source_field);
    void synchronizeSource( MEDCoupling::ParaFIELD* target_field);
    MEDCouplingFieldDouble* getField() const { return _field; }
    void setOwnSupport(bool v) const { _own_support=v; }
    DataArrayInt* returnCumulativeGlobalNumbering() const;
    DataArrayInt* returnGlobalNumbering() const;
    Topology* getTopology() const { return _topology; }
    ParaMESH* getSupport() const  { return _support; }
    int nbComponents() const;
    double getVolumeIntegral(int icomp, bool isWAbs) const;
    double getL2Norm()const { return -1; }

  private:
    MEDCouplingFieldDouble* _field;
    MEDCoupling::ComponentTopology _component_topology;
    Topology* _topology; 
    mutable bool _own_support;
    ParaMESH* _support;
  };

}

#endif
