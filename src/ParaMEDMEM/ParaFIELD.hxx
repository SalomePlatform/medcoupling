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
#ifndef __PARAFIELD_HXX__
#define __PARAFIELD_HXX__

#include "ComponentTopology.hxx"
#include "ParaMESH.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingUMesh.hxx"

namespace ParaMEDMEM
{

  class ParaSUPPORT;
  class ProcessorGroup;

  class ParaFIELD
  {
  public:

    ParaFIELD(TypeOfField type, ParaMESH* mesh, const ComponentTopology& component_topology); 


    ParaFIELD(MEDCouplingFieldDouble* field, const ProcessorGroup& group);
  
    virtual ~ParaFIELD();
    void synchronizeTarget( ParaMEDMEM::ParaFIELD* source_field);
    void synchronizeSource( ParaMEDMEM::ParaFIELD* target_field);
    MEDCouplingFieldDouble* getField() const { return _field; }
    Topology* getTopology() const { return _topology; }
    ParaMESH* getSupport() const  { return _support; }
    int nbComponents() const { return _component_topology.nbComponents(); }
    double getVolumeIntegral(int icomp) const;
    double getL2Norm()const { return -1; }
  private:
    MEDCouplingFieldDouble* _field;
    const  ParaMEDMEM::ComponentTopology& _component_topology;
    Topology* _topology; 

    ParaMESH* _support;
    bool _has_field_ownership;
    bool _has_support_ownership;
  };

}

#endif
