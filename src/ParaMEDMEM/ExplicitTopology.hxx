// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

#ifndef __EXPLICITTOPOLOGY_HXX__
#define __EXPLICITTOPOLOGY_HXX__

#include "ProcessorGroup.hxx"
#include "InterpKernelHashMap.hxx"

#include <vector>
#include <utility>
#include <iostream>

namespace MEDCoupling
{
  class ParaMESH;
  class Topology;
  class ComponentTopology;

  /*!
   * \anchor ExplicitTopology-det
   *
   * An ExplicitTopology typically represents the split of a mesh among the processors of
   * a common ProcessorGroup. Each processor gets a user-defined part of the cells in the mesh.
   * \sa BlockTopology
   */
  class ExplicitTopology : public Topology
  {
  public:
    ExplicitTopology();
    ExplicitTopology( const ExplicitTopology& topo, int nbcomponents);
    ExplicitTopology(const ParaMESH &mesh);
    virtual ~ExplicitTopology();
    
    inline mcIdType getNbElements()const;
    inline mcIdType getNbLocalElements() const;
    const ProcessorGroup* getProcGroup()const { return _proc_group; }
    mcIdType localToGlobal (const std::pair<int,mcIdType> local) const { return localToGlobal(local.second); }
    inline mcIdType localToGlobal(mcIdType) const;
    inline mcIdType globalToLocal(mcIdType) const;
    void serialize(mcIdType* & serializer, mcIdType& size) const ;
    void unserialize(const mcIdType* serializer, const CommInterface& comm_interface);
    int getNbComponents() const { return _nb_components; }
  private:
    //Processor group
    const ProcessorGroup* _proc_group;
    //nb of elements
    mcIdType _nb_elems;
    //nb of components
    int _nb_components;
    //mapping local to global
    mcIdType* _loc2glob;
    //mapping global to local
    INTERP_KERNEL::HashMap<mcIdType,mcIdType> _glob2loc;
  };

  //!converts a pair <subdomainid,local> to a global number 
  inline mcIdType ExplicitTopology::globalToLocal(const mcIdType global) const
  {
    return (_glob2loc.find(global))->second;
  }

  //!converts local number to a global number
  mcIdType ExplicitTopology::localToGlobal(mcIdType local) const
  {
    return _loc2glob[local];
  }
  
  //!Retrieves the number of elements for a given topology
  inline mcIdType ExplicitTopology::getNbElements() const
  {
    return _nb_elems;
  }

  //Retrieves the local number of elements 
  inline mcIdType ExplicitTopology::getNbLocalElements()const 
  {
    return ToIdType(_glob2loc.size());
  }
}


#endif
