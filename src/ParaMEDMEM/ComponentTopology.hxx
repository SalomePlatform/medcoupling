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

#ifndef __COMPONENTTOPOLOGY_HXX__
#define __COMPONENTTOPOLOGY_HXX__

#include "Topology.hxx"

#include <vector>

namespace MEDCoupling
{
  class ProcessorGroup;

  /*!
   * \anchor ComponentTopology-det
   *
   * The ComponentTopology can be used when building a ParaFIELD. It allows the splitting of the components
   * of the field among different processors within a single processor group.
   *
   * \sa ParaFIELD::ParaFIELD(TypeOfField , TypeOfTimeDiscretization , ParaMESH* , const ComponentTopology& )
   */
  class ComponentTopology
  {
  public:
    ComponentTopology(int nb_comp, ProcessorGroup* group);
    ComponentTopology(int nb_comp, int nb_blocks);
    ComponentTopology(int nb_comp);
    ComponentTopology();
    virtual ~ComponentTopology();
    //!returns the number of MED components in the topology
    int nbComponents() const { return _component_array.back(); }
    //!returns the number of MED components on local processor
    int nbLocalComponents() const ;
    //!returns the number of the first MED component on local processor
    int firstLocalComponent() const ;
    //!returns the number of blocks in the topology
    int nbBlocks()const {return _component_array.size()-1;}
    //!returns the block structure
    const std::vector<int>* getBlockIndices() const { return &_component_array; }
    const ProcessorGroup* getProcGroup()const { return _proc_group; } 
  private:
    std::vector<int> _component_array;
    ProcessorGroup* _proc_group;
  };
}

#endif /*COMPONENTTOPOLOGY_HXX_*/
