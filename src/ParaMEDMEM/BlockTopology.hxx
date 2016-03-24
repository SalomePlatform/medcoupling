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

#ifndef __BLOCKTOPOLOGY_HXX__
#define __BLOCKTOPOLOGY_HXX__

#include "Topology.hxx"
#include "ProcessorGroup.hxx"

#include <vector>

namespace MEDCoupling
{
  class ComponentTopology;
  class MEDCouplingCMesh;

  typedef enum{Block,Cycle} CYCLE_TYPE; 

  /*!
   * \anchor BlockTopology-det
   *
   * A BlockTopology typically represents the split of a *structured* mesh among the processors of
   * a common ProcessorGroup. Each processor gets a contiguous part of the cells in the mesh (a block).
   *
   * A BlockTopology can also be used to split a structured domain among the various components of a field.
   *
   * \sa ExplicitTopology
   */
  class BlockTopology : public Topology
  {
  public:
    BlockTopology();
    BlockTopology(const ProcessorGroup& group, MEDCouplingCMesh *grid); 
    BlockTopology(const BlockTopology& geom_topo, const ComponentTopology& comp_topo);
    BlockTopology(const ProcessorGroup& group, int nb_elem);
    virtual ~BlockTopology();
    //!Retrieves the number of elements for a given topology
    int getNbElements()const { return _nb_elems; }
    int getNbLocalElements() const;
    const ProcessorGroup* getProcGroup()const { return _proc_group; }
    std::pair<int,int> globalToLocal (const int) const ;
    int localToGlobal (const std::pair<int,int>) const;
    std::vector<std::pair<int,int> > getLocalArrayMinMax() const ;
    int getDimension() const { return _dimension; }
    void serialize(int* & serializer, int& size) const ;
    void unserialize(const int* serializer, const CommInterface& comm_interface);
  private:
    //dimension : 2 or 3
    int _dimension;
    //proc array
    std::vector<int> _nb_procs_per_dim;
    //stores the offsets vector  
    std::vector<std::vector<int> > _local_array_indices;
    //stores the cycle type (block or cyclic)
    std::vector<CYCLE_TYPE> _cycle_type;
    //Processor group
    const ProcessorGroup* _proc_group;
    //nb of elements
    int _nb_elems;
    bool _owns_processor_group;
  };
}

#endif
