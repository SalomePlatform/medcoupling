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

#ifndef __PARAGRID_HXX__
#define __PARAGRID_HXX__

#include "InterpolationUtils.hxx"

#include <vector>

namespace MEDCoupling
{
  class Topology;
  class BlockTopology;
  class MEDCouplingCMesh;

  /*!
   * This class
   * Equivalent of a ParaMESH for a structured mesh
   */
  class ParaGRID
  {
  public:
    ParaGRID(MEDCouplingCMesh* global_grid, Topology* topology) throw(INTERP_KERNEL::Exception);
    BlockTopology * getBlockTopology() const { return _block_topology; }
    virtual ~ParaGRID();
    MEDCouplingCMesh* getGrid() const { return _grid; }
  private:
    MEDCouplingCMesh* _grid;
    // structured grid topology
    MEDCoupling::BlockTopology* _block_topology;
    // stores the x,y,z axes on the global grid
    std::vector<std::vector<double> > _global_axis;
    //id of the local grid
    int _my_domain_id;
  };
}

#endif
