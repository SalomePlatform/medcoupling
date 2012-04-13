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

#ifndef __MEDPARTITIONER_JOINTFINDER_HXX__
#define __MEDPARTITIONER_JOINTFINDER_HXX__

#include <map>
#include <vector>

namespace MEDPARTITIONER
{
  class Topology;
  class MeshCollection;
  class ParaDomainSelector;
  
  class JointFinder
  {
  public:
    JointFinder(const MeshCollection& mc);
    ~JointFinder();
    void findCommonDistantNodes();
    void print();
    std::vector<std::vector<std::multimap<int,int> > >& getDistantNodeCell();
    std::vector<std::vector<std::vector<std::pair<int,int> > > >& getNodeNode();
  private:
    const MeshCollection& _mesh_collection;
    const ParaDomainSelector *_domain_selector;
    const Topology *_topology;
    std::vector<std::vector<std::multimap<int,int> > > _distant_node_cell;
    std::vector<std::vector<std::vector<std::pair<int,int> > > > _node_node;
   
  };
}
#endif
