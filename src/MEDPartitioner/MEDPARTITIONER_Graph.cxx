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

#include "MEDPARTITIONER_Graph.hxx"

#include "MEDCouplingSkyLineArray.hxx"

#include <set>

MEDPARTITIONER::Graph::Graph(MEDCoupling::MEDCouplingSkyLineArray *array, int *edgeweight):_graph(array),_partition(0),_edge_weight(edgeweight),_cell_weight(0)
{
}

MEDPARTITIONER::Graph::~Graph()
{
  delete _partition;
  delete _graph;
}

int MEDPARTITIONER::Graph::nbDomains() const
{
  std::set<int> domains;
  if ( _partition )
    if ( MEDCoupling::DataArrayInt* array = _partition->getValueArray() )
    {
      for ( const int * dom = array->begin(); dom != array->end(); ++dom )
        domains.insert( *dom );
    }
  return domains.size();
}

const int *MEDPARTITIONER::Graph::getPart() const
{
  return _partition->getValue();
}

int MEDPARTITIONER::Graph::nbVertices() const
{
  return _graph->getNumberOf();
}
