// Copyright (C) 2007-2024  CEA, EDF
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

namespace MEDPARTITIONER
{
  Graph::Graph():
    _graph(0),_partition(0),
    _edge_weight(0),_cell_weight(0)
  {
  }

  Graph::Graph(MEDCoupling::MEDCouplingSkyLineArray *array, int *edgeweight):
    _graph(array),_partition(0),
    _edge_weight(edgeweight),_cell_weight(0)
  {
  }

  Graph::~Graph()
  {
  }

  int Graph::nbDomains() const
  {
    std::set<mcIdType> domains;
    if ( _partition.isNotNull() )
      if ( MEDCoupling::DataArrayIdType* array = _partition->getValuesArray() )
      {
        for ( const mcIdType * dom = array->begin(); dom != array->end(); ++dom )
          domains.insert( *dom );
      }
    return (int)domains.size();
  }

  const mcIdType *Graph::getPart() const
  {
    return _partition->getValues();
  }

  mcIdType Graph::nbVertices() const
  {
    return _graph->getNumberOf();
  }

};
