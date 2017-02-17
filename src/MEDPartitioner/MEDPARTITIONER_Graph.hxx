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

#ifndef __MEDPARTITIONER_GRAPH_HXX__
#define __MEDPARTITIONER_GRAPH_HXX__

#include "MEDPARTITIONER.hxx"
#include "MCAuto.hxx"

#include <string>

namespace MEDCoupling
{
  class MEDCouplingSkyLineArray;
}

using namespace MEDCoupling;

namespace MEDPARTITIONER 
{
  class ParaDomainSelector;
  class MEDPARTITIONER_EXPORT Graph
  {
  public:
    typedef enum {METIS,SCOTCH} splitter_type;

    Graph();
    //creates a graph from a SKYLINEARRAY- WARNING!! Graph takes ownership of the array.
    Graph(MEDCouplingSkyLineArray* graph, int* edgeweight=0);
    virtual ~Graph();

    void setEdgesWeights(int *edgeweight) { _edge_weight=edgeweight; }
    void setVerticesWeights(int *cellweight) { _cell_weight=cellweight; }
    
    //computes partitioning of the graph
    virtual void partGraph(int ndomain, const std::string& options_string="", ParaDomainSelector *sel=0) = 0;
    
    //returns the partitioning
    const int *getPart() const;
    
    //returns the number of graph vertices (which can correspond to the cells in the mesh!)
    int nbVertices() const;

    // returns nb of domains in _partition
    int nbDomains() const;
    
    const MEDCouplingSkyLineArray *getGraph() const { return (const MEDCouplingSkyLineArray*)_graph; }
    const MEDCouplingSkyLineArray *getPartition() const { return (const MEDCouplingSkyLineArray*)_partition; }

  protected:
    MCAuto<MEDCouplingSkyLineArray> _graph;
    MCAuto<MEDCouplingSkyLineArray> _partition;
    int* _edge_weight;  
    int* _cell_weight;
  };
}
#endif
