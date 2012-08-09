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

#ifndef __MEDPARTITIONER_GRAPH_HXX__
#define __MEDPARTITIONER_GRAPH_HXX__

#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_SkyLineArray.hxx"

#include <string>

namespace MEDPARTITIONER 
{
  class ParaDomainSelector;
  class MEDPARTITIONER_EXPORT Graph
  {
  public:
    typedef enum {METIS,SCOTCH} splitter_type;

    Graph() { }
    //creates a graph from a SKYLINEARRAY
    Graph(MEDPARTITIONER::SkyLineArray* graph, int* edgeweight=0);
    virtual ~Graph();

    void setEdgesWeights(int *edgeweight) { _edge_weight=edgeweight; }
    void setVerticesWeights(int *cellweight) { _cell_weight=cellweight; }
    
    //computes partitioning of the graph
    virtual void partGraph(int ndomain, const std::string&, ParaDomainSelector *sel=0) = 0;
    
    //returns the partitioning
    const int *getPart() const { return _partition->getValue(); }
    
    //returns the number of graph vertices (which can correspond to the cells in the mesh!)
    int nbVertices() const { return _graph->getNumberOf(); }
    
    const SkyLineArray *getGraph() const { return _graph; }

  protected:
    SkyLineArray* _graph;
    SkyLineArray* _partition;
    int* _edge_weight;  
    int* _cell_weight;
  };
}
#endif
