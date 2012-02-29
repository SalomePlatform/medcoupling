//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#ifndef MEDPARTITIONER_GRAPH_HXX_
#define MEDPARTITIONER_GRAPH_HXX_

#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_SkyLineArray.hxx"
#include <string>
namespace MEDPARTITIONER {
        
class ParaDomainSelector;
class MEDPARTITIONER_EXPORT Graph
{
  public:
    typedef enum {METIS,SCOTCH} splitter_type;

    Graph(){};

    //creates a graph from a SKYLINEARRAY
    Graph(MEDPARTITIONER::MEDSKYLINEARRAY* graph, int* edgeweight=0);

    virtual ~Graph();

    void setEdgesWeights(int* edgeweight){_edgeweight=edgeweight;}
    void setVerticesWeights(int* cellweight){_cellweight=cellweight;}

    //computes partitioning of the graph
    virtual void partGraph(int ndomain, const std::string&, ParaDomainSelector* sel=0)=0;

    //! returns the partitioning
    const int* getPart() const {return _partition->getValue();}

    //! returns the number of graph vertices (which can correspond to the cells in the mesh!)
    int nbVertices() const {return _graph->getNumberOf();}

    const MEDPARTITIONER::MEDSKYLINEARRAY* getGraph() const {return _graph;}

  protected:

    MEDPARTITIONER::MEDSKYLINEARRAY* _graph;

    MEDPARTITIONER::MEDSKYLINEARRAY* _partition;

    int* _edgeweight;  

    int* _cellweight;
  };

}
#endif /*GRAPH_HXX_*/
