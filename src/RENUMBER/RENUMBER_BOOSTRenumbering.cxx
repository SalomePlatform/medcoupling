//  Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include "RENUMBER_BOOSTRenumbering.hxx"

void BOOSTRenumbering::renumber(const int* graph,const int* index_graph,int nb_cell,std::vector<int>& iperm,std::vector<int>& perm)
{
  iperm.resize(nb_cell,0);
  perm.resize(nb_cell,0);

  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
     boost::property<boost::vertex_color_t, boost::default_color_type,
       boost::property<boost::vertex_degree_t,int> > > Graph;
  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<Graph>::vertices_size_type size_type;
  Graph G(nb_cell);
  for (int i=0;i<nb_cell;++i)
    for (int j=index_graph[i]-1;j<index_graph[i+1]-1;++j)
      add_edge(i,graph[j]-1,G);
  boost::property_map<Graph, boost::vertex_index_t>::type
    index_map = boost::get(boost::vertex_index, G);
  boost::cuthill_mckee_ordering(G, iperm.rbegin(), boost::get(boost::vertex_color, G),
                           boost::make_degree_map(G));
  for (size_type c = 0; c != iperm.size(); ++c)
    perm[index_map[iperm[c]]] = c;
  for(int i=0;i<nb_cell;++i)
    {
      perm[i]+=1;
      iperm[i]+=1;
    }
}
