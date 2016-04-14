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

#include "RENUMBER_BOOSTRenumbering.hxx"

#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

void BOOSTRenumbering::renumber(const int *graph, const int *index_graph, int nbCell, MEDCoupling::DataArrayInt *&iperm, MEDCoupling::DataArrayInt *&perm)
{
  MEDCoupling::MCAuto<MEDCoupling::DataArrayInt> out0(MEDCoupling::DataArrayInt::New()),out1(MEDCoupling::DataArrayInt::New());
  out0->alloc(nbCell,1); out1->alloc(nbCell,1);
  out0->fillWithZero(); out1->fillWithZero();
  //
  typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
     boost::property<boost::vertex_color_t, boost::default_color_type,
       boost::property<boost::vertex_degree_t,int> > > Graph;
  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<Graph>::vertices_size_type size_type;
  Graph G(nbCell);
  for (int i=0;i<nbCell;++i)
    for (int j=index_graph[i];j<index_graph[i+1];++j)
      add_edge(i,graph[j],G);
  boost::property_map<Graph, boost::vertex_index_t>::type
    index_map = boost::get(boost::vertex_index, G);
  boost::cuthill_mckee_ordering(G, out0->getPointer(), boost::get(boost::vertex_color, G),
                           boost::make_degree_map(G));
  int *out0Ptr(out0->getPointer()),*out1Ptr(out1->getPointer());
  for(int c=0;c!=nbCell;++c)
    out1Ptr[index_map[out0Ptr[nbCell-c-1]]]=c;
  out0->reverse();
  iperm=out0.retn(); perm=out1.retn();
}
