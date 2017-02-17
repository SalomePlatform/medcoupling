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
#include "MEDPARTITIONER_UserGraph.hxx"

#include "MEDCouplingSkyLineArray.hxx"

#include <iostream>
#include <vector>

using namespace MEDPARTITIONER;

/*! constructor that allows to specify the desired partition
 * \param partition as table giving the domain number for each cell 
 *        (domain numbers range from 0 to ndomain-1
 * \param n number of cells in the mesh
 */
UserGraph::UserGraph(MEDCoupling::MEDCouplingSkyLineArray *array, const int *partition, int n):Graph(array,0)
{

  std::vector<int> index(n+1),value(n);

  index[0]=0;
  for (int i=0; i<n; i++)
    {
      index[i+1]=index[i]+1;
      value[i]=partition[i];
    }

  _partition = MEDCoupling::MEDCouplingSkyLineArray::New(index,value);

}

UserGraph::~UserGraph()
{
}

void UserGraph::partGraph(int ndomain, const std::string& options, ParaDomainSelector* sel)
{
  std::cerr << "MEDPARTITIONER::UserGraph::partGraph() should not have to be used" << std::endl;
}

