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

#include "ParaGRID.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingCMesh.hxx"
#include "InterpKernelUtilities.hxx"

#include <iostream>

using namespace std;

namespace ParaMEDMEM
{
  
  ParaGRID::ParaGRID(MEDCouplingCMesh* global_grid, Topology* topology) throw(INTERP_KERNEL::Exception)
  {
  
    _block_topology = dynamic_cast<BlockTopology*>(topology);
    if(_block_topology==0)
      throw INTERP_KERNEL::Exception(LOCALIZED("ParaGRID::ParaGRID topology must be block topology"));
    
    if (!_block_topology->getProcGroup()->containsMyRank())
      return;
    
    int dimension=_block_topology->getDimension() ;
    if (dimension != global_grid->getSpaceDimension())
      throw INTERP_KERNEL::Exception(LOCALIZED("ParaGrid::ParaGrid incompatible topology"));
    _grid=global_grid;
    _grid->incrRef();
    /*vector<vector<double> > xyz_array(dimension);
      vector<pair<int,int> > local_indices = _block_topology->getLocalArrayMinMax();
      vector <string> coordinates_names;
      vector <string> coordinates_units;
      for (int idim=0; idim<dimension ; idim++)
      {
      DataArrayDouble *array=global_grid->getCoordsAt(idim);
      double *arrayC=array->getPointer();
      cout << " Indices "<< local_indices[idim].first <<" "<<local_indices[idim].second<<endl;
      for (int i=(local_indices)[idim].first; i<(local_indices)[idim].second; i++)
      xyz_array[idim].push_back(arrayC[i]);
      coordinates_names.push_back(array->getName());
      coordinates_units.push_back(array->getInfoOnComponentAt(0));
      }
      _grid=MEDCouplingCMesh::New();
      _grid->set(xyz_array, coordinates_names,coordinates_units);
      _grid->setName(global_grid->getName());
      _grid->setDescription(global_grid->getDescription());*/
  }

  ParaGRID::~ParaGRID()
  {
    if(_grid)
      _grid->decrRef();
  }
}
