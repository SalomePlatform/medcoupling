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

#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "ParaMESH.hxx"
#include "Topology.hxx"
#include "ExplicitTopology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"

#include <vector>
#include <algorithm>

using namespace std;
namespace MEDCoupling
{

ExplicitTopology::ExplicitTopology():
   _proc_group(NULL), _nb_elems(0), _nb_components(0),
   _loc2glob(NULL), _glob2loc()
  {}

ExplicitTopology::ExplicitTopology(const ParaMESH& paramesh ):
_proc_group(paramesh.getBlockTopology()->getProcGroup()),
_nb_components(1)
{
  _nb_elems=paramesh.getCellMesh()->getNumberOfCells();
  const int* global=paramesh.getGlobalNumberingCell();
  _loc2glob=new int[_nb_elems]; 
  
    for (int i=0; i<_nb_elems; i++)
    {
      _loc2glob[i]=global[i];
      _glob2loc[global[i]]=i;
    }
}

ExplicitTopology::ExplicitTopology(const ExplicitTopology& topo, int nb_components)
{
  _proc_group = topo._proc_group;
  _nb_elems = topo._nb_elems;
  _nb_components = nb_components;
  _loc2glob=new int[_nb_elems];
  for (int i=0; i<_nb_elems; i++)
    {
      _loc2glob[i]=topo._loc2glob[i];
    }
  _glob2loc=topo._glob2loc;
}


ExplicitTopology::~ExplicitTopology()
{
  if (_loc2glob != 0) delete[] _loc2glob;
}


/*! Serializes the data contained in the Explicit Topology
 * for communication purposes*/
void ExplicitTopology::serialize(int* & serializer, int& size) const 
{
  vector <int> buffer;
  
  buffer.push_back(_nb_elems);
  for (int i=0; i<_nb_elems; i++)
  {
    buffer.push_back(_loc2glob[i]);
  }
    
  serializer=new int[buffer.size()];
  size=  buffer.size();
  copy(buffer.begin(), buffer.end(), serializer);
  
}
/*! Unserializes the data contained in the Explicit Topology
 * after communication. Uses the same structure as the one used for serialize()
 * 
 * */
void ExplicitTopology::unserialize(const int* serializer,const CommInterface& comm_interface)
{
  const int* ptr_serializer=serializer;
  cout << "unserialize..."<<endl;
  _nb_elems=*ptr_serializer++;
  cout << "nbelems "<<_nb_elems<<endl;
  _loc2glob=new int[_nb_elems];
  for (int i=0; i<_nb_elems; i++)
  {
    _loc2glob[i]=*ptr_serializer;
    _glob2loc[*ptr_serializer]=i;
    ptr_serializer++;
    
  }

}

}
