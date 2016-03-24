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

#include "ComponentTopology.hxx"
#include "ProcessorGroup.hxx"
#include "InterpolationUtils.hxx"

namespace MEDCoupling
{
  /* Generic constructor for \a nb_comp components equally parted
   * in \a nb_blocks blocks
   */
  ComponentTopology::ComponentTopology(int nb_comp, ProcessorGroup* group):_proc_group(group)
  {
    int nb_blocks=group->size();
  
    if (nb_blocks>nb_comp)
      throw INTERP_KERNEL::Exception("ComponentTopology Number of components must be larger than number of blocks");

    _component_array.resize(nb_blocks+1);
    _component_array[0]=0;
    for (int i=1; i<=nb_blocks; i++)
      {
        _component_array[i]=_component_array[i-1]+nb_comp/nb_blocks;
        if (i<=nb_comp%nb_blocks)
          _component_array[i]++;
      }
  }
  
  /* Generic constructor for \a nb_comp components equally parted
   * in \a nb_blocks blocks
   */
  ComponentTopology::ComponentTopology(int nb_comp, int nb_blocks):_proc_group(0)
  {
    if (nb_blocks>nb_comp)
      throw INTERP_KERNEL::Exception("ComponentTopology Number of components must be larger than number of blocks");
    
    _component_array.resize(nb_blocks+1);
    _component_array[0]=0;
    for (int i=1; i<=nb_blocks; i++)
      {
        _component_array[i]=_component_array[i-1]+nb_comp/nb_blocks;
        if (i<=nb_comp%nb_blocks)
          _component_array[i]++;
      }
  
  }
  
  //!Constructor for one block of \a nb_comp components
  ComponentTopology::ComponentTopology(int nb_comp):_proc_group(0)
  {
    
    _component_array.resize(2);
    _component_array[0]=0;
    _component_array[1]=nb_comp;
  
  }

  //! Constructor for one component
  ComponentTopology::ComponentTopology():_proc_group(0)
  {
    _component_array.resize(2);
    _component_array[0]=0;
    _component_array[1]=1;
  
  }
  
  ComponentTopology::~ComponentTopology()
  {
  }

  int ComponentTopology::nbLocalComponents() const
  {
    if (_proc_group==0)
      return nbComponents();
  
    int nbcomp;
    int myrank = _proc_group->myRank();
    if (myrank!=-1)
      nbcomp = _component_array[myrank+1]-_component_array[myrank];
    else 
      nbcomp=0;
    return nbcomp;
  }

  int ComponentTopology::firstLocalComponent() const
  {
    if (_proc_group==0)
      return 0;
  
    int icomp;
    int myrank = _proc_group->myRank();
    if (myrank!=-1)
      icomp = _component_array[myrank];
    else 
      icomp=-1;
    return icomp;
  }
}
