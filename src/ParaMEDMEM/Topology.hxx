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

#ifndef __TOPOLOGY_HXX__
#define __TOPOLOGY_HXX__

#include <utility>

namespace ParaMEDMEM
{
  class ProcessorGroup;

  class Topology
  {
  public:
    Topology() { }
    virtual ~Topology() { }
    virtual int getNbElements() const = 0;
    virtual int getNbLocalElements() const  = 0;
    virtual const ProcessorGroup* getProcGroup()const  = 0;
  };
}

#endif
