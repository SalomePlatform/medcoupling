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

#include "ProcessorGroup.hxx"
#include "InterpolationUtils.hxx"

namespace ParaMEDMEM
{
  ProcessorGroup::ProcessorGroup (const CommInterface& interface, int start, int end):_comm_interface(interface)
  {
    if (start>end)
      throw INTERP_KERNEL::Exception("wrong call to Processor group constructor");
    for (int i=start; i<=end;i++)
      _proc_ids.insert(i);
  }
}
