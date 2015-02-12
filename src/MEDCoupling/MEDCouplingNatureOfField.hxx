// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGNATUREOFFIELD_HXX__
#define __PARAMEDMEM_MEDCOUPLINGNATUREOFFIELD_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingNatureOfFieldEnum"

#include "InterpKernelException.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingNatureOfField
  {
  public:
    MEDCOUPLING_EXPORT static const char *GetRepr(NatureOfField nat);
    MEDCOUPLING_EXPORT static std::string GetReprNoThrow(NatureOfField nat);
    MEDCOUPLING_EXPORT static std::string GetAllPossibilitiesStr();
  private:
    static const int NB_OF_POSSIBILITIES=5;
    static const char *REPR_OF_NATUREOFFIELD[NB_OF_POSSIBILITIES];
    static const int POS_OF_NATUREOFFIELD[NB_OF_POSSIBILITIES];
  };
}

#endif
