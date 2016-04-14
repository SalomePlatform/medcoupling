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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingNatureOfField.hxx"

#include <algorithm>
#include <sstream>

namespace MEDCoupling
{
  const char *MEDCouplingNatureOfField::REPR_OF_NATUREOFFIELD[NB_OF_POSSIBILITIES]=
  { "NoNature",
    "IntensiveMaximum",
    "ExtensiveMaximum",
    "ExtensiveConservation",
    "IntensiveConservation"};

  const int MEDCouplingNatureOfField::POS_OF_NATUREOFFIELD[NB_OF_POSSIBILITIES]={17,26,32,35,37};

  const char *MEDCouplingNatureOfField::GetRepr(NatureOfField nat)
  {
    const int *pos=std::find(POS_OF_NATUREOFFIELD,POS_OF_NATUREOFFIELD+NB_OF_POSSIBILITIES,(int)nat);
    if(pos==POS_OF_NATUREOFFIELD+NB_OF_POSSIBILITIES)
      {
        std::ostringstream oss; oss << "MEDCouplingNatureOfField::getRepr : Unrecognized nature of field ! ";
        oss << GetAllPossibilitiesStr() << " !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    std::size_t pos2=std::distance(POS_OF_NATUREOFFIELD,pos);
    return REPR_OF_NATUREOFFIELD[pos2];
  }

  std::string MEDCouplingNatureOfField::GetReprNoThrow(NatureOfField nat)
  {
    const int *pos=std::find(POS_OF_NATUREOFFIELD,POS_OF_NATUREOFFIELD+NB_OF_POSSIBILITIES,(int)nat);
    if(pos==POS_OF_NATUREOFFIELD+NB_OF_POSSIBILITIES)
      return std::string("Unrecognized nature of field !");
    std::size_t pos2=std::distance(POS_OF_NATUREOFFIELD,pos);
    return std::string(REPR_OF_NATUREOFFIELD[pos2]);
  }

  std::string MEDCouplingNatureOfField::GetAllPossibilitiesStr()
  {
    std::ostringstream oss; oss << "Possibilities are : ";
    for(int i=0;i<NB_OF_POSSIBILITIES;i++)
      {
        oss << REPR_OF_NATUREOFFIELD[i] << "(value=" << POS_OF_NATUREOFFIELD[i] << ")";
        if(i!=NB_OF_POSSIBILITIES-1)
          oss << ", ";
      }
    return oss.str();
  }
}
