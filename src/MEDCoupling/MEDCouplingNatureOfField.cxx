// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#include "MEDCouplingNatureOfField.hxx"

#include <algorithm>

namespace ParaMEDMEM
{
  const char *MEDCouplingNatureOfField::REPR_OF_NATUREOFFIELD[NB_OF_POSSIBILITIES]=
    { "NoNature",
      "ConservativeVolumic",
      "Integral",
      "IntegralGlobConstraint",
      "RevIntegral"};
  
  const int MEDCouplingNatureOfField::POS_OF_NATUREOFFIELD[NB_OF_POSSIBILITIES]={17,26,32,35,37};

  const char *MEDCouplingNatureOfField::getRepr(NatureOfField nat) throw(INTERP_KERNEL::Exception)
  {
    const int *pos=std::find(POS_OF_NATUREOFFIELD,POS_OF_NATUREOFFIELD+NB_OF_POSSIBILITIES,(int)nat);
    if(pos==POS_OF_NATUREOFFIELD+NB_OF_POSSIBILITIES)
      throw INTERP_KERNEL::Exception("MEDCouplingNatureOfField::getRepr : Unrecognized nature of field !");
    int pos2=std::distance(POS_OF_NATUREOFFIELD,pos);
    return REPR_OF_NATUREOFFIELD[pos2];
  }
}
 
