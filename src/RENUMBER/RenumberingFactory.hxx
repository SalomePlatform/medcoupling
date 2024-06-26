// Copyright (C) 2007-2024  CEA, EDF
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

#ifndef __RENUMBERINGFACTORY_HXX__
#define __RENUMBERINGFACTORY_HXX__

#include "RENUMBERDefines.hxx"
#include "RENUMBER_Renumbering.hxx"

#include <string>

namespace MED_RENUMBER
{
  RENUMBER_EXPORT Renumbering* RenumberingFactory(const std::string& s);
  RENUMBER_EXPORT std::vector<std::string> RenumberAvailableMethods();
  RENUMBER_EXPORT std::vector<std::string> AllRenumberMethods();
  const char METIS_ALG[]="Metis";
  const char BOOST_ALG[]="Boost";
}

#endif /*RENUMBERINGFACTORY_HXX_*/
