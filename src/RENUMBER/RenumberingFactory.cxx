// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

#include "RenumberingFactory.hxx"
#include "RENUMBER_Renumbering.hxx"
#ifdef MED_ENABLE_METIS
#include "RENUMBER_METISRenumbering.hxx"
#endif
#ifdef ENABLE_BOOST
#include "RENUMBER_BOOSTRenumbering.hxx"
#endif

#include <iostream>
#include <algorithm>

namespace MED_RENUMBER
{
  bool CompareRenumMeth(const std::string& s1, const char *s2)
  {
    std::string ss1(s1),ss2(s2);
    std::transform(ss1.begin(), ss1.end(), ss1.begin(), ::tolower);
    std::transform(ss2.begin(), ss2.end(), ss2.begin(), ::tolower);
    return ss1==ss2;
  }
  
  Renumbering* RenumberingFactory(const std::string &s)
  {
#ifdef MED_ENABLE_METIS
#ifdef ENABLE_BOOST
    if ( CompareRenumMeth(s,METIS_ALG) )
      {
        return new METISRenumbering;
      }
    else if( CompareRenumMeth(s,BOOST_ALG) )
      {
        return new BOOSTRenumbering;
      }
    else 
      {
        std::cerr << "The method has to be METIS or BOOST" << std::endl;
        return 0;
      }
#endif
#ifndef ENABLE_BOOST
    if ( CompareRenumMeth(s,METIS_ALG) )
      {
        return new METISRenumbering;
      }
    else
      {
        std::cerr << "The method has to be METIS!" << std::endl;
        return 0;
      }
#endif
#endif
#ifndef MED_ENABLE_METIS
#ifdef ENABLE_BOOST
    if ( CompareRenumMeth(s,BOOST_ALG) )
      {
        return new BOOSTRenumbering;
      }
    else
      {
        std::cerr << "The method has to be BOOST!" << std::endl;
        return 0;
      }
#endif
#ifndef ENABLE_BOOST
    std::cerr << "Error, no method compiled" << std::endl;
    return 0;
#endif
#endif
  }

  std::vector<std::string> AllRenumberMethods()
  {
    std::vector<std::string> ret;
    ret.push_back(std::string(BOOST_ALG));
    ret.push_back(std::string(METIS_ALG));
    return ret;
  }
  
  std::vector<std::string> RenumberAvailableMethods()
  {
    std::vector<std::string> ret;
#ifdef ENABLE_BOOST
    ret.push_back(std::string(BOOST_ALG));
#endif
#ifdef MED_ENABLE_METIS
    ret.push_back(std::string(METIS_ALG));
#endif
    return ret;
  }
}
