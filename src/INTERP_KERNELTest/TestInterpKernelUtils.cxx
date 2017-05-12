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

#include "TestInterpKernelUtils.hxx"
#include "InterpKernelException.hxx"

#include <cstdlib>
#ifdef WIN32
#include<direct.h>
#define getcwd _getcwd
#else
#include <unistd.h>
#endif
#include <sstream>
#include <fstream>

namespace INTERP_TEST
{
  std::string getResourceFile( const std::string& filename, int levelUp)
  {
    std::string resourceFile = "";
    if ( getenv("MEDCOUPLING_ROOT_DIR") ) {
      // use MEDCOUPLING_ROOT_DIR env.var
      resourceFile = getenv("MEDCOUPLING_ROOT_DIR");
      resourceFile += "/share/resources/med/";
      resourceFile += filename;
      std::ifstream my_file(resourceFile.c_str());
      if (my_file.good())
        return resourceFile;
    }
    // else
    char * tmp_c = getcwd(NULL, 0);
    resourceFile = tmp_c;
    free(tmp_c);
    resourceFile += "/";
    for(int i=0; i<levelUp; i++)
      resourceFile += "../";
    resourceFile += "resources/";
    resourceFile += filename;
    std::ifstream my_file(resourceFile.c_str());
    if (!my_file.good())
      {
        std::stringstream ss;
        ss << "INTERP_TEST::getResourceFile(): could not open resource test file: " << filename << "\n";
        throw INTERP_KERNEL::Exception(ss.str().c_str());
      }

    return resourceFile;
  }

} // namespace INTERP_TEST
