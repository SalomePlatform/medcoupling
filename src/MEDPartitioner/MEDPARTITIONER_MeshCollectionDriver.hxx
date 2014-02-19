// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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

#ifndef __MEDPARTITIONER_MESHCOLLECTIONDRIVER_HXX__
#define __MEDPARTITIONER_MESHCOLLECTIONDRIVER_HXX__

#include "MEDPARTITIONER.hxx"

#include <vector>
#include <string>

namespace MEDPARTITIONER
{
  class MeshCollection;
  class ParaDomainSelector;

  class MEDPARTITIONER_EXPORT MeshCollectionDriver
  {
  public:
    MeshCollectionDriver(MeshCollection*);
    virtual ~MeshCollectionDriver() { }
    virtual int read(const char*, ParaDomainSelector* sel=0) = 0;
    int readSeq(const char*,const char*);
    virtual void write(const char*, ParaDomainSelector* sel=0) const = 0;
  protected:
    void readSubdomain(std::vector<int*>& cellglobal,
                       std::vector<int*>& faceglobal,
                       std::vector<int*>& nodeglobal, int idomain);
    void readSubdomain(int idomain);
    void writeMedFile(int idomain, const std::string& distfilename) const;
  protected:
    MeshCollection* _collection;
  };
}
#endif
