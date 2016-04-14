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

#ifndef __MEDPARTITIONER_MESHCOLLECTIONMEDASCIIDRIVER_HXX__
#define __MEDPARTITIONER_MESHCOLLECTIONMEDASCIIDRIVER_HXX__

#include "MEDPARTITIONER.hxx"
#include "MEDPARTITIONER_MeshCollectionDriver.hxx"

namespace MEDPARTITIONER
{
  class MeshCollection;
  class MEDPARTITIONER_EXPORT MeshCollectionMedAsciiDriver : public MeshCollectionDriver
  {
  public:
    MeshCollectionMedAsciiDriver(MeshCollection*);
    virtual ~MeshCollectionMedAsciiDriver() { }
    int read(const char*, ParaDomainSelector* sel=0);
    int read(MEDCoupling::MEDFileData*);
    void write(const char*, ParaDomainSelector* sel=0) const;
  private:
    std::string _master_filename;
  };
}
#endif
