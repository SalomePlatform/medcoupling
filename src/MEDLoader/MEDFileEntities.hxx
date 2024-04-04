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
// Author : Anthony Geay (EDF R&D)

#ifndef __MEDFILEENTITIES_HXX__
#define __MEDFILEENTITIES_HXX__

#include "MCAuto.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDLoaderDefines.hxx"

#include "MEDFileStructureElement.hxx"
#include "NormalizedGeometricTypes"

#include <string>
#include <vector>
#include <utility>

namespace MEDCoupling
{
  class MEDLOADER_EXPORT MEDFileEntities
  {
  public:
    static MEDFileEntities *BuildFrom(const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities);
    static MEDFileEntities *BuildFrom(const MEDFileStructureElements& se);
    virtual std::vector<int> getDynGTAvail() const = 0;
    virtual bool areAllStaticTypesPresent() const = 0;
    virtual bool areAllStaticPresentAndNoDyn() const = 0;
    virtual ~MEDFileEntities();
  };

  class MEDLOADER_EXPORT MEDFileStaticEntities : public MEDFileEntities
  {
  public:
    MEDFileStaticEntities(const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& entities):_entities(entities) { }
    const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& getEntries() const { return _entities; }
    std::vector<int> getDynGTAvail() const override;
    bool areAllStaticTypesPresent() const override;
    bool areAllStaticPresentAndNoDyn() const override;
  private:
    std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > _entities;
  };

  class MEDLOADER_EXPORT MEDFileAllStaticEntites : public MEDFileEntities
  {
  public:
    MEDFileAllStaticEntites() = default;
    std::vector<int> getDynGTAvail() const override;
    bool areAllStaticTypesPresent() const override;
    bool areAllStaticPresentAndNoDyn() const override;
  };

  class MEDLOADER_EXPORT MEDFileAllStaticEntitiesPlusDyn : public MEDFileEntities
  {
  public:
    MEDFileAllStaticEntitiesPlusDyn(const MEDFileStructureElements *se);
    std::vector<int> getDynGTAvail() const override;
    bool areAllStaticTypesPresent() const override;
    bool areAllStaticPresentAndNoDyn() const override;
    const MEDFileStructureElement *getWithGT(int idGT) const;
    const MEDFileUMesh *getSupMeshWithName(const std::string& name) const;
  private:
    MCConstAuto<MEDFileStructureElements> _se;
  };
}

#endif
