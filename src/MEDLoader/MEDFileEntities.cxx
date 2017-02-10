// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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

#include "MEDFileEntities.hxx"

using namespace MEDCoupling;

MEDFileEntities *MEDFileEntities::BuildFrom(const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> > *entities)
{
  if(!entities)
    return new MEDFileAllStaticEntites;
  else
    return new MEDFileStaticEntities(*entities);
}

MEDFileEntities *MEDFileEntities::BuildFrom(const MEDFileStructureElements& se)
{
  if(se.getNumberOf()==0)
    return new MEDFileAllStaticEntites;
  else
    return new MEDFileAllStaticEntitiesPlusDyn(&se);
}

MEDFileEntities::~MEDFileEntities()
{
}

//////////////

std::vector<int> MEDFileStaticEntities::getDynGTAvail() const
{
  return std::vector<int>();
}

bool MEDFileStaticEntities::areAllStaticTypesPresent() const
{
  return false;
}

//////////////


std::vector<int> MEDFileAllStaticEntites::getDynGTAvail() const
{
  return std::vector<int>();
}

bool MEDFileAllStaticEntites::areAllStaticTypesPresent() const
{
  return true;
}

//////////////

MEDFileAllStaticEntitiesPlusDyn::MEDFileAllStaticEntitiesPlusDyn(const MEDFileStructureElements *se):_se(se)
{
  if(se)
    se->incrRef();
}

std::vector<int> MEDFileAllStaticEntitiesPlusDyn::getDynGTAvail() const
{
  return _se->getDynGTAvail();
}

bool MEDFileAllStaticEntitiesPlusDyn::areAllStaticTypesPresent() const
{
  return true;
}

const MEDFileStructureElement *MEDFileAllStaticEntitiesPlusDyn::getWithGT(int idGT) const
{
  return _se->getWithGT(idGT);
}

const MEDFileUMesh *MEDFileAllStaticEntitiesPlusDyn::getSupMeshWithName(const std::string& name) const
{
  return _se->getSupMeshWithName(name);
}
