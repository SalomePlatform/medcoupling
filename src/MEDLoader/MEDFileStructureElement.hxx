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

#ifndef __MEDFILESTRUCTUREELEMENT_HXX__
#define __MEDFILESTRUCTUREELEMENT_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileUtilities.txx"
#include "MEDFileMesh.hxx"

#include "MEDCouplingRefCountObject.hxx"

namespace MEDCoupling
{
  class MEDFileStructureElement : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileStructureElement *New(med_idt fid, int idSE);
  public:
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    void writeLL(med_idt fid) const;
  private:
    MEDFileStructureElement(med_idt fid, int idSE);
  private:
    int _id_type;
    std::string _model_name;
    INTERP_KERNEL::NormalizedCellType _geo_type;
    int _dim;
  };
  
  class MEDFileStructureElements : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileStructureElements *New(med_idt fid);
    MEDLOADER_EXPORT static MEDFileStructureElements *New();
  public:
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    void writeLL(med_idt fid) const;
  private:
    MEDFileStructureElements(med_idt fid);
    MEDFileStructureElements();
    ~MEDFileStructureElements();
  private:
    std::vector< MCAuto<MEDFileStructureElement> > _elems;
  };
}

#endif
