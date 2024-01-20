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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELD_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELD_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "MCAuto.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "InterpKernelException.hxx"

#include <string>
#include <vector>

namespace MEDCoupling
{
  class DataArrayIdType;
  class DataArrayDouble;
  class MEDCouplingMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingGaussLocalization;

  class MEDCouplingField : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT virtual void checkConsistencyLight() const;
    MEDCOUPLING_EXPORT virtual bool areCompatibleForMerge(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT virtual bool areStrictlyCompatible(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT virtual bool areStrictlyCompatibleForMulDiv(const MEDCouplingField *other) const;
    MEDCOUPLING_EXPORT virtual void copyTinyStringsFrom(const MEDCouplingField *other);
    MEDCOUPLING_EXPORT void setMesh(const MEDCoupling::MEDCouplingMesh *mesh);
    MEDCOUPLING_EXPORT const MEDCoupling::MEDCouplingMesh *getMesh() const { return _mesh; }
    MEDCOUPLING_EXPORT MEDCoupling::MEDCouplingMesh *getMesh() { return const_cast<MEDCoupling::MEDCouplingMesh *>(_mesh); }
    MEDCOUPLING_EXPORT void setName(const std::string& name) { _name=name; }
    MEDCOUPLING_EXPORT std::string getDescription() const { return _desc; }
    MEDCOUPLING_EXPORT void setDescription(const std::string& desc) { _desc=desc; }
    MEDCOUPLING_EXPORT std::string getName() const { return _name; }
    MEDCOUPLING_EXPORT TypeOfField getTypeOfField() const;
    MEDCOUPLING_EXPORT NatureOfField getNature() const;
    MEDCOUPLING_EXPORT virtual void setNature(NatureOfField nat);
    MEDCOUPLING_EXPORT DataArrayDouble *getLocalizationOfDiscr() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildMeasureField(bool isAbs) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshData(const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshDataRange(mcIdType begin, mcIdType end, mcIdType step, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const;
    MEDCOUPLING_EXPORT DataArrayIdType *computeTupleIdsToSelectFromCellIds(const mcIdType *startCellIds, const mcIdType *endCellIds) const;
    MEDCOUPLING_EXPORT const MEDCouplingFieldDiscretization *getDiscretization() const { return _type; }
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *getDiscretization() { return _type; }
    MEDCOUPLING_EXPORT void setDiscretization(MEDCouplingFieldDiscretization *newDisc);
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuplesExpected() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfMeshPlacesExpected() const;
    // Gauss point specific methods
    MEDCOUPLING_EXPORT void setGaussLocalizationOnType(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                       const std::vector<double>& gsCoo, const std::vector<double>& wg);
    MEDCOUPLING_EXPORT void setGaussLocalizationOnCells(const mcIdType *begin, const mcIdType *end, const std::vector<double>& refCoo,
                                                        const std::vector<double>& gsCoo, const std::vector<double>& wg);
    MEDCOUPLING_EXPORT void clearGaussLocalizations();
    MEDCOUPLING_EXPORT MEDCouplingGaussLocalization& getGaussLocalization(int locId);
    MEDCOUPLING_EXPORT mcIdType getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT std::set<mcIdType> getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT mcIdType getNbOfGaussLocalization() const;
    MEDCOUPLING_EXPORT mcIdType getGaussLocalizationIdOfOneCell(mcIdType cellId) const;
    MEDCOUPLING_EXPORT void getCellIdsHavingGaussLocalization(int locId, std::vector<mcIdType>& cellIds) const;
    MEDCOUPLING_EXPORT const MEDCouplingGaussLocalization& getGaussLocalization(int locId) const;
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    // for MED file RW
    MEDCOUPLING_EXPORT mcIdType getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const;
    MEDCOUPLING_EXPORT virtual void reprQuickOverview(std::ostream& stream) const = 0;
  protected:
    MEDCOUPLING_EXPORT MEDCouplingField(TypeOfField type);
    MEDCOUPLING_EXPORT MEDCouplingField(const MEDCouplingField& other, bool deepCopy=true);
    MEDCOUPLING_EXPORT MEDCouplingField(MEDCouplingFieldDiscretization *type, NatureOfField nature=NoNature);
    MEDCOUPLING_EXPORT virtual ~MEDCouplingField();
    MEDCOUPLING_EXPORT bool isEqualIfNotWhyProtected(const MEDCouplingField *other, double meshPrec, std::string& reason) const;
    MEDCOUPLING_EXPORT bool isEqualWithoutConsideringStrProtected(const MEDCouplingField *other, double meshPrec) const;
  protected:
    std::string _name;
    std::string _desc;
    NatureOfField _nature;
    const MEDCouplingMesh *_mesh;
    MCAuto<MEDCouplingFieldDiscretization> _type;
  };
}

#endif
