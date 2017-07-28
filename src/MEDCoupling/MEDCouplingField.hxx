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
  class DataArrayInt;
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
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshData(const int *start, const int *end, DataArrayInt *&di) const;
    MEDCOUPLING_EXPORT MEDCouplingMesh *buildSubMeshDataRange(int begin, int end, int step, int& beginOut, int& endOut, int& stepOut, DataArrayInt *&di) const;
    MEDCOUPLING_EXPORT DataArrayInt *computeTupleIdsToSelectFromCellIds(const int *startCellIds, const int *endCellIds) const;
    MEDCOUPLING_EXPORT const MEDCouplingFieldDiscretization *getDiscretization() const { return _type; }
    MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization *getDiscretization() { return _type; }
    MEDCOUPLING_EXPORT void setDiscretization(MEDCouplingFieldDiscretization *newDisc);
    MEDCOUPLING_EXPORT int getNumberOfTuplesExpected() const;
    MEDCOUPLING_EXPORT int getNumberOfMeshPlacesExpected() const;
    // Gauss point specific methods
    MEDCOUPLING_EXPORT void setGaussLocalizationOnType(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                       const std::vector<double>& gsCoo, const std::vector<double>& wg);
    MEDCOUPLING_EXPORT void setGaussLocalizationOnCells(const int *begin, const int *end, const std::vector<double>& refCoo,
                                                        const std::vector<double>& gsCoo, const std::vector<double>& wg);
    MEDCOUPLING_EXPORT void clearGaussLocalizations();
    MEDCOUPLING_EXPORT MEDCouplingGaussLocalization& getGaussLocalization(int locId);
    MEDCOUPLING_EXPORT int getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT std::set<int> getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const;
    MEDCOUPLING_EXPORT int getNbOfGaussLocalization() const;
    MEDCOUPLING_EXPORT int getGaussLocalizationIdOfOneCell(int cellId) const;
    MEDCOUPLING_EXPORT void getCellIdsHavingGaussLocalization(int locId, std::vector<int>& cellIds) const;
    MEDCOUPLING_EXPORT const MEDCouplingGaussLocalization& getGaussLocalization(int locId) const;
    MEDCOUPLING_EXPORT void updateTime() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    // for MED file RW
    MEDCOUPLING_EXPORT int getNumberOfTuplesExpectedRegardingCode(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const;
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
