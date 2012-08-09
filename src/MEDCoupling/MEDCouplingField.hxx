// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELD_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELD_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "InterpKernelException.hxx"

#include <string>
#include <vector>

namespace ParaMEDMEM
{
  class DataArrayInt;
  class DataArrayDouble;
  class MEDCouplingMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingGaussLocalization;

  class MEDCOUPLING_EXPORT MEDCouplingField : public RefCountObject, public TimeLabel
  {
  public:
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool areCompatibleForMerge(const MEDCouplingField *other) const;
    virtual bool areStrictlyCompatible(const MEDCouplingField *other) const;
    virtual bool isEqualIfNotWhy(const MEDCouplingField *other, double meshPrec, double valsPrec, std::string& reason) const throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    void setMesh(const ParaMEDMEM::MEDCouplingMesh *mesh);
    const ParaMEDMEM::MEDCouplingMesh *getMesh() const { return _mesh; }
    void setName(const char *name) { _name=name; }
    const char *getDescription() const { return _desc.c_str(); }
    void setDescription(const char *desc) { _desc=desc; }
    const char *getName() const { return _name.c_str(); }
    TypeOfField getTypeOfField() const;
    NatureOfField getNature() const;
    virtual void setNature(NatureOfField nat) throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getLocalizationOfDiscr() const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *buildMeasureField(bool isAbs) const throw(INTERP_KERNEL::Exception);
    MEDCouplingMesh *buildSubMeshData(const int *start, const int *end, DataArrayInt *&di) const;
    DataArrayInt *computeTupleIdsToSelectFromCellIds(const int *startCellIds, const int *endCellIds) const;
    const MEDCouplingFieldDiscretization *getDiscretization() const { return _type; }
    MEDCouplingFieldDiscretization *getDiscretization() { return _type; }
    int getNumberOfTuplesExpected() const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshPlacesExpected() const throw(INTERP_KERNEL::Exception);
    // Gauss point specific methods
    void setGaussLocalizationOnType(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                    const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    void setGaussLocalizationOnCells(const int *begin, const int *end, const std::vector<double>& refCoo,
                                     const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception);
    void clearGaussLocalizations();
    MEDCouplingGaussLocalization& getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception);
    int getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception);
    int getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception);
    int getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception);
    void getCellIdsHavingGaussLocalization(int locId, std::vector<int>& cellIds) const throw(INTERP_KERNEL::Exception);
    const MEDCouplingGaussLocalization& getGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception);
    void updateTime() const;
  protected:
    MEDCouplingField(TypeOfField type);
    MEDCouplingField(const MEDCouplingField& other);
    MEDCouplingField(MEDCouplingFieldDiscretization *type, NatureOfField nature=NoNature);
    virtual ~MEDCouplingField();
  protected:
    std::string _name;
    std::string _desc;
    NatureOfField _nature;
    const MEDCouplingMesh *_mesh;
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDiscretization> _type;
  };
}

#endif
