//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELD_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELD_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <string>
#include <vector>

namespace ParaMEDMEM
{
  class DataArrayInt;
  class MEDCouplingMesh;
  class MEDCouplingFieldDiscretization;
  class MEDCouplingGaussLocalization;

  class MEDCOUPLING_EXPORT MEDCouplingField : public RefCountObject, public TimeLabel
  {
  public:
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool areCompatibleForMerge(const MEDCouplingField *other) const;
    virtual bool areStrictlyCompatible(const MEDCouplingField *other) const;
    virtual bool isEqual(const MEDCouplingField *other, double meshPrec, double valsPrec) const;
    void setMesh(const ParaMEDMEM::MEDCouplingMesh *mesh);
    const ParaMEDMEM::MEDCouplingMesh *getMesh() const { return _mesh; }
    void setName(const char *name) { _name=name; }
    const char *getDescription() const { return _desc.c_str(); }
    void setDescription(const char *desc) { _desc=desc; }
    const char *getName() const { return _name.c_str(); }
    TypeOfField getTypeOfField() const;
    MEDCouplingMesh *buildSubMeshData(const int *start, const int *end, DataArrayInt *&di) const;
    MEDCouplingFieldDiscretization *getDiscretization() const { return _type; }
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
  protected:
    void updateTime();
  protected:
    MEDCouplingField(TypeOfField type);
    MEDCouplingField(const MEDCouplingField& other);
    virtual ~MEDCouplingField();
  protected:
    std::string _name;
    std::string _desc;
    const MEDCouplingMesh *_mesh;
    MEDCouplingFieldDiscretization *_type;
  };
}

#endif
