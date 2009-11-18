//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#ifndef __MEDCOUPLINGFIELDDISCRETIZATION_HXX__
#define __MEDCOUPLINGFIELDDISCRETIZATION_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"
#include "MEDCouplingNatureOfField.hxx"

namespace ParaMEDMEM
{
  class DataArrayInt;
  class MEDCouplingMesh;
  class DataArrayDouble;
  class MEDCouplingFieldDouble;

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization
  {
  public:
    static MEDCouplingFieldDiscretization *New(TypeOfField type);
    static TypeOfField getTypeOfFieldFromStringRepr(const char *repr) throw(INTERP_KERNEL::Exception);
    virtual TypeOfField getEnum() const = 0;
    virtual bool isEqual(const MEDCouplingFieldDiscretization *other) const = 0;
    virtual MEDCouplingFieldDiscretization *clone() const = 0;
    virtual const char *getStringRepr() const = 0;
    virtual int getNumberOfTuples(const MEDCouplingMesh *mesh) const = 0;
    virtual void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingFieldDouble *getWeightingField(const MEDCouplingMesh *mesh, bool isAbs) const = 0;
    virtual MEDCouplingMesh *buildSubMeshData(const int *start, const int *end, const MEDCouplingMesh *mesh, DataArrayInt *&di) const = 0;
  };

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationP0 : public MEDCouplingFieldDiscretization
  {
  public:
    TypeOfField getEnum() const;
    MEDCouplingFieldDiscretization *clone() const;
    const char *getStringRepr() const;
    bool isEqual(const MEDCouplingFieldDiscretization *other) const;
    int getNumberOfTuples(const MEDCouplingMesh *mesh) const;
    void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception);
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getWeightingField(const MEDCouplingMesh *mesh, bool isAbs) const;
    MEDCouplingMesh *buildSubMeshData(const int *start, const int *end, const MEDCouplingMesh *mesh, DataArrayInt *&di) const;
  public:
    static const char REPR[];
    static const TypeOfField TYPE;
  };

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationP1 : public MEDCouplingFieldDiscretization
  {
  public:
    TypeOfField getEnum() const;
    MEDCouplingFieldDiscretization *clone() const;
    const char *getStringRepr() const;
    bool isEqual(const MEDCouplingFieldDiscretization *other) const;
    int getNumberOfTuples(const MEDCouplingMesh *mesh) const;
    void checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception);
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getWeightingField(const MEDCouplingMesh *mesh, bool isAbs) const;
    MEDCouplingMesh *buildSubMeshData(const int *start, const int *end, const MEDCouplingMesh *mesh, DataArrayInt *&di) const;
    static DataArrayInt *invertArrayO2N2N2O(const MEDCouplingMesh *mesh, const DataArrayInt *di);
  public:
    static const char REPR[];
    static const TypeOfField TYPE;
  };
}

#endif
