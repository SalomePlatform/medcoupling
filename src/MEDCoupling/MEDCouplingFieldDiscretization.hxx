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
#include "RefCountObject.hxx"
#include "InterpKernelException.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingMesh;
  class DataArrayDouble;
  class MEDCouplingFieldDouble;

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretization
  {
  public:
    static MEDCouplingFieldDiscretization *New(TypeOfField type);
    virtual TypeOfField getEnum() const = 0;
    virtual MEDCouplingFieldDiscretization *clone() const = 0;
    virtual const char *getStringRepr() const = 0;
    virtual int getNumberOfTuples(const MEDCouplingMesh *mesh) const = 0;
    virtual void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingFieldDouble *getWeightingField(const MEDCouplingMesh *mesh) const = 0;
  };

  class MEDCOUPLING_EXPORT MEDCouplingFieldDiscretizationP0 : public MEDCouplingFieldDiscretization
  {
  public:
    TypeOfField getEnum() const;
    MEDCouplingFieldDiscretization *clone() const;
    const char *getStringRepr() const;
    int getNumberOfTuples(const MEDCouplingMesh *mesh) const;
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getWeightingField(const MEDCouplingMesh *mesh) const;
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
    int getNumberOfTuples(const MEDCouplingMesh *mesh) const;
    void checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getWeightingField(const MEDCouplingMesh *mesh) const;
  public:
    static const char REPR[];
    static const TypeOfField TYPE;
  };
}

#endif
