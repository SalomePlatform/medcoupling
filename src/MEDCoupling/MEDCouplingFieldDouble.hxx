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
#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDDOUBLE_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDDOUBLE_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingField.hxx"
#include "MemArray.hxx"

namespace ParaMEDMEM
{
  class MEDCOUPLING_EXPORT MEDCouplingFieldDouble : public MEDCouplingField
  {
  public:
    static MEDCouplingFieldDouble *New(TypeOfField type);
    MEDCouplingFieldDouble *clone(bool recDeepCpy) const;
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    double getIJ(int tupleId, int compoId) const { return _array->getIJ(tupleId,compoId); }
    void setArray(DataArrayDouble *array);
    DataArrayDouble *getArray() const { return _array; }
    //! \b temporary
    void applyLin(double a, double b, int compoId);
    int getNumberOfComponents() const;
    int getNumberOfTuples() const throw(INTERP_KERNEL::Exception);
    void updateTime();
  private:
    MEDCouplingFieldDouble(TypeOfField type);
    MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCpy);
    ~MEDCouplingFieldDouble();
  private:
    DataArrayDouble *_array;
  };
}

#endif
