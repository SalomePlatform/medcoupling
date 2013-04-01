// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGFIELDOVERTIME_HXX__
#define __PARAMEDMEM_MEDCOUPLINGFIELDOVERTIME_HXX__

#include "MEDCouplingMultiFields.hxx"
#include "MEDCouplingDefinitionTime.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <vector>

namespace ParaMEDMEM
{
  class MEDCOUPLING_EXPORT MEDCouplingFieldOverTime : public MEDCouplingMultiFields
  {
  public:
    static MEDCouplingFieldOverTime *New(const std::vector<MEDCouplingFieldDouble *>& fs) throw(INTERP_KERNEL::Exception);
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    double getTimeTolerance() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    bool isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    bool isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    //void getIdsToFetch(double time, int& fieldId, int& arrId, int& meshId) const;
    //void setFieldOnId(int fieldId, MEDCouplingFieldDouble *f);
    //void dispatchPointers();
    std::vector<MEDCouplingMesh *> getMeshes() const throw(INTERP_KERNEL::Exception);
    std::vector<MEDCouplingMesh *> getDifferentMeshes(std::vector<int>& refs) const throw(INTERP_KERNEL::Exception);
    std::vector<DataArrayDouble *> getArrays() const throw(INTERP_KERNEL::Exception);
    std::vector<DataArrayDouble *> getDifferentArrays(std::vector< std::vector<int> >& refs) const throw(INTERP_KERNEL::Exception);
    MEDCouplingDefinitionTime getDefinitionTimeZone() const;
  protected:
    MEDCouplingFieldOverTime();
  private:
    MEDCouplingFieldOverTime(const std::vector<MEDCouplingFieldDouble *>& fs) throw(INTERP_KERNEL::Exception);
  };
}

#endif
