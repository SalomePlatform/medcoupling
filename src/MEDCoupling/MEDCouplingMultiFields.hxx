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
// Author : Anthony Geay (CEA/DEN)

#ifndef __PARAMEDMEM_MEDCOUPLINGMULTIFIELDS_HXX__
#define __PARAMEDMEM_MEDCOUPLINGMULTIFIELDS_HXX__

#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "InterpKernelException.hxx"

#include <vector>

namespace ParaMEDMEM
{
  class MEDCouplingMesh;
  class DataArrayDouble;
  class MEDCouplingFieldDouble;
  class MEDCouplingFieldTemplate;

  class MEDCOUPLING_EXPORT MEDCouplingMultiFields : public RefCountObject, public TimeLabel
  {
  public:
    static MEDCouplingMultiFields *New(const std::vector<MEDCouplingFieldDouble *>& fs) throw(INTERP_KERNEL::Exception);
    static MEDCouplingMultiFields *New();
    MEDCouplingMultiFields *deepCpy() const;
    std::string getName() const;
    std::string getDescription() const;
    std::string getTimeUnit() const;
    double getTimeResolution() const throw(INTERP_KERNEL::Exception);
    virtual std::string simpleRepr() const;
    virtual std::string advancedRepr() const;
    virtual bool isEqual(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    virtual bool isEqualWithoutConsideringStr(const MEDCouplingMultiFields *other, double meshPrec, double valsPrec) const;
    const MEDCouplingFieldDouble *getFieldWithId(int id) const throw(INTERP_KERNEL::Exception);
    std::vector<const MEDCouplingFieldDouble *> getFields() const;
    int getNumberOfFields() const;
    const MEDCouplingFieldDouble *getFieldAtPos(int id) const throw(INTERP_KERNEL::Exception);
    virtual std::vector<MEDCouplingMesh *> getMeshes() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<MEDCouplingMesh *> getDifferentMeshes(std::vector<int>& refs) const throw(INTERP_KERNEL::Exception);
    virtual std::vector<DataArrayDouble *> getArrays() const throw(INTERP_KERNEL::Exception);
    virtual std::vector<DataArrayDouble *> getDifferentArrays(std::vector< std::vector<int> >& refs) const throw(INTERP_KERNEL::Exception);
    void updateTime() const;
    void getTinySerializationInformation(std::vector<int>& tinyInfo, std::vector<double>& tinyInfo2, int& nbOfDiffMeshes, int& nbOfDiffArr) const;
    void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD,
                               const std::vector<MEDCouplingFieldTemplate *>& ft, const std::vector<MEDCouplingMesh *>& ms,
                               const std::vector<DataArrayDouble *>& das);
    virtual void checkCoherency() const throw(INTERP_KERNEL::Exception);
  protected:
    MEDCouplingMultiFields(const std::vector<MEDCouplingFieldDouble *>& fs) throw(INTERP_KERNEL::Exception);
    MEDCouplingMultiFields(const MEDCouplingMultiFields& other);
    MEDCouplingMultiFields();
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> > _fs;
  };
}

#endif

