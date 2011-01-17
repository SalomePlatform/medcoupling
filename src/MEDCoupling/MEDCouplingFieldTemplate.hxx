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

#ifndef __MEDCOUPLINGFIELDTEMPLATE_HXX__
#define __MEDCOUPLINGFIELDTEMPLATE_HXX__

#include "MEDCouplingField.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingFieldDouble;

  class MEDCouplingFieldTemplate : public MEDCouplingField
  {
  public:
    static MEDCouplingFieldTemplate *New(const MEDCouplingFieldDouble *f) throw(INTERP_KERNEL::Exception);
    static MEDCouplingFieldTemplate *New(TypeOfField type);
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    //
    void getTinySerializationIntInformation(std::vector<int>& tinyInfo) const;
    void getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const;
    void getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const;
    void resizeForUnserialization(const std::vector<int>& tinyInfoI, DataArrayInt *&dataInt);
    void finishUnserialization(const std::vector<int>& tinyInfoI, const std::vector<double>& tinyInfoD, const std::vector<std::string>& tinyInfoS);
    void serialize(DataArrayInt *&dataInt) const;
    //
  private:
    MEDCouplingFieldTemplate(const MEDCouplingFieldDouble *f) throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldTemplate(TypeOfField type);
  };
}

#endif
