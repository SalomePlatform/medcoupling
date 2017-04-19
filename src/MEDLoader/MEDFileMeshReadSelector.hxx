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

#ifndef __MEDFILEMESHREADSELECTOR_HXX__
#define __MEDFILEMESHREADSELECTOR_HXX__

#include "MEDLoaderDefines.hxx"

#include <sstream>
#include <string>

namespace MEDCoupling
{
  class MEDLOADER_EXPORT MEDFileMeshReadSelector
  {
  public:
    MEDFileMeshReadSelector();
    MEDFileMeshReadSelector(unsigned int code);
    unsigned int getCode() const;
    void setCode(unsigned int newCode);
    bool isCellFamilyFieldReading() const;
    bool isNodeFamilyFieldReading() const;
    bool isCellNameFieldReading() const;
    bool isNodeNameFieldReading() const;
    bool isCellNumFieldReading() const;
    bool isNodeNumFieldReading() const;
    bool isGlobalNodeNumFieldReading() const;
    void setCellFamilyFieldReading(bool b);
    void setNodeFamilyFieldReading(bool b);
    void setCellNameFieldReading(bool b);
    void setNodeNameFieldReading(bool b);
    void setCellNumFieldReading(bool b);
    void setNodeNumFieldReading(bool b);
    void setGlobalNodeNumFieldReading(bool b);
    void reprAll(std::ostream& str) const;
  private:
    static std::string ReprStatus(bool v);
  private:
    //bit #0 cell family field
    //bit #1 node family field
    //bit #2 cell name field
    //bit #3 node name field
    //bit #4 cell num field
    //bit #5 node num field
    unsigned int _code;
  };
}

#endif
