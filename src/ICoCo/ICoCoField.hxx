// Copyright (C) 2007-2021  CEA/DEN, EDF R&D
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

// WARNING: this file is part of the official ICoCo API and should not be modified.
// The official version can be found at the following URL:
//
//    https://github.com/cea-trust-platform/icoco-coupling

#ifndef ICoCoField_included
#define ICoCoField_included
#include <string>

#include <ICoCo_DeclSpec.hxx>

namespace ICoCo
{
  /*! @brief Top abstract class defining field objects that can be exchanged via the ICoCo interface.
   *
   * The Field class holds the name of the field.
   */
  class ICOCO_EXPORT Field
  {
  public:
    /*! @brief Set the name of the field.
     * @param name name of the field
     */
    void setName(const std::string& name);

    /*! @brief Retrieves the name of the field.
     * @return name of the field.
     */
    const std::string& getName() const;

    /*!
     * @brief Retrieves the name of the field as a char *
     * @return name of the field.
     */
    const char* getCharName() const;

  protected:
    Field();
    virtual ~Field();

  private:
    std::string* _name;
  };
} // namespace ICoCo
#endif
