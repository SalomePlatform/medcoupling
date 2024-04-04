// Copyright (C) 2007-2024  CEA, EDF
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

#ifndef ICoCoMEDDoubleField_included
#define ICoCoMEDDoubleField_included

#include "ICoCoField.hxx"
#include "ICoCo_DeclSpec.hxx"

namespace MEDCoupling
{
  class MEDCouplingFieldDouble;
}

namespace ICoCo
{
  /*! @brief Field data stored internally as a MEDCoupling object.
   *
   * This class is a wrapper around a MEDCoupling::MEDCouplingFieldDouble object, which holds the field data.
   * In version 2 of ICoCo, MEDCoupling::MEDCouplingFieldDouble objects are not anymore exposed directly into the API
   * of ICoCo::Problem. The rationale is to make the interface ICoCo::Problem free of external dependencies
   * (MEDCoupling particularly).
   *
   * @sa the MEDCoupling documentation, notably the reference counter mechanism used to manage the lifecycle of
   * MEDCoupling objects.
   */
  class ICOCO_EXPORT MEDDoubleField : public ICoCo::Field
  {
  public:
    /*! Builds an empty field (internal MEDCoupling object not set).
     */
    MEDDoubleField();

    /*! @brief Builds a field and assign its internal MEDCouplingField object.
     *
     * The name of 'this' is also automatically set the name of the MEDCouplingField object. If 'field' is null
     * the name is set to the empty string.
     *
     * @param field MEDCoupling field instance to use for field data. The field reference counter is incremented.
     */
    MEDDoubleField(MEDCoupling::MEDCouplingFieldDouble* field);

    /*! @brief Copy construcotr.
     */
    MEDDoubleField(const MEDDoubleField& field);

    /*! @brief Assignement operator.
     * @param field another MEDDoubleField instance. The previous internal MEDCoupling field reference (if any) has
     * its counter decremented.
     */
    MEDDoubleField& operator=(const MEDDoubleField& field);

    /*! @brief Destructor.
     */
    ~MEDDoubleField() override;

    /*! @brief Get the internal MEDCoupling field object.
     * @return a pointer to the MEDCouplingField object detained by this instance. Note that the corresponding
     * object should not be deleted, or its reference counter should not be decreased! Doing so will result in an
     * invalid instance of the current MEDField.
     */
    MEDCoupling::MEDCouplingFieldDouble *getMCField() const;

    /*! @brief Set the internal MEDCoupling field object.
     *
     * Any previously set field is discarded (its reference counter is decreased) and the reference counter of the
     * field being set is increased.
     *
     * The name of 'this' is also automatically set the name of the MEDCouplingField object. If 'field' is nullptr
     * the name is reset to the empty string.
     *
     * @param field MEDCouplingFieldDouble object to be used.
     */
    void setMCField(MEDCoupling::MEDCouplingFieldDouble * field);

  private:
    MEDCoupling::MEDCouplingFieldDouble *_field;
  };
} // namespace ICoCo

#endif
