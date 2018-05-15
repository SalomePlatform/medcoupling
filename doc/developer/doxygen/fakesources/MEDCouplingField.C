// Copyright (C) 2013-2016  CEA/DEN, EDF R&D, OPEN CASCADE
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

// This file contains some code used only for
// generation of documentation for inline methods.


namespace MEDCoupling
{
  /*!
   * Checks if \a this field is correctly defined, else an exception is thrown.
   *  \throw If the mesh is not set.
   *  \throw If the data array is not set.
   *  \throw If the spatial discretization of \a this field is NULL.
   *  \throw If \a this->getTimeTolerance() < 0.
   *  \throw If the temporal discretization data is incorrect.
   *  \throw If mesh data does not correspond to field data.
   */
  void MEDCouplingField::checkCoherency() const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns the underlying mesh of \a this field.
   *  \return const MEDCoupling::MEDCouplingMesh * - a const pointer to the underlying mesh.
   */
  const MEDCoupling::MEDCouplingMesh *MEDCouplingField::getMesh() const {}
  /*!
   * Returns the description of \a this field.
   *  \return const char * - a string containing the field description.
   */
  const char *MEDCouplingField::getDescription() const {}
  /*!
   * Sets the description of \a this field.
   *  \param [in] desc - a string containing the field description.
   */
  void MEDCouplingField::setDescription(const char *desc) {}
  /*!
   * Returns the name of \a this field.
   *  \return const char * - a string containing the field name.
   */
  const char *MEDCouplingField::getName() const {}
}
