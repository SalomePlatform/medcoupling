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
// * generation of documentation for inline methods of MEDCouplingPointSet

namespace MEDCoupling
{

  /*!
   * Tries to use a coordinates array of \a other mesh for \a this one. If all nodes
   * of \a this mesh coincide, within a specified precision, with some nodes of \a
   * other mesh, then \a this mesh refers to the coordinates array of the \a other mesh,
   * i.e. \a this->_coords = \a other._coords. Otherwise an exception is thrown and \a
   * this remains unchanged.
   *  \param [in] other - the other mesh.
   *  \param [in] epsilon - the precision to compare node coordinates of the two meshes.
   *  \throw If the coordinates array of \a this is not set.
   *  \throw If the coordinates array of \a other is not set.
   *  \throw If not all nodes of \a this mesh are present in the \a other mesh.
   */
  void MEDCouplingPointSet::tryToShareSameCoordsPermute(const MEDCouplingPointSet& other, double epsilon) throw(INTERP_KERNEL::Exception) {}

  /*!
   * Returns a const pointer to the node coordinates array of \a this mesh \b without
   * incrementing its reference counter, thus there is no need to decrRef() it by the caller.
   */
  const DataArrayDouble *MEDCouplingPointSet::getCoords() const { return _coords; }

  /*!
   * Returns a pointer to the node coordinates array of \a this mesh \b without
   * incrementing its reference counter, thus there is no need to decrRef() it by the caller.
   */
  DataArrayDouble *MEDCouplingPointSet::getCoords() { return _coords; }


  //! size of returned tinyInfo must be always the same.
  void MEDCouplingPointSet::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const {}
}
