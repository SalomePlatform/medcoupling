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


namespace MEDCoupling // inline methods of MEDFileField1TSWithoutSDA
{
  /*!
   * Returns the number of iteration where \a this field has been calculated.
   *  \return int - the iteration number.
   */
//  int MEDFileField1TSWithoutSDA::getIteration() const {}
  /*!
   * Returns the order number of iteration where \a this field has been calculated.
   *  \return int - the order number.
   */
//  int MEDFileField1TSWithoutSDA::getOrder() const {}
  /*!
   * Returns time, number of iteration and order number of iteration when
   * \a this field has been calculated.
   *  \param [out] iteration - the iteration number.
   *  \param [out] order - the order number.
   *  \return double - the time value.
   */
//  double MEDFileField1TSWithoutSDA::getTime(int& iteration, int& order) const {}
  /*!
   * Sets time, number of iteration and order number of iteration when
   * \a this field has been calculated.
   *  \param [in] val - the time value.
   *  \param [in] iteration - the iteration number.
   *  \param [in] order - the order number.
   */
//  void MEDFileField1TSWithoutSDA::setTime(int iteration, int order, double val) {}
  /*!
   * Returns units in which the time is measured.
   *  \return const char * - the time unit name.
   */
//  const std::string& MEDFileField1TSWithoutSDA::getDtUnit() const {}
}

namespace MEDCoupling // inline methods of MEDFileFieldGlobsReal
{
  /*!
   * Returns non empty names of all used profiles. To get all profiles call getPfls().
   *  \warning If a profile is used several times, its name is returned **only once**.
   *           To have a profile name in the result each time it is used, call
   *           getPflsReallyUsedMulti().
   *  \return std::vector<std::string> - a sequence of names of used profiles.
   */
  std::vector<std::string> MEDFileFieldGlobsReal::getPflsReallyUsed() const {}
  /*!
   * Returns non empty names of all used localizations. To get all localizations call getLocs().
   *  \warning If a localization is used several times, its name is returned **only once**.
   *           To have a localization name in the result each time it is used, call
   *           getLocsReallyUsedMulti().
   *  \return std::vector<std::string> - a sequence of names of used localizations.
   */
  std::vector<std::string> MEDFileFieldGlobsReal::getLocsReallyUsed() const {}
  /*!
   * Returns non empty names of all used profiles as many times as they are used.
   *  \return std::vector<std::string> - a sequence of names of used profiles.
   */
  std::vector<std::string> MEDFileFieldGlobsReal::getPflsReallyUsedMulti() const {}
  /*!
   * Returns non empty names of all used localizations as many times as they are used.
   *  \return std::vector<std::string> - a sequence of names of used localizations.
   */
  std::vector<std::string> MEDFileFieldGlobsReal::getLocsReallyUsedMulti() const {}
  /*!
   * Replaces references to some profiles (a reference is a profile name) by references
   * to other profiles.
   * \param [in] mapOfModif - a sequence describing required replacements. Each element of
   *        this sequence is a pair whose 
   *        - the first item is a vector of profile names to replace by the second item,
   *        - the second item is a profile name to replace every profile of the first item.
   */
  void MEDFileFieldGlobsReal::changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception) {}
  /*!
   * Replaces references to some localizations (a reference is a localization name) by references
   * to other localizations.
   * \param [in] mapOfModif - a sequence describing required replacements. Each element of
   *        this sequence is a pair whose 
   *        - the first item is a vector of localization names to replace by the second item,
   *        - the second item is a localization name to replace every localization of the first
   *          item.
   */
  void MEDFileFieldGlobsReal::changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif) throw(INTERP_KERNEL::Exception) {}
}
