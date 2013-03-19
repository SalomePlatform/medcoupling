// This file contains some code used only for
// * generation of documentation for inline methods of MEDCouplingPointSet

namespace ParaMEDMEM
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

  //! This method returns directly the array in \a this \b without incrementing ref counter. The pointer is dealed by the mesh. The caller should not deal (decrRef) with this pointer
  const DataArrayDouble *MEDCouplingPointSet::getCoords() const { return _coords; }

  //! This method returns directly the array in \a this \b without incrementing ref counter. The pointer is dealed by the mesh. The caller should not deal (decrRef) with this pointer
  DataArrayDouble *MEDCouplingPointSet::getCoords() { return _coords; }


  //! size of returned tinyInfo must be always the same.
  void MEDCouplingPointSet::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const {}
}
