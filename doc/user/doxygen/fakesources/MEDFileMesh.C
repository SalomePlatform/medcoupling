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
// * generation of documentation for inline methods,
// * groupping methods into "Basic API", "Advanced" and "Others..." sections


namespace MEDCoupling
{
  /*!
   * Sets the name of \a this mesh.
   *  \param [in] name - the new mesh name.
   */
  void MEDFileMesh::setName(const std::string& name) {}
  /*!
   * Returns the name of \a this mesh.
   *  \return const char* name - the mesh name.
   */
  const std::string& MEDFileMesh::getName() const {}
  /*!
   * Sets the universal name of \a this mesh. The universal name uniquely identifies the mesh.
   *  \param [in] name - the new universal mesh name.
   */
  void MEDFileMesh::setUnivName(const std::string& name) {}
  /*!
   * Returns the universal name of \a this mesh. The universal name uniquely identifies the mesh.
   *  \return const std::string&  - the universal mesh name.
   */
  const std::string& MEDFileMesh::getUnivName() const {}
  /*!
   * Sets the description of \a this mesh.
   *  \param [in] name - the new mesh description.
   */
  void MEDFileMesh::setDescription(const std::string& name) {}
  /*!
   * Returns the description of \a this mesh.
   *  \return const char* - the mesh description.
   */
  const std::string& MEDFileMesh::getDescription() const {}
  /*!
   * Sets the order number of iteration of \a this mesh state.
   *  \param [in] order - the order number.
   */
  void MEDFileMesh::setOrder(int order) {}
  /*!
   * Returns the order number of iteration of \a this mesh state.
   *  \return int - the order number.
   */
  int MEDFileMesh::getOrder() const {}
  /*!
   * Sets the number of iteration of \a this mesh state.
   *  \param [in] it - the iteration number.
   */
  void MEDFileMesh::setIteration(int it) {}
  /*!
   * Returns the number of iteration of \a this mesh state.
   *  \return int - the iteration number.
   */
  int MEDFileMesh::getIteration() const {}
  /*!
   * Sets the time of \a this mesh state.
   *  \param [in] val - the time value.
   */
  void MEDFileMesh::setTimeValue(double time) {}
  /*!
   * Sets time, the number of iteration and the order number of iteration 
   * of \a this mesh state.
   *  \param [in] val - the time value.
   *  \param [in] iteration - the iteration number.
   *  \param [in] order - the order number.
   */
  void MEDFileMesh::setTime(int dt, int it, double time) {}
  /*!
   * Returns time, the number of iteration and the order number of iteration 
   * of \a this mesh state.
   *  \param [out] iteration - the iteration number.
   *  \param [out] order - the order number.
   *  \return double - the time value.
   */
  double MEDFileMesh::getTime(int& dt, int& it) {}
  /*!
   * Returns the time of \a this mesh state.
   *  \return double - the time value.
   */
  double MEDFileMesh::getTimeValue() const {}
  /*!
   * Sets units in which the time is measured.
   *  \param [in] unit - the time unit name.
   */
  void MEDFileMesh::setTimeUnit(const std::string& unit) {}
  /*!
   * Returns units in which the time is measured.
   *  \return const std::string&  - the time unit name.
   */
  const std::string& MEDFileMesh::getTimeUnit() const {}
  /*!
   * Returns names and ids of all families in \a this mesh.
   *  \return const std::map<std::string,int>& - a map of a family name to a family id.
   */
  const std::map<std::string,int>& MEDFileMesh::getFamilyInfo() const {}
  /*!
   * Returns names of all groups and families constituting them in \a this mesh.
   *  \return const std::map<std::string, std::vector<std::string> >& - 
   *  a map of a group name to a vector of names of families constituting the group.
   */
  const std::map<std::string, std::vector<std::string> >& MEDFileMesh::getGroupInfo() const {}
  /*!
   * Returns relative dimensions of mesh entities (excluding nodes) present in \a this mesh.
   *  \return std::vector<int> - a sequence of the relative dimensions.
   */
  std::vector<int> MEDFileMesh::getNonEmptyLevels() const {}
  /*!
   * Returns relative dimensions of mesh entities (including nodes) present in \a this mesh.
   *  \return std::vector<int> - a sequence of the relative dimensions.
   */
  std::vector<int> MEDFileMesh::getNonEmptyLevelsExt() const {}
  /*!
   * Returns number of mesh entities of a given relative dimension in \a this mesh.
   *  \param [in] meshDimRelToMaxExt - the relative dimension of interest.
   *  \return int - the number of entities.
   *  \throw If no mesh entities of dimension \a meshDimRelToMaxExt are available in \a this mesh.
   */
  int MEDFileMesh::getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns a MEDCouplingMesh of a given relative dimension.
   *  \param [in] meshDimRelToMax - the relative dimension of interest.
   *  \param [in] renum - if \c true, the returned mesh is permuted according to the
   *          optional numbers of mesh entities.
   *  \return MEDCouplingMesh * - a pointer to MEDCouplingMesh that the caller is to
   *          delete using decrRef() as it is no more needed. 
   *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
   *  \throw If \a renum == \c true but permutation is impossible.
   */
  MEDCouplingMesh *MEDFileMesh::getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns the dimension on cells in \a this mesh.
   *  \return int - the mesh dimension.
   *  \throw If there are no cells in this mesh.
   */
  int MEDFileMesh::getMeshDimension() const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns a full textual description of \a this mesh.
   *  \return std::string - the string holding the mesh description.
   */
  std::string MEDFileMesh::advancedRepr() const {}
  /*!
   * Sets the family field of a given relative dimension.
   *  \param [in] meshDimRelToMaxExt - the relative dimension of entities for which
   *          the family field is set.
   *  \param [in] famArr - the array of the family field.
   *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
   *  \throw If \a famArr has an invalid size.
   */
  void MEDFileMesh::setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception) {}
  /*!
   * Sets the optional numbers of mesh entities of a given dimension.
   *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
   *  \param [in] renumArr - the array of the numbers.
   *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
   *  \throw If \a renumArr has an invalid size.
   */
  void MEDFileMesh::setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns the family field for mesh entities of a given dimension.
   *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
   *  \return const DataArrayInt * - the family field. It is an array of ids of families
   *          each mesh entity belongs to. It can be NULL.
   */
  const DataArrayInt *MEDFileMesh::getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns the optional numbers of mesh entities of a given dimension.
   *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
   *  \return const DataArrayInt * - the array of the entity numbers.
   *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
   */
  const DataArrayInt *MEDFileMesh::getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns the optional numbers of mesh entities of a given dimension transformed using
   * DataArrayInt::invertArrayN2O2O2N().
   *  \param [in] meshDimRelToMaxExt - the relative dimension of mesh entities.
   *  \return const DataArrayInt * - the array of the entity numbers transformed using
   *          DataArrayInt::invertArrayN2O2O2N().
   *  \throw If there are no mesh entities of \a meshDimRelToMaxExt dimension in \a this mesh.
   */
  const DataArrayInt *MEDFileMesh::getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) {}
  /*!
   * Returns ids of mesh entities contained in given families of a given dimension.
   *  \param [in] meshDimRelToMaxExt - a relative dimension of the mesh entities whose ids
   *          are required.
   *  \param [in] fams - the names of the families of interest.
   *  \param [in] renum - if \c true, the optional numbers of entities, if available, are
   *          returned instead of ids.
   *  \return DataArrayInt * - a new instance of DataArrayInt holding either ids or
   *          numbers, if available and required, of mesh entities of the families. The caller
   *          is to delete this array using decrRef() as it is no more needed. 
   *  \throw If the family field is missing for \a meshDimRelToMaxExt.
   */
  DataArrayInt *MEDFileMesh::getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception) {}
}


namespace MEDCoupling
{
  /*! \name Basic API   */
  ///@{
  MEDFileMesh::FindOrCreateAndGiveFamilyWithId(std::map<std::string,int>& families, int id, bool& created);
//MEDFileMesh::New(const std::string& fileName);
//MEDFileMesh::New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1);
MEDFileMesh::addFamily(const std::string& familyName, int id);
MEDFileMesh::addFamilyOnGrp(const std::string& grpName, const std::string& famName);
MEDFileMesh::advancedRepr() const = 0;
MEDFileMesh::areFamsEqual(const MEDFileMesh *other, std::string& what) const;
MEDFileMesh::areGrpsEqual(const MEDFileMesh *other, std::string& what) const;
MEDFileMesh::assignFamilyNameWithGroupName();
MEDFileMesh::changeFamilyId(int oldId, int newId);
MEDFileMesh::changeFamilyName(const std::string& oldName, const std::string& newName);
MEDFileMesh::changeGroupName(const std::string& oldName, const std::string& newName);
MEDFileMesh::copyFamGrpMapsFrom(const MEDFileMesh& other);
MEDFileMesh::createGroupOnAll(int meshDimRelToMaxExt, const std::string& groupName);
MEDFileMesh::existsFamily(const std::string& familyName) const;
MEDFileMesh::existsFamily(int famId) const;
MEDFileMesh::existsGroup(const std::string& groupName) const;
MEDFileMesh::findOrCreateAndGiveFamilyWithId(int id, bool& created);
MEDFileMesh::getDescription() const;
MEDFileMesh::getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const;
MEDFileMesh::getFamiliesIds(const std::vector<std::string>& famNames) const;
MEDFileMesh::getFamiliesIdsOnGroup(const std::string& name) const;
MEDFileMesh::getFamiliesNames() const;
MEDFileMesh::getFamiliesOnGroup(const std::string& name) const;
MEDFileMesh::getFamiliesOnGroups(const std::vector<std::string>& grps) const;
MEDFileMesh::getFamilyArr(int meshDimRelToMaxExt, const std::string& fam, bool renum=false) const;
MEDFileMesh::getFamilyFieldAtLevel(int meshDimRelToMaxExt) const;
MEDFileMesh::getFamilyId(const std::string& name) const;
MEDFileMesh::getFamilyInfo() const;
MEDFileMesh::getFamilyNameGivenId(int id) const;
MEDFileMesh::getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const;
MEDFileMesh::getGroupArr(int meshDimRelToMaxExt, const std::string& grp, bool renum=false) const;
MEDFileMesh::getGroupInfo() const;
MEDFileMesh::getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const;
MEDFileMesh::getGroupsNames() const;
MEDFileMesh::getGroupsOnFamily(const std::string& name) const;
MEDFileMesh::getIteration() const;
MEDFileMesh::getMaxFamilyId() const;
MEDFileMesh::getMeshDimension() const;
MEDFileMesh::getName() const;
MEDFileMesh::getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=false) const;
MEDFileMesh::getNodeFamilyArr(const std::string& fam, bool renum=false) const;
MEDFileMesh::getNodeGroupArr(const std::string& grp, bool renum=false) const;
MEDFileMesh::getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const;
MEDFileMesh::getNonEmptyLevels() const = 0;
MEDFileMesh::getNonEmptyLevelsExt() const = 0;
MEDFileMesh::getNumberFieldAtLevel(int meshDimRelToMaxExt) const;
MEDFileMesh::getOrder() const;
MEDFileMesh::getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const;
MEDFileMesh::getSizeAtLevel(int meshDimRelToMaxExt) const;
MEDFileMesh::getTime(int& dt, int& it);
MEDFileMesh::getTimeUnit() const;
MEDFileMesh::getTimeValue() const;
MEDFileMesh::getUnivName() const;
MEDFileMesh::isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
MEDFileMesh::keepFamIdsOnlyOnLevs(const std::vector<int>& famIds, const std::vector<int>& levs);
MEDFileMesh::removeFamily(const std::string& name);
MEDFileMesh::removeGroup(const std::string& name);
MEDFileMesh::setDescription(const std::string& name);
MEDFileMesh::setFamiliesIdsOnGroup(const std::string& name, const std::vector<int>& famIds);
MEDFileMesh::setFamiliesOnGroup(const std::string& name, const std::vector<std::string>& fams);
MEDFileMesh::setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr);
MEDFileMesh::setFamilyId(const std::string& familyName, int id);
MEDFileMesh::setFamilyInfo(const std::map<std::string,int>& info);
MEDFileMesh::setGroupInfo(const std::map<std::string, std::vector<std::string> >&info);
MEDFileMesh::setGroupsAtLevel(int meshDimRelToMaxExt, const std::vector<const DataArrayInt *>& grps, bool renum=false);
MEDFileMesh::setGroupsOnFamily(const std::string& famName, const std::vector<std::string>& grps);
MEDFileMesh::setIteration(int it);
MEDFileMesh::setName(const std::string& name);
MEDFileMesh::setOrder(int order);
MEDFileMesh::setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr);
MEDFileMesh::setTime(int dt, int it, double time);
MEDFileMesh::setTimeUnit(const std::string& unit);
MEDFileMesh::setTimeValue(double time);
MEDFileMesh::setUnivName(const std::string& name);
MEDFileMesh::simpleRepr() const;
MEDFileMesh::write(const std::string& fileName, int mode) const;
MEDFileMesh::write(med_idt fid) const;
///@}

/*! \name   Advanced API   */
///@{
MEDFileMesh::clearNonDiscrAttributes() const;
///@} 

/*! \name Others... */
///@{
MEDFileMesh::ChangeAllGroupsContainingFamily(std::map<std::string, std::vector<std::string> >& groups, const std::string& familyNameToChange, const std::vector<std::string>& newFamiliesNames);
MEDFileMesh::CreateNameNotIn(const std::string& nameTry, const std::vector<std::string>& namesToAvoid);
MEDFileMesh::MEDFileMesh();
MEDFileMesh::PutInThirdComponentOfCodeOffset(std::vector<int>& code, int strt);
MEDFileMesh::TranslateFamilyIds(int offset, DataArrayInt *famArr, std::vector< std::vector<int> >& famIdsPerGrp);
MEDFileMesh::addFamilyOnAllGroupsHaving(const std::string& famName, const std::string& otherFamName);
MEDFileMesh::appendFamilyEntries(const DataArrayInt *famIds, const std::vector< std::vector<int> >& fidsOfGrps, const std::vector<std::string>& grpNames);
MEDFileMesh::changeAllGroupsContainingFamily(const std::string& familyNameToChange, const std::vector<std::string>& newFamiliesNames);
MEDFileMesh::changeFamilyIdArr(int oldId, int newId);
MEDFileMesh::changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
MEDFileMesh::dealWithTinyInfo(const MEDCouplingMesh *m);
MEDFileMesh::deepCpy() const;
MEDFileMesh::ensureDifferentFamIdsPerLevel();
MEDFileMesh::getAllFamiliesIdsReferenced() const;
MEDFileMesh::getFamilyRepr(std::ostream& oss) const;
//MEDFileMesh::getHeapMemorySize() const;
MEDFileMesh::getMaxFamilyIdInArrays() const;
MEDFileMesh::getMinFamilyId() const;
MEDFileMesh::getMinFamilyIdInArrays() const;
MEDFileMesh::getNameFieldAtLevel(int meshDimRelToMaxExt) const;
MEDFileMesh::getNumberOfNodes() const;
MEDFileMesh::getTheMaxFamilyId() const;
MEDFileMesh::getTheMinFamilyId() const;
MEDFileMesh::normalizeFamIdsMEDFile();
MEDFileMesh::normalizeFamIdsTrio();
MEDFileMesh::setFamilyIdUnique(const std::string& familyName, int id);
MEDFileMesh::setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr);
MEDFileMesh::shallowCpy() const;
MEDFileMesh::synchronizeTinyInfoOnLeaves() const = 0;
MEDFileMesh::unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell);
MEDFileMesh::writeLL(med_idt fid) const;
int MEDFileMesh::_order;
int MEDFileMesh::_iteration;
double MEDFileMesh::_time;
std::string MEDFileMesh::_dt_unit;
std::string MEDFileMesh::_name;
std::string MEDFileMesh::_univ_name;
std::string MEDFileMesh::_desc_name;
std::map<std::string, std::vector<std::string> > MEDFileMesh::_groups;
std::map<std::string,int> MEDFileMesh::_families;
static const char MEDFileMesh::DFT_FAM_NAME[];
///@} 
}

