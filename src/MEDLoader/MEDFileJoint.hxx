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

#ifndef __MEDFILEJOINT_HXX__
#define __MEDFILEJOINT_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileUtilities.txx"
#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"

#include "NormalizedGeometricTypes"

namespace MEDCoupling
{
/*!
 * \brief Joint Correspondence enumerates pairs of corresponding entities of a
 *        certain geometrical type in adjacent mesh domains.
 *        Correspondence of nodes is constructed when you specify no cell type,
 *        else Correspondence of cells is constructed.
 */
class MEDFileJointCorrespondence : public RefCountObject, public MEDFileWritable
{
public:
  MEDLOADER_EXPORT static MEDFileJointCorrespondence *New();
  MEDLOADER_EXPORT static MEDFileJointCorrespondence *New(DataArrayInt* correspondence); // nodes
  MEDLOADER_EXPORT static MEDFileJointCorrespondence *New(DataArrayInt* correspondence,  // cells
                                                          INTERP_KERNEL::NormalizedCellType loc_geo_type,
                                                          INTERP_KERNEL::NormalizedCellType rem_geo_type);
  MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
  MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
  MEDLOADER_EXPORT MEDFileJointCorrespondence *deepCopy() const;
  MEDLOADER_EXPORT MEDFileJointCorrespondence *shallowCpy() const;
  MEDLOADER_EXPORT bool isEqual(const MEDFileJointCorrespondence *other) const;
  MEDLOADER_EXPORT void setIsNodal(bool isNodal) { _is_nodal = isNodal; }
  MEDLOADER_EXPORT bool getIsNodal() const { return _is_nodal; }
  MEDLOADER_EXPORT void setLocalGeometryType(INTERP_KERNEL::NormalizedCellType type) { _loc_geo_type=type; }
  MEDLOADER_EXPORT INTERP_KERNEL::NormalizedCellType getLocalGeometryType() const { return _loc_geo_type; }
  MEDLOADER_EXPORT void setRemoteGeometryType(INTERP_KERNEL::NormalizedCellType type) { _rem_geo_type=type; }
  MEDLOADER_EXPORT INTERP_KERNEL::NormalizedCellType getRemoteGeometryType() const { return _rem_geo_type; }
  MEDLOADER_EXPORT void setCorrespondence(DataArrayInt *corr);
  MEDLOADER_EXPORT const DataArrayInt *getCorrespondence() const { return _correspondence; }
  MEDLOADER_EXPORT void write(const std::string& fileName, int mode, const std::string& localMeshName, const std::string& jointName, int order, int iteration) const;

  MEDLOADER_EXPORT std::string simpleRepr() const;
  MEDLOADER_EXPORT void writeLL(med_idt fid, const std::string& localMeshName, const std::string& jointName, int order, int iteration) const;
private:
  MEDFileJointCorrespondence();
  MEDFileJointCorrespondence(DataArrayInt*                     correspondence,
                             bool                              is_nodal = true,
                             INTERP_KERNEL::NormalizedCellType loc_geo_type = INTERP_KERNEL::NORM_ERROR,
                             INTERP_KERNEL::NormalizedCellType rem_geo_type = INTERP_KERNEL::NORM_ERROR);
private:
  bool                                           _is_nodal;
  INTERP_KERNEL::NormalizedCellType              _loc_geo_type;
  INTERP_KERNEL::NormalizedCellType              _rem_geo_type;
  MCAuto<DataArrayInt> _correspondence;
};

/*!
 * \brief Joint of one iteration holds correspondences of entities of all types
 */
class MEDFileJointOneStep : public RefCountObject, public MEDFileWritable
{
public:
  MEDLOADER_EXPORT static MEDFileJointOneStep *New(int dt=-1, int it=-1);
  MEDLOADER_EXPORT static MEDFileJointOneStep *New(const std::string& fileName, const std::string& mName, const std::string& jointName, int number=1);
  MEDLOADER_EXPORT static MEDFileJointOneStep *New(med_idt fid, const std::string& mName, const std::string& jointName, int number=1);
  MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
  MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
  MEDLOADER_EXPORT MEDFileJointOneStep *deepCopy() const;
  MEDLOADER_EXPORT MEDFileJointOneStep *shallowCpy() const;
  MEDLOADER_EXPORT bool isEqual(const MEDFileJointOneStep *other) const;
  MEDLOADER_EXPORT void setOrder(int order) { _order=order; }
  MEDLOADER_EXPORT int getOrder() const { return _order; }
  MEDLOADER_EXPORT void setIteration(int it) { _iteration=it; }
  MEDLOADER_EXPORT int getIteration() const { return _iteration; }
  MEDLOADER_EXPORT void pushCorrespondence(MEDFileJointCorrespondence* correspondence);
  MEDLOADER_EXPORT int getNumberOfCorrespondences() const;
  MEDLOADER_EXPORT MEDFileJointCorrespondence *getCorrespondenceAtPos(int i) const;

  MEDLOADER_EXPORT void write(const std::string& fileName, int mode, const std::string& localMeshName, const std::string& jointName) const;

  MEDLOADER_EXPORT std::string simpleRepr() const;
  MEDLOADER_EXPORT void writeLL(med_idt fid, const std::string& localMeshName, const std::string& jointName) const;
private:
  MEDFileJointOneStep();
  MEDFileJointOneStep(med_idt fid, const std::string& mName, const std::string& jointName, int number);
  MEDLOADER_EXPORT INTERP_KERNEL::NormalizedCellType convertGeometryType(med_geometry_type geotype);
protected:
  int _order;
  int _iteration;
private:
  std::vector<MCAuto<MEDFileJointCorrespondence> > _correspondences;
};

/*!
 * \brief Joint holds a sequence of joints of different iterations relating to
 *        a pair of mesh domains: a local one and a distant one
 */
class MEDFileJoint : public RefCountObject, public MEDFileWritableStandAlone
{
public:
    MEDLOADER_EXPORT static MEDFileJoint *New();
    MEDLOADER_EXPORT static MEDFileJoint *New(const std::string& fileName, const std::string& mName, int num);
    MEDLOADER_EXPORT static MEDFileJoint *New(med_idt fid, const std::string& mName, int num);
    MEDLOADER_EXPORT static MEDFileJoint *New(const std::string& jointName, const std::string& locMeshName, const std::string& remoteMeshName, int remoteMeshNum );
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT MEDFileJoint *deepCopy() const;
    MEDLOADER_EXPORT MEDFileJoint *shallowCpy() const;
    MEDLOADER_EXPORT bool isEqual(const MEDFileJoint *other) const;
    MEDLOADER_EXPORT void setLocalMeshName(const std::string& name) { _loc_mesh_name=name; }
    MEDLOADER_EXPORT std::string getLocalMeshName() const { return _loc_mesh_name; }
    MEDLOADER_EXPORT void setRemoteMeshName(const std::string& name) { _rem_mesh_name=name; }
    MEDLOADER_EXPORT std::string getRemoteMeshName() const { return _rem_mesh_name; }
    MEDLOADER_EXPORT void setDescription(const std::string& name) { _desc_name=name; }
    MEDLOADER_EXPORT std::string getDescription() const { return _desc_name; }
    MEDLOADER_EXPORT void setJointName(const std::string& name) { _joint_name=name; }
    MEDLOADER_EXPORT std::string getJointName() const { return _joint_name; }
    MEDLOADER_EXPORT bool changeJointNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT void setDomainNumber(const int& number) { _domain_number=number; }
    MEDLOADER_EXPORT int getDomainNumber() const { return _domain_number; }
    MEDLOADER_EXPORT void pushStep(MEDFileJointOneStep* step);
    MEDLOADER_EXPORT int getNumberOfSteps() const;
    MEDLOADER_EXPORT MEDFileJointOneStep *getStepAtPos(int i) const;

    MEDLOADER_EXPORT void writeLL(med_idt fid) const;

    MEDLOADER_EXPORT std::string simpleRepr() const;
  private:
    MEDFileJoint();
    MEDFileJoint(med_idt fid, const std::string& mName, int num);
  private:
    std::string _loc_mesh_name;
    std::string _joint_name;
    std::string _desc_name;
    int _domain_number;
    std::string _rem_mesh_name;
    std::vector< MCAuto<MEDFileJointOneStep> > _joint;
  };

  /*!
   * \brief Joints of a mesh domain relating to all other mesh domains
   */
  class MEDFileJoints : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileJoints *New();
    MEDLOADER_EXPORT static MEDFileJoints *New(const std::string& fileName, const std::string& meshName);
    MEDLOADER_EXPORT static MEDFileJoints *New(med_idt fid, const std::string& meshName);
    MEDLOADER_EXPORT MEDFileJoints *deepCopy() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT void simpleReprWithoutHeader(std::ostream& oss) const;
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT int getNumberOfJoints() const;
    MEDLOADER_EXPORT MEDFileJoint *getJointAtPos(int i) const;
    MEDLOADER_EXPORT MEDFileJoint *getJointWithName(const std::string& jname) const;
    MEDLOADER_EXPORT std::vector<std::string> getJointsNames() const;
    MEDLOADER_EXPORT bool changeJointNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    //
    MEDLOADER_EXPORT void resize(int newSize);
    MEDLOADER_EXPORT void pushJoint(MEDFileJoint *joint);
    MEDLOADER_EXPORT void setJointAtPos(int i, MEDFileJoint *joint);
    MEDLOADER_EXPORT void destroyJointAtPos(int i);
  private:
    ~MEDFileJoints() { }
    MEDFileJoints();
    MEDFileJoints(med_idt fid, const std::string& meshName);
  private:
    std::vector< MCAuto<MEDFileJoint> > _joints;
  };
}

#endif
