// Copyright (C) 2007-2014  CEA/DEN, EDF R&D
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
// Author : Anthony Geay

#ifndef __MEDCOUPLINGCARTESIANAMRMESH_HXX__
#define __MEDCOUPLINGCARTESIANAMRMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "BoxSplittingOptions.hxx"
#include "InterpKernelException.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingIMesh;
  class MEDCouplingUMesh;
  class DataArrayInt;
  class DataArrayByte;
  class DataArrayDouble;
  class MEDCoupling1SGTUMesh;
  class MEDCouplingCartesianAMRMesh;

  /// @cond INTERNAL
  class MEDCouplingCartesianAMRPatch : public RefCountObject
  {
  public:
    MEDCouplingCartesianAMRPatch(MEDCouplingCartesianAMRMesh *mesh, const std::vector< std::pair<int,int> >& bottomLeftTopRight);
    // direct forward to _mesh
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithOverlap() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithoutOverlap() const;
    MEDCOUPLING_EXPORT int getMaxNumberOfLevelsRelativeToThis() const;
    MEDCOUPLING_EXPORT void addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors);
    // end of direct forward to _mesh
    MEDCOUPLING_EXPORT int getNumberOfOverlapedCellsForFather() const;
    MEDCOUPLING_EXPORT bool isInMyNeighborhood(const MEDCouplingCartesianAMRPatch *other, int ghostLev) const;
    // basic set/get
    MEDCOUPLING_EXPORT const std::vector< std::pair<int,int> >& getBLTRRange() const { return _bl_tr; }
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMesh *getMesh() const { return _mesh; }
  private:
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
  private:
    //! bottom left/top right cell range relative to \a _father
    std::vector< std::pair<int,int> > _bl_tr;
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRMesh> _mesh;
  };
  /// @endcond

  /*!
   * This class is the base class dedicated to AMR using Adaptative Hierarchical Overlapped image Grid.
   * This class does \b NOT inherit from MEDCouplingMesh because this class overlaps image grid structured meshes to perform adaptative mesh refinement.
   * But this class aggregates MEDCouplingMesh instances !
   */
  class MEDCouplingCartesianAMRMesh : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingCartesianAMRMesh *New(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                               const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCOUPLING_EXPORT int getSpaceDimension() const;
    MEDCOUPLING_EXPORT const std::vector<int>& getFactors() const { return _factors; }
    MEDCOUPLING_EXPORT void setFactors(const std::vector<int>& newFactors);
    MEDCOUPLING_EXPORT int getMaxNumberOfLevelsRelativeToThis() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsAtCurrentLevel() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithOverlap() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithoutOverlap() const;
    MEDCOUPLING_EXPORT const MEDCouplingIMesh *getImageMesh() const { return _mesh; }
    //
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMesh *getFather() const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMesh *getGodFather() const;
    MEDCOUPLING_EXPORT void detachFromFather();
    MEDCOUPLING_EXPORT void addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors);
    MEDCOUPLING_EXPORT void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const std::vector<bool>& criterion, const std::vector<int>& factors);
    MEDCOUPLING_EXPORT void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayByte *criterion, const std::vector<int>& factors);
    MEDCOUPLING_EXPORT void removeAllPatches();
    MEDCOUPLING_EXPORT void removePatch(int patchId);
    MEDCOUPLING_EXPORT int getNumberOfPatches() const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRPatch *getPatch(int patchId) const;
    MEDCOUPLING_EXPORT bool isPatchInNeighborhoodOf(int patchId1, int patchId2, int ghostLev) const;
    MEDCOUPLING_EXPORT DataArrayDouble *createCellFieldOnPatch(int patchId, const DataArrayDouble *cellFieldOnThis) const;
    MEDCOUPLING_EXPORT void fillCellFieldOnPatch(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch) const;
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchGhost(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, int ghostLev) const;
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchGhostAdv(int patchId, const DataArrayDouble *cellFieldOnThis, int ghostLev, const std::vector<const DataArrayDouble *>& arrsOnPatches) const;
    MEDCOUPLING_EXPORT void fillCellFieldComingFromPatch(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis) const;
    MEDCOUPLING_EXPORT void fillCellFieldComingFromPatchGhost(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, int ghostLev) const;
    MEDCOUPLING_EXPORT DataArrayInt *findPatchesInTheNeighborhoodOf(int patchId, int ghostLev) const;
    //
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *buildMeshFromPatchEnvelop() const;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *buildMeshOfDirectChildrenOnly() const;
  private:
    MEDCouplingCartesianAMRMesh(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCouplingCartesianAMRMesh(MEDCouplingCartesianAMRMesh *father, MEDCouplingIMesh *mesh);
    void checkPatchId(int patchId) const;
    void checkFactorsAndIfNotSetAssign(const std::vector<int>& factors);
    static void ApplyFactorsOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, const std::vector<int>& factors);
    static void ApplyGhostOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, int ghostSize);
    static void ApplyAllGhostOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, int ghostSize);
  protected:
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDCOUPLING_EXPORT void updateTime() const;
  private:
    MEDCouplingCartesianAMRMesh *_father;
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingIMesh> _mesh;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> > _patches;
    std::vector<int> _factors;
  };
}

#endif

