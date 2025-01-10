// Copyright (C) 2007-2025  CEA, EDF
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
#include "MCAuto.hxx"
#include "MCType.hxx"

#include "BoxSplittingOptions.hxx"
#include "InterpKernelException.hxx"

namespace MEDCoupling
{
  class MEDCouplingIMesh;
  class MEDCouplingUMesh;
  class DataArrayIdType;
  class DataArrayByte;
  class DataArrayDouble;
  class MEDCoupling1SGTUMesh;
  class MEDCouplingFieldDouble;
  class MEDCouplingCartesianAMRMesh;
  class MEDCouplingCartesianAMRMeshGen;

  /// @cond INTERNAL

  /*!
   * This class does not inherit from TimeLabel so only const method should exist.
   */
  class MEDCouplingCartesianAMRPatchGen : public RefCountObject
  {
  public:
    MEDCOUPLING_EXPORT virtual MEDCouplingCartesianAMRPatchGen *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const = 0;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsRecursiveWithOverlap() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsRecursiveWithoutOverlap() const;
    MEDCOUPLING_EXPORT mcIdType getMaxNumberOfLevelsRelativeToThis() const;
    const MEDCouplingCartesianAMRMeshGen *getMesh() const { return _mesh; }
  protected:
    MEDCouplingCartesianAMRPatchGen(const MEDCouplingCartesianAMRPatchGen& other, MEDCouplingCartesianAMRMeshGen *father);
    MEDCouplingCartesianAMRPatchGen(MEDCouplingCartesianAMRMeshGen *mesh);
    const MEDCouplingCartesianAMRMeshGen *getMeshSafe() const;
    MEDCouplingCartesianAMRMeshGen *getMeshSafe();
  private:
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
  protected:
    MCAuto<MEDCouplingCartesianAMRMeshGen> _mesh;
  };

  /*!
   * This class does not inherit from TimeLabel so only const method should exist.
   */
  class MEDCouplingCartesianAMRPatch : public MEDCouplingCartesianAMRPatchGen
  {
  public:
    MEDCouplingCartesianAMRPatch(MEDCouplingCartesianAMRMeshGen *mesh, const std::vector< std::pair<mcIdType,mcIdType> >& bottomLeftTopRight);
    std::string getClassName() const override { return std::string("MEDCouplingCartesianAMRPatch"); }
    MEDCouplingCartesianAMRPatch *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const;
    // direct forward to _mesh
    MEDCOUPLING_EXPORT void addPatch(const std::vector< std::pair<mcIdType,mcIdType> >& bottomLeftTopRight, const std::vector<mcIdType>& factors);
    // end of direct forward to _mesh
    MEDCOUPLING_EXPORT mcIdType getNumberOfOverlapedCellsForFather() const;
    MEDCOUPLING_EXPORT bool isInMyNeighborhood(const MEDCouplingCartesianAMRPatch *other, mcIdType ghostLev) const;
    MEDCOUPLING_EXPORT bool isInMyNeighborhoodExt(const MEDCouplingCartesianAMRPatch *other, mcIdType ghostLev) const;
    MEDCOUPLING_EXPORT bool isInMyNeighborhoodDiffLev(const MEDCouplingCartesianAMRPatch *other, mcIdType ghostLev) const;
    // basic set/get
    const std::vector< std::pair<mcIdType,mcIdType> >& getBLTRRange() const { return _bl_tr; }
    MEDCOUPLING_EXPORT std::vector< std::pair<mcIdType,mcIdType> > getBLTRRangeRelativeToGF() const;
    MEDCOUPLING_EXPORT std::vector<mcIdType> computeCellGridSt() const;
    MEDCOUPLING_EXPORT static bool IsInMyNeighborhood(mcIdType ghostLev, const std::vector< std::pair<mcIdType,mcIdType> >& p1, const std::vector< std::pair<mcIdType,mcIdType> >& p2);
    //
    static std::vector< std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > > FindNeighborsOfSubPatchesOfSameLev(mcIdType ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2);
    static void FindNeighborsOfSubPatchesOf(mcIdType ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >& ret);
    static void UpdateNeighborsOfOneWithTwo(mcIdType ghostLev, const std::vector<mcIdType>& factors, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2);
    static void UpdateNeighborsOfOneWithTwoExt(mcIdType ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2);
    static void UpdateNeighborsOfOneWithTwoMixedLev(mcIdType ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2, bool isConservative);
  private:
    static void ComputeZonesOfTwoRelativeToOneDiffLev(mcIdType ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, std::vector< std::pair<mcIdType,mcIdType> >& p1Zone, std::vector< std::pair<mcIdType,mcIdType> >& p2Zone, std::vector<mcIdType>& factToApplyOn2);
  private:
    std::size_t getHeapMemorySizeWithoutChildren() const;
    static const MEDCouplingCartesianAMRMeshGen *FindCommonAncestor(const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, mcIdType& lev);
    static std::vector<mcIdType> ComputeOffsetFromTwoToOne(const MEDCouplingCartesianAMRMeshGen *comAncestor, mcIdType lev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2);
    static void UpdateNeighborsOfOneWithTwoInternal(mcIdType ghostLev, const std::vector<mcIdType>& factors, const std::vector< std::pair<mcIdType,mcIdType> >&p1 ,const std::vector< std::pair<mcIdType,mcIdType> >&p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2);
  public:
    static void ApplyFactorsOnCompactFrmt(std::vector< std::pair<mcIdType,mcIdType> >& partBeforeFact, const std::vector<mcIdType>& factors);
    static void ApplyAllGhostOnCompactFrmt(std::vector< std::pair<mcIdType,mcIdType> >& partBeforeFact, mcIdType ghostSize);
  private:
    MEDCouplingCartesianAMRPatch(const MEDCouplingCartesianAMRPatch& other, MEDCouplingCartesianAMRMeshGen *father);
  private:
    //! bottom left/top right cell range relative to \a _father
    std::vector< std::pair<mcIdType,mcIdType> > _bl_tr;
  };

  /*!
   * This class does not inherit from TimeLabel so only const method should exist.
   */
  class MEDCouplingCartesianAMRPatchGF : public MEDCouplingCartesianAMRPatchGen
  {
  public:
    MEDCouplingCartesianAMRPatchGF(MEDCouplingCartesianAMRMesh *mesh);
    std::string getClassName() const override { return std::string("MEDCouplingCartesianAMRPatchGF"); }
    MEDCouplingCartesianAMRPatchGF *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const;
  private:
    std::size_t getHeapMemorySizeWithoutChildren() const;
  private:
    MEDCouplingCartesianAMRPatchGF(const MEDCouplingCartesianAMRPatchGF& other, MEDCouplingCartesianAMRMeshGen *father);
  };

  /// @endcond

  /*!
   * This class is the base class dedicated to AMR using Adaptative Hierarchical Overlapped image Grid.
   * This class does \b NOT inherit from MEDCouplingMesh because this class overlaps image grid structured meshes to perform adaptative mesh refinement.
   * But this class aggregates MEDCouplingMesh instances !
   */
  class MEDCouplingCartesianAMRMeshGen : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT virtual MEDCouplingCartesianAMRMeshGen *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const = 0;
    MEDCOUPLING_EXPORT int getSpaceDimension() const;
    const std::vector<mcIdType>& getFactors() const { return _factors; }
    MEDCOUPLING_EXPORT void setFactors(const std::vector<mcIdType>& newFactors);
    MEDCOUPLING_EXPORT mcIdType getMaxNumberOfLevelsRelativeToThis() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsAtCurrentLevel() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsAtCurrentLevelGhost(mcIdType ghostLev) const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsRecursiveWithOverlap() const;
    MEDCOUPLING_EXPORT mcIdType getNumberOfCellsRecursiveWithoutOverlap() const;
    const MEDCouplingIMesh *getImageMesh() const { return _mesh; }
    //
    MEDCOUPLING_EXPORT virtual const MEDCouplingCartesianAMRMeshGen *getFather() const = 0;
    MEDCOUPLING_EXPORT virtual const MEDCouplingCartesianAMRMeshGen *getGodFather() const = 0;
    MEDCOUPLING_EXPORT virtual mcIdType getAbsoluteLevel() const = 0;
    MEDCOUPLING_EXPORT virtual void detachFromFather() = 0;
    MEDCOUPLING_EXPORT virtual std::vector< std::pair<mcIdType,mcIdType> > positionRelativeToGodFather(std::vector<mcIdType>& st) const = 0;
    MEDCOUPLING_EXPORT virtual mcIdType getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const = 0;
    MEDCOUPLING_EXPORT std::vector<mcIdType> getPositionRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRPatch *getPatchAtPosition(const std::vector<mcIdType>& pos) const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getMeshAtPosition(const std::vector<mcIdType>& pos) const;
    MEDCOUPLING_EXPORT virtual std::vector<MEDCouplingCartesianAMRPatchGen *> retrieveGridsAt(mcIdType absoluteLev) const;
    MEDCOUPLING_EXPORT void addPatch(const std::vector< std::pair<mcIdType,mcIdType> >& bottomLeftTopRight, const std::vector<mcIdType>& factors);
    MEDCOUPLING_EXPORT void removeAllPatches();
    MEDCOUPLING_EXPORT void removePatch(mcIdType patchId);
    MEDCOUPLING_EXPORT mcIdType getNumberOfPatches() const;
    MEDCOUPLING_EXPORT void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const std::vector<bool>& criterion, const std::vector<mcIdType>& factors);
    MEDCOUPLING_EXPORT void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayByte *criterion, const std::vector<mcIdType>& factors);
    MEDCOUPLING_EXPORT void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayDouble *criterion, const std::vector<mcIdType>& factors, double eps);
    MEDCOUPLING_EXPORT mcIdType getPatchIdFromChildMesh(const MEDCouplingCartesianAMRMeshGen *mesh) const;
    MEDCOUPLING_EXPORT std::vector< const MEDCouplingCartesianAMRPatch *> getPatches() const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRPatch *getPatch(mcIdType patchId) const;
    MEDCOUPLING_EXPORT bool isPatchInNeighborhoodOf(mcIdType patchId1, mcIdType patchId2, mcIdType ghostLev) const;
    MEDCOUPLING_EXPORT DataArrayDouble *createCellFieldOnPatch(mcIdType patchId, const DataArrayDouble *cellFieldOnThis) const;
    // coarse to fine
    MEDCOUPLING_EXPORT void fillCellFieldOnPatch(mcIdType patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, bool isConservative=true) const;
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchGhost(mcIdType patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, mcIdType ghostLev, bool isConservative=true) const;
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchOnlyOnGhostZone(mcIdType patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, mcIdType ghostLev) const;
    // coarse to fine + fine to fine
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchGhostAdv(mcIdType patchId, const DataArrayDouble *cellFieldOnThis, mcIdType ghostLev, const std::vector<const DataArrayDouble *>& arrsOnPatches, bool isConservative=true) const;
    // fine to fine
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchOnlyGhostAdv(mcIdType patchId, mcIdType ghostLev, const std::vector<const DataArrayDouble *>& arrsOnPatches) const;
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchOnlyOnGhostZoneWith(mcIdType ghostLev, const MEDCouplingCartesianAMRPatch *patchToBeModified, const MEDCouplingCartesianAMRPatch *neighborPatch, DataArrayDouble *cellFieldOnPatch, const DataArrayDouble *cellFieldNeighbor) const;
    // fine to coarse
    MEDCOUPLING_EXPORT void fillCellFieldComingFromPatch(mcIdType patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, bool isConservative=true) const;
    MEDCOUPLING_EXPORT void fillCellFieldComingFromPatchGhost(mcIdType patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, mcIdType ghostLev, bool isConservative=true) const;
    //
    MEDCOUPLING_EXPORT DataArrayIdType *findPatchesInTheNeighborhoodOf(mcIdType patchId, mcIdType ghostLev) const;
    //
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *buildMeshFromPatchEnvelop() const;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *buildMeshOfDirectChildrenOnly() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildCellFieldOnRecurseWithoutOverlapWithoutGhost(mcIdType ghostSz, const std::vector<const DataArrayDouble *>& recurseArrs) const;
    MEDCOUPLING_EXPORT DataArrayDouble *extractGhostFrom(mcIdType ghostSz, const DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT std::vector<mcIdType> getPatchIdsInTheNeighborhoodOf(mcIdType patchId, mcIdType ghostLev) const;
    MEDCOUPLING_EXPORT std::string buildPythonDumpOfThis() const;
  protected:
    MEDCouplingCartesianAMRMeshGen(const MEDCouplingCartesianAMRMeshGen& other);
    MEDCouplingCartesianAMRMeshGen(const std::string& meshName, int spaceDim, const mcIdType *nodeStrctStart, const mcIdType *nodeStrctStop,
                                   const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCouplingCartesianAMRMeshGen(MEDCouplingIMesh *mesh);
    void checkPatchId(mcIdType patchId) const;
    void checkFactorsAndIfNotSetAssign(const std::vector<mcIdType>& factors);
    void retrieveGridsAtInternal(mcIdType lev, std::vector< MCAuto<MEDCouplingCartesianAMRPatchGen> >& grids) const;
    static mcIdType GetGhostLevelInFineRef(mcIdType ghostLev, const std::vector<mcIdType>& factors);
    std::vector<const DataArrayDouble *> extractSubTreeFromGlobalFlatten(const MEDCouplingCartesianAMRMeshGen *head, const std::vector<const DataArrayDouble *>& all) const;
    void dumpPatchesOf(const std::string& varName, std::ostream& oss) const;
  public:
    virtual void getPositionRelativeToInternal(const MEDCouplingCartesianAMRMeshGen *ref, std::vector<mcIdType>& ret) const = 0;
  protected:
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void updateTime() const;
  protected:
    MCAuto<MEDCouplingIMesh> _mesh;
    std::vector< MCAuto<MEDCouplingCartesianAMRPatch> > _patches;
    std::vector<mcIdType> _factors;
  };

  class MEDCouplingCartesianAMRMeshSub : public MEDCouplingCartesianAMRMeshGen
  {
  public:
    MEDCouplingCartesianAMRMeshSub(MEDCouplingCartesianAMRMeshGen *father, MEDCouplingIMesh *mesh);
    std::string getClassName() const override { return std::string("MEDCouplingCartesianAMRMeshSub"); }
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getFather() const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getGodFather() const;
    MEDCOUPLING_EXPORT mcIdType getAbsoluteLevel() const;
    MEDCOUPLING_EXPORT void detachFromFather();
    MEDCOUPLING_EXPORT std::vector< std::pair<mcIdType,mcIdType> > positionRelativeToGodFather(std::vector<mcIdType>& st) const;
    MEDCOUPLING_EXPORT mcIdType getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const;
  private:
    MEDCouplingCartesianAMRMeshSub(const MEDCouplingCartesianAMRMeshSub& other, MEDCouplingCartesianAMRMeshGen *father);
    MEDCouplingCartesianAMRMeshSub *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const;
    void getPositionRelativeToInternal(const MEDCouplingCartesianAMRMeshGen *ref, std::vector<mcIdType>& ret) const;
  protected:
    MEDCouplingCartesianAMRMeshGen *_father;
  };

  class MEDCouplingCartesianAMRMesh : public MEDCouplingCartesianAMRMeshGen
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingCartesianAMRMesh *New(const std::string& meshName, int spaceDim, const mcIdType *nodeStrctStart, const mcIdType *nodeStrctStop,
                                                               const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCOUPLING_EXPORT static MEDCouplingCartesianAMRMesh *New(MEDCouplingIMesh *mesh);
    std::string getClassName() const override { return std::string("MEDCouplingCartesianAMRMesh"); }
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getFather() const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getGodFather() const;
    MEDCOUPLING_EXPORT mcIdType getAbsoluteLevel() const;
    MEDCOUPLING_EXPORT void detachFromFather();
    MEDCOUPLING_EXPORT std::vector< std::pair<mcIdType,mcIdType> > positionRelativeToGodFather(std::vector<mcIdType>& st) const;
    MEDCOUPLING_EXPORT mcIdType getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const;
    MEDCOUPLING_EXPORT std::vector<MEDCouplingCartesianAMRPatchGen *> retrieveGridsAt(mcIdType absoluteLev) const;
    MEDCouplingCartesianAMRMesh *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const;
    MEDCOUPLING_EXPORT void createPatchesFromCriterionML(const std::vector<const INTERP_KERNEL::BoxSplittingOptions *>& bso, const DataArrayDouble *criterion, const std::vector< std::vector<mcIdType> >& factors, double eps);
  private:
    void getPositionRelativeToInternal(const MEDCouplingCartesianAMRMeshGen *ref, std::vector<mcIdType>& ret) const;
    MEDCouplingCartesianAMRMesh(const MEDCouplingCartesianAMRMesh& other);
    MEDCouplingCartesianAMRMesh(const std::string& meshName, int spaceDim, const mcIdType *nodeStrctStart, const mcIdType *nodeStrctStop,
                                const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCouplingCartesianAMRMesh(MEDCouplingIMesh *mesh);
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    ~MEDCouplingCartesianAMRMesh();
  };
}

#endif

