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
  class DataArrayInt;
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
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithOverlap() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithoutOverlap() const;
    MEDCOUPLING_EXPORT int getMaxNumberOfLevelsRelativeToThis() const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getMesh() const { return _mesh; }
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
    MEDCouplingCartesianAMRPatch(MEDCouplingCartesianAMRMeshGen *mesh, const std::vector< std::pair<int,int> >& bottomLeftTopRight);
    MEDCouplingCartesianAMRPatch *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const;
    // direct forward to _mesh
    MEDCOUPLING_EXPORT void addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors);
    // end of direct forward to _mesh
    MEDCOUPLING_EXPORT int getNumberOfOverlapedCellsForFather() const;
    MEDCOUPLING_EXPORT bool isInMyNeighborhood(const MEDCouplingCartesianAMRPatch *other, int ghostLev) const;
    MEDCOUPLING_EXPORT bool isInMyNeighborhoodExt(const MEDCouplingCartesianAMRPatch *other, int ghostLev) const;
    MEDCOUPLING_EXPORT bool isInMyNeighborhoodDiffLev(const MEDCouplingCartesianAMRPatch *other, int ghostLev) const;
    // basic set/get
    MEDCOUPLING_EXPORT const std::vector< std::pair<int,int> >& getBLTRRange() const { return _bl_tr; }
    MEDCOUPLING_EXPORT std::vector< std::pair<int,int> > getBLTRRangeRelativeToGF() const;
    MEDCOUPLING_EXPORT std::vector<int> computeCellGridSt() const;
    MEDCOUPLING_EXPORT static bool IsInMyNeighborhood(int ghostLev, const std::vector< std::pair<int,int> >& p1, const std::vector< std::pair<int,int> >& p2);
    //
    static std::vector< std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > > FindNeighborsOfSubPatchesOfSameLev(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2);
    static void FindNeighborsOfSubPatchesOf(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> >& ret);
    static void UpdateNeighborsOfOneWithTwo(int ghostLev, const std::vector<int>& factors, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2);
    static void UpdateNeighborsOfOneWithTwoExt(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2);
    static void UpdateNeighborsOfOneWithTwoMixedLev(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2, bool isConservative);
  private:
    static void ComputeZonesOfTwoRelativeToOneDiffLev(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, std::vector< std::pair<int,int> >& p1Zone, std::vector< std::pair<int,int> >& p2Zone, std::vector<int>& factToApplyOn2);
  private:
    std::size_t getHeapMemorySizeWithoutChildren() const;
    static const MEDCouplingCartesianAMRMeshGen *FindCommonAncestor(const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, int& lev);
    static std::vector<int> ComputeOffsetFromTwoToOne(const MEDCouplingCartesianAMRMeshGen *comAncestor, int lev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2);
    static void UpdateNeighborsOfOneWithTwoInternal(int ghostLev, const std::vector<int>& factors, const std::vector< std::pair<int,int> >&p1 ,const std::vector< std::pair<int,int> >&p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2);
  public:
    static void ApplyFactorsOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, const std::vector<int>& factors);
    static void ApplyAllGhostOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, int ghostSize);
  private:
    MEDCouplingCartesianAMRPatch(const MEDCouplingCartesianAMRPatch& other, MEDCouplingCartesianAMRMeshGen *father);
  private:
    //! bottom left/top right cell range relative to \a _father
    std::vector< std::pair<int,int> > _bl_tr;
  };

  /*!
   * This class does not inherit from TimeLabel so only const method should exist.
   */
  class MEDCouplingCartesianAMRPatchGF : public MEDCouplingCartesianAMRPatchGen
  {
  public:
    MEDCouplingCartesianAMRPatchGF(MEDCouplingCartesianAMRMesh *mesh);
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
    MEDCOUPLING_EXPORT const std::vector<int>& getFactors() const { return _factors; }
    MEDCOUPLING_EXPORT void setFactors(const std::vector<int>& newFactors);
    MEDCOUPLING_EXPORT int getMaxNumberOfLevelsRelativeToThis() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsAtCurrentLevel() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsAtCurrentLevelGhost(int ghostLev) const;
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithOverlap() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithoutOverlap() const;
    MEDCOUPLING_EXPORT const MEDCouplingIMesh *getImageMesh() const { return _mesh; }
    //
    MEDCOUPLING_EXPORT virtual const MEDCouplingCartesianAMRMeshGen *getFather() const = 0;
    MEDCOUPLING_EXPORT virtual const MEDCouplingCartesianAMRMeshGen *getGodFather() const = 0;
    MEDCOUPLING_EXPORT virtual int getAbsoluteLevel() const = 0;
    MEDCOUPLING_EXPORT virtual void detachFromFather() = 0;
    MEDCOUPLING_EXPORT virtual std::vector< std::pair<int,int> > positionRelativeToGodFather(std::vector<int>& st) const = 0;
    MEDCOUPLING_EXPORT virtual int getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const = 0;
    MEDCOUPLING_EXPORT std::vector<int> getPositionRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRPatch *getPatchAtPosition(const std::vector<int>& pos) const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getMeshAtPosition(const std::vector<int>& pos) const;
    MEDCOUPLING_EXPORT virtual std::vector<MEDCouplingCartesianAMRPatchGen *> retrieveGridsAt(int absoluteLev) const;
    MEDCOUPLING_EXPORT void addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors);
    MEDCOUPLING_EXPORT void removeAllPatches();
    MEDCOUPLING_EXPORT void removePatch(int patchId);
    MEDCOUPLING_EXPORT int getNumberOfPatches() const;
    MEDCOUPLING_EXPORT void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const std::vector<bool>& criterion, const std::vector<int>& factors);
    MEDCOUPLING_EXPORT void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayByte *criterion, const std::vector<int>& factors);
    MEDCOUPLING_EXPORT void createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayDouble *criterion, const std::vector<int>& factors, double eps);
    MEDCOUPLING_EXPORT int getPatchIdFromChildMesh(const MEDCouplingCartesianAMRMeshGen *mesh) const;
    MEDCOUPLING_EXPORT std::vector< const MEDCouplingCartesianAMRPatch *> getPatches() const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRPatch *getPatch(int patchId) const;
    MEDCOUPLING_EXPORT bool isPatchInNeighborhoodOf(int patchId1, int patchId2, int ghostLev) const;
    MEDCOUPLING_EXPORT DataArrayDouble *createCellFieldOnPatch(int patchId, const DataArrayDouble *cellFieldOnThis) const;
    // coarse to fine
    MEDCOUPLING_EXPORT void fillCellFieldOnPatch(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, bool isConservative=true) const;
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchGhost(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, int ghostLev, bool isConservative=true) const;
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchOnlyOnGhostZone(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, int ghostLev) const;
    // coarse to fine + fine to fine
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchGhostAdv(int patchId, const DataArrayDouble *cellFieldOnThis, int ghostLev, const std::vector<const DataArrayDouble *>& arrsOnPatches, bool isConservative=true) const;
    // fine to fine
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchOnlyGhostAdv(int patchId, int ghostLev, const std::vector<const DataArrayDouble *>& arrsOnPatches) const;
    MEDCOUPLING_EXPORT void fillCellFieldOnPatchOnlyOnGhostZoneWith(int ghostLev, const MEDCouplingCartesianAMRPatch *patchToBeModified, const MEDCouplingCartesianAMRPatch *neighborPatch, DataArrayDouble *cellFieldOnPatch, const DataArrayDouble *cellFieldNeighbor) const;
    // fine to coarse
    MEDCOUPLING_EXPORT void fillCellFieldComingFromPatch(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, bool isConservative=true) const;
    MEDCOUPLING_EXPORT void fillCellFieldComingFromPatchGhost(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, int ghostLev, bool isConservative=true) const;
    //
    MEDCOUPLING_EXPORT DataArrayInt *findPatchesInTheNeighborhoodOf(int patchId, int ghostLev) const;
    //
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *buildMeshFromPatchEnvelop() const;
    MEDCOUPLING_EXPORT MEDCoupling1SGTUMesh *buildMeshOfDirectChildrenOnly() const;
    MEDCOUPLING_EXPORT MEDCouplingFieldDouble *buildCellFieldOnRecurseWithoutOverlapWithoutGhost(int ghostSz, const std::vector<const DataArrayDouble *>& recurseArrs) const;
    MEDCOUPLING_EXPORT DataArrayDouble *extractGhostFrom(int ghostSz, const DataArrayDouble *arr) const;
    MEDCOUPLING_EXPORT std::vector<int> getPatchIdsInTheNeighborhoodOf(int patchId, int ghostLev) const;
    MEDCOUPLING_EXPORT std::string buildPythonDumpOfThis() const;
  protected:
    MEDCouplingCartesianAMRMeshGen(const MEDCouplingCartesianAMRMeshGen& other);
    MEDCouplingCartesianAMRMeshGen(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                   const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCouplingCartesianAMRMeshGen(MEDCouplingIMesh *mesh);
    void checkPatchId(int patchId) const;
    void checkFactorsAndIfNotSetAssign(const std::vector<int>& factors);
    void retrieveGridsAtInternal(int lev, std::vector< MCAuto<MEDCouplingCartesianAMRPatchGen> >& grids) const;
    static int GetGhostLevelInFineRef(int ghostLev, const std::vector<int>& factors);
    std::vector<const DataArrayDouble *> extractSubTreeFromGlobalFlatten(const MEDCouplingCartesianAMRMeshGen *head, const std::vector<const DataArrayDouble *>& all) const;
    void dumpPatchesOf(const std::string& varName, std::ostream& oss) const;
  public:
    virtual void getPositionRelativeToInternal(const MEDCouplingCartesianAMRMeshGen *ref, std::vector<int>& ret) const = 0;
  protected:
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void updateTime() const;
  protected:
    MCAuto<MEDCouplingIMesh> _mesh;
    std::vector< MCAuto<MEDCouplingCartesianAMRPatch> > _patches;
    std::vector<int> _factors;
  };

  class MEDCouplingCartesianAMRMeshSub : public MEDCouplingCartesianAMRMeshGen
  {
  public:
    MEDCouplingCartesianAMRMeshSub(MEDCouplingCartesianAMRMeshGen *father, MEDCouplingIMesh *mesh);
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getFather() const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getGodFather() const;
    MEDCOUPLING_EXPORT int getAbsoluteLevel() const;
    MEDCOUPLING_EXPORT void detachFromFather();
    MEDCOUPLING_EXPORT std::vector< std::pair<int,int> > positionRelativeToGodFather(std::vector<int>& st) const;
    MEDCOUPLING_EXPORT int getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const;
  private:
    MEDCouplingCartesianAMRMeshSub(const MEDCouplingCartesianAMRMeshSub& other, MEDCouplingCartesianAMRMeshGen *father);
    MEDCouplingCartesianAMRMeshSub *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const;
    void getPositionRelativeToInternal(const MEDCouplingCartesianAMRMeshGen *ref, std::vector<int>& ret) const;
  protected:
    MEDCouplingCartesianAMRMeshGen *_father;
  };

  class MEDCouplingCartesianAMRMesh : public MEDCouplingCartesianAMRMeshGen
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingCartesianAMRMesh *New(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                               const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCOUPLING_EXPORT static MEDCouplingCartesianAMRMesh *New(MEDCouplingIMesh *mesh);
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getFather() const;
    MEDCOUPLING_EXPORT const MEDCouplingCartesianAMRMeshGen *getGodFather() const;
    MEDCOUPLING_EXPORT int getAbsoluteLevel() const;
    MEDCOUPLING_EXPORT void detachFromFather();
    MEDCOUPLING_EXPORT std::vector< std::pair<int,int> > positionRelativeToGodFather(std::vector<int>& st) const;
    MEDCOUPLING_EXPORT int getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const;
    MEDCOUPLING_EXPORT std::vector<MEDCouplingCartesianAMRPatchGen *> retrieveGridsAt(int absoluteLev) const;
    MEDCouplingCartesianAMRMesh *deepCopy(MEDCouplingCartesianAMRMeshGen *father) const;
    MEDCOUPLING_EXPORT void createPatchesFromCriterionML(const std::vector<const INTERP_KERNEL::BoxSplittingOptions *>& bso, const DataArrayDouble *criterion, const std::vector< std::vector<int> >& factors, double eps);
  private:
    void getPositionRelativeToInternal(const MEDCouplingCartesianAMRMeshGen *ref, std::vector<int>& ret) const;
    MEDCouplingCartesianAMRMesh(const MEDCouplingCartesianAMRMesh& other);
    MEDCouplingCartesianAMRMesh(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCouplingCartesianAMRMesh(MEDCouplingIMesh *mesh);
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    ~MEDCouplingCartesianAMRMesh();
  };
}

#endif

