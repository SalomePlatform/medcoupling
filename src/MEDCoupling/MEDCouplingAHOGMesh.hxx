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

#ifndef __MEDCOUPLINGAHOGMESH_HXX__
#define __MEDCOUPLINGAHOGMESH_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "InterpKernelException.hxx"

namespace ParaMEDMEM
{
  class MEDCouplingIMesh;
  class MEDCouplingUMesh;
  class MEDCouplingAHOGMesh;

  /// @cond INTERNAL
  class MEDCouplingAHOGPatch : public RefCountObject
  {
  public:
    MEDCouplingAHOGPatch(MEDCouplingAHOGMesh *mesh, const std::vector< std::pair<int,int> >& bottomLeftTopRight);
    // direct forward to _mesh
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithOverlap() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithoutOverlap() const;
    MEDCOUPLING_EXPORT int getMaxNumberOfLevelsRelativeToThis() const;
    MEDCOUPLING_EXPORT void addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, int factor);
    // end of direct forward to _mesh
    MEDCOUPLING_EXPORT int getNumberOfOverlapedCellsForFather() const;
    // basic set/get
    MEDCOUPLING_EXPORT const std::vector< std::pair<int,int> >& getBLTRRange() const { return _bl_tr; }
    MEDCOUPLING_EXPORT const MEDCouplingAHOGMesh *getMesh() const { return _mesh; }
  private:
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildren() const;
  private:
    //! bottom left/top right cell range relative to \a _father
    std::vector< std::pair<int,int> > _bl_tr;
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingAHOGMesh> _mesh;
  };
  /// @endcond

  /*!
   * This class is the base class dedicated to AMR using Adaptative Hierarchical Overlapped image Grid.
   * This class does \b NOT inherit from MEDCouplingMesh because this class overlaps image grid structured meshes to perform adaptative mesh refinement.
   * But this class aggregates MEDCouplingMesh instances !
   */
  class MEDCouplingAHOGMesh : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static MEDCouplingAHOGMesh *New(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                       const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCOUPLING_EXPORT int getSpaceDimension() const;
    MEDCOUPLING_EXPORT int getMaxNumberOfLevelsRelativeToThis() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsAtCurrentLevel() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithOverlap() const;
    MEDCOUPLING_EXPORT int getNumberOfCellsRecursiveWithoutOverlap() const;
    //
    MEDCOUPLING_EXPORT const MEDCouplingAHOGMesh *getFather() const;
    MEDCOUPLING_EXPORT const MEDCouplingAHOGMesh *getGodFather() const;
    MEDCOUPLING_EXPORT void addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, int factor);
    MEDCOUPLING_EXPORT int getNumberOfPatches() const;
    MEDCOUPLING_EXPORT const MEDCouplingAHOGPatch *getPatch(int patchId) const;
    //
    MEDCOUPLING_EXPORT MEDCouplingUMesh *buildUnstructured() const;
  private:
    MEDCouplingAHOGMesh(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                        const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop);
    MEDCouplingAHOGMesh(MEDCouplingAHOGMesh *father, MEDCouplingIMesh *mesh);
  protected:
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDCOUPLING_EXPORT void updateTime() const;
  private:
    MEDCouplingAHOGMesh *_father;
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingIMesh> _mesh;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingAHOGPatch> > _patches;
  };
}

#endif

