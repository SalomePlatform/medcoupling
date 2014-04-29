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

#include "MEDCouplingAHOGMesh.hxx"
#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingUMesh.hxx"

#include <limits>
#include <sstream>

using namespace ParaMEDMEM;

/// @cond INTERNAL

/*!
 * \param [in] mesh not null pointer of refined mesh replacing the cell range of \a father defined by the bottom left and top right just after.
 * \param [in] bottomLeftTopRight a vector equal to the space dimension of \a mesh that specifies for each dimension, the included cell start of the range for the first element of the pair,
 *                                a the end cell (\b excluded) of the range for the second element of the pair.
 */
MEDCouplingAHOGPatch::MEDCouplingAHOGPatch(MEDCouplingAHOGMesh *mesh, const std::vector< std::pair<int,int> >& bottomLeftTopRight)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("EDCouplingAHOGPatch constructor : input mesh is NULL !");
  _mesh=mesh; _mesh->incrRef();
  int dim((int)bottomLeftTopRight.size());
  if(dim!=_mesh->getSpaceDimension())
    throw INTERP_KERNEL::Exception("MEDCouplingAHOGPatch constructor : space dimension of father and input bottomLeft/topRight size mismatches !");
  _bl_tr=bottomLeftTopRight;
}

int MEDCouplingAHOGPatch::getNumberOfCellsRecursiveWithOverlap() const
{
  return _mesh->getNumberOfCellsRecursiveWithOverlap();
}

int MEDCouplingAHOGPatch::getNumberOfCellsRecursiveWithoutOverlap() const
{
  return _mesh->getNumberOfCellsRecursiveWithoutOverlap();
}

int MEDCouplingAHOGPatch::getMaxNumberOfLevelsRelativeToThis() const
{
  return _mesh->getMaxNumberOfLevelsRelativeToThis();
}

void MEDCouplingAHOGPatch::addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, int factor)
{
  return _mesh->addPatch(bottomLeftTopRight,factor);
}

int MEDCouplingAHOGPatch::getNumberOfOverlapedCellsForFather() const
{
  return MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt(_bl_tr);
}

std::size_t MEDCouplingAHOGPatch::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(sizeof(MEDCouplingAHOGPatch));
  ret+=_bl_tr.capacity()*sizeof(std::pair<int,int>);
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingAHOGPatch::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  if((const MEDCouplingAHOGMesh *)_mesh)
    ret.push_back((const MEDCouplingAHOGMesh *)_mesh);
  return ret;
}

/// @endcond


MEDCouplingAHOGMesh *MEDCouplingAHOGMesh::New(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                              const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop)
{
  return new MEDCouplingAHOGMesh(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop);
}

int MEDCouplingAHOGMesh::getSpaceDimension() const
{
  return _mesh->getSpaceDimension();
}

int MEDCouplingAHOGMesh::getMaxNumberOfLevelsRelativeToThis() const
{
  int ret(1);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingAHOGPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    ret=std::max(ret,(*it)->getMaxNumberOfLevelsRelativeToThis()+1);
  return ret;
}

int MEDCouplingAHOGMesh::getNumberOfCellsAtCurrentLevel() const
{
  return _mesh->getNumberOfCells();
}

int MEDCouplingAHOGMesh::getNumberOfCellsRecursiveWithOverlap() const
{
  int ret(_mesh->getNumberOfCells());
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingAHOGPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      ret+=(*it)->getNumberOfCellsRecursiveWithOverlap();
    }
  return ret;
}

int MEDCouplingAHOGMesh::getNumberOfCellsRecursiveWithoutOverlap() const
{
  int ret(_mesh->getNumberOfCells());
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingAHOGPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      ret-=(*it)->getNumberOfOverlapedCellsForFather();
      ret+=(*it)->getNumberOfCellsRecursiveWithoutOverlap();
    }
  return ret;
}

const MEDCouplingAHOGMesh *MEDCouplingAHOGMesh::getFather() const
{
  return _father;
}

const MEDCouplingAHOGMesh *MEDCouplingAHOGMesh::getGodFather() const
{
  if(_father==0)
    return this;
  else
    return _father->getGodFather();
}

/*!
 * \param [in] bottomLeftTopRight a vector equal to the space dimension of \a mesh that specifies for each dimension, the included cell start of the range for the first element of the pair,
 *                                a the end cell (\b excluded) of the range for the second element of the pair.
 * \param [in] factor The != 0 factor of refinement.
 */
void MEDCouplingAHOGMesh::addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, int factor)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingIMesh> mesh(static_cast<MEDCouplingIMesh *>(_mesh->buildStructuredSubPart(bottomLeftTopRight)));
  mesh->refineWithFactor(factor);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingAHOGMesh> zeMesh(new MEDCouplingAHOGMesh(this,mesh));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingAHOGPatch> elt(new MEDCouplingAHOGPatch(zeMesh,bottomLeftTopRight));
  _patches.push_back(elt);
}

int MEDCouplingAHOGMesh::getNumberOfPatches() const
{
  return (int)_patches.size();
}

const MEDCouplingAHOGPatch *MEDCouplingAHOGMesh::getPatch(int patchId) const
{
  int sz(getNumberOfPatches());
  if(patchId<0 || patchId>=sz)
    {
      std::ostringstream oss; oss << "MEDCouplingAHOGMesh::getPatch : invalid patchId (" << patchId << ") ! Must be in [0," << sz << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return _patches[patchId];
}

MEDCouplingUMesh *MEDCouplingAHOGMesh::buildUnstructured() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> part(_mesh->buildUnstructured());
  std::vector<bool> bs(_mesh->getNumberOfCells(),false);
  std::vector<int> cgs(_mesh->getCellGridStructure());
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> > msSafe(_patches.size()+1);
  std::size_t ii(0);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingAHOGPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++,ii++)
    {
      MEDCouplingStructuredMesh::SwitchOnIdsFrom(cgs,(*it)->getBLTRRange(),bs);
      msSafe[ii+1]=(*it)->getMesh()->buildUnstructured();
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> eltsOff(DataArrayInt::BuildListOfSwitchedOff(bs));
  msSafe[0]=static_cast<MEDCouplingUMesh *>(part->buildPartOfMySelf(eltsOff->begin(),eltsOff->end(),false));
  std::vector< const MEDCouplingUMesh * > ms(msSafe.size());
  for(std::size_t i=0;i<msSafe.size();i++)
    ms[i]=msSafe[i];
  return MEDCouplingUMesh::MergeUMeshes(ms);
}

MEDCouplingAHOGMesh::MEDCouplingAHOGMesh(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                         const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop):_father(0)
{
  _mesh=MEDCouplingIMesh::New(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop);
}

MEDCouplingAHOGMesh::MEDCouplingAHOGMesh(MEDCouplingAHOGMesh *father, MEDCouplingIMesh *mesh):_father(father)
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDCouplingAHOGMesh(MEDCouplingIMesh *mesh) constructor : empty father !");
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingAHOGMesh(MEDCouplingIMesh *mesh) constructor : The input mesh is null !");
  mesh->checkCoherency();
  _mesh=mesh; _mesh->incrRef();
}

std::size_t MEDCouplingAHOGMesh::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDCouplingAHOGMesh);
}

std::vector<const BigMemoryObject *> MEDCouplingAHOGMesh::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  if((const MEDCouplingIMesh *)_mesh)
    ret.push_back((const MEDCouplingIMesh *)_mesh);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingAHOGPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      if((const MEDCouplingAHOGPatch*)*it)
        ret.push_back((const MEDCouplingAHOGPatch*)*it);
    }
  return ret;
}

void MEDCouplingAHOGMesh::updateTime() const
{
  if((const MEDCouplingIMesh *)_mesh)
    updateTimeWith(*_mesh);
}
