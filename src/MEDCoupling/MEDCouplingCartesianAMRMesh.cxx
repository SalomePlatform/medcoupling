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

#include "MEDCouplingCartesianAMRMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingUMesh.hxx"

#include <limits>
#include <sstream>
#include <numeric>

using namespace ParaMEDMEM;

/// @cond INTERNAL

int MEDCouplingCartesianAMRPatchGen::getNumberOfCellsRecursiveWithOverlap() const
{
  return _mesh->getNumberOfCellsRecursiveWithOverlap();
}

int MEDCouplingCartesianAMRPatchGen::getNumberOfCellsRecursiveWithoutOverlap() const
{
  return _mesh->getNumberOfCellsRecursiveWithoutOverlap();
}

int MEDCouplingCartesianAMRPatchGen::getMaxNumberOfLevelsRelativeToThis() const
{
  return _mesh->getMaxNumberOfLevelsRelativeToThis();
}

MEDCouplingCartesianAMRPatchGen::MEDCouplingCartesianAMRPatchGen(MEDCouplingCartesianAMRMeshGen *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatchGen constructor : input mesh is NULL !");
  _mesh=mesh; _mesh->incrRef();
}

const MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRPatchGen::getMeshSafe() const
{
  const MEDCouplingCartesianAMRMeshGen *mesh(_mesh);
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatchGen::getMeshSafe const : the mesh is NULL !");
  return mesh;
}

MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRPatchGen::getMeshSafe()
{
  MEDCouplingCartesianAMRMeshGen *mesh(_mesh);
    if(!mesh)
      throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatchGen::getMeshSafe : the mesh is NULL !");
    return mesh;
}

std::vector<const BigMemoryObject *> MEDCouplingCartesianAMRPatchGen::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  if((const MEDCouplingCartesianAMRMeshGen *)_mesh)
    ret.push_back((const MEDCouplingCartesianAMRMeshGen *)_mesh);
  return ret;
}

/*!
 * \param [in] mesh not null pointer of refined mesh replacing the cell range of \a father defined by the bottom left and top right just after.
 * \param [in] bottomLeftTopRight a vector equal to the space dimension of \a mesh that specifies for each dimension, the included cell start of the range for the first element of the pair,
 *                                a the end cell (\b excluded) of the range for the second element of the pair.
 */
MEDCouplingCartesianAMRPatch::MEDCouplingCartesianAMRPatch(MEDCouplingCartesianAMRMeshGen *mesh, const std::vector< std::pair<int,int> >& bottomLeftTopRight):MEDCouplingCartesianAMRPatchGen(mesh),_bl_tr(bottomLeftTopRight)
{
  int dim((int)bottomLeftTopRight.size()),dimExp(_mesh->getSpaceDimension());
  if(dim!=dimExp)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch constructor : space dimension of father and input bottomLeft/topRight size mismatches !");
}

void MEDCouplingCartesianAMRPatch::addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors)
{
  return getMeshSafe()->addPatch(bottomLeftTopRight,factors);
}

int MEDCouplingCartesianAMRPatch::getNumberOfOverlapedCellsForFather() const
{
  return MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt(_bl_tr);
}

/*!
 * This method states if \a other patch is in the neighborhood of \a this. The neighborhood zone is defined by \a ghostLev parameter
 * the must be >= 0. \b WARNING this method only works if \a this and \a other share the same father (no check of this will be done !).
 * Call isInMyNeighborhoodExt to deal with 2 patches not sharing directly the same father.
 *
 * \param [in] other - The other patch
 * \param [in] ghostLev - The size of the neighborhood zone.
 *
 * \throw if \a this or \a other are invalid (end before start).
 * \throw if \a ghostLev is \b not >= 0 .
 * \throw if \a this and \a other have not the same space dimension.
 *
 * \sa isInMyNeighborhoodExt
 */
bool MEDCouplingCartesianAMRPatch::isInMyNeighborhood(const MEDCouplingCartesianAMRPatch *other, int ghostLev) const
{
  if(ghostLev<0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::isInMyNeighborhood : the size of the neighborhood must be >= 0 !");
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::isInMyNeighborhood : the input patch is NULL !");
  const std::vector< std::pair<int,int> >& thisp(getBLTRRange());
  const std::vector< std::pair<int,int> >& otherp(other->getBLTRRange());
  return IsInMyNeighborhood(ghostLev,thisp,otherp);
}

/*!
 * This method states if \a other patch is in the neighborhood of \a this. The neighborhood zone is defined by \a ghostLev parameter
 * the must be >= 0. This method works even if \a this and \a other does not share the same father.
 *
 * \param [in] other - The other patch
 * \param [in] ghostLev - The size of the neighborhood zone.
 *
 * \throw if \a this or \a other are invalid (end before start).
 * \throw if \a ghostLev is \b not >= 0 .
 * \throw if \a this and \a other have not the same space dimension.
 * \throw if there is not common ancestor of \a this and \a other.
 *
 * \sa isInMyNeighborhood
 */
bool MEDCouplingCartesianAMRPatch::isInMyNeighborhoodExt(const MEDCouplingCartesianAMRPatch *other, int ghostLev) const
{
  if(ghostLev<0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::isInMyNeighborhoodExt : the size of the neighborhood must be >= 0 !");
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::isInMyNeighborhoodExt : the input patch is NULL !");
  int lev;
  const MEDCouplingCartesianAMRMeshGen *com(FindCommonAncestor(this,other,lev));//check that factors are OK
  if(lev==0)
    return isInMyNeighborhood(other,ghostLev);
  std::vector<int> offset(ComputeOffsetFromTwoToOne(com,lev,this,other));
  const std::vector< std::pair<int,int> >& thisp(getBLTRRange());
  std::vector< std::pair<int,int> > otherp(other->getBLTRRange());
  std::size_t sz(offset.size());
  for(std::size_t i=0;i<sz;i++)
    {
      otherp[i].first+=offset[i];
      otherp[i].second+=offset[i];
    }
  return IsInMyNeighborhood(ghostLev,thisp,otherp);
}

bool MEDCouplingCartesianAMRPatch::IsInMyNeighborhood(int ghostLev, const std::vector< std::pair<int,int> >& p1, const std::vector< std::pair<int,int> >& p2)
{
  std::size_t thispsize(p1.size());
  if(thispsize!=p2.size())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::isInMyNeighborhood : the dimensions must be the same !");
  for(std::size_t i=0;i<thispsize;i++)
    {
      const std::pair<int,int>& thispp(p1[i]);
      const std::pair<int,int>& otherpp(p2[i]);
      if(thispp.second<thispp.first)
        throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::isInMyNeighborhood : this patch is invalid !");
      if(otherpp.second<otherpp.first)
        throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::isInMyNeighborhood : this patch is invalid !");
      if(otherpp.first==thispp.second+ghostLev-1)
        continue;
      if(otherpp.second+ghostLev-1==thispp.first)
        continue;
      int start(std::max(thispp.first,otherpp.first)),end(std::min(thispp.second,otherpp.second));
      if(end<start)
        return false;
    }
  return true;
}

std::size_t MEDCouplingCartesianAMRPatch::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(sizeof(MEDCouplingCartesianAMRPatch));
  ret+=_bl_tr.capacity()*sizeof(std::pair<int,int>);
  return ret;
}

const MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRPatch::FindCommonAncestor(const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, int& lev)
{
  const MEDCouplingCartesianAMRMeshGen *f1(p1->_mesh),*f2(p2->_mesh);
  lev=0;
  while(f1!=f2 || f1==0 || f2==0)
    {
      f1=f1->getFather(); f2=f2->getFather();
      if(f1->getFactors()!=f2->getFactors())
        throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::FindCommonAncestor : factors differ !");
      lev++;
    }
  if(f1!=f2)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::FindCommonAncestor : no common ancestor between p1 and p2 !");
  return f1;
}

std::vector<int> MEDCouplingCartesianAMRPatch::ComputeOffsetFromTwoToOne(const MEDCouplingCartesianAMRMeshGen *comAncestor, int lev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2)
{
  if(lev<=0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::ComputeOffsetFromTwoToOne : this method is useful only for lev > 0 !");
  int dim(p1->getMesh()->getSpaceDimension());
  if(p2->getMesh()->getSpaceDimension()!=dim)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::ComputeOffsetFromTwoToOne : dimension must be the same !");
  std::vector< int > ret(dim,0);

  for(int i=0;i<lev;i++)
    {
      const MEDCouplingCartesianAMRMeshGen *f1(p1->_mesh),*f2(p2->_mesh);
      const MEDCouplingCartesianAMRPatch *p1h(0),*p2h(0);
      for(int j=0;j<lev-i;j++)
        {
          const MEDCouplingCartesianAMRMeshGen *f1tmp(f1->getFather()),*f2tmp(f2->getFather());
          int pid1(f1tmp->getPatchIdFromChildMesh(f1)),pid2(f2tmp->getPatchIdFromChildMesh(f2));
          p1h=f1tmp->getPatch(pid1); p2h=f2tmp->getPatch(pid2);
          f1=f1tmp; f2=f2tmp;
        }
      std::vector< std::pair<int,int> > p2c(p2h->getBLTRRange());
      for(int k=0;k<dim;k++)
        {
          p2c[k].first+=ret[k];
          p2c[k].second+=ret[k];
        }
      for(int k=0;k<dim;k++)
        {
          ret[k]=p2c[k].first-p1h->getBLTRRange()[k].first;
          ret[k]*=f1->getFactors()[k];
        }
    }
  return ret;
}

MEDCouplingCartesianAMRPatchGF::MEDCouplingCartesianAMRPatchGF(MEDCouplingCartesianAMRMesh *mesh):MEDCouplingCartesianAMRPatchGen(mesh)
{
}

std::size_t MEDCouplingCartesianAMRPatchGF::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDCouplingCartesianAMRPatchGF);
}

MEDCouplingDataForGodFather::MEDCouplingDataForGodFather(MEDCouplingCartesianAMRMesh *gf):_gf(gf),_tlc(gf)
{
  if(!_gf)
    throw INTERP_KERNEL::Exception("MEDCouplingDataForGodFather constructor : A data has to be attached to a AMR Mesh instance !");
  _gf->setData(this);
}

void MEDCouplingDataForGodFather::checkGodFatherFrozen() const
{
  _tlc.checkConst();
}

bool MEDCouplingDataForGodFather::changeGodFather(MEDCouplingCartesianAMRMesh *gf)
{
  bool ret(_tlc.keepTrackOfNewTL(gf));
  if(ret)
    {
      _gf=gf;
    }
  return ret;
}

/// @endcond

int MEDCouplingCartesianAMRMeshGen::getSpaceDimension() const
{
  return _mesh->getSpaceDimension();
}

void MEDCouplingCartesianAMRMeshGen::setFactors(const std::vector<int>& newFactors)
{
  if(getSpaceDimension()!=(int)newFactors.size())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::setFactors : size of input factors is not equal to the space dimension !");
  if(_factors.empty())
    {
      _factors=newFactors;
      return ;
    }
  if(_factors==newFactors)
    return ;
  if(!_patches.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::setFactors : modification of factors is not allowed when presence of patches !");
  _factors=newFactors;
  declareAsNew();
}

int MEDCouplingCartesianAMRMeshGen::getMaxNumberOfLevelsRelativeToThis() const
{
  int ret(1);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    ret=std::max(ret,(*it)->getMaxNumberOfLevelsRelativeToThis()+1);
  return ret;
}

/*!
 * This method returns the number of cells of \a this with the help of the MEDCouplingIMesh instance representing \a this.
 * The patches in \a this are ignored here.
 * \sa getNumberOfCellsAtCurrentLevelGhost, getNumberOfCellsRecursiveWithOverlap
 */
int MEDCouplingCartesianAMRMeshGen::getNumberOfCellsAtCurrentLevel() const
{
  return _mesh->getNumberOfCells();
}

/*!
 * This method returns the number of cells of \a this with the help of the MEDCouplingIMesh instance representing \a this enlarged by \a ghostLev size
 * to take into account of the ghost cells for future computation.
 * The patches in \a this are ignored here.
 *
 * \sa getNumberOfCellsAtCurrentLevel
 */
int MEDCouplingCartesianAMRMeshGen::getNumberOfCellsAtCurrentLevelGhost(int ghostLev) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingIMesh> tmp(_mesh->buildWithGhost(ghostLev));
  return tmp->getNumberOfCells();
}

/*!
 * This method returns the number of cells including the current level but \b also \b including recursively all cells of other levels
 * starting from this. The set of cells which size is returned here are generally overlapping each other.
 */
int MEDCouplingCartesianAMRMeshGen::getNumberOfCellsRecursiveWithOverlap() const
{
  int ret(_mesh->getNumberOfCells());
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      ret+=(*it)->getNumberOfCellsRecursiveWithOverlap();
    }
  return ret;
}

/*!
 * This method returns the max number of cells covering all the space without overlapping.
 * It returns the number of cells of the mesh with the highest resolution.
 * The returned value is equal to the number of cells of mesh returned by buildUnstructured.
 *
 * \sa buildUnstructured
 */
int MEDCouplingCartesianAMRMeshGen::getNumberOfCellsRecursiveWithoutOverlap() const
{
  int ret(_mesh->getNumberOfCells());
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      ret-=(*it)->getNumberOfOverlapedCellsForFather();
      ret+=(*it)->getNumberOfCellsRecursiveWithoutOverlap();
    }
  return ret;
}

const MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRMeshGen::getFather() const
{
  return _father;
}

const MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRMeshGen::getGodFather() const
{
  if(_father==0)
    return this;
  else
    return _father->getGodFather();
}

/*!
 * This method returns the level of \a this. 0 for god father. -1 for children of god father ...
 */
int MEDCouplingCartesianAMRMeshGen::getAbsoluteLevel() const
{
  if(_father==0)
    return 0;
  else
    return _father->getAbsoluteLevel()-1;
}

/*!
 * This method returns grids relative to god father to specified level \a absoluteLev.
 *
 * \return std::vector<MEDCouplingCartesianAMRPatchGen *> - objects in vector are to be managed (decrRef) by the caller.
 */
std::vector<MEDCouplingCartesianAMRPatchGen *> MEDCouplingCartesianAMRMeshGen::retrieveGridsAt(int absoluteLev) const
{
  if(absoluteLev<0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::retrieveGridsAt : absolute level must be >=0 !");
  if(_father)
    return getGodFather()->retrieveGridsAt(absoluteLev);
  //
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatchGen> > rets;
  retrieveGridsAtInternal(absoluteLev,rets);
  std::vector< MEDCouplingCartesianAMRPatchGen * > ret(rets.size());
  for(std::size_t i=0;i<rets.size();i++)
    {
      ret[i]=rets[i].retn();
    }
  return ret;
}

void MEDCouplingCartesianAMRMeshGen::detachFromFather()
{
  _father=0;
  declareAsNew();
}

/*!
 * \param [in] bottomLeftTopRight a vector equal to the space dimension of \a mesh that specifies for each dimension, the included cell start of the range for the first element of the pair,
 *                                a the end cell (\b excluded) of the range for the second element of the pair.
 * \param [in] factors The factor of refinement per axis (different from 0).
 */
void MEDCouplingCartesianAMRMeshGen::addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors)
{
  checkFactorsAndIfNotSetAssign(factors);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingIMesh> mesh(static_cast<MEDCouplingIMesh *>(_mesh->buildStructuredSubPart(bottomLeftTopRight)));
  mesh->refineWithFactor(factors);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRMeshSub> zeMesh(new MEDCouplingCartesianAMRMeshSub(this,mesh));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> elt(new MEDCouplingCartesianAMRPatch(zeMesh,bottomLeftTopRight));
  _patches.push_back(elt);
  declareAsNew();
}

/// @cond INTERNAL

class InternalPatch : public RefCountObjectOnly
{
public:
  InternalPatch():_nb_of_true(0) { }
  int getDimension() const { return (int)_part.size(); }
  double getEfficiency() const { return (double)_nb_of_true/(double)_crit.size(); }
  int getNumberOfCells() const { return (int)_crit.size(); }
  void setNumberOfTrue(int nboft) { _nb_of_true=nboft; }
  std::vector<bool>& getCriterion() { return _crit; }
  const std::vector<bool>& getConstCriterion() const { return _crit; }
  void setPart(const std::vector< std::pair<int,int> >& part) { _part=part; }
  std::vector< std::pair<int,int> >& getPart() { return _part; }
  const std::vector< std::pair<int,int> >& getConstPart() const { return _part; }
  bool presenceOfTrue() const { return _nb_of_true>0; }
  std::vector<int> computeCGS() const { return MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(_part); }
  std::vector< std::vector<int> > computeSignature() const { return MEDCouplingStructuredMesh::ComputeSignaturePerAxisOf(computeCGS(),getConstCriterion()); }
  double getEfficiencyPerAxis(int axisId) const { return (double)_nb_of_true/((double)(_part[axisId].second-_part[axisId].first)); }
  void zipToFitOnCriterion();
  void updateNumberOfTrue() const;
  MEDCouplingAutoRefCountObjectPtr<InternalPatch> extractPart(const std::vector< std::pair<int,int> >&partInGlobal) const;
  MEDCouplingAutoRefCountObjectPtr<InternalPatch> deepCpy() const;
protected:
  ~InternalPatch() { }
private:
  mutable int _nb_of_true;
  std::vector<bool> _crit;
  //! _part is global
  std::vector< std::pair<int,int> > _part;
};

void InternalPatch::zipToFitOnCriterion()
{
  std::vector<int> cgs(computeCGS());
  std::vector<bool> newCrit;
  std::vector< std::pair<int,int> > newPart,newPart2;
  int newNbOfTrue(MEDCouplingStructuredMesh::FindMinimalPartOf(cgs,_crit,newCrit,newPart));
  MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt(_part,newPart,newPart2);
  if(newNbOfTrue!=_nb_of_true)
    throw INTERP_KERNEL::Exception("InternalPatch::zipToFitOnCrit : internal error !");
  _crit=newCrit; _part=newPart2;
}

void InternalPatch::updateNumberOfTrue() const
{
  _nb_of_true=(int)std::count(_crit.begin(),_crit.end(),true);
}

MEDCouplingAutoRefCountObjectPtr<InternalPatch> InternalPatch::extractPart(const std::vector< std::pair<int,int> >&partInGlobal) const
{
  MEDCouplingAutoRefCountObjectPtr<InternalPatch> ret(new InternalPatch);
  std::vector<int> cgs(computeCGS());
  std::vector< std::pair<int,int> > newPart;
  MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt(_part,partInGlobal,newPart);
  MEDCouplingStructuredMesh::ExtractFieldOfBoolFrom(cgs,_crit,newPart,ret->getCriterion());
  ret->setPart(partInGlobal);
  ret->updateNumberOfTrue();
  return ret;
}

MEDCouplingAutoRefCountObjectPtr<InternalPatch> InternalPatch::deepCpy() const
{
  MEDCouplingAutoRefCountObjectPtr<InternalPatch> ret(new InternalPatch);
  (*ret)=*this;
  return ret;
}

void DissectBigPatch(const INTERP_KERNEL::BoxSplittingOptions& bso, const InternalPatch *patchToBeSplit, int axisId, int rangeOfAxisId, bool& cutFound, int& cutPlace)
{
  cutFound=false; cutPlace=-1;
  std::vector<double> ratio(rangeOfAxisId-1);
  for(int id=0;id<rangeOfAxisId-1;id++)
    {
      double efficiency[2];
      for(int h=0;h<2;h++)
        {
          std::vector< std::pair<int,int> > rectH(patchToBeSplit->getConstPart());
          if(h==0)
            rectH[axisId].second=patchToBeSplit->getConstPart()[axisId].first+id;
          else
            rectH[axisId].first=patchToBeSplit->getConstPart()[axisId].first+id;
          MEDCouplingAutoRefCountObjectPtr<InternalPatch> p(patchToBeSplit->deepCpy());
          p->zipToFitOnCriterion();
          //anouar rectH ?
          efficiency[h]=p->getEfficiencyPerAxis(axisId);
        }
      ratio[id]=std::max(efficiency[0],efficiency[1])/std::min(efficiency[0],efficiency[1]);
    }
  int minCellDirection(bso.getMinCellDirection()),indexMin(-1);
  int dimRatioInner(rangeOfAxisId-1-2*(minCellDirection-1));
  std::vector<double> ratio_inner(dimRatioInner);
  double minRatio(1.e10);
  for(int i=0; i<dimRatioInner; i++)
    {
      if(ratio[minCellDirection-1+i]<minRatio)
        {
          minRatio=ratio[minCellDirection-1+i];
          indexMin=i+minCellDirection;
        }
    }
  cutFound=true; cutPlace=indexMin+patchToBeSplit->getConstPart()[axisId].first-1;
}

void FindHole(const INTERP_KERNEL::BoxSplittingOptions& bso, const InternalPatch *patchToBeSplit, int& axisId, bool& cutFound, int& cutPlace)
{
  cutPlace=-1; cutFound=false;
  int minCellDirection(bso.getMinCellDirection());
  const int dim(patchToBeSplit->getDimension());
  std::vector< std::vector<int> > signatures(patchToBeSplit->computeSignature());
  for(int id=0;id<dim;id++)
    {
      const std::vector<int>& signature(signatures[id]);
      std::vector<int> hole;
      std::vector<double> distance;
      int len((int)signature.size());
      for(int i=0;i<len;i++)
        if(signature[i]==0)
          if(len>= 2*minCellDirection && i >= minCellDirection-1 && i <= len-minCellDirection-1)
            hole.push_back(i);
      if(!hole.empty())
        {
          double center(((double)len/2.));
          for(std::size_t i=0;i<hole.size();i++)
            distance.push_back(fabs(hole[i]+1.-center));

          std::size_t posDistanceMin(std::distance(distance.begin(),std::min_element(distance.begin(),distance.end())));
          cutFound=true;
          axisId=id;
          cutPlace=hole[posDistanceMin]+patchToBeSplit->getConstPart()[axisId].first+1;
          return ;
        }
    }
}

void FindInflection(const INTERP_KERNEL::BoxSplittingOptions& bso, const InternalPatch *patchToBeSplit, bool& cutFound, int& cutPlace, int& axisId)
{
  cutFound=false; cutPlace=-1;// do not set axisId before to be sure that cutFound was set to true

  const std::vector< std::pair<int,int> >& part(patchToBeSplit->getConstPart());
  int sign,minCellDirection(bso.getMinCellDirection());
  const int dim(patchToBeSplit->getDimension());

  std::vector<int> zeroCrossDims(dim,-1);
  std::vector<int> zeroCrossVals(dim,-1);
  std::vector< std::vector<int> > signatures(patchToBeSplit->computeSignature());
  for (int id=0;id<dim;id++)
    {
      const std::vector<int>& signature(signatures[id]);

      std::vector<int> derivate_second_order,gradient_absolute,signe_change,zero_cross,edge,max_cross_list ;
      std::vector<double> distance ;

      for (unsigned int i=1;i<signature.size()-1;i++)
        derivate_second_order.push_back(signature[i-1]-2*signature[i]+signature[i+1]) ;

      // Gradient absolute value
      for ( unsigned int i=1;i<derivate_second_order.size();i++)
        gradient_absolute.push_back(fabs(derivate_second_order[i]-derivate_second_order[i-1])) ;
      if(derivate_second_order.empty())
        continue;
      for (unsigned int i=0;i<derivate_second_order.size()-1;i++)
        {
          if (derivate_second_order[i]*derivate_second_order[i+1] < 0 )
            sign = -1 ;
          if (derivate_second_order[i]*derivate_second_order[i+1] > 0 )
            sign = 1 ;
          if (derivate_second_order[i]*derivate_second_order[i+1] == 0 )
            sign = 0 ;
          if ( sign==0 || sign==-1 )
            if ( i >= (unsigned int)minCellDirection-2 && i <= signature.size()-minCellDirection-2 )
              {
                zero_cross.push_back(i) ;
                edge.push_back(gradient_absolute[i]) ;
              }
          signe_change.push_back(sign) ;
        }
      if ( zero_cross.size() > 0 )
        {
          int max_cross=*max_element(edge.begin(),edge.end()) ;
          for (unsigned int i=0;i<edge.size();i++)
            if (edge[i]==max_cross)
              max_cross_list.push_back(zero_cross[i]+1) ;

          double center((signature.size()/2.0));
          for (unsigned int i=0;i<max_cross_list.size();i++)
            distance.push_back(fabs(max_cross_list[i]+1-center));

          float distance_min=*min_element(distance.begin(),distance.end()) ;
          int pos_distance_min=find(distance.begin(),distance.end(),distance_min)-distance.begin();
          int best_place = max_cross_list[pos_distance_min] + part[id].first ;
          if ( max_cross >=0 )
            {
              zeroCrossDims[id] = best_place ;
              zeroCrossVals[id] = max_cross ;
            }
        }
      derivate_second_order.clear() ;
      gradient_absolute.clear() ;
      signe_change.clear() ;
      zero_cross.clear() ;
      edge.clear() ;
      max_cross_list.clear() ;
      distance.clear() ;
    }

  if ( zeroCrossDims[0]!=-1 || zeroCrossDims[1]!=-1  )
    {
      int max_cross_dims = *max_element(zeroCrossVals.begin(),zeroCrossVals.end()) ;

      if (zeroCrossVals[0]==max_cross_dims &&  zeroCrossVals[1]==max_cross_dims )
        {
          int nl_left(part[0].second-part[0].first);
          int nc_left(part[1].second-part[1].first);
          if ( nl_left >=  nc_left )
            max_cross_dims = 0 ;
          else
            max_cross_dims = 1 ;
        }
      else
        max_cross_dims=std::find(zeroCrossVals.begin(),zeroCrossVals.end(),max_cross_dims)-zeroCrossVals.begin();
      cutFound=true;
      cutPlace=zeroCrossDims[max_cross_dims];
      axisId=max_cross_dims ;
    }
}

void TryAction4(const INTERP_KERNEL::BoxSplittingOptions& bso, const InternalPatch *patchToBeSplit, int axisId, int rangeOfAxisId, bool& cutFound, int& cutPlace)
{
  cutFound=false;
  if(patchToBeSplit->getEfficiency()<=bso.getEffeciencySnd())
    {
      if(rangeOfAxisId>=2*bso.getMinCellDirection())
        {
          cutFound=true;
          cutPlace=rangeOfAxisId/2+patchToBeSplit->getConstPart()[axisId].first-1;
        }
    }
  else
    {
      if(patchToBeSplit->getNumberOfCells()>bso.getMaxCells())
        {
          DissectBigPatch(bso,patchToBeSplit,axisId,rangeOfAxisId,cutFound,cutPlace);
        }
    }
}

MEDCouplingAutoRefCountObjectPtr<InternalPatch> DealWithNoCut(const InternalPatch *patch)
{
  MEDCouplingAutoRefCountObjectPtr<InternalPatch> ret(const_cast<InternalPatch *>(patch));
  ret->incrRef();
  return ret;
}

void DealWithCut(const InternalPatch *patchToBeSplit, int axisId, int cutPlace, std::vector<MEDCouplingAutoRefCountObjectPtr<InternalPatch> >& listOfPatches)
{
  MEDCouplingAutoRefCountObjectPtr<InternalPatch> leftPart,rightPart;
  std::vector< std::pair<int,int> > rect(patchToBeSplit->getConstPart());
  std::vector< std::pair<int,int> > leftRect(rect),rightRect(rect);
  leftRect[axisId].second=cutPlace+1;
  rightRect[axisId].first=cutPlace+1;
  leftPart=patchToBeSplit->extractPart(leftRect);
  rightPart=patchToBeSplit->extractPart(rightRect);
  leftPart->zipToFitOnCriterion(); rightPart->zipToFitOnCriterion();
  listOfPatches.push_back(leftPart);
  listOfPatches.push_back(rightPart);
}

/// @endcond

/*!
 * This method creates patches in \a this (by destroying the patches if any). This method uses \a criterion array as a field on cells on this level.
 */
void MEDCouplingCartesianAMRMeshGen::createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const std::vector<bool>& criterion, const std::vector<int>& factors)
{
  int nbCells(getNumberOfCellsAtCurrentLevel());
  if(nbCells!=(int)criterion.size())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createPatchesFromCriterion : the number of tuples of criterion array must be equal to the number of cells at the current level !");
  _patches.clear();
  std::vector<int> cgs(_mesh->getCellGridStructure());
  std::vector< MEDCouplingAutoRefCountObjectPtr<InternalPatch> > listOfPatches,listOfPatchesOK;
  //
  MEDCouplingAutoRefCountObjectPtr<InternalPatch> p(new InternalPatch);
  p->setNumberOfTrue(MEDCouplingStructuredMesh::FindMinimalPartOf(cgs,criterion,p->getCriterion(),p->getPart()));
  if(p->presenceOfTrue())
    listOfPatches.push_back(p);
  while(!listOfPatches.empty())
    {
      std::vector< MEDCouplingAutoRefCountObjectPtr<InternalPatch> > listOfPatchesTmp;
      for(std::vector< MEDCouplingAutoRefCountObjectPtr<InternalPatch> >::iterator it=listOfPatches.begin();it!=listOfPatches.end();it++)
        {
          //
          int axisId,rangeOfAxisId,cutPlace;
          bool cutFound;
          MEDCouplingStructuredMesh::FindTheWidestAxisOfGivenRangeInCompactFrmt((*it)->getConstPart(),axisId,rangeOfAxisId);
          if((*it)->getEfficiency()>=bso.getEffeciency() && (*it)->getNumberOfCells()<bso.getMaxCells())
            { listOfPatchesOK.push_back(DealWithNoCut(*it)); continue; }//action 1
          FindHole(bso,*it,axisId,cutFound,cutPlace);//axisId overwritten here if FindHole equal to true !
          if(cutFound)
            { DealWithCut(*it,axisId,cutPlace,listOfPatchesTmp); continue; }//action 2
          FindInflection(bso,*it,cutFound,cutPlace,axisId);//axisId overwritten here if cutFound equal to true !
          if(cutFound)
            { DealWithCut(*it,axisId,cutPlace,listOfPatchesTmp); continue; }//action 3
          TryAction4(bso,*it,axisId,rangeOfAxisId,cutFound,cutPlace);
          if(cutFound)
            { DealWithCut(*it,axisId,cutPlace,listOfPatchesTmp); continue; }//action 4
          listOfPatchesOK.push_back(DealWithNoCut(*it));
        }
      listOfPatches=listOfPatchesTmp;
    }
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<InternalPatch> >::const_iterator it=listOfPatchesOK.begin();it!=listOfPatchesOK.end();it++)
    addPatch((*it)->getConstPart(),factors);
  declareAsNew();
}

/*!
 * This method creates patches in \a this (by destroying the patches if any). This method uses \a criterion array as a field on cells on this level.
 */
void MEDCouplingCartesianAMRMeshGen::createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayByte *criterion, const std::vector<int>& factors)
{
  if(!criterion || !criterion->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createPatchesFromCriterion : the criterion DataArrayByte instance must be allocated and not NULL !");
  std::vector<bool> crit(criterion->toVectorOfBool());//check that criterion has one component.
  createPatchesFromCriterion(bso,crit,factors);
  declareAsNew();
}

void MEDCouplingCartesianAMRMeshGen::removeAllPatches()
{
  _patches.clear();
  declareAsNew();
}

void MEDCouplingCartesianAMRMeshGen::removePatch(int patchId)
{
  checkPatchId(patchId);
  int sz((int)_patches.size()),j(0);
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> > patches(sz-1);
  for(int i=0;i<sz;i++)
    if(i!=patchId)
      patches[j++]=_patches[i];
  (const_cast<MEDCouplingCartesianAMRMeshGen *>(_patches[patchId]->getMesh()))->detachFromFather();
  _patches=patches;
  declareAsNew();
}

int MEDCouplingCartesianAMRMeshGen::getNumberOfPatches() const
{
  return (int)_patches.size();
}

int MEDCouplingCartesianAMRMeshGen::getPatchIdFromChildMesh(const MEDCouplingCartesianAMRMeshGen *mesh) const
{
  int ret(0);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++,ret++)
    {
      if((*it)->getMesh()==mesh)
        return ret;
    }
  throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::getPatchIdFromChildMesh : no such a mesh in my direct progeny !");
}

std::vector< const MEDCouplingCartesianAMRPatch *> MEDCouplingCartesianAMRMeshGen::getPatches() const
{
  std::size_t sz(_patches.size());
  std::vector< const MEDCouplingCartesianAMRPatch *> ret(sz);
  for(std::size_t i=0;i<sz;i++)
    ret[i]=_patches[i];
  return ret;
}

const MEDCouplingCartesianAMRPatch *MEDCouplingCartesianAMRMeshGen::getPatch(int patchId) const
{
  checkPatchId(patchId);
  return _patches[patchId];
}

/*!
 * This method states if patch2 (with id \a patchId2) is in the neighborhood of patch1 (with id \a patchId1).
 * The neighborhood size is defined by \a ghostLev in the reference of \a this ( \b not in the reference of patches !).
 */
bool MEDCouplingCartesianAMRMeshGen::isPatchInNeighborhoodOf(int patchId1, int patchId2, int ghostLev) const
{
  const MEDCouplingCartesianAMRPatch *p1(getPatch(patchId1)),*p2(getPatch(patchId2));
  return p1->isInMyNeighborhood(p2,GetGhostLevelInFineRef(ghostLev,_factors));
}

/*!
 * This method creates a new cell field array on given \a patchId patch in \a this starting from a coarse cell field on \a this \a cellFieldOnThis.
 * This method can be seen as a fast projection from the cell field \a cellFieldOnThis on \c this->getImageMesh() to a refined part of \a this
 * defined by the patch with id \a patchId.
 *
 * \param [in] patchId - The id of the patch \a cellFieldOnThis has to be put on.
 * \param [in] cellFieldOnThis - The array of the cell field on \c this->getImageMesh() to be projected to patch having id \a patchId.
 * \return DataArrayDouble * - The array of the cell field on the requested patch
 *
 * \throw if \a patchId is not in [ 0 , \c this->getNumberOfPatches() )
 * \throw if \a cellFieldOnThis is NULL or not allocated
 * \sa fillCellFieldOnPatch, MEDCouplingIMesh::SpreadCoarseToFine
 */
DataArrayDouble *MEDCouplingCartesianAMRMeshGen::createCellFieldOnPatch(int patchId, const DataArrayDouble *cellFieldOnThis) const
{
  if(!cellFieldOnThis || !cellFieldOnThis->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createCellFieldOnPatch : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  const MEDCouplingIMesh *fine(patch->getMesh()->getImageMesh());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(fine->getNumberOfCells(),cellFieldOnThis->getNumberOfComponents());
  ret->copyStringInfoFrom(*cellFieldOnThis);
  MEDCouplingIMesh::SpreadCoarseToFine(cellFieldOnThis,_mesh->getCellGridStructure(),ret,patch->getBLTRRange(),getFactors());
  return ret.retn();
}

/*!
 * This method is equivalent to MEDCouplingCartesianAMRMesh::createCellFieldOnPatch except that here instead of creating a new instance
 * it fills the value into the \a cellFieldOnPatch data.
 *
 * \param [in] patchId - The id of the patch \a cellFieldOnThis has to be put on.
 * \param [in] cellFieldOnThis - The array of the cell field on \c this->getImageMesh() to be projected to patch having id \a patchId.
 * \param [in,out] cellFieldOnPatch - The array of the cell field on the requested patch to be filled.
 *
 * \sa createCellFieldOnPatch, fillCellFieldComingFromPatch
 */
void MEDCouplingCartesianAMRMeshGen::fillCellFieldOnPatch(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch) const
{
  if(!cellFieldOnThis || !cellFieldOnThis->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createCellFieldOnPatch : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  MEDCouplingIMesh::SpreadCoarseToFine(cellFieldOnThis,_mesh->getCellGridStructure(),cellFieldOnPatch,patch->getBLTRRange(),getFactors());
}

/*!
 * This method is the generalization of fillCellFieldOnPatch method. This method only projects coarse to fine without considering the
 * potential neighbor patches covered by the ghost cells of patch with id \a patchId.
 *
 * \param [in] patchId - The id of the patch \a cellFieldOnThis has to be put on.
 * \param [in] cellFieldOnThis - The array of the cell field on \c this->getImageMesh() to be projected to patch having id \a patchId.
 * \param [in,out] cellFieldOnPatch - The array of the cell field on the requested patch to be filled.
 * \param [in] ghostLev - The size of the ghost zone (must be >=0 !)
 *
 * \sa fillCellFieldOnPatch, fillCellFieldOnPatchGhostAdv
 */
void MEDCouplingCartesianAMRMeshGen::fillCellFieldOnPatchGhost(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, int ghostLev) const
{
  if(!cellFieldOnThis || !cellFieldOnThis->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createCellFieldOnPatchGhost : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  MEDCouplingIMesh::SpreadCoarseToFineGhost(cellFieldOnThis,_mesh->getCellGridStructure(),cellFieldOnPatch,patch->getBLTRRange(),getFactors(),ghostLev);
}

/*!
 * This method is equivalent to  fillCellFieldOnPatchGhost except that here \b ONLY \b the \b ghost \b zone will be updated
 * in \a cellFieldOnPatch.
 *
 * \param [in] patchId - The id of the patch \a cellFieldOnThis has to be put on.
 * \param [in] cellFieldOnThis - The array of the cell field on \c this->getImageMesh() to be projected to patch having id \a patchId.
 * \param [in,out] cellFieldOnPatch - The array of the cell field on the requested patch to be filled \b only \b in \b the \b ghost \b zone.
 * \param [in] ghostLev - The size of the ghost zone (must be >=0 !)
 */
void MEDCouplingCartesianAMRMeshGen::fillCellFieldOnPatchOnlyOnGhostZone(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, int ghostLev) const
{
  if(!cellFieldOnThis || !cellFieldOnThis->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::fillCellFieldOnPatchOnlyOnGhostZone : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  MEDCouplingIMesh::SpreadCoarseToFineGhostZone(cellFieldOnThis,_mesh->getCellGridStructure(),cellFieldOnPatch,patch->getBLTRRange(),getFactors(),ghostLev);
}

/*!
 * This method is a refinement of fillCellFieldOnPatchGhost. fillCellFieldOnPatchGhost is first called.
 * Then for all other patches than those pointed by \a patchId that overlap the ghost zone of the patch impact the ghost zone adequately.
 *
 * \param [in] patchId - The id of the patch \a cellFieldOnThis has to be put on.
 * \param [in] cellFieldOnThis - The array of the cell field on \c this->getImageMesh() to be projected to patch having id \a patchId.
 * \param [in,out] cellFieldOnPatch - The array of the cell field on the requested patch to be filled.
 * \param [in] ghostLev - The size of the ghost zone (must be >=0 !)
 * \param [in] arrsOnPatches - \b WARNING arrsOnPatches[patchId] is \b NOT \b const. All others are const.
 *
 * \sa fillCellFieldOnPatchOnlyGhostAdv
 */
void MEDCouplingCartesianAMRMeshGen::fillCellFieldOnPatchGhostAdv(int patchId, const DataArrayDouble *cellFieldOnThis, int ghostLev, const std::vector<const DataArrayDouble *>& arrsOnPatches) const
{
  int nbp(getNumberOfPatches());
  if(nbp!=(int)arrsOnPatches.size())
    {
      std::ostringstream oss; oss << "MEDCouplingCartesianAMRMesh::fillCellFieldOnPatchGhostAdv : there are " << nbp << " patches in this and " << arrsOnPatches.size() << " arrays in the last parameter !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  DataArrayDouble *theFieldToFill(const_cast<DataArrayDouble *>(arrsOnPatches[patchId]));
  // first, do as usual
  fillCellFieldOnPatchGhost(patchId,cellFieldOnThis,theFieldToFill,ghostLev);
  fillCellFieldOnPatchOnlyGhostAdv(patchId,ghostLev,arrsOnPatches);
}

/*!
 * This method updates the patch with id \a patchId considering the only the all the patches in \a this to fill ghost zone.
 * So \b warning, the DataArrayDouble instance \a arrsOnPatches[patchId] is non const.
 *
 * \sa getPatchIdsInTheNeighborhoodOf
 */
void MEDCouplingCartesianAMRMeshGen::fillCellFieldOnPatchOnlyGhostAdv(int patchId, int ghostLev, const std::vector<const DataArrayDouble *>& arrsOnPatches) const
{
  int nbp(getNumberOfPatches()),dim(getSpaceDimension());
  if(nbp!=(int)arrsOnPatches.size())
    {
      std::ostringstream oss; oss << "MEDCouplingCartesianAMRMesh::fillCellFieldOnPatchOnlyGhostAdv : there are " << nbp << " patches in this and " << arrsOnPatches.size() << " arrays in the last parameter !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  DataArrayDouble *theFieldToFill(const_cast<DataArrayDouble *>(arrsOnPatches[patchId]));
  const MEDCouplingCartesianAMRPatch *refP(getPatch(patchId));
  const std::vector< std::pair<int,int> >& refBLTR(refP->getBLTRRange());//[(1,4),(2,4)]
  std::vector<int> dimsCoarse(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(refBLTR));//[3,2]
  std::transform(dimsCoarse.begin(),dimsCoarse.end(),_factors.begin(),dimsCoarse.begin(),std::multiplies<int>());//[12,8]
  std::transform(dimsCoarse.begin(),dimsCoarse.end(),dimsCoarse.begin(),std::bind2nd(std::plus<int>(),2*ghostLev));//[14,10]
  std::vector< std::pair<int,int> > rangeCoarse(MEDCouplingStructuredMesh::GetCompactFrmtFromDimensions(dimsCoarse));//[(0,14),(0,10)]
  std::vector<int> fakeFactors(dim,1),ids(getPatchIdsInTheNeighborhoodOf(patchId,ghostLev));
  //
  for(std::vector<int>::const_iterator it=ids.begin();it!=ids.end();it++)
    {
      const MEDCouplingCartesianAMRPatch *otherP(getPatch(*it));
      const std::vector< std::pair<int,int> >& otherBLTR(otherP->getBLTRRange());//[(4,5),(3,4)]
      std::vector< std::pair<int,int> > tmp0,tmp1,tmp2;
      MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt(refBLTR,otherBLTR,tmp0,false);//tmp0=[(3,4),(1,2)]
      ApplyFactorsOnCompactFrmt(tmp0,_factors);//tmp0=[(12,16),(4,8)]
      ApplyGhostOnCompactFrmt(tmp0,ghostLev);//tmp0=[(13,17),(5,9)]
      std::vector< std::pair<int,int> > interstRange(MEDCouplingStructuredMesh::IntersectRanges(tmp0,rangeCoarse));//interstRange=[(13,14),(5,9)]
      MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt(otherBLTR,refBLTR,tmp1,false);//tmp1=[(-3,0),(-1,1)]
      ApplyFactorsOnCompactFrmt(tmp1,_factors);//tmp1=[(-12,-4),(-4,0)]
      MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt(tmp1,interstRange,tmp2,false);//tmp2=[(1,2),(1,5)]
      //
      std::vector< std::pair<int,int> > dimsFine(otherBLTR);
      ApplyFactorsOnCompactFrmt(dimsFine,_factors);
      ApplyAllGhostOnCompactFrmt(dimsFine,ghostLev);
      //
      MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ghostVals(MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(dimsFine),arrsOnPatches[*it],tmp2));
      MEDCouplingIMesh::CondenseFineToCoarse(dimsCoarse,ghostVals,interstRange,fakeFactors,theFieldToFill);
    }
}

/*!
 * This method updates \a cellFieldOnThis part of values coming from the cell field \a cellFieldOnPatch lying on patch having id \a patchId.
 *
 * \param [in] patchId - The id of the patch \a cellFieldOnThis has to be put on.
 * \param [in] cellFieldOnPatch - The array of the cell field on patch with id \a patchId.
 * \param [in,out] cellFieldOnThis The array of the cell field on \a this to be updated only on the part concerning the patch with id \a patchId.
 *
 * \throw if \a patchId is not in [ 0 , \c this->getNumberOfPatches() )
 * \throw if \a cellFieldOnPatch is NULL or not allocated
 * \sa createCellFieldOnPatch, MEDCouplingIMesh::CondenseFineToCoarse,fillCellFieldComingFromPatchGhost
 */
void MEDCouplingCartesianAMRMeshGen::fillCellFieldComingFromPatch(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis) const
{
  if(!cellFieldOnPatch || !cellFieldOnPatch->isAllocated())
      throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::fillCellFieldComingFromPatch : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  MEDCouplingIMesh::CondenseFineToCoarse(_mesh->getCellGridStructure(),cellFieldOnPatch,patch->getBLTRRange(),getFactors(),cellFieldOnThis);
}

/*!
 * This method is the extension of MEDCouplingCartesianAMRMesh::fillCellFieldComingFromPatch managing the ghost cells. If this
 * method is called with \a ghostLev equal to 0 it behaves exactly as MEDCouplingCartesianAMRMesh::fillCellFieldComingFromPatch.
 *
 * \param [in] patchId - The id of the patch \a cellFieldOnThis has to be put on.
 * \param [in] cellFieldOnPatch - The array of the cell field on patch with id \a patchId.
 * \param [in,out] cellFieldOnThis The array of the cell field on \a this to be updated only on the part concerning the patch with id \a patchId.
 * \param [in] ghostLev The size of ghost zone (must be >= 0 !)
 *
 * \throw if \a patchId is not in [ 0 , \c this->getNumberOfPatches() )
 * \throw if \a cellFieldOnPatch is NULL or not allocated
 * \sa fillCellFieldComingFromPatch
 */
void MEDCouplingCartesianAMRMeshGen::fillCellFieldComingFromPatchGhost(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, int ghostLev) const
{
  if(!cellFieldOnPatch || !cellFieldOnPatch->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::fillCellFieldComingFromPatchGhost : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  MEDCouplingIMesh::CondenseFineToCoarseGhost(_mesh->getCellGridStructure(),cellFieldOnPatch,patch->getBLTRRange(),getFactors(),cellFieldOnThis,ghostLev);
}

/*!
 * This method finds all patches (located by their ids) that are in the neighborhood of patch with id \a patchId. The neighborhood size is
 * defined by ghostLev.
 *
 * \param [in] patchId - the id of the considered patch.
 * \param [in] ghostLev - the size of the neighborhood.
 * \return DataArrayInt * - the newly allocated array containing the list of patches in the neighborhood of the considered patch. This array is to be deallocated by the caller.
 */
DataArrayInt *MEDCouplingCartesianAMRMeshGen::findPatchesInTheNeighborhoodOf(int patchId, int ghostLev) const
{
  int nbp(getNumberOfPatches());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  for(int i=0;i<nbp;i++)
    {
      if(i!=patchId)
        if(isPatchInNeighborhoodOf(i,patchId,ghostLev))
          ret->pushBackSilent(i);
    }
  return ret.retn();
}

MEDCouplingUMesh *MEDCouplingCartesianAMRMeshGen::buildUnstructured() const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> part(_mesh->buildUnstructured());
  std::vector<bool> bs(_mesh->getNumberOfCells(),false);
  std::vector<int> cgs(_mesh->getCellGridStructure());
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> > msSafe(_patches.size()+1);
  std::size_t ii(0);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++,ii++)
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

/*!
 * This method returns a mesh containing as cells that there is patches at the current level.
 * The patches are seen like 'boxes' that is too say the refinement will not appear here.
 *
 * \return MEDCoupling1SGTUMesh * - A new object to be managed by the caller containing as cells as there are patches in \a this.
 */
MEDCoupling1SGTUMesh *MEDCouplingCartesianAMRMeshGen::buildMeshFromPatchEnvelop() const
{
  std::vector<const MEDCoupling1SGTUMesh *> cells;
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCoupling1SGTUMesh> > cellsSafe;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      const MEDCouplingCartesianAMRPatch *patch(*it);
      if(patch)
        {
          MEDCouplingAutoRefCountObjectPtr<MEDCouplingIMesh> cell(patch->getMesh()->getImageMesh()->asSingleCell());
          MEDCouplingAutoRefCountObjectPtr<MEDCoupling1SGTUMesh> cell1SGT(cell->build1SGTUnstructured());
          cellsSafe.push_back(cell1SGT); cells.push_back(cell1SGT);
        }
    }
  return MEDCoupling1SGTUMesh::Merge1SGTUMeshes(cells);
}

MEDCoupling1SGTUMesh *MEDCouplingCartesianAMRMeshGen::buildMeshOfDirectChildrenOnly() const
{
  std::vector<const MEDCoupling1SGTUMesh *> patches;
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCoupling1SGTUMesh> > patchesSafe;
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
      {
        const MEDCouplingCartesianAMRPatch *patch(*it);
        if(patch)
          {
            MEDCouplingAutoRefCountObjectPtr<MEDCoupling1SGTUMesh> patchMesh(patch->getMesh()->getImageMesh()->build1SGTUnstructured());
            patchesSafe.push_back(patchMesh); patches.push_back(patchMesh);
          }
      }
    return MEDCoupling1SGTUMesh::Merge1SGTUMeshes(patches);
}

/*!
 * This method works same as buildUnstructured except that arrays are given in input to build a field on cell in output.
 * \return MEDCouplingFieldDouble * - a newly created instance the caller has reponsability to deal with.
 * \sa buildUnstructured
 */
MEDCouplingFieldDouble *MEDCouplingCartesianAMRMeshGen::buildCellFieldOnRecurseWithoutOverlapWithoutGhost(int ghostSz, const std::vector<const DataArrayDouble *>& recurseArrs) const
{
  if(recurseArrs.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::buildCellFieldOnRecurseWithoutOverlapWithoutGhost : array is empty ! Should never happen !");
  //
  std::vector<bool> bs(_mesh->getNumberOfCells(),false);
  std::vector<int> cgs(_mesh->getCellGridStructure());
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> > msSafe(_patches.size()+1);
  std::size_t ii(0);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++,ii++)
    {
      MEDCouplingStructuredMesh::SwitchOnIdsFrom(cgs,(*it)->getBLTRRange(),bs);
      std::vector<const DataArrayDouble *> tmpArrs(extractSubTreeFromGlobalFlatten((*it)->getMesh(),recurseArrs));
      msSafe[ii+1]=(*it)->getMesh()->buildCellFieldOnRecurseWithoutOverlapWithoutGhost(ghostSz,tmpArrs);
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> eltsOff(DataArrayInt::BuildListOfSwitchedOff(bs));
  //
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_CELLS));
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2(extractGhostFrom(ghostSz,recurseArrs[0]));
  arr2=arr2->selectByTupleIdSafe(eltsOff->begin(),eltsOff->end());
  ret->setArray(arr2);
  ret->setName(arr2->getName());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> part(_mesh->buildUnstructured());
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingMesh> mesh(part->buildPartOfMySelf(eltsOff->begin(),eltsOff->end(),false));
  ret->setMesh(mesh);
  msSafe[0]=ret;
  //
  std::vector< const MEDCouplingFieldDouble * > ms(msSafe.size());
  for(std::size_t i=0;i<msSafe.size();i++)
    ms[i]=msSafe[i];
  //
  return MEDCouplingFieldDouble::MergeFields(ms);
}

/*!
 * This method extracts from \arr arr the part inside \a arr by cutting the \a ghostSz external part.
 * \arr is expected to be an array having a number of tuples equal to \c getImageMesh()->buildWithGhost(ghostSz).
 */
DataArrayDouble *MEDCouplingCartesianAMRMeshGen::extractGhostFrom(int ghostSz, const DataArrayDouble *arr) const
{
  std::vector<int> st(_mesh->getCellGridStructure());
  std::vector< std::pair<int,int> > p(MEDCouplingStructuredMesh::GetCompactFrmtFromDimensions(st));
  std::transform(st.begin(),st.end(),st.begin(),std::bind2nd(std::plus<int>(),2*ghostSz));
  ApplyGhostOnCompactFrmt(p,ghostSz);
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret(MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom(st,arr,p));
  return ret.retn();
}

/*!
 * This method returns all the patches in \a this not equal to \a patchId that are in neighborhood of patch with id \a patchId.
 *
 * \sa fillCellFieldOnPatchOnlyGhostAdv
 */
std::vector<int> MEDCouplingCartesianAMRMeshGen::getPatchIdsInTheNeighborhoodOf(int patchId, int ghostLev) const
{
  std::vector<int> ret;
  int nbp(getNumberOfPatches());
  //
  for(int i=0;i<nbp;i++)
    {
      if(i!=patchId)
        if(isPatchInNeighborhoodOf(i,patchId,ghostLev))
          ret.push_back(i);
    }
  return ret;
}

MEDCouplingCartesianAMRMeshGen::MEDCouplingCartesianAMRMeshGen(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                               const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop):_father(0)
{
  _mesh=MEDCouplingIMesh::New(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop);
}

MEDCouplingCartesianAMRMeshGen::MEDCouplingCartesianAMRMeshGen(MEDCouplingCartesianAMRMeshGen *father, MEDCouplingIMesh *mesh):_father(father)
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen(MEDCouplingIMesh *mesh) constructor : empty father !");
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen(MEDCouplingIMesh *mesh) constructor : The input mesh is null !");
  mesh->checkCoherency();
  _mesh=mesh; _mesh->incrRef();
}

void MEDCouplingCartesianAMRMeshGen::checkPatchId(int patchId) const
{
  int sz(getNumberOfPatches());
  if(patchId<0 || patchId>=sz)
    {
      std::ostringstream oss; oss << "MEDCouplingCartesianAMRMeshGen::checkPatchId : invalid patchId (" << patchId << ") ! Must be in [0," << sz << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void MEDCouplingCartesianAMRMeshGen::checkFactorsAndIfNotSetAssign(const std::vector<int>& factors)
{
  if(getSpaceDimension()!=(int)factors.size())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::checkFactorsAndIfNotSetAssign : invalid size of factors ! size must be equal to the spaceDimension !");
  if(_factors.empty())
    {
      _factors=factors;
    }
  else
    {
      if(_factors!=factors)
        throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::checkFactorsAndIfNotSetAssign : the factors ");
    }
}

void MEDCouplingCartesianAMRMeshGen::retrieveGridsAtInternal(int lev, std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatchGen> >& grids) const
{
  if(lev==0)
    {
      const MEDCouplingCartesianAMRMesh *thisc(dynamic_cast<const MEDCouplingCartesianAMRMesh *>(this));//tony
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatchGF> elt(new MEDCouplingCartesianAMRPatchGF(const_cast<MEDCouplingCartesianAMRMesh *>(thisc)));
      grids.push_back(DynamicCastSafe<MEDCouplingCartesianAMRPatchGF,MEDCouplingCartesianAMRPatchGen>(elt));
    }
  else if(lev==1)
    {
      for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
        {
          const MEDCouplingCartesianAMRPatch *pt(*it);
          if(pt)
            {
              MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> tmp1(*it);
              grids.push_back(DynamicCastSafe<MEDCouplingCartesianAMRPatch,MEDCouplingCartesianAMRPatchGen>(tmp1));
            }
        }
    }
  else
    {
      for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
        {
          const MEDCouplingCartesianAMRPatch *pt(*it);
          if(pt)
            pt->getMesh()->retrieveGridsAtInternal(lev-1,grids);
        }
    }
}

/*!
 * \param [in,out] partBeforeFact - the part of a image mesh in compact format that will be put in refined reference.
 * \param [in] factors - the factors per axis.
 */
void MEDCouplingCartesianAMRMeshGen::ApplyFactorsOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, const std::vector<int>& factors)
{
  std::size_t sz(factors.size());
  if(sz!=partBeforeFact.size())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::ApplyFactorsOnCompactFrmt : size of input vectors must be the same !");
  for(std::size_t i=0;i<sz;i++)
    {
      partBeforeFact[i].first*=factors[i];
      partBeforeFact[i].second*=factors[i];
    }
}

/*!
 * \param [in,out] partBeforeFact - the part of a image mesh in compact format that will be put in ghost reference.
 * \param [in] ghostSize - the ghost size of zone for all axis.
 */
void MEDCouplingCartesianAMRMeshGen::ApplyGhostOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, int ghostSize)
{
  if(ghostSize<0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::ApplyGhostOnCompactFrmt : ghost size must be >= 0 !");
  std::size_t sz(partBeforeFact.size());
  for(std::size_t i=0;i<sz;i++)
    {
      partBeforeFact[i].first+=ghostSize;
      partBeforeFact[i].second+=ghostSize;
    }
}

/*!
 * This method is different than ApplyGhostOnCompactFrmt. The \a partBeforeFact parameter is enlarger contrary to ApplyGhostOnCompactFrmt.
 *
 * \param [in,out] partBeforeFact - the part of a image mesh in compact format that will be put in ghost reference.
 * \param [in] ghostSize - the ghost size of zone for all axis.
 */
void MEDCouplingCartesianAMRMeshGen::ApplyAllGhostOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, int ghostSize)
{
  if(ghostSize<0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::ApplyAllGhostOnCompactFrmt : ghost size must be >= 0 !");
  std::size_t sz(partBeforeFact.size());
  for(std::size_t i=0;i<sz;i++)
    {
      partBeforeFact[i].first-=ghostSize;
      partBeforeFact[i].second+=ghostSize;
    }
}

int MEDCouplingCartesianAMRMeshGen::GetGhostLevelInFineRef(int ghostLev, const std::vector<int>& factors)
{
  if(ghostLev<0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::GetGhostLevelInFineRef : the ghost size must be >=0 !");
  if(factors.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::GetGhostLevelInFineRef : no factors defined !");
  int ghostLevInPatchRef;
  if(ghostLev==0)
    ghostLevInPatchRef=0;
  else
    {
      ghostLevInPatchRef=(ghostLev-1)/factors[0]+1;
      for(std::size_t i=0;i<factors.size();i++)
        ghostLevInPatchRef=std::max(ghostLevInPatchRef,(ghostLev-1)/factors[i]+1);
    }
  return ghostLevInPatchRef;
}

/*!
 * This method returns a sub set of \a all. The subset is defined by the \a head in the tree defined by \a this.
 * Elements in \a all are expected to be sorted from god father to most refined structure.
 */
std::vector<const DataArrayDouble *> MEDCouplingCartesianAMRMeshGen::extractSubTreeFromGlobalFlatten(const MEDCouplingCartesianAMRMeshGen *head, const std::vector<const DataArrayDouble *>& all) const
{
  int maxLev(getMaxNumberOfLevelsRelativeToThis());
  std::vector<const DataArrayDouble *> ret;
  std::vector<const MEDCouplingCartesianAMRMeshGen *> meshes(1,this);
  std::size_t kk(0);
  for(int i=0;i<maxLev;i++)
    {
      std::vector<const MEDCouplingCartesianAMRMeshGen *> meshesTmp;
      for(std::vector<const MEDCouplingCartesianAMRMeshGen *>::const_iterator it=meshes.begin();it!=meshes.end();it++)
        {
          if((*it)==head || head->isObjectInTheProgeny(*it))
            ret.push_back(all[kk]);
          kk++;
          std::vector< const MEDCouplingCartesianAMRPatch *> ps((*it)->getPatches());
          for(std::vector< const MEDCouplingCartesianAMRPatch *>::const_iterator it0=ps.begin();it0!=ps.end();it0++)
            {
              const MEDCouplingCartesianAMRMeshGen *mesh((*it0)->getMesh());
              meshesTmp.push_back(mesh);
            }
        }
      meshes=meshesTmp;
    }
  if(kk!=all.size())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::extractSubTreeFromGlobalFlatten : the size of input vector is not compatible with number of leaves in this !");
  return ret;
}

std::size_t MEDCouplingCartesianAMRMeshGen::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDCouplingCartesianAMRMeshGen);
}

std::vector<const BigMemoryObject *> MEDCouplingCartesianAMRMeshGen::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  if((const MEDCouplingIMesh *)_mesh)
    ret.push_back((const MEDCouplingIMesh *)_mesh);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      if((const MEDCouplingCartesianAMRPatch*)*it)
        ret.push_back((const MEDCouplingCartesianAMRPatch*)*it);
    }
  return ret;
}

void MEDCouplingCartesianAMRMeshGen::updateTime() const
{
  if((const MEDCouplingIMesh *)_mesh)
    updateTimeWith(*_mesh);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      const MEDCouplingCartesianAMRPatch *elt(*it);
      if(!elt)
        continue;
      const MEDCouplingCartesianAMRMeshGen *mesh(elt->getMesh());
      if(mesh)
        updateTimeWith(*mesh);
    }
}

MEDCouplingCartesianAMRMeshSub::MEDCouplingCartesianAMRMeshSub(MEDCouplingCartesianAMRMeshGen *father, MEDCouplingIMesh *mesh):MEDCouplingCartesianAMRMeshGen(father,mesh)
{
}

MEDCouplingCartesianAMRMesh *MEDCouplingCartesianAMRMesh::New(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                              const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop)
{
  return new MEDCouplingCartesianAMRMesh(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop);
}

void MEDCouplingCartesianAMRMesh::setData(MEDCouplingDataForGodFather *data)
{
  MEDCouplingDataForGodFather *myData(_data);
  if(myData==data)
    return ;
  if(myData)
    myData->changeGodFather(0);
  _data=data;
  if(data)
    data->incrRef();
}

void MEDCouplingCartesianAMRMesh::allocData(int ghostLev) const
{
  checkData();
  _data->alloc(ghostLev);
}

void MEDCouplingCartesianAMRMesh::deallocData() const
{
  checkData();
  _data->dealloc();
}

MEDCouplingCartesianAMRMesh::MEDCouplingCartesianAMRMesh(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                         const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop):MEDCouplingCartesianAMRMeshGen(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop)
{
}

std::vector<const BigMemoryObject *> MEDCouplingCartesianAMRMesh::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret(MEDCouplingCartesianAMRMeshGen::getDirectChildren());
  const MEDCouplingDataForGodFather *pt(_data);
  if(pt)
    ret.push_back(pt);
  return ret;
}

void MEDCouplingCartesianAMRMesh::checkData() const
{
  const MEDCouplingDataForGodFather *data(_data);
  if(!data)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::checkData : no data set !");
}

MEDCouplingCartesianAMRMesh::~MEDCouplingCartesianAMRMesh()
{
  MEDCouplingDataForGodFather *data(_data);
  if(data)
    data->changeGodFather(0);
}
