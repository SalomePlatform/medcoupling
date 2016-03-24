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

#include "MEDCouplingCartesianAMRMesh.hxx"
#include "MEDCouplingAMRAttribute.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingUMesh.hxx"

#include <limits>
#include <sstream>
#include <numeric>

using namespace MEDCoupling;

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

MEDCouplingCartesianAMRPatchGen::MEDCouplingCartesianAMRPatchGen(const MEDCouplingCartesianAMRPatchGen& other, MEDCouplingCartesianAMRMeshGen *father):RefCountObject(other),_mesh(other._mesh)
{
  const MEDCouplingCartesianAMRMeshGen *mesh(other._mesh);
  if(mesh)
    _mesh=mesh->deepCopy(father);
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

std::vector<const BigMemoryObject *> MEDCouplingCartesianAMRPatchGen::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
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

MEDCouplingCartesianAMRPatch *MEDCouplingCartesianAMRPatch::deepCopy(MEDCouplingCartesianAMRMeshGen *father) const
{
  return new MEDCouplingCartesianAMRPatch(*this,father);
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
  return IsInMyNeighborhood(ghostLev==0?0:1,thisp,otherp);//make hypothesis that nb this->_mesh->getFather->getFactors() is >= ghostLev
}

/*!
 * This method states if \a other patch is in the neighborhood of \a this. The neighborhood zone is defined by \a ghostLev parameter
 * the must be >= 0. This method works even if \a this and \a other does not share the same father. But the level between their common
 * ancestor must be the same. If they don't have the same ancestor an exception will be thrown.
 *
 * \param [in] other - The other patch
 * \param [in] ghostLev - The size of the neighborhood zone.
 *
 * \throw if \a this or \a other are invalid (end before start).
 * \throw if \a ghostLev is \b not >= 0 .
 * \throw if \a this and \a other have not the same space dimension.
 * \throw if there is not common ancestor of \a this and \a other.
 *
 * \sa isInMyNeighborhood, isInMyNeighborhoodDiffLev
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
  otherp=MEDCouplingStructuredMesh::TranslateCompactFrmt(otherp,offset);
  return IsInMyNeighborhood(ghostLev,thisp,otherp);
}

/*!
 * This method states if \a other patch is in the neighborhood of \a this. The neighborhood zone is defined by \a ghostLev parameter
 * the must be >= 0. This method works even if \a this and \a other does not share the same father.
 * \a this is expected to be more refined than \a other. That is to say lev of \a this is greater than level of \a other.
 *
 * \param [in] other - The other patch
 * \param [in] ghostLev - The size of the neighborhood zone.
 *
 * \throw if \a this or \a other are invalid (end before start).
 * \throw if \a ghostLev is \b not >= 0 .
 * \throw if \a this and \a other have not the same space dimension.
 * \throw if there is not common ancestor of \a this and \a other.
 *
 * \sa isInMyNeighborhoodExt
 */
bool MEDCouplingCartesianAMRPatch::isInMyNeighborhoodDiffLev(const MEDCouplingCartesianAMRPatch *other, int ghostLev) const
{
  std::vector< std::pair<int,int> > thispp,otherpp;
  std::vector<int> factors;
  ComputeZonesOfTwoRelativeToOneDiffLev(ghostLev,this,other,thispp,otherpp,factors);
  return IsInMyNeighborhood(ghostLev>0?1:0,thispp,otherpp);//1 not ghostLev ! It is not a bug ( I hope :) ) ! Because as \a this is a refinement of \a other ghostLev is supposed to be <= factors
}

std::vector< std::pair<int,int> > MEDCouplingCartesianAMRPatch::getBLTRRangeRelativeToGF() const
{
  std::vector< std::pair<int,int> > ret(_bl_tr);
  const MEDCouplingCartesianAMRMeshGen *mesh(getMesh());
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::getBLTRRangeRelativeToGF : not valid !");
  const MEDCouplingCartesianAMRMeshGen *fath(mesh->getFather());
  if(!fath)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::getBLTRRangeRelativeToGF : not valid 2 !");
  std::vector<int> factors(fath->getFactors());
  std::size_t sz(ret.size());
  for(std::size_t ii=0;ii<sz;ii++)
    {
      ret[ii].first*=factors[ii];
      ret[ii].second*=factors[ii];
    }
  const MEDCouplingCartesianAMRMeshGen *oldFather(fath);
  fath=oldFather->getFather();
  while(fath)
    {
      int pos(fath->getPatchIdFromChildMesh(oldFather));
      const MEDCouplingCartesianAMRPatch *p(fath->getPatch(pos));
      const std::vector< std::pair<int,int> >& tmp(p->getBLTRRange());
      const std::vector<int>& factors2(fath->getFactors());
      std::transform(factors.begin(),factors.end(),factors2.begin(),factors.begin(),std::multiplies<int>());
      for(std::size_t ii=0;ii<sz;ii++)
        {
          int delta(ret[ii].second-ret[ii].first);
          ret[ii].first+=tmp[ii].first*factors[ii];
          ret[ii].second=ret[ii].first+delta;
        }
      oldFather=fath;
      fath=oldFather->getFather();
    }
  return ret;
}

std::vector<int> MEDCouplingCartesianAMRPatch::computeCellGridSt() const
{
  const MEDCouplingCartesianAMRMeshGen *m(getMesh());
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::computeCellGridSt : no mesh held by this !");
  const MEDCouplingCartesianAMRMeshGen *father(m->getFather());
  if(!father)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::computeCellGridSt : no father help by underlying mesh !");
  const std::vector< std::pair<int,int> >& bltr(getBLTRRange());
  const std::vector<int>& factors(father->getFactors());
  std::vector<int> ret(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(bltr));
  std::transform(ret.begin(),ret.end(),factors.begin(),ret.begin(),std::multiplies<int>());
  return ret;
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

/*!
 * \sa FindNeighborsOfSubPatchesOf
 */
std::vector< std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > > MEDCouplingCartesianAMRPatch::FindNeighborsOfSubPatchesOfSameLev(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2)
{
  if(!p1 || !p2)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::FindNeighborsOfSubPatchesOfSameLev : the input pointers must be not NULL !");
  std::vector< std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > > ret;
  std::vector< const MEDCouplingCartesianAMRPatch *> p1Work(p1->getMesh()->getPatches()),p2Work(p2->getMesh()->getPatches());
  while(!p1Work.empty())
    {
      std::vector< std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *> > retTmp;
      std::vector<const MEDCouplingCartesianAMRPatch *> p1Work2,p2Work2;
      for(std::vector<const MEDCouplingCartesianAMRPatch *>::const_iterator it1=p1Work.begin();it1!=p1Work.end();it1++)
        {
          for(std::vector<const MEDCouplingCartesianAMRPatch *>::const_iterator it2=p2Work.begin();it2!=p2Work.end();it2++)
            {
              if((*it1)->isInMyNeighborhoodExt(*it2,ghostLev>0?1:0))//1 not ghostLev ! It is not a bug ( I hope :) ) ! Because as \a this is a refinement of \a other ghostLev is supposed to be <= factors
                retTmp.push_back(std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *>(*it1,*it2));
            }
          std::vector<const MEDCouplingCartesianAMRPatch *> tmp1((*it1)->getMesh()->getPatches());
          p1Work2.insert(p1Work2.end(),tmp1.begin(),tmp1.end());
        }
      for(std::vector<const MEDCouplingCartesianAMRPatch *>::const_iterator it2=p2Work.begin();it2!=p2Work.end();it2++)
        {
          std::vector<const MEDCouplingCartesianAMRPatch *> tmp2((*it2)->getMesh()->getPatches());
          p2Work2.insert(p2Work2.end(),tmp2.begin(),tmp2.end());
        }
      ret.push_back(retTmp);
      p1Work=p1Work2;
      p2Work=p2Work2;
    }
  return ret;
}

/*!
 * This method returns all pair of patches (pa,pb) so that pb is in the neighborhood of pa (size of neighborhood is \a ghostLev).
 * pa is a refinement (a child) of \b p1 and pb is equal to \a p2. So the returned pair do not have the same level as it is the case for
 * FindNeighborsOfSubPatchesOfSameLev.
 *
 * \sa FindNeighborsOfSubPatchesOfSameLev
 */
void MEDCouplingCartesianAMRPatch::FindNeighborsOfSubPatchesOf(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, std::vector< std::pair<const MEDCouplingCartesianAMRPatch *, const MEDCouplingCartesianAMRPatch *> >& ret)
{
  if(!p1 || !p2)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::FindNeighborsOfSubPatchesOf : the input pointers must be not NULL !");
  std::vector< const MEDCouplingCartesianAMRPatch *> p1Work(p1->getMesh()->getPatches());
  while(!p1Work.empty())
    {
      std::vector<const MEDCouplingCartesianAMRPatch *> p1Work2;
      for(std::vector<const MEDCouplingCartesianAMRPatch *>::const_iterator it0=p1Work.begin();it0!=p1Work.end();it0++)
        {
          if((*it0)->isInMyNeighborhoodDiffLev(p2,ghostLev))
            ret.push_back(std::pair<const MEDCouplingCartesianAMRPatch *,const MEDCouplingCartesianAMRPatch *>(*it0,p2));
          std::vector<const MEDCouplingCartesianAMRPatch *> tmp2((*it0)->getMesh()->getPatches());
          p1Work2.insert(p1Work2.end(),tmp2.begin(),tmp2.end());
        }
      p1Work=p1Work2;
    }
}

/*!
 * \a p1 and \a p2 are expected to be neighbors (inside the \a ghostLev zone). This method updates \a dataOnP1 only in the ghost part using a part of \a dataOnP2.
 *
 * \saUpdateNeighborsOfOneWithTwoExt
 */
void MEDCouplingCartesianAMRPatch::UpdateNeighborsOfOneWithTwo(int ghostLev, const std::vector<int>& factors, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2)
{
  const std::vector< std::pair<int,int> >& p1BLTR(p1->getBLTRRange());
  const std::vector< std::pair<int,int> >& p2BLTR(p2->getBLTRRange());
  UpdateNeighborsOfOneWithTwoInternal(ghostLev,factors,p1BLTR,p2BLTR,dataOnP1,dataOnP2);
}

/*!
 * Idem than UpdateNeighborsOfOneWithTwo, except that here \a p1 and \a p2 are not sharing the same direct father.
 *
 * \sa UpdateNeighborsOfOneWithTwo
 */
void MEDCouplingCartesianAMRPatch::UpdateNeighborsOfOneWithTwoExt(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2)
{
  const std::vector< std::pair<int,int> >& p1BLTR(p1->getBLTRRange());//p1BLTR=[(10,12),(5,8)]
  std::vector< std::pair<int,int> > p2BLTR(p2->getBLTRRange());//p2BLTR=[(0,1),(0,5)]
  int lev(0);
  const MEDCouplingCartesianAMRMeshGen *ca(FindCommonAncestor(p1,p2,lev));
  std::vector<int> offset(ComputeOffsetFromTwoToOne(ca,lev,p1,p2));//[12,4]
  p2BLTR=MEDCouplingStructuredMesh::TranslateCompactFrmt(p2BLTR,offset);//p2BLTR=[(12,13),(4,9)]
  UpdateNeighborsOfOneWithTwoInternal(ghostLev,p1->getMesh()->getFather()->getFactors(),p1BLTR,p2BLTR,dataOnP1,dataOnP2);
}

/*!
 * \a p1 is expected to be more refined than \a p2. \a p1 and \a p2 have to share a common ancestor. Compared to UpdateNeighborsOfOneWithTwoExt here \a p1 and \a p2 are \b not at the same level !
 */
void MEDCouplingCartesianAMRPatch::UpdateNeighborsOfOneWithTwoMixedLev(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2, bool isConservative)
{
  std::vector< std::pair<int,int> > p1pp,p2pp;
  std::vector<int> factors;
  ComputeZonesOfTwoRelativeToOneDiffLev(ghostLev,p1,p2,p1pp,p2pp,factors);
  //
  std::vector<int> dimsP2NotRefined(p2->computeCellGridSt());
  std::vector<int> dimsP2Refined(dimsP2NotRefined);
  std::transform(dimsP2NotRefined.begin(),dimsP2NotRefined.end(),factors.begin(),dimsP2Refined.begin(),std::multiplies<int>());
  std::vector< std::pair<int,int> > p2RefinedAbs(MEDCouplingStructuredMesh::GetCompactFrmtFromDimensions(dimsP2NotRefined));
  std::vector<int> dimsP2RefinedGhost(dimsP2Refined.size());
  std::transform(dimsP2Refined.begin(),dimsP2Refined.end(),dimsP2RefinedGhost.begin(),std::bind2nd(std::plus<int>(),2*ghostLev));
  MCAuto<DataArrayDouble> fineP2(DataArrayDouble::New()); fineP2->alloc(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(dimsP2RefinedGhost),dataOnP2->getNumberOfComponents());
  MEDCouplingIMesh::SpreadCoarseToFineGhost(dataOnP2,dimsP2NotRefined,fineP2,p2RefinedAbs,factors,ghostLev);
  if(isConservative)
    {
      int fact(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(factors));
      std::transform(fineP2->begin(),fineP2->end(),fineP2->getPointer(),std::bind2nd(std::multiplies<double>(),1./((double)fact)));
    }
  //
  UpdateNeighborsOfOneWithTwoInternal(ghostLev,p1->getMesh()->getFather()->getFactors(),p1pp,p2pp,dataOnP1,fineP2);
}

/*!
 * \a p1 is expected to be more refined than \a p2. \a p1 and \a p2 have to share a common ancestor. Compared to UpdateNeighborsOfOneWithTwoExt here \a p1 and \a p2 are \b not at the same level !
 * This method has 3 outputs. 2 two first are the resp the position of \a p1 and \a p2 relative to \a p1. And \a factToApplyOn2 is the coeff of refinement to be applied on \a p2 to be virtualy
 * on the same level as \a p1.
 */
void MEDCouplingCartesianAMRPatch::ComputeZonesOfTwoRelativeToOneDiffLev(int ghostLev, const MEDCouplingCartesianAMRPatch *p1, const MEDCouplingCartesianAMRPatch *p2, std::vector< std::pair<int,int> >& p1Zone, std::vector< std::pair<int,int> >& p2Zone, std::vector<int>& factToApplyOn2)
{
  std::vector<const MEDCouplingCartesianAMRMeshGen *> ancestorsOfThis;
  const MEDCouplingCartesianAMRMeshGen *work(p1->getMesh()),*work2(0);
  ancestorsOfThis.push_back(work);
  while(work)
    {
      work=work->getFather();
      if(work)
        ancestorsOfThis.push_back(work);
    }
  //
  work=p2->getMesh();
  bool found(false);
  std::size_t levThis(0),levOther(0);
  while(work && !found)
    {
      work2=work;
      work=work->getFather();
      if(work)
        {
          levOther++;
          std::vector<const MEDCouplingCartesianAMRMeshGen *>::iterator it(std::find(ancestorsOfThis.begin(),ancestorsOfThis.end(),work));
          if(it!=ancestorsOfThis.end())
            {
              levThis=std::distance(ancestorsOfThis.begin(),it);
              found=true;
            }
        }
    }
  if(!found)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::ComputeZonesOfTwoRelativeToOneDiffLev : no common ancestor found !");
  if(levThis<=levOther)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::ComputeZonesOfTwoRelativeToOneDiffLev : this method is not called correctly !");
  //
  const MEDCouplingCartesianAMRMeshGen *comAncestor(ancestorsOfThis[levThis]);
  int idThis(comAncestor->getPatchIdFromChildMesh(ancestorsOfThis[levThis-1])),idOther(comAncestor->getPatchIdFromChildMesh(work2));
  const MEDCouplingCartesianAMRPatch *thisp(comAncestor->getPatch(idThis)),*otherp(comAncestor->getPatch(idOther));
  std::vector<int> offset(ComputeOffsetFromTwoToOne(comAncestor,levOther,thisp,otherp));
  p1Zone=thisp->getBLTRRange(); p2Zone=MEDCouplingStructuredMesh::TranslateCompactFrmt(otherp->getBLTRRange(),offset);
  factToApplyOn2.resize(p1Zone.size()); std::fill(factToApplyOn2.begin(),factToApplyOn2.end(),1);
  //
  std::size_t nbOfTurn(levThis-levOther);
  for(std::size_t i=0;i<nbOfTurn;i++)
    {
      std::vector< std::pair<int,int> > tmp0;
      MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt(p1Zone,p2Zone,tmp0,false);
      p2Zone=tmp0;
      const MEDCouplingCartesianAMRMeshGen *curAncestor(ancestorsOfThis[levThis-i]);
      ApplyFactorsOnCompactFrmt(p2Zone,curAncestor->getFactors());
      curAncestor=ancestorsOfThis[levThis-1-i];
      const std::vector<int>& factors(curAncestor->getFactors());
      std::transform(factToApplyOn2.begin(),factToApplyOn2.end(),factors.begin(),factToApplyOn2.begin(),std::multiplies<int>());
      int tmpId(curAncestor->getPatchIdFromChildMesh(ancestorsOfThis[levThis-2-i]));
      p1Zone=curAncestor->getPatch(tmpId)->getBLTRRange();
    }
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
  int zeLev(lev-1);
  int dim(p1->getMesh()->getSpaceDimension());
  if(p2->getMesh()->getSpaceDimension()!=dim)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::ComputeOffsetFromTwoToOne : dimension must be the same !");
  std::vector< int > ret(dim,0);
  for(int i=0;i<zeLev;i++)
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

void MEDCouplingCartesianAMRPatch::UpdateNeighborsOfOneWithTwoInternal(int ghostLev, const std::vector<int>& factors, const std::vector< std::pair<int,int> >&p1 ,const std::vector< std::pair<int,int> >&p2, DataArrayDouble *dataOnP1, const DataArrayDouble *dataOnP2)
{//p1=[(1,4),(2,4)] p2=[(4,5),(3,4)]
  int dim((int)factors.size());
  std::vector<int> dimsCoarse(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(p1));//[3,2]
  std::transform(dimsCoarse.begin(),dimsCoarse.end(),factors.begin(),dimsCoarse.begin(),std::multiplies<int>());//[12,8]
  std::transform(dimsCoarse.begin(),dimsCoarse.end(),dimsCoarse.begin(),std::bind2nd(std::plus<int>(),2*ghostLev));//[14,10]
  std::vector< std::pair<int,int> > rangeCoarse(MEDCouplingStructuredMesh::GetCompactFrmtFromDimensions(dimsCoarse));//[(0,14),(0,10)]
  std::vector<int> fakeFactors(dim,1);
  //
  std::vector< std::pair<int,int> > tmp0,tmp1,tmp2;
  MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt(p1,p2,tmp0,false);//tmp0=[(3,4),(1,2)]
  ApplyFactorsOnCompactFrmt(tmp0,factors);//tmp0=[(12,16),(4,8)]
  MEDCouplingStructuredMesh::ApplyGhostOnCompactFrmt(tmp0,ghostLev);//tmp0=[(13,17),(5,9)]
  std::vector< std::pair<int,int> > interstRange(MEDCouplingStructuredMesh::IntersectRanges(tmp0,rangeCoarse));//interstRange=[(13,14),(5,9)]
  MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt(p2,p1,tmp1,false);//tmp1=[(-3,0),(-1,1)]
  ApplyFactorsOnCompactFrmt(tmp1,factors);//tmp1=[(-12,-4),(-4,0)]
  MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt(tmp1,interstRange,tmp2,false);//tmp2=[(1,2),(1,5)]
  //
  std::vector< std::pair<int,int> > dimsFine(p2);
  ApplyFactorsOnCompactFrmt(dimsFine,factors);
  ApplyAllGhostOnCompactFrmt(dimsFine,ghostLev);
  //
  MCAuto<DataArrayDouble> ghostVals(MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(dimsFine),dataOnP2,tmp2));
  MEDCouplingIMesh::CondenseFineToCoarse(dimsCoarse,ghostVals,interstRange,fakeFactors,dataOnP1);
}

MEDCouplingCartesianAMRPatch::MEDCouplingCartesianAMRPatch(const MEDCouplingCartesianAMRPatch& other, MEDCouplingCartesianAMRMeshGen *father):MEDCouplingCartesianAMRPatchGen(other,father),_bl_tr(other._bl_tr)
{
}

/*!
 * \param [in,out] partBeforeFact - the part of a image mesh in compact format that will be put in refined reference.
 * \param [in] factors - the factors per axis.
 */
void MEDCouplingCartesianAMRPatch::ApplyFactorsOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, const std::vector<int>& factors)
{
  std::size_t sz(factors.size());
  if(sz!=partBeforeFact.size())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::ApplyFactorsOnCompactFrmt : size of input vectors must be the same !");
  for(std::size_t i=0;i<sz;i++)
    {
      partBeforeFact[i].first*=factors[i];
      partBeforeFact[i].second*=factors[i];
    }
}

/*!
 * This method is different than ApplyGhostOnCompactFrmt. The \a partBeforeFact parameter is enlarger contrary to ApplyGhostOnCompactFrmt.
 *
 * \param [in,out] partBeforeFact - the part of a image mesh in compact format that will be put in ghost reference.
 * \param [in] ghostSize - the ghost size of zone for all axis.
 */
void MEDCouplingCartesianAMRPatch::ApplyAllGhostOnCompactFrmt(std::vector< std::pair<int,int> >& partBeforeFact, int ghostSize)
{
  if(ghostSize<0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch::ApplyAllGhostOnCompactFrmt : ghost size must be >= 0 !");
  std::size_t sz(partBeforeFact.size());
  for(std::size_t i=0;i<sz;i++)
    {
      partBeforeFact[i].first-=ghostSize;
      partBeforeFact[i].second+=ghostSize;
    }
}

MEDCouplingCartesianAMRPatchGF::MEDCouplingCartesianAMRPatchGF(MEDCouplingCartesianAMRMesh *mesh):MEDCouplingCartesianAMRPatchGen(mesh)
{
}

MEDCouplingCartesianAMRPatchGF *MEDCouplingCartesianAMRPatchGF::deepCopy(MEDCouplingCartesianAMRMeshGen *father) const
{
  return new MEDCouplingCartesianAMRPatchGF(*this,father);
}

std::size_t MEDCouplingCartesianAMRPatchGF::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDCouplingCartesianAMRPatchGF);
}

MEDCouplingCartesianAMRPatchGF::MEDCouplingCartesianAMRPatchGF(const MEDCouplingCartesianAMRPatchGF& other, MEDCouplingCartesianAMRMeshGen *father):MEDCouplingCartesianAMRPatchGen(other,father)
{
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
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
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
  MCAuto<MEDCouplingIMesh> tmp(_mesh->buildWithGhost(ghostLev));
  return tmp->getNumberOfCells();
}

/*!
 * This method returns the number of cells including the current level but \b also \b including recursively all cells of other levels
 * starting from this. The set of cells which size is returned here are generally overlapping each other.
 */
int MEDCouplingCartesianAMRMeshGen::getNumberOfCellsRecursiveWithOverlap() const
{
  int ret(_mesh->getNumberOfCells());
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
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
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      ret-=(*it)->getNumberOfOverlapedCellsForFather();
      ret+=(*it)->getNumberOfCellsRecursiveWithoutOverlap();
    }
  return ret;
}

/*!
 * This method returns a vector of size equal to getAbsoluteLevelRelativeTo. It allows to find position an absolute position of \a this
 * relative to \a ref (that is typically the god father).
 *
 * \sa getPatchAtPosition
 */
std::vector<int> MEDCouplingCartesianAMRMeshGen::getPositionRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const
{
  if(!ref)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::getPositionRelativeTo : input pointer is NULL !");
  std::vector<int> ret;
  getPositionRelativeToInternal(ref,ret);
  std::reverse(ret.begin(),ret.end());
  return ret;
}

/*!
 * \sa getPositionRelativeTo, getMeshAtPosition
 */
const MEDCouplingCartesianAMRPatch *MEDCouplingCartesianAMRMeshGen::getPatchAtPosition(const std::vector<int>& pos) const
{
  std::size_t sz(pos.size());
  if(sz==0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::getPatchAtPosition : empty input -> no patch by definition !");
  int patchId(pos[0]);
  const MEDCouplingCartesianAMRPatch *elt(getPatch(patchId));
  if(sz==1)
    return elt;
  if(!elt || !elt->getMesh())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::getPatchAtPosition : NULL element found during walk !");
  std::vector<int> pos2(pos.begin()+1,pos.end());
  return elt->getMesh()->getPatchAtPosition(pos2);
}

const MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRMeshGen::getMeshAtPosition(const std::vector<int>& pos) const
{
  std::size_t sz(pos.size());
  if(sz==0)
    return this;
  int patchId(pos[0]);
  const MEDCouplingCartesianAMRPatch *elt(getPatch(patchId));
  if(sz==1)
    {
      if(!elt)
        throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::getMeshAtPosition : NULL patch !");
      return elt->getMesh();
    }
  if(!elt || !elt->getMesh())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::getPatchAtPosition : NULL element found during walk !");
  std::vector<int> pos2(pos.begin()+1,pos.end());
  return elt->getMesh()->getMeshAtPosition(pos2);
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
  return getGodFather()->retrieveGridsAt(absoluteLev);
}

/*!
 * \param [in] bottomLeftTopRight a vector equal to the space dimension of \a mesh that specifies for each dimension, the included cell start of the range for the first element of the pair,
 *                                a the end cell (\b excluded) of the range for the second element of the pair.
 * \param [in] factors The factor of refinement per axis (different from 0).
 */
void MEDCouplingCartesianAMRMeshGen::addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors)
{
  checkFactorsAndIfNotSetAssign(factors);
  MCAuto<MEDCouplingIMesh> mesh(static_cast<MEDCouplingIMesh *>(_mesh->buildStructuredSubPart(bottomLeftTopRight)));
  mesh->refineWithFactor(factors);
  MCAuto<MEDCouplingCartesianAMRMeshSub> zeMesh(new MEDCouplingCartesianAMRMeshSub(this,mesh));
  MCAuto<MEDCouplingCartesianAMRPatch> elt(new MEDCouplingCartesianAMRPatch(zeMesh,bottomLeftTopRight));
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
  void zipToFitOnCriterion(int minPatchLgth);
  void updateNumberOfTrue() const;
  MCAuto<InternalPatch> extractPart(const std::vector< std::pair<int,int> >&partInGlobal) const;
  MCAuto<InternalPatch> deepCopy() const;
protected:
  ~InternalPatch() { }
private:
  mutable int _nb_of_true;
  std::vector<bool> _crit;
  //! _part is global
  std::vector< std::pair<int,int> > _part;
};

void InternalPatch::zipToFitOnCriterion(int minPatchLgth)
{
  std::vector<int> cgs(computeCGS());
  std::vector<bool> newCrit;
  std::vector< std::pair<int,int> > newPart,newPart2;
  int newNbOfTrue(MEDCouplingStructuredMesh::FindMinimalPartOf(minPatchLgth,cgs,_crit,newCrit,newPart));
  MEDCouplingStructuredMesh::ChangeReferenceToGlobalOfCompactFrmt(_part,newPart,newPart2);
  if(newNbOfTrue!=_nb_of_true)
    throw INTERP_KERNEL::Exception("InternalPatch::zipToFitOnCrit : internal error !");
  _crit=newCrit; _part=newPart2;
}

void InternalPatch::updateNumberOfTrue() const
{
  _nb_of_true=(int)std::count(_crit.begin(),_crit.end(),true);
}

MCAuto<InternalPatch> InternalPatch::extractPart(const std::vector< std::pair<int,int> >&partInGlobal) const
{
  MCAuto<InternalPatch> ret(new InternalPatch);
  std::vector<int> cgs(computeCGS());
  std::vector< std::pair<int,int> > newPart;
  MEDCouplingStructuredMesh::ChangeReferenceFromGlobalOfCompactFrmt(_part,partInGlobal,newPart);
  MEDCouplingStructuredMesh::ExtractFieldOfBoolFrom(cgs,_crit,newPart,ret->getCriterion());
  ret->setPart(partInGlobal);
  ret->updateNumberOfTrue();
  return ret;
}

MCAuto<InternalPatch> InternalPatch::deepCopy() const
{
  MCAuto<InternalPatch> ret(new InternalPatch);
  (*ret)=*this;
  return ret;
}

void DissectBigPatch(const INTERP_KERNEL::BoxSplittingOptions& bso, const InternalPatch *patchToBeSplit, int axisId, int largestLength, int& cutPlace)
{
  int minimumPatchLength(bso.getMinimumPatchLength());
  std::vector<double> ratio(largestLength-minimumPatchLength,std::numeric_limits<double>::max());
  int index_min = -1;
  double minSemiEfficiencyRatio(std::numeric_limits<double>::max());
  double efficiencyPerAxis[2];

  for(int i=minimumPatchLength-1;i<largestLength-minimumPatchLength;i++)
    {
      for(int h=0;h<2;h++)
        {
          std::vector< std::pair<int,int> > rectH(patchToBeSplit->getConstPart());
          if(h==0)
            rectH[axisId].second=patchToBeSplit->getConstPart()[axisId].first+i;
          else
            rectH[axisId].first=patchToBeSplit->getConstPart()[axisId].first+i;
          MCAuto<InternalPatch> p(patchToBeSplit->deepCopy());
          p->zipToFitOnCriterion(bso.getMinimumPatchLength());
          efficiencyPerAxis[h]=p->getEfficiencyPerAxis(axisId);
        }
      ratio[i]=std::max(efficiencyPerAxis[0], efficiencyPerAxis[1]) / std::min(efficiencyPerAxis[0], efficiencyPerAxis[1]);
      if(ratio[i]<minSemiEfficiencyRatio)
        {
          minSemiEfficiencyRatio = ratio[i];
          index_min = i;
        }
    }

  if(index_min==-1)
    throw INTERP_KERNEL::Exception("DissectBigPatch : just call to Arthur !");

  cutPlace=index_min+patchToBeSplit->getConstPart()[axisId].first;
}

bool FindHole(const INTERP_KERNEL::BoxSplittingOptions& bso, const InternalPatch *patchToBeSplit, int axisId, int& cutPlace)
{
  cutPlace=-1;
  int minimumPatchLength(bso.getMinimumPatchLength());
  const int dim(patchToBeSplit->getDimension());
  std::vector< std::vector<int> > signatures(patchToBeSplit->computeSignature());
  for(int id=0;id<dim;id++)
    {
      const std::vector<int>& signature(signatures[id]);
      std::vector<int> hole;
      std::vector<double> distance;
      int len((int)signature.size());
      for(int i=minimumPatchLength-1;i<len-minimumPatchLength;i++)
        if(signature[i]==0)
          hole.push_back(i);
      if(!hole.empty())
        {
          int closestHoleToMiddle(hole[0]);
          int oldDistanceToMiddle(std::abs(hole[0]-len/2));
          int newDistanceToMiddle(oldDistanceToMiddle);
          for(std::size_t i=0;i<hole.size();i++)
            {
              newDistanceToMiddle=std::abs(hole[i]-len/2);
              if(newDistanceToMiddle < oldDistanceToMiddle)
                {
                  oldDistanceToMiddle = newDistanceToMiddle;
                  closestHoleToMiddle = hole[i];
                }
            }
          cutPlace=closestHoleToMiddle+patchToBeSplit->getConstPart()[axisId].first;
          return true;
        }
    }
  return false;
}

bool FindInflection(const INTERP_KERNEL::BoxSplittingOptions& bso, const InternalPatch *patchToBeSplit, int& cutPlace, int& axisId)
{
  bool cutFound(false); cutPlace=-1;// do not set axisId before to be sure that cutFound was set to true
  const std::vector< std::pair<int,int> >& part(patchToBeSplit->getConstPart());
  int sign,minimumPatchLength(bso.getMinimumPatchLength());
  const int dim(patchToBeSplit->getDimension());

  std::vector<int> zeroCrossDims(dim,-1);
  std::vector<int> zeroCrossVals(dim,-1);
  std::vector< std::vector<int> > signatures(patchToBeSplit->computeSignature());
  for (int id=0;id<dim;id++)
    {
      const std::vector<int>& signature(signatures[id]);

      std::vector<int> derivate_second_order,gradient_absolute,zero_cross,edge,max_cross_list ;
      std::vector<double> distance ;

      for(std::size_t i=1;i<signature.size()-1;i++)
        derivate_second_order.push_back(signature[i-1]-2*signature[i]+signature[i+1]) ;

      // Gradient absolute value
      for(std::size_t i=1;i<derivate_second_order.size();i++)
        gradient_absolute.push_back(abs(derivate_second_order[i]-derivate_second_order[i-1])) ;
      if(derivate_second_order.empty())
        continue;
      for(std::size_t i=1;i<derivate_second_order.size()-1;i++)
        {
          if (derivate_second_order[i]*derivate_second_order[i+1] < 0 )
            sign = -1 ;
          if (derivate_second_order[i]*derivate_second_order[i+1] > 0 )
            sign = 1 ;
          if (derivate_second_order[i]*derivate_second_order[i+1] == 0 )
            sign = 0 ;
          if ( sign==0 || sign==-1 )
            if ( i >= (std::size_t)minimumPatchLength-2 && i <= signature.size()-minimumPatchLength-2 )
              {
                zero_cross.push_back(i) ;
                edge.push_back(gradient_absolute[i]) ;
              }
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

          double distance_min=*min_element(distance.begin(),distance.end()) ;
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
  return cutFound;
}

bool TryAction4(const INTERP_KERNEL::BoxSplittingOptions& bso, const InternalPatch *patchToBeSplit, int axisId, int rangeOfAxisId, int& cutPlace)
{
  if(patchToBeSplit->getEfficiency()<=bso.getEfficiencyGoal())
    {
      if(rangeOfAxisId>=2*bso.getMinimumPatchLength())
        {
          cutPlace=rangeOfAxisId/2+patchToBeSplit->getConstPart()[axisId].first-1;
        }
      else
        return false;
    }
  else
    {
      if(patchToBeSplit->getNumberOfCells()>bso.getMaximumNbOfCellsInPatch() || rangeOfAxisId>bso.getMaximumPatchLength())
        {
          DissectBigPatch(bso,patchToBeSplit,axisId,rangeOfAxisId,cutPlace);
        }
      else
        return false;
    }
  return true;
}

MCAuto<InternalPatch> DealWithNoCut(const InternalPatch *patch)
{
  MCAuto<InternalPatch> ret(const_cast<InternalPatch *>(patch));
  ret->incrRef();
  return ret;
}

void DealWithCut(double minPatchLgth, const InternalPatch *patchToBeSplit, int axisId, int cutPlace, std::vector<MCAuto<InternalPatch> >& listOfPatches)
{
  MCAuto<InternalPatch> leftPart,rightPart;
  std::vector< std::pair<int,int> > rect(patchToBeSplit->getConstPart());
  std::vector< std::pair<int,int> > leftRect(rect),rightRect(rect);
  leftRect[axisId].second=cutPlace+1;
  rightRect[axisId].first=cutPlace+1;
  leftPart=patchToBeSplit->extractPart(leftRect);
  rightPart=patchToBeSplit->extractPart(rightRect);
  leftPart->zipToFitOnCriterion(minPatchLgth); rightPart->zipToFitOnCriterion(minPatchLgth);
  listOfPatches.push_back(leftPart);
  listOfPatches.push_back(rightPart);
}

/// @endcond

void MEDCouplingCartesianAMRMeshGen::removeAllPatches()
{
  _patches.clear();
  declareAsNew();
}

void MEDCouplingCartesianAMRMeshGen::removePatch(int patchId)
{
  checkPatchId(patchId);
  int sz((int)_patches.size()),j(0);
  std::vector< MCAuto<MEDCouplingCartesianAMRPatch> > patches(sz-1);
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

/*!
 * This method is a generic algorithm to create patches in \a this (by destroying the patches if any).
 * This method uses \a criterion array as a field on cells on this level.
 * This method only create patches at level 0 relative to \a this.
 *
 * This generic algorithm can be degenerated into three child ones, depending on the arguments given; in particular depending
 * on whether they are equal to 0 or not.
 * 1/ If  minimumPatchLength = maximumPatchLength = maximumPatchVolume = 0, then we have the Berger-Rigoutsos algorithm.
 * This algorithm was developed in 1991 and seems appropriate for sequential programming.
 * 2/ If maximumPatchLength = 0, then we have the Livne algorithm.
 * This algorithm was developed in 2004 and is an improvement of the Berger-Rigoutsos algorithm.
 * 3/ If maximumPatchVolume = 0, the we have the lmin-lmax algorithm.
 * This algorithm was developed by Arthur TALPAERT in 2014 and is an improvement of the Livne algorithm. It is especially
 * appropriate for parallel computing, where one patch would be given to one CPU. See Arthur TALPAERT's 2014 CANUM poster
 * for more information.
 *
 */
void MEDCouplingCartesianAMRMeshGen::createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const std::vector<bool>& criterion, const std::vector<int>& factors)
{
  int nbCells(getNumberOfCellsAtCurrentLevel());
  if(nbCells!=(int)criterion.size())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::createPatchesFromCriterion : the number of tuples of criterion array must be equal to the number of cells at the current level !");
  _patches.clear();
  std::vector<int> cgs(_mesh->getCellGridStructure());
  std::vector< MCAuto<InternalPatch> > listOfPatches,listOfPatchesOK;
  //
  MCAuto<InternalPatch> p(new InternalPatch);
  p->setNumberOfTrue(MEDCouplingStructuredMesh::FindMinimalPartOf(bso.getMinimumPatchLength(),cgs,criterion,p->getCriterion(),p->getPart()));
  if(p->presenceOfTrue())
    listOfPatches.push_back(p);
  while(!listOfPatches.empty())
    {
      std::vector< MCAuto<InternalPatch> > listOfPatchesTmp;
      for(std::vector< MCAuto<InternalPatch> >::iterator it=listOfPatches.begin();it!=listOfPatches.end();it++)
        {
          //
          int axisId,largestLength,cutPlace;
          MEDCouplingStructuredMesh::FindTheWidestAxisOfGivenRangeInCompactFrmt((*it)->getConstPart(),axisId,largestLength);
          if((*it)->getEfficiency()>=bso.getEfficiencyThreshold() && ((*it)->getNumberOfCells()>bso.getMaximumNbOfCellsInPatch() || largestLength>bso.getMaximumPatchLength()))
            {
              DissectBigPatch(bso,*it,axisId,largestLength,cutPlace);
              DealWithCut(bso.getMinimumPatchLength(),*it,axisId,cutPlace,listOfPatchesTmp);
              continue;
            }//action 1
          if(FindHole(bso,*it,axisId,cutPlace))//axisId overwritten here if FindHole equal to true !
            { DealWithCut(bso.getMinimumPatchLength(),*it,axisId,cutPlace,listOfPatchesTmp); continue; }//action 2
          if(FindInflection(bso,*it,cutPlace,axisId))//axisId overwritten here if cutFound equal to true !
            { DealWithCut(bso.getMinimumPatchLength(),*it,axisId,cutPlace,listOfPatchesTmp); continue; }//action 3
          if(TryAction4(bso,*it,axisId,largestLength,cutPlace))
            { DealWithCut(bso.getMinimumPatchLength(),*it,axisId,cutPlace,listOfPatchesTmp); continue; }//action 4
          else
            listOfPatchesOK.push_back(DealWithNoCut(*it));
        }
      listOfPatches=listOfPatchesTmp;
    }
  for(std::vector< MCAuto<InternalPatch> >::const_iterator it=listOfPatchesOK.begin();it!=listOfPatchesOK.end();it++)
    addPatch((*it)->getConstPart(),factors);
  declareAsNew();
}

/*!
 * This method creates patches in \a this (by destroying the patches if any). This method uses \a criterion array as a field on cells on this level.
 * This method only create patches at level 0 relative to \a this.
 */
void MEDCouplingCartesianAMRMeshGen::createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayByte *criterion, const std::vector<int>& factors)
{
  if(!criterion || !criterion->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::createPatchesFromCriterion : the criterion DataArrayByte instance must be allocated and not NULL !");
  std::vector<bool> crit(criterion->toVectorOfBool());//check that criterion has one component.
  createPatchesFromCriterion(bso,crit,factors);
  declareAsNew();
}

void MEDCouplingCartesianAMRMeshGen::createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayDouble *criterion, const std::vector<int>& factors, double eps)
{
  if(!criterion)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen::createPatchesFromCriterion : null criterion pointer !");
  std::vector<bool> inp(criterion->toVectorOfBool(eps));
  createPatchesFromCriterion(bso,inp,factors);
}

int MEDCouplingCartesianAMRMeshGen::getPatchIdFromChildMesh(const MEDCouplingCartesianAMRMeshGen *mesh) const
{
  int ret(0);
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++,ret++)
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
  return p1->isInMyNeighborhood(p2,ghostLev);
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
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(fine->getNumberOfCells(),cellFieldOnThis->getNumberOfComponents());
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
void MEDCouplingCartesianAMRMeshGen::fillCellFieldOnPatch(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, bool isConservative) const
{
  if(!cellFieldOnThis || !cellFieldOnThis->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createCellFieldOnPatch : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  MEDCouplingIMesh::SpreadCoarseToFine(cellFieldOnThis,_mesh->getCellGridStructure(),cellFieldOnPatch,patch->getBLTRRange(),getFactors());
  if(isConservative)
    {
      int fact(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(getFactors()));
      std::transform(cellFieldOnPatch->begin(),cellFieldOnPatch->end(),cellFieldOnPatch->getPointer(),std::bind2nd(std::multiplies<double>(),1./((double)fact)));
    }
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
void MEDCouplingCartesianAMRMeshGen::fillCellFieldOnPatchGhost(int patchId, const DataArrayDouble *cellFieldOnThis, DataArrayDouble *cellFieldOnPatch, int ghostLev, bool isConservative) const
{
  if(!cellFieldOnThis || !cellFieldOnThis->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createCellFieldOnPatchGhost : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  MEDCouplingIMesh::SpreadCoarseToFineGhost(cellFieldOnThis,_mesh->getCellGridStructure(),cellFieldOnPatch,patch->getBLTRRange(),getFactors(),ghostLev);
  if(isConservative)
    {
      int fact(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(getFactors()));
      std::transform(cellFieldOnPatch->begin(),cellFieldOnPatch->end(),cellFieldOnPatch->getPointer(),std::bind2nd(std::multiplies<double>(),1./((double)fact)));
    }
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
void MEDCouplingCartesianAMRMeshGen::fillCellFieldOnPatchGhostAdv(int patchId, const DataArrayDouble *cellFieldOnThis, int ghostLev, const std::vector<const DataArrayDouble *>& arrsOnPatches, bool isConservative) const
{
  int nbp(getNumberOfPatches());
  if(nbp!=(int)arrsOnPatches.size())
    {
      std::ostringstream oss; oss << "MEDCouplingCartesianAMRMesh::fillCellFieldOnPatchGhostAdv : there are " << nbp << " patches in this and " << arrsOnPatches.size() << " arrays in the last parameter !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  DataArrayDouble *theFieldToFill(const_cast<DataArrayDouble *>(arrsOnPatches[patchId]));
  // first, do as usual
  fillCellFieldOnPatchGhost(patchId,cellFieldOnThis,theFieldToFill,ghostLev,isConservative);
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
  int nbp(getNumberOfPatches());
  if(nbp!=(int)arrsOnPatches.size())
    {
      std::ostringstream oss; oss << "MEDCouplingCartesianAMRMesh::fillCellFieldOnPatchOnlyGhostAdv : there are " << nbp << " patches in this and " << arrsOnPatches.size() << " arrays in the last parameter !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const MEDCouplingCartesianAMRPatch *refP(getPatch(patchId));
  DataArrayDouble *theFieldToFill(const_cast<DataArrayDouble *>(arrsOnPatches[patchId]));
  std::vector<int> ids(getPatchIdsInTheNeighborhoodOf(patchId,ghostLev));
  for(std::vector<int>::const_iterator it=ids.begin();it!=ids.end();it++)
    {
      const MEDCouplingCartesianAMRPatch *otherP(getPatch(*it));
      MEDCouplingCartesianAMRPatch::UpdateNeighborsOfOneWithTwo(ghostLev,_factors,refP,otherP,theFieldToFill,arrsOnPatches[*it]);
    }
}

void MEDCouplingCartesianAMRMeshGen::fillCellFieldOnPatchOnlyOnGhostZoneWith(int ghostLev, const MEDCouplingCartesianAMRPatch *patchToBeModified, const MEDCouplingCartesianAMRPatch *neighborPatch, DataArrayDouble *cellFieldOnPatch, const DataArrayDouble *cellFieldNeighbor) const
{
  MEDCouplingCartesianAMRPatch::UpdateNeighborsOfOneWithTwo(ghostLev,_factors,patchToBeModified,neighborPatch,cellFieldOnPatch,cellFieldNeighbor);
}

/*!
 * This method updates \a cellFieldOnThis part of values coming from the cell field \a cellFieldOnPatch lying on patch having id \a patchId.
 *
 * \param [in] patchId - The id of the patch \a cellFieldOnThis has to be put on.
 * \param [in] cellFieldOnPatch - The array of the cell field on patch with id \a patchId.
 * \param [in,out] cellFieldOnThis The array of the cell field on \a this to be updated only on the part concerning the patch with id \a patchId.
 * \param [in] isConservative - true if the field needs to be conserved. false if maximum principle has to be applied.
 *
 * \throw if \a patchId is not in [ 0 , \c this->getNumberOfPatches() )
 * \throw if \a cellFieldOnPatch is NULL or not allocated
 * \sa createCellFieldOnPatch, MEDCouplingIMesh::CondenseFineToCoarse,fillCellFieldComingFromPatchGhost
 */
void MEDCouplingCartesianAMRMeshGen::fillCellFieldComingFromPatch(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, bool isConservative) const
{
  if(!cellFieldOnPatch || !cellFieldOnPatch->isAllocated())
      throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::fillCellFieldComingFromPatch : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  MEDCouplingIMesh::CondenseFineToCoarse(_mesh->getCellGridStructure(),cellFieldOnPatch,patch->getBLTRRange(),getFactors(),cellFieldOnThis);
  if(!isConservative)
    {
      int fact(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(getFactors()));
      MEDCouplingStructuredMesh::MultiplyPartOf(_mesh->getCellGridStructure(),patch->getBLTRRange(),1./((double)fact),cellFieldOnThis);
    }
}

/*!
 * This method is the extension of MEDCouplingCartesianAMRMesh::fillCellFieldComingFromPatch managing the ghost cells. If this
 * method is called with \a ghostLev equal to 0 it behaves exactly as MEDCouplingCartesianAMRMesh::fillCellFieldComingFromPatch.
 *
 * \param [in] patchId - The id of the patch \a cellFieldOnThis has to be put on.
 * \param [in] cellFieldOnPatch - The array of the cell field on patch with id \a patchId.
 * \param [in,out] cellFieldOnThis The array of the cell field on \a this to be updated only on the part concerning the patch with id \a patchId.
 * \param [in] ghostLev The size of ghost zone (must be >= 0 !)
 * \param [in] isConservative - true if the field needs to be conserved. false if maximum principle has to be applied.
 *
 * \throw if \a patchId is not in [ 0 , \c this->getNumberOfPatches() )
 * \throw if \a cellFieldOnPatch is NULL or not allocated
 * \sa fillCellFieldComingFromPatch
 */
void MEDCouplingCartesianAMRMeshGen::fillCellFieldComingFromPatchGhost(int patchId, const DataArrayDouble *cellFieldOnPatch, DataArrayDouble *cellFieldOnThis, int ghostLev, bool isConservative) const
{
  if(!cellFieldOnPatch || !cellFieldOnPatch->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::fillCellFieldComingFromPatchGhost : the input cell field array is NULL or not allocated !");
  const MEDCouplingCartesianAMRPatch *patch(getPatch(patchId));
  MEDCouplingIMesh::CondenseFineToCoarseGhost(_mesh->getCellGridStructure(),cellFieldOnPatch,patch->getBLTRRange(),getFactors(),cellFieldOnThis,ghostLev);
  if(!isConservative)
    {
      int fact(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(getFactors()));
      MEDCouplingStructuredMesh::MultiplyPartOfByGhost(_mesh->getCellGridStructure(),patch->getBLTRRange(),ghostLev,1./((double)fact),cellFieldOnThis);
    }
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
  MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
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
  MCAuto<MEDCouplingUMesh> part(_mesh->buildUnstructured());
  std::vector<bool> bs(_mesh->getNumberOfCells(),false);
  std::vector<int> cgs(_mesh->getCellGridStructure());
  std::vector< MCAuto<MEDCouplingUMesh> > msSafe(_patches.size()+1);
  std::size_t ii(0);
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++,ii++)
    {
      MEDCouplingStructuredMesh::SwitchOnIdsFrom(cgs,(*it)->getBLTRRange(),bs);
      msSafe[ii+1]=(*it)->getMesh()->buildUnstructured();
    }
  MCAuto<DataArrayInt> eltsOff(DataArrayInt::BuildListOfSwitchedOff(bs));
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
  std::vector< MCAuto<MEDCoupling1SGTUMesh> > cellsSafe;
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      const MEDCouplingCartesianAMRPatch *patch(*it);
      if(patch)
        {
          MCAuto<MEDCouplingIMesh> cell(patch->getMesh()->getImageMesh()->asSingleCell());
          MCAuto<MEDCoupling1SGTUMesh> cell1SGT(cell->build1SGTUnstructured());
          cellsSafe.push_back(cell1SGT); cells.push_back(cell1SGT);
        }
    }
  return MEDCoupling1SGTUMesh::Merge1SGTUMeshes(cells);
}

MEDCoupling1SGTUMesh *MEDCouplingCartesianAMRMeshGen::buildMeshOfDirectChildrenOnly() const
{
  std::vector<const MEDCoupling1SGTUMesh *> patches;
  std::vector< MCAuto<MEDCoupling1SGTUMesh> > patchesSafe;
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
      {
        const MEDCouplingCartesianAMRPatch *patch(*it);
        if(patch)
          {
            MCAuto<MEDCoupling1SGTUMesh> patchMesh(patch->getMesh()->getImageMesh()->build1SGTUnstructured());
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
  std::vector< MCAuto<MEDCouplingFieldDouble> > msSafe(_patches.size()+1);
  std::size_t ii(0);
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++,ii++)
    {
      MEDCouplingStructuredMesh::SwitchOnIdsFrom(cgs,(*it)->getBLTRRange(),bs);
      std::vector<const DataArrayDouble *> tmpArrs(extractSubTreeFromGlobalFlatten((*it)->getMesh(),recurseArrs));
      msSafe[ii+1]=(*it)->getMesh()->buildCellFieldOnRecurseWithoutOverlapWithoutGhost(ghostSz,tmpArrs);
    }
  MCAuto<DataArrayInt> eltsOff(DataArrayInt::BuildListOfSwitchedOff(bs));
  //
  MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_CELLS));
  MCAuto<DataArrayDouble> arr2(extractGhostFrom(ghostSz,recurseArrs[0]));
  arr2=arr2->selectByTupleIdSafe(eltsOff->begin(),eltsOff->end());
  ret->setArray(arr2);
  ret->setName(arr2->getName());
  MCAuto<MEDCouplingUMesh> part(_mesh->buildUnstructured());
  MCAuto<MEDCouplingMesh> mesh(part->buildPartOfMySelf(eltsOff->begin(),eltsOff->end(),false));
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
  MEDCouplingStructuredMesh::ApplyGhostOnCompactFrmt(p,ghostSz);
  MCAuto<DataArrayDouble> ret(MEDCouplingStructuredMesh::ExtractFieldOfDoubleFrom(st,arr,p));
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

/*!
 * This method returns a dump python of \a this. It is useful for users of createPatchesFromCriterion method for debugging.
 *
 * \sa dumpPatchesOf, createPatchesFromCriterion, createPatchesFromCriterionML
 */
std::string MEDCouplingCartesianAMRMeshGen::buildPythonDumpOfThis() const
{
  std::ostringstream oss;
  oss << "amr=MEDCouplingCartesianAMRMesh(\""<< getImageMesh()->getName() << "\"," << getSpaceDimension() << ",[";
  std::vector<int> ngs(getImageMesh()->getNodeGridStructure());
  std::vector<double> orig(getImageMesh()->getOrigin()),dxyz(getImageMesh()->getDXYZ());
  std::copy(ngs.begin(),ngs.end(),std::ostream_iterator<int>(oss,","));
  oss <<  "],[";
  std::copy(orig.begin(),orig.end(),std::ostream_iterator<double>(oss,","));
  oss << "],[";
  std::copy(dxyz.begin(),dxyz.end(),std::ostream_iterator<double>(oss,","));
  oss << "])\n";
  dumpPatchesOf("amr",oss);
  return oss.str();
}

MEDCouplingCartesianAMRMeshGen::MEDCouplingCartesianAMRMeshGen(const MEDCouplingCartesianAMRMeshGen& other):RefCountObject(other),_mesh(other._mesh),_patches(other._patches),_factors(other._factors)
{
  const MEDCouplingIMesh *mesh(other._mesh);
  if(mesh)
    _mesh=static_cast<MEDCouplingIMesh *>(mesh->deepCopy());
  std::size_t sz(other._patches.size());
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDCouplingCartesianAMRPatch *patch(other._patches[i]);
      if(patch)
        _patches[i]=patch->deepCopy(this);
    }
}

MEDCouplingCartesianAMRMeshGen::MEDCouplingCartesianAMRMeshGen(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                               const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop)
{
  _mesh=MEDCouplingIMesh::New(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop);
}

MEDCouplingCartesianAMRMeshGen::MEDCouplingCartesianAMRMeshGen(MEDCouplingIMesh *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshGen(MEDCouplingIMesh *mesh) constructor : The input mesh is null !");
  mesh->checkConsistencyLight();
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

void MEDCouplingCartesianAMRMeshGen::retrieveGridsAtInternal(int lev, std::vector< MCAuto<MEDCouplingCartesianAMRPatchGen> >& grids) const
{
  if(lev==0)
    {
      const MEDCouplingCartesianAMRMesh *thisc(dynamic_cast<const MEDCouplingCartesianAMRMesh *>(this));//tony
      MCAuto<MEDCouplingCartesianAMRPatchGF> elt(new MEDCouplingCartesianAMRPatchGF(const_cast<MEDCouplingCartesianAMRMesh *>(thisc)));
      grids.push_back(DynamicCastSafe<MEDCouplingCartesianAMRPatchGF,MEDCouplingCartesianAMRPatchGen>(elt));
    }
  else if(lev==1)
    {
      for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
        {
          const MEDCouplingCartesianAMRPatch *pt(*it);
          if(pt)
            {
              MCAuto<MEDCouplingCartesianAMRPatch> tmp1(*it);
              grids.push_back(DynamicCastSafe<MEDCouplingCartesianAMRPatch,MEDCouplingCartesianAMRPatchGen>(tmp1));
            }
        }
    }
  else
    {
      for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
        {
          const MEDCouplingCartesianAMRPatch *pt(*it);
          if(pt)
            pt->getMesh()->retrieveGridsAtInternal(lev-1,grids);
        }
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

void MEDCouplingCartesianAMRMeshGen::dumpPatchesOf(const std::string& varName, std::ostream& oss) const
{
  std::size_t j(0);
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      const MEDCouplingCartesianAMRPatch *patch(*it);
      if(patch)
        {
          std::ostringstream oss2; oss2 << varName << ".addPatch([";
          const std::vector< std::pair<int,int> >& bltr(patch->getBLTRRange());
          std::size_t sz(bltr.size());
          for(std::size_t i=0;i<sz;i++)
            {
              oss2 << "(" << bltr[i].first << "," << bltr[i].second << ")";
              if(i!=sz-1)
                oss2 << ",";
            }
          oss2 << "],[";
          std::copy(_factors.begin(),_factors.end(),std::ostream_iterator<int>(oss2,","));
          oss2 << "])\n";
          oss << oss2.str();
          std::ostringstream oss3; oss3 << varName << "[" << j++ << "]";
          patch->getMesh()->dumpPatchesOf(oss3.str(),oss);
        }
    }
}

std::size_t MEDCouplingCartesianAMRMeshGen::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDCouplingCartesianAMRMeshGen);
}

std::vector<const BigMemoryObject *> MEDCouplingCartesianAMRMeshGen::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back((const MEDCouplingIMesh *)_mesh);
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    ret.push_back((const MEDCouplingCartesianAMRPatch*)*it);
  return ret;
}

void MEDCouplingCartesianAMRMeshGen::updateTime() const
{
  if((const MEDCouplingIMesh *)_mesh)
    updateTimeWith(*_mesh);
  for(std::vector< MCAuto<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      const MEDCouplingCartesianAMRPatch *elt(*it);
      if(!elt)
        continue;
      const MEDCouplingCartesianAMRMeshGen *mesh(elt->getMesh());
      if(mesh)
        updateTimeWith(*mesh);
    }
}

MEDCouplingCartesianAMRMeshSub::MEDCouplingCartesianAMRMeshSub(MEDCouplingCartesianAMRMeshGen *father, MEDCouplingIMesh *mesh):MEDCouplingCartesianAMRMeshGen(mesh),_father(father)
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshSub(MEDCouplingCartesianAMRMeshGen *father, MEDCouplingIMesh *mesh) constructor : empty father !");
}

const MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRMeshSub::getFather() const
{
  return _father;
}

const MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRMeshSub::getGodFather() const
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshSub::getGodFather : Impossible to find a god father because there is a hole in chain !");
  return _father->getGodFather();
}

/*!
 * This method returns the level of \a this. 0 for god father. 1 for children of god father ...
 */
int MEDCouplingCartesianAMRMeshSub::getAbsoluteLevel() const
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshSub::getAbsoluteLevel : Impossible to find a god father because there is a hole in chain !");
  return _father->getAbsoluteLevel()+1;
}

void MEDCouplingCartesianAMRMeshSub::detachFromFather()
{
  _father=0;
  declareAsNew();
}

std::vector< std::pair<int,int> > MEDCouplingCartesianAMRMeshSub::positionRelativeToGodFather(std::vector<int>& st) const
{
  st=_father->getFactors();
  std::size_t dim(st.size());
  std::vector<int> prev(st);
  int id(_father->getPatchIdFromChildMesh(this));
  const MEDCouplingCartesianAMRPatch *p(_father->getPatch(id));
  std::vector< std::pair<int,int> > ret(p->getBLTRRange());
  std::vector<int> delta(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(ret)),start(dim);
  std::transform(delta.begin(),delta.end(),prev.begin(),delta.begin(),std::multiplies<int>());
  for(std::size_t i=0;i<dim;i++)
    start[i]=ret[i].first;
  std::transform(start.begin(),start.end(),prev.begin(),start.begin(),std::multiplies<int>());
  const MEDCouplingCartesianAMRMeshGen *it(_father);
  while(!dynamic_cast<const MEDCouplingCartesianAMRMesh *>(it))
    {
      const MEDCouplingCartesianAMRMeshSub *itc(static_cast<const MEDCouplingCartesianAMRMeshSub *>(it));
      int id2(itc->_father->getPatchIdFromChildMesh(itc));
      const MEDCouplingCartesianAMRPatch *p2(itc->_father->getPatch(id2));
      const std::vector< std::pair<int,int> >& start2(p2->getBLTRRange());
      std::vector<int> tmp(dim);
      for(std::size_t i=0;i<dim;i++)
        tmp[i]=start2[i].first;
      //
      prev=itc->_father->getFactors();
      std::transform(st.begin(),st.end(),prev.begin(),st.begin(),std::multiplies<int>());
      std::transform(st.begin(),st.end(),tmp.begin(),tmp.begin(),std::multiplies<int>());
      std::transform(start.begin(),start.end(),tmp.begin(),start.begin(),std::plus<int>());
      it=itc->_father;
    }
  for(std::size_t i=0;i<dim;i++)
    {
      ret[i].first=start[i];
      ret[i].second=start[i]+delta[i];
    }
  return ret;
}

int MEDCouplingCartesianAMRMeshSub::getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const
{
  if(this==ref)
    return 0;
  if(_father==0)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshSub::getAbsoluteLevelRelativeTo : ref is not in the progeny of this !");
  else
    return _father->getAbsoluteLevelRelativeTo(ref)+1;
}

MEDCouplingCartesianAMRMeshSub::MEDCouplingCartesianAMRMeshSub(const MEDCouplingCartesianAMRMeshSub& other, MEDCouplingCartesianAMRMeshGen *father):MEDCouplingCartesianAMRMeshGen(other),_father(father)
{
}

MEDCouplingCartesianAMRMeshSub *MEDCouplingCartesianAMRMeshSub::deepCopy(MEDCouplingCartesianAMRMeshGen *fath) const
{
  return new MEDCouplingCartesianAMRMeshSub(*this,fath);
}

/*!
 * \sa getPositionRelativeTo
 */
void MEDCouplingCartesianAMRMeshSub::getPositionRelativeToInternal(const MEDCouplingCartesianAMRMeshGen *ref, std::vector<int>& ret) const
{
  if(this==ref)
    return ;
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMeshSub::getPositionRelativeToInternal : ref is not in the progeny of this !");
  int myId(_father->getPatchIdFromChildMesh(this));
  ret.push_back(myId);
  _father->getPositionRelativeToInternal(ref,ret);
}

MEDCouplingCartesianAMRMesh *MEDCouplingCartesianAMRMesh::New(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                              const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop)
{
  return new MEDCouplingCartesianAMRMesh(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop);
}

MEDCouplingCartesianAMRMesh *MEDCouplingCartesianAMRMesh::New(MEDCouplingIMesh *mesh)
{
  return new MEDCouplingCartesianAMRMesh(mesh);
}

const MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRMesh::getFather() const
{
  //I'm god father ! No father !
  return 0;
}

const MEDCouplingCartesianAMRMeshGen *MEDCouplingCartesianAMRMesh::getGodFather() const
{
  return this;
}

int MEDCouplingCartesianAMRMesh::getAbsoluteLevel() const
{
  return 0;
}

void MEDCouplingCartesianAMRMesh::detachFromFather()
{//not a bug - do nothing
}

std::vector< std::pair<int,int> > MEDCouplingCartesianAMRMesh::positionRelativeToGodFather(std::vector<int>& st) const
{
  st=_mesh->getCellGridStructure();
  return MEDCouplingStructuredMesh::GetCompactFrmtFromDimensions(st);
}

int MEDCouplingCartesianAMRMesh::getAbsoluteLevelRelativeTo(const MEDCouplingCartesianAMRMeshGen *ref) const
{
  if(this==ref)
    return 0;
  throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::getAbsoluteLevelRelativeTo : ref is not in the progeny of this !");
}

std::vector<MEDCouplingCartesianAMRPatchGen *> MEDCouplingCartesianAMRMesh::retrieveGridsAt(int absoluteLev) const
{
  std::vector< MCAuto<MEDCouplingCartesianAMRPatchGen> > rets;
  retrieveGridsAtInternal(absoluteLev,rets);
  std::vector< MEDCouplingCartesianAMRPatchGen * > ret(rets.size());
  for(std::size_t i=0;i<rets.size();i++)
    {
      ret[i]=rets[i].retn();
    }
  return ret;
}

MEDCouplingCartesianAMRMesh *MEDCouplingCartesianAMRMesh::deepCopy(MEDCouplingCartesianAMRMeshGen *father) const
{
  if(father)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::deepCopy : specifying a not null father for a God Father object !");
  return new MEDCouplingCartesianAMRMesh(*this);
}

/*!
 * This method creates a multi level patches split at once.
 * This method calls as times as size of \a bso createPatchesFromCriterion. Size of \a bso and size of \a factors must be the same !
 * \b WARNING, after the call the number of levels in \a this is equal to bso.size() + 1 !
 *
 * \param [in] bso
 * \param [in] criterion
 * \param [in] factors
 * \param [in] eps - See DataArrayDouble::toVectorOfBool for more information about the semantic of eps.
 *
 * \sa createPatchesFromCriterion
 */
void MEDCouplingCartesianAMRMesh::createPatchesFromCriterionML(const std::vector<const INTERP_KERNEL::BoxSplittingOptions *>& bso, const DataArrayDouble *criterion, const std::vector< std::vector<int> >& factors, double eps)
{
  std::size_t nbOfLevs(bso.size());
  if(nbOfLevs!=factors.size())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createPatchesFromCriterionML : size of vectors must be the same !");
  if(nbOfLevs==0)
    return ;
  if(!bso[0])
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createPatchesFromCriterionML : pointers in 1st arg must be not NULL !");
  createPatchesFromCriterion(*bso[0],criterion,factors[0],eps);
  for(std::size_t i=1;i<nbOfLevs;i++)
    {
      if(!bso[i])
        throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createPatchesFromCriterionML : presence of a NULL BoxSplittingOptions in input vector !");
      //
      std::vector<MEDCouplingCartesianAMRPatchGen *> elts(retrieveGridsAt((int)(i)));
      std::size_t sz(elts.size());
      std::vector< MCAuto<MEDCouplingCartesianAMRPatchGen> > elts2(sz);
      std::vector< MCAuto<DataArrayDouble> > elts3(sz);
      for(std::size_t ii=0;ii<sz;ii++)
        elts2[ii]=elts[ii];
      //
      static const char TMP_STR[]="TMP";
      std::vector< std::pair<std::string,int> > fieldNames(1); fieldNames[0].first=TMP_STR; fieldNames[0].second=1;
      MCAuto<MEDCouplingAMRAttribute> att(MEDCouplingAMRAttribute::New(this,fieldNames,0));
      att->alloc();
      DataArrayDouble *tmpDa(const_cast<DataArrayDouble *>(att->getFieldOn(this,TMP_STR)));
      tmpDa->deepCopyFrom(*criterion);
      att->synchronizeCoarseToFine();
      for(std::size_t ii=0;ii<sz;ii++)
        {
          const DataArrayDouble *critOnLeaf(att->getFieldOn(const_cast<MEDCouplingCartesianAMRMeshGen *>(elts[ii]->getMesh()),TMP_STR));
          elts3[ii]=const_cast<DataArrayDouble *>(critOnLeaf); elts3[ii]->incrRef();
        }
      att=0;
      for(std::size_t ii=0;ii<sz;ii++)
        const_cast<MEDCouplingCartesianAMRMeshGen *>(elts[ii]->getMesh())->createPatchesFromCriterion(*bso[i],elts3[ii],factors[i],eps);
    }
}

void MEDCouplingCartesianAMRMesh::getPositionRelativeToInternal(const MEDCouplingCartesianAMRMeshGen *ref, std::vector<int>& ret) const
{

}

MEDCouplingCartesianAMRMesh::MEDCouplingCartesianAMRMesh(const MEDCouplingCartesianAMRMesh& other):MEDCouplingCartesianAMRMeshGen(other)
{
}

MEDCouplingCartesianAMRMesh::MEDCouplingCartesianAMRMesh(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                         const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop):MEDCouplingCartesianAMRMeshGen(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop)
{
}

MEDCouplingCartesianAMRMesh::MEDCouplingCartesianAMRMesh(MEDCouplingIMesh *mesh):MEDCouplingCartesianAMRMeshGen(mesh)
{
}

std::vector<const BigMemoryObject *> MEDCouplingCartesianAMRMesh::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDCouplingCartesianAMRMeshGen::getDirectChildrenWithNull());
  return ret;
}

MEDCouplingCartesianAMRMesh::~MEDCouplingCartesianAMRMesh()
{
}
