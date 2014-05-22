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
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingUMesh.hxx"

#include <limits>
#include <sstream>
#include <numeric>

using namespace ParaMEDMEM;

/// @cond INTERNAL

/*!
 * \param [in] mesh not null pointer of refined mesh replacing the cell range of \a father defined by the bottom left and top right just after.
 * \param [in] bottomLeftTopRight a vector equal to the space dimension of \a mesh that specifies for each dimension, the included cell start of the range for the first element of the pair,
 *                                a the end cell (\b excluded) of the range for the second element of the pair.
 * \param [in] factors The refinement per axis relative to the father of \a this.
 */
MEDCouplingCartesianAMRPatch::MEDCouplingCartesianAMRPatch(MEDCouplingCartesianAMRMesh *mesh, const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("EDCouplingCartesianAMRPatch constructor : input mesh is NULL !");
  _mesh=mesh; _mesh->incrRef();
  int dim((int)bottomLeftTopRight.size()),dimExp(_mesh->getSpaceDimension());
  if(dim!=dimExp)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch constructor : space dimension of father and input bottomLeft/topRight size mismatches !");
  _bl_tr=bottomLeftTopRight;
  if((int)factors.size()!=dimExp)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRPatch constructor : space dimension of father and input factors per axis size mismatches !");
  _factors=factors;
}

int MEDCouplingCartesianAMRPatch::getNumberOfCellsRecursiveWithOverlap() const
{
  return _mesh->getNumberOfCellsRecursiveWithOverlap();
}

int MEDCouplingCartesianAMRPatch::getNumberOfCellsRecursiveWithoutOverlap() const
{
  return _mesh->getNumberOfCellsRecursiveWithoutOverlap();
}

int MEDCouplingCartesianAMRPatch::getMaxNumberOfLevelsRelativeToThis() const
{
  return _mesh->getMaxNumberOfLevelsRelativeToThis();
}

void MEDCouplingCartesianAMRPatch::addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors)
{
  return _mesh->addPatch(bottomLeftTopRight,factors);
}

int MEDCouplingCartesianAMRPatch::getNumberOfOverlapedCellsForFather() const
{
  return MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt(_bl_tr);
}

std::size_t MEDCouplingCartesianAMRPatch::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(sizeof(MEDCouplingCartesianAMRPatch));
  ret+=_bl_tr.capacity()*sizeof(std::pair<int,int>);
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingCartesianAMRPatch::getDirectChildren() const
{
  std::vector<const BigMemoryObject *> ret;
  if((const MEDCouplingCartesianAMRMesh *)_mesh)
    ret.push_back((const MEDCouplingCartesianAMRMesh *)_mesh);
  return ret;
}

/// @endcond


MEDCouplingCartesianAMRMesh *MEDCouplingCartesianAMRMesh::New(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                              const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop)
{
  return new MEDCouplingCartesianAMRMesh(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop);
}

int MEDCouplingCartesianAMRMesh::getSpaceDimension() const
{
  return _mesh->getSpaceDimension();
}

int MEDCouplingCartesianAMRMesh::getMaxNumberOfLevelsRelativeToThis() const
{
  int ret(1);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    ret=std::max(ret,(*it)->getMaxNumberOfLevelsRelativeToThis()+1);
  return ret;
}

int MEDCouplingCartesianAMRMesh::getNumberOfCellsAtCurrentLevel() const
{
  return _mesh->getNumberOfCells();
}

int MEDCouplingCartesianAMRMesh::getNumberOfCellsRecursiveWithOverlap() const
{
  int ret(_mesh->getNumberOfCells());
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      ret+=(*it)->getNumberOfCellsRecursiveWithOverlap();
    }
  return ret;
}

int MEDCouplingCartesianAMRMesh::getNumberOfCellsRecursiveWithoutOverlap() const
{
  int ret(_mesh->getNumberOfCells());
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      ret-=(*it)->getNumberOfOverlapedCellsForFather();
      ret+=(*it)->getNumberOfCellsRecursiveWithoutOverlap();
    }
  return ret;
}

const MEDCouplingCartesianAMRMesh *MEDCouplingCartesianAMRMesh::getFather() const
{
  return _father;
}

const MEDCouplingCartesianAMRMesh *MEDCouplingCartesianAMRMesh::getGodFather() const
{
  if(_father==0)
    return this;
  else
    return _father->getGodFather();
}

void MEDCouplingCartesianAMRMesh::detachFromFather()
{
  _father=0;
}

/*!
 * \param [in] bottomLeftTopRight a vector equal to the space dimension of \a mesh that specifies for each dimension, the included cell start of the range for the first element of the pair,
 *                                a the end cell (\b excluded) of the range for the second element of the pair.
 * \param [in] factors The != 0 factor of refinement per axis.
 */
void MEDCouplingCartesianAMRMesh::addPatch(const std::vector< std::pair<int,int> >& bottomLeftTopRight, const std::vector<int>& factors)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingIMesh> mesh(static_cast<MEDCouplingIMesh *>(_mesh->buildStructuredSubPart(bottomLeftTopRight)));
  mesh->refineWithFactor(factors);
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRMesh> zeMesh(new MEDCouplingCartesianAMRMesh(this,mesh));
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> elt(new MEDCouplingCartesianAMRPatch(zeMesh,bottomLeftTopRight,factors));
  _patches.push_back(elt);
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
      if (hole.size()>0)
        {
          double center(((double)len/2.));
          for(std::size_t i=0;i<hole.size();i++)
            distance.push_back(fabs(hole[i]+1.-center));

          double distanceMin=*std::min_element(distance.begin(),distance.end());
          int posDistanceMin=std::find(distance.begin(),distance.end(),distanceMin)-distance.begin()-1;
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
void MEDCouplingCartesianAMRMesh::createPatchesFromCriterion(const INTERP_KERNEL::BoxSplittingOptions& bso, const DataArrayByte *criterion, const std::vector<int>& factors)
{
  if(!criterion || !criterion->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createPatchesFromCriterion : the criterion DataArrayByte instance must be allocated and not NULL !");
  int nbCells(getNumberOfCellsAtCurrentLevel());
  if(nbCells!=criterion->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh::createPatchesFromCriterion : the number of tuples of criterion array must be equal to the number of cells at the current level !");
  _patches.clear();
  //
  std::vector<int> cgs(_mesh->getCellGridStructure());
  std::vector<bool> crit(criterion->toVectorOfBool());//check that criterion has one component.
  std::vector< MEDCouplingAutoRefCountObjectPtr<InternalPatch> > listOfPatches,listOfPatchesOK;
  //
  MEDCouplingAutoRefCountObjectPtr<InternalPatch> p(new InternalPatch);
  p->setNumberOfTrue(MEDCouplingStructuredMesh::FindMinimalPartOf(cgs,crit,p->getCriterion(),p->getPart()));
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
          std::cerr << axisId << std::endl;
          if(cutFound)
            { DealWithCut(*it,axisId,cutPlace,listOfPatchesTmp); continue; }//action 4
          listOfPatchesOK.push_back(DealWithNoCut(*it));
        }
      listOfPatches=listOfPatchesTmp;
    }
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<InternalPatch> >::const_iterator it=listOfPatchesOK.begin();it!=listOfPatchesOK.end();it++)
    addPatch((*it)->getConstPart(),factors);
}

void MEDCouplingCartesianAMRMesh::removePatch(int patchId)
{
  checkPatchId(patchId);
  int sz((int)_patches.size()),j(0);
  std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> > patches(sz-1);
  for(int i=0;i<sz;i++)
    if(i!=patchId)
      patches[j++]=_patches[i];
  (const_cast<MEDCouplingCartesianAMRMesh *>(_patches[patchId]->getMesh()))->detachFromFather();
  _patches=patches;
  declareAsNew();
}

int MEDCouplingCartesianAMRMesh::getNumberOfPatches() const
{
  return (int)_patches.size();
}

const MEDCouplingCartesianAMRPatch *MEDCouplingCartesianAMRMesh::getPatch(int patchId) const
{
  checkPatchId(patchId);
  return _patches[patchId];
}

MEDCouplingUMesh *MEDCouplingCartesianAMRMesh::buildUnstructured() const
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
MEDCoupling1SGTUMesh *MEDCouplingCartesianAMRMesh::buildMeshFromPatchEnvelop() const
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

MEDCouplingCartesianAMRMesh::MEDCouplingCartesianAMRMesh(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                                         const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop):_father(0)
{
  _mesh=MEDCouplingIMesh::New(meshName,spaceDim,nodeStrctStart,nodeStrctStop,originStart,originStop,dxyzStart,dxyzStop);
}

MEDCouplingCartesianAMRMesh::MEDCouplingCartesianAMRMesh(MEDCouplingCartesianAMRMesh *father, MEDCouplingIMesh *mesh):_father(father)
{
  if(!_father)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh(MEDCouplingIMesh *mesh) constructor : empty father !");
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingCartesianAMRMesh(MEDCouplingIMesh *mesh) constructor : The input mesh is null !");
  mesh->checkCoherency();
  _mesh=mesh; _mesh->incrRef();
}

void MEDCouplingCartesianAMRMesh::checkPatchId(int patchId) const
{
  int sz(getNumberOfPatches());
  if(patchId<0 || patchId>=sz)
    {
      std::ostringstream oss; oss << "MEDCouplingCartesianAMRMesh::checkPatchId : invalid patchId (" << patchId << ") ! Must be in [0," << sz << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

std::size_t MEDCouplingCartesianAMRMesh::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(MEDCouplingCartesianAMRMesh);
}

std::vector<const BigMemoryObject *> MEDCouplingCartesianAMRMesh::getDirectChildren() const
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

void MEDCouplingCartesianAMRMesh::updateTime() const
{
  if((const MEDCouplingIMesh *)_mesh)
    updateTimeWith(*_mesh);
  for(std::vector< MEDCouplingAutoRefCountObjectPtr<MEDCouplingCartesianAMRPatch> >::const_iterator it=_patches.begin();it!=_patches.end();it++)
    {
      const MEDCouplingCartesianAMRPatch *elt(*it);
      if(!elt)
        continue;
      const MEDCouplingCartesianAMRMesh *mesh(elt->getMesh());
      if(mesh)
        updateTimeWith(*mesh);
    }
}
