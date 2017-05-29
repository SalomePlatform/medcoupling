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
// Author : Anthony Geay (EDF R&D)

#include "MEDCouplingIMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <functional>
#include <algorithm>
#include <sstream>
#include <numeric>

using namespace MEDCoupling;

MEDCouplingIMesh::MEDCouplingIMesh():_space_dim(-1)
{
  _origin[0]=0.; _origin[1]=0.; _origin[2]=0.;
  _dxyz[0]=0.; _dxyz[1]=0.; _dxyz[2]=0.;
  _structure[0]=0; _structure[1]=0; _structure[2]=0;
}

MEDCouplingIMesh::MEDCouplingIMesh(const MEDCouplingIMesh& other, bool deepCpy):MEDCouplingStructuredMesh(other,deepCpy),_space_dim(other._space_dim),_axis_unit(other._axis_unit)
{
  _origin[0]=other._origin[0]; _origin[1]=other._origin[1]; _origin[2]=other._origin[2];
  _dxyz[0]=other._dxyz[0]; _dxyz[1]=other._dxyz[1]; _dxyz[2]=other._dxyz[2];
  _structure[0]=other._structure[0]; _structure[1]=other._structure[1]; _structure[2]=other._structure[2];
}

MEDCouplingIMesh::~MEDCouplingIMesh()
{
}

MEDCouplingIMesh *MEDCouplingIMesh::New()
{
  return new MEDCouplingIMesh;
}

MEDCouplingIMesh *MEDCouplingIMesh::New(const std::string& meshName, int spaceDim, const int *nodeStrctStart, const int *nodeStrctStop,
                                        const double *originStart, const double *originStop, const double *dxyzStart, const double *dxyzStop)
{
  MCAuto<MEDCouplingIMesh> ret(new MEDCouplingIMesh);
  ret->setName(meshName);
  ret->setSpaceDimension(spaceDim);
  ret->setNodeStruct(nodeStrctStart,nodeStrctStop);
  ret->setOrigin(originStart,originStop);
  ret->setDXYZ(dxyzStart,dxyzStop);
  return ret.retn();
}

MEDCouplingIMesh *MEDCouplingIMesh::deepCopy() const
{
  return clone(true);
}

MEDCouplingIMesh *MEDCouplingIMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingIMesh(*this,recDeepCpy);
}

const DataArrayDouble *MEDCouplingIMesh::getDirectAccessOfCoordsArrIfInStructure() const
{
  throw INTERP_KERNEL::Exception("MEDCouplingIMesh::getDirectAccessOfCoordsArrIfInStructure : MEDCouplingIMesh does not aggregate array of coordinates !");
}

/*!
 * This method creates a copy of \a this enlarged by \a ghostLev cells on each axis.
 * If \a ghostLev equal to 0 this method behaves as MEDCouplingIMesh::clone.
 *
 * \param [in] ghostLev - the ghost level expected
 * \return MEDCouplingIMesh * - a newly alloacted object to be managed by the caller.
 * \throw if \a ghostLev < 0.
 */
MEDCouplingIMesh *MEDCouplingIMesh::buildWithGhost(int ghostLev) const
{
  if(ghostLev<0)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::buildWithGhost : the ghostLev must be >= 0 !");
  checkConsistencyLight();
  int spaceDim(getSpaceDimension());
  double origin[3],dxyz[3];
  int structure[3];
  for(int i=0;i<spaceDim;i++)
    {
      origin[i]=_origin[i]-ghostLev*_dxyz[i];
      dxyz[i]=_dxyz[i];
      structure[i]=_structure[i]+2*ghostLev;
    }
  MCAuto<MEDCouplingIMesh> ret(MEDCouplingIMesh::New(getName(),spaceDim,structure,structure+spaceDim,origin,origin+spaceDim,dxyz,dxyz+spaceDim));
  ret->copyTinyInfoFrom(this);
  return ret.retn();
}

void MEDCouplingIMesh::setNodeStruct(const int *nodeStrctStart, const int *nodeStrctStop)
{
  checkSpaceDimension();
  int sz((int)std::distance(nodeStrctStart,nodeStrctStop));
  if(sz!=_space_dim)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::setNodeStruct : input vector of node structure has not the right size ! Or change space dimension before calling it !");
  std::copy(nodeStrctStart,nodeStrctStop,_structure);
  declareAsNew();
}

std::vector<int> MEDCouplingIMesh::getNodeStruct() const
{
  checkSpaceDimension();
  return std::vector<int>(_structure,_structure+_space_dim);
}

void MEDCouplingIMesh::setOrigin(const double *originStart, const double *originStop)
{
  checkSpaceDimension();
  int sz((int)std::distance(originStart,originStop));
  if(sz!=_space_dim)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::setOrigin : input vector of origin vector has not the right size ! Or change space dimension before calling it !");
  std::copy(originStart,originStop,_origin);
  declareAsNew();
}

std::vector<double> MEDCouplingIMesh::getOrigin() const
{
  checkSpaceDimension();
  return std::vector<double>(_origin,_origin+_space_dim);
}

void MEDCouplingIMesh::setDXYZ(const double *dxyzStart, const double *dxyzStop)
{
  checkSpaceDimension();
  int sz((int)std::distance(dxyzStart,dxyzStop));
  if(sz!=_space_dim)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::setDXYZ : input vector of dxyz vector has not the right size ! Or change space dimension before calling it !");
  std::copy(dxyzStart,dxyzStop,_dxyz);
  declareAsNew();
}

std::vector<double> MEDCouplingIMesh::getDXYZ() const
{
  checkSpaceDimension();
  return std::vector<double>(_dxyz,_dxyz+_space_dim);
}

void MEDCouplingIMesh::setAxisUnit(const std::string& unitName)
{
  _axis_unit=unitName;
  declareAsNew();
}

std::string MEDCouplingIMesh::getAxisUnit() const
{
  return _axis_unit;
}

/*!
 * This method returns the measure of any cell in \a this.
 * This specific method of image grid mesh utilizes the fact that any cell in \a this have the same measure.
 * The value returned by this method is those used to feed the returned field in the MEDCouplingIMesh::getMeasureField.
 *
 * \sa getMeasureField
 */
double MEDCouplingIMesh::getMeasureOfAnyCell() const
{
  checkConsistencyLight();
  int dim(getSpaceDimension());
  double ret(1.);
  for(int i=0;i<dim;i++)
    ret*=fabs(_dxyz[i]);
  return ret;
}

/*!
 * This method is allows to convert \a this into MEDCouplingCMesh instance.
 * This method is the middle level between MEDCouplingIMesh and the most general MEDCouplingUMesh.
 * This method is useful for MED writers that do not have still the image grid support.
 *
 * \sa MEDCouplingMesh::buildUnstructured
 */
MEDCouplingCMesh *MEDCouplingIMesh::convertToCartesian() const
{
  checkConsistencyLight();
  MCAuto<MEDCouplingCMesh> ret(MEDCouplingCMesh::New());
  try
  { ret->copyTinyInfoFrom(this); }
  catch(INTERP_KERNEL::Exception& ) { }
  int spaceDim(getSpaceDimension());
  std::vector<std::string> infos(buildInfoOnComponents());
  for(int i=0;i<spaceDim;i++)
    {
      MCAuto<DataArrayDouble> arr(DataArrayDouble::New()); arr->alloc(_structure[i],1); arr->setInfoOnComponent(0,infos[i]);
      arr->iota(); arr->applyLin(_dxyz[i],_origin[i]);
      ret->setCoordsAt(i,arr);
    }
  return ret.retn();
}

/*!
 * This method refines \a this uniformaly along all of its dimensions. In case of success the space covered by \a this will remain
 * the same before the invocation except that the number of cells will be multiplied by \a factor ^ this->getMeshDimension().
 * The origin of \a this will be not touched only spacing and node structure will be changed.
 * This method can be useful for AMR users.
 */
void MEDCouplingIMesh::refineWithFactor(const std::vector<int>& factors)
{
  if((int)factors.size()!=_space_dim)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::refineWithFactor : refinement factors must have size equal to spaceDim !");
  checkConsistencyLight();
  std::vector<int> structure(_structure,_structure+3);
  std::vector<double> dxyz(_dxyz,_dxyz+3);
  for(int i=0;i<_space_dim;i++)
    {
      if(factors[i]<=0)
        {
          std::ostringstream oss; oss << "MEDCouplingIMesh::refineWithFactor : factor for axis #" << i << " (" << factors[i] << ")is invalid ! Must be > 0 !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      int factAbs(std::abs(factors[i]));
      double fact2(1./(double)factors[i]);
      structure[i]=(_structure[i]-1)*factAbs+1;
      dxyz[i]=fact2*_dxyz[i];
    }
  std::copy(structure.begin(),structure.end(),_structure);
  std::copy(dxyz.begin(),dxyz.end(),_dxyz);
  declareAsNew();
}

/*!
 * This method returns a newly created mesh containing a single cell in it. This returned cell covers exactly the space covered by \a this.
 *
 * \return MEDCouplingIMesh * - A newly created object (to be managed by the caller with decrRef) containing simply one cell.
 *
 * \throw if \a this does not pass the \c checkConsistencyLight test.
 */
MEDCouplingIMesh *MEDCouplingIMesh::asSingleCell() const
{
  checkConsistencyLight();
  int spaceDim(getSpaceDimension()),nodeSt[3];
  double dxyz[3];
  for(int i=0;i<spaceDim;i++)
    {
      if(_structure[i]>=2)
        {
          nodeSt[i]=2;
          dxyz[i]=(_structure[i]-1)*_dxyz[i];
        }
      else
        {
          nodeSt[i]=_structure[i];
          dxyz[i]=_dxyz[i];
        }
    }
  MCAuto<MEDCouplingIMesh> ret(MEDCouplingIMesh::New(getName(),getSpaceDimension(),nodeSt,nodeSt+spaceDim,_origin,_origin+spaceDim,dxyz,dxyz+spaceDim));
  ret->copyTinyInfoFrom(this);
  return ret.retn();
}

/*!
 * This static method is useful to condense field on cells of a MEDCouplingIMesh instance coming from a refinement ( MEDCouplingIMesh::refineWithFactor for example)
 * to a coarse MEDCouplingIMesh instance. So this method can be seen as a specialization in P0P0 conservative interpolation non overlaping from fine image mesh
 * to a coarse image mesh. Only tuples ( deduced from \a fineLocInCoarse ) of \a coarseDA will be modified. Other tuples of \a coarseDA will be let unchanged.
 *
 * \param [in] coarseSt The cell structure of coarse mesh.
 * \param [in] fineDA The DataArray containing the cell field on uniformly refined mesh
 * \param [in] fineLocInCoarse The cell localization of refined mesh into the coarse one.
 * \param [in] facts The refinement coefficient per axis.
 * \param [in,out] coarseDA The DataArrayDouble corresponding to the a cell field of a coarse mesh whose cell structure is defined by \a coarseSt.
 *
 * \sa CondenseFineToCoarseGhost,SpreadCoarseToFine
 */
void MEDCouplingIMesh::CondenseFineToCoarse(const std::vector<int>& coarseSt, const DataArrayDouble *fineDA, const std::vector< std::pair<int,int> >& fineLocInCoarse, const std::vector<int>& facts, DataArrayDouble *coarseDA)
{
  if(coarseSt.size()!=fineLocInCoarse.size() || coarseSt.size()!=facts.size())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarse : All input vectors (dimension) must have the same size !");
  if(!coarseDA || !coarseDA->isAllocated() || !fineDA || !fineDA->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarse : the parameters 1 or 3 are NULL or not allocated !");
  int meshDim((int)coarseSt.size()),nbOfTuplesInCoarseExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(coarseSt)),nbOfTuplesInFineExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt(fineLocInCoarse));
  int nbCompo(fineDA->getNumberOfComponents());
  if((int)coarseDA->getNumberOfComponents()!=nbCompo)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarse : the number of components of fine DA and coarse one mismatches !");
  if(meshDim!=(int)fineLocInCoarse.size() || meshDim!=(int)facts.size())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarse : the size of fineLocInCoarse (4th param) and facts (5th param) must be equal to the sier of coarseSt (2nd param) !");
  if(coarseDA->getNumberOfTuples()!=nbOfTuplesInCoarseExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::CondenseFineToCoarse : Expecting " << nbOfTuplesInCoarseExp << " tuples having " << coarseDA->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbTuplesFine(fineDA->getNumberOfTuples());
  if(nbOfTuplesInFineExp==0)
    {
      if(nbTuplesFine==0)
        return ;
      else
        throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarse : Nothing to condense considering the range specified ! But DataArray is not empty !");
    }
  if(nbTuplesFine%nbOfTuplesInFineExp!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarse : Invalid nb of tuples in fine DataArray regarding its structure !");
  int fact(std::accumulate(facts.begin(),facts.end(),1,std::multiplies<int>()));
  if(nbTuplesFine!=fact*nbOfTuplesInFineExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::CondenseFineToCoarse : Invalid number of tuples ("  << nbTuplesFine << ") of fine dataarray is invalid ! Must be " << fact*nbOfTuplesInFineExp << "!";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  // to improve use jump-iterator. Factorizes with SwitchOnIdsFrom BuildExplicitIdsFrom
  double *outPtr(coarseDA->getPointer());
  const double *inPtr(fineDA->begin());
  //
  std::vector<int> dims(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
  switch(meshDim)
  {
    case 1:
      {
        int offset(fineLocInCoarse[0].first),fact0(facts[0]);
        for(int i=0;i<dims[0];i++)
          {
            double *loc(outPtr+(offset+i)*nbCompo);
            for(int ifact=0;ifact<fact0;ifact++,inPtr+=nbCompo)
              {
                if(ifact!=0)
                  std::transform(inPtr,inPtr+nbCompo,loc,loc,std::plus<double>());
                else
                  std::copy(inPtr,inPtr+nbCompo,loc);
              }
          }
        break;
      }
    case 2:
      {
        int kk(fineLocInCoarse[0].first+coarseSt[0]*fineLocInCoarse[1].first),fact1(facts[1]),fact0(facts[0]);
        for(int j=0;j<dims[1];j++)
          {
            for(int jfact=0;jfact<fact1;jfact++)
              {
                for(int i=0;i<dims[0];i++)
                  {
                    double *loc(outPtr+(kk+i)*nbCompo);
                    for(int ifact=0;ifact<fact0;ifact++,inPtr+=nbCompo)
                      {
                        if(jfact!=0 || ifact!=0)
                          std::transform(inPtr,inPtr+nbCompo,loc,loc,std::plus<double>());
                        else
                          std::copy(inPtr,inPtr+nbCompo,loc);
                      }
                  }
              }
            kk+=coarseSt[0];
          }
        break;
      }
    case 3:
      {
        int kk(fineLocInCoarse[0].first+coarseSt[0]*fineLocInCoarse[1].first+coarseSt[0]*coarseSt[1]*fineLocInCoarse[2].first),fact2(facts[2]),fact1(facts[1]),fact0(facts[0]);
        for(int k=0;k<dims[2];k++)
          {
            for(int kfact=0;kfact<fact2;kfact++)
              {
                for(int j=0;j<dims[1];j++)
                  {
                    for(int jfact=0;jfact<fact1;jfact++)
                      {
                        for(int i=0;i<dims[0];i++)
                          {
                            double *loc(outPtr+(kk+i+j*coarseSt[0])*nbCompo);
                            for(int ifact=0;ifact<fact0;ifact++,inPtr+=nbCompo)
                              {
                                if(kfact!=0 || jfact!=0 || ifact!=0)
                                  std::transform(inPtr,inPtr+nbCompo,loc,loc,std::plus<double>());
                                else
                                  std::copy(inPtr,inPtr+nbCompo,loc);
                              }
                          }
                      }
                  }
              }
            kk+=coarseSt[0]*coarseSt[1];
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarse : only dimensions 1, 2 and 3 supported !");
  }
}

/*!
 * This static method is useful to condense field on cells of a MEDCouplingIMesh instance coming from a refinement ( MEDCouplingIMesh::refineWithFactor for example)
 * to a coarse MEDCouplingIMesh instance. So this method can be seen as a specialization in P0P0 conservative interpolation non overlaping from fine image mesh
 * to a coarse image mesh. Only tuples ( deduced from \a fineLocInCoarse ) of \a coarseDA will be modified. Other tuples of \a coarseDA will be let unchanged.
 *
 * \param [in] coarseSt The cell structure of coarse mesh.
 * \param [in] fineDA The DataArray containing the cell field on uniformly refined mesh
 * \param [in] fineLocInCoarse The cell localization of refined mesh into the coarse one.
 * \param [in] facts The refinement coefficient per axis.
 * \param [in,out] coarseDA The DataArrayDouble corresponding to the a cell field of a coarse mesh whose cell structure is defined by \a coarseSt.
 * \param [in] ghostSize - The size of the ghost zone. The ghost zone is expected to be the same for all axis and both for coarse and fine meshes.
 *
 * \sa CondenseFineToCoarse,SpreadCoarseToFineGhost
 */
void MEDCouplingIMesh::CondenseFineToCoarseGhost(const std::vector<int>& coarseSt, const DataArrayDouble *fineDA, const std::vector< std::pair<int,int> >& fineLocInCoarse, const std::vector<int>& facts, DataArrayDouble *coarseDA, int ghostSize)
{
  if(ghostSize<0)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarseGhost : ghost level has to be >= 0 !");
  if(coarseSt.size()!=fineLocInCoarse.size() || coarseSt.size()!=facts.size())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarseGhost : All input vectors (dimension) must have the same size !");
  if(!coarseDA || !coarseDA->isAllocated() || !fineDA || !fineDA->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarseGhost : the parameters 1 or 3 are NULL or not allocated !");
  std::vector<int> coarseStG(coarseSt.size()); std::transform(coarseSt.begin(),coarseSt.end(),coarseStG.begin(),std::bind2nd(std::plus<int>(),2*ghostSize));
  int meshDim((int)coarseSt.size()),nbOfTuplesInCoarseExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(coarseStG));
  int nbCompo(fineDA->getNumberOfComponents());
  if((int)coarseDA->getNumberOfComponents()!=nbCompo)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarseGhost : the number of components of fine DA and coarse one mismatches !");
  if(meshDim!=(int)fineLocInCoarse.size() || meshDim!=(int)facts.size())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarseGhost : the size of fineLocInCoarse (4th param) and facts (5th param) must be equal to the sier of coarseSt (2nd param) !");
  if(coarseDA->getNumberOfTuples()!=nbOfTuplesInCoarseExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::CondenseFineToCoarseGhost : Expecting " << nbOfTuplesInCoarseExp << " tuples in coarse DataArray having " << coarseDA->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  std::vector<int> fineStG(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
  std::transform(fineStG.begin(),fineStG.end(),facts.begin(),fineStG.begin(),std::multiplies<int>());
  std::transform(fineStG.begin(),fineStG.end(),fineStG.begin(),std::bind2nd(std::plus<int>(),2*ghostSize));
  int nbTuplesFine(fineDA->getNumberOfTuples()),nbTuplesFineExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(fineStG));
  if(fineDA->getNumberOfTuples()!=nbTuplesFineExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::CondenseFineToCoarseGhost : Expecting " << nbTuplesFineExp << " tuples in fine DataArray having " << nbTuplesFine << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  double *outPtr(coarseDA->getPointer());
  const double *inPtr(fineDA->begin());
  //
  std::vector<int> dims(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
  switch(meshDim)
  {
    case 1:
      {
        int offset(fineLocInCoarse[0].first+ghostSize),fact0(facts[0]);
        inPtr+=ghostSize*nbCompo;
        for(int i=0;i<dims[0];i++)
          {
            double *loc(outPtr+(offset+i)*nbCompo);
            for(int ifact=0;ifact<fact0;ifact++,inPtr+=nbCompo)
              {
                if(ifact!=0)
                  std::transform(inPtr,inPtr+nbCompo,loc,loc,std::plus<double>());
                else
                  std::copy(inPtr,inPtr+nbCompo,loc);
              }
          }
        break;
      }
    case 2:
      {
        int nxwg(coarseSt[0]+2*ghostSize);
        int kk(fineLocInCoarse[0].first+ghostSize+nxwg*(fineLocInCoarse[1].first+ghostSize)),fact1(facts[1]),fact0(facts[0]);
        inPtr+=(dims[0]*fact0+2*ghostSize)*ghostSize*nbCompo;
        for(int j=0;j<dims[1];j++)
          {
             for(int jfact=0;jfact<fact1;jfact++)
              {
                inPtr+=ghostSize*nbCompo;
                for(int i=0;i<dims[0];i++)
                  {
                    double *loc(outPtr+(kk+i)*nbCompo);
                    for(int ifact=0;ifact<fact0;ifact++,inPtr+=nbCompo)
                      {
                        if(jfact!=0 || ifact!=0)
                          std::transform(inPtr,inPtr+nbCompo,loc,loc,std::plus<double>());
                        else
                          std::copy(inPtr,inPtr+nbCompo,loc);
                      }
                  }
                inPtr+=ghostSize*nbCompo;
              }
            kk+=nxwg;
          }
        break;
      }
    case 3:
      {
        int nxwg(coarseSt[0]+2*ghostSize),nxywg((coarseSt[0]+2*ghostSize)*(coarseSt[1]+2*ghostSize));
        int kk(fineLocInCoarse[0].first+ghostSize+nxwg*(fineLocInCoarse[1].first+ghostSize)+nxywg*(fineLocInCoarse[2].first+ghostSize)),fact2(facts[2]),fact1(facts[1]),fact0(facts[0]);
        inPtr+=(dims[0]*fact0+2*ghostSize)*(dims[1]*fact1+2*ghostSize)*ghostSize*nbCompo;
        for(int k=0;k<dims[2];k++)
          {
            for(int kfact=0;kfact<fact2;kfact++)
              {
                inPtr+=ghostSize*(dims[0]*fact0+2*ghostSize)*nbCompo;
                for(int j=0;j<dims[1];j++)
                  {
                    int kky(j*nxwg);
                    for(int jfact=0;jfact<fact1;jfact++)
                      {
                        inPtr+=ghostSize*nbCompo;
                        for(int i=0;i<dims[0];i++)
                          {
                            double *loc(outPtr+(kky+kk+i)*nbCompo);
                            for(int ifact=0;ifact<fact0;ifact++,inPtr+=nbCompo)
                              {
                                if(kfact!=0 || jfact!=0 || ifact!=0)
                                  std::transform(inPtr,inPtr+nbCompo,loc,loc,std::plus<double>());
                                else
                                  std::copy(inPtr,inPtr+nbCompo,loc);
                              }
                          }
                        inPtr+=ghostSize*nbCompo;
                      }
                  }
                inPtr+=ghostSize*(dims[0]*fact0+2*ghostSize)*nbCompo;
              }
            kk+=nxywg;
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CondenseFineToCoarseGhost : only dimensions 1, 2, 3 supported !");
  }
}

/*!
 * This method spreads the values of coarse data \a coarseDA into \a fineDA.
 *
 * \param [in] coarseDA The DataArrayDouble corresponding to the a cell field of a coarse mesh whose cell structure is defined by \a coarseSt.
 * \param [in] coarseSt The cell structure of coarse mesh.
 * \param [in,out] fineDA The DataArray containing the cell field on uniformly refined mesh
 * \param [in] fineLocInCoarse The cell localization of refined mesh into the coarse one.
 * \param [in] facts The refinement coefficient per axis.
 * \sa SpreadCoarseToFineGhost, CondenseFineToCoarse
 */
void MEDCouplingIMesh::SpreadCoarseToFine(const DataArrayDouble *coarseDA, const std::vector<int>& coarseSt, DataArrayDouble *fineDA, const std::vector< std::pair<int,int> >& fineLocInCoarse, const std::vector<int>& facts)
{
  if(coarseSt.size()!=fineLocInCoarse.size() || coarseSt.size()!=facts.size())
      throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFine : All input vectors (dimension) must have the same size !");
  if(!coarseDA || !coarseDA->isAllocated() || !fineDA || !fineDA->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFine : the parameters 1 or 3 are NULL or not allocated !");
  int meshDim((int)coarseSt.size()),nbOfTuplesInCoarseExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(coarseSt)),nbOfTuplesInFineExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenRangeInCompactFrmt(fineLocInCoarse));
  int nbCompo(fineDA->getNumberOfComponents());
  if((int)coarseDA->getNumberOfComponents()!=nbCompo)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFine : the number of components of fine DA and coarse one mismatches !");
  if(meshDim!=(int)fineLocInCoarse.size() || meshDim!=(int)facts.size())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFine : the size of fineLocInCoarse (4th param) and facts (5th param) must be equal to the sier of coarseSt (2nd param) !");
  if(coarseDA->getNumberOfTuples()!=nbOfTuplesInCoarseExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::SpreadCoarseToFine : Expecting " << nbOfTuplesInCoarseExp << " tuples having " << coarseDA->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  int nbTuplesFine(fineDA->getNumberOfTuples());
  if(nbTuplesFine%nbOfTuplesInFineExp!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFine : Invalid nb of tuples in fine DataArray regarding its structure !");
  int fact(std::accumulate(facts.begin(),facts.end(),1,std::multiplies<int>()));
  if(nbTuplesFine!=fact*nbOfTuplesInFineExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::SpreadCoarseToFine : Invalid number of tuples ("  << nbTuplesFine << ") of fine dataarray is invalid ! Must be " << fact*nbOfTuplesInFineExp << "!";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  // to improve use jump-iterator. Factorizes with SwitchOnIdsFrom BuildExplicitIdsFrom
  double *outPtr(fineDA->getPointer());
  const double *inPtr(coarseDA->begin());
  //
  std::vector<int> dims(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
  switch(meshDim)
  {
    case 1:
      {
        int offset(fineLocInCoarse[0].first),fact0(facts[0]);
        for(int i=0;i<dims[0];i++)
          {
            const double *loc(inPtr+(offset+i)*nbCompo);
            for(int ifact=0;ifact<fact0;ifact++)
              outPtr=std::copy(loc,loc+nbCompo,outPtr);
          }
        break;
      }
    case 2:
      {
        int kk(fineLocInCoarse[0].first+coarseSt[0]*fineLocInCoarse[1].first),fact0(facts[0]),fact1(facts[1]);
        for(int j=0;j<dims[1];j++)
          {
            for(int jfact=0;jfact<fact1;jfact++)
              {
                for(int i=0;i<dims[0];i++)
                  {
                    const double *loc(inPtr+(kk+i)*nbCompo);
                    for(int ifact=0;ifact<fact0;ifact++)
                      outPtr=std::copy(loc,loc+nbCompo,outPtr);
                  }
              }
            kk+=coarseSt[0];
          }
        break;
      }
    case 3:
      {
        int kk(fineLocInCoarse[0].first+coarseSt[0]*fineLocInCoarse[1].first+coarseSt[0]*coarseSt[1]*fineLocInCoarse[2].first),fact0(facts[0]),fact1(facts[2]),fact2(facts[2]);
        for(int k=0;k<dims[2];k++)
          {
            for(int kfact=0;kfact<fact2;kfact++)
              {
                for(int j=0;j<dims[1];j++)
                  {
                    for(int jfact=0;jfact<fact1;jfact++)
                      {
                        for(int i=0;i<dims[0];i++)
                          {
                            const double *loc(inPtr+(kk+i+j*coarseSt[0])*nbCompo);
                            for(int ifact=0;ifact<fact0;ifact++)
                              outPtr=std::copy(loc,loc+nbCompo,outPtr);
                          }
                      }
                  }
              }
            kk+=coarseSt[0]*coarseSt[1];
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFine : only dimensions 1, 2 and 3 supported !");
  }
}

/*!
 * This method spreads the values of coarse data \a coarseDA into \a fineDA.
 *
 * \param [in] coarseDA The DataArrayDouble corresponding to the a cell field of a coarse mesh whose cell structure is defined by \a coarseSt.
 * \param [in] coarseSt The cell structure of coarse mesh.
 * \param [in,out] fineDA The DataArray containing the cell field on uniformly refined mesh
 * \param [in] fineLocInCoarse The cell localization of refined mesh into the coarse one.
 * \param [in] facts The refinement coefficient per axis.
 * \param [in] ghostSize - The size of the ghost zone. The ghost zone is expected to be the same for all axis and both for coarse and fine meshes.
 *
 * \sa CondenseFineToCoarse, SpreadCoarseToFineGhostZone
 */
void MEDCouplingIMesh::SpreadCoarseToFineGhost(const DataArrayDouble *coarseDA, const std::vector<int>& coarseSt, DataArrayDouble *fineDA, const std::vector< std::pair<int,int> >& fineLocInCoarse, const std::vector<int>& facts, int ghostSize)
{
  if(ghostSize<0)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhost : ghost level has to be >= 0 !");
  if(coarseSt.size()!=fineLocInCoarse.size() || coarseSt.size()!=facts.size())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhost : All input vectors (dimension) must have the same size !");
  if(!coarseDA || !coarseDA->isAllocated() || !fineDA || !fineDA->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhost : the parameters 1 or 3 are NULL or not allocated !");
  std::vector<int> coarseStG(coarseSt.size()); std::transform(coarseSt.begin(),coarseSt.end(),coarseStG.begin(),std::bind2nd(std::plus<int>(),2*ghostSize));
  int meshDim((int)coarseSt.size()),nbOfTuplesInCoarseExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(coarseStG));
  int nbCompo(fineDA->getNumberOfComponents());
  if((int)coarseDA->getNumberOfComponents()!=nbCompo)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhost : the number of components of fine DA and coarse one mismatches !");
  if(meshDim!=(int)fineLocInCoarse.size() || meshDim!=(int)facts.size())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhost : the size of fineLocInCoarse (4th param) and facts (5th param) must be equal to the sier of coarseSt (2nd param) !");
  if(coarseDA->getNumberOfTuples()!=nbOfTuplesInCoarseExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::SpreadCoarseToFineGhost : Expecting " << nbOfTuplesInCoarseExp << " tuples having " << coarseDA->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  std::vector<int> fineStG(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
  std::transform(fineStG.begin(),fineStG.end(),facts.begin(),fineStG.begin(),std::multiplies<int>());
  std::transform(fineStG.begin(),fineStG.end(),fineStG.begin(),std::bind2nd(std::plus<int>(),2*ghostSize));
  int nbTuplesFine(fineDA->getNumberOfTuples()),nbTuplesFineExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(fineStG));
  if(fineDA->getNumberOfTuples()!=nbTuplesFineExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::SpreadCoarseToFineGhost : Expecting " << nbTuplesFineExp << " tuples in fine DataArray having " << nbTuplesFine << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  double *outPtr(fineDA->getPointer());
  const double *inPtr(coarseDA->begin());
  //
  switch(meshDim)
  {
    case 1:
      {
        std::vector<int> dims(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
        int offset(fineLocInCoarse[0].first+ghostSize-1),fact0(facts[0]);//offset is always >=0 thanks to the fact that ghostSize>=1 !
        for(int i=0;i<ghostSize;i++)
          outPtr=std::copy(inPtr+offset*nbCompo,inPtr+(offset+1)*nbCompo,outPtr);
        offset=fineLocInCoarse[0].first+ghostSize;
        for(int i=0;i<dims[0];i++)
          {
            const double *loc(inPtr+(offset+i)*nbCompo);
            for(int ifact=0;ifact<fact0;ifact++)
              outPtr=std::copy(loc,loc+nbCompo,outPtr);
          }
        offset=fineLocInCoarse[0].second+ghostSize;
        for(int i=0;i<ghostSize;i++)
          outPtr=std::copy(inPtr+offset*nbCompo,inPtr+(offset+1)*nbCompo,outPtr);
        break;
      }
    case 2:
      {
        SpreadCoarseToFineGhost2D(inPtr,outPtr,nbCompo,coarseSt,fineLocInCoarse,facts,ghostSize);
        break;
      }
    case 3:
      {
        std::vector<int> dims(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
        int fact0(facts[0]),fact1(facts[1]),fact2(facts[2]);
        int nxyWgCoarse((coarseSt[0]+2*ghostSize)*(coarseSt[1]+2*ghostSize)),nxyWgFine((dims[0]*fact0+2*ghostSize)*(dims[1]*fact1+2*ghostSize));
        int offset((fineLocInCoarse[2].first+ghostSize-1)*nxyWgCoarse);//offset is always >=0 thanks to the fact that ghostSize>=1 !
        for(int i=0;i<ghostSize;i++,outPtr+=nxyWgFine*nbCompo)
          SpreadCoarseToFineGhost2D(inPtr+offset*nbCompo,outPtr,nbCompo,coarseSt,fineLocInCoarse,facts,ghostSize);
        offset+=nxyWgCoarse;
        for(int i=0;i<dims[2];i++,offset+=nxyWgCoarse)
          for(int j=0;j<fact2;j++,outPtr+=nxyWgFine*nbCompo)
            SpreadCoarseToFineGhost2D(inPtr+offset*nbCompo,outPtr,nbCompo,coarseSt,fineLocInCoarse,facts,ghostSize);
        for(int i=0;i<ghostSize;i++,outPtr+=nxyWgFine*nbCompo)
          SpreadCoarseToFineGhost2D(inPtr+offset*nbCompo,outPtr,nbCompo,coarseSt,fineLocInCoarse,facts,ghostSize);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhost : only dimensions 1, 2, 3 supported !");
  }
}

/*!
 * This method spreads the values of coarse data \a coarseDA into \a fineDA \b ONLY \b in \b the \b ghost \b zone (contrary to SpreadCoarseToFineGhost that spread the values everywhere).
 *
 * \param [in] coarseDA The DataArrayDouble corresponding to the a cell field of a coarse mesh whose cell structure is defined by \a coarseSt.
 * \param [in] coarseSt The cell structure of coarse mesh.
 * \param [in,out] fineDA The DataArray containing the cell field on uniformly refined mesh
 * \param [in] fineLocInCoarse The cell localization of refined mesh into the coarse one.
 * \param [in] facts The refinement coefficient per axis.
 * \param [in] ghostSize - The size of the ghost zone. The ghost zone is expected to be the same for all axis and both for coarse and fine meshes.
 *
 * \sa SpreadCoarseToFineGhost
 */
void MEDCouplingIMesh::SpreadCoarseToFineGhostZone(const DataArrayDouble *coarseDA, const std::vector<int>& coarseSt, DataArrayDouble *fineDA, const std::vector< std::pair<int,int> >& fineLocInCoarse, const std::vector<int>& facts, int ghostSize)
{
  if(ghostSize<0)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhostZone : ghost level has to be >= 0 !");
  if(coarseSt.size()!=fineLocInCoarse.size() || coarseSt.size()!=facts.size())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhostZone : All input vectors (dimension) must have the same size !");
  if(!coarseDA || !coarseDA->isAllocated() || !fineDA || !fineDA->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhostZone : the parameters 1 or 3 are NULL or not allocated !");
  std::vector<int> coarseStG(coarseSt.size()); std::transform(coarseSt.begin(),coarseSt.end(),coarseStG.begin(),std::bind2nd(std::plus<int>(),2*ghostSize));
  int meshDim((int)coarseSt.size()),nbOfTuplesInCoarseExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(coarseStG));
  int nbCompo(fineDA->getNumberOfComponents());
  if((int)coarseDA->getNumberOfComponents()!=nbCompo)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhostZone : the number of components of fine DA and coarse one mismatches !");
  if(meshDim!=(int)fineLocInCoarse.size() || meshDim!=(int)facts.size())
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhostZone : the size of fineLocInCoarse (4th param) and facts (5th param) must be equal to the sier of coarseSt (2nd param) !");
  if(coarseDA->getNumberOfTuples()!=nbOfTuplesInCoarseExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::SpreadCoarseToFineGhostZone : Expecting " << nbOfTuplesInCoarseExp << " tuples having " << coarseDA->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  std::vector<int> fineStG(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
  std::transform(fineStG.begin(),fineStG.end(),facts.begin(),fineStG.begin(),std::multiplies<int>());
  std::transform(fineStG.begin(),fineStG.end(),fineStG.begin(),std::bind2nd(std::plus<int>(),2*ghostSize));
  int nbTuplesFine(fineDA->getNumberOfTuples()),nbTuplesFineExp(MEDCouplingStructuredMesh::DeduceNumberOfGivenStructure(fineStG));
  if(fineDA->getNumberOfTuples()!=nbTuplesFineExp)
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::SpreadCoarseToFineGhostZone : Expecting " << nbTuplesFineExp << " tuples in fine DataArray having " << nbTuplesFine << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  //
  double *outPtr(fineDA->getPointer());
  const double *inPtr(coarseDA->begin());
  //
  std::vector<int> dims(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
  switch(meshDim)
  {
    case 1:
      {
        int offset(fineLocInCoarse[0].first+ghostSize-1),fact0(facts[0]);//offset is always >=0 thanks to the fact that ghostSize>=1 !
        for(int i=0;i<ghostSize;i++)
          outPtr=std::copy(inPtr+offset*nbCompo,inPtr+(offset+1)*nbCompo,outPtr);
        outPtr+=nbCompo*fact0*dims[0];
        offset=fineLocInCoarse[0].second+ghostSize;
        for(int i=0;i<ghostSize;i++)
          outPtr=std::copy(inPtr+offset*nbCompo,inPtr+(offset+1)*nbCompo,outPtr);
        break;
      }
    case 2:
      {
        SpreadCoarseToFineGhostZone2D(inPtr,outPtr,nbCompo,coarseSt,fineLocInCoarse,facts,ghostSize);
        break;
      }
    case 3:
      {
        int fact0(facts[0]),fact1(facts[1]),fact2(facts[2]);
        int nxyWgCoarse((coarseSt[0]+2*ghostSize)*(coarseSt[1]+2*ghostSize)),nxyWgFine((dims[0]*fact0+2*ghostSize)*(dims[1]*fact1+2*ghostSize));
        int offset((fineLocInCoarse[2].first+ghostSize-1)*nxyWgCoarse);//offset is always >=0 thanks to the fact that ghostSize>=1 !
        for(int i=0;i<ghostSize;i++,outPtr+=nxyWgFine*nbCompo)
          SpreadCoarseToFineGhost2D(inPtr+offset*nbCompo,outPtr,nbCompo,coarseSt,fineLocInCoarse,facts,ghostSize);
        offset+=nxyWgCoarse;
        for(int i=0;i<dims[2];i++,offset+=nxyWgCoarse)
          for(int j=0;j<fact2;j++,outPtr+=nxyWgFine*nbCompo)
            SpreadCoarseToFineGhostZone2D(inPtr+offset*nbCompo,outPtr,nbCompo,coarseSt,fineLocInCoarse,facts,ghostSize);
        for(int i=0;i<ghostSize;i++,outPtr+=nxyWgFine*nbCompo)
          SpreadCoarseToFineGhost2D(inPtr+offset*nbCompo,outPtr,nbCompo,coarseSt,fineLocInCoarse,facts,ghostSize);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingIMesh::SpreadCoarseToFineGhostZone : only dimensions 1, 2, 3 supported !");
  }
}

void MEDCouplingIMesh::setSpaceDimension(int spaceDim)
{
  if(spaceDim==_space_dim)
    return ;
  CheckSpaceDimension(spaceDim);
  _space_dim=spaceDim;
  declareAsNew();
}

void MEDCouplingIMesh::updateTime() const
{
}

std::size_t MEDCouplingIMesh::getHeapMemorySizeWithoutChildren() const
{
  return MEDCouplingStructuredMesh::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDCouplingIMesh::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

/*!
 * This method copyies all tiny strings from other (name and components name).
 * @throw if other and this have not same mesh type.
 */
void MEDCouplingIMesh::copyTinyStringsFrom(const MEDCouplingMesh *other)
{ 
  const MEDCouplingIMesh *otherC=dynamic_cast<const MEDCouplingIMesh *>(other);
  if(!otherC)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::copyTinyStringsFrom : meshes have not same type !");
  MEDCouplingStructuredMesh::copyTinyStringsFrom(other);
  declareAsNew();
}

bool MEDCouplingIMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::isEqualIfNotWhy : input other pointer is null !");
  const MEDCouplingIMesh *otherC(dynamic_cast<const MEDCouplingIMesh *>(other));
  if(!otherC)
    {
      reason="mesh given in input is not castable in MEDCouplingIMesh !";
      return false;
    }
  if(!MEDCouplingStructuredMesh::isEqualIfNotWhy(other,prec,reason))
    return false;
  if(!isEqualWithoutConsideringStrInternal(otherC,prec,reason))
    return false;
  if(_axis_unit!=otherC->_axis_unit)
    {
      reason="The units of axis are not the same !";
      return false;
    }
  return true;
}

bool MEDCouplingIMesh::isEqualWithoutConsideringStr(const MEDCouplingMesh *other, double prec) const
{
  const MEDCouplingIMesh *otherC=dynamic_cast<const MEDCouplingIMesh *>(other);
  if(!otherC)
    return false;
  std::string tmp;
  return isEqualWithoutConsideringStrInternal(other,prec,tmp);
}

bool MEDCouplingIMesh::isEqualWithoutConsideringStrInternal(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  const MEDCouplingIMesh *otherC=dynamic_cast<const MEDCouplingIMesh *>(other);
  if(!otherC)
    return false;
  if(_space_dim!=otherC->_space_dim)
    {
      std::ostringstream oss;
      oss << "The spaceDimension of this (" << _space_dim << ") is not equal to those of other (" << otherC->_space_dim << ") !";
      return false;
    }
  checkSpaceDimension();
  for(int i=0;i<_space_dim;i++)
    {
      if(fabs(_origin[i]-otherC->_origin[i])>prec)
        {
          std::ostringstream oss;
          oss << "The origin of this and other differs at " << i << " !";
          reason=oss.str();
          return false;
        }
    }
  for(int i=0;i<_space_dim;i++)
    {
      if(fabs(_dxyz[i]-otherC->_dxyz[i])>prec)
        {
          std::ostringstream oss;
          oss << "The delta of this and other differs at " << i << " !";
          reason=oss.str();
          return false;
        }
    }
  for(int i=0;i<_space_dim;i++)
    {
      if(_structure[i]!=otherC->_structure[i])
        {
          std::ostringstream oss;
          oss << "The structure of this and other differs at " << i << " !";
          reason=oss.str();
          return false;
        }
    }
  return true;
}

void MEDCouplingIMesh::checkDeepEquivalWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                            DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const
{
  if(!isEqualWithoutConsideringStr(other,prec))
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::checkDeepEquivalWith : Meshes are not the same !");
}

/*!
 * Nothing is done here (except to check that the other is a MEDCoupling::MEDCouplingIMesh instance too).
 * The user intend that the nodes are the same, so by construction of MEDCoupling::MEDCouplingIMesh, \a this and \a other are the same !
 */
void MEDCouplingIMesh::checkDeepEquivalOnSameNodesWith(const MEDCouplingMesh *other, int cellCompPol, double prec,
                                                       DataArrayInt *&cellCor) const
{
  if(!isEqualWithoutConsideringStr(other,prec))
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::checkDeepEquivalOnSameNodesWith : Meshes are not the same !");
}

void MEDCouplingIMesh::checkConsistencyLight() const
{
  checkSpaceDimension();
  for(int i=0;i<_space_dim;i++)
    if(_structure[i]<1)
      {
        std::ostringstream oss; oss << "MEDCouplingIMesh::checkConsistencyLight : On axis " << i << "/" << _space_dim << ", number of nodes is equal to " << _structure[i] << " ! must be >=1 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
}

void MEDCouplingIMesh::checkConsistency(double eps) const
{
  checkConsistencyLight();
}

void MEDCouplingIMesh::getNodeGridStructure(int *res) const
{
  checkSpaceDimension();
  std::copy(_structure,_structure+_space_dim,res);
}

std::vector<int> MEDCouplingIMesh::getNodeGridStructure() const
{
  checkSpaceDimension();
  std::vector<int> ret(_structure,_structure+_space_dim);
  return ret;
}

MEDCouplingStructuredMesh *MEDCouplingIMesh::buildStructuredSubPart(const std::vector< std::pair<int,int> >& cellPart) const
{
  checkConsistencyLight();
  int dim(getSpaceDimension());
  if(dim!=(int)cellPart.size())
    {
      std::ostringstream oss; oss << "MEDCouplingIMesh::buildStructuredSubPart : the space dimension is " << dim << " and cell part size is " << cellPart.size() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  double retOrigin[3]={0.,0.,0.};
  int retStruct[3]={0,0,0};
  MCAuto<MEDCouplingIMesh> ret(dynamic_cast<MEDCouplingIMesh *>(deepCopy()));
  for(int i=0;i<dim;i++)
    {
      int startNode(cellPart[i].first),endNode(cellPart[i].second+1);
      int myDelta(endNode-startNode);
      if(startNode<0 || startNode>=_structure[i])
        {
          std::ostringstream oss; oss << "MEDCouplingIMesh::buildStructuredSubPart : At dimension #" << i << " the start node id is " << startNode << " it should be in [0," << _structure[i] << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(myDelta<0 || myDelta>_structure[i])
        {
          std::ostringstream oss; oss << "MEDCouplingIMesh::buildStructuredSubPart : Along dimension #" << i << " the number of nodes is " << _structure[i] << ", and you are requesting for " << myDelta << " nodes wide range !" << std::endl;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      retOrigin[i]=_origin[i]+startNode*_dxyz[i];
      retStruct[i]=myDelta;
    }
  ret->setNodeStruct(retStruct,retStruct+dim);
  ret->setOrigin(retOrigin,retOrigin+dim);
  ret->checkConsistencyLight();
  return ret.retn();
}

/*!
 * Return the space dimension of \a this.
 */
int MEDCouplingIMesh::getSpaceDimension() const
{
  return _space_dim;
}

void MEDCouplingIMesh::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const
{
  int tmp[3];
  int spaceDim(getSpaceDimension());
  getSplitNodeValues(tmp);
  int tmp2[3];
  GetPosFromId(nodeId,spaceDim,tmp,tmp2);
  for(int j=0;j<spaceDim;j++)
    coo.push_back(_origin[j]+_dxyz[j]*tmp2[j]);
}

std::string MEDCouplingIMesh::simpleRepr() const
{
  std::ostringstream ret;
  ret << "Image grid with name : \"" << getName() << "\"\n";
  ret << "Description of mesh : \"" << getDescription() << "\"\n";
  int tmpp1,tmpp2;
  double tt(getTime(tmpp1,tmpp2));
  int spaceDim(_space_dim);
  ret << "Time attached to the mesh [unit] : " << tt << " [" << getTimeUnit() << "]\n";
  ret << "Iteration : " << tmpp1  << " Order : " << tmpp2 << "\n";
  ret << "Space dimension : " << spaceDim << "\n";
  if(spaceDim<0 || spaceDim>3)
    return ret.str();
  ret << "The nodal structure is : "; std::copy(_structure,_structure+spaceDim,std::ostream_iterator<int>(ret," ")); ret << "\n";
  ret << "The origin position is [" << _axis_unit << "]: ";
  std::copy(_origin,_origin+spaceDim,std::ostream_iterator<double>(ret," ")); ret << "\n";
  ret << "The intervals along axis are : ";
  std::copy(_dxyz,_dxyz+spaceDim,std::ostream_iterator<double>(ret," ")); ret << "\n";
  return ret.str();
}

std::string MEDCouplingIMesh::advancedRepr() const
{
  return simpleRepr();
}

void MEDCouplingIMesh::getBoundingBox(double *bbox) const
{
  checkConsistencyLight();
  int dim(getSpaceDimension());
  for(int idim=0; idim<dim; idim++)
    {
      bbox[2*idim]=_origin[idim];
      int coeff(_structure[idim]);
      if(_structure[idim]<0)
        {
          std::ostringstream oss; oss << "MEDCouplingIMesh::getBoundingBox : on axis #" << idim << " number of nodes in structure is < 0 !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(_structure[idim]>1)
        coeff=_structure[idim]-1;
      bbox[2*idim+1]=_origin[idim]+_dxyz[idim]*coeff;
    }
}

/*!
 * Returns a new MEDCouplingFieldDouble containing volumes of cells constituting \a this
 * mesh.<br>
 * For 1D cells, the returned field contains lengths.<br>
 * For 2D cells, the returned field contains areas.<br>
 * For 3D cells, the returned field contains volumes.
 *  \param [in] isAbs - a not used parameter.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on cells
 *         and one time . The caller is to delete this field using decrRef() as it is no
 *         more needed.
 */
MEDCouplingFieldDouble *MEDCouplingIMesh::getMeasureField(bool isAbs) const
{
  checkConsistencyLight();
  std::string name="MeasureOfMesh_";
  name+=getName();
  int nbelem(getNumberOfCells());
  MEDCouplingFieldDouble *field(MEDCouplingFieldDouble::New(ON_CELLS,ONE_TIME));
  field->setName(name);
  DataArrayDouble* array(DataArrayDouble::New());
  array->alloc(nbelem,1);
  array->fillWithValue(getMeasureOfAnyCell());
  field->setArray(array) ;
  array->decrRef();
  field->setMesh(const_cast<MEDCouplingIMesh *>(this));
  field->synchronizeTimeWithMesh();
  return field;
}

/*!
 * not implemented yet !
 */
MEDCouplingFieldDouble *MEDCouplingIMesh::getMeasureFieldOnNode(bool isAbs) const
{
  throw INTERP_KERNEL::Exception("MEDCouplingIMesh::getMeasureFieldOnNode : not implemented yet !");
  //return 0;
}

int MEDCouplingIMesh::getCellContainingPoint(const double *pos, double eps) const
{
  int dim(getSpaceDimension()),ret(0),coeff(1);
  for(int i=0;i<dim;i++)
    {
      int nbOfCells(_structure[i]-1);
      double ref(pos[i]);
      int tmp((int)((ref-_origin[i])/_dxyz[i]));
      if(tmp>=0 && tmp<nbOfCells)
        {
          ret+=coeff*tmp;
          coeff*=nbOfCells;
        }
      else
        return -1;
    }
  return ret;
}

void MEDCouplingIMesh::getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const
{
  int ret(getCellContainingPoint(pos,eps));
  elts.push_back(ret);
}

void MEDCouplingIMesh::rotate(const double *center, const double *vector, double angle)
{
  throw INTERP_KERNEL::Exception("No rotation available on IMesh : Traduce it to unstructured mesh to apply it !");
}

/*!
 * Translates all nodes of \a this mesh by a given vector. Actually, it adds each
 * component of the \a vector to all node coordinates of a corresponding axis.
 *  \param [in] vector - the translation vector whose size must be not less than \a
 *         this->getSpaceDimension().
 */
void MEDCouplingIMesh::translate(const double *vector)
{
  checkSpaceDimension();
  int dim(getSpaceDimension());
  std::transform(_origin,_origin+dim,vector,_origin,std::plus<double>());
  declareAsNew();
}

/*!
 * Applies scaling transformation to all nodes of \a this mesh.
 *  \param [in] point - coordinates of a scaling center. This array is to be of
 *         size \a this->getSpaceDimension() at least.
 *  \param [in] factor - a scale factor.
 */
void MEDCouplingIMesh::scale(const double *point, double factor)
{
  checkSpaceDimension();
  int dim(getSpaceDimension());
  std::transform(_origin,_origin+dim,point,_origin,std::minus<double>());
  std::transform(_origin,_origin+dim,_origin,std::bind2nd(std::multiplies<double>(),factor));
  std::transform(_dxyz,_dxyz+dim,_dxyz,std::bind2nd(std::multiplies<double>(),factor));
  std::transform(_origin,_origin+dim,point,_origin,std::plus<double>());
  declareAsNew();
}

MEDCouplingMesh *MEDCouplingIMesh::mergeMyselfWith(const MEDCouplingMesh *other) const
{
  //not implemented yet !
  return 0;
}

/*!
 * Returns a new DataArrayDouble holding coordinates of all nodes of \a this mesh.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble, of size \a
 *          this->getNumberOfNodes() tuples per \a this->getSpaceDimension()
 *          components. The caller is to delete this array using decrRef() as it is
 *          no more needed.
 */
DataArrayDouble *MEDCouplingIMesh::getCoordinatesAndOwner() const
{
  checkConsistencyLight();
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  int spaceDim(getSpaceDimension()),nbNodes(getNumberOfNodes());
  ret->alloc(nbNodes,spaceDim);
  double *pt(ret->getPointer());
  ret->setInfoOnComponents(buildInfoOnComponents());
  int tmp2[3],tmp[3];
  getSplitNodeValues(tmp);
  for(int i=0;i<nbNodes;i++)
    {
      GetPosFromId(i,spaceDim,tmp,tmp2);
      for(int j=0;j<spaceDim;j++)
        pt[i*spaceDim+j]=_dxyz[j]*tmp2[j]+_origin[j];
    }
  return ret.retn();
}

/*!
 * Returns a new DataArrayDouble holding barycenters of all cells. The barycenter is
 * computed by averaging coordinates of cell nodes.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble, of size \a
 *          this->getNumberOfCells() tuples per \a this->getSpaceDimension()
 *          components. The caller is to delete this array using decrRef() as it is
 *          no more needed.
 */
DataArrayDouble *MEDCouplingIMesh::computeCellCenterOfMass() const
{
  checkConsistencyLight();
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  int spaceDim(getSpaceDimension()),nbCells(getNumberOfCells()),tmp[3],tmp2[3];
  ret->alloc(nbCells,spaceDim);
  double *pt(ret->getPointer()),shiftOrigin[3];
  std::transform(_dxyz,_dxyz+spaceDim,shiftOrigin,std::bind2nd(std::multiplies<double>(),0.5));
  std::transform(_origin,_origin+spaceDim,shiftOrigin,shiftOrigin,std::plus<double>());
  getSplitCellValues(tmp);
  ret->setInfoOnComponents(buildInfoOnComponents());
  for(int i=0;i<nbCells;i++)
    {
      GetPosFromId(i,spaceDim,tmp,tmp2);
      for(int j=0;j<spaceDim;j++)
        pt[i*spaceDim+j]=_dxyz[j]*tmp2[j]+shiftOrigin[j];
    }
  return ret.retn();
}

DataArrayDouble *MEDCouplingIMesh::computeIsoBarycenterOfNodesPerCell() const
{
  return MEDCouplingIMesh::computeCellCenterOfMass();
}

void MEDCouplingIMesh::renumberCells(const int *old2NewBg, bool check)
{
  throw INTERP_KERNEL::Exception("Functionnality of renumbering cell not available for IMesh !");
}

void MEDCouplingIMesh::getTinySerializationInformation(std::vector<double>& tinyInfoD, std::vector<int>& tinyInfo, std::vector<std::string>& littleStrings) const
{
  int it,order;
  double time(getTime(it,order));
  tinyInfo.clear();
  tinyInfoD.clear();
  littleStrings.clear();
  littleStrings.push_back(getName());
  littleStrings.push_back(getDescription());
  littleStrings.push_back(getTimeUnit());
  littleStrings.push_back(getAxisUnit());
  tinyInfo.push_back(it);
  tinyInfo.push_back(order);
  tinyInfo.push_back(_space_dim);
  tinyInfo.insert(tinyInfo.end(),_structure,_structure+3);
  tinyInfoD.push_back(time);
  tinyInfoD.insert(tinyInfoD.end(),_dxyz,_dxyz+3);
  tinyInfoD.insert(tinyInfoD.end(),_origin,_origin+3);
}

void MEDCouplingIMesh::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2, std::vector<std::string>& littleStrings) const
{
  a1->alloc(0,1);
  a2->alloc(0,1);
}

void MEDCouplingIMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2) const
{
  a1=DataArrayInt::New();
  a1->alloc(0,1);
  a2=DataArrayDouble::New();
  a2->alloc(0,1);
}

void MEDCouplingIMesh::unserialization(const std::vector<double>& tinyInfoD, const std::vector<int>& tinyInfo, const DataArrayInt *a1, DataArrayDouble *a2,
                                       const std::vector<std::string>& littleStrings)
{
  setName(littleStrings[0]);
  setDescription(littleStrings[1]);
  setTimeUnit(littleStrings[2]);
  setAxisUnit(littleStrings[3]);
  setTime(tinyInfoD[0],tinyInfo[0],tinyInfo[1]);
  _space_dim=tinyInfo[2];
  _structure[0]=tinyInfo[3]; _structure[1]=tinyInfo[4]; _structure[2]=tinyInfo[5];
  _dxyz[0]=tinyInfoD[1]; _dxyz[1]=tinyInfoD[2]; _dxyz[2]=tinyInfoD[3];
  _origin[0]=tinyInfoD[4]; _origin[1]=tinyInfoD[5]; _origin[2]=tinyInfoD[6];
  declareAsNew();
}

void MEDCouplingIMesh::writeVTKLL(std::ostream& ofs, const std::string& cellData, const std::string& pointData, DataArrayByte *byteData) const
{
  checkConsistencyLight();
  std::ostringstream extent,origin,spacing;
  for(int i=0;i<3;i++)
    {
      if(i<_space_dim)
        { extent << "0 " <<  _structure[i]-1 << " "; origin << _origin[i] << " "; spacing << _dxyz[i] << " "; }
      else
        { extent << "0 0 "; origin << "0 "; spacing << "0 "; }
    }
  ofs << "  <" << getVTKDataSetType() << " WholeExtent=\"" << extent.str() << "\" Origin=\"" << origin.str() << "\" Spacing=\"" << spacing.str() << "\">\n";
  ofs << "    <Piece Extent=\"" << extent.str() << "\">\n";
  ofs << "      <PointData>\n" << pointData << std::endl;
  ofs << "      </PointData>\n";
  ofs << "      <CellData>\n" << cellData << std::endl;
  ofs << "      </CellData>\n";
  ofs << "      <Coordinates>\n";
  ofs << "      </Coordinates>\n";
  ofs << "    </Piece>\n";
  ofs << "  </" << getVTKDataSetType() << ">\n";
}

void MEDCouplingIMesh::reprQuickOverview(std::ostream& stream) const
{
  stream << "MEDCouplingIMesh C++ instance at " << this << ". Name : \"" << getName() << "\". Space dimension : " << _space_dim << ".";
  if(_space_dim<0 || _space_dim>3)
    return ;
  stream << "\n";
  std::ostringstream stream0,stream1;
  int nbNodes(1),nbCells(0);
  bool isPb(false);
  for(int i=0;i<_space_dim;i++)
    {
      char tmp('X'+i);
      int tmpNodes(_structure[i]);
      stream1 << "- Axis " << tmp << " : " << tmpNodes << " nodes (orig=" << _origin[i] << ", inter=" << _dxyz[i] << ").";
      if(i!=_space_dim-1)
        stream1 << std::endl;
      if(tmpNodes>=1)
        nbNodes*=tmpNodes;
      else
        isPb=true;
      if(tmpNodes>=2)
        nbCells=nbCells==0?tmpNodes-1:nbCells*(tmpNodes-1);
    }
  if(!isPb)
    {
      stream0 << "Number of cells : " << nbCells << ", Number of nodes : " << nbNodes;
      stream << stream0.str();
      if(_space_dim>0)
        stream << std::endl;
    }
  stream << stream1.str();
}

std::string MEDCouplingIMesh::getVTKFileExtension() const
{
  return std::string("vti");
}

std::string MEDCouplingIMesh::getVTKDataSetType() const
{
  return std::string("ImageData");
}

std::vector<std::string> MEDCouplingIMesh::buildInfoOnComponents() const
{
  checkSpaceDimension();
  int dim(getSpaceDimension());
  std::vector<std::string> ret(dim);
  for(int i=0;i<dim;i++)
    {
      std::ostringstream oss;
      char tmp('X'+i); oss << tmp;
      ret[i]=DataArray::BuildInfoFromVarAndUnit(oss.str(),_axis_unit);
    }
  return ret;
}

void MEDCouplingIMesh::checkSpaceDimension() const
{
  CheckSpaceDimension(_space_dim);
}

void MEDCouplingIMesh::CheckSpaceDimension(int spaceDim)
{
  if(spaceDim<0 || spaceDim>3)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::CheckSpaceDimension : input spaceDim must be in [0,1,2,3] !");
}

int MEDCouplingIMesh::FindIntRoot(int val, int order)
{
  if(order==0)
    return 1;
  if(val<0)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::FindIntRoot : input val is < 0 ! Not possible to compute a root !");
  if(order==1)
    return val;
  if(order!=2 && order!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingIMesh::FindIntRoot : the order available are 0,1,2 or 3 !");
  double valf((double)val);
  if(order==2)
    {
      double retf(sqrt(valf));
      int ret((int)retf);
      if(ret*ret!=val)
        throw INTERP_KERNEL::Exception("MEDCouplingIMesh::FindIntRoot : the input val is not a perfect square root !");
      return ret;
    }
  else//order==3
    {
      double retf(std::pow(val,0.3333333333333333));
      int ret((int)retf),ret2(ret+1);
      if(ret*ret*ret!=val && ret2*ret2*ret2!=val)
        throw INTERP_KERNEL::Exception("MEDCouplingIMesh::FindIntRoot : the input val is not a perfect cublic root !");
      if(ret*ret*ret==val)
        return ret;
      else
        return ret2;
    }
}

void MEDCouplingIMesh::SpreadCoarseToFineGhost2D(const double *inPtr, double *outPtr, int nbCompo, const std::vector<int>& coarseSt, const std::vector< std::pair<int,int> >& fineLocInCoarse, const std::vector<int>& facts, int ghostSize)
{
  double *outPtrWork(outPtr);
  std::vector<int> dims(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
  int nxwg(coarseSt[0]+2*ghostSize),fact0(facts[0]),fact1(facts[1]);
  int kk(fineLocInCoarse[0].first+ghostSize-1+nxwg*(fineLocInCoarse[1].first+ghostSize-1));//kk is always >=0 thanks to the fact that ghostSize>=1 !
  for(int jg=0;jg<ghostSize;jg++)
    {
      for(int ig=0;ig<ghostSize;ig++)
        outPtrWork=std::copy(inPtr+kk*nbCompo,inPtr+(kk+1)*nbCompo,outPtrWork);
      int kk0(kk+1);
      for(int ig=0;ig<dims[0];ig++,kk0++)
        for(int ifact=0;ifact<fact0;ifact++)
          outPtrWork=std::copy(inPtr+(kk0)*nbCompo,inPtr+(kk0+1)*nbCompo,outPtrWork);
      for(int ik=0;ik<ghostSize;ik++)
        outPtrWork=std::copy(inPtr+kk0*nbCompo,inPtr+(kk0+1)*nbCompo,outPtrWork);
    }
  for(int j=0;j<dims[1];j++)
    {
      kk=fineLocInCoarse[0].first-1+ghostSize+nxwg*(fineLocInCoarse[1].first+ghostSize+j);
      for(int jfact=0;jfact<fact1;jfact++)
        {
          for(int ig=0;ig<ghostSize;ig++)
            outPtrWork=std::copy(inPtr+kk*nbCompo,inPtr+(kk+1)*nbCompo,outPtrWork);
          int kk0(kk+1);//1 not ghost. We make the hypothesis that factors is >= ghostlev
          for(int i=0;i<dims[0];i++,kk0++)
            {
              const double *loc(inPtr+kk0*nbCompo);
              for(int ifact=0;ifact<fact0;ifact++)
                outPtrWork=std::copy(loc,loc+nbCompo,outPtrWork);
            }
          for(int ig=0;ig<ghostSize;ig++)
            outPtrWork=std::copy(inPtr+kk0*nbCompo,inPtr+(kk0+1)*nbCompo,outPtrWork);
        }
    }
  kk=fineLocInCoarse[0].first+ghostSize-1+nxwg*(fineLocInCoarse[1].second+ghostSize);
  for(int jg=0;jg<ghostSize;jg++)
    {
      for(int ig=0;ig<ghostSize;ig++)
        outPtrWork=std::copy(inPtr+kk*nbCompo,inPtr+(kk+1)*nbCompo,outPtrWork);
      int kk0(kk+1);
      for(int ig=0;ig<dims[0];ig++,kk0++)
        for(int ifact=0;ifact<fact0;ifact++)
          outPtrWork=std::copy(inPtr+(kk0)*nbCompo,inPtr+(kk0+1)*nbCompo,outPtrWork);
      for(int ik=0;ik<ghostSize;ik++)
        outPtrWork=std::copy(inPtr+kk0*nbCompo,inPtr+(kk0+1)*nbCompo,outPtrWork);
    }
}

void MEDCouplingIMesh::SpreadCoarseToFineGhostZone2D(const double *inPtr, double *outPtr, int nbCompo, const std::vector<int>& coarseSt, const std::vector< std::pair<int,int> >& fineLocInCoarse, const std::vector<int>& facts, int ghostSize)
{
  double *outPtr2(outPtr);
  std::vector<int> dims(MEDCouplingStructuredMesh::GetDimensionsFromCompactFrmt(fineLocInCoarse));
  int nxwg(coarseSt[0]+2*ghostSize),fact0(facts[0]),fact1(facts[1]);
  int kk(fineLocInCoarse[0].first+ghostSize-1+nxwg*(fineLocInCoarse[1].first+ghostSize-1));//kk is always >=0 thanks to the fact that ghostSize>=1 !
  for(int jg=0;jg<ghostSize;jg++)
    {
      for(int ig=0;ig<ghostSize;ig++)
        outPtr2=std::copy(inPtr+kk*nbCompo,inPtr+(kk+1)*nbCompo,outPtr2);
      int kk0(kk+1);
      for(int ig=0;ig<dims[0];ig++,kk0++)
        for(int ifact=0;ifact<fact0;ifact++)
          outPtr2=std::copy(inPtr+(kk0)*nbCompo,inPtr+(kk0+1)*nbCompo,outPtr2);
      for(int ik=0;ik<ghostSize;ik++)
        outPtr2=std::copy(inPtr+kk0*nbCompo,inPtr+(kk0+1)*nbCompo,outPtr2);
    }
  for(int j=0;j<dims[1];j++)
    {
      kk=fineLocInCoarse[0].first-1+ghostSize+nxwg*(fineLocInCoarse[1].first+ghostSize+j);
      for(int jfact=0;jfact<fact1;jfact++)
        {
          for(int ig=0;ig<ghostSize;ig++)
            outPtr2=std::copy(inPtr+kk*nbCompo,inPtr+(kk+1)*nbCompo,outPtr2);
          int kk0(kk+1+dims[0]);//1 not ghost. We make the hypothesis that factors is >= ghostlev
          outPtr2+=fact0*nbCompo*dims[0];
          for(int ig=0;ig<ghostSize;ig++)
            outPtr2=std::copy(inPtr+kk0*nbCompo,inPtr+(kk0+1)*nbCompo,outPtr2);
        }
    }
  kk=fineLocInCoarse[0].first+ghostSize-1+nxwg*(fineLocInCoarse[1].second+ghostSize);
  for(int jg=0;jg<ghostSize;jg++)
    {
      for(int ig=0;ig<ghostSize;ig++)
        outPtr2=std::copy(inPtr+kk*nbCompo,inPtr+(kk+1)*nbCompo,outPtr2);
      int kk0(kk+1);
      for(int ig=0;ig<dims[0];ig++,kk0++)
        for(int ifact=0;ifact<fact0;ifact++)
          outPtr2=std::copy(inPtr+(kk0)*nbCompo,inPtr+(kk0+1)*nbCompo,outPtr2);
      for(int ik=0;ik<ghostSize;ik++)
        outPtr2=std::copy(inPtr+kk0*nbCompo,inPtr+(kk0+1)*nbCompo,outPtr2);
    }
}
