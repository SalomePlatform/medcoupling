// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include "CellModel.hxx"
#include "InterpolationUtils.hxx"
#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGaussCoords.hxx"
#include "InterpKernelMatrixTools.hxx"

#include <set>
#include <list>
#include <limits>
#include <sstream>
#include <algorithm>
#include <functional>

using namespace ParaMEDMEM;

const double MEDCouplingFieldDiscretization::DFLT_PRECISION=1.e-12;

const char MEDCouplingFieldDiscretizationP0::REPR[]="P0";

const TypeOfField MEDCouplingFieldDiscretizationP0::TYPE=ON_CELLS;

const char MEDCouplingFieldDiscretizationP1::REPR[]="P1";

const TypeOfField MEDCouplingFieldDiscretizationP1::TYPE=ON_NODES;

const int MEDCouplingFieldDiscretizationPerCell::DFT_INVALID_LOCID_VALUE=-1;

const char MEDCouplingFieldDiscretizationGauss::REPR[]="GAUSS";

const TypeOfField MEDCouplingFieldDiscretizationGauss::TYPE=ON_GAUSS_PT;

const char MEDCouplingFieldDiscretizationGaussNE::REPR[]="GSSNE";

const TypeOfField MEDCouplingFieldDiscretizationGaussNE::TYPE=ON_GAUSS_NE;

const char MEDCouplingFieldDiscretizationKriging::REPR[]="KRIGING";

const TypeOfField MEDCouplingFieldDiscretizationKriging::TYPE=ON_NODES_KR;

MEDCouplingFieldDiscretization::MEDCouplingFieldDiscretization():_precision(DFLT_PRECISION)
{
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretization::New(TypeOfField type)
{
  switch(type)
    {
    case MEDCouplingFieldDiscretizationP0::TYPE:
      return new MEDCouplingFieldDiscretizationP0;
    case MEDCouplingFieldDiscretizationP1::TYPE:
      return new MEDCouplingFieldDiscretizationP1;
    case MEDCouplingFieldDiscretizationGauss::TYPE:
      return new MEDCouplingFieldDiscretizationGauss;
    case MEDCouplingFieldDiscretizationGaussNE::TYPE:
      return new MEDCouplingFieldDiscretizationGaussNE;
    case MEDCouplingFieldDiscretizationKriging::TYPE:
      return new MEDCouplingFieldDiscretizationKriging;
    default:
      throw INTERP_KERNEL::Exception("Choosen discretization is not implemented yet.");
    }
}

TypeOfField MEDCouplingFieldDiscretization::getTypeOfFieldFromStringRepr(const char *repr) throw(INTERP_KERNEL::Exception)
{
  std::string reprCpp(repr);
  if(reprCpp==MEDCouplingFieldDiscretizationP0::REPR)
    return MEDCouplingFieldDiscretizationP0::TYPE;
  if(reprCpp==MEDCouplingFieldDiscretizationP1::REPR)
    return MEDCouplingFieldDiscretizationP1::TYPE;
  if(reprCpp==MEDCouplingFieldDiscretizationGauss::REPR)
    return MEDCouplingFieldDiscretizationGauss::TYPE;
  if(reprCpp==MEDCouplingFieldDiscretizationGaussNE::REPR)
    return MEDCouplingFieldDiscretizationGaussNE::TYPE;
  if(reprCpp==MEDCouplingFieldDiscretizationKriging::REPR)
    return MEDCouplingFieldDiscretizationKriging::TYPE;
  throw INTERP_KERNEL::Exception("Representation does not match with any field discretization !");
}

bool MEDCouplingFieldDiscretization::isEqual(const MEDCouplingFieldDiscretization *other, double eps) const
{
  std::string reason;
  return isEqualIfNotWhy(other,eps,reason);
}

bool MEDCouplingFieldDiscretization::isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const
{
  return isEqual(other,eps);
}

/*!
 * For all field discretization excepted GaussPts the [ \a startCellIds, \a endCellIds ) has no impact on the cloned instance.
 */
MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretization::clonePart(const int *startCellIds, const int *endCellIds) const
{
  return clone();
}

/*!
 * Excepted for MEDCouplingFieldDiscretizationPerCell no underlying TimeLabel object : nothing to do in generally.
 */
void MEDCouplingFieldDiscretization::updateTime() const
{
}

/*!
 * Computes normL1 of DataArrayDouble instance arr.
 * @param res output parameter expected to be of size arr->getNumberOfComponents();
 * @throw when the field discretization fails on getMeasure fields (gauss points for example)
 */
void MEDCouplingFieldDiscretization::normL1(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, double *res) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingFieldDouble *vol=getMeasureField(mesh,true);
  int nbOfCompo=arr->getNumberOfComponents();
  int nbOfElems=getNumberOfTuples(mesh);
  std::fill(res,res+nbOfCompo,0.);
  const double *arrPtr=arr->getConstPointer();
  const double *volPtr=vol->getArray()->getConstPointer();
  double deno=0.;
  for(int i=0;i<nbOfElems;i++)
    {
      double v=fabs(volPtr[i]);
      for(int j=0;j<nbOfCompo;j++)
        res[j]+=fabs(arrPtr[i*nbOfCompo+j])*v;
      deno+=v;
    }
  std::transform(res,res+nbOfCompo,res,std::bind2nd(std::multiplies<double>(),1./deno));
  vol->decrRef();
}

/*!
 * Computes normL2 of DataArrayDouble instance arr.
 * @param res output parameter expected to be of size arr->getNumberOfComponents();
 * @throw when the field discretization fails on getMeasure fields (gauss points for example)
 */
void MEDCouplingFieldDiscretization::normL2(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, double *res) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingFieldDouble *vol=getMeasureField(mesh,true);
  int nbOfCompo=arr->getNumberOfComponents();
  int nbOfElems=getNumberOfTuples(mesh);
  std::fill(res,res+nbOfCompo,0.);
  const double *arrPtr=arr->getConstPointer();
  const double *volPtr=vol->getArray()->getConstPointer();
  double deno=0.;
  for(int i=0;i<nbOfElems;i++)
    {
      double v=fabs(volPtr[i]);
      for(int j=0;j<nbOfCompo;j++)
        res[j]+=arrPtr[i*nbOfCompo+j]*arrPtr[i*nbOfCompo+j]*v;
      deno+=v;
    }
  std::transform(res,res+nbOfCompo,res,std::bind2nd(std::multiplies<double>(),1./deno));
  std::transform(res,res+nbOfCompo,res,std::ptr_fun<double,double>(std::sqrt));
  vol->decrRef();
}

/*!
 * Computes integral of DataArrayDouble instance arr.
 * @param res output parameter expected to be of size arr->getNumberOfComponents();
 * @throw when the field discretization fails on getMeasure fields (gauss points for example)
 */
void MEDCouplingFieldDiscretization::integral(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, bool isWAbs, double *res) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingFieldDouble *vol=getMeasureField(mesh,isWAbs);
  int nbOfCompo=arr->getNumberOfComponents();
  int nbOfElems=getNumberOfTuples(mesh);
  std::fill(res,res+nbOfCompo,0.);
  const double *arrPtr=arr->getConstPointer();
  const double *volPtr=vol->getArray()->getConstPointer();
  double *tmp=new double[nbOfCompo];
  for (int i=0;i<nbOfElems;i++)
    {
      std::transform(arrPtr+i*nbOfCompo,arrPtr+(i+1)*nbOfCompo,tmp,std::bind2nd(std::multiplies<double>(),volPtr[i]));
      std::transform(tmp,tmp+nbOfCompo,res,res,std::plus<double>());
    }
  delete [] tmp;
  vol->decrRef();
}

void MEDCouplingFieldDiscretization::getSerializationIntArray(DataArrayInt *& arr) const
{
  arr=0;
}

/*!
 * Empty : Not a bug
 */
void MEDCouplingFieldDiscretization::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
}

/*!
 * Empty : Not a bug
 */
void MEDCouplingFieldDiscretization::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
}

void MEDCouplingFieldDiscretization::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *& arr)
{
  arr=0;
}

/*!
 * Empty : Not a bug
 */
void MEDCouplingFieldDiscretization::finishUnserialization(const std::vector<double>& tinyInfo)
{
}

/*!
 * This method is typically the first step of renumbering. The implementation is empty it is not a bug only gauss is impacted
 * virtualy by this method.
 */
void MEDCouplingFieldDiscretization::renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
}

double MEDCouplingFieldDiscretization::getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da,
                                              int cellId, int nodeIdInCell, int compoId) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("getIJK Invalid ! only for GaussPoint and GaussNE discretizations !");
}

void MEDCouplingFieldDiscretization::setGaussLocalizationOnType(const MEDCouplingMesh *m, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                                const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

void MEDCouplingFieldDiscretization::setGaussLocalizationOnCells(const MEDCouplingMesh *m, const int *begin, const int *end, const std::vector<double>& refCoo,
                                                                 const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

void MEDCouplingFieldDiscretization::clearGaussLocalizations() throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

MEDCouplingGaussLocalization& MEDCouplingFieldDiscretization::getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

const MEDCouplingGaussLocalization& MEDCouplingFieldDiscretization::getGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

int MEDCouplingFieldDiscretization::getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

int MEDCouplingFieldDiscretization::getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

int MEDCouplingFieldDiscretization::getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

std::set<int> MEDCouplingFieldDiscretization::getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

void MEDCouplingFieldDiscretization::getCellIdsHavingGaussLocalization(int locId, std::vector<int>& cellIds) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

void MEDCouplingFieldDiscretization::RenumberEntitiesFromO2NArr(double eps, const int *old2NewPtr, int newNbOfEntity, DataArrayDouble *arr, const char *msg)
{
  int oldNbOfElems=arr->getNumberOfTuples();
  int nbOfComp=arr->getNumberOfComponents();
  int newNbOfTuples=newNbOfEntity;
  DataArrayDouble *arrCpy=arr->deepCpy();
  const double *ptSrc=arrCpy->getConstPointer();
  arr->reAlloc(newNbOfTuples);
  double *ptToFill=arr->getPointer();
  std::fill(ptToFill,ptToFill+nbOfComp*newNbOfTuples,std::numeric_limits<double>::max());
  INTERP_KERNEL::AutoPtr<double> tmp=new double[nbOfComp];
  for(int i=0;i<oldNbOfElems;i++)
    {
      int newNb=old2NewPtr[i];
      if(newNb>=0)//if newNb<0 the node is considered as out.
        {
          if(std::find_if(ptToFill+newNb*nbOfComp,ptToFill+(newNb+1)*nbOfComp,std::bind2nd(std::not_equal_to<double>(),std::numeric_limits<double>::max()))
             ==ptToFill+(newNb+1)*nbOfComp)
            std::copy(ptSrc+i*nbOfComp,ptSrc+(i+1)*nbOfComp,ptToFill+newNb*nbOfComp);
          else
            {
              std::transform(ptSrc+i*nbOfComp,ptSrc+(i+1)*nbOfComp,ptToFill+newNb*nbOfComp,(double *)tmp,std::minus<double>());
              std::transform((double *)tmp,((double *)tmp)+nbOfComp,(double *)tmp,std::ptr_fun<double,double>(fabs));
              //if(!std::equal(ptSrc+i*nbOfComp,ptSrc+(i+1)*nbOfComp,ptToFill+newNb*nbOfComp))
              if(*std::max_element((double *)tmp,((double *)tmp)+nbOfComp)>eps)
                {
                  arrCpy->decrRef();
                  std::ostringstream oss;
                  oss << msg << " " << i << " and " << std::find(old2NewPtr,old2NewPtr+i,newNb)-old2NewPtr
                      << " have been merged and " << msg << " field on them are different !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
        }
    }
  arrCpy->decrRef();
}

void MEDCouplingFieldDiscretization::RenumberEntitiesFromN2OArr(const int *new2OldPtr, int new2OldSz, DataArrayDouble *arr, const char *msg)
{
  int nbOfComp=arr->getNumberOfComponents();
  DataArrayDouble *arrCpy=arr->deepCpy();
  const double *ptSrc=arrCpy->getConstPointer();
  arr->reAlloc(new2OldSz);
  double *ptToFill=arr->getPointer();
  for(int i=0;i<new2OldSz;i++)
    {
      int oldNb=new2OldPtr[i];
      std::copy(ptSrc+oldNb*nbOfComp,ptSrc+(oldNb+1)*nbOfComp,ptToFill+i*nbOfComp);
    }
  arrCpy->decrRef();
}

MEDCouplingFieldDiscretization::~MEDCouplingFieldDiscretization()
{
}

TypeOfField MEDCouplingFieldDiscretizationP0::getEnum() const
{
  return TYPE;
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationP0::clone() const
{
  return new MEDCouplingFieldDiscretizationP0;
}

std::string MEDCouplingFieldDiscretizationP0::getStringRepr() const
{
  return std::string(REPR);
}

const char *MEDCouplingFieldDiscretizationP0::getRepr() const
{
  return REPR;
}

bool MEDCouplingFieldDiscretizationP0::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  const MEDCouplingFieldDiscretizationP0 *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationP0 *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason="Spatial discrtization of this is ON_CELLS, which is not the case of other.";
  return ret;
}

int MEDCouplingFieldDiscretizationP0::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfCells();
}

int MEDCouplingFieldDiscretizationP0::getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfCells();
}

DataArrayInt *MEDCouplingFieldDiscretizationP0::getOffsetArr(const MEDCouplingMesh *mesh) const
{
  int nbOfTuples=mesh->getNumberOfCells();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfTuples+1,1);
  ret->iota(0);
  return ret;
}

void MEDCouplingFieldDiscretizationP0::renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                                                             const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  const int *array=old2NewBg;
  if(check)
    array=DataArrayInt::CheckAndPreparePermutation(old2NewBg,old2NewBg+mesh->getNumberOfCells());
  for(std::vector<DataArrayDouble *>::const_iterator it=arrays.begin();it!=arrays.end();it++)
    {
      if(*it)
        (*it)->renumberInPlace(array);
    }
  if(check)
    delete [] array;
}

DataArrayDouble *MEDCouplingFieldDiscretizationP0::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  return mesh->getBarycenterAndOwner();
}

void MEDCouplingFieldDiscretizationP0::computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *partBg, const int *partEnd,
                                                                          DataArrayInt *&cellRest)
{
  cellRest=DataArrayInt::New();
  cellRest->alloc((int)std::distance(partBg,partEnd),1);
  std::copy(partBg,partEnd,cellRest->getPointer());
}

void MEDCouplingFieldDiscretizationP0::checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception)
{
}

void MEDCouplingFieldDiscretizationP0::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception)
{
  if(mesh->getNumberOfCells()!=da->getNumberOfTuples())
    {
      std::ostringstream message;
      message << "Field on cells invalid because there are " << mesh->getNumberOfCells();
      message << " cells in mesh and " << da->getNumberOfTuples() << " tuples in field !";
      throw INTERP_KERNEL::Exception(message.str().c_str());
    }
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationP0::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  return mesh->getMeasureField(isAbs);
}

void MEDCouplingFieldDiscretizationP0::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  int id=mesh->getCellContainingPoint(loc,_precision);
  if(id==-1)
    throw INTERP_KERNEL::Exception("Specified point is detected outside of mesh : unable to apply P0::getValueOn !");
  arr->getTuple(id,res);
}

void MEDCouplingFieldDiscretizationP0::getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const
{
  const MEDCouplingCMesh *meshC=dynamic_cast<const MEDCouplingCMesh *>(mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("P0::getValueOnPos is only accessible for structured meshes !");
  int id=meshC->getCellIdFromPos(i,j,k);
  arr->getTuple(id,res);
}

DataArrayDouble *MEDCouplingFieldDiscretizationP0::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const
{
  std::vector<int> elts,eltsIndex;
  mesh->getCellsContainingPoints(loc,nbOfPoints,_precision,elts,eltsIndex);
  int spaceDim=mesh->getSpaceDimension();
  int nbOfComponents=arr->getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(nbOfPoints,nbOfComponents);
  double *ptToFill=ret->getPointer();
  for(int i=0;i<nbOfPoints;i++,ptToFill+=nbOfComponents)
    if(eltsIndex[i+1]-eltsIndex[i]>=1)
      arr->getTuple(elts[eltsIndex[i]],ptToFill);
    else
      {
        std::ostringstream oss; oss << "Point #" << i << " with coordinates : (";
        std::copy(loc+i*spaceDim,loc+(i+1)*spaceDim,std::ostream_iterator<double>(oss,", "));
        oss << ") detected outside mesh : unable to apply P0::getValueOnMulti ! ";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  ret->incrRef();
  return ret;
}

/*!
 * Nothing to do. It's not a bug.
 */
void MEDCouplingFieldDiscretizationP0::renumberValuesOnNodes(double , const int *, int newNbOfNodes, DataArrayDouble *) const
{
}

void MEDCouplingFieldDiscretizationP0::renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const int *old2New, int newSz, DataArrayDouble *arr) const
{
  RenumberEntitiesFromO2NArr(epsOnVals,old2New,newSz,arr,"Cell");
}

void MEDCouplingFieldDiscretizationP0::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const
{
  RenumberEntitiesFromN2OArr(new2old,newSz,arr,"Cell");
}

/*!
 * This method returns a tuple ids selection from cell ids selection [start;end).
 * This method is called by MEDCouplingFieldDiscretizationP0::buildSubMeshData to return parameter \b di.
 * Here for P0 it's very simple !
 *
 * \return a newly allocated array containing ids to select into the DataArrayDouble of the field.
 * 
 */
DataArrayInt *MEDCouplingFieldDiscretizationP0::computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const int *startCellIds, const int *endCellIds) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc((int)std::distance(startCellIds,endCellIds),1);
  std::copy(startCellIds,endCellIds,ret->getPointer());
  ret->incrRef(); return ret;
}

/*!
 * This method returns a submesh of 'mesh' instance constituting cell ids contained in array defined as an interval [start;end).
 * @param di is an array returned that specifies entity ids (here cells ids) in mesh 'mesh' of entity in returned submesh.
 * Example : The first cell id of returned mesh has the (*di)[0] id in 'mesh'
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationP0::buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const
{
  MEDCouplingMesh *ret=mesh->buildPart(start,end);
  di=DataArrayInt::New();
  di->alloc((int)std::distance(start,end),1);
  int *pt=di->getPointer();
  std::copy(start,end,pt);
  return ret;
}

int MEDCouplingFieldDiscretizationOnNodes::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfNodes();
}

int MEDCouplingFieldDiscretizationOnNodes::getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfNodes();
}

/*!
 * Nothing to do here.
 */
void MEDCouplingFieldDiscretizationOnNodes::renumberArraysForCell(const MEDCouplingMesh *, const std::vector<DataArrayDouble *>& arrays,
                                                                  const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
}

DataArrayInt *MEDCouplingFieldDiscretizationOnNodes::getOffsetArr(const MEDCouplingMesh *mesh) const
{
  int nbOfTuples=mesh->getNumberOfNodes();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfTuples+1,1);
  ret->iota(0);
  return ret;
}

DataArrayDouble *MEDCouplingFieldDiscretizationOnNodes::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  return mesh->getCoordinatesAndOwner();
}

void MEDCouplingFieldDiscretizationOnNodes::computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *partBg, const int *partEnd,
                                                                               DataArrayInt *&cellRest)
{
  cellRest=mesh->getCellIdsFullyIncludedInNodeIds(partBg,partEnd);
}

void MEDCouplingFieldDiscretizationOnNodes::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception)
{
  if(mesh->getNumberOfNodes()!=da->getNumberOfTuples())
    {
      std::ostringstream message;
      message << "Field on nodes invalid because there are " << mesh->getNumberOfNodes();
      message << " nodes in mesh and " << da->getNumberOfTuples() << " tuples in field !";
      throw INTERP_KERNEL::Exception(message.str().c_str());
    }
}

/*!
 * This method returns a submesh of 'mesh' instance constituting cell ids contained in array defined as an interval [start;end).
* @param di is an array returned that specifies entity ids (here nodes ids) in mesh 'mesh' of entity in returned submesh.
 * Example : The first node id of returned mesh has the (*di)[0] id in 'mesh'
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationOnNodes::buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const
{
  MEDCouplingMesh *ret=mesh->buildPartAndReduceNodes(start,end,di);
  DataArrayInt *di2=di->invertArrayO2N2N2O(ret->getNumberOfNodes());
  di->decrRef();
  di=di2;
  return ret;
}

/*!
 * This method returns a tuple ids selection from cell ids selection [start;end).
 * This method is called by MEDCouplingFieldDiscretizationP0::buildSubMeshData to return parameter \b di.
 * Here for P1 only nodes fetched by submesh of mesh[startCellIds:endCellIds) is returned !
 *
 * \return a newly allocated array containing ids to select into the DataArrayDouble of the field.
 * 
 */
DataArrayInt *MEDCouplingFieldDiscretizationOnNodes::computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const int *startCellIds, const int *endCellIds) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP1::computeTupleIdsToSelectFromCellIds : null mesh !");
  const MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> umesh=mesh->buildUnstructured();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> umesh2=static_cast<MEDCouplingUMesh *>(umesh->buildPartOfMySelf(startCellIds,endCellIds,true));
  return umesh2->computeFetchedNodeIds();
}

void MEDCouplingFieldDiscretizationOnNodes::renumberValuesOnNodes(double epsOnVals, const int *old2NewPtr, int newNbOfNodes, DataArrayDouble *arr) const
{
  RenumberEntitiesFromO2NArr(epsOnVals,old2NewPtr,newNbOfNodes,arr,"Node");
}

/*!
 * Nothing to do it's not a bug.
 */
void MEDCouplingFieldDiscretizationOnNodes::renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const int *old2New, int newSz, DataArrayDouble *arr) const
{
}

/*!
 * Nothing to do it's not a bug.
 */
void MEDCouplingFieldDiscretizationOnNodes::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const
{
}

void MEDCouplingFieldDiscretizationOnNodes::getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const
{
  const MEDCouplingCMesh *meshC=dynamic_cast<const MEDCouplingCMesh *>(mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("OnNodes::getValueOnPos(i,j,k) is only accessible for structured meshes !");
  int id=meshC->getNodeIdFromPos(i,j,k);
  arr->getTuple(id,res);
}

TypeOfField MEDCouplingFieldDiscretizationP1::getEnum() const
{
  return TYPE;
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationP1::clone() const
{
  return new MEDCouplingFieldDiscretizationP1;
}

std::string MEDCouplingFieldDiscretizationP1::getStringRepr() const
{
  return std::string(REPR);
}

const char *MEDCouplingFieldDiscretizationP1::getRepr() const
{
  return REPR;
}

bool MEDCouplingFieldDiscretizationP1::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  const MEDCouplingFieldDiscretizationP1 *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationP1 *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason="Spatial discrtization of this is ON_NODES, which is not the case of other.";
  return ret;
}

void MEDCouplingFieldDiscretizationP1::checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception)
{
  if(nat!=ConservativeVolumic)
    throw INTERP_KERNEL::Exception("Invalid nature for P1 field  : expected ConservativeVolumic !");
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationP1::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  return mesh->getMeasureFieldOnNode(isAbs);
}

void MEDCouplingFieldDiscretizationP1::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  int id=mesh->getCellContainingPoint(loc,_precision);
  if(id==-1)
    throw INTERP_KERNEL::Exception("Specified point is detected outside of mesh : unable to apply P1::getValueOn !");
  INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(id);
  if(type!=INTERP_KERNEL::NORM_SEG2 && type!=INTERP_KERNEL::NORM_TRI3 && type!=INTERP_KERNEL::NORM_TETRA4)
    throw INTERP_KERNEL::Exception("P1 getValueOn is not specified for not simplex cells !");
  getValueInCell(mesh,id,arr,loc,res);
}

/*!
 * This method localizes a point defined by 'loc' in a cell with id 'cellId' into mesh 'mesh'.
 * The result is put into res expected to be of size at least arr->getNumberOfComponents()
 */
void MEDCouplingFieldDiscretizationP1::getValueInCell(const MEDCouplingMesh *mesh, int cellId, const DataArrayDouble *arr, const double *loc, double *res) const
{
  std::vector<int> conn;
  std::vector<double> coo;
  mesh->getNodeIdsOfCell(cellId,conn);
  for(std::vector<int>::const_iterator iter=conn.begin();iter!=conn.end();iter++)
    mesh->getCoordinatesOfNode(*iter,coo);
  int spaceDim=mesh->getSpaceDimension();
  std::size_t nbOfNodes=conn.size();
  std::vector<const double *> vec(nbOfNodes);
  for(std::size_t i=0;i<nbOfNodes;i++)
    vec[i]=&coo[i*spaceDim];
  INTERP_KERNEL::AutoPtr<double> tmp=new double[nbOfNodes];
  INTERP_KERNEL::barycentric_coords(vec,loc,tmp);
  int sz=arr->getNumberOfComponents();
  INTERP_KERNEL::AutoPtr<double> tmp2=new double[sz];
  std::fill(res,res+sz,0.);
  for(std::size_t i=0;i<nbOfNodes;i++)
    {
      arr->getTuple(conn[i],(double *)tmp2);
      std::transform((double *)tmp2,((double *)tmp2)+sz,(double *)tmp2,std::bind2nd(std::multiplies<double>(),tmp[i]));
      std::transform(res,res+sz,(double *)tmp2,res,std::plus<double>());
    }
}

DataArrayDouble *MEDCouplingFieldDiscretizationP1::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const
{
  std::vector<int> elts,eltsIndex;
  mesh->getCellsContainingPoints(loc,nbOfPoints,_precision,elts,eltsIndex);
  int spaceDim=mesh->getSpaceDimension();
  int nbOfComponents=arr->getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(nbOfPoints,nbOfComponents);
  double *ptToFill=ret->getPointer();
  for(int i=0;i<nbOfPoints;i++)
    if(eltsIndex[i+1]-eltsIndex[i]>=1)
      getValueInCell(mesh,elts[eltsIndex[i]],arr,loc+i*spaceDim,ptToFill+i*nbOfComponents);
    else
      {
        std::ostringstream oss; oss << "Point #" << i << " with coordinates : (";
        std::copy(loc+i*spaceDim,loc+(i+1)*spaceDim,std::ostream_iterator<double>(oss,", "));
        oss << ") detected outside mesh : unable to apply P1::getValueOnMulti ! ";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  ret->incrRef();
  return ret;
}

MEDCouplingFieldDiscretizationPerCell::MEDCouplingFieldDiscretizationPerCell():_discr_per_cell(0)
{
}

MEDCouplingFieldDiscretizationPerCell::~MEDCouplingFieldDiscretizationPerCell()
{
  if(_discr_per_cell)
    _discr_per_cell->decrRef();
}

MEDCouplingFieldDiscretizationPerCell::MEDCouplingFieldDiscretizationPerCell(const MEDCouplingFieldDiscretizationPerCell& other, const int *startCellIds, const int *endCellIds):_discr_per_cell(0)
{
  DataArrayInt *arr=other._discr_per_cell;
  if(arr)
    {
      if(startCellIds==0 && endCellIds==0)
        _discr_per_cell=arr->deepCpy();
      else
        _discr_per_cell=arr->selectByTupleIdSafe(startCellIds,endCellIds);
    }
}

void MEDCouplingFieldDiscretizationPerCell::updateTime() const
{
  if(_discr_per_cell)
    updateTimeWith(*_discr_per_cell);
}

void MEDCouplingFieldDiscretizationPerCell::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception)
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell has no discretization per cell !");
  int nbOfTuples=_discr_per_cell->getNumberOfTuples();
  if(nbOfTuples!=mesh->getNumberOfCells())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell has a discretization per cell but it's not matching the underlying mesh !");
}

bool MEDCouplingFieldDiscretizationPerCell::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  const MEDCouplingFieldDiscretizationPerCell *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationPerCell *>(other);
  if(!otherC)
    {
      reason="Spatial discrtization of this is ON_GAUSS, which is not the case of other.";
      return false;
    }
  if(_discr_per_cell==0)
    return otherC->_discr_per_cell==0;
  if(otherC->_discr_per_cell==0)
    return false;
  bool ret=_discr_per_cell->isEqualIfNotWhy(*otherC->_discr_per_cell,reason);
  if(!ret)
    reason.insert(0,"Field discretization per cell DataArrayInt given the discid per cell :");
  return ret;
}

bool MEDCouplingFieldDiscretizationPerCell::isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const
{
  const MEDCouplingFieldDiscretizationPerCell *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationPerCell *>(other);
  if(!otherC)
    return false;
  if(_discr_per_cell==0)
    return otherC->_discr_per_cell==0;
  if(otherC->_discr_per_cell==0)
    return false;
  return _discr_per_cell->isEqualWithoutConsideringStr(*otherC->_discr_per_cell);
}

/*!
 * This method is typically the first step of renumbering. The impact on _discr_per_cell is necessary here.
 * virtualy by this method.
 */
void MEDCouplingFieldDiscretizationPerCell::renumberCells(const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  int nbCells=_discr_per_cell->getNumberOfTuples();
  const int *array=old2NewBg;
  if(check)
    array=DataArrayInt::CheckAndPreparePermutation(old2NewBg,old2NewBg+nbCells);
  //
  DataArrayInt *dpc=_discr_per_cell->renumber(array);
  _discr_per_cell->decrRef();
  _discr_per_cell=dpc;
  //
  if(check)
    delete [] const_cast<int *>(array);
}

void MEDCouplingFieldDiscretizationPerCell::buildDiscrPerCellIfNecessary(const MEDCouplingMesh *m)
{
  if(!_discr_per_cell)
    {
      _discr_per_cell=DataArrayInt::New();
      int nbTuples=m->getNumberOfCells();
      _discr_per_cell->alloc(nbTuples,1);
      int *ptr=_discr_per_cell->getPointer();
      std::fill(ptr,ptr+nbTuples,DFT_INVALID_LOCID_VALUE);
    }
}

void MEDCouplingFieldDiscretizationPerCell::checkNoOrphanCells() const throw(INTERP_KERNEL::Exception)
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell::checkNoOrphanCells : no discretization defined !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> test=_discr_per_cell->getIdsEqual(DFT_INVALID_LOCID_VALUE);
  if(test->getNumberOfTuples()!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell::checkNoOrphanCells : presence of orphan cells !");
}

const DataArrayInt *MEDCouplingFieldDiscretizationPerCell::getArrayOfDiscIds() const
{
  return _discr_per_cell;
}

MEDCouplingFieldDiscretizationGauss::MEDCouplingFieldDiscretizationGauss()
{
}

MEDCouplingFieldDiscretizationGauss::MEDCouplingFieldDiscretizationGauss(const MEDCouplingFieldDiscretizationGauss& other, const int *startCellIds, const int *endCellIds):MEDCouplingFieldDiscretizationPerCell(other,startCellIds,endCellIds),_loc(other._loc)
{
}

TypeOfField MEDCouplingFieldDiscretizationGauss::getEnum() const
{
  return TYPE;
}

bool MEDCouplingFieldDiscretizationGauss::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  const MEDCouplingFieldDiscretizationGauss *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(other);
  if(!otherC)
    {
      reason="Spatial discrtization of this is ON_GAUSS, which is not the case of other.";
      return false;
    }
  if(!MEDCouplingFieldDiscretizationPerCell::isEqualIfNotWhy(other,eps,reason))
    return false;
  if(_loc.size()!=otherC->_loc.size())
    {
      reason="Gauss spatial discretization : localization sizes differ";
      return false;
    }
  std::size_t sz=_loc.size();
  for(std::size_t i=0;i<sz;i++)
    if(!_loc[i].isEqual(otherC->_loc[i],eps))
      {
        std::ostringstream oss; oss << "Gauss spatial discretization : Localization #" << i << " differ from this to other.";
        reason=oss.str();
        return false;
      }
  return true;
}

bool MEDCouplingFieldDiscretizationGauss::isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const
{
  const MEDCouplingFieldDiscretizationGauss *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(other);
  if(!otherC)
    return false;
  if(!MEDCouplingFieldDiscretizationPerCell::isEqualWithoutConsideringStr(other,eps))
    return false;
  if(_loc.size()!=otherC->_loc.size())
    return false;
  std::size_t sz=_loc.size();
  for(std::size_t i=0;i<sz;i++)
    if(!_loc[i].isEqual(otherC->_loc[i],eps))
      return false;
  return true;
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationGauss::clone() const
{
  return new MEDCouplingFieldDiscretizationGauss(*this);
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationGauss::clonePart(const int *startCellIds, const int *endCellIds) const
{
  return new MEDCouplingFieldDiscretizationGauss(*this,startCellIds,endCellIds);
}

std::string MEDCouplingFieldDiscretizationGauss::getStringRepr() const
{
  std::ostringstream oss; oss << REPR << "." << std::endl;
  if(_discr_per_cell)
    {
      if(_discr_per_cell->isAllocated())
        {
          oss << "Discretization per cell : ";
          std::copy(_discr_per_cell->begin(),_discr_per_cell->end(),std::ostream_iterator<int>(oss,", "));
          oss << std::endl;
        }
    }
  oss << "Presence of " << _loc.size() << " localizations." << std::endl;
  int i=0;
  for(std::vector<MEDCouplingGaussLocalization>::const_iterator it=_loc.begin();it!=_loc.end();it++,i++)
    {
      oss << "+++++ Localization #" << i << " +++++" << std::endl;
      oss << (*it).getStringRepr();
      oss << "++++++++++" << std::endl;
    }
  return oss.str();
}

const char *MEDCouplingFieldDiscretizationGauss::getRepr() const
{
  return REPR;
}

int MEDCouplingFieldDiscretizationGauss::getNumberOfTuples(const MEDCouplingMesh *) const
{
  int ret=0;
  const int *dcPtr=_discr_per_cell->getConstPointer();
  int nbOfTuples=_discr_per_cell->getNumberOfTuples();
  for(const int *w=dcPtr;w!=dcPtr+nbOfTuples;w++)
    ret+=_loc[*w].getNumberOfGaussPt();
  return ret;
}

int MEDCouplingFieldDiscretizationGauss::getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfCells();
}

DataArrayInt *MEDCouplingFieldDiscretizationGauss::getOffsetArr(const MEDCouplingMesh *mesh) const
{
  int nbOfTuples=mesh->getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbOfTuples+1,1);
  int *retPtr=ret->getPointer();
  const int *start=_discr_per_cell->getConstPointer();
  if(_discr_per_cell->getNumberOfTuples()!=nbOfTuples)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getOffsetArr : mismatch between the mesh and the discretization ids array length !");
  int maxPossible=(int)_loc.size();
  retPtr[0]=0;
  for(int i=0;i<nbOfTuples;i++,start++)
    {
      if(*start>=0 && *start<maxPossible)
        retPtr[i+1]=retPtr[i]+_loc[*start].getNumberOfGaussPt();
      else
        {
          std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::getOffsetArr : At position #" << i << " the locid = " << *start << " whereas it should be in [0," << maxPossible << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->incrRef();
  return ret;
}

void MEDCouplingFieldDiscretizationGauss::renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                                                                const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  const int *array=old2NewBg;
  if(check)
    array=DataArrayInt::CheckAndPreparePermutation(old2NewBg,old2NewBg+mesh->getNumberOfCells());
  int nbOfCells=_discr_per_cell->getNumberOfTuples();
  int nbOfTuples=getNumberOfTuples(0);
  const int *dcPtr=_discr_per_cell->getConstPointer();
  int *array2=new int[nbOfTuples];//stores the final conversion array old2New to give to arrays in renumberInPlace.
  int *array3=new int[nbOfCells];//store for each cell in present dcp array (already renumbered) the offset needed by each cell in new numbering.
  array3[0]=0;
  for(int i=1;i<nbOfCells;i++)
    array3[i]=array3[i-1]+_loc[dcPtr[i-1]].getNumberOfGaussPt();
  int j=0;
  for(int i=0;i<nbOfCells;i++)
    {
      int nbOfGaussPt=_loc[dcPtr[array[i]]].getNumberOfGaussPt();
      for(int k=0;k<nbOfGaussPt;k++,j++)
        array2[j]=array3[array[i]]+k;
    }
  delete [] array3;
  for(std::vector<DataArrayDouble *>::const_iterator it=arrays.begin();it!=arrays.end();it++)
    if(*it)
      (*it)->renumberInPlace(array2);
  delete [] array2;
  if(check)
    delete [] const_cast<int*>(array);
}

DataArrayDouble *MEDCouplingFieldDiscretizationGauss::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  checkNoOrphanCells();
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> umesh=mesh->buildUnstructured();//in general do nothing
  int nbOfTuples=getNumberOfTuples(mesh);
  DataArrayDouble *ret=DataArrayDouble::New();
  int spaceDim=mesh->getSpaceDimension();
  ret->alloc(nbOfTuples,spaceDim);
  std::vector< int > locIds;
  std::vector<DataArrayInt *> parts=splitIntoSingleGaussDicrPerCellType(locIds);
  std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > parts2(parts.size());
  std::copy(parts.begin(),parts.end(),parts2.begin());
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> offsets=buildNbOfGaussPointPerCellField();
  offsets->computeOffsets();
  const int *ptrOffsets=offsets->getConstPointer();
  const double *coords=umesh->getCoords()->getConstPointer();
  const int *connI=umesh->getNodalConnectivityIndex()->getConstPointer();
  const int *conn=umesh->getNodalConnectivity()->getConstPointer();
  double *valsToFill=ret->getPointer();
  for(std::size_t i=0;i<parts2.size();i++)
    {
      INTERP_KERNEL::GaussCoords calculator;
      //
      const MEDCouplingGaussLocalization& cli=_loc[locIds[i]];//curLocInfo
      INTERP_KERNEL::NormalizedCellType typ=cli.getType();
      const std::vector<double>& wg=cli.getWeights();
      calculator.addGaussInfo(typ,INTERP_KERNEL::CellModel::GetCellModel(typ).getDimension(),
                                  &cli.getGaussCoords()[0],(int)wg.size(),&cli.getRefCoords()[0],
                                  INTERP_KERNEL::CellModel::GetCellModel(typ).getNumberOfNodes());
      //
      int nbt=parts2[i]->getNumberOfTuples();
      for(const int *w=parts2[i]->getConstPointer();w!=parts2[i]->getConstPointer()+nbt;w++)
        calculator.calculateCoords(cli.getType(),coords,spaceDim,conn+connI[*w]+1,valsToFill+spaceDim*(ptrOffsets[*w]));
    }
  ret->copyStringInfoFrom(*umesh->getCoords());
  return ret;
}

void MEDCouplingFieldDiscretizationGauss::computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *partBg, const int *partEnd,
                                                                             DataArrayInt *&cellRest)
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

/*!
 * Empty : not a bug
 */
void MEDCouplingFieldDiscretizationGauss::checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception)
{
}

void MEDCouplingFieldDiscretizationGauss::getTinySerializationIntInformation(std::vector<int>& tinyInfo) const
{
  int val=-1;
  if(_discr_per_cell)
    val=_discr_per_cell->getNumberOfTuples();
  tinyInfo.push_back(val);
  tinyInfo.push_back((int)_loc.size());
  if(_loc.empty())
    tinyInfo.push_back(-1);
  else
    tinyInfo.push_back(_loc[0].getDimension());
  for(std::vector<MEDCouplingGaussLocalization>::const_iterator iter=_loc.begin();iter!=_loc.end();iter++)
    (*iter).pushTinySerializationIntInfo(tinyInfo);
}

void MEDCouplingFieldDiscretizationGauss::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
  for(std::vector<MEDCouplingGaussLocalization>::const_iterator iter=_loc.begin();iter!=_loc.end();iter++)
    (*iter).pushTinySerializationDblInfo(tinyInfo);
}

void MEDCouplingFieldDiscretizationGauss::getSerializationIntArray(DataArrayInt *& arr) const
{
  arr=0;
  if(_discr_per_cell)
    arr=_discr_per_cell;
}

void MEDCouplingFieldDiscretizationGauss::resizeForUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *& arr)
{
  int val=tinyInfo[0];
  if(val>=0)
    {
      _discr_per_cell=DataArrayInt::New();
      _discr_per_cell->alloc(val,1);
    }
  else
    _discr_per_cell=0;
  arr=_discr_per_cell;
  int nbOfLoc=tinyInfo[1];
  _loc.clear();
  int dim=tinyInfo[2];
  int delta=-1;
  if(nbOfLoc>0)
    delta=((int)tinyInfo.size()-3)/nbOfLoc;
  for(int i=0;i<nbOfLoc;i++)
    {
      std::vector<int> tmp(tinyInfo.begin()+3+i*delta,tinyInfo.begin()+3+(i+1)*delta);
      MEDCouplingGaussLocalization elt=MEDCouplingGaussLocalization::BuildNewInstanceFromTinyInfo(dim,tmp);
      _loc.push_back(elt);
    }
}

void MEDCouplingFieldDiscretizationGauss::finishUnserialization(const std::vector<double>& tinyInfo)
{
  double *tmp=new double[tinyInfo.size()];
  std::copy(tinyInfo.begin(),tinyInfo.end(),tmp);
  const double *work=tmp;
  for(std::vector<MEDCouplingGaussLocalization>::iterator iter=_loc.begin();iter!=_loc.end();iter++)
    work=(*iter).fillWithValues(work);
  delete [] tmp;
}

double MEDCouplingFieldDiscretizationGauss::getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da,
                                                   int cellId, int nodeIdInCell, int compoId) const throw(INTERP_KERNEL::Exception)
{
  int offset=getOffsetOfCell(cellId);
  return da->getIJ(offset+nodeIdInCell,compoId);
}

void MEDCouplingFieldDiscretizationGauss::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingFieldDiscretizationPerCell::checkCoherencyBetween(mesh,da);
  for(std::vector<MEDCouplingGaussLocalization>::const_iterator iter=_loc.begin();iter!=_loc.end();iter++)
    (*iter).checkCoherency();
  int nbOfDesc=(int)_loc.size();
  int nbOfCells=mesh->getNumberOfCells();
  const int *dc=_discr_per_cell->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      if(dc[i]>=nbOfDesc)
        {
          std::ostringstream oss; oss << "Cell # " << i << " of mesh \"" << mesh->getName() << "\" has an undefined gauss location ! Should never happend !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(dc[i]<0)
        {
          std::ostringstream oss; oss << "Cell # " << i << " of mesh \"" << mesh->getName() << "\" has no gauss location !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      if(mesh->getTypeOfCell(i)!=_loc[dc[i]].getType())
        {
          std::ostringstream oss; oss << "Types of mesh and gauss location mismatch for cell # " << i;
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  int nbOfTuples=getNumberOfTuples(mesh);
  if(nbOfTuples!=da->getNumberOfTuples())
    {
      std::ostringstream oss; oss << "Invalid number of tuples in the array : expecting " << nbOfTuples << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationGauss::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGauss::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGauss::getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const
{
  throw INTERP_KERNEL::Exception("getValueOnPos(i,j,k) : Not applyable for Gauss points !");
}

DataArrayDouble *MEDCouplingFieldDiscretizationGauss::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const
{
  throw INTERP_KERNEL::Exception("getValueOnMulti : Not implemented yet for gauss points !");
}

MEDCouplingMesh *MEDCouplingFieldDiscretizationGauss::buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const
{
  di=computeTupleIdsToSelectFromCellIds(mesh,start,end);
  return mesh->buildPart(start,end);
}

/*!
 * This method returns a tuple ids selection from cell ids selection [start;end).
 * This method is called by MEDCouplingFieldDiscretizationGauss::buildSubMeshData to return parameter \b di.
 *
 * \return a newly allocated array containing ids to select into the DataArrayDouble of the field.
 * 
 */
DataArrayInt *MEDCouplingFieldDiscretizationGauss::computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const int *startCellIds, const int *endCellIds) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::computeTupleIdsToSelectFromCellIds : null mesh !");
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::computeTupleIdsToSelectFromCellIds : null discretization ids !");
  int nbOfCells=mesh->getNumberOfCells();
  if(_discr_per_cell->getNumberOfTuples()!=nbOfCells)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::computeTupleIdsToSelectFromCellIds : mismatch of nb of tuples of cell ids array and number of cells !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> nbOfNodesPerCell=DataArrayInt::New(); nbOfNodesPerCell->alloc(nbOfCells,1);
  int *retPtr=nbOfNodesPerCell->getPointer();
  const int *pt=_discr_per_cell->getConstPointer();
  int nbMaxOfLocId=(int)_loc.size();
  for(int i=0;i<nbOfCells;i++,retPtr++,pt++)
    {
      if(*pt>=0 && *pt<nbMaxOfLocId)
        *retPtr=_loc[*pt].getNumberOfGaussPt();
    }
  nbOfNodesPerCell->computeOffsets2();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> sel=DataArrayInt::New(); sel->useArray(startCellIds,false,CPP_DEALLOC,(int)std::distance(startCellIds,endCellIds),1);
  return sel->buildExplicitArrByRanges(nbOfNodesPerCell);
}

/*!
 * No implementation needed !
 */
void MEDCouplingFieldDiscretizationGauss::renumberValuesOnNodes(double , const int *, int newNbOfNodes, DataArrayDouble *) const
{
}

void MEDCouplingFieldDiscretizationGauss::renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const int *old2New, int newSz, DataArrayDouble *arr) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGauss::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const
{
  throw INTERP_KERNEL::Exception("Number of cells has changed and becomes higher with some cells that have been split ! Unable to conserve the Gauss field !");
}

void MEDCouplingFieldDiscretizationGauss::setGaussLocalizationOnType(const MEDCouplingMesh *m, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                                     const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if((int)cm.getDimension()!=m->getMeshDimension())
    {
      std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::setGaussLocalizationOnType : mismatch of dimensions ! MeshDim==" << m->getMeshDimension();
      oss << " whereas Type '" << cm.getRepr() << "' has dimension " << cm.getDimension() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  buildDiscrPerCellIfNecessary(m);
  int id=(int)_loc.size();
  MEDCouplingGaussLocalization elt(type,refCoo,gsCoo,wg);
  _loc.push_back(elt);
  int *ptr=_discr_per_cell->getPointer();
  int nbCells=m->getNumberOfCells();
  for(int i=0;i<nbCells;i++)
    if(m->getTypeOfCell(i)==type)
      ptr[i]=id;
  zipGaussLocalizations();
}

void MEDCouplingFieldDiscretizationGauss::setGaussLocalizationOnCells(const MEDCouplingMesh *m, const int *begin, const int *end, const std::vector<double>& refCoo,
                                                                      const std::vector<double>& gsCoo, const std::vector<double>& wg) throw(INTERP_KERNEL::Exception)
{
  buildDiscrPerCellIfNecessary(m);
  if(std::distance(begin,end)<1)
    throw INTERP_KERNEL::Exception("Size of [begin,end) must be equal or greater than 1 !");
  INTERP_KERNEL::NormalizedCellType type=m->getTypeOfCell(*begin);
  MEDCouplingGaussLocalization elt(type,refCoo,gsCoo,wg);
  int id=(int)_loc.size();
  int *ptr=_discr_per_cell->getPointer();
  for(const int *w=begin+1;w!=end;w++)
    {
      if(m->getTypeOfCell(*w)!=type)
        {
          std::ostringstream oss; oss << "The cell with id " << *w << " has been detected to be incompatible in the [begin,end) array specified !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  //
  for(const int *w2=begin;w2!=end;w2++)
    ptr[*w2]=id;
  //
  _loc.push_back(elt);
  zipGaussLocalizations();
}

void MEDCouplingFieldDiscretizationGauss::clearGaussLocalizations() throw(INTERP_KERNEL::Exception)
{
  if(_discr_per_cell)
    {
      _discr_per_cell->decrRef();
      _discr_per_cell=0;
    }
  _loc.clear();
}

MEDCouplingGaussLocalization& MEDCouplingFieldDiscretizationGauss::getGaussLocalization(int locId) throw(INTERP_KERNEL::Exception)
{
  checkLocalizationId(locId);
  return _loc[locId];
}

int MEDCouplingFieldDiscretizationGauss::getNbOfGaussLocalization() const throw(INTERP_KERNEL::Exception)
{
  return (int)_loc.size();
}

int MEDCouplingFieldDiscretizationGauss::getGaussLocalizationIdOfOneCell(int cellId) const throw(INTERP_KERNEL::Exception)
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("No Gauss localization still set !");
  int locId=_discr_per_cell->getConstPointer()[cellId];
  if(locId<0)
    throw INTERP_KERNEL::Exception("No Gauss localization set for the specified cell !");
  return locId;
}

int MEDCouplingFieldDiscretizationGauss::getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
{
  std::set<int> ret=getGaussLocalizationIdsOfOneType(type);
  if(ret.empty())
    throw INTERP_KERNEL::Exception("No gauss discretization found for the specified type !");
  if(ret.size()>1)
    throw INTERP_KERNEL::Exception("Several gauss discretizations have been found for the specified type !");
  return *ret.begin();
}

std::set<int> MEDCouplingFieldDiscretizationGauss::getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("No Gauss localization still set !");
  std::set<int> ret;
  int id=0;
  for(std::vector<MEDCouplingGaussLocalization>::const_iterator iter=_loc.begin();iter!=_loc.end();iter++,id++)
    if((*iter).getType()==type)
      ret.insert(id);
  return ret;
}

void MEDCouplingFieldDiscretizationGauss::getCellIdsHavingGaussLocalization(int locId, std::vector<int>& cellIds) const throw(INTERP_KERNEL::Exception)
{
  if(locId<0 || locId>=(int)_loc.size())
    throw INTERP_KERNEL::Exception("Invalid locId given : must be in range [0:getNbOfGaussLocalization()) !");
  int nbOfTuples=_discr_per_cell->getNumberOfTuples();
  const int *ptr=_discr_per_cell->getConstPointer();
  for(int i=0;i<nbOfTuples;i++)
    if(ptr[i]==locId)
      cellIds.push_back(i);
}

const MEDCouplingGaussLocalization& MEDCouplingFieldDiscretizationGauss::getGaussLocalization(int locId) const throw(INTERP_KERNEL::Exception)
{
  checkLocalizationId(locId);
  return _loc[locId];
}

void MEDCouplingFieldDiscretizationGauss::checkLocalizationId(int locId) const throw(INTERP_KERNEL::Exception)
{
  if(locId<0 || locId>=(int)_loc.size())
    throw INTERP_KERNEL::Exception("Invalid locId given : must be in range [0:getNbOfGaussLocalization()) !");
}

int MEDCouplingFieldDiscretizationGauss::getOffsetOfCell(int cellId) const throw(INTERP_KERNEL::Exception)
{
  int ret=0;
  const int *start=_discr_per_cell->getConstPointer();
  for(const int *w=start;w!=start+cellId;w++)
    ret+=_loc[*w].getNumberOfGaussPt();
  return ret;
}

/*!
 * This method do the assumption that there is no orphan cell. If there is an exception is thrown.
 * This method makes the assumption too that '_discr_per_cell' is defined. If not an exception is thrown.
 * This method returns a newly created array with number of tuples equals to '_discr_per_cell->getNumberOfTuples' and number of components equal to 1.
 * The i_th tuple in returned array is the number of gauss point if the corresponding cell.
 */
DataArrayInt *MEDCouplingFieldDiscretizationGauss::buildNbOfGaussPointPerCellField() const throw(INTERP_KERNEL::Exception)
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::buildNbOfGaussPointPerCellField : no discretization array set !");
  int nbOfTuples=_discr_per_cell->getNumberOfTuples();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  const int *w=_discr_per_cell->getConstPointer();
  ret->alloc(nbOfTuples,1);
  int *valsToFill=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++,w++)
    if(*w!=DFT_INVALID_LOCID_VALUE)
      valsToFill[i]=_loc[*w].getNumberOfGaussPt();
    else
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::buildNbOfGaussPointPerCellField : orphan cell detected !");
  ret->incrRef();
  return ret;
}

/*!
 * This method makes the assumption that _discr_per_cell is set.
 * This method reduces as much as possible number size of _loc.
 * This method is usefull when several set on same cells has been done and that some Gauss Localization are no more used.
 */
void MEDCouplingFieldDiscretizationGauss::zipGaussLocalizations()
{
  const int *start=_discr_per_cell->getConstPointer();
  int nbOfTuples=_discr_per_cell->getNumberOfTuples();
  INTERP_KERNEL::AutoPtr<int> tmp=new int[_loc.size()];
  std::fill((int *)tmp,(int *)tmp+_loc.size(),-2);
  for(const int *w=start;w!=start+nbOfTuples;w++)
    if(*w>=0)
      tmp[*w]=1;
  int fid=0;
  for(int i=0;i<(int)_loc.size();i++)
    if(tmp[i]!=-2)
      tmp[i]=fid++;
  if(fid==(int)_loc.size())
    return;
  // zip needed
  int *start2=_discr_per_cell->getPointer();
  for(int *w2=start2;w2!=start2+nbOfTuples;w2++)
    if(*w2>=0)
      *w2=tmp[*w2];
  std::vector<MEDCouplingGaussLocalization> tmpLoc;
  for(int i=0;i<(int)_loc.size();i++)
    if(tmp[i]!=-2)
      tmpLoc.push_back(_loc[tmp[i]]);
  _loc=tmpLoc;
}

/*!
 * This method is usefull when 'this' describes a field discretization with several gauss discretization on a \b same cell type.
 * For example same NORM_TRI3 cells having 6 gauss points and others with 12 gauss points.
 * This method returns 2 arrays with same size : the return value and 'locIds' output parameter.
 * For a given i into [0,locIds.size) ret[i] represents the set of cell ids of i_th set an locIds[i] represents the set of discretisation of the set.
 * The return vector contains a set of newly created instance to deal with.
 * The returned vector represents a \b partition of cells ids with a gauss discretization set.
 * 
 * If no descretization is set in 'this' and exception will be thrown.
 */
std::vector<DataArrayInt *> MEDCouplingFieldDiscretizationGauss::splitIntoSingleGaussDicrPerCellType(std::vector<int>& locIds) const throw(INTERP_KERNEL::Exception)
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::splitIntoSingleGaussDicrPerCellType : no descretization set !");
  return _discr_per_cell->partitionByDifferentValues(locIds);
}

MEDCouplingFieldDiscretizationGaussNE::MEDCouplingFieldDiscretizationGaussNE()
{
}

TypeOfField MEDCouplingFieldDiscretizationGaussNE::getEnum() const
{
  return TYPE;
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationGaussNE::clone() const
{
  return new MEDCouplingFieldDiscretizationGaussNE(*this);
}

std::string MEDCouplingFieldDiscretizationGaussNE::getStringRepr() const
{
  return std::string(REPR);
}

const char *MEDCouplingFieldDiscretizationGaussNE::getRepr() const
{
  return REPR;
}

bool MEDCouplingFieldDiscretizationGaussNE::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  const MEDCouplingFieldDiscretizationGaussNE *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationGaussNE *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason="Spatial discrtization of this is ON_GAUSS_NE, which is not the case of other.";
  return ret;
}

int MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  int ret=0;
  int nbOfCells=mesh->getNumberOfCells();
  for(int i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      if(cm.isDynamic())
        throw INTERP_KERNEL::Exception("Not implemented yet Gauss node on elements for polygons and polyedrons !");
      ret+=cm.getNumberOfNodes();
    }
  return ret;
}

int MEDCouplingFieldDiscretizationGaussNE::getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfCells();
}

DataArrayInt *MEDCouplingFieldDiscretizationGaussNE::getOffsetArr(const MEDCouplingMesh *mesh) const
{
  int nbOfTuples=mesh->getNumberOfCells();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(nbOfTuples+1,1);
  int *retPtr=ret->getPointer();
  retPtr[0]=0;
  for(int i=0;i<nbOfTuples;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      if(cm.isDynamic())
        throw INTERP_KERNEL::Exception("Not implemented yet Gauss node on elements for polygons and polyedrons !");
      retPtr[i+1]=retPtr[i]+cm.getNumberOfNodes();
    }
  return ret;
}

void MEDCouplingFieldDiscretizationGaussNE::renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                                                                  const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  const int *array=old2NewBg;
  if(check)
    array=DataArrayInt::CheckAndPreparePermutation(old2NewBg,old2NewBg+mesh->getNumberOfCells());
  int nbOfCells=mesh->getNumberOfCells();
  int nbOfTuples=getNumberOfTuples(mesh);
  int *array2=new int[nbOfTuples];//stores the final conversion array old2New to give to arrays in renumberInPlace.
  int *array3=new int[nbOfCells];//store for each cell in after renumbering the offset needed by each cell in new numbering.
  array3[0]=0;
  for(int i=1;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell((int)std::distance(array,std::find(array,array+nbOfCells,i-1)));
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      array3[i]=array3[i-1]+cm.getNumberOfNodes();
    }
  int j=0;
  for(int i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      for(int k=0;k<(int)cm.getNumberOfNodes();k++,j++)
        array2[j]=array3[array[i]]+k;
    }
  delete [] array3;
  for(std::vector<DataArrayDouble *>::const_iterator it=arrays.begin();it!=arrays.end();it++)
    if(*it)
      (*it)->renumberInPlace(array2);
  delete [] array2;
  if(check)
    delete [] const_cast<int *>(array);
}

DataArrayDouble *MEDCouplingFieldDiscretizationGaussNE::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGaussNE::computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *partBg, const int *partEnd,
                                                                               DataArrayInt *&cellRest)
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGaussNE::checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception)
{
}

double MEDCouplingFieldDiscretizationGaussNE::getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da,
                                                     int cellId, int nodeIdInCell, int compoId) const throw(INTERP_KERNEL::Exception)
{
  int offset=0;
  for(int i=0;i<cellId;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      offset+=cm.getNumberOfNodes();
    }
  return da->getIJ(offset+nodeIdInCell,compoId);
}

void MEDCouplingFieldDiscretizationGaussNE::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception)
{
  int nbOfTuples=getNumberOfTuples(mesh);
  if(nbOfTuples!=da->getNumberOfTuples())
    {
      std::ostringstream oss; oss << "Invalid number of tuples in the array : expecting " << nbOfTuples << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationGaussNE::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGaussNE::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGaussNE::getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const
{
  throw INTERP_KERNEL::Exception("getValueOnPos(i,j,k) : Not applyable for Gauss points !");
}

DataArrayDouble *MEDCouplingFieldDiscretizationGaussNE::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfPoints) const
{
  throw INTERP_KERNEL::Exception("getValueOnMulti : Not implemented for Gauss NE !");
}

MEDCouplingMesh *MEDCouplingFieldDiscretizationGaussNE::buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const
{
  di=computeTupleIdsToSelectFromCellIds(mesh,start,end);
  return mesh->buildPart(start,end);
}

/*!
 * This method returns a tuple ids selection from cell ids selection [start;end).
 * This method is called by MEDCouplingFieldDiscretizationGaussNE::buildSubMeshData to return parameter \b di.
 *
 * \return a newly allocated array containing ids to select into the DataArrayDouble of the field.
 * 
 */
DataArrayInt *MEDCouplingFieldDiscretizationGaussNE::computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const int *startCellIds, const int *endCellIds) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::computeTupleIdsToSelectFromCellIds : null mesh !");
  const MEDCouplingAutoRefCountObjectPtr<MEDCouplingUMesh> umesh=mesh->buildUnstructured();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> nbOfNodesPerCell=umesh->computeNbOfNodesPerCell();
  nbOfNodesPerCell->computeOffsets2();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> sel=DataArrayInt::New(); sel->useArray(startCellIds,false,CPP_DEALLOC,(int)std::distance(startCellIds,endCellIds),1);
  return sel->buildExplicitArrByRanges(nbOfNodesPerCell);
}

/*!
 * No implementation needed !
 */
void MEDCouplingFieldDiscretizationGaussNE::renumberValuesOnNodes(double , const int *, int newNbOfNodes, DataArrayDouble *) const
{
}

void MEDCouplingFieldDiscretizationGaussNE::renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const int *old2New, int newSz, DataArrayDouble *arr) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGaussNE::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

MEDCouplingFieldDiscretizationGaussNE::MEDCouplingFieldDiscretizationGaussNE(const MEDCouplingFieldDiscretizationGaussNE& other):MEDCouplingFieldDiscretization(other)
{
}

TypeOfField MEDCouplingFieldDiscretizationKriging::getEnum() const
{
  return TYPE;
}

const char *MEDCouplingFieldDiscretizationKriging::getRepr() const
{
  return REPR;
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationKriging::clone() const
{
  return new MEDCouplingFieldDiscretizationKriging;
}

std::string MEDCouplingFieldDiscretizationKriging::getStringRepr() const
{
  return std::string(REPR);
}

void MEDCouplingFieldDiscretizationKriging::checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception)
{
  if(nat!=ConservativeVolumic)
    throw INTERP_KERNEL::Exception("Invalid nature for Kriging field : expected ConservativeVolumic !");
}

bool MEDCouplingFieldDiscretizationKriging::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  const MEDCouplingFieldDiscretizationKriging *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationKriging *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason="Spatial discrtization of this is ON_NODES_KR, which is not the case of other.";
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationKriging::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  throw INTERP_KERNEL::Exception("getMeasureField on FieldDiscretizationKriging : not implemented yet !");
}

void MEDCouplingFieldDiscretizationKriging::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> res2=MEDCouplingFieldDiscretizationKriging::getValueOnMulti(arr,mesh,loc,1);
  std::copy(res2->begin(),res2->end(),res);
}

DataArrayDouble *MEDCouplingFieldDiscretizationKriging::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, int nbOfTargetPoints) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coords=getLocalizationOfDiscValues(mesh);
  int nbOfPts=coords->getNumberOfTuples();
  int dimension=coords->getNumberOfComponents();
  //
  int delta=0;
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> KnewiK=computeVectorOfCoefficients(mesh,arr,delta);
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> locArr=DataArrayDouble::New();
  locArr->useArray(loc,false,CPP_DEALLOC,nbOfTargetPoints,dimension);
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> matrix2=coords->buildEuclidianDistanceDenseMatrixWith(locArr);
  operateOnDenseMatrix(mesh->getSpaceDimension(),nbOfPts*nbOfTargetPoints,matrix2->getPointer());
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> matrix3=DataArrayDouble::New();
  matrix3->alloc((nbOfPts+delta)*nbOfTargetPoints,1);
  double *work=matrix3->getPointer();
  const double *workCst=matrix2->getConstPointer();
  const double *workCst2=loc;
  for(int i=0;i<nbOfTargetPoints;i++,workCst+=nbOfPts,workCst2+=delta-1)
    {
      for(int j=0;j<nbOfPts;j++)
        work[j*nbOfTargetPoints+i]=workCst[j];
      work[nbOfPts*nbOfTargetPoints+i]=1.0;
      for(int j=0;j<delta-1;j++)
        work[(nbOfPts+1+j)*nbOfTargetPoints+i]=workCst2[j];
    }
  //
  int nbOfCompo=arr->getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(nbOfTargetPoints,nbOfCompo);
  INTERP_KERNEL::matrixProduct(KnewiK->getConstPointer(),1,nbOfPts+delta,matrix3->getConstPointer(),nbOfPts+delta,nbOfTargetPoints*nbOfCompo,ret->getPointer());
  ret->incrRef();
  return ret;
}

/*!
 * This method computes coefficients to apply to each representing points of \a mesh, that is to say the nodes of \a mesh given a field array \a arr whose
 * number of tuples should be equal to the number of representing points in \a mesh.
 * 
 * \param [in] mesh is the sources of nodes on which kriging will be done regarding the parameters and the value of \c this->getSpaceDimension()
 * \param [in] arr input field DataArrayDouble whose number of tuples must be equal to the number of nodes in \a mesh
 * \param [out] isDrift return if drift coefficients are present in the returned vector of coefficients, and if. If different from 0 there is presence of drift coefficients.
 *              Whatever the value of \a isDrift the number of tuples of returned DataArrayDouble  will be equal to \c arr->getNumberOfTuples() + \a isDrift.
 * \return a newly allocated array containing coefficients including or not drift coefficient at the end depending the value of \a isDrift parameter.
 */
DataArrayDouble *MEDCouplingFieldDiscretizationKriging::computeVectorOfCoefficients(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, int& isDrift) const
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coords=getLocalizationOfDiscValues(mesh);
  int nbOfPts=coords->getNumberOfTuples();
  //int dimension=coords->getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> matrix=coords->buildEuclidianDistanceDenseMatrix();
  operateOnDenseMatrix(mesh->getSpaceDimension(),nbOfPts*nbOfPts,matrix->getPointer());
  // Drift
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> matrixWithDrift=performDrift(matrix,coords,isDrift);
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> matrixInv=DataArrayDouble::New();
  matrixInv->alloc((nbOfPts+isDrift)*(nbOfPts+isDrift),1);
  INTERP_KERNEL::inverseMatrix(matrixWithDrift->getConstPointer(),nbOfPts+isDrift,matrixInv->getPointer());
  //
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> KnewiK=DataArrayDouble::New();
  KnewiK->alloc((nbOfPts+isDrift)*1,1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arr2=DataArrayDouble::New();
  arr2->alloc((nbOfPts+isDrift)*1,1);
  double *work=std::copy(arr->begin(),arr->end(),arr2->getPointer());
  std::fill(work,work+isDrift,0.);
  INTERP_KERNEL::matrixProduct(matrixInv->getConstPointer(),nbOfPts+isDrift,nbOfPts+isDrift,arr2->getConstPointer(),nbOfPts+isDrift,1,KnewiK->getPointer());
  KnewiK->incrRef();
  return KnewiK;
}

/*!
 * Apply \f f(x) on each element x in \a matrixPtr. \a matrixPtr is expected to be a dense matrix represented by a chunck of memory of size at least equal to \a nbOfElems.
 *
 * \param [in] spaceDimension space dimension of the input mesh on which the Kriging has to be performed
 * \param [in] nbOfElems is the result of the product of nb of rows and the nb of columns of matrix \a matrixPtr
 * \param [in,out] matrixPtr is the dense matrix whose on each values the operation will be applied
 */
void MEDCouplingFieldDiscretizationKriging::operateOnDenseMatrix(int spaceDimension, int nbOfElems, double *matrixPtr) const
{
  switch(spaceDimension)
    {
    case 1:
      {
        for(int i=0;i<nbOfElems;i++)
          {
            double val=matrixPtr[i];
            matrixPtr[i]=val*val*val;
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::operateOnDenseMatrix : only dimension 1 implemented !");
    }
}

/*!
 * Starting from a square matrix \a matr, this method returns a newly allocated dense square matrix whose \a matr is included in returned matrix
 * in the top left corner, and in the remaining returned matrix the parameters to take into account about the kriging drift.
 * For the moment only linear srift is implemented.
 *
 * \param [in] arr the position of points were input mesh geometry is considered for Kriging
 * \param [in] matr input matrix whose drift part will be added
 * \param [out] delta the difference between the size of the output matrix and the input matrix \a matr.
 * \return a newly allocated matrix bigger than input matrix \a matr.
 */
DataArrayDouble *MEDCouplingFieldDiscretizationKriging::performDrift(const DataArrayDouble *matr, const DataArrayDouble *arr, int& delta) const
{
  int spaceDimension=arr->getNumberOfComponents();
  delta=spaceDimension+1;
  int szOfMatrix=arr->getNumberOfTuples();
  if(szOfMatrix*szOfMatrix!=matr->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::performDrift : invalid size");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc((szOfMatrix+delta)*(szOfMatrix+delta),1);
  const double *srcWork=matr->getConstPointer();
  const double *srcWork2=arr->getConstPointer();
  double *destWork=ret->getPointer();
  for(int i=0;i<szOfMatrix;i++)
    {
      destWork=std::copy(srcWork,srcWork+szOfMatrix,destWork);
      srcWork+=szOfMatrix;
      *destWork++=1.;
      destWork=std::copy(srcWork2,srcWork2+spaceDimension,destWork);
      srcWork2+=spaceDimension;
    }
  std::fill(destWork,destWork+szOfMatrix,1.); destWork+=szOfMatrix;
  std::fill(destWork,destWork+spaceDimension+1,0.); destWork+=spaceDimension+1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> arrNoI=arr->toNoInterlace();
  srcWork2=arrNoI->getConstPointer();
  for(int i=0;i<spaceDimension;i++)
    {
      destWork=std::copy(srcWork2,srcWork2+szOfMatrix,destWork);
      srcWork2+=szOfMatrix;
      std::fill(destWork,destWork+spaceDimension+1,0.);
      destWork+=spaceDimension;
    }
  //
  ret->incrRef();
  return ret;
}

