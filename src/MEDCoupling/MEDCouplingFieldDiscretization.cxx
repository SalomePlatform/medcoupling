//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "CellModel.hxx"

#include "InterpolationUtils.hxx"

#include <set>
#include <limits>
#include <algorithm>
#include <functional>

using namespace ParaMEDMEM;

const double MEDCouplingFieldDiscretization::DFLT_PRECISION=1.e-12;

const char MEDCouplingFieldDiscretizationP0::REPR[]="P0";

const TypeOfField MEDCouplingFieldDiscretizationP0::TYPE=ON_CELLS;

const char MEDCouplingFieldDiscretizationP1::REPR[]="P1";

const TypeOfField MEDCouplingFieldDiscretizationP1::TYPE=ON_NODES;

const char MEDCouplingFieldDiscretizationGauss::REPR[]="GAUSS";

const TypeOfField MEDCouplingFieldDiscretizationGauss::TYPE=ON_GAUSS_PT;

const char MEDCouplingFieldDiscretizationGaussNE::REPR[]="GSSNE";

const TypeOfField MEDCouplingFieldDiscretizationGaussNE::TYPE=ON_GAUSS_NE;

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
  throw INTERP_KERNEL::Exception("Representation does not match with any field discretization !");
}

bool MEDCouplingFieldDiscretization::isEqualWithoutConsideringStr(const MEDCouplingFieldDiscretization *other, double eps) const
{
  return isEqual(other,eps);
}

/*!
 * Excepted for MEDCouplingFieldDiscretizationPerCell no underlying TimeLabel object : nothing to do in generally.
 */
void MEDCouplingFieldDiscretization::updateTime()
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
  for(int i=0;i<nbOfElems;i++)
    {
      for(int j=0;j<nbOfCompo;j++)
        res[j]+=fabs(arrPtr[i*nbOfCompo+j])*volPtr[i];
    }
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
  for(int i=0;i<nbOfElems;i++)
    {
      for(int j=0;j<nbOfCompo;j++)
        res[j]+=arrPtr[i*nbOfCompo+j]*arrPtr[i*nbOfCompo+j]*fabs(volPtr[i]);
    }
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

void MEDCouplingFieldDiscretization::getCellIdsHavingGaussLocalization(int locId, std::vector<int>& cellIds) const throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

void MEDCouplingFieldDiscretization::renumberEntitiesFromO2NArr(const int *old2NewPtr, DataArrayDouble *arr, const char *msg)
{
  int oldNbOfElems=arr->getNumberOfTuples();
  int nbOfComp=arr->getNumberOfComponents();
  int newNbOfTuples=(*std::max_element(old2NewPtr,old2NewPtr+oldNbOfElems))+1;
  DataArrayDouble *arrCpy=arr->deepCpy();
  const double *ptSrc=arrCpy->getConstPointer();
  arr->reAlloc(newNbOfTuples);
  double *ptToFill=arr->getPointer();
  std::fill(ptToFill,ptToFill+nbOfComp*newNbOfTuples,std::numeric_limits<double>::max());
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
              if(!std::equal(ptSrc+i*nbOfComp,ptSrc+(i+1)*nbOfComp,ptToFill+newNb*nbOfComp))
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

void MEDCouplingFieldDiscretization::renumberEntitiesFromN2OArr(const int *new2OldPtr, int new2OldSz, DataArrayDouble *arr, const char *msg)
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

const char *MEDCouplingFieldDiscretizationP0::getStringRepr() const
{
  return REPR;
}

bool MEDCouplingFieldDiscretizationP0::isEqual(const MEDCouplingFieldDiscretization *other, double eps) const
{
  const MEDCouplingFieldDiscretizationP0 *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationP0 *>(other);
  return otherC!=0;
}

int MEDCouplingFieldDiscretizationP0::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfCells();
}

void MEDCouplingFieldDiscretizationP0::renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                                                             const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  const int *array=old2NewBg;
  if(check)
    array=DataArrayInt::checkAndPreparePermutation(old2NewBg,old2NewBg+mesh->getNumberOfCells());
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
  cellRest->alloc(std::distance(partBg,partEnd),1);
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

/*!
 * Nothing to do. It's not a bug.
 */
void MEDCouplingFieldDiscretizationP0::renumberValuesOnNodes(const int *, DataArrayDouble *) const
{
}

void MEDCouplingFieldDiscretizationP0::renumberValuesOnCells(const MEDCouplingMesh *mesh, const int *old2New, DataArrayDouble *arr) const
{
  renumberEntitiesFromO2NArr(old2New,arr,"Cell");
}

void MEDCouplingFieldDiscretizationP0::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const
{
  renumberEntitiesFromN2OArr(new2old,newSz,arr,"Cell");
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
  di->alloc(std::distance(start,end),1);
  int *pt=di->getPointer();
  std::copy(start,end,pt);
  return ret;
}

TypeOfField MEDCouplingFieldDiscretizationP1::getEnum() const
{
  return TYPE;
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationP1::clone() const
{
  return new MEDCouplingFieldDiscretizationP1;
}

const char *MEDCouplingFieldDiscretizationP1::getStringRepr() const
{
  return REPR;
}

bool MEDCouplingFieldDiscretizationP1::isEqual(const MEDCouplingFieldDiscretization *other, double eps) const
{
  const MEDCouplingFieldDiscretizationP1 *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationP1 *>(other);
  return otherC!=0;
}

/*!
 * Nothing to do here.
 */
void MEDCouplingFieldDiscretizationP1::renumberArraysForCell(const MEDCouplingMesh *, const std::vector<DataArrayDouble *>& arrays,
                                                             const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
}

int MEDCouplingFieldDiscretizationP1::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfNodes();
}

DataArrayDouble *MEDCouplingFieldDiscretizationP1::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  return mesh->getCoordinatesAndOwner();
}

void MEDCouplingFieldDiscretizationP1::computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const int *partBg, const int *partEnd,
                                                                          DataArrayInt *&cellRest)
{
  cellRest=mesh->getCellIdsFullyIncludedInNodeIds(partBg,partEnd);
}

void MEDCouplingFieldDiscretizationP1::checkCompatibilityWithNature(NatureOfField nat) const throw(INTERP_KERNEL::Exception)
{
  if(nat!=ConservativeVolumic)
    throw INTERP_KERNEL::Exception("Invalid nature for P1 field !");
}

void MEDCouplingFieldDiscretizationP1::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArrayDouble *da) const throw(INTERP_KERNEL::Exception)
{
  if(mesh->getNumberOfNodes()!=da->getNumberOfTuples())
    {
      std::ostringstream message;
      message << "Field on nodes invalid because there are " << mesh->getNumberOfNodes();
      message << " nodes in mesh and " << da->getNumberOfTuples() << " tuples in field !";
      throw INTERP_KERNEL::Exception(message.str().c_str());
    }
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
  std::vector<int> conn;
  std::vector<double> coo;
  mesh->getNodeIdsOfCell(id,conn);
  for(std::vector<int>::const_iterator iter=conn.begin();iter!=conn.end();iter++)
    mesh->getCoordinatesOfNode(*iter,coo);
  int spaceDim=mesh->getSpaceDimension();
  int nbOfNodes=conn.size();
  std::vector<const double *> vec(nbOfNodes);
  for(int i=0;i<nbOfNodes;i++)
    vec[i]=&coo[i*spaceDim];
  double *tmp=new double[nbOfNodes];
  INTERP_KERNEL::barycentric_coords(vec,loc,tmp);
  int sz=arr->getNumberOfComponents();
  double *tmp2=new double[sz];
  std::fill(res,res+sz,0.);
  for(int i=0;i<nbOfNodes;i++)
    {
      arr->getTuple(conn[i],tmp2);
      std::transform(tmp2,tmp2+sz,tmp2,std::bind2nd(std::multiplies<double>(),tmp[i]));
      std::transform(res,res+sz,tmp2,res,std::plus<double>());
    }
  delete [] tmp;
  delete [] tmp2;
}

void MEDCouplingFieldDiscretizationP1::getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, int i, int j, int k, double *res) const
{
  const MEDCouplingCMesh *meshC=dynamic_cast<const MEDCouplingCMesh *>(mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("P1::getValueOnPos is only accessible for structured meshes !");
  int id=meshC->getNodeIdFromPos(i,j,k);
  arr->getTuple(id,res);
}

void MEDCouplingFieldDiscretizationP1::renumberValuesOnNodes(const int *old2NewPtr, DataArrayDouble *arr) const
{
  renumberEntitiesFromO2NArr(old2NewPtr,arr,"Node");
}

/*!
 * Nothing to do it's not a bug.
 */
void MEDCouplingFieldDiscretizationP1::renumberValuesOnCells(const MEDCouplingMesh *mesh, const int *old2New, DataArrayDouble *arr) const
{
}

/*!
 * Nothing to do it's not a bug.
 */
void MEDCouplingFieldDiscretizationP1::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const int *new2old, int newSz, DataArrayDouble *arr) const
{
}

/*!
 * This method returns a submesh of 'mesh' instance constituting cell ids contained in array defined as an interval [start;end).
* @param di is an array returned that specifies entity ids (here nodes ids) in mesh 'mesh' of entity in returned submesh.
 * Example : The first node id of returned mesh has the (*di)[0] id in 'mesh'
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationP1::buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const
{
  MEDCouplingMesh *ret=mesh->buildPartAndReduceNodes(start,end,di);
  DataArrayInt *di2=di->invertArrayO2N2N2O(ret->getNumberOfNodes());
  di->decrRef();
  di=di2;
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

MEDCouplingFieldDiscretizationPerCell::MEDCouplingFieldDiscretizationPerCell(const MEDCouplingFieldDiscretizationPerCell& other):_discr_per_cell(0)
{
  DataArrayInt *arr=other._discr_per_cell;
  if(arr)
    _discr_per_cell=arr->deepCpy();
}

void MEDCouplingFieldDiscretizationPerCell::updateTime()
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

bool MEDCouplingFieldDiscretizationPerCell::isEqual(const MEDCouplingFieldDiscretization *other, double eps) const
{
  const MEDCouplingFieldDiscretizationPerCell *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationPerCell *>(other);
  if(!otherC)
    return false;
  if(_discr_per_cell==0)
    return otherC->_discr_per_cell==0;
  if(otherC->_discr_per_cell==0)
    return false;
  return _discr_per_cell->isEqual(*otherC->_discr_per_cell);
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
    array=DataArrayInt::checkAndPreparePermutation(old2NewBg,old2NewBg+nbCells);
  //
  DataArrayInt *dpc=_discr_per_cell->renumber(array);
  _discr_per_cell->decrRef();
  _discr_per_cell=dpc;
  //
  if(check)
    delete [] (int *)array;
}

void MEDCouplingFieldDiscretizationPerCell::buildDiscrPerCellIfNecessary(const MEDCouplingMesh *m)
{
  if(!_discr_per_cell)
    {
      _discr_per_cell=DataArrayInt::New();
      int nbTuples=m->getNumberOfCells();
      _discr_per_cell->alloc(nbTuples,1);
      int *ptr=_discr_per_cell->getPointer();
      std::fill(ptr,ptr+nbTuples,-1);
    }
}

MEDCouplingFieldDiscretizationGauss::MEDCouplingFieldDiscretizationGauss()
{
}

MEDCouplingFieldDiscretizationGauss::MEDCouplingFieldDiscretizationGauss(const MEDCouplingFieldDiscretizationGauss& other):MEDCouplingFieldDiscretizationPerCell(other),_loc(other._loc)
{
}

TypeOfField MEDCouplingFieldDiscretizationGauss::getEnum() const
{
  return TYPE;
}

bool MEDCouplingFieldDiscretizationGauss::isEqual(const MEDCouplingFieldDiscretization *other, double eps) const
{
  const MEDCouplingFieldDiscretizationGauss *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(other);
  if(!otherC)
    return false;
  if(!MEDCouplingFieldDiscretizationPerCell::isEqual(other,eps))
    return false;
  if(_loc.size()!=otherC->_loc.size())
    return false;
  int sz=_loc.size();
  for(int i=0;i<sz;i++)
    if(!_loc[i].isEqual(otherC->_loc[i],eps))
      return false;
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
  int sz=_loc.size();
  for(int i=0;i<sz;i++)
    if(!_loc[i].isEqual(otherC->_loc[i],eps))
      return false;
  return true;
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationGauss::clone() const
{
  return new MEDCouplingFieldDiscretizationGauss(*this);
}

const char *MEDCouplingFieldDiscretizationGauss::getStringRepr() const
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

void MEDCouplingFieldDiscretizationGauss::renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                                                                const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  const int *array=old2NewBg;
  if(check)
    array=DataArrayInt::checkAndPreparePermutation(old2NewBg,old2NewBg+mesh->getNumberOfCells());
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
    delete [] (int*)array;
}

DataArrayDouble *MEDCouplingFieldDiscretizationGauss::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
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
  tinyInfo.push_back(_loc.size());
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
    delta=(tinyInfo.size()-3)/nbOfLoc;
  for(int i=0;i<nbOfLoc;i++)
    {
      std::vector<int> tmp(tinyInfo.begin()+3+i*delta,tinyInfo.begin()+3+(i+1)*delta);
      MEDCouplingGaussLocalization elt=MEDCouplingGaussLocalization::buildNewInstanceFromTinyInfo(dim,tmp);
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
  int nbOfDesc=_loc.size();
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

MEDCouplingMesh *MEDCouplingFieldDiscretizationGauss::buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

/*!
 * No implementation needed !
 */
void MEDCouplingFieldDiscretizationGauss::renumberValuesOnNodes(const int *, DataArrayDouble *) const
{
}

void MEDCouplingFieldDiscretizationGauss::renumberValuesOnCells(const MEDCouplingMesh *mesh, const int *old2New, DataArrayDouble *arr) const
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
  buildDiscrPerCellIfNecessary(m);
  int id=_loc.size();
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
  int id=_loc.size();
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
  return _loc.size();
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
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("No Gauss localization still set !");
  std::set<int> ret;
  int id=0;
  for(std::vector<MEDCouplingGaussLocalization>::const_iterator iter=_loc.begin();iter!=_loc.end();iter++,id++)
    if((*iter).getType()==type)
      ret.insert(id);
  if(ret.empty())
    throw INTERP_KERNEL::Exception("No gauss discretization found for the specified type !");
  if(ret.size()>1)
    throw INTERP_KERNEL::Exception("Several gauss discretizations have been found for the specified type !");
  return *ret.begin();
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
 * This method makes the assumption that _discr_per_cell is set.
 * This method reduces as much as possible number size of _loc.
 * This method is usefull when several set on same cells has been done and that some Gauss Localization are no more used.
 */
void MEDCouplingFieldDiscretizationGauss::zipGaussLocalizations()
{
  const int *start=_discr_per_cell->getConstPointer();
  int nbOfTuples=_discr_per_cell->getNumberOfTuples();
  int *tmp=new int[_loc.size()];
  std::fill(tmp,tmp+_loc.size(),-2);
  for(const int *w=start;w!=start+nbOfTuples;w++)
    if(*w>=0)
      tmp[*w]=1;
  int fid=0;
  for(int i=0;i<(int)_loc.size();i++)
    if(tmp[i]!=-2)
      tmp[i]=fid++;
  if(fid==(int)_loc.size())
    {//no zip needed
      delete [] tmp;
      return;
    }
  // zip needed
  int *start2=_discr_per_cell->getPointer();
  for(int *w2=start2;w2!=start2+nbOfTuples;w2++)
    *w2=tmp[*w2];
  std::vector<MEDCouplingGaussLocalization> tmpLoc;
  for(int i=0;i<(int)_loc.size();i++)
    if(tmp[i]!=-2)
      tmpLoc.push_back(_loc[tmp[i]]);
  delete [] tmp;
  _loc=tmpLoc;
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

const char *MEDCouplingFieldDiscretizationGaussNE::getStringRepr() const
{
  return REPR;
}

bool MEDCouplingFieldDiscretizationGaussNE::isEqual(const MEDCouplingFieldDiscretization *other, double eps) const
{
  const MEDCouplingFieldDiscretizationGaussNE *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationGaussNE *>(other);
  return otherC!=0;
}

int MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  int ret=0;
  int nbOfCells=mesh->getNumberOfCells();
  for(int i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::getCellModel(type);
      if(cm.isDynamic())
        throw INTERP_KERNEL::Exception("Not implemented yet Gauss node on elements for polygons and polyedrons !");
      ret+=cm.getNumberOfNodes();
    }
  return ret;
}

void MEDCouplingFieldDiscretizationGaussNE::renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArrayDouble *>& arrays,
                                                                  const int *old2NewBg, bool check) throw(INTERP_KERNEL::Exception)
{
  const int *array=old2NewBg;
  if(check)
    array=DataArrayInt::checkAndPreparePermutation(old2NewBg,old2NewBg+mesh->getNumberOfCells());
  int nbOfCells=mesh->getNumberOfCells();
  int nbOfTuples=getNumberOfTuples(mesh);
  int *array2=new int[nbOfTuples];//stores the final conversion array old2New to give to arrays in renumberInPlace.
  int *array3=new int[nbOfCells];//store for each cell in after renumbering the offset needed by each cell in new numbering.
  array3[0]=0;
  for(int i=1;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(std::distance(array,std::find(array,array+nbOfCells,i-1)));
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::getCellModel(type);
      array3[i]=array3[i-1]+cm.getNumberOfNodes();
    }
  int j=0;
  for(int i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::getCellModel(type);
      for(int k=0;k<(int)cm.getNumberOfNodes();k++,j++)
        array2[j]=array3[array[i]]+k;
    }
  delete [] array3;
  for(std::vector<DataArrayDouble *>::const_iterator it=arrays.begin();it!=arrays.end();it++)
    if(*it)
      (*it)->renumberInPlace(array2);
  delete [] array2;
  if(check)
    delete [] (int*)array;
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
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::getCellModel(type);
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

MEDCouplingMesh *MEDCouplingFieldDiscretizationGaussNE::buildSubMeshData(const MEDCouplingMesh *mesh, const int *start, const int *end, DataArrayInt *&di) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

/*!
 * No implementation needed !
 */
void MEDCouplingFieldDiscretizationGaussNE::renumberValuesOnNodes(const int *, DataArrayDouble *) const
{
}

void MEDCouplingFieldDiscretizationGaussNE::renumberValuesOnCells(const MEDCouplingMesh *mesh, const int *old2New, DataArrayDouble *arr) const
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

