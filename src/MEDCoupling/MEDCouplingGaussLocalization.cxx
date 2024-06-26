// Copyright (C) 2007-2024  CEA, EDF
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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingGaussLocalization.hxx"
#include "InterpKernelGaussCoords.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "CellModel.hxx"

#include <cmath>
#include <numeric>
#include <sstream>
#include <iterator>
#include <algorithm>

using namespace MEDCoupling;

MEDCouplingGaussLocalization::MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                           const std::vector<double>& gsCoo, const std::vector<double>& w)
:_type(type),_ref_coord(refCoo),_gauss_coord(gsCoo),_weight(w)
{
  // Will potentially throw (and then release memory for above objects _ref_coord etc ...)
  checkConsistencyLight();
}

MEDCouplingGaussLocalization::MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType typ)
:_type(typ)
{
  // Will potentially throw
  INTERP_KERNEL::CellModel::GetCellModel(_type);
}

void MEDCouplingGaussLocalization::setType(INTERP_KERNEL::NormalizedCellType typ)
{
  INTERP_KERNEL::CellModel::GetCellModel(typ);//throws if not found. This is a check
  _type=typ;
}

void MEDCouplingGaussLocalization::checkConsistencyLight() const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_type);
  int nbNodes=cm.getNumberOfNodes();
  int dim=cm.getDimension();
  if(!cm.isDynamic())
    {
      if(ToIdType(_ref_coord.size())!=nbNodes*dim)
        {
          std::ostringstream oss; oss << "Invalid size of refCoo : expecting to be : " << nbNodes << " (nbNodePerCell) * " << dim << " (dim) !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(_gauss_coord.size()!=dim*_weight.size())
    {
      std::ostringstream oss; oss << "Invalid gsCoo size and weight size : gsCoo.size() must be equal to _weight.size() * " << dim << " (dim) !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

int MEDCouplingGaussLocalization::getDimension() const
{
  if(_weight.empty())
    THROW_IK_EXCEPTION("getDimension : weight is empty !");
  return (int)_gauss_coord.size()/(int)_weight.size();
}

int MEDCouplingGaussLocalization::getNumberOfPtsInRefCell() const
{
  if(_gauss_coord.empty())
  {
    if(!_weight.empty())
      THROW_IK_EXCEPTION("getNumberOfPtsInRefCell : gauss_coords are empty whereas weights are not empty !");
    const INTERP_KERNEL::CellModel& cm = INTERP_KERNEL::CellModel::GetCellModel(_type);
    return ((int)_ref_coord.size()) / ((int)cm.getDimension());
  }
  int dim( getDimension() );
  return (int)_ref_coord.size()/dim;
}

std::string MEDCouplingGaussLocalization::getStringRepr() const
{
  std::ostringstream oss;
  oss << "CellType : " << INTERP_KERNEL::CellModel::GetCellModel(_type).getRepr() << std::endl;
  oss << "Ref coords : "; std::copy(_ref_coord.begin(),_ref_coord.end(),std::ostream_iterator<double>(oss,", ")); oss << std::endl;
  oss << "Localization coords : "; std::copy(_gauss_coord.begin(),_gauss_coord.end(),std::ostream_iterator<double>(oss,", ")); oss << std::endl;
  oss << "Weight : "; std::copy(_weight.begin(),_weight.end(),std::ostream_iterator<double>(oss,", ")); oss << std::endl;
  return oss.str();
}

std::size_t MEDCouplingGaussLocalization::getMemorySize() const
{
  std::size_t ret=0;
  ret+=_ref_coord.capacity()*sizeof(double);
  ret+=_gauss_coord.capacity()*sizeof(double);
  ret+=_weight.capacity()*sizeof(double);
  return ret;
}

bool MEDCouplingGaussLocalization::isEqual(const MEDCouplingGaussLocalization& other, double eps) const
{
  if(_type!=other._type)
    return false;
  if(!AreAlmostEqual(_ref_coord,other._ref_coord,eps))
    return false;
  if(!AreAlmostEqual(_gauss_coord,other._gauss_coord,eps))
    return false;
  if(!AreAlmostEqual(_weight,other._weight,eps))
    return false;
  return true;
}

double MEDCouplingGaussLocalization::getRefCoord(int ptIdInCell, int comp) const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_type);
  int nbNodes=cm.getNumberOfNodes();
  int dim=cm.getDimension();
  if(ptIdInCell<0 || ptIdInCell>=nbNodes)
    throw INTERP_KERNEL::Exception("ptIdInCell specified is invalid : must be in [0;nbNodesPerCell) !");
  if(comp<0 || comp>=dim)
    throw INTERP_KERNEL::Exception("comp specified is invalid : must be in [0:dimOfCell) !");
  return _ref_coord[ptIdInCell*dim+comp];
}

double MEDCouplingGaussLocalization::getGaussCoord(int gaussPtIdInCell, int comp) const
{
  int dim=checkCoherencyOfRequest(gaussPtIdInCell,comp);
  return _gauss_coord[gaussPtIdInCell*dim+comp];
}

double MEDCouplingGaussLocalization::getWeight(int gaussPtIdInCell) const
{
  checkCoherencyOfRequest(gaussPtIdInCell,0);
  return _weight[gaussPtIdInCell];
}

/*!
 * Completely useless method for end user. Only for CORBA MPI serialization/unserialization.
 * push at the end of tinyInfo its basic serialization info. The size of pushed data is always the same.
 * @param tinyInfo inout parameter.
 */
void MEDCouplingGaussLocalization::pushTinySerializationIntInfo(std::vector<mcIdType>& tinyInfo) const
{
  tinyInfo.push_back(ToIdType(_type));
  tinyInfo.push_back(getNumberOfPtsInRefCell());
  tinyInfo.push_back(getNumberOfGaussPt());
}

/*!
 * Completely useless method for end user. Only for CORBA MPI serialization/unserialization.
 * push at the end of tinyInfo its basic serialization info. The size of pushed data is \b NOT always the same contrary to pushTinySerializationIntInfo.
 * @param tinyInfo inout parameter.
 */
void MEDCouplingGaussLocalization::pushTinySerializationDblInfo(std::vector<double>& tinyInfo) const
{
  tinyInfo.insert(tinyInfo.end(),_ref_coord.begin(),_ref_coord.end());
  tinyInfo.insert(tinyInfo.end(),_gauss_coord.begin(),_gauss_coord.end());
  tinyInfo.insert(tinyInfo.end(),_weight.begin(),_weight.end());
}

/*!
 * This method operates the exact inverse operation than MEDCouplingGaussLocalization::pushTinySerializationDblInfo method. This is one of the last step of unserialization process.
 * This method should be called on an object resized by buildNewInstanceFromTinyInfo static method.
 * This method takes in argument a pointer 'vals' that point to the begin of double data pushed remotely by pushTinySerializationDblInfo method.
 * This method returns the pointer 'vals' with an offset of size what it has been read in this method.
 */
const double *MEDCouplingGaussLocalization::fillWithValues(const double *vals)
{
  const double *work=vals;
  std::copy(work,work+_ref_coord.size(),_ref_coord.begin());
  work+=_ref_coord.size();
  std::copy(work,work+_gauss_coord.size(),_gauss_coord.begin());
  work+=_gauss_coord.size();
  std::copy(work,work+_weight.size(),_weight.begin());
  work+=_weight.size();
  return work;
}

/*!
 * Given points in \a ptsInRefCoo in the reference coordinates of \a this (in _ref_coord attribute)
 * this method computes their coordinates in real world for each cells in \a mesh.
 * So the returned array will contain nbCells* \a ptsInRefCoo->getNumberOfTuples() tuples and the number of component
 * will be equal to the dimension of \a this.
 * 
 * This method ignores Gauss points in \a this and only those in \a ptsInRefCoo are considered here.
 */
MCAuto<DataArrayDouble> MEDCouplingGaussLocalization::localizePtsInRefCooForEachCell(const DataArrayDouble *ptsInRefCoo, const MEDCouplingUMesh *mesh) const
{
  if(!ptsInRefCoo || !mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingGaussLocalization::localizePtsInRefCooForEachCell : null pointer !");
  ptsInRefCoo->checkAllocated();
  mesh->checkConsistencyLight();
  //
  mcIdType nbCells=mesh->getNumberOfCells();
  const double *coords(mesh->getCoords()->begin());
  const mcIdType *connI(mesh->getNodalConnectivityIndex()->begin()),*conn(mesh->getNodalConnectivity()->begin());
  //
  mcIdType nbPts(ptsInRefCoo->getNumberOfTuples());
  INTERP_KERNEL::NormalizedCellType typ(getType());
  int dim(INTERP_KERNEL::CellModel::GetCellModel(typ).getDimension()),outDim(mesh->getSpaceDimension());
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbPts*nbCells,outDim);
  double *retPtr(ret->getPointer());
  if(dim!=ToIdType(ptsInRefCoo->getNumberOfComponents()))
    throw INTERP_KERNEL::Exception("MEDCouplingGaussLocalization::localizePtsInRefCooForEachCell : number of components of input coo is not equal to dim of element !");
  INTERP_KERNEL::GaussCoords calculator;
  calculator.addGaussInfo(typ,dim, ptsInRefCoo->begin(),nbPts,&_ref_coord[0],getNumberOfPtsInRefCell());
  //
  for(mcIdType i=0;i<nbCells;i++,retPtr+=nbPts*outDim)
    calculator.calculateCoords(getType(),coords,outDim,conn+connI[i]+1,retPtr);
  return ret;
}

/*!
 * This method returns an unstructured mesh that represents the reference cell.
 */
MCAuto<MEDCouplingUMesh> MEDCouplingGaussLocalization::buildRefCell() const
{
  MCAuto<DataArrayDouble> coo(DataArrayDouble::New());
  const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(getType()));
  if(getDimension()!=ToIdType(cm.getDimension()))
    throw INTERP_KERNEL::Exception("BuildRefCell : dimension mistmatch !");
  coo->alloc(cm.getNumberOfNodes(),getDimension());
  std::copy(_ref_coord.begin(),_ref_coord.end(),coo->getPointer());
  MCAuto<MEDCoupling1SGTUMesh> ret(MEDCoupling1SGTUMesh::New("",getType()));
  ret->setCoords(coo);
  MCAuto<DataArrayIdType> conn(DataArrayIdType::New());
  conn->alloc(cm.getNumberOfNodes(),1);
  conn->iota();
  ret->setNodalConnectivity(conn);
  return MCAuto<MEDCouplingUMesh>(ret->buildUnstructured());
}

/*!
 * This method returns an array containing for each Gauss Points in \a this, function values relative to the points of the
 * reference cell. Number of components of returned array is equal to the number of points of the reference cell.
 */
MCAuto<DataArrayDouble> MEDCouplingGaussLocalization::getShapeFunctionValues() const
{
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  int nbGaussPt(getNumberOfGaussPt()),nbPtsRefCell(getNumberOfPtsInRefCell()),dim(getDimension());
  ret->alloc(nbGaussPt,nbPtsRefCell);
  double *retPtr(ret->getPointer());
  for(int iGaussPt = 0 ; iGaussPt < nbGaussPt ; ++iGaussPt)
  {
    std::vector<double> curGaussPt(_gauss_coord.begin()+iGaussPt*dim,_gauss_coord.begin()+(iGaussPt+1)*dim);
    INTERP_KERNEL::GaussInfo gi(_type,curGaussPt,1,_ref_coord,nbPtsRefCell);
    gi.initLocalInfo();
    const double *funcVal( gi.getFunctionValues(0) );
    std::copy(funcVal,funcVal+nbPtsRefCell,retPtr);
    retPtr += nbPtsRefCell;
  }
  return ret;
}

MCAuto<DataArrayDouble> MEDCouplingGaussLocalization::getDerivativeOfShapeFunctionValues() const
{
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  int nbGaussPt(getNumberOfGaussPt()),nbPtsRefCell(getNumberOfPtsInRefCell()),dim(getDimension());
  ret->alloc(nbGaussPt,nbPtsRefCell*dim);
  double *retPtr(ret->getPointer());
  for(int iGaussPt = 0 ; iGaussPt < nbGaussPt ; ++iGaussPt)
  {
    std::vector<double> curGaussPt(_gauss_coord.begin()+iGaussPt*dim,_gauss_coord.begin()+(iGaussPt+1)*dim);
    INTERP_KERNEL::GaussInfo gi(_type,curGaussPt,1,_ref_coord,nbPtsRefCell);
    gi.initLocalInfo();
    const double *devOfFuncVal( gi.getDerivativeOfShapeFunctionAt(0) );
    std::copy(devOfFuncVal,devOfFuncVal+nbPtsRefCell*dim,retPtr);
    retPtr += nbPtsRefCell*dim;
  }
  return ret;
}

/*!
 * This method sets the comp_th component of ptIdInCell_th point coordinate of reference element of type this->_type.
 * @throw if not 0<=ptIdInCell<nbOfNodePerCell or if not 0<=comp<dim
 */
void MEDCouplingGaussLocalization::setRefCoord(int ptIdInCell, int comp, double newVal)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_type);
  int nbNodes=cm.getNumberOfNodes();
  int dim=cm.getDimension();
  if(ptIdInCell<0 || ptIdInCell>=nbNodes)
    throw INTERP_KERNEL::Exception("ptIdInCell specified is invalid : must be in [0;nbNodesPerCell) !");
  if(comp<0 || comp>=dim)
    throw INTERP_KERNEL::Exception("comp specified is invalid : must be in [0:dimOfCell) !");
  _ref_coord[ptIdInCell*dim+comp]=newVal;
}

void MEDCouplingGaussLocalization::setGaussCoord(int gaussPtIdInCell, int comp, double newVal)
{
  int dim=checkCoherencyOfRequest(gaussPtIdInCell,comp);
  _gauss_coord[gaussPtIdInCell*dim+comp]=newVal;
}

void MEDCouplingGaussLocalization::setWeight(int gaussPtIdInCell, double newVal)
{
  checkCoherencyOfRequest(gaussPtIdInCell,0);
  _weight[gaussPtIdInCell]=newVal;
}

void MEDCouplingGaussLocalization::setRefCoords(const std::vector<double>& refCoo)
{
  _ref_coord=refCoo;
}

void MEDCouplingGaussLocalization::setGaussCoords(const std::vector<double>& gsCoo)
{
  _gauss_coord=gsCoo;
}

void MEDCouplingGaussLocalization::setWeights(const std::vector<double>& w)
{
  _weight=w;
}

/*!
 * The format of 'tinyData' parameter is the same than pushed in method MEDCouplingGaussLocalization::pushTinySerializationIntInfo.
 */
MEDCouplingGaussLocalization MEDCouplingGaussLocalization::BuildNewInstanceFromTinyInfo(mcIdType dim, const std::vector<mcIdType>& tinyData)
{
  std::vector<double> v1(dim*tinyData[1]),v2(dim*tinyData[2]),v3(tinyData[2]);
  return MEDCouplingGaussLocalization((INTERP_KERNEL::NormalizedCellType)tinyData[0],v1,v2,v3);
}

int MEDCouplingGaussLocalization::checkCoherencyOfRequest(mcIdType gaussPtIdInCell, int comp) const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_type);
  int dim=cm.getDimension();
  mcIdType nbGsPts=getNumberOfGaussPt();
  if(gaussPtIdInCell<0 || gaussPtIdInCell>=nbGsPts)
    throw INTERP_KERNEL::Exception("gaussPtIdInCell specified is invalid : must be in [0:nbGsPts) !");
  if(comp<0 || comp>=dim)
    throw INTERP_KERNEL::Exception("comp specified is invalid : must be in [0:dimOfCell) !");
  return dim;
}

bool MEDCouplingGaussLocalization::AreAlmostEqual(const std::vector<double>& v1, const std::vector<double>& v2, double eps)
{
  std::size_t sz=v1.size();
  if(sz!=v2.size())
    return false;
  std::vector<double> tmp(sz);
  std::transform(v1.begin(),v1.end(),v2.begin(),tmp.begin(),std::minus<double>());
  std::transform(tmp.begin(),tmp.end(),tmp.begin(),[](double c){return fabs(c);});
  return *std::max_element(tmp.begin(),tmp.end())<eps;
}

MCAuto<DataArrayDouble> MEDCouplingGaussLocalization::GetDefaultReferenceCoordinatesOf(INTERP_KERNEL::NormalizedCellType type)
{
  std::vector<double> retCpp(INTERP_KERNEL::GaussInfo::GetDefaultReferenceCoordinatesOf(type));
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  auto nbDim(cm.getDimension());
  std::size_t sz(retCpp.size());
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  if( sz%std::size_t(nbDim) != 0 )
    THROW_IK_EXCEPTION("GetDefaultReferenceCoordinatesOf : unexpected size of defaut array : " << sz << " % " << nbDim << " != 0 !");
  ret->alloc(sz/size_t(nbDim),nbDim);
  std::copy(retCpp.begin(),retCpp.end(),ret->getPointer());
  return ret;
}
