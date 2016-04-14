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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingGaussLocalization.hxx"
#include "CellModel.hxx"

#include <cmath>
#include <numeric>
#include <sstream>
#include <iterator>
#include <algorithm>

MEDCoupling::MEDCouplingGaussLocalization::MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                                       const std::vector<double>& gsCoo, const std::vector<double>& w)
try:_type(type),_ref_coord(refCoo),_gauss_coord(gsCoo),_weight(w)
{
  checkConsistencyLight();
}
catch(INTERP_KERNEL::Exception& e)
{
    _type=INTERP_KERNEL::NORM_ERROR;
    _ref_coord.clear();
    _gauss_coord.clear();
    _weight.clear();
    throw e;
}

MEDCoupling::MEDCouplingGaussLocalization::MEDCouplingGaussLocalization(INTERP_KERNEL::NormalizedCellType typ)
try:_type(typ)
{
  INTERP_KERNEL::CellModel::GetCellModel(_type);
}
catch(INTERP_KERNEL::Exception& e)
{
    _type=INTERP_KERNEL::NORM_ERROR;
    throw e;
}

void MEDCoupling::MEDCouplingGaussLocalization::setType(INTERP_KERNEL::NormalizedCellType typ)
{
  INTERP_KERNEL::CellModel::GetCellModel(typ);//throws if not found. This is a check
  _type=typ;
}

void MEDCoupling::MEDCouplingGaussLocalization::checkConsistencyLight() const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_type);
  int nbNodes=cm.getNumberOfNodes();
  int dim=cm.getDimension();
  if(!cm.isDynamic())
    {
      if((int)_ref_coord.size()!=nbNodes*dim)
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

int MEDCoupling::MEDCouplingGaussLocalization::getDimension() const
{
  if(_weight.empty())
    return -1;
  return (int)_gauss_coord.size()/(int)_weight.size();
}

int MEDCoupling::MEDCouplingGaussLocalization::getNumberOfPtsInRefCell() const
{
  int dim=getDimension();
  if(dim==0)
    return -1;
  return (int)_ref_coord.size()/dim;
}

std::string MEDCoupling::MEDCouplingGaussLocalization::getStringRepr() const
{
  std::ostringstream oss;
  oss << "CellType : " << INTERP_KERNEL::CellModel::GetCellModel(_type).getRepr() << std::endl;
  oss << "Ref coords : "; std::copy(_ref_coord.begin(),_ref_coord.end(),std::ostream_iterator<double>(oss,", ")); oss << std::endl;
  oss << "Localization coords : "; std::copy(_gauss_coord.begin(),_gauss_coord.end(),std::ostream_iterator<double>(oss,", ")); oss << std::endl;
  oss << "Weight : "; std::copy(_weight.begin(),_weight.end(),std::ostream_iterator<double>(oss,", ")); oss << std::endl;
  return oss.str();
}

std::size_t MEDCoupling::MEDCouplingGaussLocalization::getMemorySize() const
{
  std::size_t ret=0;
  ret+=_ref_coord.capacity()*sizeof(double);
  ret+=_gauss_coord.capacity()*sizeof(double);
  ret+=_weight.capacity()*sizeof(double);
  return ret;
}

bool MEDCoupling::MEDCouplingGaussLocalization::isEqual(const MEDCouplingGaussLocalization& other, double eps) const
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

double MEDCoupling::MEDCouplingGaussLocalization::getRefCoord(int ptIdInCell, int comp) const
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

double MEDCoupling::MEDCouplingGaussLocalization::getGaussCoord(int gaussPtIdInCell, int comp) const
{
  int dim=checkCoherencyOfRequest(gaussPtIdInCell,comp);
  return _gauss_coord[gaussPtIdInCell*dim+comp];
}

double MEDCoupling::MEDCouplingGaussLocalization::getWeight(int gaussPtIdInCell, double newVal) const
{
  checkCoherencyOfRequest(gaussPtIdInCell,0);
  return _weight[gaussPtIdInCell];
}

/*!
 * Completely useless method for end user. Only for CORBA MPI serialization/unserialization.
 * push at the end of tinyInfo its basic serialization info. The size of pushed data is always the same.
 * @param tinyInfo inout parameter.
 */
void MEDCoupling::MEDCouplingGaussLocalization::pushTinySerializationIntInfo(std::vector<int>& tinyInfo) const
{
  tinyInfo.push_back((int)_type);
  tinyInfo.push_back(getNumberOfPtsInRefCell());
  tinyInfo.push_back(getNumberOfGaussPt());
}

/*!
 * Completely useless method for end user. Only for CORBA MPI serialization/unserialization.
 * push at the end of tinyInfo its basic serialization info. The size of pushed data is \b NOT always the same contrary to pushTinySerializationIntInfo.
 * @param tinyInfo inout parameter.
 */
void MEDCoupling::MEDCouplingGaussLocalization::pushTinySerializationDblInfo(std::vector<double>& tinyInfo) const
{
  tinyInfo.insert(tinyInfo.end(),_ref_coord.begin(),_ref_coord.end());
  tinyInfo.insert(tinyInfo.end(),_gauss_coord.begin(),_gauss_coord.end());
  tinyInfo.insert(tinyInfo.end(),_weight.begin(),_weight.end());
}

/*!
 * This method operates the exact inverse operation than MEDCoupling::MEDCouplingGaussLocalization::pushTinySerializationDblInfo method. This is one of the last step of unserialization process.
 * This method should be called on an object resized by buildNewInstanceFromTinyInfo static method.
 * This method takes in argument a pointer 'vals' that point to the begin of double data pushed remotely by pushTinySerializationDblInfo method.
 * This method returns the pointer 'vals' with an offset of size what it has been read in this method.
 */
const double *MEDCoupling::MEDCouplingGaussLocalization::fillWithValues(const double *vals)
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
 * This method sets the comp_th component of ptIdInCell_th point coordinate of reference element of type this->_type.
 * @throw if not 0<=ptIdInCell<nbOfNodePerCell or if not 0<=comp<dim
 */
void MEDCoupling::MEDCouplingGaussLocalization::setRefCoord(int ptIdInCell, int comp, double newVal)
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

void MEDCoupling::MEDCouplingGaussLocalization::setGaussCoord(int gaussPtIdInCell, int comp, double newVal)
{
  int dim=checkCoherencyOfRequest(gaussPtIdInCell,comp);
  _gauss_coord[gaussPtIdInCell*dim+comp]=newVal;
}

void MEDCoupling::MEDCouplingGaussLocalization::setWeight(int gaussPtIdInCell, double newVal)
{
  checkCoherencyOfRequest(gaussPtIdInCell,0);
  _weight[gaussPtIdInCell]=newVal;
}

void MEDCoupling::MEDCouplingGaussLocalization::setRefCoords(const std::vector<double>& refCoo)
{
  _ref_coord=refCoo;
}

void MEDCoupling::MEDCouplingGaussLocalization::setGaussCoords(const std::vector<double>& gsCoo)
{
  _gauss_coord=gsCoo;
}

void MEDCoupling::MEDCouplingGaussLocalization::setWeights(const std::vector<double>& w)
{
  _weight=w;
}

/*!
 * The format of 'tinyData' parameter is the same than pushed in method MEDCoupling::MEDCouplingGaussLocalization::pushTinySerializationIntInfo.
 */
MEDCoupling::MEDCouplingGaussLocalization MEDCoupling::MEDCouplingGaussLocalization::BuildNewInstanceFromTinyInfo(int dim, const std::vector<int>& tinyData)
{
  std::vector<double> v1(dim*tinyData[1]),v2(dim*tinyData[2]),v3(tinyData[2]);
  return MEDCoupling::MEDCouplingGaussLocalization((INTERP_KERNEL::NormalizedCellType)tinyData[0],v1,v2,v3);
}

int MEDCoupling::MEDCouplingGaussLocalization::checkCoherencyOfRequest(int gaussPtIdInCell, int comp) const
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(_type);
  int dim=cm.getDimension();
  int nbGsPts=getNumberOfGaussPt();
  if(gaussPtIdInCell<0 || gaussPtIdInCell>=nbGsPts)
    throw INTERP_KERNEL::Exception("gaussPtIdInCell specified is invalid : must be in [0:nbGsPts) !");
  if(comp<0 || comp>=dim)
    throw INTERP_KERNEL::Exception("comp specified is invalid : must be in [0:dimOfCell) !");
  return dim;
}

bool MEDCoupling::MEDCouplingGaussLocalization::AreAlmostEqual(const std::vector<double>& v1, const std::vector<double>& v2, double eps)
{
  std::size_t sz=v1.size();
  if(sz!=v2.size())
    return false;
  std::vector<double> tmp(sz);
  std::transform(v1.begin(),v1.end(),v2.begin(),tmp.begin(),std::minus<double>());
  std::transform(tmp.begin(),tmp.end(),tmp.begin(),std::ptr_fun<double,double>(fabs));
  return *std::max_element(tmp.begin(),tmp.end())<eps;
}
