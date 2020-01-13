// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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

#include "MEDCouplingMemArray.txx"

#include "BBTree.txx"
#include "GenMathFormulae.hxx"
#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelExprParser.hxx"

#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"

#include <set>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <functional>

typedef double (*MYFUNCPTR)(double);

using namespace MEDCoupling;

template class MEDCOUPLING_EXPORT MEDCoupling::MemArray<mcIdType>;
template class MEDCOUPLING_EXPORT MEDCoupling::MemArray<double>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayTemplate<mcIdType>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayTemplate<double>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayTemplateClassic<Int32>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayTemplateClassic<Int64>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayTemplateClassic<double>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayTemplateFP<double>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayIterator<double>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayIterator<mcIdType>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayDiscrete<Int32>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayDiscreteSigned<Int32>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayDiscrete<Int64>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayDiscreteSigned<Int64>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayTuple<mcIdType>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayTuple<double>;
template class MEDCOUPLING_EXPORT MEDCoupling::DataArrayTuple<float>;

template<mcIdType SPACEDIM>
void DataArrayDouble::findCommonTuplesAlg(const double *bbox, mcIdType nbNodes, mcIdType limitNodeId, double prec, DataArrayIdType *c, DataArrayIdType *cI) const
{
  const double *coordsPtr=getConstPointer();
  BBTreePts<SPACEDIM,mcIdType> myTree(bbox,0,0,nbNodes,prec);
  std::vector<bool> isDone(nbNodes);
  for(mcIdType i=0;i<nbNodes;i++)
    {
      if(!isDone[i])
        {
          std::vector<mcIdType> intersectingElems;
          myTree.getElementsAroundPoint(coordsPtr+i*SPACEDIM,intersectingElems);
          if(intersectingElems.size()>1)
            {
              std::vector<mcIdType> commonNodes;
              for(std::vector<mcIdType>::const_iterator it=intersectingElems.begin();it!=intersectingElems.end();it++)
                if(*it!=i)
                  if(*it>=limitNodeId)
                    {
                      commonNodes.push_back(*it);
                      isDone[*it]=true;
                    }
              if(!commonNodes.empty())
                {
                  cI->pushBackSilent(cI->back()+ToIdType(commonNodes.size())+1);
                  c->pushBackSilent(i);
                  c->insertAtTheEnd(commonNodes.begin(),commonNodes.end());
                }
            }
        }
    }
}

template<mcIdType SPACEDIM>
void DataArrayDouble::FindTupleIdsNearTuplesAlg(const BBTreePts<SPACEDIM,mcIdType>& myTree, const double *pos, mcIdType nbOfTuples, double eps,
                                                DataArrayIdType *c, DataArrayIdType *cI)
{
  for(mcIdType i=0;i<nbOfTuples;i++)
    {
      std::vector<mcIdType> intersectingElems;
      myTree.getElementsAroundPoint(pos+i*SPACEDIM,intersectingElems);
      std::vector<mcIdType> commonNodes;
      for(std::vector<mcIdType>::const_iterator it=intersectingElems.begin();it!=intersectingElems.end();it++)
        commonNodes.push_back(*it);
      cI->pushBackSilent(cI->back()+ToIdType(commonNodes.size()));
      c->insertAtTheEnd(commonNodes.begin(),commonNodes.end());
    }
}

template<mcIdType SPACEDIM>
void DataArrayDouble::FindClosestTupleIdAlg(const BBTreePts<SPACEDIM,mcIdType>& myTree, double dist, const double *pos, mcIdType nbOfTuples, const double *thisPt, mcIdType thisNbOfTuples, mcIdType *res)
{
  double distOpt(dist);
  const double *p(pos);
  mcIdType *r(res);
  for(mcIdType i=0;i<nbOfTuples;i++,p+=SPACEDIM,r++)
    {
      while(true)
        {
          mcIdType elem=-1;
          double ret=myTree.getElementsAroundPoint2(p,distOpt,elem);
          if(ret!=std::numeric_limits<double>::max())
            {
              distOpt=std::max(ret,1e-4);
              *r=elem;
              break;
            }
          else
            { distOpt=2*distOpt; continue; }
        }
    }
}

mcIdType DataArray::EffectiveCircPerm(mcIdType nbOfShift, mcIdType nbOfTuples)
{
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArray::EffectiveCircPerm : number of tuples is expected to be > 0 !");
  if(nbOfShift>=0)
    {
      return nbOfShift%nbOfTuples;
    }
  else
    {
      mcIdType tmp(-nbOfShift);
      tmp=tmp%nbOfTuples;
      return nbOfTuples-tmp;
    }
}

std::size_t DataArray::getHeapMemorySizeWithoutChildren() const
{
  std::size_t sz1=_name.capacity();
  std::size_t sz2=_info_on_compo.capacity();
  std::size_t sz3=0;
  for(std::vector<std::string>::const_iterator it=_info_on_compo.begin();it!=_info_on_compo.end();it++)
    sz3+=(*it).capacity();
  return sz1+sz2+sz3;
}

std::vector<const BigMemoryObject *> DataArray::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

/*!
 * Sets the attribute \a _name of \a this array.
 * See \ref MEDCouplingArrayBasicsName "DataArrays infos" for more information.
 *  \param [in] name - new array name
 */
void DataArray::setName(const std::string& name)
{
  _name=name;
}

/*!
 * Copies textual data from an \a other DataArray. The copied data are
 * - the name attribute,
 * - the information of components.
 *
 * For more information on these data see \ref MEDCouplingArrayBasicsName "DataArrays infos".
 *
 *  \param [in] other - another instance of DataArray to copy the textual data from.
 *  \throw If number of components of \a this array differs from that of the \a other.
 */
void DataArray::copyStringInfoFrom(const DataArray& other)
{
  if(_info_on_compo.size()!=other._info_on_compo.size())
    throw INTERP_KERNEL::Exception("Size of arrays mismatches on copyStringInfoFrom !");
  _name=other._name;
  _info_on_compo=other._info_on_compo;
}

void DataArray::copyPartOfStringInfoFrom(const DataArray& other, const std::vector<std::size_t>& compoIds)
{
  std::size_t nbOfCompoOth=other.getNumberOfComponents();
  std::size_t newNbOfCompo=compoIds.size();
  for(std::size_t i=0;i<newNbOfCompo;i++)
    if(compoIds[i]>=nbOfCompoOth || compoIds[i]<0)
      {
        std::ostringstream oss; oss << "Specified component id is out of range (" << compoIds[i] << ") compared with nb of actual components (" << nbOfCompoOth << ")";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  for(std::size_t i=0;i<newNbOfCompo;i++)
    setInfoOnComponent(i,other.getInfoOnComponent(compoIds[i]));
}

void DataArray::copyPartOfStringInfoFrom2(const std::vector<std::size_t>& compoIds, const DataArray& other)
{
  if(compoIds.size()!=other.getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Given compoIds has not the same size as number of components of given array !");
  std::size_t partOfCompoToSet=compoIds.size();
  std::size_t nbOfCompo=getNumberOfComponents();
  for(std::size_t i=0;i<partOfCompoToSet;i++)
    if(compoIds[i]>=nbOfCompo || compoIds[i]<0)
      {
        std::ostringstream oss; oss << "Specified component id is out of range (" << compoIds[i] << ") compared with nb of actual components (" << nbOfCompo << ")";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  for(std::size_t i=0;i<partOfCompoToSet;i++)
    setInfoOnComponent(compoIds[i],other.getInfoOnComponent(i));
}

bool DataArray::areInfoEqualsIfNotWhy(const DataArray& other, std::string& reason) const
{
  std::ostringstream oss;
  if(_name!=other._name)
    {
      oss << "Names DataArray mismatch : this name=\"" << _name << " other name=\"" << other._name << "\" !";
      reason=oss.str();
      return false;
    }
  if(_info_on_compo!=other._info_on_compo)
    {
      oss << "Components DataArray mismatch : \nThis components=";
      for(std::vector<std::string>::const_iterator it=_info_on_compo.begin();it!=_info_on_compo.end();it++)
        oss << "\"" << *it << "\",";
      oss << "\nOther components=";
      for(std::vector<std::string>::const_iterator it=other._info_on_compo.begin();it!=other._info_on_compo.end();it++)
        oss << "\"" << *it << "\",";
      reason=oss.str();
      return false;
    }
  return true;
}

/*!
 * Compares textual information of \a this DataArray with that of an \a other one.
 * The compared data are
 * - the name attribute,
 * - the information of components.
 *
 * For more information on these data see \ref MEDCouplingArrayBasicsName "DataArrays infos".
 *  \param [in] other - another instance of DataArray to compare the textual data of.
 *  \return bool - \a true if the textual information is same, \a false else.
 */
bool DataArray::areInfoEquals(const DataArray& other) const
{
  std::string tmp;
  return areInfoEqualsIfNotWhy(other,tmp);
}

void DataArray::reprWithoutNameStream(std::ostream& stream) const
{
  stream << "Number of components : "<< getNumberOfComponents() << "\n";
  stream << "Info of these components : ";
  for(std::vector<std::string>::const_iterator iter=_info_on_compo.begin();iter!=_info_on_compo.end();iter++)
    stream << "\"" << *iter << "\"   ";
  stream << "\n";
}

std::string DataArray::cppRepr(const std::string& varName) const
{
  std::ostringstream ret;
  reprCppStream(varName,ret);
  return ret.str();
}

/*!
 * Sets information on all components. To know more on format of this information
 * see \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] info - a vector of strings.
 *  \throw If size of \a info differs from the number of components of \a this.
 */
void DataArray::setInfoOnComponents(const std::vector<std::string>& info)
{
  if(getNumberOfComponents()!=info.size())
    {
      std::ostringstream oss; oss << "DataArray::setInfoOnComponents : input is of size " << info.size() << " whereas number of components is equal to " << getNumberOfComponents() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  _info_on_compo=info;
}

/*!
 * This method is only a dispatcher towards DataArrayDouble::setPartOfValues3, DataArrayInt::setPartOfValues3, DataArrayChar::setPartOfValues3 depending on the true
 * type of \a this and \a aBase.
 *
 * \throw If \a aBase and \a this do not have the same type.
 *
 * \sa DataArrayDouble::setPartOfValues3, DataArrayInt::setPartOfValues3, DataArrayChar::setPartOfValues3.
 */
void DataArray::setPartOfValuesBase3(const DataArray *aBase, const mcIdType *bgTuples, const mcIdType *endTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare)
{
  if(!aBase)
    throw INTERP_KERNEL::Exception("DataArray::setPartOfValuesBase3 : input aBase object is NULL !");
  DataArrayDouble *this1(dynamic_cast<DataArrayDouble *>(this));
  DataArrayIdType *this2(dynamic_cast<DataArrayIdType *>(this));
  DataArrayChar *this3(dynamic_cast<DataArrayChar *>(this));
  const DataArrayDouble *a1(dynamic_cast<const DataArrayDouble *>(aBase));
  const DataArrayIdType *a2(dynamic_cast<const DataArrayIdType *>(aBase));
  const DataArrayChar *a3(dynamic_cast<const DataArrayChar *>(aBase));
  if(this1 && a1)
    {
      this1->setPartOfValues3(a1,bgTuples,endTuples,bgComp,endComp,stepComp,strictCompoCompare);
      return ;
    }
  if(this2 && a2)
    {
      this2->setPartOfValues3(a2,bgTuples,endTuples,bgComp,endComp,stepComp,strictCompoCompare);
      return ;
    }
  if(this3 && a3)
    {
      this3->setPartOfValues3(a3,bgTuples,endTuples,bgComp,endComp,stepComp,strictCompoCompare);
      return ;
    }
  throw INTERP_KERNEL::Exception("DataArray::setPartOfValuesBase3 : input aBase object and this do not have the same type !");
}

std::vector<std::string> DataArray::getVarsOnComponent() const
{
  std::size_t nbOfCompo=_info_on_compo.size();
  std::vector<std::string> ret(nbOfCompo);
  for(std::size_t i=0;i<nbOfCompo;i++)
    ret[i]=getVarOnComponent(i);
  return ret;
}

std::vector<std::string> DataArray::getUnitsOnComponent() const
{
  std::size_t nbOfCompo=_info_on_compo.size();
  std::vector<std::string> ret(nbOfCompo);
  for(std::size_t i=0;i<nbOfCompo;i++)
    ret[i]=getUnitOnComponent(i);
  return ret;
}

/*!
 * Returns information on a component specified by an index.
 * To know more on format of this information
 * see \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] i - the index (zero based) of the component of interest.
 *  \return std::string - a string containing the information on \a i-th component.
 *  \throw If \a i is not a valid component index.
 */
std::string DataArray::getInfoOnComponent(std::size_t i) const
{
  if(i<_info_on_compo.size())
    return _info_on_compo[i];
  else
    {
      std::ostringstream oss; oss << "DataArray::getInfoOnComponent : Specified component id is out of range (" << i << ") compared with nb of actual components (" << _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * Returns the var part of the full information of the \a i-th component.
 * For example, if \c getInfoOnComponent(0) returns "SIGXY [N/m^2]", then
 * \c getVarOnComponent(0) returns "SIGXY".
 * If a unit part of information is not detected by presence of
 * two square brackets, then the full information is returned.
 * To read more about the component information format, see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] i - the index (zero based) of the component of interest.
 *  \return std::string - a string containing the var information, or the full info.
 *  \throw If \a i is not a valid component index.
 */
std::string DataArray::getVarOnComponent(std::size_t i) const
{
  if(i<_info_on_compo.size())
    {
      return GetVarNameFromInfo(_info_on_compo[i]);
    }
  else
    {
      std::ostringstream oss; oss << "DataArray::getVarOnComponent : Specified component id is out of range  (" << i << ") compared with nb of actual components (" << _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * Returns the unit part of the full information of the \a i-th component.
 * For example, if \c getInfoOnComponent(0) returns "SIGXY [ N/m^2]", then
 * \c getUnitOnComponent(0) returns " N/m^2".
 * If a unit part of information is not detected by presence of
 * two square brackets, then an empty string is returned.
 * To read more about the component information format, see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] i - the index (zero based) of the component of interest.
 *  \return std::string - a string containing the unit information, if any, or "".
 *  \throw If \a i is not a valid component index.
 */
std::string DataArray::getUnitOnComponent(std::size_t i) const
{
  if(i<_info_on_compo.size())
    {
      return GetUnitFromInfo(_info_on_compo[i]);
    }
  else
    {
      std::ostringstream oss; oss << "DataArray::getUnitOnComponent : Specified component id is out of range  (" << i << ") compared with nb of actual components (" << _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * Returns the var part of the full component information.
 * For example, if \a info == "SIGXY [N/m^2]", then this method returns "SIGXY".
 * If a unit part of information is not detected by presence of
 * two square brackets, then the whole \a info is returned.
 * To read more about the component information format, see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] info - the full component information.
 *  \return std::string - a string containing only var information, or the \a info.
 */
std::string DataArray::GetVarNameFromInfo(const std::string& info)
{
  std::size_t p1=info.find_last_of('[');
  std::size_t p2=info.find_last_of(']');
  if(p1==std::string::npos || p2==std::string::npos)
    return info;
  if(p1>p2)
    return info;
  if(p1==0)
    return std::string();
  std::size_t p3=info.find_last_not_of(' ',p1-1);
  return info.substr(0,p3+1);
}

/*!
 * Returns the unit part of the full component information.
 * For example, if \a info == "SIGXY [ N/m^2]", then this method returns " N/m^2".
 * If a unit part of information is not detected by presence of
 * two square brackets, then an empty string is returned.
 * To read more about the component information format, see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] info - the full component information.
 *  \return std::string - a string containing only unit information, if any, or "".
 */
std::string DataArray::GetUnitFromInfo(const std::string& info)
{
  std::size_t p1=info.find_last_of('[');
  std::size_t p2=info.find_last_of(']');
  if(p1==std::string::npos || p2==std::string::npos)
    return std::string();
  if(p1>p2)
    return std::string();
  return info.substr(p1+1,p2-p1-1);
}

/*!
 * This method put in info format the result of the merge of \a var and \a unit.
 * The standard format for that is "var [unit]".
 * Inversely you can retrieve the var part or the unit part of info string using resp. GetVarNameFromInfo and GetUnitFromInfo.
 */
std::string DataArray::BuildInfoFromVarAndUnit(const std::string& var, const std::string& unit)
{
  std::ostringstream oss;
  oss << var << " [" << unit << "]";
  return oss.str();
}

std::string DataArray::GetAxisTypeRepr(MEDCouplingAxisType at)
{
  switch(at)
    {
    case AX_CART:
      return std::string("AX_CART");
    case AX_CYL:
      return std::string("AX_CYL");
    case AX_SPHER:
      return std::string("AX_SPHER");
    default:
      throw INTERP_KERNEL::Exception("DataArray::GetAxisTypeRepr : unrecognized axis type enum !");
    }
}

/*!
 * Returns a new DataArray by concatenating all given arrays, so that (1) the number
 * of tuples in the result array is a sum of the number of tuples of given arrays and (2)
 * the number of component in the result array is same as that of each of given arrays.
 * Info on components is copied from the first of the given arrays. Number of components
 * in the given arrays must be  the same.
 *  \param [in] arrs - a sequence of arrays to include in the result array. All arrays must have the same type.
 *  \return DataArray * - the new instance of DataArray (that can be either DataArrayInt, DataArrayDouble, DataArrayChar).
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If all arrays within \a arrs are NULL.
 *  \throw If all not null arrays in \a arrs have not the same type.
 *  \throw If getNumberOfComponents() of arrays within \a arrs.
 */
DataArray *DataArray::Aggregate(const std::vector<const DataArray *>& arrs)
{
  std::vector<const DataArray *> arr2;
  for(std::vector<const DataArray *>::const_iterator it=arrs.begin();it!=arrs.end();it++)
    if(*it)
      arr2.push_back(*it);
  if(arr2.empty())
    throw INTERP_KERNEL::Exception("DataArray::Aggregate : only null instance in input vector !");
  std::vector<const DataArrayDouble *> arrd;
  std::vector<const DataArrayIdType *> arri;
  std::vector<const DataArrayChar *> arrc;
  for(std::vector<const DataArray *>::const_iterator it=arr2.begin();it!=arr2.end();it++)
    {
      const DataArrayDouble *a=dynamic_cast<const DataArrayDouble *>(*it);
      if(a)
        { arrd.push_back(a); continue; }
      const DataArrayIdType *b=dynamic_cast<const DataArrayIdType *>(*it);
      if(b)
        { arri.push_back(b); continue; }
      const DataArrayChar *c=dynamic_cast<const DataArrayChar *>(*it);
      if(c)
        { arrc.push_back(c); continue; }
      throw INTERP_KERNEL::Exception("DataArray::Aggregate : presence of not null instance in inuput that is not in [DataArrayDouble, DataArrayInt, DataArrayChar] !");
    }
  if(arr2.size()==arrd.size())
    return DataArrayDouble::Aggregate(arrd);
  if(arr2.size()==arri.size())
    return DataArrayIdType::Aggregate(arri);
  if(arr2.size()==arrc.size())
    return DataArrayChar::Aggregate(arrc);
  throw INTERP_KERNEL::Exception("DataArray::Aggregate : all input arrays must have the same type !");
}

/*!
 * Sets information on a component specified by an index.
 * To know more on format of this information
 * see \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \warning Don't pass NULL as \a info!
 *  \param [in] i - the index (zero based) of the component of interest.
 *  \param [in] info - the string containing the information.
 *  \throw If \a i is not a valid component index.
 */
void DataArray::setInfoOnComponent(std::size_t i, const std::string& info)
{
  if(i<_info_on_compo.size())
    _info_on_compo[i]=info;
  else
    {
      std::ostringstream oss; oss << "DataArray::setInfoOnComponent : Specified component id is out of range  (" << i << ") compared with nb of actual components (" << _info_on_compo.size();
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * Sets information on all components. This method can change number of components
 * at certain conditions; if the conditions are not respected, an exception is thrown.
 * The number of components can be changed in \a this only if \a this is not allocated.
 * The condition of number of components must not be changed.
 *
 * To know more on format of the component information see
 * \ref MEDCouplingArrayBasicsCompoName "DataArrays infos".
 *  \param [in] info - a vector of component infos.
 *  \throw If \a this->getNumberOfComponents() != \a info.size() && \a this->isAllocated()
 */
void DataArray::setInfoAndChangeNbOfCompo(const std::vector<std::string>& info)
{
  if(getNumberOfComponents()!=info.size())
    {
      if(!isAllocated())
        _info_on_compo=info;
      else
        {
          std::ostringstream oss; oss << "DataArray::setInfoAndChangeNbOfCompo : input is of size " << info.size() << " whereas number of components is equal to " << getNumberOfComponents() << "  and this is already allocated !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  else
    _info_on_compo=info;
}

void DataArray::checkNbOfTuples(mcIdType nbOfTuples, const std::string& msg) const
{
  if(getNumberOfTuples()!=nbOfTuples)
    {
      std::ostringstream oss; oss << msg << " : mismatch number of tuples : expected " <<  nbOfTuples << " having " << getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::checkNbOfComps(std::size_t nbOfCompo, const std::string& msg) const
{
  if (getNumberOfComponents()!=nbOfCompo)
    {
      std::ostringstream oss; oss << msg << " : mismatch number of components : expected " << nbOfCompo << " having " << getNumberOfComponents() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::checkNbOfElems(mcIdType nbOfElems, const std::string& msg) const
{
  if(getNbOfElems()!=nbOfElems)
    {
      std::ostringstream oss; oss << msg << " : mismatch number of elems : Expected " << nbOfElems << " having " << getNbOfElems() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::checkNbOfTuplesAndComp(const DataArray& other, const std::string& msg) const
{
  if(getNumberOfTuples()!=other.getNumberOfTuples())
    {
      std::ostringstream oss; oss << msg << " : mismatch number of tuples : expected " <<  other.getNumberOfTuples() << " having " << getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(getNumberOfComponents()!=other.getNumberOfComponents())
    {
      std::ostringstream oss; oss << msg << " : mismatch number of components : expected " << other.getNumberOfComponents() << " having " << getNumberOfComponents() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::checkNbOfTuplesAndComp(mcIdType nbOfTuples, std::size_t nbOfCompo, const std::string& msg) const
{
  checkNbOfTuples(nbOfTuples,msg);
  checkNbOfComps(nbOfCompo,msg);
}

/*!
 * Simply this method checks that \b value is in [0,\b ref).
 */
void DataArray::CheckValueInRange(mcIdType ref, mcIdType value, const std::string& msg)
{
  if(value<0 || value>=ref)
    {
      std::ostringstream oss; oss << "DataArray::CheckValueInRange : " << msg  << " ! Expected in range [0," << ref << "[ having " << value << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * This method checks that [\b start, \b end) is compliant with ref length \b value.
 * typically start in [0,\b value) and end in [0,\b value). If value==start and start==end, it is supported.
 */
void DataArray::CheckValueInRangeEx(mcIdType value, mcIdType start, mcIdType end, const std::string& msg)
{
  if(start<0 || start>=value)
    {
      if(value!=start || end!=start)
        {
          std::ostringstream oss; oss << "DataArray::CheckValueInRangeEx : " << msg  << " ! Expected start " << start << " of input range, in [0," << value << "[ !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(end<0 || end>value)
    {
      std::ostringstream oss; oss << "DataArray::CheckValueInRangeEx : " << msg  << " ! Expected end " << end << " of input range, in [0," << value << "] !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

void DataArray::CheckClosingParInRange(mcIdType ref, mcIdType value, const std::string& msg)
{
  if(value<0 || value>ref)
    {
      std::ostringstream oss; oss << "DataArray::CheckClosingParInRange : " << msg  << " ! Expected input range in [0," << ref << "] having closing open parenthesis " << value << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * This method is useful to slice work among a pool of threads or processes. \a begin, \a end \a step is the input whole slice of work to perform,
 * typically it is a whole slice of tuples of DataArray or cells, nodes of a mesh...
 *
 * The input \a sliceId should be an id in [0, \a nbOfSlices) that specifies the slice of work.
 *
 * \param [in] start - the start of the input slice of the whole work to perform split into slices.
 * \param [in] stop - the stop of the input slice of the whole work to perform split into slices.
 * \param [in] step - the step (that can be <0) of the input slice of the whole work to perform split into slices.
 * \param [in] sliceId - the slice id considered
 * \param [in] nbOfSlices - the number of slices (typically the number of cores on which the work is expected to be sliced)
 * \param [out] startSlice - the start of the slice considered
 * \param [out] stopSlice - the stop of the slice consided
 *
 * \throw If \a step == 0
 * \throw If \a nbOfSlices not > 0
 * \throw If \a sliceId not in [0,nbOfSlices)
 */
void DataArray::GetSlice(mcIdType start, mcIdType stop, mcIdType step, mcIdType sliceId, mcIdType nbOfSlices, mcIdType& startSlice, mcIdType& stopSlice)
{
  DataArrayTools<mcIdType>::GetSlice(start, stop, step, sliceId, nbOfSlices, startSlice, stopSlice);
}

mcIdType DataArray::GetNumberOfItemGivenBES(mcIdType begin, mcIdType end, mcIdType step, const std::string& msg)
{
  return DataArrayTools<mcIdType>::GetNumberOfItemGivenBES(begin, end, step, msg);
}

mcIdType DataArray::GetNumberOfItemGivenBESRelative(mcIdType begin, mcIdType end, mcIdType step, const std::string& msg)
{
  return DataArrayTools<mcIdType>::GetNumberOfItemGivenBESRelative(begin, end, step, msg);
}

mcIdType DataArray::GetPosOfItemGivenBESRelativeNoThrow(mcIdType value, mcIdType begin, mcIdType end, mcIdType step)
{
  return DataArrayTools<mcIdType>::GetPosOfItemGivenBESRelativeNoThrow(value, begin, end, step);
}

/*!
 * Returns a new instance of DataArrayDouble. The caller is to delete this array
 * using decrRef() as it is no more needed.
 */
DataArrayDouble *DataArrayDouble::New()
{
  return new DataArrayDouble;
}

/*!
 * Returns the only one value in \a this, if and only if number of elements
 * (nb of tuples * nb of components) is equal to 1, and that \a this is allocated.
 *  \return double - the sole value stored in \a this array.
 *  \throw If at least one of conditions stated above is not fulfilled.
 */
double DataArrayDouble::doubleValue() const
{
  if(isAllocated())
    {
      if(getNbOfElems()==1)
        {
          return *getConstPointer();
        }
      else
        throw INTERP_KERNEL::Exception("DataArrayDouble::doubleValue : DataArrayDouble instance is allocated but number of elements is not equal to 1 !");
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayDouble::doubleValue : DataArrayDouble instance is not allocated !");
}

/*!
 * Returns a full copy of \a this. For more info on copying data arrays see
 * \ref MEDCouplingArrayBasicsCopyDeep.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble. The caller is to
 *          delete this array using decrRef() as it is no more needed.
 */
DataArrayDouble *DataArrayDouble::deepCopy() const
{
  return new DataArrayDouble(*this);
}

/*!
 * Checks that \a this array is consistently **increasing** or **decreasing** in value,
 * with at least absolute difference value of |\a eps| at each step.
 * If not an exception is thrown.
 *  \param [in] increasing - if \a true, the array values should be increasing.
 *  \param [in] eps - minimal absolute difference between the neighbor values at which
 *                    the values are considered different.
 *  \throw If sequence of values is not strictly monotonic in agreement with \a
 *         increasing arg.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this is not allocated.
 */
void DataArrayDouble::checkMonotonic(bool increasing, double eps) const
{
  if(!isMonotonic(increasing,eps))
    {
      if (increasing)
        throw INTERP_KERNEL::Exception("DataArrayDouble::checkMonotonic : 'this' is not INCREASING monotonic !");
      else
        throw INTERP_KERNEL::Exception("DataArrayDouble::checkMonotonic : 'this' is not DECREASING monotonic !");
    }
}

/*!
 * Checks that \a this array is consistently **increasing** or **decreasing** in value,
 * with at least absolute difference value of |\a eps| at each step.
 *  \param [in] increasing - if \a true, array values should be increasing.
 *  \param [in] eps - minimal absolute difference between the neighbor values at which
 *                    the values are considered different.
 *  \return bool - \a true if values change in accordance with \a increasing arg.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this is not allocated.
 */
bool DataArrayDouble::isMonotonic(bool increasing, double eps) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::isMonotonic : only supported with 'this' array with ONE component !");
  mcIdType nbOfElements(getNumberOfTuples());
  const double *ptr=getConstPointer();
  if(nbOfElements==0)
    return true;
  double ref=ptr[0];
  double absEps=fabs(eps);
  if(increasing)
    {
      for(mcIdType i=1;i<nbOfElements;i++)
        {
          if(ptr[i]<(ref+absEps))
            return false;
          ref=ptr[i];
        }
      return true;
    }
  else
    {
      for(mcIdType i=1;i<nbOfElements;i++)
        {
          if(ptr[i]>(ref-absEps))
            return false;
          ref=ptr[i];
        }
      return true;
    }
}

void DataArrayDouble::writeVTK(std::ostream& ofs, mcIdType indent, const std::string& nameInFile, DataArrayByte *byteArr) const
{
  static const char SPACE[4]={' ',' ',' ',' '};
  checkAllocated();
  std::string idt(indent,' ');
  ofs.precision(17);
  ofs << idt << "<DataArray type=\"Float32\" Name=\"" << nameInFile << "\" NumberOfComponents=\"" << getNumberOfComponents() << "\"";
  //
  bool areAllEmpty(true);
  for(std::vector<std::string>::const_iterator it=_info_on_compo.begin();it!=_info_on_compo.end();it++)
    if(!(*it).empty())
      areAllEmpty=false;
  if(!areAllEmpty)
    for(std::size_t i=0;i<_info_on_compo.size();i++)
      ofs << " ComponentName" << i << "=\"" << _info_on_compo[i] << "\"";
  //
  if(byteArr)
    {
      ofs << " format=\"appended\" offset=\"" << byteArr->getNumberOfTuples() << "\">";
      INTERP_KERNEL::AutoPtr<float> tmp(new float[getNbOfElems()]);
      float *pt(tmp);
      // to make Visual C++ happy : instead of std::copy(begin(),end(),(float *)tmp);
      for(const double *src=begin();src!=end();src++,pt++)
        *pt=float(*src);
      const char *data(reinterpret_cast<const char *>((float *)tmp));
      std::size_t sz(getNbOfElems()*sizeof(float));
      byteArr->insertAtTheEnd(data,data+sz);
      byteArr->insertAtTheEnd(SPACE,SPACE+4);
    }
  else
    {
      ofs << " RangeMin=\"" << getMinValueInArray() << "\" RangeMax=\"" << getMaxValueInArray() << "\" format=\"ascii\">\n" << idt;
      std::copy(begin(),end(),std::ostream_iterator<double>(ofs," "));
    }
  ofs << std::endl << idt << "</DataArray>\n";
}

void DataArrayDouble::reprCppStream(const std::string& varName, std::ostream& stream) const
{
  mcIdType nbTuples=getNumberOfTuples();
  std::size_t nbComp=getNumberOfComponents();
  const double *data(getConstPointer());
  stream.precision(17);
  stream << "DataArrayDouble *" << varName << "=DataArrayDouble::New();" << std::endl;
  if(nbTuples*nbComp>=1)
    {
      stream << "const double " << varName << "Data[" << nbTuples*nbComp << "]={";
      std::copy(data,data+nbTuples*nbComp-1,std::ostream_iterator<double>(stream,","));
      stream << data[nbTuples*nbComp-1] << "};" << std::endl;
      stream << varName << "->useArray(" << varName << "Data,false,CPP_DEALLOC," << nbTuples << "," << nbComp << ");" << std::endl;
    }
  else
    stream << varName << "->alloc(" << nbTuples << "," << nbComp << ");" << std::endl;
  stream << varName << "->setName(\"" << getName() << "\");" << std::endl;
}

/*!
 * Method that gives a quick overvien of \a this for python.
 */
void DataArrayDouble::reprQuickOverview(std::ostream& stream) const
{
  static const std::size_t MAX_NB_OF_BYTE_IN_REPR=300;
  stream << "DataArrayDouble C++ instance at " << this << ". ";
  if(isAllocated())
    {
      std::size_t nbOfCompo(_info_on_compo.size());
      if(nbOfCompo>=1)
        {
          mcIdType nbOfTuples(getNumberOfTuples());
          stream << "Number of tuples : " << nbOfTuples << ". Number of components : " << nbOfCompo << "." << std::endl;
          reprQuickOverviewData(stream,MAX_NB_OF_BYTE_IN_REPR);
        }
      else
        stream << "Number of components : 0.";
    }
  else
    stream << "*** No data allocated ****";
}

void DataArrayDouble::reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const
{
  const double *data=begin();
  mcIdType nbOfTuples(getNumberOfTuples());
  std::size_t nbOfCompo(_info_on_compo.size());
  std::ostringstream oss2; oss2 << "[";
  oss2.precision(17);
  std::string oss2Str(oss2.str());
  bool isFinished=true;
  for(mcIdType i=0;i<nbOfTuples && isFinished;i++)
    {
      if(nbOfCompo>1)
        {
          oss2 << "(";
          for(std::size_t j=0;j<nbOfCompo;j++,data++)
            {
              oss2 << *data;
              if(j!=nbOfCompo-1) oss2 << ", ";
            }
          oss2 << ")";
        }
      else
        oss2 << *data++;
      if(i!=nbOfTuples-1) oss2 << ", ";
      std::string oss3Str(oss2.str());
      if(oss3Str.length()<maxNbOfByteInRepr)
        oss2Str=oss3Str;
      else
        isFinished=false;
    }
  stream << oss2Str;
  if(!isFinished)
    stream << "... ";
  stream << "]";
}

/*!
 * Equivalent to DataArrayDouble::isEqual except that if false the reason of
 * mismatch is given.
 *
 * \param [in] other the instance to be compared with \a this
 * \param [in] prec the precision to compare numeric data of the arrays.
 * \param [out] reason In case of inequality returns the reason.
 * \sa DataArrayDouble::isEqual
 */
bool DataArrayDouble::isEqualIfNotWhy(const DataArrayDouble& other, double prec, std::string& reason) const
{
  if(!areInfoEqualsIfNotWhy(other,reason))
    return false;
  return _mem.isEqual(other._mem,prec,reason);
}

/*!
 * Checks if \a this and another DataArrayDouble are fully equal. For more info see
 * \ref MEDCouplingArrayBasicsCompare.
 *  \param [in] other - an instance of DataArrayDouble to compare with \a this one.
 *  \param [in] prec - precision value to compare numeric data of the arrays.
 *  \return bool - \a true if the two arrays are equal, \a false else.
 */
bool DataArrayDouble::isEqual(const DataArrayDouble& other, double prec) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,prec,tmp);
}

/*!
 * Checks if values of \a this and another DataArrayDouble are equal. For more info see
 * \ref MEDCouplingArrayBasicsCompare.
 *  \param [in] other - an instance of DataArrayDouble to compare with \a this one.
 *  \param [in] prec - precision value to compare numeric data of the arrays.
 *  \return bool - \a true if the values of two arrays are equal, \a false else.
 */
bool DataArrayDouble::isEqualWithoutConsideringStr(const DataArrayDouble& other, double prec) const
{
  std::string tmp;
  return _mem.isEqual(other._mem,prec,tmp);
}

/*!
 * This method checks that all tuples in \a other are in \a this.
 * If true, the output param \a tupleIds contains the tuples ids of \a this that correspond to tupes in \a this.
 * For each i in [ 0 , other->getNumberOfTuples() ) tuple #i in \a other is equal ( regarding input precision \a prec ) to tuple tupleIds[i] in \a this.
 *
 * \param [in] other - the array having the same number of components than \a this.
 * \param [out] tupleIds - the tuple ids containing the same number of tuples than \a other has.
 * \sa DataArrayDouble::findCommonTuples
 */
bool DataArrayDouble::areIncludedInMe(const DataArrayDouble *other, double prec, DataArrayIdType *&tupleIds) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::areIncludedInMe : input array is NULL !");
  checkAllocated(); other->checkAllocated();
  if(getNumberOfComponents()!=other->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayDouble::areIncludedInMe : the number of components does not match !");
  MCAuto<DataArrayDouble> a=DataArrayDouble::Aggregate(this,other);
  DataArrayIdType *c=0,*ci=0;
  a->findCommonTuples(prec,getNumberOfTuples(),c,ci);
  MCAuto<DataArrayIdType> cSafe(c),ciSafe(ci);
  mcIdType newNbOfTuples=-1;
  MCAuto<DataArrayIdType> ids=DataArrayIdType::ConvertIndexArrayToO2N(a->getNumberOfTuples(),c->begin(),ci->begin(),ci->end(),newNbOfTuples);
  MCAuto<DataArrayIdType> ret1=ids->selectByTupleIdSafeSlice(getNumberOfTuples(),a->getNumberOfTuples(),1);
  tupleIds=ret1.retn();
  return newNbOfTuples==getNumberOfTuples();
}

/*!
 * Searches for tuples coincident within \a prec tolerance. Each tuple is considered
 * as coordinates of a point in getNumberOfComponents()-dimensional space. The
 * distance separating two points is computed with the infinite norm.
 *
 * Indices of coincident tuples are stored in output arrays.
 * A pair of arrays (\a comm, \a commIndex) is called "Surjective Format 2".
 *
 * This method is typically used by MEDCouplingPointSet::findCommonNodes() and
 * MEDCouplingUMesh::mergeNodes().
 *  \param [in] prec - minimal absolute distance between two tuples (infinite norm) at which they are
 *              considered not coincident.
 *  \param [in] limitTupleId - limit tuple id. If all tuples within a group of coincident
 *              tuples have id strictly lower than \a limitTupleId then they are not returned.
 *  \param [out] comm - the array holding ids (== indices) of coincident tuples.
 *               \a comm->getNumberOfComponents() == 1.
 *               \a comm->getNumberOfTuples() == \a commIndex->back().
 *  \param [out] commIndex - the array dividing all indices stored in \a comm into
 *               groups of (indices of) coincident tuples. Its every value is a tuple
 *               index where a next group of tuples begins. For example the second
 *               group of tuples in \a comm is described by following range of indices:
 *               [ \a commIndex[1], \a commIndex[2] ). \a commIndex->getNumberOfTuples()-1
 *               gives the number of groups of coincident tuples.
 *  \throw If \a this is not allocated.
 *  \throw If the number of components is not in [1,2,3,4].
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcdataarraydouble_findcommontuples "Here is a C++ example".
 *
 *  \ref py_mcdataarraydouble_findcommontuples  "Here is a Python example".
 *  \endif
 *  \sa DataArrayInt::ConvertIndexArrayToO2N(), DataArrayDouble::areIncludedInMe
 */
void DataArrayDouble::findCommonTuples(double prec, mcIdType limitTupleId, DataArrayIdType *&comm, DataArrayIdType *&commIndex) const
{
  checkAllocated();
  std::size_t nbOfCompo=getNumberOfComponents();
  if ((nbOfCompo<1) || (nbOfCompo>4)) //test before work
    throw INTERP_KERNEL::Exception("DataArrayDouble::findCommonTuples : Unexpected spacedim of coords. Must be 1, 2, 3 or 4.");

  mcIdType nbOfTuples(getNumberOfTuples());
  //
  MCAuto<DataArrayIdType> c(DataArrayIdType::New()),cI(DataArrayIdType::New()); c->alloc(0,1); cI->pushBackSilent(0);
  switch(nbOfCompo)
  {
    case 4:
      findCommonTuplesAlg<4>(begin(),nbOfTuples,limitTupleId,prec,c,cI);
      break;
    case 3:
      findCommonTuplesAlg<3>(begin(),nbOfTuples,limitTupleId,prec,c,cI);
      break;
    case 2:
      findCommonTuplesAlg<2>(begin(),nbOfTuples,limitTupleId,prec,c,cI);
      break;
    case 1:
      findCommonTuplesAlg<1>(begin(),nbOfTuples,limitTupleId,prec,c,cI);
      break;
    default:
      throw INTERP_KERNEL::Exception("DataArrayDouble::findCommonTuples : nb of components managed are 1,2,3 and 4 ! not implemented for other number of components !");
  }
  comm=c.retn();
  commIndex=cI.retn();
}

/*!
 * This methods returns the minimal distance between the two set of points \a this and \a other.
 * So \a this and \a other have to have the same number of components. If not an INTERP_KERNEL::Exception will be thrown.
 * This method works only if number of components of \a this (equal to those of \a other) is in 1, 2 or 3.
 *
 * \param [out] thisTupleId the tuple id in \a this corresponding to the returned minimal distance
 * \param [out] otherTupleId the tuple id in \a other corresponding to the returned minimal distance
 * \return the minimal distance between the two set of points \a this and \a other.
 * \sa DataArrayDouble::findClosestTupleId
 */
double DataArrayDouble::minimalDistanceTo(const DataArrayDouble *other, mcIdType& thisTupleId, mcIdType& otherTupleId) const
{
  MCAuto<DataArrayIdType> part1=findClosestTupleId(other);
  std::size_t nbOfCompo=getNumberOfComponents();
  mcIdType otherNbTuples=other->getNumberOfTuples();
  const double *thisPt(begin()),*otherPt(other->begin());
  const mcIdType *part1Pt(part1->begin());
  double ret=std::numeric_limits<double>::max();
  for(mcIdType i=0;i<otherNbTuples;i++,part1Pt++,otherPt+=nbOfCompo)
    {
      double tmp(0.);
      for(std::size_t j=0;j<nbOfCompo;j++)
        tmp+=(otherPt[j]-thisPt[nbOfCompo*(*part1Pt)+j])*(otherPt[j]-thisPt[nbOfCompo*(*part1Pt)+j]);
      if(tmp<ret)
        { ret=tmp; thisTupleId=*part1Pt; otherTupleId=i; }
    }
  return sqrt(ret);
}

/*!
 * This methods returns for each tuple in \a other which tuple in \a this is the closest.
 * So \a this and \a other have to have the same number of components. If not an INTERP_KERNEL::Exception will be thrown.
 * This method works only if number of components of \a this (equal to those of \a other) is in 1, 2 or 3.
 *
 * \return a newly allocated (new object to be dealt by the caller) DataArrayInt having \c other->getNumberOfTuples() tuples and one components.
 * \sa DataArrayDouble::minimalDistanceTo
 */
DataArrayIdType *DataArrayDouble::findClosestTupleId(const DataArrayDouble *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::findClosestTupleId : other instance is NULL !");
  checkAllocated(); other->checkAllocated();
  std::size_t nbOfCompo(getNumberOfComponents());
  if(nbOfCompo!=other->getNumberOfComponents())
    {
      std::ostringstream oss; oss << "DataArrayDouble::findClosestTupleId : number of components in this is " << nbOfCompo;
      oss << ", whereas number of components in other is " << other->getNumberOfComponents() << "! Should be equal !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  mcIdType nbOfTuples(other->getNumberOfTuples());
  mcIdType thisNbOfTuples(getNumberOfTuples());
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New(); ret->alloc(nbOfTuples,1);
  double bounds[6];
  getMinMaxPerComponent(bounds);
  switch(nbOfCompo)
  {
    case 3:
      {
        double xDelta(fabs(bounds[1]-bounds[0])),yDelta(fabs(bounds[3]-bounds[2])),zDelta(fabs(bounds[5]-bounds[4]));
        double delta=std::max(xDelta,yDelta); delta=std::max(delta,zDelta);
        double characSize=pow((delta*delta*delta)/((double)thisNbOfTuples),1./3.);
        BBTreePts<3,mcIdType> myTree(begin(),0,0,getNumberOfTuples(),characSize*1e-12);
        FindClosestTupleIdAlg<3>(myTree,3.*characSize*characSize,other->begin(),nbOfTuples,begin(),thisNbOfTuples,ret->getPointer());
        break;
      }
    case 2:
      {
        double xDelta(fabs(bounds[1]-bounds[0])),yDelta(fabs(bounds[3]-bounds[2]));
        double delta=std::max(xDelta,yDelta);
        double characSize=sqrt(delta/(double)thisNbOfTuples);
        BBTreePts<2,mcIdType> myTree(begin(),0,0,getNumberOfTuples(),characSize*1e-12);
        FindClosestTupleIdAlg<2>(myTree,2.*characSize*characSize,other->begin(),nbOfTuples,begin(),thisNbOfTuples,ret->getPointer());
        break;
      }
    case 1:
      {
        double characSize=fabs(bounds[1]-bounds[0])/FromIdType<double>(thisNbOfTuples);
        BBTreePts<1,mcIdType> myTree(begin(),0,0,getNumberOfTuples(),characSize*1e-12);
        FindClosestTupleIdAlg<1>(myTree,1.*characSize*characSize,other->begin(),nbOfTuples,begin(),thisNbOfTuples,ret->getPointer());
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Unexpected spacedim of coords for findClosestTupleId. Must be 1, 2 or 3.");
  }
  return ret.retn();
}

/*!
 * This method expects that \a this and \a otherBBoxFrmt arrays are bounding box arrays ( as the output of MEDCouplingPointSet::getBoundingBoxForBBTree method ).
 * This method will return a DataArrayInt array having the same number of tuples than \a this. This returned array tells for each cell in \a this
 * how many bounding boxes in \a otherBBoxFrmt.
 * So, this method expects that \a this and \a otherBBoxFrmt have the same number of components.
 *
 * \param [in] otherBBoxFrmt - It is an array .
 * \param [in] eps - the absolute precision of the detection. when eps < 0 the bboxes are enlarged so more interactions are detected. Inversely when > 0 the bboxes are stretched.
 * \sa MEDCouplingPointSet::getBoundingBoxForBBTree
 * \throw If \a this and \a otherBBoxFrmt have not the same number of components.
 * \throw If \a this and \a otherBBoxFrmt number of components is not even (BBox format).
 */
DataArrayIdType *DataArrayDouble::computeNbOfInteractionsWith(const DataArrayDouble *otherBBoxFrmt, double eps) const
{
  if(!otherBBoxFrmt)
    throw INTERP_KERNEL::Exception("DataArrayDouble::computeNbOfInteractionsWith : input array is NULL !");
  if(!isAllocated() || !otherBBoxFrmt->isAllocated())
    throw INTERP_KERNEL::Exception("DataArrayDouble::computeNbOfInteractionsWith : this and input array must be allocated !");
  std::size_t nbOfComp(getNumberOfComponents());
  mcIdType nbOfTuples(getNumberOfTuples());
  if(nbOfComp!=otherBBoxFrmt->getNumberOfComponents())
    {
      std::ostringstream oss; oss << "DataArrayDouble::computeNbOfInteractionsWith : this number of components (" << nbOfComp << ") must be equal to the number of components of input array (" << otherBBoxFrmt->getNumberOfComponents() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(nbOfComp%2!=0)
    {
      std::ostringstream oss; oss << "DataArrayDouble::computeNbOfInteractionsWith : Number of components (" << nbOfComp << ") is not even ! It should be to be compatible with bbox format !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(nbOfTuples,1);
  const double *thisBBPtr(begin());
  mcIdType *retPtr(ret->getPointer());
  switch(nbOfComp/2)
  {
    case 3:
      {
        BBTree<3,mcIdType> bbt(otherBBoxFrmt->begin(),0,0,otherBBoxFrmt->getNumberOfTuples(),eps);
        for(mcIdType i=0;i<nbOfTuples;i++,retPtr++,thisBBPtr+=nbOfComp)
          *retPtr=bbt.getNbOfIntersectingElems(thisBBPtr);
        break;
      }
    case 2:
      {
        BBTree<2,mcIdType> bbt(otherBBoxFrmt->begin(),0,0,otherBBoxFrmt->getNumberOfTuples(),eps);
        for(mcIdType i=0;i<nbOfTuples;i++,retPtr++,thisBBPtr+=nbOfComp)
          *retPtr=bbt.getNbOfIntersectingElems(thisBBPtr);
        break;
      }
    case 1:
      {
        BBTree<1,mcIdType> bbt(otherBBoxFrmt->begin(),0,0,otherBBoxFrmt->getNumberOfTuples(),eps);
        for(mcIdType i=0;i<nbOfTuples;i++,retPtr++,thisBBPtr+=nbOfComp)
          *retPtr=bbt.getNbOfIntersectingElems(thisBBPtr);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("DataArrayDouble::computeNbOfInteractionsWith : space dimension supported are [1,2,3] !");
  }

  return ret.retn();
}

/*!
 * Returns a copy of \a this array by excluding coincident tuples. Each tuple is
 * considered as coordinates of a point in getNumberOfComponents()-dimensional
 * space. The distance between tuples is computed using norm2. If several tuples are
 * not far each from other than \a prec, only one of them remains in the result
 * array. The order of tuples in the result array is same as in \a this one except
 * that coincident tuples are excluded.
 *  \param [in] prec - minimal absolute distance between two tuples at which they are
 *              considered not coincident.
 *  \param [in] limitTupleId - limit tuple id. If all tuples within a group of coincident
 *              tuples have id strictly lower than \a limitTupleId then they are not excluded.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If the number of components is not in [1,2,3,4].
 *
 *  \if ENABLE_EXAMPLES
 *  \ref py_mcdataarraydouble_getdifferentvalues "Here is a Python example".
 *  \endif
 */
DataArrayDouble *DataArrayDouble::getDifferentValues(double prec, mcIdType limitTupleId) const
{
  checkAllocated();
  DataArrayIdType *c0=0,*cI0=0;
  findCommonTuples(prec,limitTupleId,c0,cI0);
  MCAuto<DataArrayIdType> c(c0),cI(cI0);
  mcIdType newNbOfTuples=-1;
  MCAuto<DataArrayIdType> o2n=DataArrayIdType::ConvertIndexArrayToO2N(getNumberOfTuples(),c0->begin(),cI0->begin(),cI0->end(),newNbOfTuples);
  return renumberAndReduce(o2n->getConstPointer(),newNbOfTuples);
}

/*!
 * Copy all components in a specified order from another DataArrayDouble.
 * Both numerical and textual data is copied. The number of tuples in \a this and
 * the other array can be different.
 *  \param [in] a - the array to copy data from.
 *  \param [in] compoIds - sequence of zero based indices of components, data of which is
 *              to be copied.
 *  \throw If \a a is NULL.
 *  \throw If \a compoIds.size() != \a a->getNumberOfComponents().
 *  \throw If \a compoIds[i] < 0 or \a compoIds[i] > \a this->getNumberOfComponents().
 *
 *  \if ENABLE_EXAMPLES
 *  \ref py_mcdataarraydouble_setselectedcomponents "Here is a Python example".
 *  \endif
 */
void DataArrayDouble::setSelectedComponents(const DataArrayDouble *a, const std::vector<std::size_t>& compoIds)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayDouble::setSelectedComponents : input DataArrayDouble is NULL !");
  checkAllocated();
  copyPartOfStringInfoFrom2(compoIds,*a);
  std::size_t partOfCompoSz=compoIds.size();
  std::size_t nbOfCompo=getNumberOfComponents();
  mcIdType nbOfTuples=std::min(getNumberOfTuples(),a->getNumberOfTuples());
  const double *ac=a->getConstPointer();
  double *nc=getPointer();
  for(mcIdType i=0;i<nbOfTuples;i++)
    for(std::size_t j=0;j<partOfCompoSz;j++,ac++)
      nc[nbOfCompo*i+compoIds[j]]=*ac;
}
/*!
 * Checks if 0.0 value is present in \a this array. If it is the case, an exception
 * is thrown.
 * \throw If zero is found in \a this array.
 */
void DataArrayDouble::checkNoNullValues() const
{
  const double *tmp=getConstPointer();
  mcIdType nbOfElems=getNbOfElems();
  const double *where=std::find(tmp,tmp+nbOfElems,0.);
  if(where!=tmp+nbOfElems)
    throw INTERP_KERNEL::Exception("A value 0.0 have been detected !");
}

/*!
 * Computes minimal and maximal value in each component. An output array is filled
 * with \c 2 * \a this->getNumberOfComponents() values, so the caller is to allocate
 * enough memory before calling this method.
 *  \param [out] bounds - array of size at least 2 *\a this->getNumberOfComponents().
 *               It is filled as follows:<br>
 *               \a bounds[0] = \c min_of_component_0 <br>
 *               \a bounds[1] = \c max_of_component_0 <br>
 *               \a bounds[2] = \c min_of_component_1 <br>
 *               \a bounds[3] = \c max_of_component_1 <br>
 *               ...
 */
void DataArrayDouble::getMinMaxPerComponent(double *bounds) const
{
  checkAllocated();
  std::size_t dim=getNumberOfComponents();
  for (std::size_t idim=0; idim<dim; idim++)
    {
      bounds[idim*2]=std::numeric_limits<double>::max();
      bounds[idim*2+1]=-std::numeric_limits<double>::max();
    }
  const double *ptr=getConstPointer();
  mcIdType nbOfTuples=getNumberOfTuples();
  for(mcIdType i=0;i<nbOfTuples;i++)
    {
      for(std::size_t idim=0;idim<dim;idim++)
        {
          if(bounds[idim*2]>ptr[i*dim+idim])
            {
              bounds[idim*2]=ptr[i*dim+idim];
            }
          if(bounds[idim*2+1]<ptr[i*dim+idim])
            {
              bounds[idim*2+1]=ptr[i*dim+idim];
            }
        }
    }
}

/*!
 * This method retrieves a newly allocated DataArrayDouble instance having same number of tuples than \a this and twice number of components than \a this
 * to store both the min and max per component of each tuples.
 * \param [in] epsilon the width of the bbox (identical in each direction) - 0.0 by default
 *
 * \return a newly created DataArrayDouble instance having \c this->getNumberOfTuples() tuples and 2 * \c this->getNumberOfComponent() components
 *
 * \throw If \a this is not allocated yet.
 */
DataArrayDouble *DataArrayDouble::computeBBoxPerTuple(double epsilon) const
{
  checkAllocated();
  const double *dataPtr=getConstPointer();
  std::size_t nbOfCompo=getNumberOfComponents();
  mcIdType nbTuples=getNumberOfTuples();
  MCAuto<DataArrayDouble> bbox=DataArrayDouble::New();
  bbox->alloc(nbTuples,2*nbOfCompo);
  double *bboxPtr=bbox->getPointer();
  for(mcIdType i=0;i<nbTuples;i++)
    {
      for(std::size_t j=0;j<nbOfCompo;j++)
        {
          bboxPtr[2*nbOfCompo*i+2*j]=dataPtr[nbOfCompo*i+j]-epsilon;
          bboxPtr[2*nbOfCompo*i+2*j+1]=dataPtr[nbOfCompo*i+j]+epsilon;
        }
    }
  return bbox.retn();
}

/*!
 * For each tuples **t** in \a other, this method retrieves tuples in \a this that are equal to **t**.
 * Two tuples are considered equal if the euclidian distance between the two tuples is lower than \a eps.
 *
 * \param [in] other a DataArrayDouble having same number of components than \a this.
 * \param [in] eps absolute precision representing distance (using infinite norm) between 2 tuples behind which 2 tuples are considered equal.
 * \param [out] c will contain the set of tuple ids in \a this that are equal to to the tuple ids in \a other contiguously.
 *             \a cI allows to extract information in \a c.
 * \param [out] cI is an indirection array that allows to extract the data contained in \a c.
 *
 * \throw In case of:
 *  - \a this is not allocated
 *  - \a other is not allocated or null
 *  - \a this and \a other do not have the same number of components
 *  - if number of components of \a this is not in [1,2,3]
 *
 * \sa MEDCouplingPointSet::getNodeIdsNearPoints, DataArrayDouble::getDifferentValues
 */
void DataArrayDouble::computeTupleIdsNearTuples(const DataArrayDouble *other, double eps, DataArrayIdType *& c, DataArrayIdType *& cI) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::computeTupleIdsNearTuples : input pointer other is null !");
  checkAllocated();
  other->checkAllocated();
  std::size_t nbOfCompo=getNumberOfComponents();
  std::size_t otherNbOfCompo=other->getNumberOfComponents();
  if(nbOfCompo!=otherNbOfCompo)
    throw INTERP_KERNEL::Exception("DataArrayDouble::computeTupleIdsNearTuples : number of components should be equal between this and other !");
  mcIdType nbOfTuplesOther=other->getNumberOfTuples();
  mcIdType nbOfTuples=getNumberOfTuples();
  MCAuto<DataArrayIdType> cArr(DataArrayIdType::New()),cIArr(DataArrayIdType::New()); cArr->alloc(0,1); cIArr->pushBackSilent(0);
  switch(nbOfCompo)
  {
    case 3:
      {
        BBTreePts<3,mcIdType> myTree(begin(),0,0,nbOfTuples,eps);
        FindTupleIdsNearTuplesAlg<3>(myTree,other->getConstPointer(),nbOfTuplesOther,eps,cArr,cIArr);
        break;
      }
    case 2:
      {
        BBTreePts<2,mcIdType> myTree(begin(),0,0,nbOfTuples,eps);
        FindTupleIdsNearTuplesAlg<2>(myTree,other->getConstPointer(),nbOfTuplesOther,eps,cArr,cIArr);
        break;
      }
    case 1:
      {
        BBTreePts<1,mcIdType> myTree(begin(),0,0,nbOfTuples,eps);
        FindTupleIdsNearTuplesAlg<1>(myTree,other->getConstPointer(),nbOfTuplesOther,eps,cArr,cIArr);
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("Unexpected spacedim of coords for computeTupleIdsNearTuples. Must be 1, 2 or 3.");
  }
  c=cArr.retn(); cI=cIArr.retn();
}

/*!
 * This method recenter tuples in \b this in order to be centered at the origin to benefit about the advantages of maximal precision to be around the box
 * around origin of 'radius' 1.
 *
 * \param [in] eps absolute epsilon. under that value of delta between max and min no scale is performed.
 */
void DataArrayDouble::recenterForMaxPrecision(double eps)
{
  checkAllocated();
  std::size_t dim=getNumberOfComponents();
  std::vector<double> bounds(2*dim);
  getMinMaxPerComponent(&bounds[0]);
  for(std::size_t i=0;i<dim;i++)
    {
      double delta=bounds[2*i+1]-bounds[2*i];
      double offset=(bounds[2*i]+bounds[2*i+1])/2.;
      if(delta>eps)
        applyLin(1./delta,-offset/delta,i);
      else
        applyLin(1.,-offset,i);
    }
}

/*!
 * Returns the maximal value and all its locations within \a this one-dimensional array.
 *  \param [out] tupleIds - a new instance of DataArrayInt containing indices of
 *               tuples holding the maximal value. The caller is to delete it using
 *               decrRef() as it is no more needed.
 *  \return double - the maximal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
double DataArrayDouble::getMaxValue2(DataArrayIdType*& tupleIds) const
{
  mcIdType tmp;
  tupleIds=0;
  double ret=getMaxValue(tmp);
  tupleIds=findIdsInRange(ret,ret);
  return ret;
}

/*!
 * Returns the minimal value and all its locations within \a this one-dimensional array.
 *  \param [out] tupleIds - a new instance of DataArrayInt containing indices of
 *               tuples holding the minimal value. The caller is to delete it using
 *               decrRef() as it is no more needed.
 *  \return double - the minimal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
double DataArrayDouble::getMinValue2(DataArrayIdType*& tupleIds) const
{
  mcIdType tmp;
  tupleIds=0;
  double ret=getMinValue(tmp);
  tupleIds=findIdsInRange(ret,ret);
  return ret;
}

/*!
 * This method returns the number of values in \a this that are equals ( within an absolute precision of \a eps ) to input parameter \a value.
 * This method only works for single component array.
 *
 * \return a value in [ 0, \c this->getNumberOfTuples() )
 *
 * \throw If \a this is not allocated
 *
 */
mcIdType DataArrayDouble::count(double value, double eps) const
{
  mcIdType ret=0;
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::count : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before !");
  const double *vals=begin();
  mcIdType nbOfTuples=getNumberOfTuples();
  for(mcIdType i=0;i<nbOfTuples;i++,vals++)
    if(fabs(*vals-value)<=eps)
      ret++;
  return ret;
}

/*!
 * Returns the average value of \a this one-dimensional array.
 *  \return double - the average value over all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
double DataArrayDouble::getAverageValue() const
{
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getAverageValue : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before !");
  mcIdType nbOfTuples(getNumberOfTuples());
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::getAverageValue : array exists but number of tuples must be > 0 !");
  const double *vals=getConstPointer();
  double ret=std::accumulate(vals,vals+nbOfTuples,0.);
  return ret/FromIdType<double>(nbOfTuples);
}

/*!
 * Returns the Euclidean norm of the vector defined by \a this array.
 *  \return double - the value of the Euclidean norm, i.e.
 *          the square root of the inner product of vector.
 *  \throw If \a this is not allocated.
 */
double DataArrayDouble::norm2() const
{
  checkAllocated();
  double ret=0.;
  std::size_t nbOfElems=getNbOfElems();
  const double *pt=getConstPointer();
  for(std::size_t i=0;i<nbOfElems;i++,pt++)
    ret+=(*pt)*(*pt);
  return sqrt(ret);
}

/*!
 * Returns the maximum norm of the vector defined by \a this array.
 * This method works even if the number of components is different from one.
 * If the number of elements in \a this is 0, -1. is returned.
 *  \return double - the value of the maximum norm, i.e.
 *          the maximal absolute value among values of \a this array (whatever its number of components).
 *  \throw If \a this is not allocated.
 */
double DataArrayDouble::normMax() const
{
  checkAllocated();
  double ret(-1.);
  std::size_t nbOfElems(getNbOfElems());
  const double *pt(getConstPointer());
  for(std::size_t i=0;i<nbOfElems;i++,pt++)
    {
      double val(std::abs(*pt));
      if(val>ret)
        ret=val;
    }
  return ret;
}

/*!
 * Returns the maximum norm of for each component of \a this array.
 * If the number of elements in \a this is 0, -1. is returned.
*  \param [out] res - pointer to an array of result values, of size at least \a
 *         this->getNumberOfComponents(), that is to be allocated by the caller.
 *  \throw If \a this is not allocated.
 */
void DataArrayDouble::normMaxPerComponent(double * res) const
{
  checkAllocated();
  mcIdType nbOfTuples(getNumberOfTuples());
  std::size_t nbOfCompos(getNumberOfComponents());
  std::fill(res, res+nbOfCompos, -1.0);
  const double *pt(getConstPointer());
  for(mcIdType i=0;i<nbOfTuples;i++)
    for (std::size_t j=0; j<nbOfCompos; j++, pt++)
      {
        double val(std::abs(*pt));
        if(val>res[j])
          res[j]=val;
      }
}


/*!
 * Returns the minimum norm (absolute value) of the vector defined by \a this array.
 * This method works even if the number of components is different from one.
 * If the number of elements in \a this is 0, std::numeric_limits<double>::max() is returned.
 *  \return double - the value of the minimum norm, i.e.
 *          the minimal absolute value among values of \a this array (whatever its number of components).
 *  \throw If \a this is not allocated.
 */
double DataArrayDouble::normMin() const
{
  checkAllocated();
  double ret(std::numeric_limits<double>::max());
  std::size_t nbOfElems(getNbOfElems());
  const double *pt(getConstPointer());
  for(std::size_t i=0;i<nbOfElems;i++,pt++)
    {
      double val(std::abs(*pt));
      if(val<ret)
        ret=val;
    }
  return ret;
}

/*!
 * Accumulates values of each component of \a this array.
 *  \param [out] res - an array of length \a this->getNumberOfComponents(), allocated
 *         by the caller, that is filled by this method with sum value for each
 *         component.
 *  \throw If \a this is not allocated.
 */
void DataArrayDouble::accumulate(double *res) const
{
  checkAllocated();
  const double *ptr=getConstPointer();
  mcIdType nbTuple(getNumberOfTuples());
  std::size_t nbComps(getNumberOfComponents());
  std::fill(res,res+nbComps,0.);
  for(mcIdType i=0;i<nbTuple;i++)
    std::transform(ptr+i*nbComps,ptr+(i+1)*nbComps,res,res,std::plus<double>());
}

/*!
 * This method returns the min distance from an external tuple defined by [ \a tupleBg , \a tupleEnd ) to \a this and
 * the first tuple in \a this that matches the returned distance. If there is no tuples in \a this an exception will be thrown.
 *
 *
 * \a this is expected to be allocated and expected to have a number of components equal to the distance from \a tupleBg to
 * \a tupleEnd. If not an exception will be thrown.
 *
 * \param [in] tupleBg start pointer (included) of input external tuple
 * \param [in] tupleEnd end pointer (not included) of input external tuple
 * \param [out] tupleId the tuple id in \a this that matches the min of distance between \a this and input external tuple
 * \return the min distance.
 * \sa MEDCouplingUMesh::distanceToPoint
 */
double DataArrayDouble::distanceToTuple(const double *tupleBg, const double *tupleEnd, mcIdType& tupleId) const
{
  checkAllocated();
  mcIdType nbTuple(getNumberOfTuples());
  std::size_t nbComps(getNumberOfComponents());
  if(nbComps!=(std::size_t)std::distance(tupleBg,tupleEnd))
    { std::ostringstream oss; oss << "DataArrayDouble::distanceToTuple : size of input tuple is " << std::distance(tupleBg,tupleEnd) << " should be equal to the number of components in this : " << nbComps << " !"; throw INTERP_KERNEL::Exception(oss.str().c_str()); }
  if(nbTuple==0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::distanceToTuple : no tuple in this ! No distance to compute !");
  double ret0=std::numeric_limits<double>::max();
  tupleId=-1;
  const double *work=getConstPointer();
  for(mcIdType i=0;i<nbTuple;i++)
    {
      double val=0.;
      for(std::size_t j=0;j<nbComps;j++,work++)
        val+=(*work-tupleBg[j])*((*work-tupleBg[j]));
      if(val>=ret0)
        continue;
      else
        { ret0=val; tupleId=i; }
    }
  return sqrt(ret0);
}

/*!
 * Accumulate values of the given component of \a this array.
 *  \param [in] compId - the index of the component of interest.
 *  \return double - a sum value of \a compId-th component.
 *  \throw If \a this is not allocated.
 *  \throw If \a the condition ( 0 <= \a compId < \a this->getNumberOfComponents() ) is
 *         not respected.
 */
double DataArrayDouble::accumulate(std::size_t compId) const
{
  checkAllocated();
  const double *ptr=getConstPointer();
  mcIdType nbTuple(getNumberOfTuples());
  std::size_t nbComps(getNumberOfComponents());
  if(compId>=nbComps)
    throw INTERP_KERNEL::Exception("DataArrayDouble::accumulate : Invalid compId specified : No such nb of components !");
  double ret=0.;
  for(mcIdType i=0;i<nbTuple;i++)
    ret+=ptr[i*nbComps+compId];
  return ret;
}

/*!
 * This method accumulate using addition tuples in \a this using input index array [ \a bgOfIndex, \a endOfIndex ).
 * The returned array will have same number of components than \a this and number of tuples equal to
 * \c std::distance(bgOfIndex,endOfIndex) \b minus \b one.
 *
 * The input index array is expected to be ascendingly sorted in which the all referenced ids should be in [0, \c this->getNumberOfTuples).
 * This method is quite useful for users that need to put a field on cells to field on nodes on the same mesh without a need of conservation.
 *
 * \param [in] bgOfIndex - begin (included) of the input index array.
 * \param [in] endOfIndex - end (excluded) of the input index array.
 * \return DataArrayDouble * - the new instance having the same number of components than \a this.
 *
 * \throw If bgOfIndex or end is NULL.
 * \throw If input index array is not ascendingly sorted.
 * \throw If there is an id in [ \a bgOfIndex, \a endOfIndex ) not in [0, \c this->getNumberOfTuples).
 * \throw If std::distance(bgOfIndex,endOfIndex)==0.
 */
DataArrayDouble *DataArrayDouble::accumulatePerChunck(const mcIdType *bgOfIndex, const mcIdType *endOfIndex) const
{
  if(!bgOfIndex || !endOfIndex)
    throw INTERP_KERNEL::Exception("DataArrayDouble::accumulatePerChunck : input pointer NULL !");
  checkAllocated();
  std::size_t nbCompo(getNumberOfComponents());
  mcIdType nbOfTuples(getNumberOfTuples());
  std::size_t sz=std::distance(bgOfIndex,endOfIndex);
  if(sz<1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::accumulatePerChunck : invalid size of input index array !");
  sz--;
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New(); ret->alloc(sz,nbCompo);
  const mcIdType *w=bgOfIndex;
  if(*w<0 || *w>=nbOfTuples)
    throw INTERP_KERNEL::Exception("DataArrayDouble::accumulatePerChunck : The first element of the input index not in [0,nbOfTuples) !");
  const double *srcPt=begin()+(*w)*nbCompo;
  double *tmp=ret->getPointer();
  for(std::size_t i=0;i<sz;i++,tmp+=nbCompo,w++)
    {
      std::fill(tmp,tmp+nbCompo,0.);
      if(w[1]>=w[0])
        {
          for(mcIdType j=w[0];j<w[1];j++,srcPt+=nbCompo)
            {
              if(j>=0 && j<nbOfTuples)
                std::transform(srcPt,srcPt+nbCompo,tmp,tmp,std::plus<double>());
              else
                {
                  std::ostringstream oss; oss << "DataArrayDouble::accumulatePerChunck : At rank #" << i << " the input index array points to id " << j << " should be in [0," << nbOfTuples << ") !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::accumulatePerChunck : At rank #" << i << " the input index array is not in ascendingly sorted.";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * This method is close to numpy cumSum except that number of element is equal to \a this->getNumberOfTuples()+1. First element of DataArray returned is equal to 0.
 * This method expects that \a this as only one component. The returned array will have \a this->getNumberOfTuples()+1 tuple with also one component.
 * The ith element of returned array is equal to the sum of elements in \a this with rank strictly lower than i.
 *
 * \return DataArrayDouble - A newly built array containing cum sum of \a this.
 */
MCAuto<DataArrayDouble> DataArrayDouble::cumSum() const
{
  checkAllocated();
  checkNbOfComps(1,"DataArrayDouble::cumSum : this is expected to be single component");
  mcIdType nbOfTuple(getNumberOfTuples());
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(nbOfTuple+1,1);
  double *ptr(ret->getPointer());
  ptr[0]=0.;
  const double *thisPtr(begin());
  for(mcIdType i=0;i<nbOfTuple;i++)
    ptr[i+1]=ptr[i]+thisPtr[i];
  return ret;
}

/*!
 * Converts each 2D point defined by the tuple of \a this array from the Polar to the
 * Cartesian coordinate system. The two components of the tuple of \a this array are
 * considered to contain (1) radius and (2) angle of the point in the Polar CS.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          contains X and Y coordinates of the point in the Cartesian CS. The caller
 *          is to delete this array using decrRef() as it is no more needed. The array
 *          does not contain any textual info on components.
 *  \throw If \a this->getNumberOfComponents() != 2.
 * \sa fromCartToPolar
 */
DataArrayDouble *DataArrayDouble::fromPolarToCart() const
{
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  if(nbOfComp!=2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromPolarToCart : must be an array with exactly 2 components !");
  mcIdType nbOfTuple(getNumberOfTuples());
  DataArrayDouble *ret(DataArrayDouble::New());
  ret->alloc(nbOfTuple,2);
  double *w(ret->getPointer());
  const double *wIn(getConstPointer());
  for(mcIdType i=0;i<nbOfTuple;i++,w+=2,wIn+=2)
    {
      w[0]=wIn[0]*cos(wIn[1]);
      w[1]=wIn[0]*sin(wIn[1]);
    }
  return ret;
}

/*!
 * Converts each 3D point defined by the tuple of \a this array from the Cylindrical to
 * the Cartesian coordinate system. The three components of the tuple of \a this array
 * are considered to contain (1) radius, (2) azimuth and (3) altitude of the point in
 * the Cylindrical CS.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          contains X, Y and Z coordinates of the point in the Cartesian CS. The info
 *          on the third component is copied from \a this array. The caller
 *          is to delete this array using decrRef() as it is no more needed.
 *  \throw If \a this->getNumberOfComponents() != 3.
 * \sa fromCartToCyl
 */
DataArrayDouble *DataArrayDouble::fromCylToCart() const
{
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  if(nbOfComp!=3)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCylToCart : must be an array with exactly 3 components !");
  mcIdType nbOfTuple(getNumberOfTuples());
  DataArrayDouble *ret(DataArrayDouble::New());
  ret->alloc(getNumberOfTuples(),3);
  double *w(ret->getPointer());
  const double *wIn(getConstPointer());
  for(mcIdType i=0;i<nbOfTuple;i++,w+=3,wIn+=3)
    {
      w[0]=wIn[0]*cos(wIn[1]);
      w[1]=wIn[0]*sin(wIn[1]);
      w[2]=wIn[2];
    }
  ret->setInfoOnComponent(2,getInfoOnComponent(2));
  return ret;
}

/*!
 * Converts each 3D point defined by the tuple of \a this array from the Spherical to
 * the Cartesian coordinate system. The three components of the tuple of \a this array
 * are considered to contain (1) radius, (2) polar angle and (3) azimuthal angle of the
 * point in the Cylindrical CS.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          contains X, Y and Z coordinates of the point in the Cartesian CS. The info
 *          on the third component is copied from \a this array. The caller
 *          is to delete this array using decrRef() as it is no more needed.
 *  \throw If \a this->getNumberOfComponents() != 3.
 * \sa fromCartToSpher
 */
DataArrayDouble *DataArrayDouble::fromSpherToCart() const
{
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  if(nbOfComp!=3)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromSpherToCart : must be an array with exactly 3 components !");
  mcIdType nbOfTuple(getNumberOfTuples());
  DataArrayDouble *ret(DataArrayDouble::New());
  ret->alloc(getNumberOfTuples(),3);
  double *w(ret->getPointer());
  const double *wIn(getConstPointer());
  for(mcIdType i=0;i<nbOfTuple;i++,w+=3,wIn+=3)
    {
      w[0]=wIn[0]*cos(wIn[2])*sin(wIn[1]);
      w[1]=wIn[0]*sin(wIn[2])*sin(wIn[1]);
      w[2]=wIn[0]*cos(wIn[1]);
    }
  return ret;
}

/*!
 * This method returns a new array containing the same number of tuples than \a this. To do this, this method needs \a at parameter to specify the convention of \a this.
 * All the tuples of the returned array will be in cartesian sense. So if \a at equals to AX_CART the returned array is basically a deep copy of \a this.
 * If \a at equals to AX_CYL the returned array will be the result of operation cylindric to cartesian of \a this...
 *
 * \param [in] atOfThis - The axis type of \a this.
 * \return DataArrayDouble * - the new instance of DataArrayDouble (that must be dealed by caller) containing the result of the cartesianizification of \a this.
 */
DataArrayDouble *DataArrayDouble::cartesianize(MEDCouplingAxisType atOfThis) const
{
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  MCAuto<DataArrayDouble> ret;
  switch(atOfThis)
    {
    case AX_CART:
      ret=deepCopy();
    case AX_CYL:
      if(nbOfComp==3)
        {
          ret=fromCylToCart();
          break;
        }
      if(nbOfComp==2)
        {
          ret=fromPolarToCart();
          break;
        }
      else
        throw INTERP_KERNEL::Exception("DataArrayDouble::cartesianize : For AX_CYL, number of components must be in [2,3] !");
    case AX_SPHER:
      if(nbOfComp==3)
        {
          ret=fromSpherToCart();
          break;
        }
      if(nbOfComp==2)
        {
          ret=fromPolarToCart();
          break;
        }
      else
        throw INTERP_KERNEL::Exception("DataArrayDouble::cartesianize : For AX_CYL, number of components must be in [2,3] !");
    default:
      throw INTERP_KERNEL::Exception("DataArrayDouble::cartesianize : not recognized axis type ! Only AX_CART, AX_CYL and AX_SPHER supported !");
    }
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * This method returns a newly created array to be deallocated that contains the result of conversion from cartesian to polar.
 * This method expects that \a this has exactly 2 components.
 * \sa fromPolarToCart
 */
DataArrayDouble *DataArrayDouble::fromCartToPolar() const
{
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  mcIdType nbTuples(getNumberOfTuples());
  if(nbOfComp!=2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCartToPolar : must be an array with exactly 2 components !");
  ret->alloc(nbTuples,2);
  double *retPtr(ret->getPointer());
  const double *ptr(begin());
  for(mcIdType i=0;i<nbTuples;i++,ptr+=2,retPtr+=2)
    {
      retPtr[0]=sqrt(ptr[0]*ptr[0]+ptr[1]*ptr[1]);
      retPtr[1]=atan2(ptr[1],ptr[0]);
    }
  return ret.retn();
}

/*!
 * This method returns a newly created array to be deallocated that contains the result of conversion from cartesian to cylindrical.
 * This method expects that \a this has exactly 3 components.
 * \sa fromCylToCart
 */
DataArrayDouble *DataArrayDouble::fromCartToCyl() const
{
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  mcIdType nbTuples(getNumberOfTuples());
  if(nbOfComp!=3)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCartToCyl : must be an array with exactly 3 components !");
  ret->alloc(nbTuples,3);
  double *retPtr(ret->getPointer());
  const double *ptr(begin());
  for(mcIdType i=0;i<nbTuples;i++,ptr+=3,retPtr+=3)
    {
      retPtr[0]=sqrt(ptr[0]*ptr[0]+ptr[1]*ptr[1]);
      retPtr[1]=atan2(ptr[1],ptr[0]);
      retPtr[2]=ptr[2];
    }
  return ret.retn();
}

/*!
 * This method returns a newly created array to be deallocated that contains the result of conversion from cartesian to spherical coordinates.
 * \sa fromSpherToCart
 */
DataArrayDouble *DataArrayDouble::fromCartToSpher() const
{
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  mcIdType nbTuples(getNumberOfTuples());
  if(nbOfComp!=3)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCartToSpher : must be an array with exactly 3 components !");
  ret->alloc(nbTuples,3);
  double *retPtr(ret->getPointer());
  const double *ptr(begin());
  for(mcIdType i=0;i<nbTuples;i++,ptr+=3,retPtr+=3)
    {
      retPtr[0]=sqrt(ptr[0]*ptr[0]+ptr[1]*ptr[1]+ptr[2]*ptr[2]);
      retPtr[1]=acos(ptr[2]/retPtr[0]);
      retPtr[2]=atan2(ptr[1],ptr[0]);
    }
  return ret.retn();
}

/*!
 * This method returns a newly created array to be deallocated that contains the result of conversion from cartesian to cylindrical relative to the given \a center and a \a vector.
 * This method expects that \a this has exactly 3 components.
 * \sa MEDCouplingFieldDouble::computeVectorFieldCyl
 */
DataArrayDouble *DataArrayDouble::fromCartToCylGiven(const DataArrayDouble *coords, const double center[3], const double vect[3]) const
{
  if(!coords)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCartToCylGiven : input coords are NULL !");
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  checkAllocated(); coords->checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  mcIdType nbTuples(getNumberOfTuples());
  if(nbOfComp!=3)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCartToCylGiven : must be an array with exactly 3 components !");
  if(coords->getNumberOfComponents()!=3)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCartToCylGiven : coords array must have exactly 3 components !");
  if(coords->getNumberOfTuples()!=nbTuples)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCartToCylGiven : coords array must have the same number of tuples !");
  ret->alloc(nbTuples,nbOfComp);
  double magOfVect(sqrt(vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]));
  if(magOfVect<1e-12)
    throw INTERP_KERNEL::Exception("DataArrayDouble::fromCartToCylGiven : magnitude of vect is too low !");
  double Ur[3],Uteta[3],Uz[3],*retPtr(ret->getPointer());
  const double *coo(coords->begin()),*vectField(begin());
  std::transform(vect,vect+3,Uz,std::bind2nd(std::multiplies<double>(),1./magOfVect));
  for(mcIdType i=0;i<nbTuples;i++,vectField+=3,retPtr+=3,coo+=3)
    {
      std::transform(coo,coo+3,center,Ur,std::minus<double>());
      Uteta[0]=Uz[1]*Ur[2]-Uz[2]*Ur[1]; Uteta[1]=Uz[2]*Ur[0]-Uz[0]*Ur[2]; Uteta[2]=Uz[0]*Ur[1]-Uz[1]*Ur[0];
      double magOfTeta(sqrt(Uteta[0]*Uteta[0]+Uteta[1]*Uteta[1]+Uteta[2]*Uteta[2]));
      std::transform(Uteta,Uteta+3,Uteta,std::bind2nd(std::multiplies<double>(),1./magOfTeta));
      Ur[0]=Uteta[1]*Uz[2]-Uteta[2]*Uz[1]; Ur[1]=Uteta[2]*Uz[0]-Uteta[0]*Uz[2]; Ur[2]=Uteta[0]*Uz[1]-Uteta[1]*Uz[0];
      retPtr[0]=Ur[0]*vectField[0]+Ur[1]*vectField[1]+Ur[2]*vectField[2];
      retPtr[1]=Uteta[0]*vectField[0]+Uteta[1]*vectField[1]+Uteta[2]*vectField[2];
      retPtr[2]=Uz[0]*vectField[0]+Uz[1]*vectField[1]+Uz[2]*vectField[2];
    }
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Computes the doubly contracted product of every tensor defined by the tuple of \a this
 * array containing 6 components.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          is calculated from the tuple <em>(t)</em> of \a this array as follows:
 *          \f$ t[0]^2+t[1]^2+t[2]^2+2*t[3]^2+2*t[4]^2+2*t[5]^2\f$.
 *         The caller is to delete this result array using decrRef() as it is no more needed. 
 *  \throw If \a this->getNumberOfComponents() != 6.
 */
DataArrayDouble *DataArrayDouble::doublyContractedProduct() const
{
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  if(nbOfComp!=6)
    throw INTERP_KERNEL::Exception("DataArrayDouble::doublyContractedProduct : must be an array with exactly 6 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  mcIdType nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(mcIdType i=0;i<nbOfTuple;i++,dest++,src+=6)
    *dest=src[0]*src[0]+src[1]*src[1]+src[2]*src[2]+2.*src[3]*src[3]+2.*src[4]*src[4]+2.*src[5]*src[5];
  return ret;
}

/*!
 * Computes the determinant of every square matrix defined by the tuple of \a this
 * array, which contains either 4, 6 or 9 components. The case of 6 components
 * corresponds to that of the upper triangular matrix.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble, whose each tuple
 *          is the determinant of matrix of the corresponding tuple of \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() is not in [4,6,9].
 */
DataArrayDouble *DataArrayDouble::determinant() const
{
  checkAllocated();
  DataArrayDouble *ret=DataArrayDouble::New();
  mcIdType nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  switch(getNumberOfComponents())
  {
    case 6:
      for(mcIdType i=0;i<nbOfTuple;i++,dest++,src+=6)
        *dest=src[0]*src[1]*src[2]+2.*src[4]*src[5]*src[3]-src[0]*src[4]*src[4]-src[2]*src[3]*src[3]-src[1]*src[5]*src[5];
      return ret;
    case 4:
      for(mcIdType i=0;i<nbOfTuple;i++,dest++,src+=4)
        *dest=src[0]*src[3]-src[1]*src[2];
      return ret;
    case 9:
      for(mcIdType i=0;i<nbOfTuple;i++,dest++,src+=9)
        *dest=src[0]*src[4]*src[8]+src[1]*src[5]*src[6]+src[2]*src[3]*src[7]-src[0]*src[5]*src[7]-src[1]*src[3]*src[8]-src[2]*src[4]*src[6];
      return ret;
    default:
      ret->decrRef();
      throw INTERP_KERNEL::Exception("DataArrayDouble::determinant : Invalid number of components ! must be in 4,6,9 !");
  }
}

/*!
 * Computes 3 eigenvalues of every upper triangular matrix defined by the tuple of
 * \a this array, which contains 6 components.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing 3
 *          components, whose each tuple contains the eigenvalues of the matrix of
 *          corresponding tuple of \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() != 6.
 */
DataArrayDouble *DataArrayDouble::eigenValues() const
{
  checkAllocated();
  std::size_t nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6)
    throw INTERP_KERNEL::Exception("DataArrayDouble::eigenValues : must be an array with exactly 6 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  mcIdType nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,3);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(mcIdType i=0;i<nbOfTuple;i++,dest+=3,src+=6)
    INTERP_KERNEL::computeEigenValues6(src,dest);
  return ret;
}

/*!
 * Computes 3 eigenvectors of every upper triangular matrix defined by the tuple of
 * \a this array, which contains 6 components.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing 9
 *          components, whose each tuple contains 3 eigenvectors of the matrix of
 *          corresponding tuple of \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() != 6.
 */
DataArrayDouble *DataArrayDouble::eigenVectors() const
{
  checkAllocated();
  std::size_t nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6)
    throw INTERP_KERNEL::Exception("DataArrayDouble::eigenVectors : must be an array with exactly 6 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  mcIdType nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,9);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(mcIdType i=0;i<nbOfTuple;i++,src+=6)
    {
      double tmp[3];
      INTERP_KERNEL::computeEigenValues6(src,tmp);
      for(mcIdType j=0;j<3;j++,dest+=3)
        INTERP_KERNEL::computeEigenVectorForEigenValue6(src,tmp[j],1e-12,dest);
    }
  return ret;
}

/*!
 * Computes the inverse matrix of every matrix defined by the tuple of \a this
 * array, which contains either 4, 6 or 9 components. The case of 6 components
 * corresponds to that of the upper triangular matrix.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of components as \a this one, whose each tuple is the inverse
 *          matrix of the matrix of corresponding tuple of \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() is not in [4,6,9].
 */
DataArrayDouble *DataArrayDouble::inverse() const
{
  checkAllocated();
  std::size_t nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6 && nbOfComp!=9 && nbOfComp!=4)
    throw INTERP_KERNEL::Exception("DataArrayDouble::inversion : must be an array with 4,6 or 9 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  mcIdType nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,nbOfComp);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  if(nbOfComp==6)
    for(mcIdType i=0;i<nbOfTuple;i++,dest+=6,src+=6)
      {
        double det=src[0]*src[1]*src[2]+2.*src[4]*src[5]*src[3]-src[0]*src[4]*src[4]-src[2]*src[3]*src[3]-src[1]*src[5]*src[5];
        dest[0]=(src[1]*src[2]-src[4]*src[4])/det;
        dest[1]=(src[0]*src[2]-src[5]*src[5])/det;
        dest[2]=(src[0]*src[1]-src[3]*src[3])/det;
        dest[3]=(src[5]*src[4]-src[3]*src[2])/det;
        dest[4]=(src[5]*src[3]-src[0]*src[4])/det;
        dest[5]=(src[3]*src[4]-src[1]*src[5])/det;
      }
  else if(nbOfComp==4)
    for(mcIdType i=0;i<nbOfTuple;i++,dest+=4,src+=4)
      {
        double det=src[0]*src[3]-src[1]*src[2];
        dest[0]=src[3]/det;
        dest[1]=-src[1]/det;
        dest[2]=-src[2]/det;
        dest[3]=src[0]/det;
      }
  else
    for(mcIdType i=0;i<nbOfTuple;i++,dest+=9,src+=9)
      {
        double det=src[0]*src[4]*src[8]+src[1]*src[5]*src[6]+src[2]*src[3]*src[7]-src[0]*src[5]*src[7]-src[1]*src[3]*src[8]-src[2]*src[4]*src[6];
        dest[0]=(src[4]*src[8]-src[7]*src[5])/det;
        dest[1]=(src[7]*src[2]-src[1]*src[8])/det;
        dest[2]=(src[1]*src[5]-src[4]*src[2])/det;
        dest[3]=(src[6]*src[5]-src[3]*src[8])/det;
        dest[4]=(src[0]*src[8]-src[6]*src[2])/det;
        dest[5]=(src[2]*src[3]-src[0]*src[5])/det;
        dest[6]=(src[3]*src[7]-src[6]*src[4])/det;
        dest[7]=(src[6]*src[1]-src[0]*src[7])/det;
        dest[8]=(src[0]*src[4]-src[1]*src[3])/det;
      }
  return ret;
}

/*!
 * Computes the trace of every matrix defined by the tuple of \a this
 * array, which contains either 4, 6 or 9 components. The case of 6 components
 * corresponds to that of the upper triangular matrix.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing
 *          1 component, whose each tuple is the trace of
 *          the matrix of corresponding tuple of \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() is not in [4,6,9].
 */
DataArrayDouble *DataArrayDouble::trace() const
{
  checkAllocated();
  std::size_t nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6 && nbOfComp!=9 && nbOfComp!=4)
    throw INTERP_KERNEL::Exception("DataArrayDouble::trace : must be an array with 4,6 or 9 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  mcIdType nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  if(nbOfComp==6)
    for(mcIdType i=0;i<nbOfTuple;i++,dest++,src+=6)
      *dest=src[0]+src[1]+src[2];
  else if(nbOfComp==4)
    for(mcIdType i=0;i<nbOfTuple;i++,dest++,src+=4)
      *dest=src[0]+src[3];
  else
    for(mcIdType i=0;i<nbOfTuple;i++,dest++,src+=9)
      *dest=src[0]+src[4]+src[8];
  return ret;
}

/*!
 * Computes the stress deviator tensor of every stress tensor defined by the tuple of
 * \a this array, which contains 6 components.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of components and tuples as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() != 6.
 */
DataArrayDouble *DataArrayDouble::deviator() const
{
  checkAllocated();
  std::size_t nbOfComp=getNumberOfComponents();
  if(nbOfComp!=6)
    throw INTERP_KERNEL::Exception("DataArrayDouble::deviator : must be an array with exactly 6 components !");
  DataArrayDouble *ret=DataArrayDouble::New();
  mcIdType nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,6);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(mcIdType i=0;i<nbOfTuple;i++,dest+=6,src+=6)
    {
      double tr=(src[0]+src[1]+src[2])/3.;
      dest[0]=src[0]-tr;
      dest[1]=src[1]-tr;
      dest[2]=src[2]-tr;
      dest[3]=src[3];
      dest[4]=src[4];
      dest[5]=src[5];
    }
  return ret;
}

/*!
 * Computes the magnitude of every vector defined by the tuple of
 * \a this array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array and one component.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 */
DataArrayDouble *DataArrayDouble::magnitude() const
{
  checkAllocated();
  std::size_t nbOfComp=getNumberOfComponents();
  DataArrayDouble *ret=DataArrayDouble::New();
  mcIdType nbOfTuple=getNumberOfTuples();
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(mcIdType i=0;i<nbOfTuple;i++,dest++)
    {
      double sum=0.;
      for(std::size_t j=0;j<nbOfComp;j++,src++)
        sum+=(*src)*(*src);
      *dest=sqrt(sum);
    }
  return ret;
}

/*!
 * Computes the maximal value within every tuple of \a this array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array and one component.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \sa DataArrayDouble::maxPerTupleWithCompoId
 */
DataArrayDouble *DataArrayDouble::maxPerTuple() const
{
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  mcIdType nbOfTuple(getNumberOfTuples());
  ret->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret->getPointer();
  for(mcIdType i=0;i<nbOfTuple;i++,dest++,src+=nbOfComp)
    *dest=*std::max_element(src,src+nbOfComp);
  return ret.retn();
}

/*!
 * Computes the maximal value within every tuple of \a this array and it returns the first component
 * id for each tuple that corresponds to the maximal value within the tuple.
 *
 *  \param [out] compoIdOfMaxPerTuple - the new new instance of DataArrayInt containing the
 *          same number of tuples and only one component.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array and one component.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \sa DataArrayDouble::maxPerTuple
 */
DataArrayDouble *DataArrayDouble::maxPerTupleWithCompoId(DataArrayIdType* &compoIdOfMaxPerTuple) const
{
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  MCAuto<DataArrayDouble> ret0=DataArrayDouble::New();
  MCAuto<DataArrayIdType> ret1=DataArrayIdType::New();
  mcIdType nbOfTuple=getNumberOfTuples();
  ret0->alloc(nbOfTuple,1); ret1->alloc(nbOfTuple,1);
  const double *src=getConstPointer();
  double *dest=ret0->getPointer(); mcIdType *dest1=ret1->getPointer();
  for(mcIdType i=0;i<nbOfTuple;i++,dest++,dest1++,src+=nbOfComp)
    {
      const double *loc=std::max_element(src,src+nbOfComp);
      *dest=*loc;
      *dest1=ToIdType(std::distance(src,loc));
    }
  compoIdOfMaxPerTuple=ret1.retn();
  return ret0.retn();
}

/*!
 * This method returns a newly allocated DataArrayDouble instance having one component and \c this->getNumberOfTuples() * \c this->getNumberOfTuples() tuples.
 * \n This returned array contains the euclidian distance for each tuple in \a this.
 * \n So the returned array can be seen as a dense symmetrical matrix whose diagonal elements are equal to 0.
 * \n The returned array has only one component (and **not** \c this->getNumberOfTuples() components to avoid the useless memory consumption due to components info in returned DataArrayDouble)
 *
 * \warning use this method with care because it can leads to big amount of consumed memory !
 *
 * \return A newly allocated (huge) MEDCoupling::DataArrayDouble instance that the caller should deal with.
 *
 * \throw If \a this is not allocated.
 *
 * \sa DataArrayDouble::buildEuclidianDistanceDenseMatrixWith
 */
DataArrayDouble *DataArrayDouble::buildEuclidianDistanceDenseMatrix() const
{
  checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  mcIdType nbOfTuples(getNumberOfTuples());
  const double *inData=getConstPointer();
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(nbOfTuples*nbOfTuples,1);
  double *outData=ret->getPointer();
  for(mcIdType i=0;i<nbOfTuples;i++)
    {
      outData[i*nbOfTuples+i]=0.;
      for(mcIdType j=i+1;j<nbOfTuples;j++)
        {
          double dist=0.;
          for(std::size_t k=0;k<nbOfComp;k++)
            { double delta=inData[i*nbOfComp+k]-inData[j*nbOfComp+k]; dist+=delta*delta; }
          dist=sqrt(dist);
          outData[i*nbOfTuples+j]=dist;
          outData[j*nbOfTuples+i]=dist;
        }
    }
  return ret.retn();
}

/*!
 * This method returns a newly allocated DataArrayDouble instance having one component and \c this->getNumberOfTuples() * \c other->getNumberOfTuples() tuples.
 * \n This returned array contains the euclidian distance for each tuple in \a other with each tuple in \a this.
 * \n So the returned array can be seen as a dense rectangular matrix with \c other->getNumberOfTuples() rows and \c this->getNumberOfTuples() columns.
 * \n Output rectangular matrix is sorted along rows.
 * \n The returned array has only one component (and **not** \c this->getNumberOfTuples() components to avoid the useless memory consumption due to components info in returned DataArrayDouble)
 *
 * \warning use this method with care because it can leads to big amount of consumed memory !
 *
 * \param [in] other DataArrayDouble instance having same number of components than \a this.
 * \return A newly allocated (huge) MEDCoupling::DataArrayDouble instance that the caller should deal with.
 *
 * \throw If \a this is not allocated, or if \a other is null or if \a other is not allocated, or if number of components of \a other and \a this differs.
 *
 * \sa DataArrayDouble::buildEuclidianDistanceDenseMatrix
 */
DataArrayDouble *DataArrayDouble::buildEuclidianDistanceDenseMatrixWith(const DataArrayDouble *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::buildEuclidianDistanceDenseMatrixWith : input parameter is null !");
  checkAllocated();
  other->checkAllocated();
  std::size_t nbOfComp(getNumberOfComponents());
  std::size_t otherNbOfComp(other->getNumberOfComponents());
  if(nbOfComp!=otherNbOfComp)
    {
      std::ostringstream oss; oss << "DataArrayDouble::buildEuclidianDistanceDenseMatrixWith : this nb of compo=" << nbOfComp << " and other nb of compo=" << otherNbOfComp << ". It should match !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  mcIdType nbOfTuples(getNumberOfTuples());
  mcIdType otherNbOfTuples(other->getNumberOfTuples());
  const double *inData=getConstPointer();
  const double *inDataOther=other->getConstPointer();
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(otherNbOfTuples*nbOfTuples,1);
  double *outData=ret->getPointer();
  for(mcIdType i=0;i<otherNbOfTuples;i++,inDataOther+=nbOfComp)
    {
      for(mcIdType j=0;j<nbOfTuples;j++)
        {
          double dist=0.;
          for(std::size_t k=0;k<nbOfComp;k++)
            { double delta=inDataOther[k]-inData[j*nbOfComp+k]; dist+=delta*delta; }
          dist=sqrt(dist);
          outData[i*nbOfTuples+j]=dist;
        }
    }
  return ret.retn();
}

/*!
 * This method expects that \a this stores 3 tuples containing 2 components each.
 * Each of this tuples represent a point into 2D space.
 * This method tries to find an arc of circle starting from first point (tuple) to 2nd and middle point (tuple) along 3nd and last point (tuple).
 * If such arc of circle exists, the corresponding center, radius of circle is returned. And additionnaly the length of arc expressed with an \a ang output variable in ]0,2*pi[.
 *
 *  \throw If \a this is not allocated.
 *  \throw If \a this has not 3 tuples of 2 components
 *  \throw If tuples/points in \a this are aligned
 */
void DataArrayDouble::asArcOfCircle(double center[2], double& radius, double& ang) const
{
  checkAllocated();
  INTERP_KERNEL::QuadraticPlanarPrecision arcPrec(1e-14);
  if(getNumberOfTuples()!=3 && getNumberOfComponents()!=2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::asArcCircle : this method expects");
  const double *pt(begin());
  MCAuto<INTERP_KERNEL::Node> n0(new INTERP_KERNEL::Node(pt[0],pt[1])),n1(new INTERP_KERNEL::Node(pt[2],pt[3])),n2(new INTERP_KERNEL::Node(pt[4],pt[5]));
  {
    INTERP_KERNEL::AutoCppPtr<INTERP_KERNEL::EdgeLin> e1(new INTERP_KERNEL::EdgeLin(n0,n2)),e2(new INTERP_KERNEL::EdgeLin(n2,n1));
    INTERP_KERNEL::SegSegIntersector inters(*e1,*e2);
    bool colinearity(inters.areColinears());
    if(colinearity)
      throw INTERP_KERNEL::Exception("DataArrayDouble::asArcOfCircle : 3 points in this have been detected as colinear !");
  }
  INTERP_KERNEL::AutoCppPtr<INTERP_KERNEL::EdgeArcCircle> ret(new INTERP_KERNEL::EdgeArcCircle(n0,n2,n1));
  const double *c(ret->getCenter());
  center[0]=c[0]; center[1]=c[1];
  radius=ret->getRadius();
  ang=ret->getAngle();
}

/*!
 * Sorts value within every tuple of \a this array.
 *  \param [in] asc - if \a true, the values are sorted in ascending order, else,
 *              in descending order.
 *  \throw If \a this is not allocated.
 */
void DataArrayDouble::sortPerTuple(bool asc)
{
  checkAllocated();
  double *pt=getPointer();
  mcIdType nbOfTuple(getNumberOfTuples());
  std::size_t nbOfComp(getNumberOfComponents());
  if(asc)
    for(mcIdType i=0;i<nbOfTuple;i++,pt+=nbOfComp)
      std::sort(pt,pt+nbOfComp);
  else
    for(mcIdType i=0;i<nbOfTuple;i++,pt+=nbOfComp)
      std::sort(pt,pt+nbOfComp,std::greater<double>());
  declareAsNew();
}

/*!
 * Modify all elements of \a this array, so that
 * an element _x_ becomes \f$ numerator / x \f$.
 *  \warning If an exception is thrown because of presence of 0.0 element in \a this
 *           array, all elements processed before detection of the zero element remain
 *           modified.
 *  \param [in] numerator - the numerator used to modify array elements.
 *  \throw If \a this is not allocated.
 *  \throw If there is an element equal to 0.0 in \a this array.
 */
void DataArrayDouble::applyInv(double numerator)
{
  checkAllocated();
  double *ptr=getPointer();
  std::size_t nbOfElems=getNbOfElems();
  for(std::size_t i=0;i<nbOfElems;i++,ptr++)
    {
      if(std::abs(*ptr)>std::numeric_limits<double>::min())
        {
          *ptr=numerator/(*ptr);
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::applyInv : presence of null value in tuple #" << i/getNumberOfComponents() << " component #" << i%getNumberOfComponents();
          oss << " !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  declareAsNew();
}

/*!
 * Modify all elements of \a this array, so that
 * an element _x_ becomes <em> val ^ x </em>. Contrary to DataArrayInt::applyPow
 * all values in \a this have to be >= 0 if val is \b not integer.
 *  \param [in] val - the value used to apply pow on all array elements.
 *  \throw If \a this is not allocated.
 *  \warning If an exception is thrown because of presence of 0 element in \a this
 *           array and \a val is \b not integer, all elements processed before detection of the zero element remain
 *           modified.
 */
void DataArrayDouble::applyPow(double val)
{
  checkAllocated();
  double *ptr=getPointer();
  std::size_t nbOfElems=getNbOfElems();
  int val2=(int)val;
  bool isInt=((double)val2)==val;
  if(!isInt)
    {
      for(std::size_t i=0;i<nbOfElems;i++,ptr++)
        {
          if(*ptr>=0)
            *ptr=pow(*ptr,val);
          else
            {
              std::ostringstream oss; oss << "DataArrayDouble::applyPow (double) : At elem # " << i << " value is " << *ptr << " ! must be >=0. !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
  else
    {
      for(std::size_t i=0;i<nbOfElems;i++,ptr++)
        *ptr=pow(*ptr,val2);
    }
  declareAsNew();
}

/*!
 * Modify all elements of \a this array, so that
 * an element _x_ becomes \f$ val ^ x \f$.
 *  \param [in] val - the value used to apply pow on all array elements.
 *  \throw If \a this is not allocated.
 *  \throw If \a val < 0.
 *  \warning If an exception is thrown because of presence of 0 element in \a this
 *           array, all elements processed before detection of the zero element remain
 *           modified.
 */
void DataArrayDouble::applyRPow(double val)
{
  checkAllocated();
  if(val<0.)
    throw INTERP_KERNEL::Exception("DataArrayDouble::applyRPow : the input value has to be >= 0 !");
  double *ptr=getPointer();
  std::size_t nbOfElems=getNbOfElems();
  for(std::size_t i=0;i<nbOfElems;i++,ptr++)
    *ptr=pow(val,*ptr);
  declareAsNew();
}

/*!
 * Returns a new DataArrayDouble created from \a this one by applying \a
 * FunctionToEvaluate to every tuple of \a this array. Textual data is not copied.
 * For more info see \ref MEDCouplingArrayApplyFunc
 *  \param [in] nbOfComp - number of components in the result array.
 *  \param [in] func - the \a FunctionToEvaluate declared as
 *              \c bool (*\a func)(\c const \c double *\a pos, \c double *\a res),
 *              where \a pos points to the first component of a tuple of \a this array
 *              and \a res points to the first component of a tuple of the result array.
 *              Note that length (number of components) of \a pos can differ from
 *              that of \a res.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a func returns \a false.
 */
DataArrayDouble *DataArrayDouble::applyFunc(std::size_t nbOfComp, FunctionToEvaluate func) const
{
  checkAllocated();
  DataArrayDouble *newArr=DataArrayDouble::New();
  mcIdType nbOfTuples(getNumberOfTuples());
  std::size_t oldNbOfComp(getNumberOfComponents());
  newArr->alloc(nbOfTuples,nbOfComp);
  const double *ptr=getConstPointer();
  double *ptrToFill=newArr->getPointer();
  for(mcIdType i=0;i<nbOfTuples;i++)
    {
      if(!func(ptr+i*oldNbOfComp,ptrToFill+i*nbOfComp))
        {
          std::ostringstream oss; oss << "For tuple # " << i << " with value (";
          std::copy(ptr+oldNbOfComp*i,ptr+oldNbOfComp*(i+1),std::ostream_iterator<double>(oss,", "));
          oss << ") : Evaluation of function failed !";
          newArr->decrRef();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return newArr;
}

/*!
 * Returns a new DataArrayDouble created from \a this one by applying a function to every
 * tuple of \a this array. Textual data is not copied.
 * For more info see \ref MEDCouplingArrayApplyFunc1.
 *  \param [in] nbOfComp - number of components in the result array.
 *  \param [in] func - the expression defining how to transform a tuple of \a this array.
 *              Supported expressions are described \ref MEDCouplingArrayApplyFuncExpr "here".
 *  \param [in] isSafe - By default true. If true invalid operation (division by 0. acos of value > 1. ...) leads to a throw of an exception.
 *              If false the computation is carried on without any notification. When false the evaluation is a little faster.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array and \a nbOfComp components.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \throw If computing \a func fails.
 */
DataArrayDouble *DataArrayDouble::applyFunc(std::size_t nbOfComp, const std::string& func, bool isSafe) const
{
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  std::set<std::string> vars;
  expr.getTrueSetOfVars(vars);
  std::vector<std::string> varsV(vars.begin(),vars.end());
  return applyFuncNamedCompo(nbOfComp,varsV,func,isSafe);
}

/*!
 * Returns a new DataArrayDouble created from \a this one by applying a function to every
 * tuple of \a this array. Textual data is not copied. This method works by tuples (whatever its size).
 * If \a this is a one component array, call applyFuncOnThis instead that performs the same work faster.
 *
 * For more info see \ref MEDCouplingArrayApplyFunc0.
 *  \param [in] func - the expression defining how to transform a tuple of \a this array.
 *              Supported expressions are described \ref MEDCouplingArrayApplyFuncExpr "here".
 *  \param [in] isSafe - By default true. If true invalid operation (division by 0. acos of value > 1. ...) leads to a throw of an exception.
 *                       If false the computation is carried on without any notification. When false the evaluation is a little faster.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples and components as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \sa applyFuncOnThis
 *  \throw If \a this is not allocated.
 *  \throw If computing \a func fails.
 */
DataArrayDouble *DataArrayDouble::applyFunc(const std::string& func, bool isSafe) const
{
  std::size_t nbOfComp(getNumberOfComponents());
  if(nbOfComp<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::applyFunc : output number of component must be > 0 !");
  checkAllocated();
  mcIdType nbOfTuples(getNumberOfTuples());
  MCAuto<DataArrayDouble> newArr(DataArrayDouble::New());
  newArr->alloc(nbOfTuples,nbOfComp);
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  std::set<std::string> vars;
  expr.getTrueSetOfVars(vars);
  if(vars.size()>1)
    {
      std::ostringstream oss; oss << "DataArrayDouble::applyFunc : this method works only with at most one var func expression ! If you need to map comps on variables please use applyFuncCompo or applyFuncNamedCompo instead ! Vars in expr are : ";
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(vars.empty())
    {
      expr.prepareFastEvaluator();
      newArr->rearrange(1);
      newArr->fillWithValue(expr.evaluateDouble());
      newArr->rearrange(nbOfComp);
      return newArr.retn();
    }
  std::vector<std::string> vars2(vars.begin(),vars.end());
  double buff,*ptrToFill(newArr->getPointer());
  const double *ptr(begin());
  std::vector<double> stck;
  expr.prepareExprEvaluationDouble(vars2,1,1,0,&buff,&buff+1);
  expr.prepareFastEvaluator();
  if(!isSafe)
    {
      for(mcIdType i=0;i<nbOfTuples;i++)
        {
          for(std::size_t iComp=0;iComp<nbOfComp;iComp++,ptr++,ptrToFill++)
            {
              buff=*ptr;
              expr.evaluateDoubleInternal(stck);
              *ptrToFill=stck.back();
              stck.pop_back();
            }
        }
    }
  else
    {
      for(mcIdType i=0;i<nbOfTuples;i++)
        {
          for(std::size_t iComp=0;iComp<nbOfComp;iComp++,ptr++,ptrToFill++)
            {
              buff=*ptr;
              try
              {
                  expr.evaluateDoubleInternalSafe(stck);
              }
              catch(INTERP_KERNEL::Exception& e)
              {
                  std::ostringstream oss; oss << "For tuple # " << i << " component # " << iComp << " with value (";
                  oss << buff;
                  oss << ") : Evaluation of function failed !" << e.what();
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
              *ptrToFill=stck.back();
              stck.pop_back();
            }
        }
    }
  return newArr.retn();
}

/*!
 * This method is a non const method that modify the array in \a this.
 * This method only works on one component array. It means that function \a func must
 * contain at most one variable.
 * This method is a specialization of applyFunc method with one parameter on one component array.
 *
 *  \param [in] func - the expression defining how to transform a tuple of \a this array.
 *              Supported expressions are described \ref MEDCouplingArrayApplyFuncExpr "here".
 *  \param [in] isSafe - By default true. If true invalid operation (division by 0. acos of value > 1. ...) leads to a throw of an exception.
 *              If false the computation is carried on without any notification. When false the evaluation is a little faster.
 *
 * \sa applyFunc
 */
void DataArrayDouble::applyFuncOnThis(const std::string& func, bool isSafe)
{
  std::size_t nbOfComp(getNumberOfComponents());
  if(nbOfComp<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::applyFuncOnThis : output number of component must be > 0 !");
  checkAllocated();
  mcIdType nbOfTuples(getNumberOfTuples());
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  std::set<std::string> vars;
  expr.getTrueSetOfVars(vars);
  if(vars.size()>1)
    {
      std::ostringstream oss; oss << "DataArrayDouble::applyFuncOnThis : this method works only with at most one var func expression ! If you need to map comps on variables please use applyFuncCompo or applyFuncNamedCompo instead ! Vars in expr are : ";
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(vars.empty())
    {
      expr.prepareFastEvaluator();
      std::vector<std::string> compInfo(getInfoOnComponents());
      rearrange(1);
      fillWithValue(expr.evaluateDouble());
      rearrange(nbOfComp);
      setInfoOnComponents(compInfo);
      return ;
    }
  std::vector<std::string> vars2(vars.begin(),vars.end());
  double buff,*ptrToFill(getPointer());
  const double *ptr(begin());
  std::vector<double> stck;
  expr.prepareExprEvaluationDouble(vars2,1,1,0,&buff,&buff+1);
  expr.prepareFastEvaluator();
  if(!isSafe)
    {
      for(mcIdType i=0;i<nbOfTuples;i++)
        {
          for(std::size_t iComp=0;iComp<nbOfComp;iComp++,ptr++,ptrToFill++)
            {
              buff=*ptr;
              expr.evaluateDoubleInternal(stck);
              *ptrToFill=stck.back();
              stck.pop_back();
            }
        }
    }
  else
    {
      for(mcIdType i=0;i<nbOfTuples;i++)
        {
          for(std::size_t iComp=0;iComp<nbOfComp;iComp++,ptr++,ptrToFill++)
            {
              buff=*ptr;
              try
              {
                  expr.evaluateDoubleInternalSafe(stck);
              }
              catch(INTERP_KERNEL::Exception& e)
              {
                  std::ostringstream oss; oss << "For tuple # " << i << " component # " << iComp << " with value (";
                  oss << buff;
                  oss << ") : Evaluation of function failed !" << e.what();
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
              *ptrToFill=stck.back();
              stck.pop_back();
            }
        }
    }
}

/*!
 * Returns a new DataArrayDouble created from \a this one by applying a function to every
 * tuple of \a this array. Textual data is not copied.
 * For more info see \ref MEDCouplingArrayApplyFunc2.
 *  \param [in] nbOfComp - number of components in the result array.
 *  \param [in] func - the expression defining how to transform a tuple of \a this array.
 *              Supported expressions are described \ref MEDCouplingArrayApplyFuncExpr "here".
 *  \param [in] isSafe - By default true. If true invalid operation (division by 0. acos of value > 1. ...) leads to a throw of an exception.
 *              If false the computation is carried on without any notification. When false the evaluation is a little faster.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a func contains vars that are not in \a this->getInfoOnComponent().
 *  \throw If computing \a func fails.
 */
DataArrayDouble *DataArrayDouble::applyFuncCompo(std::size_t nbOfComp, const std::string& func, bool isSafe) const
{
  return applyFuncNamedCompo(nbOfComp,getVarsOnComponent(),func,isSafe);
}

/*!
 * Returns a new DataArrayDouble created from \a this one by applying a function to every
 * tuple of \a this array. Textual data is not copied.
 * For more info see \ref MEDCouplingArrayApplyFunc3.
 *  \param [in] nbOfComp - number of components in the result array.
 *  \param [in] varsOrder - sequence of vars defining their order.
 *  \param [in] func - the expression defining how to transform a tuple of \a this array.
 *              Supported expressions are described \ref MEDCouplingArrayApplyFuncExpr "here".
 *  \param [in] isSafe - By default true. If true invalid operation (division by 0. acos of value > 1. ...) leads to a throw of an exception.
 *              If false the computation is carried on without any notification. When false the evaluation is a little faster.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
 *          same number of tuples as \a this array.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a func contains vars not in \a varsOrder.
 *  \throw If computing \a func fails.
 */
DataArrayDouble *DataArrayDouble::applyFuncNamedCompo(std::size_t nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func, bool isSafe) const
{
  if(nbOfComp<=0)
    throw INTERP_KERNEL::Exception("DataArrayDouble::applyFuncNamedCompo : output number of component must be > 0 !");
  std::vector<std::string> varsOrder2(varsOrder);
  std::size_t oldNbOfComp(getNumberOfComponents());
  for(std::size_t i=varsOrder.size();i<oldNbOfComp;i++)
    varsOrder2.push_back(std::string());
  checkAllocated();
  mcIdType nbOfTuples(getNumberOfTuples());
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  std::set<std::string> vars;
  expr.getTrueSetOfVars(vars);
  if(vars.size()>oldNbOfComp)
    {
      std::ostringstream oss; oss << "The field has " << oldNbOfComp << " components and there are ";
      oss << vars.size() << " variables : ";
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  MCAuto<DataArrayDouble> newArr(DataArrayDouble::New());
  newArr->alloc(nbOfTuples,nbOfComp);
  INTERP_KERNEL::AutoPtr<double> buff(new double[oldNbOfComp]);
  double *buffPtr(buff),*ptrToFill;
  std::vector<double> stck;
  for(std::size_t iComp=0;iComp<nbOfComp;iComp++)
    {
      expr.prepareExprEvaluationDouble(varsOrder2,(int)oldNbOfComp,(int)nbOfComp,(int)iComp,buffPtr,buffPtr+oldNbOfComp);
      expr.prepareFastEvaluator();
      const double *ptr(getConstPointer());
      ptrToFill=newArr->getPointer()+iComp;
      if(!isSafe)
        {
          for(mcIdType i=0;i<nbOfTuples;i++,ptrToFill+=nbOfComp,ptr+=oldNbOfComp)
            {
              std::copy(ptr,ptr+oldNbOfComp,buffPtr);
              expr.evaluateDoubleInternal(stck);
              *ptrToFill=stck.back();
              stck.pop_back();
            }
        }
      else
        {
          for(mcIdType i=0;i<nbOfTuples;i++,ptrToFill+=nbOfComp,ptr+=oldNbOfComp)
            {
              std::copy(ptr,ptr+oldNbOfComp,buffPtr);
              try
              {
                  expr.evaluateDoubleInternalSafe(stck);
                  *ptrToFill=stck.back();
                  stck.pop_back();
              }
              catch(INTERP_KERNEL::Exception& e)
              {
                  std::ostringstream oss; oss << "For tuple # " << i << " with value (";
                  std::copy(ptr+oldNbOfComp*i,ptr+oldNbOfComp*(i+1),std::ostream_iterator<double>(oss,", "));
                  oss << ") : Evaluation of function failed !" << e.what();
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            }
        }
    }
  return newArr.retn();
}

void DataArrayDouble::applyFuncFast32(const std::string& func)
{
  checkAllocated();
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  char *funcStr=expr.compileX86();
  MYFUNCPTR funcPtr;
  *((void **)&funcPtr)=funcStr;//he he...
  //
  double *ptr=getPointer();
  std::size_t nbOfComp=getNumberOfComponents();
  mcIdType nbOfTuples=getNumberOfTuples();
  std::size_t nbOfElems=nbOfTuples*nbOfComp;
  for(std::size_t i=0;i<nbOfElems;i++,ptr++)
    *ptr=funcPtr(*ptr);
  declareAsNew();
}

void DataArrayDouble::applyFuncFast64(const std::string& func)
{
  checkAllocated();
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  char *funcStr=expr.compileX86_64();
  MYFUNCPTR funcPtr;
  *((void **)&funcPtr)=funcStr;//he he...
  //
  double *ptr=getPointer();
  std::size_t nbOfComp=getNumberOfComponents();
  mcIdType nbOfTuples=getNumberOfTuples();
  std::size_t nbOfElems=nbOfTuples*nbOfComp;
  for(std::size_t i=0;i<nbOfElems;i++,ptr++)
    *ptr=funcPtr(*ptr);
  declareAsNew();
}

/*!
 * \return a new object that is the result of the symmetry along 3D plane defined by its normal vector \a normalVector and a point \a point.
 */
MCAuto<DataArrayDouble> DataArrayDouble::symmetry3DPlane(const double point[3], const double normalVector[3]) const
{
  checkAllocated();
  if(getNumberOfComponents()!=3)
    throw INTERP_KERNEL::Exception("DataArrayDouble::symmetry3DPlane : this is excepted to have 3 components !");
  mcIdType nbTuples(getNumberOfTuples());
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbTuples,3);
  Symmetry3DPlane(point,normalVector,nbTuples,begin(),ret->getPointer());
  return ret;
}

DataArrayDoubleIterator *DataArrayDouble::iterator()
{
  return new DataArrayDoubleIterator(this);
}

/*!
 * Returns a new DataArrayInt containing indices of tuples of \a this one-dimensional
 * array whose values are within a given range. Textual data is not copied.
 *  \param [in] vmin - a lowest acceptable value (included).
 *  \param [in] vmax - a greatest acceptable value (included).
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *
 *  \sa DataArrayDouble::findIdsNotInRange
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcdataarraydouble_getidsinrange "Here is a C++ example".<br>
 *  \ref py_mcdataarraydouble_getidsinrange "Here is a Python example".
 *  \endif
 */
DataArrayIdType *DataArrayDouble::findIdsInRange(double vmin, double vmax) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::findIdsInRange : this must have exactly one component !");
  const double *cptr(begin());
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);
  mcIdType nbOfTuples(getNumberOfTuples());
  for(mcIdType i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr>=vmin && *cptr<=vmax)
      ret->pushBackSilent(i);
  return ret.retn();
}

/*!
 * Returns a new DataArrayInt containing indices of tuples of \a this one-dimensional
 * array whose values are not within a given range. Textual data is not copied.
 *  \param [in] vmin - a lowest not acceptable value (excluded).
 *  \param [in] vmax - a greatest not acceptable value (excluded).
 *  \return DataArrayInt * - the new instance of DataArrayInt.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *
 *  \sa DataArrayDouble::findIdsInRange
 */
DataArrayIdType *DataArrayDouble::findIdsNotInRange(double vmin, double vmax) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::findIdsNotInRange : this must have exactly one component !");
  const double *cptr(begin());
  MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);
  mcIdType nbOfTuples(getNumberOfTuples());
  for(mcIdType i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr<vmin || *cptr>vmax)
      ret->pushBackSilent(i);
  return ret.retn();
}

/*!
 * Returns a new DataArrayDouble by concatenating two given arrays, so that (1) the number
 * of tuples in the result array is a sum of the number of tuples of given arrays and (2)
 * the number of component in the result array is same as that of each of given arrays.
 * Info on components is copied from the first of the given arrays. Number of components
 * in the given arrays must be  the same.
 *  \param [in] a1 - an array to include in the result array.
 *  \param [in] a2 - another array to include in the result array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If both \a a1 and \a a2 are NULL.
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents().
 */
DataArrayDouble *DataArrayDouble::Aggregate(const DataArrayDouble *a1, const DataArrayDouble *a2)
{
  std::vector<const DataArrayDouble *> tmp(2);
  tmp[0]=a1; tmp[1]=a2;
  return Aggregate(tmp);
}

/*!
 * Returns a new DataArrayDouble by concatenating all given arrays, so that (1) the number
 * of tuples in the result array is a sum of the number of tuples of given arrays and (2)
 * the number of component in the result array is same as that of each of given arrays.
 * Info on components is copied from the first of the given arrays. Number of components
 * in the given arrays must be  the same.
 * If the number of non null of elements in \a arr is equal to one the returned object is a copy of it
 * not the object itself.
 *  \param [in] arr - a sequence of arrays to include in the result array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If all arrays within \a arr are NULL.
 *  \throw If getNumberOfComponents() of arrays within \a arr.
 */
DataArrayDouble *DataArrayDouble::Aggregate(const std::vector<const DataArrayDouble *>& arr)
{
  std::vector<const DataArrayDouble *> a;
  for(std::vector<const DataArrayDouble *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
    if(*it4)
      a.push_back(*it4);
  if(a.empty())
    throw INTERP_KERNEL::Exception("DataArrayDouble::Aggregate : input list must contain at least one NON EMPTY DataArrayDouble !");
  std::vector<const DataArrayDouble *>::const_iterator it=a.begin();
  std::size_t nbOfComp((*it)->getNumberOfComponents());
  mcIdType nbt=(*it++)->getNumberOfTuples();
  for(mcIdType i=1;it!=a.end();it++,i++)
    {
      if((*it)->getNumberOfComponents()!=nbOfComp)
        throw INTERP_KERNEL::Exception("DataArrayDouble::Aggregate : Nb of components mismatch for array aggregation !");
      nbt+=(*it)->getNumberOfTuples();
    }
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(nbt,nbOfComp);
  double *pt=ret->getPointer();
  for(it=a.begin();it!=a.end();it++)
    pt=std::copy((*it)->getConstPointer(),(*it)->getConstPointer()+(*it)->getNbOfElems(),pt);
  ret->copyStringInfoFrom(*(a[0]));
  return ret.retn();
}

/*!
 * Returns a new DataArrayDouble containing a dot product of two given arrays, so that
 * the i-th tuple of the result array is a sum of products of j-th components of i-th
 * tuples of given arrays (\f$ a_i = \sum_{j=1}^n a1_j * a2_j \f$).
 * Info on components and name is copied from the first of the given arrays.
 * Number of tuples and components in the given arrays must be the same.
 *  \param [in] a1 - a given array.
 *  \param [in] a2 - another given array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If any given array is not allocated.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents()
 */
DataArrayDouble *DataArrayDouble::Dot(const DataArrayDouble *a1, const DataArrayDouble *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Dot : input DataArrayDouble instance is NULL !");
  a1->checkAllocated();
  a2->checkAllocated();
  std::size_t nbOfComp(a1->getNumberOfComponents());
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array Dot !");
  mcIdType nbOfTuple(a1->getNumberOfTuples());
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array Dot !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,1);
  double *retPtr=ret->getPointer();
  const double *a1Ptr=a1->begin(),*a2Ptr(a2->begin());
  for(mcIdType i=0;i<nbOfTuple;i++)
    {
      double sum=0.;
      for(std::size_t j=0;j<nbOfComp;j++)
        sum+=a1Ptr[i*nbOfComp+j]*a2Ptr[i*nbOfComp+j];
      retPtr[i]=sum;
    }
  ret->setInfoOnComponent(0,a1->getInfoOnComponent(0));
  ret->setName(a1->getName());
  return ret;
}

/*!
 * Returns a new DataArrayDouble containing a cross product of two given arrays, so that
 * the i-th tuple of the result array contains 3 components of a vector which is a cross
 * product of two vectors defined by the i-th tuples of given arrays.
 * Info on components is copied from the first of the given arrays.
 * Number of tuples in the given arrays must be the same.
 * Number of components in the given arrays must be 3.
 *  \param [in] a1 - a given array.
 *  \param [in] a2 - another given array.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 *  \throw If \a a1->getNumberOfComponents() != 3
 *  \throw If \a a2->getNumberOfComponents() != 3
 */
DataArrayDouble *DataArrayDouble::CrossProduct(const DataArrayDouble *a1, const DataArrayDouble *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::CrossProduct : input DataArrayDouble instance is NULL !");
  std::size_t nbOfComp(a1->getNumberOfComponents());
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array crossProduct !");
  if(nbOfComp!=3)
    throw INTERP_KERNEL::Exception("Nb of components must be equal to 3 for array crossProduct !");
  mcIdType nbOfTuple(a1->getNumberOfTuples());
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array crossProduct !");
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(nbOfTuple,3);
  double *retPtr=ret->getPointer();
  const double *a1Ptr(a1->begin()),*a2Ptr(a2->begin());
  for(mcIdType i=0;i<nbOfTuple;i++)
    {
      retPtr[3*i]=a1Ptr[3*i+1]*a2Ptr[3*i+2]-a1Ptr[3*i+2]*a2Ptr[3*i+1];
      retPtr[3*i+1]=a1Ptr[3*i+2]*a2Ptr[3*i]-a1Ptr[3*i]*a2Ptr[3*i+2];
      retPtr[3*i+2]=a1Ptr[3*i]*a2Ptr[3*i+1]-a1Ptr[3*i+1]*a2Ptr[3*i];
    }
  ret->copyStringInfoFrom(*a1);
  return ret;
}

/*!
 * Returns a new DataArrayDouble containing maximal values of two given arrays.
 * Info on components is copied from the first of the given arrays.
 * Number of tuples and components in the given arrays must be the same.
 *  \param [in] a1 - an array to compare values with another one.
 *  \param [in] a2 - another array to compare values with the first one.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents()
 */
DataArrayDouble *DataArrayDouble::Max(const DataArrayDouble *a1, const DataArrayDouble *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Max : input DataArrayDouble instance is NULL !");
  std::size_t nbOfComp(a1->getNumberOfComponents());
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array Max !");
  mcIdType nbOfTuple(a1->getNumberOfTuples());
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array Max !");
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbOfTuple,nbOfComp);
  double *retPtr(ret->getPointer());
  const double *a1Ptr(a1->begin()),*a2Ptr(a2->begin());
  std::size_t nbElem(nbOfTuple*nbOfComp);
  for(std::size_t i=0;i<nbElem;i++)
    retPtr[i]=std::max(a1Ptr[i],a2Ptr[i]);
  ret->copyStringInfoFrom(*a1);
  return ret.retn();
}

/*!
 * Returns a new DataArrayDouble containing minimal values of two given arrays.
 * Info on components is copied from the first of the given arrays.
 * Number of tuples and components in the given arrays must be the same.
 *  \param [in] a1 - an array to compare values with another one.
 *  \param [in] a2 - another array to compare values with the first one.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents()
 */
DataArrayDouble *DataArrayDouble::Min(const DataArrayDouble *a1, const DataArrayDouble *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Min : input DataArrayDouble instance is NULL !");
  std::size_t nbOfComp(a1->getNumberOfComponents());
  if(nbOfComp!=a2->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("Nb of components mismatch for array min !");
  mcIdType nbOfTuple(a1->getNumberOfTuples());
  if(nbOfTuple!=a2->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array min !");
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbOfTuple,nbOfComp);
  double *retPtr(ret->getPointer());
  const double *a1Ptr(a1->begin()),*a2Ptr(a2->begin());
  std::size_t nbElem(nbOfTuple*nbOfComp);
  for(std::size_t i=0;i<nbElem;i++)
    retPtr[i]=std::min(a1Ptr[i],a2Ptr[i]);
  ret->copyStringInfoFrom(*a1);
  return ret.retn();
}

/*!
 * Returns a new DataArrayDouble that is the result of pow of two given arrays. There are 3
 * valid cases.
 *
 *  \param [in] a1 - an array to pow up.
 *  \param [in] a2 - another array to sum up.
 *  \return DataArrayDouble * - the new instance of DataArrayDouble.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 *  \throw If \a a1->getNumberOfComponents() != 1 or \a a2->getNumberOfComponents() != 1.
 *  \throw If there is a negative value in \a a1.
 */
DataArrayDouble *DataArrayDouble::Pow(const DataArrayDouble *a1, const DataArrayDouble *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Pow : at least one of input instances is null !");
  mcIdType nbOfTuple=a1->getNumberOfTuples();
  mcIdType nbOfTuple2=a2->getNumberOfTuples();
  std::size_t nbOfComp=a1->getNumberOfComponents();
  std::size_t nbOfComp2=a2->getNumberOfComponents();
  if(nbOfTuple!=nbOfTuple2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Pow : number of tuples mismatches !");
  if(nbOfComp!=1 || nbOfComp2!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Pow : number of components of both arrays must be equal to 1 !");
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New(); ret->alloc(nbOfTuple,1);
  const double *ptr1(a1->begin()),*ptr2(a2->begin());
  double *ptr=ret->getPointer();
  for(mcIdType i=0;i<nbOfTuple;i++,ptr1++,ptr2++,ptr++)
    {
      if(*ptr1>=0)
        {
          *ptr=pow(*ptr1,*ptr2);
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::Pow : on tuple #" << i << " of a1 value is < 0 (" << *ptr1 << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret.retn();
}

/*!
 * Apply pow on values of another DataArrayDouble to values of \a this one.
 *
 *  \param [in] other - an array to pow to \a this one.
 *  \throw If \a other is NULL.
 *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples()
 *  \throw If \a this->getNumberOfComponents() != 1 or \a other->getNumberOfComponents() != 1
 *  \throw If there is a negative value in \a this.
 */
void DataArrayDouble::powEqual(const DataArrayDouble *other)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayDouble::powEqual : input instance is null !");
  mcIdType nbOfTuple=getNumberOfTuples();
  mcIdType nbOfTuple2=other->getNumberOfTuples();
  std::size_t nbOfComp=getNumberOfComponents();
  std::size_t nbOfComp2=other->getNumberOfComponents();
  if(nbOfTuple!=nbOfTuple2)
    throw INTERP_KERNEL::Exception("DataArrayDouble::powEqual : number of tuples mismatches !");
  if(nbOfComp!=1 || nbOfComp2!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::powEqual : number of components of both arrays must be equal to 1 !");
  double *ptr=getPointer();
  const double *ptrc=other->begin();
  for(mcIdType i=0;i<nbOfTuple;i++,ptrc++,ptr++)
    {
      if(*ptr>=0)
        *ptr=pow(*ptr,*ptrc);
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::powEqual : on tuple #" << i << " of this value is < 0 (" << *ptr << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  declareAsNew();
}

/*!
 * This method is \b NOT wrapped into python because it can be useful only for performance reasons in C++ context.
 * All values in \a this must be 0. or 1. within eps error. 0 means false, 1 means true.
 * If an another value than 0 or 1 appear (within eps precision) an INTERP_KERNEL::Exception will be thrown.
 *
 * \throw if \a this is not allocated.
 * \throw if \a this has not exactly one component.
 */
std::vector<bool> DataArrayDouble::toVectorOfBool(double eps) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayDouble::toVectorOfBool : must be applied on single component array !");
  mcIdType nbt(getNumberOfTuples());
  std::vector<bool> ret(nbt);
  const double *pt(begin());
  for(mcIdType i=0;i<nbt;i++)
    {
      if(fabs(pt[i])<eps)
        ret[i]=false;
      else if(fabs(pt[i]-1.)<eps)
        ret[i]=true;
      else
        {
          std::ostringstream oss; oss << "DataArrayDouble::toVectorOfBool : the tuple #" << i << " has value " << pt[i] << " is invalid ! must be 0. or 1. !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret;
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * Server side.
 */
void DataArrayDouble::getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const
{
  tinyInfo.resize(2);
  if(isAllocated())
    {
      tinyInfo[0]=getNumberOfTuples();
      tinyInfo[1]=ToIdType(getNumberOfComponents());
    }
  else
    {
      tinyInfo[0]=-1;
      tinyInfo[1]=-1;
    }
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * Server side.
 */
void DataArrayDouble::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
{
  if(isAllocated())
    {
      std::size_t nbOfCompo(getNumberOfComponents());
      tinyInfo.resize(nbOfCompo+1);
      tinyInfo[0]=getName();
      for(std::size_t i=0;i<nbOfCompo;i++)
        tinyInfo[i+1]=getInfoOnComponent(i);
    }
  else
    {
      tinyInfo.resize(1);
      tinyInfo[0]=getName();
    }
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 * This method returns if a feeding is needed.
 */
bool DataArrayDouble::resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI)
{
  mcIdType nbOfTuple=tinyInfoI[0];
  mcIdType nbOfComp=tinyInfoI[1];
  if(nbOfTuple!=-1 || nbOfComp!=-1)
    {
      alloc(nbOfTuple,nbOfComp);
      return true;
    }
  return false;
}

/*!
 * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
 */
void DataArrayDouble::finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<std::string>& tinyInfoS)
{
  setName(tinyInfoS[0]);
  if(isAllocated())
    {
      std::size_t nbOfCompo(getNumberOfComponents());
      for(std::size_t i=0;i<nbOfCompo;i++)
        setInfoOnComponent(i,tinyInfoS[i+1]);
    }
}

/*!
 * Low static method that operates 3D rotation of 'nbNodes' 3D nodes whose coordinates are arranged in \a coordsIn
 * around an axe ( \a center, \a vect) and with angle \a angle.
 */
void DataArrayDouble::Rotate3DAlg(const double *center, const double *vect, double angle, mcIdType nbNodes, const double *coordsIn, double *coordsOut)
{
  if(!center || !vect)
    throw INTERP_KERNEL::Exception("DataArrayDouble::Rotate3DAlg : null vector in input !");
  double sina(sin(angle));
  double cosa(cos(angle));
  double vectorNorm[3];
  double matrix[9];
  double matrixTmp[9];
  double norm(sqrt(vect[0]*vect[0]+vect[1]*vect[1]+vect[2]*vect[2]));
  if(norm<std::numeric_limits<double>::min())
    throw INTERP_KERNEL::Exception("DataArrayDouble::Rotate3DAlg : magnitude of input vector is too close of 0. !");
  std::transform(vect,vect+3,vectorNorm,std::bind2nd(std::multiplies<double>(),1/norm));
  //rotation matrix computation
  matrix[0]=cosa; matrix[1]=0.; matrix[2]=0.; matrix[3]=0.; matrix[4]=cosa; matrix[5]=0.; matrix[6]=0.; matrix[7]=0.; matrix[8]=cosa;
  matrixTmp[0]=vectorNorm[0]*vectorNorm[0]; matrixTmp[1]=vectorNorm[0]*vectorNorm[1]; matrixTmp[2]=vectorNorm[0]*vectorNorm[2];
  matrixTmp[3]=vectorNorm[1]*vectorNorm[0]; matrixTmp[4]=vectorNorm[1]*vectorNorm[1]; matrixTmp[5]=vectorNorm[1]*vectorNorm[2];
  matrixTmp[6]=vectorNorm[2]*vectorNorm[0]; matrixTmp[7]=vectorNorm[2]*vectorNorm[1]; matrixTmp[8]=vectorNorm[2]*vectorNorm[2];
  std::transform(matrixTmp,matrixTmp+9,matrixTmp,std::bind2nd(std::multiplies<double>(),1-cosa));
  std::transform(matrix,matrix+9,matrixTmp,matrix,std::plus<double>());
  matrixTmp[0]=0.; matrixTmp[1]=-vectorNorm[2]; matrixTmp[2]=vectorNorm[1];
  matrixTmp[3]=vectorNorm[2]; matrixTmp[4]=0.; matrixTmp[5]=-vectorNorm[0];
  matrixTmp[6]=-vectorNorm[1]; matrixTmp[7]=vectorNorm[0]; matrixTmp[8]=0.;
  std::transform(matrixTmp,matrixTmp+9,matrixTmp,std::bind2nd(std::multiplies<double>(),sina));
  std::transform(matrix,matrix+9,matrixTmp,matrix,std::plus<double>());
  //rotation matrix computed.
  double tmp[3];
  for(mcIdType i=0; i<nbNodes; i++)
    {
      std::transform(coordsIn+i*3,coordsIn+(i+1)*3,center,tmp,std::minus<double>());
      coordsOut[i*3]=matrix[0]*tmp[0]+matrix[1]*tmp[1]+matrix[2]*tmp[2]+center[0];
      coordsOut[i*3+1]=matrix[3]*tmp[0]+matrix[4]*tmp[1]+matrix[5]*tmp[2]+center[1];
      coordsOut[i*3+2]=matrix[6]*tmp[0]+matrix[7]*tmp[1]+matrix[8]*tmp[2]+center[2];
    }
}

void DataArrayDouble::Symmetry3DPlane(const double point[3], const double normalVector[3], mcIdType nbNodes, const double *coordsIn, double *coordsOut)
{
  double matrix[9],matrix2[9],matrix3[9];
  double vect[3],crossVect[3];
  INTERP_KERNEL::orthogonalVect3(normalVector,vect);
  crossVect[0]=normalVector[1]*vect[2]-normalVector[2]*vect[1];
  crossVect[1]=normalVector[2]*vect[0]-normalVector[0]*vect[2];
  crossVect[2]=normalVector[0]*vect[1]-normalVector[1]*vect[0];
  double nv(INTERP_KERNEL::norm<3>(vect)),ni(INTERP_KERNEL::norm<3>(normalVector)),nc(INTERP_KERNEL::norm<3>(crossVect));
  matrix[0]=vect[0]/nv; matrix[1]=crossVect[0]/nc; matrix[2]=-normalVector[0]/ni;
  matrix[3]=vect[1]/nv; matrix[4]=crossVect[1]/nc; matrix[5]=-normalVector[1]/ni;
  matrix[6]=vect[2]/nv; matrix[7]=crossVect[2]/nc; matrix[8]=-normalVector[2]/ni;
  matrix2[0]=vect[0]/nv; matrix2[1]=vect[1]/nv; matrix2[2]=vect[2]/nv;
  matrix2[3]=crossVect[0]/nc; matrix2[4]=crossVect[1]/nc; matrix2[5]=crossVect[2]/nc;
  matrix2[6]=normalVector[0]/ni; matrix2[7]=normalVector[1]/ni; matrix2[8]=normalVector[2]/ni;
  for(mcIdType i=0;i<3;i++)
    for(mcIdType j=0;j<3;j++)
      {
        double val(0.);
        for(mcIdType k=0;k<3;k++)
          val+=matrix[3*i+k]*matrix2[3*k+j];
        matrix3[3*i+j]=val;
      }
  //rotation matrix computed.
  double tmp[3];
  for(mcIdType i=0; i<nbNodes; i++)
    {
      std::transform(coordsIn+i*3,coordsIn+(i+1)*3,point,tmp,std::minus<double>());
      coordsOut[i*3]=matrix3[0]*tmp[0]+matrix3[1]*tmp[1]+matrix3[2]*tmp[2]+point[0];
      coordsOut[i*3+1]=matrix3[3]*tmp[0]+matrix3[4]*tmp[1]+matrix3[5]*tmp[2]+point[1];
      coordsOut[i*3+2]=matrix3[6]*tmp[0]+matrix3[7]*tmp[1]+matrix3[8]*tmp[2]+point[2];
    }
}

void DataArrayDouble::GiveBaseForPlane(const double normalVector[3], double baseOfPlane[9])
{
  double vect[3],crossVect[3];
  INTERP_KERNEL::orthogonalVect3(normalVector,vect);
  crossVect[0]=normalVector[1]*vect[2]-normalVector[2]*vect[1];
  crossVect[1]=normalVector[2]*vect[0]-normalVector[0]*vect[2];
  crossVect[2]=normalVector[0]*vect[1]-normalVector[1]*vect[0];
  double nv(INTERP_KERNEL::norm<3>(vect)),ni(INTERP_KERNEL::norm<3>(normalVector)),nc(INTERP_KERNEL::norm<3>(crossVect));
  baseOfPlane[0]=vect[0]/nv; baseOfPlane[1]=vect[1]/nv; baseOfPlane[2]=vect[2]/nv;
  baseOfPlane[3]=crossVect[0]/nc; baseOfPlane[4]=crossVect[1]/nc; baseOfPlane[5]=crossVect[2]/nc;
  baseOfPlane[6]=normalVector[0]/ni; baseOfPlane[7]=normalVector[1]/ni; baseOfPlane[8]=normalVector[2]/ni;
}

/*!
 * \param [in] seg2 : coordinates of input seg2 expected to have spacedim==2
 * \param [in] tri3 : coordinates of input tri3 also expected to have spacedim==2
 * \param [out] coeffs : the result of integration normalized to 1. along \a seg2 inside tri3 sorted by the node id of \a tri3
 * \param [out] length : the length of seg2. That is too say the length of integration
 */
void DataArrayDouble::ComputeIntegralOfSeg2IntoTri3(const double seg2[4], const double tri3[6], double coeffs[3], double& length)
{
  length=INTERP_KERNEL::norme_vecteur(seg2,seg2+2);
  double mid[2];
  INTERP_KERNEL::mid_of_seg2(seg2,seg2+2,mid);
  INTERP_KERNEL::barycentric_coords<2>(tri3,mid,coeffs); // integral along seg2 is equal to value at the center of SEG2 !
}

/*!
 * Low static method that operates 3D rotation of \a nbNodes 3D nodes whose coordinates are arranged in \a coords
 * around the center point \a center and with angle \a angle.
 */
void DataArrayDouble::Rotate2DAlg(const double *center, double angle, mcIdType nbNodes, const double *coordsIn, double *coordsOut)
{
  double cosa=cos(angle);
  double sina=sin(angle);
  double matrix[4];
  matrix[0]=cosa; matrix[1]=-sina; matrix[2]=sina; matrix[3]=cosa;
  double tmp[2];
  for(mcIdType i=0; i<nbNodes; i++)
    {
      std::transform(coordsIn+i*2,coordsIn+(i+1)*2,center,tmp,std::minus<double>());
      coordsOut[i*2]=matrix[0]*tmp[0]+matrix[1]*tmp[1]+center[0];
      coordsOut[i*2+1]=matrix[2]*tmp[0]+matrix[3]*tmp[1]+center[1];
    }
}

DataArrayDoubleIterator::DataArrayDoubleIterator(DataArrayDouble *da):DataArrayIterator<double>(da)
{
}

DataArrayDoubleTuple::DataArrayDoubleTuple(double *pt, std::size_t nbOfComp):DataArrayTuple<double>(pt,nbOfComp)
{
}


std::string DataArrayDoubleTuple::repr() const
{
  std::ostringstream oss; oss.precision(17); oss << "(";
  for(std::size_t i=0;i<_nb_of_compo-1;i++)
    oss << _pt[i] << ", ";
  oss << _pt[_nb_of_compo-1] << ")";
  return oss.str();
}

double DataArrayDoubleTuple::doubleValue() const
{
  return this->zeValue();
}

/*!
 * This method returns a newly allocated instance the caller should dealed with by a MEDCoupling::DataArrayDouble::decrRef.
 * This method performs \b no copy of data. The content is only referenced using MEDCoupling::DataArrayDouble::useArray with ownership set to \b false.
 * This method throws an INTERP_KERNEL::Exception is it is impossible to match sizes of \b this that is too say \b nbOfCompo=this->_nb_of_elem and \bnbOfTuples==1 or
 * \b nbOfCompo=1 and \bnbOfTuples==this->_nb_of_elem.
 */
DataArrayDouble *DataArrayDoubleTuple::buildDADouble(std::size_t nbOfTuples, std::size_t nbOfCompo) const
{
  return this->buildDA(nbOfTuples,nbOfCompo);
}

/*!
 * Returns a full copy of \a this. For more info on copying data arrays see
 * \ref MEDCouplingArrayBasicsCopyDeep.
 *  \return DataArrayInt * - a new instance of DataArrayInt.
 */
DataArrayInt32 *DataArrayInt32::deepCopy() const
{
  return new DataArrayInt32(*this);
}

DataArrayInt32Iterator *DataArrayInt32::iterator()
{
  return new DataArrayInt32Iterator(this);
}


DataArrayInt32Iterator::DataArrayInt32Iterator(DataArrayInt32 *da):DataArrayIterator<Int32>(da)
{
}

DataArrayInt32Tuple::DataArrayInt32Tuple(Int32 *pt, std::size_t nbOfComp):DataArrayTuple<Int32>(pt,nbOfComp)
{
}

std::string DataArrayInt32Tuple::repr() const
{
  std::ostringstream oss; oss << "(";
  for(std::size_t i=0;i<_nb_of_compo-1;i++)
    oss << _pt[i] << ", ";
  oss << _pt[_nb_of_compo-1] << ")";
  return oss.str();
}

Int32 DataArrayInt32Tuple::intValue() const
{
  return this->zeValue();
}

/*!
 * This method returns a newly allocated instance the caller should dealed with by a MEDCoupling::DataArrayInt::decrRef.
 * This method performs \b no copy of data. The content is only referenced using MEDCoupling::DataArrayInt::useArray with ownership set to \b false.
 * This method throws an INTERP_KERNEL::Exception is it is impossible to match sizes of \b this that is too say \b nbOfCompo=this->_nb_of_elem and \bnbOfTuples==1 or
 * \b nbOfCompo=1 and \bnbOfTuples==this->_nb_of_elem.
 */
DataArrayInt32 *DataArrayInt32Tuple::buildDAInt(std::size_t nbOfTuples, std::size_t nbOfCompo) const
{
  return this->buildDA(nbOfTuples,nbOfCompo);
}

DataArrayInt64Iterator *DataArrayInt64::iterator()
{
  return new DataArrayInt64Iterator(this);
}


DataArrayInt64Iterator::DataArrayInt64Iterator(DataArrayInt64 *da):DataArrayIterator<Int64>(da)
{
}

DataArrayInt64Tuple::DataArrayInt64Tuple(Int64 *pt, std::size_t nbOfComp):DataArrayTuple<Int64>(pt,nbOfComp)
{
}

std::string DataArrayInt64Tuple::repr() const
{
  std::ostringstream oss; oss << "(";
  for(std::size_t i=0;i<_nb_of_compo-1;i++)
    oss << _pt[i] << ", ";
  oss << _pt[_nb_of_compo-1] << ")";
  return oss.str();
}

Int64 DataArrayInt64Tuple::intValue() const
{
  return this->zeValue();
}

DataArrayInt64 *DataArrayInt64Tuple::buildDAInt(std::size_t nbOfTuples, std::size_t nbOfCompo) const
{
  return this->buildDA(nbOfTuples,nbOfCompo);
}


DataArrayInt64 *DataArrayInt64::deepCopy() const
{
  return new DataArrayInt64(*this);
}
