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

#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldTemplate.hxx"
#include "MEDCouplingFieldT.txx"
#include "MEDCouplingFieldInt.hxx"
#include "MEDCouplingFieldFloat.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MCAuto.txx"
#include "MEDCouplingVoronoi.hxx"
#include "MEDCouplingNatureOfField.hxx"
#include "MEDCouplingMemArray.txx"

#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGaussCoords.hxx"

#include <sstream>
#include <limits>
#include <algorithm>
#include <functional>

using namespace MEDCoupling;

template class MEDCoupling::MEDCouplingFieldT<double>;

/*!
 * Creates a new MEDCouplingFieldDouble, of given spatial type and time discretization.
 * For more info, see \ref MEDCouplingFirstSteps3.
 * \param [in] type - the type of spatial discretization of the created field, one of
 *        (\ref MEDCoupling::ON_CELLS "ON_CELLS", 
 *         \ref MEDCoupling::ON_NODES "ON_NODES",
 *         \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT", 
 *         \ref MEDCoupling::ON_GAUSS_NE "ON_GAUSS_NE",
 *         \ref MEDCoupling::ON_NODES_KR "ON_NODES_KR").
 * \param [in] td - the type of time discretization of the created field, one of
 *        (\ref MEDCoupling::NO_TIME "NO_TIME", 
 *         \ref MEDCoupling::ONE_TIME "ONE_TIME", 
 *         \ref MEDCoupling::LINEAR_TIME "LINEAR_TIME", 
 *         \ref MEDCoupling::CONST_ON_TIME_INTERVAL "CONST_ON_TIME_INTERVAL").
 * \return MEDCouplingFieldDouble* - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed. 
 */
MEDCouplingFieldDouble* MEDCouplingFieldDouble::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldDouble(type,td);
}

/*!
 * Creates a new MEDCouplingFieldDouble, of a given time discretization and with a
 * spatial type and supporting mesh copied from a given 
 * \ref MEDCouplingFieldTemplatesPage "field template".
 * For more info, see \ref MEDCouplingFirstSteps3.
 * \warning This method does not deeply copy neither the mesh nor the spatial
 * discretization. Only a shallow copy (reference) is done for the mesh and the spatial
 * discretization!
 * \param [in] ft - the \ref MEDCouplingFieldTemplatesPage "field template" defining
 *        the spatial discretization and the supporting mesh.
 * \param [in] td - the type of time discretization of the created field, one of
 *        (\ref MEDCoupling::NO_TIME "NO_TIME", 
 *         \ref MEDCoupling::ONE_TIME "ONE_TIME", 
 *         \ref MEDCoupling::LINEAR_TIME "LINEAR_TIME", 
 *         \ref MEDCoupling::CONST_ON_TIME_INTERVAL "CONST_ON_TIME_INTERVAL").
 * \return MEDCouplingFieldDouble* - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed. 
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::New(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldDouble(ft,td);
}

/*!
 * Sets a time \a unit of \a this field. For more info, see \ref MEDCouplingFirstSteps3.
 * \param [in] unit \a unit (string) in which time is measured.
 */
//void MEDCouplingFieldDouble::setTimeUnit(const std::string& unit)

/*!
 * Returns a time unit of \a this field.
 * \return a string describing units in which time is measured.
 */
//std::string MEDCouplingFieldDouble::getTimeUnit() const


/*!
 * This method if possible the time information (time unit, time iteration, time unit and time value) with its support
 * that is to say its mesh.
 * 
 * \throw  If \c this->_mesh is null an exception will be thrown. An exception will also be throw if the spatial discretization is
 *         NO_TIME.
 */
void MEDCouplingFieldDouble::synchronizeTimeWithSupport()
{
  timeDiscr()->synchronizeTimeWith(_mesh);
}

/*!
 * Returns a new MEDCouplingFieldDouble which is a copy of \a this one. The data
 * of \a this field is copied either deep or shallow depending on \a recDeepCpy
 * parameter. But the underlying mesh is always shallow copied.
 * Data that can be copied either deeply or shallow are:
 * - \ref MEDCouplingTemporalDisc "temporal discretization" data that holds array(s)
 * of field values,
 * - \ref MEDCouplingSpatialDisc "a spatial discretization".
 * 
 * \c clone(false) is rather dedicated for advanced users that want to limit the amount
 * of memory. It allows the user to perform methods like operator+(), operator*()
 * etc. with \a this and the returned field. If the user wants to duplicate deeply the
 * underlying mesh he should call cloneWithMesh() method or deepCopy() instead. 
 * \warning The underlying \b mesh of the returned field is **always the same**
 *         (pointer) as \a this one **whatever the value** of \a recDeepCpy parameter.
 *  \param [in] recDeepCpy - if \c true, the copy of the underlying data arrays is
 *         deep, else all data arrays of \a this field are shared by the new field.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed.
 * \sa cloneWithMesh()
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::clone(bool recDeepCpy) const
{
  return new MEDCouplingFieldDouble(*this,recDeepCpy);
}

/*!
 * Returns a new MEDCouplingFieldDouble which is a deep copy of \a this one **including
 * the mesh**.
 * The result of this method is exactly the same as that of \c cloneWithMesh(true).
 * So the resulting field can not be used together with \a this one in the methods
 * like operator+(), operator*() etc. To avoid deep copying the underlying mesh,
 * the user can call clone().
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed.
 * \sa cloneWithMesh()
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::deepCopy() const
{
  return cloneWithMesh(true);
}

/*!
 * Creates a new MEDCouplingFieldDouble of given
 * \ref MEDCouplingTemporalDisc "temporal discretization". The result field either
 * shares the data array(s) with \a this field, or holds a deep copy of it, depending on
 * \a deepCopy parameter. But the underlying \b mesh is always **shallow copied**.
 * \param [in] td - the type of time discretization of the created field, one of
 *        (\ref MEDCoupling::NO_TIME "NO_TIME", 
 *         \ref MEDCoupling::ONE_TIME "ONE_TIME", 
 *         \ref MEDCoupling::LINEAR_TIME "LINEAR_TIME", 
 *         \ref MEDCoupling::CONST_ON_TIME_INTERVAL "CONST_ON_TIME_INTERVAL").
 * \param [in] deepCopy - if \c true, the copy of the underlying data arrays is
 *         deep, else all data arrays of \a this field are shared by the new field.
 * \return MEDCouplingFieldDouble* - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed. 
 * 
 * \if ENABLE_EXAMPLES
 * \ref cpp_mcfielddouble_buildNewTimeReprFromThis "Here is a C++ example."<br>
 * \ref py_mcfielddouble_buildNewTimeReprFromThis "Here is a Python example."
 * \endif
 * \sa clone()
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::buildNewTimeReprFromThis(TypeOfTimeDiscretization td, bool deepCpy) const
{
  MEDCouplingTimeDiscretization *tdo=timeDiscr()->buildNewTimeReprFromThis(td,deepCpy);
  MCAuto<MEDCouplingFieldDiscretization> disc;
  if(_type)
    disc=_type->clone();
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),tdo,disc.retn()));
  ret->setMesh(getMesh());
  ret->setName(getName());
  ret->setDescription(getDescription());
  return ret.retn();
}

/*!
 * This method converts a field on nodes (\a this) to a cell field (returned field). The convertion is a \b non \b conservative remapping !
 * This method is useful only for users that need a fast convertion from node to cell spatial discretization. The algorithm applied is simply to attach
 * to each cell the average of values on nodes constituting this cell.
 *
 * \return MEDCouplingFieldDouble* - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed. The returned field will share the same mesh object object than those in \a this.
 * \throw If \a this spatial discretization is empty or not ON_NODES.
 * \throw If \a this is not coherent (see MEDCouplingFieldDouble::checkConsistencyLight).
 * 
 * \warning This method is a \b non \b conservative method of remapping from node spatial discretization to cell spatial discretization.
 * If a conservative method of interpolation is required MEDCoupling::MEDCouplingRemapper class should be used instead with "P1P0" method.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::nodeToCellDiscretization() const
{
  checkConsistencyLight();
  TypeOfField tf(getTypeOfField());
  if(tf!=ON_NODES)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::nodeToCellDiscretization : this field is expected to be on ON_NODES !");
  MCAuto<MEDCouplingFieldDouble> ret(clone(false));
  MCAuto<MEDCouplingFieldDiscretizationP0> nsp(new MEDCouplingFieldDiscretizationP0);
  ret->setDiscretization(nsp);
  const MEDCouplingMesh *m(getMesh());//m is non empty thanks to checkConsistencyLight call
  int nbCells(m->getNumberOfCells());
  std::vector<DataArrayDouble *> arrs(getArrays());
  std::size_t sz(arrs.size());
  std::vector< MCAuto<DataArrayDouble> > outArrsSafe(sz); std::vector<DataArrayDouble *> outArrs(sz);
  for(std::size_t j=0;j<sz;j++)
    {
      int nbCompo(arrs[j]->getNumberOfComponents());
      outArrsSafe[j]=DataArrayDouble::New(); outArrsSafe[j]->alloc(nbCells,nbCompo);
      outArrsSafe[j]->copyStringInfoFrom(*arrs[j]);
      outArrs[j]=outArrsSafe[j];
      double *pt(outArrsSafe[j]->getPointer());
      const double *srcPt(arrs[j]->begin());
      for(int i=0;i<nbCells;i++,pt+=nbCompo)
        {
          std::vector<int> nodeIds;
          m->getNodeIdsOfCell(i,nodeIds);
          std::fill(pt,pt+nbCompo,0.);
          std::size_t nbNodesInCell(nodeIds.size());
          for(std::size_t k=0;k<nbNodesInCell;k++)
            std::transform(srcPt+nodeIds[k]*nbCompo,srcPt+(nodeIds[k]+1)*nbCompo,pt,pt,std::plus<double>());
          if(nbNodesInCell!=0)
            std::transform(pt,pt+nbCompo,pt,std::bind2nd(std::multiplies<double>(),1./((double)nbNodesInCell)));
          else
            {
              std::ostringstream oss; oss << "MEDCouplingFieldDouble::nodeToCellDiscretization : Cell id #" << i << " has been detected to have no nodes !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
        }
    }
  ret->setArrays(outArrs);
  return ret.retn();
}

/*!
 * This method converts a field on cell (\a this) to a node field (returned field). The convertion is a \b non \b conservative remapping !
 * This method is useful only for users that need a fast convertion from cell to node spatial discretization. The algorithm applied is simply to attach
 * to each node the average of values on cell sharing this node. If \a this lies on a mesh having orphan nodes the values applied on them will be NaN (division by 0.).
 *
 * \return MEDCouplingFieldDouble* - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed. The returned field will share the same mesh object object than those in \a this.
 * \throw If \a this spatial discretization is empty or not ON_CELLS.
 * \throw If \a this is not coherent (see MEDCouplingFieldDouble::checkConsistencyLight).
 *
 * \warning This method is a \b non \b conservative method of remapping from cell spatial discretization to node spatial discretization.
 * If a conservative method of interpolation is required MEDCoupling::MEDCouplingRemapper class should be used instead with "P0P1" method.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::cellToNodeDiscretization() const
{
  checkConsistencyLight();
  TypeOfField tf(getTypeOfField());
  if(tf!=ON_CELLS)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::cellToNodeDiscretization : this field is expected to be on ON_CELLS !");
  MCAuto<MEDCouplingFieldDouble> ret(clone(false));
  MCAuto<MEDCouplingFieldDiscretizationP1> nsp(new MEDCouplingFieldDiscretizationP1);
  ret->setDiscretization(nsp);
  const MEDCouplingMesh *m(getMesh());//m is non empty thanks to checkConsistencyLight call
  MCAuto<DataArrayInt> rn(DataArrayInt::New()),rni(DataArrayInt::New());
  m->getReverseNodalConnectivity(rn,rni);
  MCAuto<DataArrayInt> rni2(rni->deltaShiftIndex());
  MCAuto<DataArrayDouble> rni3(rni2->convertToDblArr()); rni2=0;
  std::vector<DataArrayDouble *> arrs(getArrays());
  std::size_t sz(arrs.size());
  std::vector< MCAuto<DataArrayDouble> > outArrsSafe(sz); std::vector<DataArrayDouble *> outArrs(sz);
  for(std::size_t j=0;j<sz;j++)
    {
      MCAuto<DataArrayDouble> tmp(arrs[j]->selectByTupleIdSafe(rn->begin(),rn->end()));
      outArrsSafe[j]=(tmp->accumulatePerChunck(rni->begin(),rni->end())); tmp=0;
      outArrsSafe[j]->divideEqual(rni3);
      outArrsSafe[j]->copyStringInfoFrom(*arrs[j]);
      outArrs[j]=outArrsSafe[j];
    }
  ret->setArrays(outArrs);
  return ret.retn();
}

/*!
 * Returns a string describing \a this field. The string includes info on
 * - name,
 * - description,
 * - \ref MEDCouplingSpatialDisc "spatial discretization",
 * - \ref MEDCouplingTemporalDisc "time discretization",
 * - components,
 * - mesh,
 * - contents of data arrays.
 *
 *  \return std::string - the string describing \a this field.
 */
std::string MEDCouplingFieldDouble::advancedRepr() const
{
  std::ostringstream ret;
  ret << "FieldDouble with name : \"" << getName() << "\"\n";
  ret << "Description of field is : \"" << getDescription() << "\"\n";
  if(_type)
    { ret << "FieldDouble space discretization is : " << _type->getStringRepr() << "\n"; }
  else
    { ret << "FieldDouble has no space discretization set !\n"; }
  if(timeDiscr())
    { ret << "FieldDouble time discretization is : " << timeDiscr()->getStringRepr() << "\n"; }
  else
    { ret << "FieldDouble has no time discretization set !\n"; }
  if(getArray())
    ret << "FieldDouble default array has " << getArray()->getNumberOfComponents() << " components and " << getArray()->getNumberOfTuples() << " tuples.\n";
  if(_mesh)
    ret << "Mesh support information :\n__________________________\n" << _mesh->advancedRepr();
  else
    ret << "Mesh support information : No mesh set !\n";
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  int arrayId=0;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++,arrayId++)
    {
      ret << "Array #" << arrayId << " :\n__________\n";
      if(*iter)
        (*iter)->reprWithoutNameStream(ret);
      else
        ret << "Array empty !";
      ret << "\n";
    }
  return ret.str();
}

std::string MEDCouplingFieldDouble::writeVTK(const std::string& fileName, bool isBinary) const
{
  std::vector<const MEDCouplingFieldDouble *> fs(1,this);
  return MEDCouplingFieldDouble::WriteVTK(fileName,fs,isBinary);
}

/*!
 * This method states if \a this and 'other' are compatibles each other before performing any treatment.
 * This method is good for methods like : mergeFields.
 * This method is not very demanding compared to areStrictlyCompatible that is better for operation on fields.
 */
bool MEDCouplingFieldDouble::areCompatibleForMerge(const MEDCouplingField *other) const
{
  if(!MEDCouplingField::areCompatibleForMerge(other))
    return false;
  const MEDCouplingFieldDouble *otherC(dynamic_cast<const MEDCouplingFieldDouble *>(other));
  if(!otherC)
    return false;
  if(!timeDiscr()->areCompatible(otherC->timeDiscr()))
    return false;
  return true;
}

/*!
 * This method is invocated before any attempt of melding. This method is very close to areStrictlyCompatible,
 * except that \a this and other can have different number of components.
 */
bool MEDCouplingFieldDouble::areCompatibleForMeld(const MEDCouplingFieldDouble *other) const
{
  if(!MEDCouplingField::areStrictlyCompatible(other))
    return false;
  if(!timeDiscr()->areCompatibleForMeld(other->timeDiscr()))
    return false;
  return true;
}

/*!
 * Permutes values of \a this field according to a given permutation array for node
 * renumbering. The underlying mesh is deeply copied and its nodes are also permuted. 
 * The number of nodes can change, contrary to renumberCells().
 *  \param [in] old2NewBg - the permutation array in "Old to New" mode. Its length is
 *         to be equal to \a this->getMesh()->getNumberOfNodes().
 *  \param [in] eps - a precision used to compare field values at merged nodes. If
 *         the values differ more than \a eps, an exception is thrown.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If \a check == \c true and \a old2NewBg contains equal ids.
 *  \throw If mesh nature does not allow renumbering (e.g. structured mesh).
 *  \throw If values at merged nodes deffer more than \a eps.
 * 
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_renumberNodes "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_renumberNodes "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::renumberNodes(const int *old2NewBg, double eps)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("Invalid mesh to apply renumberNodes on it !");
  int nbOfNodes=meshC->getNumberOfNodes();
  MCAuto<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCopy());
  int newNbOfNodes=*std::max_element(old2NewBg,old2NewBg+nbOfNodes)+1;
  renumberNodesWithoutMesh(old2NewBg,newNbOfNodes,eps);
  meshC2->renumberNodes(old2NewBg,newNbOfNodes);
  setMesh(meshC2);
}

/*!
 * Permutes values of \a this field according to a given permutation array for nodes
 * renumbering. The underlying mesh is \b not permuted. 
 * The number of nodes can change, contrary to renumberCells().
 * A given epsilon specifies a threshold of error in case of two nodes are merged but
 * the difference of values on these nodes are higher than \a eps.
 * This method performs a part of job of renumberNodes(), excluding node renumbering
 * in mesh. The reasonable use of this
 * method is only for multi-field instances lying on the same mesh to avoid a
 * systematic duplication and renumbering of _mesh attribute. 
 * \warning Use this method with a lot of care!
 * \warning In case of an exception thrown, the contents of the data array can be
 *         partially modified until the exception occurs. 
 *  \param [in] old2NewBg - the permutation array in "Old to New" mode. Its length is
 *         to be equal to \a this->getMesh()->getNumberOfNodes().
 *  \param [in] newNbOfNodes - a number of nodes in the mesh after renumbering.
 *  \param [in] eps - a precision used to compare field values at merged nodes. If
 *         the values differ more than \a eps, an exception is thrown.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If values at merged nodes deffer more than \a eps.
 */
void MEDCouplingFieldDouble::renumberNodesWithoutMesh(const int *old2NewBg, int newNbOfNodes, double eps)
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("Expecting a spatial discretization to be able to operate a renumbering !");
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    if(*iter)
      _type->renumberValuesOnNodes(eps,old2NewBg,newNbOfNodes,*iter);
}

/*!
 * Returns all tuple ids of \a this scalar field that fit the range [\a vmin,
 * \a vmax]. This method calls DataArrayDouble::findIdsInRange().
 *  \param [in] vmin - a lower boundary of the range. Tuples with values less than \a
 *         vmin are not included in the result array.
 *  \param [in] vmax - an upper boundary of the range. Tuples with values more than \a
 *         vmax are not included in the result array.
 *  \return DataArrayInt * - a new instance of DataArrayInt holding ids of selected
 *          tuples. The caller is to delete this array using decrRef() as it is no
 *          more needed.
 *  \throw If the data array is not set.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
DataArrayInt *MEDCouplingFieldDouble::findIdsInRange(double vmin, double vmax) const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::findIdsInRange : no default array set !");
  return getArray()->findIdsInRange(vmin,vmax);
}

template<class U>
typename Traits<U>::FieldType *ConvertToUField(const MEDCouplingFieldDouble *self)
{
  MCAuto<MEDCouplingFieldTemplate> tmp(MEDCouplingFieldTemplate::New(*self));
  int t1,t2;
  double t0(self->getTime(t1,t2));
  MCAuto<typename Traits<U>::FieldType > ret(Traits<U>::FieldType::New(*tmp,self->getTimeDiscretization()));
  ret->setTime(t0,t1,t2);
  if(self->getArray())
    {
      MCAuto<typename Traits<U>::ArrayType> arr(self->getArray()->convertToOtherTypeOfArr<U>());
      ret->setArray(arr);
    }
  return ret.retn();
}

MEDCouplingFieldInt *MEDCouplingFieldDouble::convertToIntField() const
{
  return ConvertToUField<int>(this);
}

MEDCouplingFieldFloat *MEDCouplingFieldDouble::convertToFloatField() const
{
  return ConvertToUField<float>(this);
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingFieldT<double>(type,MEDCouplingTimeDiscretization::New(td))
{
}

/*!
 * ** WARINING : This method do not deeply copy neither mesh nor spatial discretization. Only a shallow copy (reference) is done for mesh and spatial discretization ! **
 */
MEDCouplingFieldDouble::MEDCouplingFieldDouble(const MEDCouplingFieldTemplate& ft, TypeOfTimeDiscretization td):MEDCouplingFieldT<double>(ft,MEDCouplingTimeDiscretization::New(td),false)
{
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCpy):MEDCouplingFieldT<double>(other,deepCpy)
{
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(NatureOfField n, MEDCouplingTimeDiscretization *td, MEDCouplingFieldDiscretization *type):MEDCouplingFieldT<double>(type,n,td)
{
}

/*!
 * Accumulate values of a given component of \a this field.
 *  \param [in] compId - the index of the component of interest.
 *  \return double - a sum value of *compId*-th component.
 *  \throw If the data array is not set.
 *  \throw If \a the condition ( 0 <= \a compId < \a this->getNumberOfComponents() ) is
 *         not respected.
 */
double MEDCouplingFieldDouble::accumulate(int compId) const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::accumulate : no default array defined !");
  return getArray()->accumulate(compId);
}

/*!
 * Accumulates values of each component of \a this array.
 *  \param [out] res - an array of length \a this->getNumberOfComponents(), allocated 
 *         by the caller, that is filled by this method with sum value for each
 *         component.
 *  \throw If the data array is not set.
 */
void MEDCouplingFieldDouble::accumulate(double *res) const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::accumulate : no default array defined !");
  getArray()->accumulate(res);
}

/*!
 * Returns the maximal value within \a this scalar field. Values of all arrays stored
 * in \a this->_time_discr are checked.
 *  \return double - the maximal value among all values of \a this field.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If the data array is not set.
 *  \throw If there is an empty data array in \a this field.
 */
double MEDCouplingFieldDouble::getMaxValue() const
{
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  double ret(-std::numeric_limits<double>::max());
  bool isExistingArr=false;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    {
      if(*iter)
        {
          isExistingArr=true;
          int loc;
          ret=std::max(ret,(*iter)->getMaxValue(loc));
        }
    }
  if(!isExistingArr)
    throw INTERP_KERNEL::Exception("getMaxValue : No arrays defined !");
  return ret;
}

/*!
 * Returns the maximal value and all its locations within \a this scalar field.
 * Only the first of available data arrays is checked.
 *  \param [out] tupleIds - a new instance of DataArrayInt containg indices of
 *               tuples holding the maximal value. The caller is to delete it using
 *               decrRef() as it is no more needed.
 *  \return double - the maximal value among all values of the first array of \a this filed.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If there is an empty data array in \a this field.
 */
double MEDCouplingFieldDouble::getMaxValue2(DataArrayInt*& tupleIds) const
{
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  double ret(-std::numeric_limits<double>::max());
  bool isExistingArr=false;
  tupleIds=0;
  MCAuto<DataArrayInt> ret1;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    {
      if(*iter)
        {
          isExistingArr=true;
          DataArrayInt *tmp;
          ret=std::max(ret,(*iter)->getMaxValue2(tmp));
          MCAuto<DataArrayInt> tmpSafe(tmp);
          if(!((const DataArrayInt *)ret1))
            ret1=tmpSafe;
        }
    }
  if(!isExistingArr)
    throw INTERP_KERNEL::Exception("getMaxValue2 : No arrays defined !");
  tupleIds=ret1.retn();
  return ret;
}

/*!
 * Returns the minimal value within \a this scalar field. Values of all arrays stored
 * in \a this->_time_discr are checked.
 *  \return double - the minimal value among all values of \a this field.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If the data array is not set.
 *  \throw If there is an empty data array in \a this field.
 */
double MEDCouplingFieldDouble::getMinValue() const
{
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  double ret(std::numeric_limits<double>::max());
  bool isExistingArr=false;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    {
      if(*iter)
        {
          isExistingArr=true;
          int loc;
          ret=std::min(ret,(*iter)->getMinValue(loc));
        }
    }
  if(!isExistingArr)
    throw INTERP_KERNEL::Exception("getMinValue : No arrays defined !");
  return ret;
}

/*!
 * Returns the minimal value and all its locations within \a this scalar field.
 * Only the first of available data arrays is checked.
 *  \param [out] tupleIds - a new instance of DataArrayInt containg indices of
 *               tuples holding the minimal value. The caller is to delete it using
 *               decrRef() as it is no more needed.
 *  \return double - the minimal value among all values of the first array of \a this filed.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If there is an empty data array in \a this field.
 */
double MEDCouplingFieldDouble::getMinValue2(DataArrayInt*& tupleIds) const
{
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  double ret(-std::numeric_limits<double>::max());
  bool isExistingArr=false;
  tupleIds=0;
  MCAuto<DataArrayInt> ret1;
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    {
      if(*iter)
        {
          isExistingArr=true;
          DataArrayInt *tmp;
          ret=std::max(ret,(*iter)->getMinValue2(tmp));
          MCAuto<DataArrayInt> tmpSafe(tmp);
          if(!((const DataArrayInt *)ret1))
            ret1=tmpSafe;
        }
    }
  if(!isExistingArr)
    throw INTERP_KERNEL::Exception("getMinValue2 : No arrays defined !");
  tupleIds=ret1.retn();
  return ret;
}

/*!
 * Returns the average value of \a this scalar field.
 *  \return double - the average value over all values of the data array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If the data array is not set or it is empty.
 */
double MEDCouplingFieldDouble::getAverageValue() const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::getAverageValue : no default array defined !");
  return getArray()->getAverageValue();
}

/*!
 * This method returns the euclidean norm of \a this field.
 * \f[
 * \sqrt{\sum_{0 \leq i < nbOfEntity}val[i]*val[i]}
 * \f]
 *  \throw If the data array is not set.
 */
double MEDCouplingFieldDouble::norm2() const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::norm2 : no default array defined !");
  return getArray()->norm2();
}

/*!
 * This method returns the max norm of \a this field.
 * \f[
 * \max_{0 \leq i < nbOfEntity}{abs(val[i])}
 * \f]
 *  \throw If the data array is not set.
 */
double MEDCouplingFieldDouble::normMax() const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::normMax : no default array defined !");
  return getArray()->normMax();
}

/*!
 * Computes the weighted average of values of each component of \a this field, the weights being the
 * values returned by buildMeasureField().  
 *  \param [out] res - pointer to an array of result sum values, of size at least \a
 *         this->getNumberOfComponents(), that is to be allocated by the caller.
 *  \param [in] isWAbs - if \c true (default), \c abs() is applied to the weights computed by
 *         buildMeasureField(). It makes this method slower. If you are sure that all
 *         the cells of the underlying mesh have a correct orientation (no negative volume), you can put \a isWAbs ==
 *         \c false to speed up the method.
 *  \throw If the mesh is not set.
 *  \throw If the data array is not set.
 */
void MEDCouplingFieldDouble::getWeightedAverageValue(double *res, bool isWAbs) const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::getWeightedAverageValue : no default array defined !");
  MCAuto<MEDCouplingFieldDouble> w=buildMeasureField(isWAbs);
  double deno=w->getArray()->accumulate(0);
  MCAuto<DataArrayDouble> arr=getArray()->deepCopy();
  arr->multiplyEqual(w->getArray());
  arr->accumulate(res);
  int nCompo = getArray()->getNumberOfComponents();
  std::transform(res,res+nCompo,res,std::bind2nd(std::multiplies<double>(),1./deno));
}

/*!
 * Computes the weighted average of values of a given component of \a this field, the weights being the
 * values returned by buildMeasureField().
 *  \param [in] compId - an index of the component of interest.
 *  \param [in] isWAbs - if \c true (default), \c abs() is applied to the weights computed by
 *         buildMeasureField(). It makes this method slower. If you are sure that all
 *         the cells of the underlying mesh have a correct orientation (no negative volume), you can put \a isWAbs ==
 *         \c false to speed up the method.
 *  \throw If the mesh is not set.
 *  \throw If the data array is not set.
 *  \throw If \a compId is not valid.
           A valid range is ( 0 <= \a compId < \a this->getNumberOfComponents() ).
 */
double MEDCouplingFieldDouble::getWeightedAverageValue(int compId, bool isWAbs) const
{
  int nbComps=getArray()->getNumberOfComponents();
  if(compId<0 || compId>=nbComps)
    {
      std::ostringstream oss; oss << "MEDCouplingFieldDouble::getWeightedAverageValue : Invalid compId specified : No such nb of components ! Should be in [0," << nbComps << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  INTERP_KERNEL::AutoPtr<double> res=new double[nbComps];
  getWeightedAverageValue(res,isWAbs);
  return res[compId];
}

/*!
 * Returns the \c normL1 of values of a given component of \a this field:
 * \f[
 * \frac{\sum_{0 \leq i < nbOfEntity}|val[i]*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}
 * \f]
 *  \param [in] compId - an index of the component of interest.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If \a compId is not valid.
           A valid range is ( 0 <= \a compId < \a this->getNumberOfComponents() ).
 */
double MEDCouplingFieldDouble::normL1(int compId) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL1 !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform normL1 !");
  int nbComps=getArray()->getNumberOfComponents();
  if(compId<0 || compId>=nbComps)
    {
      std::ostringstream oss; oss << "MEDCouplingFieldDouble::normL1 : Invalid compId specified : No such nb of components ! Should be in [0," << nbComps << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  INTERP_KERNEL::AutoPtr<double> res=new double[nbComps];
  _type->normL1(_mesh,getArray(),res);
  return res[compId];
}

/*!
 * Returns the \c normL1 of values of each component of \a this field:
 * \f[
 * \frac{\sum_{0 \leq i < nbOfEntity}|val[i]*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}
 * \f]
 *  \param [out] res - pointer to an array of result values, of size at least \a
 *         this->getNumberOfComponents(), that is to be allocated by the caller.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
void MEDCouplingFieldDouble::normL1(double *res) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL1");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform normL1 !");
  _type->normL1(_mesh,getArray(),res);
}

/*!
 * Returns the \c normL2 of values of a given component of \a this field:
 * \f[
 * \sqrt{\frac{\sum_{0 \leq i < nbOfEntity}|val[i]^{2}*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}}
 * \f]
 *  \param [in] compId - an index of the component of interest.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If \a compId is not valid.
           A valid range is ( 0 <= \a compId < \a this->getNumberOfComponents() ).
 */
double MEDCouplingFieldDouble::normL2(int compId) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL2");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform normL2 !");
  int nbComps=getArray()->getNumberOfComponents();
  if(compId<0 || compId>=nbComps)
    {
      std::ostringstream oss; oss << "MEDCouplingFieldDouble::normL2 : Invalid compId specified : No such nb of components ! Should be in [0," << nbComps << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  INTERP_KERNEL::AutoPtr<double> res=new double[nbComps];
  _type->normL2(_mesh,getArray(),res);
  return res[compId];
}

/*!
 * Returns the \c normL2 of values of each component of \a this field:
 * \f[
 * \sqrt{\frac{\sum_{0 \leq i < nbOfEntity}|val[i]^{2}*Vol[i]|}{\sum_{0 \leq i < nbOfEntity}|Vol[i]|}}
 * \f]
 *  \param [out] res - pointer to an array of result values, of size at least \a
 *         this->getNumberOfComponents(), that is to be allocated by the caller.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
void MEDCouplingFieldDouble::normL2(double *res) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform normL2");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform normL2 !");
  _type->normL2(_mesh,getArray(),res);
}

/*!
 * Computes a sum of values of a given component of \a this field multiplied by
 * values returned by buildMeasureField().
 * This method is useful to check the conservativity of interpolation method.
 *  \param [in] compId - an index of the component of interest.
 *  \param [in] isWAbs - if \c true (default), \c abs() is applied to the weighs computed by
 *         buildMeasureField() that makes this method slower. If a user is sure that all
 *         cells of the underlying mesh have correct orientation, he can put \a isWAbs ==
 *         \c false that speeds up this method.
 *  \throw If the mesh is not set.
 *  \throw If the data array is not set.
 *  \throw If \a compId is not valid.
           A valid range is ( 0 <= \a compId < \a this->getNumberOfComponents() ).
 */
double MEDCouplingFieldDouble::integral(int compId, bool isWAbs) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform integral");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform integral !");
  int nbComps=getArray()->getNumberOfComponents();
  if(compId<0 || compId>=nbComps)
    {
      std::ostringstream oss; oss << "MEDCouplingFieldDouble::integral : Invalid compId specified : No such nb of components ! Should be in [0," << nbComps << ") !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  INTERP_KERNEL::AutoPtr<double> res=new double[nbComps];
  _type->integral(_mesh,getArray(),isWAbs,res);
  return res[compId];
}

/*!
 * Computes a sum of values of each component of \a this field multiplied by
 * values returned by buildMeasureField().
 * This method is useful to check the conservativity of interpolation method.
 *  \param [in] isWAbs - if \c true (default), \c abs() is applied to the weighs computed by
 *         buildMeasureField() that makes this method slower. If a user is sure that all
 *         cells of the underlying mesh have correct orientation, he can put \a isWAbs ==
 *         \c false that speeds up this method.
 *  \param [out] res - pointer to an array of result sum values, of size at least \a
 *         this->getNumberOfComponents(), that is to be allocated by the caller.
 *  \throw If the mesh is not set.
 *  \throw If the data array is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
void MEDCouplingFieldDouble::integral(bool isWAbs, double *res) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform integral2");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform integral2 !");
  _type->integral(_mesh,getArray(),isWAbs,res);
}

/*!
 * Returns a value at a given cell of a structured mesh. The cell is specified by its
 * (i,j,k) index.
 *  \param [in] i - a index of node coordinates array along X axis. The cell is
 *         located between the i-th and ( i + 1 )-th nodes along X axis.
 *  \param [in] j - a index of node coordinates array along Y axis. The cell is
 *         located between the j-th and ( j + 1 )-th nodes along Y axis.
 *  \param [in] k - a index of node coordinates array along Z axis. The cell is
 *         located between the k-th and ( k + 1 )-th nodes along Z axis.
 *  \param [out] res - pointer to an array returning a feild value, of size at least
 *         \a this->getNumberOfComponents(), that is to be allocated by the caller.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 *  \throw If the mesh is not a structured one.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_getValueOnPos "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_getValueOnPos "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::getValueOnPos(int i, int j, int k, double *res) const
{
  const DataArrayDouble *arr=timeDiscr()->getArray();
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform getValueOnPos");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getValueOnPos !");
  _type->getValueOnPos(arr,_mesh,i,j,k,res);
}

/*!
 * Returns a value of \a this at a given point using spatial discretization.
 *  \param [in] spaceLoc - the point of interest.
 *  \param [out] res - pointer to an array returning a feild value, of size at least
 *         \a this->getNumberOfComponents(), that is to be allocated by the caller.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 *  \throw If \a spaceLoc is out of the spatial discretization.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_getValueOn "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_getValueOn "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double *res) const
{
  const DataArrayDouble *arr=timeDiscr()->getArray();
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform getValueOn");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getValueOnPos !");
  _type->getValueOn(arr,_mesh,spaceLoc,res);
}

/*!
 * Returns values of \a this at given points using spatial discretization.
 *  \param [in] spaceLoc - coordinates of points of interest in full-interlace
 *          mode. This array is to be of size ( \a nbOfPoints * \a this->getNumberOfComponents() ).
 *  \param [in] nbOfPoints - number of points of interest.
 *  \return DataArrayDouble * - a new instance of DataArrayDouble holding field
 *         values relating to the input points. This array is of size \a nbOfPoints
 *         tuples per \a this->getNumberOfComponents() components. The caller is to 
 *         delete this array using decrRef() as it is no more needed.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 *  \throw If any point in \a spaceLoc is out of the spatial discretization.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_getValueOnMulti "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_getValueOnMulti "Here is a Python example".
 *  \endif
 */
DataArrayDouble *MEDCouplingFieldDouble::getValueOnMulti(const double *spaceLoc, int nbOfPoints) const
{
  const DataArrayDouble *arr=timeDiscr()->getArray();
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform getValueOnMulti");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getValueOnMulti !");
  return _type->getValueOnMulti(arr,_mesh,spaceLoc,nbOfPoints);
}

/*!
 * Returns a value of \a this field at a given point at a given time using spatial discretization.
 * If the time is not covered by \a this->_time_discr, an exception is thrown.
 *  \param [in] spaceLoc - the point of interest.
 *  \param [in] time - the time of interest.
 *  \param [out] res - pointer to an array returning a feild value, of size at least
 *         \a this->getNumberOfComponents(), that is to be allocated by the caller.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 *  \throw If \a spaceLoc is out of the spatial discretization.
 *  \throw If \a time is not covered by \a this->_time_discr.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_getValueOn_time "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_getValueOn_time "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double time, double *res) const
{
  std::vector< const DataArrayDouble *> arrs=timeDiscr()->getArraysForTime(time);
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform getValueOn");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getValueOn !");
  std::vector<double> res2;
  for(std::vector< const DataArrayDouble *>::const_iterator iter=arrs.begin();iter!=arrs.end();iter++)
    {
      int sz=(int)res2.size();
      res2.resize(sz+(*iter)->getNumberOfComponents());
      _type->getValueOn(*iter,_mesh,spaceLoc,&res2[sz]);
    }
  timeDiscr()->getValueForTime(time,res2,res);
}

/*!
 * Apply a linear function to a given component of \a this field, so that
 * a component value <em>(x)</em> becomes \f$ a * x + b \f$.
 *  \param [in] a - the first coefficient of the function.
 *  \param [in] b - the second coefficient of the function.
 *  \param [in] compoId - the index of component to modify.
 *  \throw If the data array(s) is(are) not set.
 */
void MEDCouplingFieldDouble::applyLin(double a, double b, int compoId)
{
  timeDiscr()->applyLin(a,b,compoId);
}

/*!
 * Apply a linear function to all components of \a this field, so that
 * values <em>(x)</em> becomes \f$ a * x + b \f$.
 *  \param [in] a - the first coefficient of the function.
 *  \param [in] b - the second coefficient of the function.
 *  \throw If the data array(s) is(are) not set.
 */
void MEDCouplingFieldDouble::applyLin(double a, double b)
{
  timeDiscr()->applyLin(a,b);
}

/*!
 * This method sets \a this to a uniform scalar field with one component.
 * All tuples will have the same value 'value'.
 * An exception is thrown if no underlying mesh is defined.
 */
MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator=(double value)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::operator= : no mesh defined !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform operator = !");
  int nbOfTuple=_type->getNumberOfTuples(_mesh);
  timeDiscr()->setOrCreateUniformValueOnAllComponents(nbOfTuple,value);
  return *this;
}

/*!
 * Creates data array(s) of \a this field by using a C function for value generation.
 *  \param [in] nbOfComp - the number of components for \a this field to have.
 *  \param [in] func - the function used to compute values of \a this field.
 *         This function is to compute a field value basing on coordinates of value
 *         location point.
 *  \throw If the mesh is not set.
 *  \throw If \a func returns \c false.
 *  \throw If the spatial discretization of \a this field is NULL.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_fillFromAnalytic_c_func "Here is a C++ example".
 *  \endif
 */
void MEDCouplingFieldDouble::fillFromAnalytic(int nbOfComp, FunctionToEvaluate func)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::fillFromAnalytic : no mesh defined !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform fillFromAnalytic !");
  MCAuto<DataArrayDouble> loc=_type->getLocalizationOfDiscValues(_mesh);
  timeDiscr()->fillFromAnalytic(loc,nbOfComp,func);
}

/*!
 * Creates data array(s) of \a this field by using a function for value generation.<br>
 * The function is applied to coordinates of value location points. For example, if
 * \a this field is on cells, the function is applied to cell barycenters.
 * For more info on supported expressions that can be used in the function, see \ref
 * MEDCouplingArrayApplyFuncExpr. <br>
 * The function can include arbitrary named variables
 * (e.g. "x","y" or "va44") to refer to components of point coordinates. Names of
 * variables are sorted in \b alphabetical \b order to associate a variable name with a
 * component. For example, in the expression "2*x+z", "x" stands for the component #0
 * and "z" stands for the component #1 (\b not #2)!<br>
 * In a general case, a value resulting from the function evaluation is assigned to all
 * components of a field value. But there is a possibility to have its own expression for
 * each component within one function. For this purpose, there are predefined variable
 * names (IVec, JVec, KVec, LVec etc) each dedicated to a certain component (IVec, to
 * the component #0 etc). A factor of such a variable is added to the
 * corresponding component only.<br>
 * For example, \a nbOfComp == 4, coordinates of a 3D point are (1.,3.,7.), then
 *   - "2*x + z"               produces (5.,5.,5.,5.)
 *   - "2*x + 0*y + z"         produces (9.,9.,9.,9.)
 *   - "2*x*IVec + (x+z)*LVec" produces (2.,0.,0.,4.)
 *   - "2*y*IVec + z*KVec + x" produces (7.,1.,1.,4.)
 *
 *  \param [in] nbOfComp - the number of components for \a this field to have.
 *  \param [in] func - the function used to compute values of \a this field.
 *         This function is used to compute a field value basing on coordinates of value
 *         location point. For example, if \a this field is on cells, the function
 *         is applied to cell barycenters.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_fillFromAnalytic "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_fillFromAnalytic "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::fillFromAnalytic(int nbOfComp, const std::string& func)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::fillFromAnalytic : no mesh defined !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform fillFromAnalytic !");
  MCAuto<DataArrayDouble> loc=_type->getLocalizationOfDiscValues(_mesh);
  timeDiscr()->fillFromAnalytic(loc,nbOfComp,func);
}

/*!
 * Creates data array(s) of \a this field by using a function for value generation.<br>
 * The function is applied to coordinates of value location points. For example, if
 * \a this field is on cells, the function is applied to cell barycenters.<br>
 * This method differs from
 * \ref MEDCoupling::MEDCouplingFieldDouble::fillFromAnalytic(int nbOfComp, const std::string& func) "fillFromAnalytic()"
 * by the way how variable
 * names, used in the function, are associated with components of coordinates of field
 * location points; here, a variable name corresponding to a component is retrieved from
 * a corresponding node coordinates array (where it is set via
 * DataArrayDouble::setInfoOnComponent()).<br>
 * For more info on supported expressions that can be used in the function, see \ref
 * MEDCouplingArrayApplyFuncExpr. <br> 
 * In a general case, a value resulting from the function evaluation is assigned to all
 * components of a field value. But there is a possibility to have its own expression for
 * each component within one function. For this purpose, there are predefined variable
 * names (IVec, JVec, KVec, LVec etc) each dedicated to a certain component (IVec, to
 * the component #0 etc). A factor of such a variable is added to the
 * corresponding component only.<br>
 * For example, \a nbOfComp == 4, names of spatial components are "x", "y" and "z",
 * coordinates of a 3D point are (1.,3.,7.), then
 *   - "2*x + z"               produces (9.,9.,9.,9.)
 *   - "2*x*IVec + (x+z)*LVec" produces (2.,0.,0.,8.)
 *   - "2*y*IVec + z*KVec + x" produces (7.,1.,1.,8.)
 *
 *  \param [in] nbOfComp - the number of components for \a this field to have.
 *  \param [in] func - the function used to compute values of \a this field.
 *         This function is used to compute a field value basing on coordinates of value
 *         location point. For example, if \a this field is on cells, the function
 *         is applied to cell barycenters.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_fillFromAnalytic2 "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_fillFromAnalytic2 "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::fillFromAnalyticCompo(int nbOfComp, const std::string& func)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::fillFromAnalyticCompo : no mesh defined !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform fillFromAnalyticCompo !");
  MCAuto<DataArrayDouble> loc=_type->getLocalizationOfDiscValues(_mesh);
  timeDiscr()->fillFromAnalyticCompo(loc,nbOfComp,func);
}

/*!
 * Creates data array(s) of \a this field by using a function for value generation.<br>
 * The function is applied to coordinates of value location points. For example, if
 * \a this field is on cells, the function is applied to cell barycenters.<br>
 * This method differs from
 * \ref MEDCoupling::MEDCouplingFieldDouble::fillFromAnalytic(int nbOfComp, const std::string& func) "fillFromAnalytic()"
 * by the way how variable
 * names, used in the function, are associated with components of coordinates of field
 * location points; here, a component index of a variable is defined by a
 * rank of the variable within the input array \a varsOrder.<br>
 * For more info on supported expressions that can be used in the function, see \ref
 * MEDCouplingArrayApplyFuncExpr.
 * In a general case, a value resulting from the function evaluation is assigned to all
 * components of a field value. But there is a possibility to have its own expression for
 * each component within one function. For this purpose, there are predefined variable
 * names (IVec, JVec, KVec, LVec etc) each dedicated to a certain component (IVec, to
 * the component #0 etc). A factor of such a variable is added to the
 * corresponding component only.<br>
 * For example, \a nbOfComp == 4, names of
 * spatial components are given in \a varsOrder: ["x", "y","z"], coordinates of a
 * 3D point are (1.,3.,7.), then
 *   - "2*x + z"               produces (9.,9.,9.,9.)
 *   - "2*x*IVec + (x+z)*LVec" produces (2.,0.,0.,8.)
 *   - "2*y*IVec + z*KVec + x" produces (7.,1.,1.,8.)
 *
 *  \param [in] nbOfComp - the number of components for \a this field to have.
 *  \param [in] func - the function used to compute values of \a this field.
 *         This function is used to compute a field value basing on coordinates of value
 *         location point. For example, if \a this field is on cells, the function
 *         is applied to cell barycenters.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_fillFromAnalytic3 "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_fillFromAnalytic3 "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::fillFromAnalyticNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::fillFromAnalyticCompo : no mesh defined !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform fillFromAnalyticNamedCompo !");
  MCAuto<DataArrayDouble> loc=_type->getLocalizationOfDiscValues(_mesh);
  timeDiscr()->fillFromAnalyticNamedCompo(loc,nbOfComp,varsOrder,func);
}

/*!
 * Modifies values of \a this field by applying a C function to each tuple of all
 * data arrays.
 *  \param [in] nbOfComp - the number of components for \a this field to have.
 *  \param [in] func - the function used to compute values of \a this field.
 *         This function is to compute a field value basing on a current field value.
 *  \throw If \a func returns \c false.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_applyFunc_c_func "Here is a C++ example".
 *  \endif
 */
void MEDCouplingFieldDouble::applyFunc(int nbOfComp, FunctionToEvaluate func)
{
  timeDiscr()->applyFunc(nbOfComp,func);
}

/*!
 * Fill \a this field with a given value.<br>
 * This method is a specialization of other overloaded methods. When \a nbOfComp == 1
 * this method is equivalent to MEDCoupling::MEDCouplingFieldDouble::operator=().
 *  \param [in] nbOfComp - the number of components for \a this field to have.
 *  \param [in] val - the value to assign to every atomic value of \a this field.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the mesh is not set.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_applyFunc_val "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_applyFunc_val "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::applyFunc(int nbOfComp, double val)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::applyFunc : no mesh defined !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform applyFunc !");
  int nbOfTuple=_type->getNumberOfTuples(_mesh);
  timeDiscr()->setUniformValue(nbOfTuple,nbOfComp,val);
}

/*!
 * Modifies values of \a this field by applying a function to each tuple of all
 * data arrays.
 * For more info on supported expressions that can be used in the function, see \ref
 * MEDCouplingArrayApplyFuncExpr. <br>
 * The function can include arbitrary named variables
 * (e.g. "x","y" or "va44") to refer to components of a field value. Names of
 * variables are sorted in \b alphabetical \b order to associate a variable name with a
 * component. For example, in the expression "2*x+z", "x" stands for the component #0
 * and "z" stands for the component #1 (\b not #2)!<br>
 * In a general case, a value resulting from the function evaluation is assigned to all
 * components of a field value. But there is a possibility to have its own expression for
 * each component within one function. For this purpose, there are predefined variable
 * names (IVec, JVec, KVec, LVec etc) each dedicated to a certain component (IVec, to
 * the component #0 etc). A factor of such a variable is added to the
 * corresponding component only.<br>
 * For example, \a nbOfComp == 4, components of a field value are (1.,3.,7.), then
 *   - "2*x + z"               produces (5.,5.,5.,5.)
 *   - "2*x + 0*y + z"         produces (9.,9.,9.,9.)
 *   - "2*x*IVec + (x+z)*LVec" produces (2.,0.,0.,4.)
 *   - "2*y*IVec + z*KVec + x" produces (7.,1.,1.,4.)
 *
 *  \param [in] nbOfComp - the number of components for \a this field to have.
 *  \param [in] func - the function used to compute values of \a this field.
 *         This function is to compute a field value basing on a current field value.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_applyFunc "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_applyFunc "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::applyFunc(int nbOfComp, const std::string& func)
{
  timeDiscr()->applyFunc(nbOfComp,func);
}


/*!
 * Modifies values of \a this field by applying a function to each tuple of all
 * data arrays.
 * For more info on supported expressions that can be used in the function, see \ref
 * MEDCouplingArrayApplyFuncExpr. <br>
 * This method differs from
 * \ref MEDCoupling::MEDCouplingFieldDouble::applyFunc(int nbOfComp, const std::string& func) "applyFunc()"
 * by the way how variable
 * names, used in the function, are associated with components of field values;
 * here, a variable name corresponding to a component is retrieved from
 * component information of an array (where it is set via
 * DataArrayDouble::setInfoOnComponent()).<br>
 * In a general case, a value resulting from the function evaluation is assigned to all
 * components of a field value. But there is a possibility to have its own expression for
 * each component within one function. For this purpose, there are predefined variable
 * names (IVec, JVec, KVec, LVec etc) each dedicated to a certain component (IVec, to
 * the component #0 etc). A factor of such a variable is added to the
 * corresponding component only.<br>
 * For example, \a nbOfComp == 4, components of a field value are (1.,3.,7.), then
 *   - "2*x + z"               produces (5.,5.,5.,5.)
 *   - "2*x + 0*y + z"         produces (9.,9.,9.,9.)
 *   - "2*x*IVec + (x+z)*LVec" produces (2.,0.,0.,4.)
 *   - "2*y*IVec + z*KVec + x" produces (7.,1.,1.,4.)
 *
 *  \param [in] nbOfComp - the number of components for \a this field to have.
 *  \param [in] func - the function used to compute values of \a this field.
 *         This function is to compute a new field value basing on a current field value.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_applyFunc2 "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_applyFunc2 "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::applyFuncCompo(int nbOfComp, const std::string& func)
{
  timeDiscr()->applyFuncCompo(nbOfComp,func);
}

/*!
 * Modifies values of \a this field by applying a function to each tuple of all
 * data arrays.
 * This method differs from
 * \ref MEDCoupling::MEDCouplingFieldDouble::applyFunc(int nbOfComp, const std::string& func) "applyFunc()"
 * by the way how variable
 * names, used in the function, are associated with components of field values;
 * here, a component index of a variable is defined by a
 * rank of the variable within the input array \a varsOrder.<br>
 * For more info on supported expressions that can be used in the function, see \ref
 * MEDCouplingArrayApplyFuncExpr.
 * In a general case, a value resulting from the function evaluation is assigned to all
 * components of a field value. But there is a possibility to have its own expression for
 * each component within one function. For this purpose, there are predefined variable
 * names (IVec, JVec, KVec, LVec etc) each dedicated to a certain component (IVec, to
 * the component #0 etc). A factor of such a variable is added to the
 * corresponding component only.<br>
 * For example, \a nbOfComp == 4, names of
 * components are given in \a varsOrder: ["x", "y","z"], components of a
 * 3D vector are (1.,3.,7.), then
 *   - "2*x + z"               produces (9.,9.,9.,9.)
 *   - "2*x*IVec + (x+z)*LVec" produces (2.,0.,0.,8.)
 *   - "2*y*IVec + z*KVec + x" produces (7.,1.,1.,8.)
 *
 *  \param [in] nbOfComp - the number of components for \a this field to have.
 *  \param [in] func - the function used to compute values of \a this field.
 *         This function is to compute a new field value basing on a current field value.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_applyFunc3 "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_applyFunc3 "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::applyFuncNamedCompo(int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func)
{
  timeDiscr()->applyFuncNamedCompo(nbOfComp,varsOrder,func);
}

/*!
 * Modifies values of \a this field by applying a function to each atomic value of all
 * data arrays. The function computes a new single value basing on an old single value.
 * For more info on supported expressions that can be used in the function, see \ref
 * MEDCouplingArrayApplyFuncExpr. <br>
 * The function can include **only one** arbitrary named variable
 * (e.g. "x","y" or "va44") to refer to a field atomic value. <br>
 * In a general case, a value resulting from the function evaluation is assigned to 
 * a single field value. But there is a possibility to have its own expression for
 * each component within one function. For this purpose, there are predefined variable
 * names (IVec, JVec, KVec, LVec etc) each dedicated to a certain component (IVec, to
 * the component #0 etc). A factor of such a variable is added to the
 * corresponding component only.<br>
 * For example, components of a field value are (1.,3.,7.), then
 *   - "2*x - 1"               produces (1.,5.,13.)
 *   - "2*x*IVec + (x+3)*KVec" produces (2.,0.,10.)
 *   - "2*x*IVec + (x+3)*KVec + 1" produces (3.,1.,11.)
 *
 *  \param [in] func - the function used to compute values of \a this field.
 *         This function is to compute a field value basing on a current field value.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_applyFunc_same_nb_comp "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_applyFunc_same_nb_comp "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::applyFunc(const std::string& func)
{
  timeDiscr()->applyFunc(func);
}

/*!
 * Applyies the function specified by the string repr 'func' on each tuples on all arrays contained in _time_discr.
 * The field will contain exactly the same number of components after the call.
 * Use is not warranted for the moment !
 */
void MEDCouplingFieldDouble::applyFuncFast32(const std::string& func)
{
  timeDiscr()->applyFuncFast32(func);
}

/*!
 * Applyies the function specified by the string repr 'func' on each tuples on all arrays contained in _time_discr.
 * The field will contain exactly the same number of components after the call.
 * Use is not warranted for the moment !
 */
void MEDCouplingFieldDouble::applyFuncFast64(const std::string& func)
{
  timeDiscr()->applyFuncFast64(func);
}

/*!
 * Returns number of components in the data array. For more info on the data arrays,
 * see \ref arrays.
 *  \return int - the number of components in the data array.
 *  \throw If the data array is not set.
 */
std::size_t MEDCouplingFieldDouble::getNumberOfComponents() const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::getNumberOfComponents : No array specified !");
  return getArray()->getNumberOfComponents();
}

/*!
 * Use MEDCouplingField::getNumberOfTuplesExpected instead of this method. This method will be removed soon, because it is
 * confusing compared to getNumberOfComponents() and getNumberOfValues() behaviour.
 *
 * Returns number of tuples in \a this field, that depends on 
 * - the number of entities in the underlying mesh
 * - \ref MEDCouplingSpatialDisc "spatial discretization" of \a this field (e.g. number
 * of Gauss points if \a this->getTypeOfField() == 
 * \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT").
 *
 * The returned value does \b not \b depend on the number of tuples in the data array
 * (which has to be equal to the returned value), \b contrary to
 * getNumberOfComponents() and getNumberOfValues() that retrieve information from the
 * data array (Sorry, it is confusing !).
 * So \b this \b method \b behaves \b exactly \b as MEDCouplingField::getNumberOfTuplesExpected \b method.
 *
 * \warning No checkConsistencyLight() is done here.
 * For more info on the data arrays, see \ref arrays.
 *  \return int - the number of tuples.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the spatial discretization is not fully defined.
 *  \sa MEDCouplingField::getNumberOfTuplesExpected
 */
std::size_t MEDCouplingFieldDouble::getNumberOfTuples() const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Impossible to retrieve number of tuples because no mesh specified !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getNumberOfTuples !");
  return _type->getNumberOfTuples(_mesh);
}

/*!
 * Returns number of atomic double values in the data array of \a this field.
 * For more info on the data arrays, see \ref arrays.
 *  \return int - (number of tuples) * (number of components) of the
 *  data array.
 *  \throw If the data array is not set.
 */
std::size_t MEDCouplingFieldDouble::getNumberOfValues() const
{
  if(getArray()==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::getNumberOfValues : No array specified !");
  return getArray()->getNbOfElems();
}

/*!
 * Sets own modification time by the most recently modified element of data (the mesh,
 * the data array etc). For more info, see \ref MEDCouplingTimeLabelPage.
 */
void MEDCouplingFieldDouble::updateTime() const
{
  MEDCouplingField::updateTime();
  updateTimeWith(*timeDiscr());
}

std::size_t MEDCouplingFieldDouble::getHeapMemorySizeWithoutChildren() const
{
  return MEDCouplingField::getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDCouplingFieldDouble::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDCouplingField::getDirectChildrenWithNull());
  if(timeDiscr())
    {
      std::vector<const BigMemoryObject *> ret2(timeDiscr()->getDirectChildrenWithNull());
      ret.insert(ret.end(),ret2.begin(),ret2.end());
    }
  return ret;
}

/*!
 * Returns a value of \a this field of type either
 * \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT" or
 * \ref MEDCoupling::ON_GAUSS_NE "ON_GAUSS_NE".
 *  \param [in] cellId - an id of cell of interest.
 *  \param [in] nodeIdInCell - a node index within the cell.
 *  \param [in] compoId - an index of component.
 *  \return double - the field value corresponding to the specified parameters.
 *  \throw If the data array is not set.
 *  \throw If the mesh is not set.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If \a this field if of type other than 
 *         \ref MEDCoupling::ON_GAUSS_PT "ON_GAUSS_PT" or
 *         \ref MEDCoupling::ON_GAUSS_NE "ON_GAUSS_NE".
 */
double MEDCouplingFieldDouble::getIJK(int cellId, int nodeIdInCell, int compoId) const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform getIJK !");
  return _type->getIJK(_mesh,getArray(),cellId,nodeIdInCell,compoId);
}

/*!
 * Sets the data array. 
 *  \param [in] array - the data array holding values of \a this field. It's size
 *         should correspond to the mesh and
 *         \ref MEDCouplingSpatialDisc "spatial discretization" of \a this field
 *         (see getNumberOfTuples()), but this size is not checked here.
 */
//void MEDCouplingFieldDouble::setArray(DataArrayDouble *array)

/*!
 * Sets the data array holding values corresponding to an end of a time interval
 * for which \a this field is defined.
 *  \param [in] array - the data array holding values of \a this field. It's size
 *         should correspond to the mesh and
 *         \ref MEDCouplingSpatialDisc "spatial discretization" of \a this field
 *         (see getNumberOfTuples()), but this size is not checked here.
 */
//void MEDCouplingFieldDouble::setEndArray(DataArrayDouble *array)

/*!
 * Sets all data arrays needed to define the field values.
 *  \param [in] arrs - a vector of DataArrayDouble's holding values of \a this
 *         field. Size of each array should correspond to the mesh and
 *         \ref MEDCouplingSpatialDisc "spatial discretization" of \a this field
 *         (see getNumberOfTuples()), but this size is not checked here.
 *  \throw If number of arrays in \a arrs does not correspond to type of
 *         \ref MEDCouplingTemporalDisc "temporal discretization" of \a this field.
 */
//void MEDCouplingFieldDouble::setArrays(const std::vector<DataArrayDouble *>& arrs)

/*!
 * Tries to set an \a other mesh as the support of \a this field. An attempt fails, if
 * a current and the \a other meshes are different with use of specified equality
 * criteria, and then an exception is thrown.
 *  \param [in] other - the mesh to use as the field support if this mesh can be
 *         considered equal to the current mesh.
 *  \param [in] levOfCheck - defines equality criteria used for mesh comparison. For
 *         it's meaning explanation, see MEDCouplingMesh::checkGeoEquivalWith() which
 *         is used for mesh comparison.
 *  \param [in] precOnMesh - a precision used to compare nodes of the two meshes.
 *         It is used as \a prec parameter of MEDCouplingMesh::checkGeoEquivalWith().
 *  \param [in] eps - a precision used at node renumbering (if needed) to compare field
 *         values at merged nodes. If the values differ more than \a eps, an
 *         exception is thrown.
 *  \throw If the mesh is not set.
 *  \throw If \a other == NULL.
 *  \throw If any of the meshes is not well defined.
 *  \throw If the two meshes do not match.
 *  \throw If field values at merged nodes (if any) deffer more than \a eps.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_changeUnderlyingMesh "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_changeUnderlyingMesh "Here is a Python example".
 *  \endif
 */
void MEDCouplingFieldDouble::changeUnderlyingMesh(const MEDCouplingMesh *other, int levOfCheck, double precOnMesh, double eps)
{
  if(_mesh==0 || other==0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::changeUnderlyingMesh : is expected to operate on not null meshes !");
  DataArrayInt *cellCor=0,*nodeCor=0;
  other->checkGeoEquivalWith(_mesh,levOfCheck,precOnMesh,cellCor,nodeCor);
  MCAuto<DataArrayInt> cellCor2(cellCor),nodeCor2(nodeCor);
  if(cellCor)
    renumberCellsWithoutMesh(cellCor->getConstPointer(),false);
  if(nodeCor)
    renumberNodesWithoutMesh(nodeCor->getConstPointer(),nodeCor->getMaxValueInArray()+1,eps);
  setMesh(other);
}

/*!
 * Subtracts another field from \a this one in case when the two fields have different
 * supporting meshes. The subtraction is performed provided that the tho meshes can be
 * considered equal with use of specified equality criteria, else an exception is thrown.
 * If the meshes match, the mesh of \a f is set to \a this field (\a this is permuted if 
 * necessary) and field values are subtracted. No interpolation is done here, only an
 * analysis of two underlying mesh is done to see if the meshes are geometrically
 * equivalent.<br>
 * The job of this method consists in calling
 * \a this->changeUnderlyingMesh() with \a f->getMesh() as the first parameter, and then
 * \a this -= \a f.<br>
 * This method requires that \a f and \a this are coherent (checkConsistencyLight()) and that \a f
 * and \a this are coherent for a merge.<br>
 * "DM" in the method name stands for "different meshes".
 *  \param [in] f - the field to subtract from this.
 *  \param [in] levOfCheck - defines equality criteria used for mesh comparison. For
 *         it's meaning explanation, see MEDCouplingMesh::checkGeoEquivalWith() which
 *         is used for mesh comparison.
 *  \param [in] precOnMesh - a precision used to compare nodes of the two meshes.
 *         It is used as \a prec parameter of MEDCouplingMesh::checkGeoEquivalWith().
 *  \param [in] eps - a precision used at node renumbering (if needed) to compare field
 *         values at merged nodes. If the values differ more than \a eps, an
 *         exception is thrown.
 *  \throw If \a f == NULL.
 *  \throw If any of the meshes is not set or is not well defined.
 *  \throw If the two meshes do not match.
 *  \throw If the two fields are not coherent for merge.
 *  \throw If field values at merged nodes (if any) deffer more than \a eps.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_substractInPlaceDM "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_substractInPlaceDM "Here is a Python example".
 *  \endif
 *  \sa changeUnderlyingMesh().
 */
void MEDCouplingFieldDouble::substractInPlaceDM(const MEDCouplingFieldDouble *f, int levOfCheck, double precOnMesh, double eps)
{
  checkConsistencyLight();
  if(!f)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::substractInPlaceDM : input field is NULL !");
  f->checkConsistencyLight();
  if(!areCompatibleForMerge(f))
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::substractInPlaceDM : Fields are not compatible ; unable to apply mergeFields on them !");
  changeUnderlyingMesh(f->getMesh(),levOfCheck,precOnMesh,eps);
  operator-=(*f);
}

/*!
 * Merges coincident nodes of the underlying mesh. If some nodes are coincident, the
 * underlying mesh is replaced by a new mesh instance where the coincident nodes are merged.
 *  \param [in] eps - a precision used to compare nodes of the two meshes.
 *  \param [in] epsOnVals - a precision used to compare field
 *         values at merged nodes. If the values differ more than \a epsOnVals, an
 *         exception is thrown.
 *  \return bool - \c true if some nodes have been merged and hence \a this field lies
 *         on another mesh.
 *  \throw If the mesh is of type not inheriting from MEDCouplingPointSet.
 *  \throw If the mesh is not well defined.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the data array is not set.
 *  \throw If field values at merged nodes (if any) deffer more than \a epsOnVals.
 */
bool MEDCouplingFieldDouble::mergeNodes(double eps, double epsOnVals)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("Invalid support mesh to apply mergeNodes on it : must be a MEDCouplingPointSet one !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform mergeNodes !");
  MCAuto<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCopy());
  bool ret;
  int ret2;
  MCAuto<DataArrayInt> arr=meshC2->mergeNodes(eps,ret,ret2);
  if(!ret)//no nodes have been merged.
    return ret;
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    if(*iter)
      _type->renumberValuesOnNodes(epsOnVals,arr->getConstPointer(),meshC2->getNumberOfNodes(),*iter);
  setMesh(meshC2);
  return true;
}

/*!
 * Merges coincident nodes of the underlying mesh. If some nodes are coincident, the
 * underlying mesh is replaced by a new mesh instance where the coincident nodes are
 * merged.<br>
 * In contrast to mergeNodes(), location of merged nodes is changed to be at their barycenter.
 *  \param [in] eps - a precision used to compare nodes of the two meshes.
 *  \param [in] epsOnVals - a precision used to compare field
 *         values at merged nodes. If the values differ more than \a epsOnVals, an
 *         exception is thrown.
 *  \return bool - \c true if some nodes have been merged and hence \a this field lies
 *         on another mesh.
 *  \throw If the mesh is of type not inheriting from MEDCouplingPointSet.
 *  \throw If the mesh is not well defined.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the data array is not set.
 *  \throw If field values at merged nodes (if any) deffer more than \a epsOnVals.
 */
bool MEDCouplingFieldDouble::mergeNodesCenter(double eps, double epsOnVals)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("Invalid support mesh to apply mergeNodes on it : must be a MEDCouplingPointSet one !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform mergeNodesCenter !");
  MCAuto<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCopy());
  bool ret;
  int ret2;
  MCAuto<DataArrayInt> arr=meshC2->mergeNodesCenter(eps,ret,ret2);
  if(!ret)//no nodes have been merged.
    return ret;
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    if(*iter)
      _type->renumberValuesOnNodes(epsOnVals,arr->getConstPointer(),meshC2->getNumberOfNodes(),*iter);
  setMesh(meshC2);
  return true;
}

/*!
 * Removes from the underlying mesh nodes not used in any cell. If some nodes are
 * removed, the underlying mesh is replaced by a new mesh instance where the unused
 * nodes are removed.<br>
 *  \param [in] epsOnVals - a precision used to compare field
 *         values at merged nodes. If the values differ more than \a epsOnVals, an
 *         exception is thrown.
 *  \return bool - \c true if some nodes have been removed and hence \a this field lies
 *         on another mesh.
 *  \throw If the mesh is of type not inheriting from MEDCouplingPointSet.
 *  \throw If the mesh is not well defined.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the data array is not set.
 *  \throw If field values at merged nodes (if any) deffer more than \a epsOnVals.
 */
bool MEDCouplingFieldDouble::zipCoords(double epsOnVals)
{
  const MEDCouplingPointSet *meshC=dynamic_cast<const MEDCouplingPointSet *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::zipCoords : Invalid support mesh to apply zipCoords on it : must be a MEDCouplingPointSet one !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform zipCoords !");
  MCAuto<MEDCouplingPointSet> meshC2((MEDCouplingPointSet *)meshC->deepCopy());
  int oldNbOfNodes=meshC2->getNumberOfNodes();
  MCAuto<DataArrayInt> arr=meshC2->zipCoordsTraducer();
  if(meshC2->getNumberOfNodes()!=oldNbOfNodes)
    {
      std::vector<DataArrayDouble *> arrays;
      timeDiscr()->getArrays(arrays);
      for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
        if(*iter)
          _type->renumberValuesOnNodes(epsOnVals,arr->getConstPointer(),meshC2->getNumberOfNodes(),*iter);
      setMesh(meshC2);
      return true;
    }
  return false;
}

/*!
 * Removes duplicates of cells from the understanding mesh. If some cells are
 * removed, the underlying mesh is replaced by a new mesh instance where the cells
 * duplicates are removed.<br>
 *  \param [in] compType - specifies a cell comparison technique. Meaning of its
 *          valid values [0,1,2] is explained in the description of
 *          MEDCouplingPointSet::zipConnectivityTraducer() which is called by this method.
 *  \param [in] epsOnVals - a precision used to compare field
 *         values at merged cells. If the values differ more than \a epsOnVals, an
 *         exception is thrown.
 *  \return bool - \c true if some cells have been removed and hence \a this field lies
 *         on another mesh.
 *  \throw If the mesh is not an instance of MEDCouplingUMesh.
 *  \throw If the mesh is not well defined.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the data array is not set.
 *  \throw If field values at merged cells (if any) deffer more than \a epsOnVals.
 */
bool MEDCouplingFieldDouble::zipConnectivity(int compType, double epsOnVals)
{
  const MEDCouplingUMesh *meshC=dynamic_cast<const MEDCouplingUMesh *>(_mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::zipConnectivity : Invalid support mesh to apply zipCoords on it : must be a MEDCouplingPointSet one !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform zipConnectivity !");
  MCAuto<MEDCouplingUMesh> meshC2((MEDCouplingUMesh *)meshC->deepCopy());
  std::size_t oldNbOfCells(meshC2->getNumberOfCells());
  MCAuto<DataArrayInt> arr=meshC2->zipConnectivityTraducer(compType);
  if(meshC2->getNumberOfCells()!=oldNbOfCells)
    {
      std::vector<DataArrayDouble *> arrays;
      timeDiscr()->getArrays(arrays);
      for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
        if(*iter)
          _type->renumberValuesOnCells(epsOnVals,meshC,arr->getConstPointer(),meshC2->getNumberOfCells(),*iter);
      setMesh(meshC2);
      return true;
    }
  return false;
}

/*!
 * This method calls MEDCouplingUMesh::buildSlice3D method. So this method makes the assumption that underlying mesh exists.
 * For the moment, this method is implemented for fields on cells.
 * 
 * \return a newly allocated field double containing the result that the user should deallocate.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::extractSlice3D(const double *origin, const double *vec, double eps) const
{
  const MEDCouplingMesh *mesh=getMesh();
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::extractSlice3D : underlying mesh is null !");
  if(getTypeOfField()!=ON_CELLS)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::extractSlice3D : only implemented for fields on cells !");
  const MCAuto<MEDCouplingUMesh> umesh(mesh->buildUnstructured());
  MCAuto<MEDCouplingFieldDouble> ret(clone(false));
  ret->setMesh(umesh);
  DataArrayInt *cellIds=0;
  MCAuto<MEDCouplingUMesh> mesh2=umesh->buildSlice3D(origin,vec,eps,cellIds);
  MCAuto<DataArrayInt> cellIds2=cellIds;
  ret->setMesh(mesh2);
  MCAuto<DataArrayInt> tupleIds=computeTupleIdsToSelectFromCellIds(cellIds->begin(),cellIds->end());
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  int i=0;
  std::vector<DataArrayDouble *> newArr(arrays.size());
  std::vector< MCAuto<DataArrayDouble> > newArr2(arrays.size());
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++,i++)
    {
      if(*iter)
        {
          newArr2[i]=(*iter)->selectByTupleIdSafe(cellIds->begin(),cellIds->end());
          newArr[i]=newArr2[i];
        }
    }
  ret->setArrays(newArr);
  return ret.retn();
}

/*!
 * Divides every cell of the underlying mesh into simplices (triangles in 2D and
 * tetrahedra in 3D). If some cells are divided, the underlying mesh is replaced by a new
 * mesh instance containing the simplices.<br> 
 *  \param [in] policy - specifies a pattern used for splitting. For its description, see
 *          MEDCouplingUMesh::simplexize().
 *  \return bool - \c true if some cells have been divided and hence \a this field lies
 *         on another mesh.
 *  \throw If \a policy has an invalid value. For valid values, see the description of 
 *         MEDCouplingUMesh::simplexize().
 *  \throw If MEDCouplingMesh::simplexize() is not applicable to the underlying mesh.
 *  \throw If the mesh is not well defined.
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If the data array is not set.
 */
bool MEDCouplingFieldDouble::simplexize(int policy)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No underlying mesh on this field to perform simplexize !");
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform simplexize !");
  int oldNbOfCells=_mesh->getNumberOfCells();
  MCAuto<MEDCouplingMesh> meshC2(_mesh->deepCopy());
  MCAuto<DataArrayInt> arr=meshC2->simplexize(policy);
  int newNbOfCells=meshC2->getNumberOfCells();
  if(oldNbOfCells==newNbOfCells)
    return false;
  std::vector<DataArrayDouble *> arrays;
  timeDiscr()->getArrays(arrays);
  for(std::vector<DataArrayDouble *>::const_iterator iter=arrays.begin();iter!=arrays.end();iter++)
    if(*iter)
      _type->renumberValuesOnCellsR(_mesh,arr->getConstPointer(),arr->getNbOfElems(),*iter);
  setMesh(meshC2);
  return true;
}

/*!
 * This method makes the hypothesis that \a this is a Gauss field. This method returns a newly created field on cells with same number of tuples than \a this.
 * Each Gauss points in \a this is replaced by a polygon or polyhedron cell with associated region following Voronoi algorithm.
 */
MCAuto<MEDCouplingFieldDouble> MEDCouplingFieldDouble::voronoize(double eps) const
{
  checkConsistencyLight();
  const MEDCouplingMesh *mesh(getMesh());
  INTERP_KERNEL::AutoCppPtr<Voronizer> vor;
  int meshDim(mesh->getMeshDimension()),spaceDim(mesh->getSpaceDimension());
  if(meshDim==1 && (spaceDim==1 || spaceDim==2 || spaceDim==3))
    vor=new Voronizer1D;
  else if(meshDim==2 && (spaceDim==2 || spaceDim==3))
    vor=new Voronizer2D;
  else if(meshDim==3 && spaceDim==3)
    vor=new Voronizer3D;
  else
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::voronoize : only 2D, 3D surf, and 3D are supported for the moment !");
  return voronoizeGen(vor,eps);
}

/*!
 * \sa MEDCouplingUMesh::convertQuadraticCellsToLinear
 */
MCAuto<MEDCouplingFieldDouble> MEDCouplingFieldDouble::convertQuadraticCellsToLinear() const
{
  checkConsistencyLight();
  switch(getTypeOfField())
    {
    case ON_NODES:
      {
        const MEDCouplingMesh *mesh(getMesh());
        if(!mesh)
          throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::convertQuadraticCellsToLinear : null mesh !");
        MCAuto<MEDCouplingUMesh> umesh(mesh->buildUnstructured());
        umesh=umesh->clone(false);
        umesh->convertQuadraticCellsToLinear();
        MCAuto<DataArrayInt> o2n(umesh->zipCoordsTraducer());
        MCAuto<DataArrayInt> n2o(o2n->invertArrayO2N2N2O(umesh->getNumberOfNodes()));
        MCAuto<DataArrayDouble> arr(getArray()->selectByTupleIdSafe(n2o->begin(),n2o->end()));
        MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_NODES));
        ret->setArray(arr);
        ret->setMesh(umesh);
        ret->copyAllTinyAttrFrom(this);
        return ret;
      }
    case ON_CELLS:
      {
        const MEDCouplingMesh *mesh(getMesh());
        if(!mesh)
          throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::convertQuadraticCellsToLinear : null mesh !");
        MCAuto<MEDCouplingUMesh> umesh(mesh->buildUnstructured());
        umesh=umesh->clone(false);
        umesh->convertQuadraticCellsToLinear();
        umesh->zipCoords();
        MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_CELLS));
        ret->setArray(const_cast<DataArrayDouble *>(getArray()));
        ret->setMesh(umesh);
        ret->copyAllTinyAttrFrom(this);
        return ret;
      }
    case ON_GAUSS_PT:
      {
        const MEDCouplingMesh *mesh(getMesh());
        if(!mesh)
          throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::convertQuadraticCellsToLinear : null mesh !");
        MCAuto<MEDCouplingUMesh> umesh(mesh->buildUnstructured());
        std::set<INTERP_KERNEL::NormalizedCellType> gt(umesh->getAllGeoTypes());
        MCAuto<MEDCouplingFieldDouble> ret(MEDCouplingFieldDouble::New(ON_GAUSS_PT));
        //
        const MEDCouplingFieldDiscretization *disc(getDiscretization());
        const MEDCouplingFieldDiscretizationGauss *disc2(dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(disc));
        if(!disc2)
          throw INTERP_KERNEL::Exception("convertQuadraticCellsToLinear : Not a ON_GAUSS_PT field");
        std::set<INTERP_KERNEL::NormalizedCellType> gt2(umesh->getAllGeoTypes());
        std::vector< MCAuto<DataArrayInt> > cellIdsV;
        std::vector< MCAuto<MEDCouplingUMesh> > meshesV;
        std::vector< MEDCouplingGaussLocalization > glV;
        bool isZipReq(false);
        for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator it=gt.begin();it!=gt.end();it++)
          {
            const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel(*it));
            MCAuto<DataArrayInt> cellIds(umesh->giveCellsWithType(*it));
            cellIdsV.push_back(cellIds);
            MCAuto<MEDCouplingUMesh> part(umesh->buildPartOfMySelf(cellIds->begin(),cellIds->end()));
            int id(disc2->getGaussLocalizationIdOfOneType(*it));
            const MEDCouplingGaussLocalization& gl(disc2->getGaussLocalization(id));
            if(!cm.isQuadratic())
              {
                glV.push_back(gl);
              }
            else
              {
                isZipReq=true;
                part->convertQuadraticCellsToLinear();
                INTERP_KERNEL::GaussInfo gi(*it,gl.getGaussCoords(),gl.getNumberOfGaussPt(),gl.getRefCoords(),gl.getNumberOfPtsInRefCell());
                INTERP_KERNEL::GaussInfo gi2(gi.convertToLinear());
                MEDCouplingGaussLocalization gl2(gi2.getGeoType(),gi2.getRefCoords(),gi2.getGaussCoords(),gl.getWeights());
                glV.push_back(gl2);
              }
            meshesV.push_back(part);
          }
        //
        {
          std::vector< const MEDCouplingUMesh * > meshesPtr(VecAutoToVecOfCstPt(meshesV));
          umesh=MEDCouplingUMesh::MergeUMeshesOnSameCoords(meshesPtr);
          std::vector< const DataArrayInt * > zeCellIds(VecAutoToVecOfCstPt(cellIdsV));
          MCAuto<DataArrayInt> zeIds(DataArrayInt::Aggregate(zeCellIds));
          umesh->renumberCells(zeIds->begin());
          umesh->setName(mesh->getName());
        }
        //
        if(isZipReq)
          umesh->zipCoords();
        ret->setArray(const_cast<DataArrayDouble *>(getArray()));
        ret->setMesh(umesh);
        for(std::vector< MEDCouplingGaussLocalization >::const_iterator it=glV.begin();it!=glV.end();it++)
          ret->setGaussLocalizationOnType((*it).getType(),(*it).getRefCoords(),(*it).getGaussCoords(),(*it).getWeights());
        ret->copyAllTinyAttrFrom(this);
        ret->checkConsistencyLight();
        return ret;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::convertQuadraticCellsToLinear : Only available for fields on nodes and on cells !");
    }
}

/*!
 * This is expected to be a 3 components vector field on nodes (if not an exception will be thrown). \a this is also expected to lie on a MEDCouplingPointSet mesh.
 * Finaly \a this is also expected to be consistent.
 * In these conditions this method returns a newly created field (to be dealed by the caller).
 * The returned field will also 3 compo vector field be on nodes lying on the same mesh than \a this.
 * 
 * For each 3 compo tuple \a tup in \a this the returned tuple is the result of the transformation of \a tup in the new referential. This referential is defined by \a Ur, \a Uteta, \a Uz.
 * \a Ur is the vector between \a center point and the associated node with \a tuple. \a Uz is \a vect normalized. And Uteta is the cross product of \a Uz with \a Ur.
 *
 * \sa DataArrayDouble::fromCartToCylGiven
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::computeVectorFieldCyl(const double center[3], const double vect[3]) const
{
  checkConsistencyLight();
  const DataArrayDouble *coo(getMesh()->getDirectAccessOfCoordsArrIfInStructure());
  MEDCouplingTimeDiscretization *td(timeDiscr()->computeVectorFieldCyl(coo,center,vect));
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setMesh(getMesh());
  ret->setName(getName());
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble filled with the doubly contracted product of
 * every tensor of \a this 6-componental field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, whose
 *          each tuple is calculated from the tuple <em>(t)</em> of \a this field as
 *          follows: \f$ t[0]^2+t[1]^2+t[2]^2+2*t[3]^2+2*t[4]^2+2*t[5]^2\f$. 
 *          This new field lies on the same mesh as \a this one. The caller is to delete
 *          this field using decrRef() as it is no more needed.
 *  \throw If \a this->getNumberOfComponents() != 6.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::doublyContractedProduct() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform doublyContractedProduct !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->doublyContractedProduct());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setName("DoublyContractedProduct");
  ret->setMesh(getMesh());
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble filled with the determinant of a square
 * matrix defined by every tuple of \a this field, having either 4, 6 or 9 components.
 * The case of 6 components corresponds to that of the upper triangular matrix. 
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, whose
 *          each tuple is the determinant of matrix of the corresponding tuple of \a this 
 *          field. This new field lies on the same mesh as \a this one. The caller is to 
 *          delete this field using decrRef() as it is no more needed.
 *  \throw If \a this->getNumberOfComponents() is not in [4,6,9].
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::determinant() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform determinant !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->determinant());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setName("Determinant");
  ret->setMesh(getMesh());
  return ret.retn();
}


/*!
 * Creates a new MEDCouplingFieldDouble with 3 components filled with 3 eigenvalues of
 * an upper triangular matrix defined by every tuple of \a this 6-componental field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, 
 *          having 3 components, whose each tuple contains the eigenvalues of the matrix of
 *          corresponding tuple of \a this field. This new field lies on the same mesh as
 *          \a this one. The caller is to delete this field using decrRef() as it is no
 *          more needed.  
 *  \throw If \a this->getNumberOfComponents() != 6.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::eigenValues() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform eigenValues !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->eigenValues());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setName("EigenValues");
  ret->setMesh(getMesh());
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble with 9 components filled with 3 eigenvectors of
 * an upper triangular matrix defined by every tuple of \a this 6-componental field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, 
 *          having 9 components, whose each tuple contains the eigenvectors of the matrix of
 *          corresponding tuple of \a this field. This new field lies on the same mesh as
 *          \a this one. The caller is to delete this field using decrRef() as it is no
 *          more needed.  
 *  \throw If \a this->getNumberOfComponents() != 6.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::eigenVectors() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform eigenVectors !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->eigenVectors());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setName("EigenVectors");
  ret->setMesh(getMesh());
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble filled with the inverse matrix of
 * a matrix defined by every tuple of \a this field having either 4, 6 or 9
 * components. The case of 6 components corresponds to that of the upper triangular
 * matrix.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, 
 *          having the same number of components as \a this one, whose each tuple
 *          contains the inverse matrix of the matrix of corresponding tuple of \a this
 *          field. This new field lies on the same mesh as \a this one. The caller is to
 *          delete this field using decrRef() as it is no more needed.  
 *  \throw If \a this->getNumberOfComponents() is not in [4,6,9].
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::inverse() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform inverse !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->inverse());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setName("Inversion");
  ret->setMesh(getMesh());
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble filled with the trace of
 * a matrix defined by every tuple of \a this field having either 4, 6 or 9
 * components. The case of 6 components corresponds to that of the upper triangular
 * matrix.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, 
 *          having 1 component, whose each tuple is the trace of the matrix of
 *          corresponding tuple of \a this field.
 *          This new field lies on the same mesh as \a this one. The caller is to
 *          delete this field using decrRef() as it is no more needed.  
 *  \throw If \a this->getNumberOfComponents() is not in [4,6,9].
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::trace() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform trace !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->trace());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setName("Trace");
  ret->setMesh(getMesh());
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble filled with the stress deviator tensor of
 * a stress tensor defined by every tuple of \a this 6-componental field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, 
 *          having same number of components and tuples as \a this field,
 *          whose each tuple contains the stress deviator tensor of the stress tensor of
 *          corresponding tuple of \a this field. This new field lies on the same mesh as
 *          \a this one. The caller is to delete this field using decrRef() as it is no
 *          more needed.  
 *  \throw If \a this->getNumberOfComponents() != 6.
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::deviator() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform deviator !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->deviator());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setName("Deviator");
  ret->setMesh(getMesh());
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble filled with the magnitude of
 * every vector of \a this field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, 
 *          having one component, whose each tuple is the magnitude of the vector
 *          of corresponding tuple of \a this field. This new field lies on the
 *          same mesh as \a this one. The caller is to
 *          delete this field using decrRef() as it is no more needed.  
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::magnitude() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform magnitude !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->magnitude());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setName("Magnitude");
  ret->setMesh(getMesh());
  return ret.retn();
}

/*!
 * Creates a new scalar MEDCouplingFieldDouble filled with the maximal value among
 * values of every tuple of \a this field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
 *          This new field lies on the same mesh as \a this one. The caller is to
 *          delete this field using decrRef() as it is no more needed.  
 *  \throw If the spatial discretization of \a this field is NULL.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::maxPerTuple() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform maxPerTuple !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->maxPerTuple());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  std::ostringstream oss;
  oss << "Max_" << getName();
  ret->setName(oss.str());
  ret->setMesh(getMesh());
  return ret.retn();
}

/*!
 * Changes number of components in \a this field. If \a newNbOfComp is less
 * than \a this->getNumberOfComponents() then each tuple
 * is truncated to have \a newNbOfComp components, keeping first components. If \a
 * newNbOfComp is more than \a this->getNumberOfComponents() then 
 * each tuple is populated with \a dftValue to have \a newNbOfComp components.  
 *  \param [in] newNbOfComp - number of components for the new field to have.
 *  \param [in] dftValue - value assigned to new values added to \a this field.
 *  \throw If \a this is not allocated.
 */
void MEDCouplingFieldDouble::changeNbOfComponents(int newNbOfComp, double dftValue)
{
  timeDiscr()->changeNbOfComponents(newNbOfComp,dftValue);
}

/*!
 * Creates a new MEDCouplingFieldDouble composed of selected components of \a this field.
 * The new MEDCouplingFieldDouble has the same number of tuples but includes components
 * specified by \a compoIds parameter. So that getNbOfElems() of the result field
 * can be either less, same or more than \a this->getNumberOfValues().
 *  \param [in] compoIds - sequence of zero based indices of components to include
 *              into the new field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If a component index (\a i) is not valid: 
 *         \a i < 0 || \a i >= \a this->getNumberOfComponents().
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::keepSelectedComponents(const std::vector<int>& compoIds) const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform keepSelectedComponents !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->keepSelectedComponents(compoIds));
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setName(getName());
  ret->setMesh(getMesh());
  return ret.retn();
}


/*!
 * Copy all components in a specified order from another field.
 * The number of tuples in \a this and the other field can be different.
 *  \param [in] f - the field to copy data from.
 *  \param [in] compoIds - sequence of zero based indices of components, data of which is
 *              to be copied.
 *  \throw If the two fields have different number of data arrays.
 *  \throw If a data array is set in one of fields and is not set in the other.
 *  \throw If \a compoIds.size() != \a a->getNumberOfComponents().
 *  \throw If \a compoIds[i] < 0 or \a compoIds[i] > \a this->getNumberOfComponents().
 */
void MEDCouplingFieldDouble::setSelectedComponents(const MEDCouplingFieldDouble *f, const std::vector<int>& compoIds)
{
  timeDiscr()->setSelectedComponents(f->timeDiscr(),compoIds);
}

/*!
 * Sorts value within every tuple of \a this field.
 *  \param [in] asc - if \a true, the values are sorted in ascending order, else,
 *              in descending order.
 *  \throw If a data array is not allocated.
 */
void MEDCouplingFieldDouble::sortPerTuple(bool asc)
{
  timeDiscr()->sortPerTuple(asc);
}

/*!
 * Creates a new MEDCouplingFieldDouble by concatenating two given fields.
 * Values of
 * the first field precede values of the second field within the result field.
 *  \param [in] f1 - the first field.
 *  \param [in] f2 - the second field.
 *  \return MEDCouplingFieldDouble * - the result field. It is a new instance of
 *          MEDCouplingFieldDouble. The caller is to delete this mesh using decrRef() 
 *          as it is no more needed.
 *  \throw If the fields are not compatible for the merge.
 *  \throw If the spatial discretization of \a f1 is NULL.
 *  \throw If the time discretization of \a f1 is NULL.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_MergeFields "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_MergeFields "Here is a Python example".
 *  \endif
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::MergeFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1->areCompatibleForMerge(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply MergeFields on them ! Check support mesh, field nature, and spatial and time discretisation.");
  const MEDCouplingMesh *m1(f1->getMesh()),*m2(f2->getMesh());
  if(!f1->timeDiscr())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::MergeFields : no time discr of f1 !");
  if(!f1->_type)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::MergeFields : no spatial discr of f1 !");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->aggregate(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone()));
  ret->setName(f1->getName());
  ret->setDescription(f1->getDescription());
  if(m1)
    {
      MCAuto<MEDCouplingMesh> m=m1->mergeMyselfWith(m2);
      ret->setMesh(m);
    }
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble by concatenating all given fields.
 * Values of the *i*-th field precede values of the (*i*+1)-th field within the result.
 * If there is only one field in \a a, a deepCopy() (except time information of mesh and
 * field) of the field is returned. 
 * Generally speaking the first field in \a a is used to assign tiny attributes of the
 * returned field. 
 *  \param [in] a - a vector of fields (MEDCouplingFieldDouble) to concatenate.
 *  \return MEDCouplingFieldDouble * - the result field. It is a new instance of
 *          MEDCouplingFieldDouble. The caller is to delete this mesh using decrRef() 
 *          as it is no more needed.
 *  \throw If \a a is empty.
 *  \throw If the fields are not compatible for the merge.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_MergeFields "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_MergeFields "Here is a Python example".
 *  \endif
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::MergeFields(const std::vector<const MEDCouplingFieldDouble *>& a)
{
  if(a.size()<1)
    throw INTERP_KERNEL::Exception("FieldDouble::MergeFields : size of array must be >= 1 !");
  std::vector< MCAuto<MEDCouplingUMesh> > ms(a.size());
  std::vector< const MEDCouplingUMesh *> ms2(a.size());
  std::vector< const MEDCouplingTimeDiscretization *> tds(a.size());
  std::vector<const MEDCouplingFieldDouble *>::const_iterator it=a.begin();
  const MEDCouplingFieldDouble *ref=(*it++);
  if(!ref)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::MergeFields : presence of NULL instance in first place of input vector !");
  for(;it!=a.end();it++)
    if(!ref->areCompatibleForMerge(*it))
      throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply MergeFields on them! Check support mesh, field nature, and spatial and time discretisation.");
  for(int i=0;i<(int)a.size();i++)
    {
      if(a[i]->getMesh())
        { ms[i]=a[i]->getMesh()->buildUnstructured(); ms2[i]=ms[i]; }
      else
        { ms[i]=0; ms2[i]=0; }
      tds[i]=a[i]->timeDiscr();
    }
  MEDCouplingTimeDiscretization *td(tds[0]->aggregate(tds));
  td->copyTinyAttrFrom(*(a[0]->timeDiscr()));
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(a[0]->getNature(),td,a[0]->_type->clone()));
  ret->setName(a[0]->getName());
  ret->setDescription(a[0]->getDescription());
  if(ms2[0])
    {
      MCAuto<MEDCouplingUMesh> m(MEDCouplingUMesh::MergeUMeshes(ms2));
      m->copyTinyInfoFrom(ms2[0]);
      ret->setMesh(m);
    }
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble by concatenating components of two given fields.
 * The number of components in the result field is a sum of the number of components of
 * given fields. The number of tuples in the result field is same as that of each of given
 * arrays.
 * Number of tuples in the given fields must be the same.
 *  \param [in] f1 - a field to include in the result field.
 *  \param [in] f2 - another field to include in the result field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If the fields are not compatible for a meld (areCompatibleForMeld()).
 *  \throw If any of data arrays is not allocated.
 *  \throw If \a f1->getNumberOfTuples() != \a f2->getNumberOfTuples()
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::MeldFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1 || !f2)
    throw INTERP_KERNEL::Exception("MeldFields : null input pointer !");
  if(!f1->areCompatibleForMeld(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply MeldFields on them ! Check support mesh, field nature, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->meld(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble containing a dot product of two given fields, 
 * so that the i-th tuple of the result field is a sum of products of j-th components of
 * i-th tuples of given fields (\f$ f_i = \sum_{j=1}^n f1_j * f2_j \f$). 
 * Number of tuples and components in the given fields must be the same.
 *  \param [in] f1 - a given field.
 *  \param [in] f2 - another given field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If either \a f1 or \a f2 is NULL.
 *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
 *         differ not only in values.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::DotFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::DotFields : input field is NULL !");
  if(!f1->areStrictlyCompatibleForMulDiv(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply DotFields on them!  Check support mesh, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->dot(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MEDCouplingFieldDouble *ret(new MEDCouplingFieldDouble(NoNature,td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret;
}

/*!
 * Returns a new MEDCouplingFieldDouble containing a cross product of two given fields, 
 * so that
 * the i-th tuple of the result field is a 3D vector which is a cross
 * product of two vectors defined by the i-th tuples of given fields.
 * Number of tuples in the given fields must be the same.
 * Number of components in the given fields must be 3.
 *  \param [in] f1 - a given field.
 *  \param [in] f2 - another given field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If either \a f1 or \a f2 is NULL.
 *  \throw If \a f1->getNumberOfComponents() != 3
 *  \throw If \a f2->getNumberOfComponents() != 3
 *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
 *         differ not only in values.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::CrossProductFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::CrossProductFields : input field is NULL !");
  if(!f1->areStrictlyCompatibleForMulDiv(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply CrossProductFields on them! Check support mesh, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->crossProduct(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(NoNature,td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble containing maximal values of two given fields.
 * Number of tuples and components in the given fields must be the same.
 *  \param [in] f1 - a field to compare values with another one.
 *  \param [in] f2 - another field to compare values with the first one.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If either \a f1 or \a f2 is NULL.
 *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
 *         differ not only in values.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_MaxFields "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_MaxFields "Here is a Python example".
 *  \endif
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::MaxFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::MaxFields : input field is NULL !");
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply MaxFields on them! Check support mesh, field nature, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->max(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble containing minimal values of two given fields.
 * Number of tuples and components in the given fields must be the same.
 *  \param [in] f1 - a field to compare values with another one.
 *  \param [in] f2 - another field to compare values with the first one.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If either \a f1 or \a f2 is NULL.
 *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
 *         differ not only in values.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_MaxFields "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_MaxFields "Here is a Python example".
 *  \endif
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::MinFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::MinFields : input field is NULL !");
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply MinFields on them! Check support mesh, field nature, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->min(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret.retn();
}

/*!
 * Returns a copy of \a this field in which sign of all values is reversed.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble
 *         containing the same number of tuples and components as \a this field. 
 *         The caller is to delete this result field using decrRef() as it is no more
 *         needed. 
 *  \throw If the spatial discretization of \a this field is NULL.
 *  \throw If a data array is not allocated.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::negate() const
{
  if(_type.isNull())
    throw INTERP_KERNEL::Exception("No spatial discretization underlying this field to perform negate !");
  MEDCouplingTimeDiscretization *td(timeDiscr()->negate());
  td->copyTinyAttrFrom(*timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(getNature(),td,_type->clone()));
  ret->setMesh(getMesh());
  return ret.retn();
}

/*!
 * Returns a new MEDCouplingFieldDouble containing sum values of corresponding values of
 * two given fields ( _f_ [ i, j ] = _f1_ [ i, j ] + _f2_ [ i, j ] ).
 * Number of tuples and components in the given fields must be the same.
 *  \param [in] f1 - a field to sum up.
 *  \param [in] f2 - another field to sum up.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If either \a f1 or \a f2 is NULL.
 *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
 *         differ not only in values.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::AddFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::AddFields : input field is NULL !");
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply AddFields on them! Check support mesh, field nature, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->add(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret.retn();
}

/*!
 * Adds values of another MEDCouplingFieldDouble to values of \a this one
 * ( _this_ [ i, j ] += _other_ [ i, j ] ) using DataArrayDouble::addEqual().
 * The two fields must have same number of tuples, components and same underlying mesh.
 *  \param [in] other - the field to add to \a this one.
 *  \return const MEDCouplingFieldDouble & - a reference to \a this field.
 *  \throw If \a other is NULL.
 *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
 *         differ not only in values.
 */
const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator+=(const MEDCouplingFieldDouble& other)
{
  if(!areStrictlyCompatible(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply += on them! Check support mesh, field nature, and spatial and time discretisation.");
  timeDiscr()->addEqual(other.timeDiscr());
  return *this;
}

/*!
 * Returns a new MEDCouplingFieldDouble containing subtraction of corresponding values of
 * two given fields ( _f_ [ i, j ] = _f1_ [ i, j ] - _f2_ [ i, j ] ).
 * Number of tuples and components in the given fields must be the same.
 *  \param [in] f1 - a field to subtract from.
 *  \param [in] f2 - a field to subtract.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If either \a f1 or \a f2 is NULL.
 *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
 *         differ not only in values.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::SubstractFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::SubstractFields : input field is NULL !");
  if(!f1->areStrictlyCompatible(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply SubstractFields on them! Check support mesh, field nature, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->substract(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(f1->getNature(),td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret.retn();
}

/*!
 * Subtract values of another MEDCouplingFieldDouble from values of \a this one
 * ( _this_ [ i, j ] -= _other_ [ i, j ] ) using DataArrayDouble::substractEqual().
 * The two fields must have same number of tuples, components and same underlying mesh.
 *  \param [in] other - the field to subtract from \a this one.
 *  \return const MEDCouplingFieldDouble & - a reference to \a this field.
 *  \throw If \a other is NULL.
 *  \throw If the fields are not strictly compatible (areStrictlyCompatible()), i.e. they
 *         differ not only in values.
 */
const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator-=(const MEDCouplingFieldDouble& other)
{
  if(!areStrictlyCompatible(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply -= on them! Check support mesh, field nature, and spatial and time discretisation.");
  timeDiscr()->substractEqual(other.timeDiscr());
  return *this;
}

/*!
 * Returns a new MEDCouplingFieldDouble containing product values of
 * two given fields. There are 2 valid cases.
 * 1.  The fields have same number of tuples and components. Then each value of
 *   the result field (_f_) is a product of the corresponding values of _f1_ and
 *   _f2_, i.e. _f_ [ i, j ] = _f1_ [ i, j ] * _f2_ [ i, j ].
 * 2.  The fields have same number of tuples and one field, say _f2_, has one
 *   component. Then
 *   _f_ [ i, j ] = _f1_ [ i, j ] * _f2_ [ i, 0 ].
 *
 * The two fields must have same number of tuples and same underlying mesh.
 *  \param [in] f1 - a factor field.
 *  \param [in] f2 - another factor field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, with no nature set.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If either \a f1 or \a f2 is NULL.
 *  \throw If the fields are not compatible for multiplication (areCompatibleForMul()),
 *         i.e. they differ not only in values and possibly number of components.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::MultiplyFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::MultiplyFields : input field is NULL !");
  if(!f1->areCompatibleForMul(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply MultiplyFields on them! Check support mesh, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->multiply(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(NoNature,td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret.retn();
}

/*!
 * Multiply values of another MEDCouplingFieldDouble to values of \a this one
 * using DataArrayDouble::multiplyEqual().
 * The two fields must have same number of tuples and same underlying mesh.
 * There are 2 valid cases.
 * 1.  The fields have same number of components. Then each value of
 *   \a other is multiplied to the corresponding value of \a this field, i.e.
 *   _this_ [ i, j ] *= _other_ [ i, j ].
 * 2. The _other_ field has one component. Then
 *   _this_ [ i, j ] *= _other_ [ i, 0 ].
 *
 * The two fields must have same number of tuples and same underlying mesh.
 *  \param [in] other - an field to multiply to \a this one.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, with no nature set.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If \a other is NULL.
 *  \throw If the fields are not strictly compatible for multiplication
 *         (areCompatibleForMul()),
 *         i.e. they differ not only in values and possibly in number of components.
 */
const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator*=(const MEDCouplingFieldDouble& other)
{
  if(!areCompatibleForMul(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply *= on them! Check support mesh, and spatial and time discretisation.");
  timeDiscr()->multiplyEqual(other.timeDiscr());
  _nature = NoNature;
  return *this;
}

/*!
 * Returns a new MEDCouplingFieldDouble containing division of two given fields.
 * There are 2 valid cases.
 * 1.  The fields have same number of tuples and components. Then each value of
 *   the result field (_f_) is a division of the corresponding values of \a f1 and
 *   \a f2, i.e. _f_ [ i, j ] = _f1_ [ i, j ] / _f2_ [ i, j ].
 * 2.  The fields have same number of tuples and _f2_ has one component. Then
 *   _f_ [ i, j ] = _f1_ [ i, j ] / _f2_ [ i, 0 ].
 *
 *  \param [in] f1 - a numerator field.
 *  \param [in] f2 - a denominator field.
 *  \return MEDCouplingFieldDouble * - the new instance of MEDCouplingFieldDouble, with no nature set.
 *          The caller is to delete this result field using decrRef() as it is no more
 *          needed.
 *  \throw If either \a f1 or \a f2 is NULL.
 *  \throw If the fields are not compatible for division (areCompatibleForDiv()),
 *         i.e. they differ not only in values and possibly in number of components.
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::DivideFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::DivideFields : input field is NULL !");
  if(!f1->areCompatibleForDiv(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply DivideFields on them! Check support mesh, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->divide(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(NoNature,td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret.retn();
}

/*!
 * Divide values of \a this field by values of another MEDCouplingFieldDouble
 * using DataArrayDouble::divideEqual().
 * The two fields must have same number of tuples and same underlying mesh.
 * There are 2 valid cases.
 * 1.  The fields have same number of components. Then each value of
 *    \a this field is divided by the corresponding value of \a other one, i.e.
 *   _this_ [ i, j ] /= _other_ [ i, j ].
 * 2.  The \a other field has one component. Then
 *   _this_ [ i, j ] /= _other_ [ i, 0 ].
 *
 *  \warning No check of division by zero is performed!
 *  \param [in] other - an field to divide \a this one by.
 *  \throw If \a other is NULL.
 *  \throw If the fields are not compatible for division (areCompatibleForDiv()),
 *         i.e. they differ not only in values and possibly in number of components.
 */
const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator/=(const MEDCouplingFieldDouble& other)
{
  if(!areCompatibleForDiv(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply /= on them! Check support mesh, and spatial and time discretisation.");
  timeDiscr()->divideEqual(other.timeDiscr());
  _nature = NoNature;
  return *this;
}

/*!
 * Directly called by MEDCouplingFieldDouble::operator^.
 * 
 * \sa MEDCouplingFieldDouble::operator^
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::PowFields(const MEDCouplingFieldDouble *f1, const MEDCouplingFieldDouble *f2)
{
  if(!f1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::PowFields : input field is NULL !");
  if(!f1->areCompatibleForMul(f2))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply PowFields on them! Check support mesh, and spatial and time discretisation.");
  MEDCouplingTimeDiscretization *td(f1->timeDiscr()->pow(f2->timeDiscr()));
  td->copyTinyAttrFrom(*f1->timeDiscr());
  MCAuto<MEDCouplingFieldDouble> ret(new MEDCouplingFieldDouble(NoNature,td,f1->_type->clone()));
  ret->setMesh(f1->getMesh());
  return ret.retn();
}

/*!
 * Directly call MEDCouplingFieldDouble::PowFields static method.
 * 
 * \sa MEDCouplingFieldDouble::PowFields
 */
MEDCouplingFieldDouble *MEDCouplingFieldDouble::operator^(const MEDCouplingFieldDouble& other) const
{
  return PowFields(this,&other);
}

const MEDCouplingFieldDouble &MEDCouplingFieldDouble::operator^=(const MEDCouplingFieldDouble& other)
{
  if(!areCompatibleForDiv(&other))
    throw INTERP_KERNEL::Exception("Fields are not compatible. Unable to apply ^= on them!  Check support mesh, and spatial and time discretisation.");
  timeDiscr()->powEqual(other.timeDiscr());
  _nature = NoNature;
  return *this;
}

/*!
 * Writes the field series \a fs and the mesh the fields lie on in the VTK file \a fileName.
 * If \a fs is empty no file is written.
 * The result file is valid provided that no exception is thrown.
 * \warning All the fields must be named and lie on the same non NULL mesh.
 *  \param [in] fileName - the name of a VTK file to write in.
 *  \param [in] fs - the fields to write.
 *  \param [in] isBinary - specifies the VTK format of the written file. By default true (Binary mode)
 *  \throw If \a fs[ 0 ] == NULL.
 *  \throw If the fields lie not on the same mesh.
 *  \throw If the mesh is not set.
 *  \throw If any of the fields has no name.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcfielddouble_WriteVTK "Here is a C++ example".<br>
 *  \ref  py_mcfielddouble_WriteVTK "Here is a Python example".
 *  \endif
 */
std::string MEDCouplingFieldDouble::WriteVTK(const std::string& fileName, const std::vector<const MEDCouplingFieldDouble *>& fs, bool isBinary)
{
  if(fs.empty())
    return std::string();
  std::size_t nfs=fs.size();
  if(!fs[0])
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::WriteVTK : 1st instance of field is NULL !");
  const MEDCouplingMesh *m=fs[0]->getMesh();
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::WriteVTK : 1st instance of field lies on NULL mesh !");
  for(std::size_t i=1;i<nfs;i++)
    if(fs[i]->getMesh()!=m)
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::WriteVTK : Fields are not lying on a same mesh ! Expected by VTK ! MEDCouplingFieldDouble::setMesh or MEDCouplingFieldDouble::changeUnderlyingMesh can help to that.");
  if(!m)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::WriteVTK : Fields are lying on a same mesh but it is empty !");
  std::string ret(m->getVTKFileNameOf(fileName));
  MCAuto<DataArrayByte> byteArr;
  if(isBinary)
    { byteArr=DataArrayByte::New(); byteArr->alloc(0,1); }
  std::ostringstream coss,noss;
  for(std::size_t i=0;i<nfs;i++)
    {
      const MEDCouplingFieldDouble *cur=fs[i];
      std::string name(cur->getName());
      if(name.empty())
        {
          std::ostringstream oss; oss << "MEDCouplingFieldDouble::WriteVTK : Field in pos #" << i << " has no name !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      TypeOfField typ=cur->getTypeOfField();
      if(typ==ON_CELLS)
        cur->getArray()->writeVTK(coss,8,cur->getName(),byteArr);
      else if(typ==ON_NODES)
        cur->getArray()->writeVTK(noss,8,cur->getName(),byteArr);
      else
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::WriteVTK : only node and cell fields supported for the moment !");
    }
  m->writeVTKAdvanced(ret,coss.str(),noss.str(),byteArr);
  return ret;
}

MCAuto<MEDCouplingFieldDouble> MEDCouplingFieldDouble::voronoizeGen(const Voronizer *vor, double eps) const
{
  checkConsistencyLight();
  if(!vor)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::voronoizeGen : null pointer !");
  MCAuto<MEDCouplingFieldDouble> fieldToWO;
  const MEDCouplingMesh *inpMeshBase(getMesh());
  MCAuto<MEDCouplingUMesh> inpMesh(inpMeshBase->buildUnstructured());
  std::string meshName(inpMesh->getName());
  if(!inpMesh->isPresenceOfQuadratic())
    fieldToWO=clone(false);
  else
    {
      fieldToWO=convertQuadraticCellsToLinear();
      inpMeshBase=fieldToWO->getMesh();
      inpMesh=inpMeshBase->buildUnstructured();
    }
  int nbCells(inpMesh->getNumberOfCells());
  const MEDCouplingFieldDiscretization *disc(fieldToWO->getDiscretization());
  const MEDCouplingFieldDiscretizationGauss *disc2(dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(disc));
  if(!disc2)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::voronoize2D : Not a ON_GAUSS_PT field");
  int nbLocs(disc2->getNbOfGaussLocalization());
  std::vector< MCAuto<MEDCouplingUMesh> > cells(nbCells);
  for(int i=0;i<nbLocs;i++)
    {
      const MEDCouplingGaussLocalization& gl(disc2->getGaussLocalization(i));
      if(gl.getDimension()!=vor->getDimension())
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDouble::voronoize2D : not a 2D one !");
      MCAuto<MEDCouplingUMesh> mesh(gl.buildRefCell());
      const std::vector<double>& coo(gl.getGaussCoords());
      MCAuto<DataArrayDouble> coo2(DataArrayDouble::NewFromStdVector(coo));
      coo2->rearrange(vor->getDimension());
      //
      MCAuto<MEDCouplingUMesh> coo3(MEDCouplingUMesh::Build0DMeshFromCoords(coo2));
      //
      MCAuto<MEDCouplingUMesh> vorCellsForCurDisc(vor->doIt(mesh,coo2,eps));
      std::vector<int> ids;
      MCAuto<DataArrayDouble> ptsInReal;
      disc2->getCellIdsHavingGaussLocalization(i,ids);
      {
        MCAuto<MEDCouplingUMesh> subMesh(inpMesh->buildPartOfMySelf(&ids[0],&ids[0]+ids.size()));
        ptsInReal=gl.localizePtsInRefCooForEachCell(vorCellsForCurDisc->getCoords(),subMesh);
      }
      int nbPtsPerCell(vorCellsForCurDisc->getNumberOfNodes());
      for(std::size_t j=0;j<ids.size();j++)
        {
          MCAuto<MEDCouplingUMesh> elt(vorCellsForCurDisc->clone(false));
          MCAuto<DataArrayDouble> coo4(ptsInReal->selectByTupleIdSafeSlice(j*nbPtsPerCell,(j+1)*nbPtsPerCell,1));
          elt->setCoords(coo4);
          cells[ids[j]]=elt;
        }
    }
  std::vector< const MEDCouplingUMesh * > cellsPtr(VecAutoToVecOfCstPt(cells));
  MCAuto<MEDCouplingUMesh> outMesh(MEDCouplingUMesh::MergeUMeshes(cellsPtr));
  outMesh->setName(meshName);
  MCAuto<MEDCouplingFieldDouble> onCells(MEDCouplingFieldDouble::New(ON_CELLS));
  onCells->setMesh(outMesh);
  {
    MCAuto<DataArrayDouble> arr(fieldToWO->getArray()->deepCopy());
    onCells->setArray(arr);
  }
  onCells->setTimeUnit(getTimeUnit());
  {
    int b,c;
    double a(getTime(b,c));
    onCells->setTime(a,b,c);
  }
  onCells->setName(getName());
  return onCells;
}

MEDCouplingTimeDiscretization *MEDCouplingFieldDouble::timeDiscr()
{
  MEDCouplingTimeDiscretizationTemplate<double> *ret(_time_discr);
  if(!ret)
    return 0;
  MEDCouplingTimeDiscretization *retc(dynamic_cast<MEDCouplingTimeDiscretization *>(ret));
  if(!retc)
    throw INTERP_KERNEL::Exception("Field Double Null invalid type of time discr !");
  return retc;
}

const MEDCouplingTimeDiscretization *MEDCouplingFieldDouble::timeDiscr() const
{
  const MEDCouplingTimeDiscretizationTemplate<double> *ret(_time_discr);
  if(!ret)
    return 0;
  const MEDCouplingTimeDiscretization *retc(dynamic_cast<const MEDCouplingTimeDiscretization *>(ret));
  if(!retc)
    throw INTERP_KERNEL::Exception("Field Double Null invalid type of time discr !");
  return retc;
}
