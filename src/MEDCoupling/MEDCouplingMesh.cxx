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

#include "MEDCouplingMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMemArray.txx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MCAuto.hxx"

#include <set>
#include <cmath>
#include <sstream>
#include <fstream>
#include <iterator>

using namespace MEDCoupling;

MEDCouplingMesh::MEDCouplingMesh():_time(0.),_iteration(-1),_order(-1)
{
}

MEDCouplingMesh::MEDCouplingMesh(const MEDCouplingMesh& other):RefCountObject(other),_name(other._name),_description(other._description),
                                                               _time(other._time),_iteration(other._iteration),
                                                               _order(other._order),_time_unit(other._time_unit)
{
}

std::size_t MEDCouplingMesh::getHeapMemorySizeWithoutChildren() const
{
  return _name.capacity()+_description.capacity()+_time_unit.capacity();
}

/*!
 * This method is only for ParaMEDMEM in ParaFIELD constructor.
 */
bool MEDCouplingMesh::isStructured() const
{
  return getType()==CARTESIAN;
}

bool MEDCouplingMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingMesh::isEqualIfNotWhy : other instance is NULL !");
  std::ostringstream oss; oss.precision(15);
  if(_name!=other->_name)
    {
      oss << "Mesh names differ : this name = \"" << _name << "\" and other name = \"" << other->_name << "\" !";
      reason=oss.str();
      return false;
    }
  if(_description!=other->_description)
    {
      oss << "Mesh descriptions differ : this description = \"" << _description << "\" and other description = \"" << other->_description << "\" !";
      reason=oss.str();
      return false;
    }
  if(_iteration!=other->_iteration)
    {
      oss << "Mesh iterations differ : this iteration = \"" << _iteration << "\" and other iteration = \"" << other->_iteration << "\" !";
      reason=oss.str();
      return false;
    }
  if(_order!=other->_order)
    {
      oss << "Mesh orders differ : this order = \"" << _order << "\" and other order = \"" << other->_order << "\" !";
      reason=oss.str();
      return false;
    }
  if(_time_unit!=other->_time_unit)
    {
      oss << "Mesh time units differ : this time unit = \"" << _time_unit << "\" and other time unit = \"" << other->_time_unit << "\" !";
      reason=oss.str();
      return false;
    }
  if(fabs(_time-other->_time)>=1e-12)
    {
      oss << "Mesh times differ : this time = \"" << _time << "\" and other time = \"" << other->_time << "\" !";
      reason=oss.str();
      return false;
    }
  return true;
}

/*!
 * Checks if \a this and another MEDCouplingMesh are fully equal.
 *  \param [in] other - an instance of MEDCouplingMesh to compare with \a this one.
 *  \param [in] prec - precision value used to compare node coordinates.
 *  \return bool - \c true if the two meshes are equal, \c false else.
 */
bool MEDCouplingMesh::isEqual(const MEDCouplingMesh *other, double prec) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,prec,tmp);
}

/*!
 * This method checks geo equivalence between two meshes : \a this and \a other.
 * If no exception is thrown \a this and \a other are geometrically equivalent regarding \a levOfCheck level.
 * This method is typically used to change the mesh of a field "safely" depending the \a levOfCheck level considered.
 * 
 * In case of success cell \c other[i] is equal to the cell \c this[cellCor[i]].
 * In case of success node \c other->getCoords()[i] is equal to the node \c this->getCoords()[nodeCor[i]].
 *
 * If \a cellCor is null (or Py_None) it means that for all #i cell in \a other is equal to cell # i in \a this.
 *
 * If \a nodeCor is null (or Py_None) it means that for all #i node in \a other is equal to node # i in \a this.
 *
 * So null (or Py_None) returned in \a cellCor and/or \a nodeCor means identity array. This is for optimization reason to avoid to build useless arrays
 * for some \a levOfCheck (for example 0).
 *
 * **Warning a not null output does not mean that it is not identity !**
 *
 * \param [in] other - the mesh to be compared with \a this.
 * \param [in] levOfCheck - input that specifies the level of check specified. The possible values are listed below.
 * \param [in] prec - input that specifies precision for double float data used for comparison in meshes.
 * \param [out] cellCor - output array not always informed (depending \a levOfCheck param) that gives the corresponding array for cells from \a other to \a this.
 * \param [out] nodeCor - output array not always informed (depending \a levOfCheck param) that gives the corresponding array for nodes from \a other to \a this.
 *
 * Possible values for levOfCheck :
 *   - 0 for strict equality. This is the strongest level. \a cellCor and \a nodeCor params are never informed.
 *   - 10,11,12 (10+x) for less strict equality. Two meshes are compared geometrically. In case of success \a cellCor and \a nodeCor are informed. Warning ! These equivalences are CPU/Mem costly. The 3 values correspond respectively to policy used for cell comparison (see MEDCouplingUMesh::zipConnectivityTraducer to have more details)
 *   - 20,21,22 (20+x), for less strict equality. Two meshes are compared geometrically. The difference with the previous version is that nodes(coordinates) are expected to be the same between this and other. In case of success \a cellCor is informed. Warning ! These equivalences are CPU/Mem costly. The 3 values correspond respectively to policy used for cell comparison (see MEDCouplingUMesh::zipConnectivityTraducer to have more details)
 *   - 1 for fast 'equality'. This is a lazy level. Just number of cells and number of nodes are considered here and 3 cells (begin,middle,end)
 *   - 2 for deep 'equality' as 0 option except that no control is done on all strings in mesh.
 *
 * So the most strict level of check is 0 (equality). The least strict is 12. If the level of check 12 throws, the 2 meshes \a this and \a other are not similar enough
 * to be compared. An interpolation using MEDCouplingRemapper class should be then used.
 */
void MEDCouplingMesh::checkGeoEquivalWith(const MEDCouplingMesh *other, int levOfCheck, double prec,
                                          DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const
{
  cellCor=0;
  nodeCor=0;
  if(this==other)
    return ;
  switch(levOfCheck)
    {
    case 0:
      {
        if(!isEqual(other,prec))
          throw INTERP_KERNEL::Exception("checkGeoFitWith : Meshes are not equal !");
        return ;
      }
    case 10:
    case 11:
    case 12:
      {
        checkDeepEquivalWith(other,levOfCheck-10,prec,cellCor,nodeCor);
        return ;
      }
    case 20:
    case 21:
    case 22:
      {
        checkDeepEquivalOnSameNodesWith(other,levOfCheck-20,prec,cellCor);
        return ;
      }
    case 1:
      {
        checkFastEquivalWith(other,prec);
        return;
      }
    case 2:
      {
        if(!isEqualWithoutConsideringStr(other,prec))
          throw INTERP_KERNEL::Exception("checkGeoFitWith : Meshes are not equal without considering strings !");
        return ;
      }
    default:
      throw INTERP_KERNEL::Exception("checkGeoFitWith : Invalid levOfCheck specified ! Value must be in 0,1,2,10,11 or 12.");
    }
}

/*!
 * Finds cells whose all nodes are in a given array of node ids.
 *  \param [in] partBg - the array of node ids.
 *  \param [in] partEnd - end of \a partBg, i.e. a pointer to a (last+1)-th element
 *          of \a partBg.
 *  \return DataArrayInt * - a new instance of DataArrayInt holding ids of found
 *          cells. The caller is to delete this array using decrRef() as it is no
 *          more needed.
 */
DataArrayInt *MEDCouplingMesh::getCellIdsFullyIncludedInNodeIds(const int *partBg, const int *partEnd) const
{
  std::vector<int> crest;
  std::set<int> p(partBg,partEnd);
  int nbOfCells=getNumberOfCells();
  for(int i=0;i<nbOfCells;i++)
    {
      std::vector<int> conn;
      getNodeIdsOfCell(i,conn);
      bool cont=true;
      for(std::vector<int>::const_iterator iter=conn.begin();iter!=conn.end() && cont;iter++)
        if(p.find(*iter)==p.end())
          cont=false;
      if(cont)
        crest.push_back(i);
    }
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc((int)crest.size(),1);
  std::copy(crest.begin(),crest.end(),ret->getPointer());
  return ret;
}

/*!
 * This method checks fastly that \a this and \a other are equal. All common checks are done here.
 */
void MEDCouplingMesh::checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingMesh::checkFastEquivalWith : input mesh is null !");
  if(getMeshDimension()!=other->getMeshDimension())
    throw INTERP_KERNEL::Exception("checkFastEquivalWith : Mesh dimensions are not equal !");
  if(getSpaceDimension()!=other->getSpaceDimension())
    throw INTERP_KERNEL::Exception("checkFastEquivalWith : Space dimensions are not equal !");
  if(getNumberOfCells()!=other->getNumberOfCells())
    throw INTERP_KERNEL::Exception("checkFastEquivalWith : number of cells are not equal !");
}

/*!
 * This method is very poor and looks only if \a this and \a other are candidate for merge of fields lying repectively on them.
 */
bool MEDCouplingMesh::areCompatibleForMerge(const MEDCouplingMesh *other) const
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingMesh::areCompatibleForMerge : input mesh is null !");
  if(getMeshDimension()!=other->getMeshDimension())
    return false;
  if(getSpaceDimension()!=other->getSpaceDimension())
    return false;
  return true;
}

/*!
 * This method is equivalent to MEDCouplingMesh::buildPart method except that here the cell ids are specified using slice \a beginCellIds \a endCellIds and \a stepCellIds.
 * \b WARNING , there is a big difference compared to MEDCouplingMesh::buildPart method.
 * If the input range is equal all cells in \a this, \a this is returned !
 *
 * \return a new ref to be managed by the caller. Warning this ref can be equal to \a this if input slice is exactly equal to the whole cells in the same order.
 *
 * \sa MEDCouplingMesh::buildPart
 */
MEDCouplingMesh *MEDCouplingMesh::buildPartRange(int beginCellIds, int endCellIds, int stepCellIds) const
{
  if(beginCellIds==0 && endCellIds==getNumberOfCells() && stepCellIds==1)
    {
      MEDCouplingMesh *ret(const_cast<MEDCouplingMesh *>(this));
      ret->incrRef();
      return ret;
    }
  else
    {
      MCAuto<DataArrayInt> cellIds=DataArrayInt::Range(beginCellIds,endCellIds,stepCellIds);
      return buildPart(cellIds->begin(),cellIds->end());
    }
}

/*!
 * This method is equivalent to MEDCouplingMesh::buildPartAndReduceNodes method except that here the cell ids are specified using slice \a beginCellIds \a endCellIds and \a stepCellIds.
 *
 * \sa MEDCouplingMesh::buildPartAndReduceNodes
 */
MEDCouplingMesh *MEDCouplingMesh::buildPartRangeAndReduceNodes(int beginCellIds, int endCellIds, int stepCellIds, int& beginOut, int& endOut, int& stepOut, DataArrayInt*& arr) const
{
  MCAuto<DataArrayInt> cellIds=DataArrayInt::Range(beginCellIds,endCellIds,stepCellIds);
  return buildPartAndReduceNodes(cellIds->begin(),cellIds->end(),arr);
}

/*!
 * This method builds a field lying on \a this with 'nbOfComp' components.
 * 'func' is a pointer that points to a function that takes 2 arrays in parameter and returns a boolean.
 * The first array is a in-param of size this->getSpaceDimension and the second an out param of size 'nbOfComp'.
 * The return field will have type specified by 't'. 't' is also used to determine where values of field will be
 * evaluate.
 * Contrary to other fillFromAnalytic methods this method requests a C++ function pointer as input.
 * The 'func' is a callback that takes as first parameter an input array of size 'this->getSpaceDimension()',
 * the second parameter is a pointer on a valid zone of size at least equal to 'nbOfComp' values. And too finish
 * the returned value is a boolean that is equal to False in case of invalid evaluation (log(0) for example...)
 * 
 * \param t type of field returned and specifies where the evaluation of func will be done.
 * \param nbOfComp number of components of returned field.
 * \param func pointer to a function that should return false if the evaluation failed. (division by 0. for example)
 * \return field with counter = 1.
 */
MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalytic(TypeOfField t, int nbOfComp, FunctionToEvaluate func) const
{
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(t,ONE_TIME);
  ret->setMesh(this);
  ret->fillFromAnalytic(nbOfComp,func);
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * This method copyies all tiny strings from other (name and components name).
 * @throw if other and this have not same mesh type.
 */
void MEDCouplingMesh::copyTinyStringsFrom(const MEDCouplingMesh *other)
{
  if(!other)
    throw INTERP_KERNEL::Exception("MEDCouplingMesh::copyTinyStringsFrom : input mesh is null !");
  _name=other->_name;
  _description=other->_description;
  _time_unit=other->_time_unit;
}

/*!
 * This method copies all attributes that are \b NOT arrays in this.
 * All tiny attributes not usefully for state of \a this are ignored.
 */
void MEDCouplingMesh::copyTinyInfoFrom(const MEDCouplingMesh *other)
{
  _time=other->_time;
  _iteration=other->_iteration;
  _order=other->_order;
  copyTinyStringsFrom(other);
}

/*!
 * \anchor mcmesh_fillFromAnalytic
 * Creates a new MEDCouplingFieldDouble of a given type, one time, with given number of
 * components, lying on \a this mesh, with contents got by applying a specified
 * function to coordinates of field location points (defined by the given field type).
 * For example, if \a t == MEDCoupling::ON_CELLS, the function is applied to cell
 * barycenters.<br>
 * For more info on supported expressions that can be used in the function, see \ref
 * MEDCouplingArrayApplyFuncExpr. The function can include arbitrary named variables
 * (e.g. "x","y" or "va44") to refer to components of point coordinates. Names of
 * variables are sorted in \b alphabetical \b order to associate a variable name with a
 * component. For example, in the expression "2*x+z", "x" stands for the component #0
 * and "z" stands for the component #1 (\b not #2)!<br>
 * In a general case, a value resulting from the function evaluation is assigned to all
 * components of the field. But there is a possibility to have its own expression for
 * each component within one function. For this purpose, there are predefined variable
 * names (IVec, JVec, KVec, LVec etc) each dedicated to a certain component (IVec, to
 * the component #0 etc). A factor of such a variable is added to the
 * corresponding component only.<br>
 * For example, \a nbOfComp == 4, \a this->getSpaceDimension() == 3, coordinates of a
 * point are (1.,3.,7.), then
 *   - "2*x + z"               produces (5.,5.,5.,5.)
 *   - "2*x + 0*y + z"         produces (9.,9.,9.,9.)
 *   - "2*x*IVec + (x+z)*LVec" produces (2.,0.,0.,4.)
 *   - "2*y*IVec + z*KVec + x" produces (7.,1.,1.,4.)
 *
 *  \param [in] t - the field type. It defines, apart from other things, points to
 *         coordinates of which the function is applied to get field values.
 *  \param [in] nbOfComp - the number of components in the result field.
 *  \param [in] func - a string defining the expression which is evaluated to get
 *         field values.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcmesh_fillFromAnalytic "Here is a C++ example".<br>
 *  \ref  py_mcmesh_fillFromAnalytic "Here is a Python example".
 *  \endif
 */
MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalytic(TypeOfField t, int nbOfComp, const std::string& func) const
{
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(t,ONE_TIME);
  ret->setMesh(this);
  ret->fillFromAnalytic(nbOfComp,func);
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble of a given type, one time, with given number of
 * components, lying on \a this mesh, with contents got by applying a specified
 * function to coordinates of field location points (defined by the given field type).
 * For example, if \a t == MEDCoupling::ON_CELLS, the function is applied to cell
 * barycenters. This method differs from
 * \ref MEDCouplingMesh::fillFromAnalytic(TypeOfField t, int nbOfComp, const std::string& func) const "fillFromAnalytic()"
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
 * For example, \a nbOfComp == 4, \a this->getSpaceDimension() == 3, names of
 * spatial components are "x", "y" and "z", coordinates of a
 * point are (1.,3.,7.), then
 *   - "2*x + z"               produces (9.,9.,9.,9.)
 *   - "2*x*IVec + (x+z)*LVec" produces (2.,0.,0.,8.)
 *   - "2*y*IVec + z*KVec + x" produces (7.,1.,1.,8.)
 *
 *  \param [in] t - the field type. It defines, apart from other things, the points to
 *         coordinates of which the function is applied to get field values.
 *  \param [in] nbOfComp - the number of components in the result field.
 *  \param [in] func - a string defining the expression which is evaluated to get
 *         field values.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the node coordinates are not defined.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcmesh_fillFromAnalytic2 "Here is a C++ example".<br>
 *  \ref  py_mcmesh_fillFromAnalytic2 "Here is a Python example".
 *  \endif
 */
MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalyticCompo(TypeOfField t, int nbOfComp, const std::string& func) const
{
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(t,ONE_TIME);
  ret->setMesh(this);
  ret->fillFromAnalyticCompo(nbOfComp,func);
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingFieldDouble of a given type, one time, with given number of
 * components, lying on \a this mesh, with contents got by applying a specified
 * function to coordinates of field location points (defined by the given field type).
 * For example, if \a t == MEDCoupling::ON_CELLS, the function is applied to cell
 * barycenters. This method differs from \ref  \ref mcmesh_fillFromAnalytic
 * "fillFromAnalytic()" by the way how variable
 * names, used in the function, are associated with components of coordinates of field
 * location points; here, a component index of a variable is defined by a
 * rank of the variable within the input array \a varsOrder.<br>
 * For more info on supported expressions that can be used in the function, see \ref
 * MEDCouplingArrayApplyFuncExpr.
 * In a general case, a value resulting from the function evaluation is assigned to all
 * components of the field. But there is a possibility to have its own expression for
 * each component within one function. For this purpose, there are predefined variable
 * names (IVec, JVec, KVec, LVec etc) each dedicated to a certain component (IVec, to
 * the component #0 etc). A factor of such a variable is added to the
 * corresponding component only.<br>
 * For example, \a nbOfComp == 4, \a this->getSpaceDimension() == 3, names of
 * spatial components are given in \a varsOrder: ["x", "y","z"], coordinates of a
 * point are (1.,3.,7.), then
 *   - "2*x + z"               produces (9.,9.,9.,9.)
 *   - "2*x*IVec + (x+z)*LVec" produces (2.,0.,0.,8.)
 *   - "2*y*IVec + z*KVec + x" produces (7.,1.,1.,8.)
 *
 *  \param [in] t - the field type. It defines, apart from other things, the points to
 *         coordinates of which the function is applied to get field values.
 *  \param [in] nbOfComp - the number of components in the result field.
 *  \param [in] varsOrder - the vector defining names of variables used to refer to
 *         components of coordinates of field location points. A variable named
 *         varsOrder[0] refers to the component #0 etc.
 *  \param [in] func - a string defining the expression which is evaluated to get
 *         field values.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble. The
 *         caller is to delete this field using decrRef() as it is no more needed. 
 *  \throw If the node coordinates are not defined.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If computing \a func fails.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcmesh_fillFromAnalytic3 "Here is a C++ example".<br>
 *  \ref  py_mcmesh_fillFromAnalytic3 "Here is a Python example".
 *  \endif
 */
MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalyticNamedCompo(TypeOfField t, int nbOfComp, const std::vector<std::string>& varsOrder, const std::string& func) const
{
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(t,ONE_TIME);
  ret->setMesh(this);
  ret->fillFromAnalyticNamedCompo(nbOfComp,varsOrder,func);
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

/*!
 * Creates a new MEDCouplingMesh by concatenating two given meshes, if possible.
 * Cells and nodes of
 * the first mesh precede cells and nodes of the second mesh within the result mesh.
 * The meshes must be of the same mesh type, else, an exception is thrown. The method
 * MergeMeshes(), accepting a vector of input meshes, has no such a limitation.
 *  \param [in] mesh1 - the first mesh.
 *  \param [in] mesh2 - the second mesh.
 *  \return MEDCouplingMesh * - the result mesh. It is a new instance of
 *          MEDCouplingMesh. The caller is to delete this mesh using decrRef() as it
 *          is no more needed.
 *  \throw If the meshes are of different mesh type.
 */
MEDCouplingMesh *MEDCouplingMesh::MergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2)
{
  if(!mesh1)
    throw INTERP_KERNEL::Exception("MEDCouplingMesh::MergeMeshes : first parameter is an empty mesh !");
  if(!mesh2)
    throw INTERP_KERNEL::Exception("MEDCouplingMesh::MergeMeshes : second parameter is an empty mesh !");
  return mesh1->mergeMyselfWith(mesh2);
}

/*!
 * Creates a new MEDCouplingMesh by concatenating all given meshes, if possible.
 * Cells and nodes of
 * the *i*-th mesh precede cells and nodes of the (*i*+1)-th mesh within the result mesh.
 * This method performs a systematic conversion to unstructured meshes before
 * performing aggregation contrary to the other MergeMeshes()
 * with two parameters that works only on the same type of meshes. So here it is possible
 * to mix different type of meshes. 
 *  \param [in] meshes - a vector of meshes to concatenate.
 *  \return MEDCouplingMesh * - the result mesh. It is a new instance of
 *          MEDCouplingUMesh. The caller is to delete this mesh using decrRef() as it
 *          is no more needed.
 *  \throw If \a meshes.size() == 0.
 *  \throw If \a size[ *i* ] == NULL.
 *  \throw If the coordinates is not set in none of the meshes.
 *  \throw If \a meshes[ *i* ]->getMeshDimension() < 0.
 *  \throw If the \a meshes are of different dimension (getMeshDimension()).
 */
MEDCouplingMesh *MEDCouplingMesh::MergeMeshes(std::vector<const MEDCouplingMesh *>& meshes)
{
  std::vector< MCAuto<MEDCouplingUMesh> > ms1(meshes.size());
  std::vector< const MEDCouplingUMesh * > ms2(meshes.size());
  for(std::size_t i=0;i<meshes.size();i++)
    {
      if(meshes[i])
        {
          MEDCouplingUMesh *cur=meshes[i]->buildUnstructured();
          ms1[i]=cur;  ms2[i]=cur;
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingMesh::MergeMeshes(std::vector<const MEDCouplingMesh *>& meshes) : mesh at pos #" << i << " of input vector of size " << meshes.size() << " is empty !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return MEDCouplingUMesh::MergeUMeshes(ms2);
}

/*!
 * For example if \a type is INTERP_KERNEL::NORM_TRI3 , INTERP_KERNEL::NORM_POLYGON is returned.
 * If \a type is INTERP_KERNEL::NORM_HEXA8 , INTERP_KERNEL::NORM_POLYHED is returned.
 * 
 * \param [in] type the geometric type for which the corresponding dynamic type, is asked.
 * \return the corresponding dynamic type, able to store the input \a type.
 * 
 * \throw if type is equal to \c INTERP_KERNEL::NORM_ERROR or to an unexisting geometric type.
 */
INTERP_KERNEL::NormalizedCellType MEDCouplingMesh::GetCorrespondingPolyType(INTERP_KERNEL::NormalizedCellType type)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  return cm.getCorrespondingPolyType();
}

/*!
 * \param [in] type the geometric type for which the number of nodes consituting it, is asked.
 * \return number of nodes consituting the input geometric type \a type.
 * 
 * \throw if type is dynamic as \c INTERP_KERNEL::NORM_POLYHED , \c INTERP_KERNEL::NORM_POLYGON , \c INTERP_KERNEL::NORM_QPOLYG
 * \throw if type is equal to \c INTERP_KERNEL::NORM_ERROR or to an unexisting geometric type.
 */
int MEDCouplingMesh::GetNumberOfNodesOfGeometricType(INTERP_KERNEL::NormalizedCellType type)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(cm.isDynamic())
    throw INTERP_KERNEL::Exception("MEDCouplingMesh::GetNumberOfNodesOfGeometricType : the input geometric type is dynamic ! Impossible to return a fixed number of nodes constituting it !");
  return (int) cm.getNumberOfNodes();
}

/*!
 * \param [in] type the geometric type for which the status static/dynamic is asked.
 * \return true for static geometric type, false for dynamic geometric type.
 * 
 * \throw if type is equal to \c INTERP_KERNEL::NORM_ERROR or to an unexisting geometric type.
 */
bool MEDCouplingMesh::IsStaticGeometricType(INTERP_KERNEL::NormalizedCellType type)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  return !cm.isDynamic();
}

bool MEDCouplingMesh::IsLinearGeometricType(INTERP_KERNEL::NormalizedCellType type)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  return !cm.isQuadratic();
}

/*!
 * \param [in] type the geometric type for which the dimension is asked.
 * \return the dimension associated to the input geometric type \a type.
 * 
 * \throw if type is equal to \c INTERP_KERNEL::NORM_ERROR or to an unexisting geometric type.
 */
int MEDCouplingMesh::GetDimensionOfGeometricType(INTERP_KERNEL::NormalizedCellType type)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  return (int) cm.getDimension();
}

/*!
 * \param [in] type the geometric type for which the representation is asked.
 * \return the string representation corresponding to the input geometric type \a type.
 * 
 * \throw if type is equal to \c INTERP_KERNEL::NORM_ERROR or to an unexisting geometric type.
 */
const char *MEDCouplingMesh::GetReprOfGeometricType(INTERP_KERNEL::NormalizedCellType type)
{
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  return cm.getRepr();
}

/*!
 * Finds cells in contact with several balls (i.e. points with precision).
 * This method is an extension of getCellContainingPoint() and
 * getCellsContainingPoint() for the case of multiple points.
 *  \param [in] pos - an array of coordinates of points in full interlace mode :
 *         X0,Y0,Z0,X1,Y1,Z1,... Size of the array must be \a
 *         this->getSpaceDimension() * \a nbOfPoints 
 *  \param [in] nbOfPoints - number of points to locate within \a this mesh.
 *  \param [in] eps - radius of balls (i.e. the precision).
 *  \param [out] elts - vector returning ids of found cells.
 *  \param [out] eltsIndex - an array, of length \a nbOfPoints + 1,
 *         dividing cell ids in \a elts into groups each referring to one
 *         point. Its every element (except the last one) is an index pointing to the
 *         first id of a group of cells. For example cells in contact with the *i*-th
 *         point are described by following range of indices:
 *         [ \a eltsIndex[ *i* ], \a eltsIndex[ *i*+1 ] ) and the cell ids are
 *         \a elts[ \a eltsIndex[ *i* ]], \a elts[ \a eltsIndex[ *i* ] + 1 ], ...
 *         Number of cells in contact with the *i*-th point is
 *         \a eltsIndex[ *i*+1 ] - \a eltsIndex[ *i* ].
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcumesh_getCellsContainingPoints "Here is a C++ example".<br>
 *  \ref  py_mcumesh_getCellsContainingPoints "Here is a Python example".
 *  \endif
 */
void MEDCouplingMesh::getCellsContainingPoints(const double *pos, int nbOfPoints, double eps, MCAuto<DataArrayInt>& elts, MCAuto<DataArrayInt>& eltsIndex) const
{
  eltsIndex=DataArrayInt::New(); elts=DataArrayInt::New(); eltsIndex->alloc(nbOfPoints+1,1); eltsIndex->setIJ(0,0,0); elts->alloc(0,1);
  int *eltsIndexPtr(eltsIndex->getPointer());
  int spaceDim(getSpaceDimension());
  const double *work(pos);
  for(int i=0;i<nbOfPoints;i++,work+=spaceDim)
    {
      std::vector<int> ret;
      getCellsContainingPoint(work,eps,ret);
      elts->insertAtTheEnd(ret.begin(),ret.end());
      eltsIndexPtr[i+1]=elts->getNumberOfTuples();
    }
}

/*!
 * Writes \a this mesh into a VTK format file named as specified.
 *  \param [in] fileName - the name of the file to write in. If the extension is OK the fileName will be used directly.
 *                         If extension is invalid or no extension the right extension will be appended.
 *  \return - the real fileName
 *  \throw If \a fileName is not a writable file.
 *  \sa getVTKFileNameOf
 */
std::string MEDCouplingMesh::writeVTK(const std::string& fileName, bool isBinary) const
{
  std::string ret(getVTKFileNameOf(fileName));
  //
  std::string cda,pda;
  MCAuto<DataArrayByte> byteArr;
  if(isBinary)
    { byteArr=DataArrayByte::New(); byteArr->alloc(0,1); }
  writeVTKAdvanced(ret,cda,pda,byteArr);
  return ret;
}

/*!
 * This method takes in input a file name \a fileName and considering the VTK extension of \a this (depending on the type of \a this)
 * returns a right file name. If the input \a fileName has a valid extension the returned string is equal to \a fileName.
 *
 * \sa  getVTKFileExtension
 */
std::string MEDCouplingMesh::getVTKFileNameOf(const std::string& fileName) const
{
  std::string ret;
  std::string part0,part1;
  SplitExtension(fileName,part0,part1);
  std::string ext("."); ext+=getVTKFileExtension();
  if(part1==ext)
    ret=fileName;
  else
    ret=fileName+ext;
  return ret;
}

/// @cond INTERNAL
void MEDCouplingMesh::writeVTKAdvanced(const std::string& fileName, const std::string& cda, const std::string& pda, DataArrayByte *byteData) const
{
  std::ofstream ofs(fileName.c_str());
  ofs << "<VTKFile type=\""  << getVTKDataSetType() << "\" version=\"0.1\" byte_order=\"" << MEDCouplingByteOrderStr() << "\">\n";
  writeVTKLL(ofs,cda,pda,byteData);
  if(byteData)
    {
      ofs << "<AppendedData encoding=\"raw\">\n_1234";
      ofs << std::flush; ofs.close();
      std::ofstream ofs2(fileName.c_str(),std::ios_base::binary | std::ios_base::app);
      ofs2.write(byteData->begin(),byteData->getNbOfElems()); ofs2 << std::flush; ofs2.close();
      std::ofstream ofs3(fileName.c_str(),std::ios_base::app); ofs3 << "\n</AppendedData>\n</VTKFile>\n"; ofs3.close();
    }
  else
    {
      ofs << "</VTKFile>\n";
      ofs.close();
    }
}

void MEDCouplingMesh::SplitExtension(const std::string& fileName, std::string& baseName, std::string& extension)
{
  std::size_t pos(fileName.find_last_of('.'));
  if(pos==std::string::npos)
    {
      baseName=fileName;
      extension.clear();
      return ;
    }
  baseName=fileName.substr(0,pos);
  extension=fileName.substr(pos);
}
/// @endcond
