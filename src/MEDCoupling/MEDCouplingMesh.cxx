// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#include "MEDCouplingMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include <set>
#include <cmath>
#include <sstream>
#include <iterator>

using namespace ParaMEDMEM;

MEDCouplingMesh::MEDCouplingMesh():_time(0.),_iteration(-1),_order(-1)
{
}

MEDCouplingMesh::MEDCouplingMesh(const MEDCouplingMesh& other):_name(other._name),_description(other._description),
                                                               _time(other._time),_iteration(other._iteration),
                                                               _order(other._order),_time_unit(other._time_unit)
{
}

/*!
 * This method is only for ParaMEDMEM in ParaFIELD constructor.
 */
bool MEDCouplingMesh::isStructured() const
{
  return getType()==CARTESIAN;
}

bool MEDCouplingMesh::isEqual(const MEDCouplingMesh *other, double prec) const
{
  return _name==other->_name && _description==other->_description && _iteration==other->_iteration
    && _order==other->_order && _time_unit==other->_time_unit && fabs(_time-other->_time)<1e-12;
}

/*!
 * This method checks geo equivalence between two meshes : 'this' and 'other'.
 * If no exception is throw 'this' and 'other' are geometrically equivalent regarding 'levOfCheck' level.
 * This method is typically used to change the mesh of a field "safely" depending the 'levOfCheck' level considered.
 * 
 * @param levOfCheck input that specifies the level of check specified. The possible values are listed below.
 * @param prec input that specifies precision for double float data used for comparison in meshes.
 * @param cellCor output array not always informed (depending 'levOfCheck' param) that gives the corresponding array for cells from 'other' to 'this'.
 * @param nodeCor output array not always informed (depending 'levOfCheck' param) that gives the corresponding array for nodes from 'other' to 'this'.
 *
 * Possible values for levOfCheck :
 *   - 0 for strict equality. This is the strongest level. 'cellCor' and 'nodeCor' params are never informed.
 *   - 10,11,12 for less strict equality. Two meshes are compared geometrically. In case of success 'cellCor' and 'nodeCor' are informed. Warning ! These equivalences are CPU/Mem costly. The 3 values correspond respectively to policy used for cell comparison (see MEDCouplingUMesh::zipConnectivityTraducer to have more details)
 *   - 20,21,22, for less strict equality. Two meshes are compared geometrically. The difference with the previous version is that nodes(coordinates) are expected to be the same between this and other. In case of success 'cellCor' is informed. Warning ! These equivalences are CPU/Mem costly. The 3 values correspond respectively to policy used for cell comparison (see MEDCouplingUMesh::zipConnectivityTraducer to have more details)
 *   - 1 for fast 'equality'. This is a lazy level. Just number of cells and number of nodes are considered here and 3 cells (begin,middle,end)
 *   - 2 for deep 'equality' as 0 option except that no control is done on all strings in mesh.
 */
void MEDCouplingMesh::checkGeoEquivalWith(const MEDCouplingMesh *other, int levOfCheck, double prec,
                                          DataArrayInt *&cellCor, DataArrayInt *&nodeCor) const throw(INTERP_KERNEL::Exception)
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
 * Given a nodeIds range ['partBg','partEnd'), this method returns the set of cell ids in ascendant order whose connectivity of
 * these cells are fully included in the range. As a consequence the returned set of cell ids does \b not \b always fit the nodes in ['partBg','partEnd')
 * This method returns the corresponding cells in a newly created array that the caller has the responsability.
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
  ret->alloc(crest.size(),1);
  std::copy(crest.begin(),crest.end(),ret->getPointer());
  return ret;
}

/*!
 * This method checks fastly that 'this' and 'other' are equal. All common checks are done here.
 */
void MEDCouplingMesh::checkFastEquivalWith(const MEDCouplingMesh *other, double prec) const throw(INTERP_KERNEL::Exception)
{
  if(getMeshDimension()!=other->getMeshDimension())
    throw INTERP_KERNEL::Exception("checkFastEquivalWith : Mesh dimensions are not equal !");
  if(getSpaceDimension()!=other->getSpaceDimension())
    throw INTERP_KERNEL::Exception("checkFastEquivalWith : Space dimensions are not equal !");
  if(getNumberOfCells()!=other->getNumberOfCells())
    throw INTERP_KERNEL::Exception("checkFastEquivalWith : number of cells are not equal !");
}

/*!
 * This method is very poor and looks only if 'this' and 'other' are candidate for merge of fields lying repectively on them.
 */
bool MEDCouplingMesh::areCompatibleForMerge(const MEDCouplingMesh *other) const
{
  if(getMeshDimension()!=other->getMeshDimension())
    return false;
  if(getSpaceDimension()!=other->getSpaceDimension())
    return false;
  return true;
}

/*!
 * This method builds a field lying on 'this' with 'nbOfComp' components.
 * 'func' is a pointer that points to a function that takes 2 arrays in parameter and returns a boolean.
 * The first array is a in-param of size this->getSpaceDimension and the second an out param of size 'nbOfComp'.
 * The return field will have type specified by 't'. 't' is also used to determine where values of field will be
 * evaluate.
 * @param t type of field returned and specifies where the evaluation of func will be done.
 * @param nbOfComp number of components of returned field.
 * @param func pointer to a function that should return false if the evaluation failed. (division by 0. for example)
 * @return field with counter = 1.
 */
MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalytic(TypeOfField t, int nbOfComp, FunctionToEvaluate func) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(t,NO_TIME);
  ret->setMesh(this);
  ret->fillFromAnalytic(nbOfComp,func);
  ret->incrRef();
  return ret;
}

/*!
 * This method copyies all tiny strings from other (name and components name).
 * @throw if other and this have not same mesh type.
 */
void MEDCouplingMesh::copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception)
{
  _name=other->_name;
  _description=other->_description;
  _time_unit=other->_time_unit;
}

/*!
 * This method copies all attributes that are \b NOT arrays in this.
 * All tiny attributes not usefully for state of 'this' are ignored.
 */
void MEDCouplingMesh::copyTinyInfoFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception)
{
  copyTinyStringsFrom(other);
  _time=other->_time;
  _iteration=other->_iteration;
  _order=other->_order;
}

/*!
 * This method builds a field lying on 'this' with 'nbOfComp' components.
 * 'func' is a string that is the expression to evaluate.
 * The return field will have type specified by 't'. 't' is also used to determine where values of field will be
 * evaluate.
 * @param t type of field returned and specifies where the evaluation of func will be done.
 * @param nbOfComp number of components of returned field.
 * @param func expression.
 * @return field with counter = 1.
 */
MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalytic(TypeOfField t, int nbOfComp, const char *func) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(t,NO_TIME);
  ret->setMesh(this);
  ret->fillFromAnalytic(nbOfComp,func);
  ret->incrRef();
  return ret;
}

/*!
 * This method builds a field lying on 'this' with 'nbOfComp' components.
 * 'func' is a string that is the expression to evaluate.
 * The return field will have type specified by 't'. 't' is also used to determine where values of field will be
 * evaluate. This method is different than MEDCouplingMesh::fillFromAnalytic, because the info on components are used here to determine vars pos in 'func'.
 *
 * @param t type of field returned and specifies where the evaluation of func will be done.
 * @param nbOfComp number of components of returned field.
 * @param func expression.
 * @return field with counter = 1.
 */
MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalytic2(TypeOfField t, int nbOfComp, const char *func) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(t,NO_TIME);
  ret->setMesh(this);
  ret->fillFromAnalytic2(nbOfComp,func);
  ret->incrRef();
  return ret;
}

/*!
 * This method builds a field lying on 'this' with 'nbOfComp' components.
 * 'func' is a string that is the expression to evaluate.
 * The return field will have type specified by 't'. 't' is also used to determine where values of field will be
 * evaluate. This method is different than MEDCouplingMesh::fillFromAnalytic, because 'varsOrder' specifies the pos to assign of vars in 'func'.
 *
 * @param t type of field returned and specifies where the evaluation of func will be done.
 * @param nbOfComp number of components of returned field.
 * @param func expression.
 * @return field with counter = 1.
 */
MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalytic3(TypeOfField t, int nbOfComp, const std::vector<std::string>& varsOrder, const char *func) const
{
  MEDCouplingAutoRefCountObjectPtr<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(t,NO_TIME);
  ret->setMesh(this);
  ret->fillFromAnalytic3(nbOfComp,varsOrder,func);
  ret->incrRef();
  return ret;
}

/*!
 * retruns a newly created mesh with counter=1 
 * that is the union of mesh1 and mesh2 if possible. The cells of mesh2 will appear after cells of 'mesh1'. Idem for nodes.
 */
MEDCouplingMesh *MEDCouplingMesh::MergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2)
{
  return mesh1->mergeMyselfWith(mesh2);
}

void MEDCouplingMesh::getCellsContainingPoint(const double *pos, double eps, std::vector<int>& elts) const
{
  int ret=getCellContainingPoint(pos,eps);
  elts.push_back(ret);
}

void MEDCouplingMesh::getCellsContainingPoints(const double *pos, int nbOfPoints, double eps, std::vector<int>& elts, std::vector<int>& eltsIndex) const
{
  eltsIndex.resize(nbOfPoints+1);
  eltsIndex[0]=0;
  elts.clear();
  int spaceDim=getSpaceDimension();
  const double *work=pos;
  for(int i=0;i<nbOfPoints;i++,work+=spaceDim)
    {
      int ret=getCellContainingPoint(work,eps);
      if(ret>=0)
        {
          elts.push_back(ret);
          eltsIndex[i+1]=eltsIndex[i]+1;
        }
      else
        eltsIndex[i+1]=eltsIndex[i];
    }
}
