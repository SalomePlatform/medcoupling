//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingFieldDiscretization.hxx"

#include "InterpKernelExprParser.hxx"

#include <sstream>
#include <iterator>

using namespace ParaMEDMEM;

bool MEDCouplingMesh::areCompatible(const MEDCouplingMesh *other) const
{
  if(getMeshDimension()!=other->getMeshDimension())
    return false;
  if(getSpaceDimension()!=other->getMeshDimension())
    return false;
  return true;
}

MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalytic(TypeOfField t, int nbOfComp, FunctionToEvaluate func) const
{
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(t);
  ret->setMesh(this);
  DataArrayDouble *loc=ret->getDiscretization()->getLocalizationOfDiscValues(this);
  DataArrayDouble *array=DataArrayDouble::New();
  int nbOfTuple=loc->getNumberOfTuples();
  int nbCompIn=loc->getNumberOfComponents();
  const double *locPtr=loc->getConstPointer();
  array->alloc(nbOfTuple,nbOfComp);
  double *ptToFill=array->getPointer();
  for(int i=0;i<nbOfTuple;i++)
    {
      if(!func(locPtr+nbCompIn*i,ptToFill))
        {
          std::ostringstream oss; oss << "For tuple # " << i << " localized on (";
          std::copy(locPtr+nbCompIn*i,locPtr+nbCompIn*(i+1),std::ostream_iterator<double>(oss,", "));
          oss << ") : Evaluation of function failed !";
          loc->decrRef();
          array->decrRef();
          ret->decrRef();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      ptToFill+=nbOfComp;
    }
  loc->decrRef();
  ret->setArray(array);
  array->decrRef();
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingMesh::fillFromAnalytic(TypeOfField t, int nbOfComp, const char *func) const
{
  INTERP_KERNEL::ExprParser expr(func);
  expr.parse();
  std::set<std::string> vars;
  expr.getTrueSetOfVars(vars);
  if(vars.size()>getSpaceDimension())
    {
      std::ostringstream oss; oss << "The mesh has a spaceDim==" << getSpaceDimension() << " and there are ";
      oss << vars.size() << " variables : ";
      std::copy(vars.begin(),vars.end(),std::ostream_iterator<std::string>(oss," "));
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::vector<std::string> varsV(vars.begin(),vars.end());
  expr.prepareExprEvaluation(varsV);
  //
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(t);
  ret->setMesh(this);
  DataArrayDouble *loc=ret->getDiscretization()->getLocalizationOfDiscValues(this);
  DataArrayDouble *array=DataArrayDouble::New();
  int nbOfTuple=loc->getNumberOfTuples();
  int nbCompIn=loc->getNumberOfComponents();
  const double *locPtr=loc->getConstPointer();
  array->alloc(nbOfTuple,nbOfComp);
  double *ptToFill=array->getPointer();
  for(int i=0;i<nbOfTuple;i++)
    {
      try
        {
          expr.evaluateExpr(nbOfComp,ptToFill,locPtr+nbCompIn*i);
        }
      catch(INTERP_KERNEL::Exception& e)
        {
          std::ostringstream oss; oss << "For tuple # " << i << " localized on (";
          std::copy(locPtr+nbCompIn*i,locPtr+nbCompIn*(i+1),std::ostream_iterator<double>(oss,", "));
          oss << ") : Evaluation of function failed ! " << e.what();
          loc->decrRef();
          array->decrRef();
          ret->decrRef();
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      ptToFill+=nbOfComp;
    }
  loc->decrRef();
  ret->setArray(array);
  array->decrRef();
  return ret;
}

MEDCouplingMesh *MEDCouplingMesh::mergeMeshes(const MEDCouplingMesh *mesh1, const MEDCouplingMesh *mesh2)
{
  return mesh1->mergeMyselfWith(mesh2);
}
