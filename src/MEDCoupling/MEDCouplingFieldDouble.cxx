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
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMesh.hxx"
#include "MEDCouplingTimeDiscretization.hxx"
#include "MEDCouplingFieldDiscretization.hxx"

#include <sstream>

using namespace ParaMEDMEM;

MEDCouplingFieldDouble *MEDCouplingFieldDouble::New(TypeOfField type, TypeOfTimeDiscretization td)
{
  return new MEDCouplingFieldDouble(type,td);
}

MEDCouplingFieldDouble *MEDCouplingFieldDouble::clone(bool recDeepCpy) const
{
  return new MEDCouplingFieldDouble(*this,recDeepCpy);
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(TypeOfField type, TypeOfTimeDiscretization td):MEDCouplingField(type),
                                                                                              _time_discr(MEDCouplingTimeDiscretization::New(td))
{
}

MEDCouplingFieldDouble::MEDCouplingFieldDouble(const MEDCouplingFieldDouble& other, bool deepCpy):MEDCouplingField(other),
                                                                                                  _time_discr(other._time_discr->performCpy(deepCpy))
{
}

MEDCouplingFieldDouble::~MEDCouplingFieldDouble()
{
  delete _time_discr;
}

void MEDCouplingFieldDouble::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Field invalid because no mesh specified !");
  if(!getArray())
    throw INTERP_KERNEL::Exception("Field invalid because no values set !");
  _type->checkCoherencyBetween(_mesh,getArray());
}

double MEDCouplingFieldDouble::accumulate(int compId) const
{
  const double *ptr=getArray()->getConstPointer();
  int nbTuple=getArray()->getNumberOfTuples();
  int nbComps=getArray()->getNumberOfComponents();
  double ret=0.;
  for(int i=0;i<nbTuple;i++)
    ret+=ptr[i*nbComps+compId];
  return ret;
}

double MEDCouplingFieldDouble::measureAccumulate(int compId) const
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("No mesh underlying this field to perform measureAccumulate");
  MEDCouplingFieldDouble *weight=_type->getWeightingField(_mesh);
  const double *ptr=weight->getArray()->getConstPointer();
  int nbOfValues=weight->getArray()->getNbOfElems();
  double ret=0.;
  for (int i=0; i<nbOfValues; i++)
    ret+=getIJ(i,compId)*ptr[i];
  weight->decrRef();
  return ret;
}

void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double *res) const throw(INTERP_KERNEL::Exception)
{
  _time_discr->checkNoTimePresence();
}

void MEDCouplingFieldDouble::getValueOn(const double *spaceLoc, double time, double *res) const throw(INTERP_KERNEL::Exception)
{
  _time_discr->checkTimePresence(time);
}

void MEDCouplingFieldDouble::applyLin(double a, double b, int compoId)
{
  double *ptr=getArray()->getPointer();
  ptr+=compoId;
  int nbOfComp=getArray()->getNumberOfComponents();
  int nbOfTuple=getArray()->getNumberOfTuples();
  for(int i=0;i<nbOfTuple;i++,ptr+=nbOfComp)
    *ptr=a*(*ptr)+b;
}

int MEDCouplingFieldDouble::getNumberOfComponents() const
{
  return getArray()->getNumberOfComponents();
}

int MEDCouplingFieldDouble::getNumberOfTuples() const throw(INTERP_KERNEL::Exception)
{
  if(!_mesh)
    throw INTERP_KERNEL::Exception("Impossible to retrieve number of tuples because no mesh specified !");
  return _type->getNumberOfTuples(_mesh);
}

void MEDCouplingFieldDouble::updateTime()
{
  MEDCouplingField::updateTime();
  if(getArray())
    updateTimeWith(*getArray());
}

void MEDCouplingFieldDouble::setArray(DataArrayDouble *array)
{
  _time_discr->setArray(array,this);
}
