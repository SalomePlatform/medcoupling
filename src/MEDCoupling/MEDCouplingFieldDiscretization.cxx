// Copyright (C) 2007-2023  CEA/DEN, EDF R&D
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

#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingFieldDiscretizationOnNodesFE.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MCAuto.hxx"

#include "CellModel.hxx"
#include "InterpolationUtils.hxx"
#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGaussCoords.hxx"
#include "InterpKernelMatrixTools.hxx"
#include "InterpKernelDenseMatrix.hxx"

#include <set>
#include <list>
#include <limits>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <functional>

using namespace MEDCoupling;

const double MEDCouplingFieldDiscretization::DFLT_PRECISION=1.e-12;

const char MEDCouplingFieldDiscretizationP0::REPR[]="P0";

const char MEDCouplingFieldDiscretizationP1::REPR[]="P1";

const mcIdType MEDCouplingFieldDiscretizationPerCell::DFT_INVALID_LOCID_VALUE=-1;

const char MEDCouplingFieldDiscretizationGauss::REPR[]="GAUSS";

const char MEDCouplingFieldDiscretizationGaussNE::REPR[]="GSSNE";

const char MEDCouplingFieldDiscretizationKriging::REPR[]="KRIGING";

// doc is here http://www.code-aster.org/V2/doc/default/fr/man_r/r3/r3.01.01.pdf
const double MEDCouplingFieldDiscretizationGaussNE::FGP_POINT1[1]={0.};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_SEG2[2]={1.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_SEG3[3]={0.5555555555555556,0.8888888888888888,0.5555555555555556};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_SEG4[4]={0.347854845137454,0.347854845137454,0.652145154862546,0.652145154862546};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_TRI3[3]={0.16666666666666666,0.16666666666666666,0.16666666666666666};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_TRI6[6]={0.0549758718227661,0.0549758718227661,0.0549758718227661,0.11169079483905,0.11169079483905,0.11169079483905};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_TRI7[7]={0.062969590272413,0.062969590272413,0.062969590272413,0.066197076394253,0.066197076394253,0.066197076394253,0.1125};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_QUAD4[4]={1.,1.,1.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_QUAD8[8]={1.,1.,1.,1.,1.,1.,1.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_QUAD9[9]={0.30864197530864196,0.30864197530864196,0.30864197530864196,0.30864197530864196,0.49382716049382713,0.49382716049382713,0.49382716049382713,0.49382716049382713,0.7901234567901234};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_TETRA4[4]={0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_TETRA10[10]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};//to check
const double MEDCouplingFieldDiscretizationGaussNE::FGP_PENTA6[6]={0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_PENTA15[15]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};//to check
const double MEDCouplingFieldDiscretizationGaussNE::FGP_PENTA18[18]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};//to check
const double MEDCouplingFieldDiscretizationGaussNE::FGP_HEXA8[8]={1.,1.,1.,1.,1.,1.,1.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_HEXA20[20]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_HEXA27[27]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_PYRA5[5]={0.13333333333333333,0.13333333333333333,0.13333333333333333,0.13333333333333333,0.13333333333333333};
const double MEDCouplingFieldDiscretizationGaussNE::FGP_PYRA13[13]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};//to check
const double MEDCouplingFieldDiscretizationGaussNE::REF_SEG2[2]={-1.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_SEG3[3]={-1.,1.,0.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_SEG4[4]={-1.,1.,-0.3333333333333333,0.3333333333333333};
const double MEDCouplingFieldDiscretizationGaussNE::REF_TRI3[6]={0.,0.,1.,0.,0.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_TRI6[12]={0.,0.,1.,0.,0.,1.,0.5,0.,0.5,0.5,0.,0.5};
const double MEDCouplingFieldDiscretizationGaussNE::REF_TRI7[14]={0.,0.,1.,0.,0.,1.,0.5,0.,0.5,0.5,0.,0.5,0.3333333333333333,0.3333333333333333};
const double MEDCouplingFieldDiscretizationGaussNE::REF_QUAD4[8]={-1.,-1.,1.,-1.,1.,1.,-1.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_QUAD8[16]={-1.,-1.,1.,-1.,1.,1.,-1.,1.,0.,-1.,1.,0.,0.,1.,-1.,0.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_QUAD9[18]={-1.,-1.,1.,-1.,1.,1.,-1.,1.,0.,-1.,1.,0.,0.,1.,-1.,0.,0.,0.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_TETRA4[12]={0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_TETRA10[30]={0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,0.5,0.5,0.,0.,0.5,0.,0.5,0.,0.5,0.5,0.,0.5,0.,0.5,0.5,0.,0.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_PENTA6[18]={-1.,1.,0.,-1.,0.,1.,-1.,0.,0.,1.,1.,0.,1.,0.,1.,1.,0.,0.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_PENTA15[45]={-1.,1.,0.,-1.,0.,1.,-1.,0.,0.,1.,1.,0.,1.,0.,1.,1.,0.,0.,-1.,0.5,0.5,-1.,0.,0.5,-1.,0.5,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,0.5,0.5,1.,0.,0.5,1.,0.5,0.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_PENTA18[54]={-1.,1.,0.,-1.,0.,1.,-1.,0.,0.,1.,1.,0.,1.,0.,1.,1.,0.,0.,-1.,0.5,0.5,-1.,0.,0.5,-1.,0.5,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,0.5,0.5,1.,0.,0.5,1.,0.5,0.,0.,0.5,0.5,0.,0.,0.5,0.,0.5,0.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_HEXA8[24]={-1.,-1.,-1.,1.,-1.,-1.,1.,1.,-1.,-1.,1.,-1.,-1.,-1.,1.,1.,-1.,1.,1.,1.,1.,-1.,1.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_HEXA20[60]={-1.,-1.,-1.,1.,-1.,-1.,1.,1.,-1.,-1.,1.,-1.,-1.,-1.,1.,1.,-1.,1.,1.,1.,1.,-1.,1.,1.,0.,-1.,-1.,1.,0.,-1.,0.,1.,-1.,-1.,0.,-1.,-1.,-1.,0.,1.,-1.,0.,1.,1.,0.,-1.,1.,0.,0.,-1.,1.,1.,0.,1.,0.,1.,1.,-1.,0.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_HEXA27[81]={-1.,-1.,-1.,-1.,1.,-1.,1.,1.,-1.,1.,-1.,-1.,-1.,-1.,1.,-1.,1.,1.,1.,1.,1.,1.,-1.,1.,-1.,0.,-1.,0.,1.,-1.,1.,0.,-1.,0.,-1.,-1.,-1.,0.,1.,0.,1.,1.,1.,0.,1.,0.,-1.,1.,-1.,-1.,0.,-1.,1.,0.,1.,1.,0.,1.,-1.,0.,0.,0.,-1.,-1.,0.,0.,0.,1.,0.,1.,0.,0.,0.,-1.,0.,0.,0.,1.,0.,0.,0.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_PYRA5[15]={1.,0.,0.,0.,1.,0.,-1.,0.,0.,0.,-1.,0.,0.,0.,1.};
const double MEDCouplingFieldDiscretizationGaussNE::REF_PYRA13[39]={1.,0.,0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,0.,0.,1.,0.5,-0.5,0.,-0.5,-0.5,0.,-0.5,0.5,0.,0.5,0.5,0.,0.5,0.,0.5,0.,-0.5,0.5,-0.5,0.,0.5,0.,0.5,0.5};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_SEG2[2]={0.577350269189626,-0.577350269189626};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_SEG3[3]={-0.774596669241,0.,0.774596669241};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_SEG4[4]={0.339981043584856,-0.339981043584856,0.861136311594053,-0.861136311594053};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_TRI3[6]={0.16666666666666667,0.16666666666666667,0.6666666666666667,0.16666666666666667,0.16666666666666667,0.6666666666666667};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_TRI6[12]={0.091576213509771,0.091576213509771,0.816847572980458,0.091576213509771,0.091576213509771,0.816847572980458,0.445948490915965,0.10810301816807,0.445948490915965,0.445948490915965,0.10810301816807,0.445948490915965};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_TRI7[14]={0.3333333333333333,0.3333333333333333,0.470142064105115,0.470142064105115,0.05971587178977,0.470142064105115,0.470142064105115,0.05971587178977,0.101286507323456,0.101286507323456,0.797426985353088,0.101286507323456,0.101286507323456,0.797426985353088};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_QUAD4[8]={-0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,0.774596669241483,0.774596669241483,-0.774596669241483,0.774596669241483};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_QUAD8[16]={-0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,0.774596669241483,0.774596669241483,-0.774596669241483,0.774596669241483,0.,-0.774596669241483,0.774596669241483,0.,0.,0.774596669241483,-0.774596669241483,0.};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_QUAD9[18]={-0.774596669241483,-0.774596669241483,0.774596669241483,-0.774596669241483,0.774596669241483,0.774596669241483,-0.774596669241483,0.774596669241483,0.,-0.774596669241483,0.774596669241483,0.,0.,0.774596669241483,-0.774596669241483,0.,0.,0.};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_TETRA4[12]={0.1381966011250105,0.1381966011250105,0.1381966011250105,0.1381966011250105,0.1381966011250105,0.5854101966249685,0.1381966011250105,0.5854101966249685,0.1381966011250105,0.5854101966249685,0.1381966011250105,0.1381966011250105};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_TETRA10[30]={0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,0.5,0.5,0.,0.,0.5,0.,0.5,0.,0.5,0.5,0.,0.5,0.,0.5,0.5,0.,0.};//to check
const double MEDCouplingFieldDiscretizationGaussNE::LOC_PENTA6[18]={-0.5773502691896258,0.5,0.5,-0.5773502691896258,0.,0.5,-0.5773502691896258,0.5,0.,0.5773502691896258,0.5,0.5,0.5773502691896258,0.,0.5,0.5773502691896258,0.5,0.};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_PENTA15[45]={-1.,1.,0.,-1.,0.,1.,-1.,0.,0.,1.,1.,0.,1.,0.,1.,1.,0.,0.,-1.,0.5,0.5,-1.,0.,0.5,-1.,0.5,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,0.5,0.5,1.,0.,0.5,1.,0.5,0.};//to check
const double MEDCouplingFieldDiscretizationGaussNE::LOC_PENTA18[54]={-1.,1.,0.,-1.,0.,1.,-1.,0.,0.,1.,1.,0.,1.,0.,1.,1.,0.,0.,-1.,0.5,0.5,-1.,0.,0.5,-1.,0.5,0.,0.,1.,0.,0.,0.,1.,0.,0.,0.,1.,0.5,0.5,1.,0.,0.5,1.,0.5,0.,0.,0.5,0.5,0.,0.,0.5,0.,0.5,0.};//to check
const double MEDCouplingFieldDiscretizationGaussNE::LOC_HEXA8[24]={-0.5773502691896258,-0.5773502691896258,-0.5773502691896258,-0.5773502691896258,-0.5773502691896258,0.5773502691896258,-0.5773502691896258,0.5773502691896258,-0.5773502691896258,-0.5773502691896258,0.5773502691896258,0.5773502691896258,0.5773502691896258,-0.5773502691896258,-0.5773502691896258,0.5773502691896258,-0.5773502691896258,0.5773502691896258,0.5773502691896258,0.5773502691896258,-0.5773502691896258,0.5773502691896258,0.5773502691896258,0.5773502691896258};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_HEXA20[60]={-1.,-1.,-1.,1.,-1.,-1.,1.,1.,-1.,-1.,1.,-1.,-1.,-1.,1.,1.,-1.,1.,1.,1.,1.,-1.,1.,1.,0.,-1.,-1.,1.,0.,-1.,0.,1.,-1.,-1.,0.,-1.,-1.,-1.,0.,1.,-1.,0.,1.,1.,0.,-1.,1.,0.,0.,-1.,1.,1.,0.,1.,0.,1.,1.,-1.,0.,1.};//to check
const double MEDCouplingFieldDiscretizationGaussNE::LOC_HEXA27[81]={-1.,-1.,-1.,-1.,1.,-1.,1.,1.,-1.,1.,-1.,-1.,-1.,-1.,1.,-1.,1.,1.,1.,1.,1.,1.,-1.,1.,-1.,0.,-1.,0.,1.,-1.,1.,0.,-1.,0.,-1.,-1.,-1.,0.,1.,0.,1.,1.,1.,0.,1.,0.,-1.,1.,-1.,-1.,0.,-1.,1.,0.,1.,1.,0.,1.,-1.,0.,0.,0.,-1.,-1.,0.,0.,0.,1.,0.,1.,0.,0.,0.,-1.,0.,0.,0.,1.,0.,0.,0.};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_PYRA5[15]={0.5,0.,0.1531754163448146,0.,0.5,0.1531754163448146,-0.5,0.,0.1531754163448146,0.,-0.5,0.1531754163448146,0.,0.,0.6372983346207416};
const double MEDCouplingFieldDiscretizationGaussNE::LOC_PYRA13[39]={1.,0.,0.,0.,-1.,0.,-1.,0.,0.,0.,1.,0.,0.,0.,0.999999999999,0.5,-0.5,0.,-0.5,-0.5,0.,-0.5,0.5,0.,0.5,0.5,0.,0.5,0.,0.5,0.,-0.5,0.5,-0.5,0.,0.5,0.,0.5,0.5};//to check 0.99999... to avoid nan ! on node #4 of PYRA13

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
    case MEDCouplingFieldDiscretizationOnNodesFE::TYPE:
      return new MEDCouplingFieldDiscretizationOnNodesFE;
    default:
      throw INTERP_KERNEL::Exception("Chosen discretization is not implemented yet.");
  }
}

TypeOfField MEDCouplingFieldDiscretization::GetTypeOfFieldFromStringRepr(const std::string& repr)
{
  if(repr==MEDCouplingFieldDiscretizationP0::REPR)
    return MEDCouplingFieldDiscretizationP0::TYPE;
  if(repr==MEDCouplingFieldDiscretizationP1::REPR)
    return MEDCouplingFieldDiscretizationP1::TYPE;
  if(repr==MEDCouplingFieldDiscretizationGauss::REPR)
    return MEDCouplingFieldDiscretizationGauss::TYPE;
  if(repr==MEDCouplingFieldDiscretizationGaussNE::REPR)
    return MEDCouplingFieldDiscretizationGaussNE::TYPE;
  if(repr==MEDCouplingFieldDiscretizationKriging::REPR)
    return MEDCouplingFieldDiscretizationKriging::TYPE;
  if(repr==MEDCouplingFieldDiscretizationOnNodesFE::REPR)
    return MEDCouplingFieldDiscretizationOnNodesFE::TYPE;
  throw INTERP_KERNEL::Exception("Representation does not match with any field discretization !");
}

std::string MEDCouplingFieldDiscretization::GetTypeOfFieldRepr(TypeOfField type)
{
  switch(type)
  {
    case MEDCouplingFieldDiscretizationP0::TYPE:
      return MEDCouplingFieldDiscretizationP0::REPR;
    case MEDCouplingFieldDiscretizationP1::TYPE:
      return MEDCouplingFieldDiscretizationP1::REPR;
    case MEDCouplingFieldDiscretizationGauss::TYPE:
      return MEDCouplingFieldDiscretizationGauss::REPR;
    case MEDCouplingFieldDiscretizationGaussNE::TYPE:
      return MEDCouplingFieldDiscretizationGaussNE::REPR;
    case MEDCouplingFieldDiscretizationKriging::TYPE:
      return MEDCouplingFieldDiscretizationKriging::REPR;
    case MEDCouplingFieldDiscretizationOnNodesFE::TYPE:
      return MEDCouplingFieldDiscretizationOnNodesFE::REPR;
    default:
      throw INTERP_KERNEL::Exception("GetTypeOfFieldRepr : Representation does not match with any field discretization !");
  }
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
 * This method is an alias of MEDCouplingFieldDiscretization::clone. It is only here for coherency with all the remaining of MEDCoupling.
 * \sa MEDCouplingFieldDiscretization::clone.
 */
MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretization::deepCopy() const
{
  return clone();
}

/*!
 * For all field discretization excepted GaussPts the [ \a startCellIds, \a endCellIds ) has no impact on the cloned instance.
 */
MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretization::clonePart(const mcIdType *startCellIds, const mcIdType *endCellIds) const
{
  return clone();
}

/*!
 * For all field discretization excepted GaussPts the slice( \a beginCellId, \a endCellIds, \a stepCellId ) has no impact on the cloned instance.
 */
MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretization::clonePartRange(mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds) const
{
  return clone();
}

/*!
 * Excepted for MEDCouplingFieldDiscretizationPerCell no underlying TimeLabel object : nothing to do in generally.
 */
void MEDCouplingFieldDiscretization::updateTime() const
{
}

std::size_t MEDCouplingFieldDiscretization::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> MEDCouplingFieldDiscretization::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

/*!
 * Computes normL1 of DataArrayDouble instance arr.
 * @param res output parameter expected to be of size arr->getNumberOfComponents();
 * @throw when the field discretization fails on getMeasure fields (gauss points for example)
 */
void MEDCouplingFieldDiscretization::normL1(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, double *res) const
{
  MCAuto<MEDCouplingFieldDouble> vol=getMeasureField(mesh,true);
  std::size_t nbOfCompo=arr->getNumberOfComponents();
  mcIdType nbOfElems=getNumberOfTuples(mesh);
  std::fill(res,res+nbOfCompo,0.);
  const double *arrPtr=arr->getConstPointer();
  const double *volPtr=vol->getArray()->getConstPointer();
  double deno=0.;
  for(mcIdType i=0;i<nbOfElems;i++)
    {
      double v=fabs(volPtr[i]);
      for(std::size_t j=0;j<nbOfCompo;j++)
        res[j]+=fabs(arrPtr[i*nbOfCompo+j])*v;
      deno+=v;
    }
  std::transform(res,res+nbOfCompo,res,std::bind(std::multiplies<double>(),std::placeholders::_1,1./deno));
}

/*!
 * Computes normL2 of DataArrayDouble instance arr.
 * @param res output parameter expected to be of size arr->getNumberOfComponents();
 * @throw when the field discretization fails on getMeasure fields (gauss points for example)
 */
void MEDCouplingFieldDiscretization::normL2(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, double *res) const
{
  MCAuto<MEDCouplingFieldDouble> vol=getMeasureField(mesh,true);
  std::size_t nbOfCompo=arr->getNumberOfComponents();
  mcIdType nbOfElems=getNumberOfTuples(mesh);
  std::fill(res,res+nbOfCompo,0.);
  const double *arrPtr=arr->getConstPointer();
  const double *volPtr=vol->getArray()->getConstPointer();
  double deno=0.;
  for(mcIdType i=0;i<nbOfElems;i++)
    {
      double v=fabs(volPtr[i]);
      for(std::size_t j=0;j<nbOfCompo;j++)
        res[j]+=arrPtr[i*nbOfCompo+j]*arrPtr[i*nbOfCompo+j]*v;
      deno+=v;
    }
  std::transform(res,res+nbOfCompo,res,std::bind(std::multiplies<double>(),std::placeholders::_1,1./deno));
  std::transform(res,res+nbOfCompo,res,[](double c){return sqrt(c);});
}

/*!
 * Computes integral of DataArrayDouble instance arr.
 * @param res output parameter expected to be of size arr->getNumberOfComponents();
 * @throw when the field discretization fails on getMeasure fields (gauss points for example)
 */
void MEDCouplingFieldDiscretization::integral(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, bool isWAbs, double *res) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretization::integral : mesh is NULL !");
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretization::integral : input array is NULL !");
  MCAuto<MEDCouplingFieldDouble> vol=getMeasureField(mesh,isWAbs);
  std::size_t nbOfCompo(arr->getNumberOfComponents());
  mcIdType nbOfElems(getNumberOfTuples(mesh));
  if(nbOfElems!=arr->getNumberOfTuples())
    {
      std::ostringstream oss; oss << "MEDCouplingFieldDiscretization::integral : field is not correct ! number of tuples in array is " << arr->getNumberOfTuples();
      oss << " whereas number of tuples expected is " << nbOfElems << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  std::fill(res,res+nbOfCompo,0.);
  const double *arrPtr(arr->begin()),*volPtr(vol->getArray()->begin());
  INTERP_KERNEL::AutoPtr<double> tmp=new double[nbOfCompo];
  for(mcIdType i=0;i<nbOfElems;i++)
    {
      std::transform(arrPtr+i*nbOfCompo,arrPtr+(i+1)*nbOfCompo,(double *)tmp,std::bind(std::multiplies<double>(),std::placeholders::_1,volPtr[i]));
      std::transform((double *)tmp,(double *)tmp+nbOfCompo,res,res,std::plus<double>());
    }
}

/*!
 * This method is strictly equivalent to MEDCouplingFieldDiscretization::buildSubMeshData except that it is optimized for input defined as a range of cell ids.
 * 
 * \param [out] beginOut Valid only if \a di is NULL
 * \param [out] endOut Valid only if \a di is NULL
 * \param [out] stepOut Valid only if \a di is NULL
 * \param [out] di is an array returned that specifies entity ids (nodes, cells, Gauss points... ) in array if no output range is foundable.
 *
 * \sa MEDCouplingFieldDiscretization::buildSubMeshData
 */
MEDCouplingMesh *MEDCouplingFieldDiscretization::buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const
{
  MCAuto<DataArrayIdType> da=DataArrayIdType::Range(beginCellIds,endCellIds,stepCellIds);
  return buildSubMeshData(mesh,da->begin(),da->end(),di);
}

void MEDCouplingFieldDiscretization::getSerializationIntArray(DataArrayIdType *& arr) const
{
  arr=0;
}

/*!
 * Empty : Not a bug
 */
void MEDCouplingFieldDiscretization::getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const
{
}

/*!
 * Empty : Not a bug
 */
void MEDCouplingFieldDiscretization::getTinySerializationDbleInformation(std::vector<double>& tinyInfo) const
{
}

void MEDCouplingFieldDiscretization::resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *& arr)
{
  arr=0;
}

/*!
 * Empty : Not a bug
 */
void MEDCouplingFieldDiscretization::checkForUnserialization(const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *arr)
{
}

/*!
 * Empty : Not a bug
 */
void MEDCouplingFieldDiscretization::finishUnserialization(const std::vector<double>& tinyInfo)
{
}

/*!
 * This method is typically the first step of renumbering. The implementation is empty it is not a bug only gauss is impacted
 * virtually by this method.
 */
void MEDCouplingFieldDiscretization::renumberCells(const mcIdType *old2NewBg, bool check)
{
}

double MEDCouplingFieldDiscretization::getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, mcIdType cellId, mcIdType nodeIdInCell, int compoId) const
{
  throw INTERP_KERNEL::Exception("getIJK Invalid ! only for GaussPoint and GaussNE discretizations !");
}

void MEDCouplingFieldDiscretization::setGaussLocalizationOnType(const MEDCouplingMesh *m, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                                const std::vector<double>& gsCoo, const std::vector<double>& wg)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

void MEDCouplingFieldDiscretization::setGaussLocalizationOnCells(const MEDCouplingMesh *m, const mcIdType *begin, const mcIdType *end, const std::vector<double>& refCoo,
                                                                 const std::vector<double>& gsCoo, const std::vector<double>& wg)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

void MEDCouplingFieldDiscretization::clearGaussLocalizations()
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

MEDCouplingGaussLocalization& MEDCouplingFieldDiscretization::getGaussLocalization(mcIdType locId)
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

const MEDCouplingGaussLocalization& MEDCouplingFieldDiscretization::getGaussLocalization(mcIdType locId) const
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

mcIdType MEDCouplingFieldDiscretization::getNbOfGaussLocalization() const
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

mcIdType MEDCouplingFieldDiscretization::getGaussLocalizationIdOfOneCell(mcIdType cellId) const
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

mcIdType MEDCouplingFieldDiscretization::getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

std::set<mcIdType> MEDCouplingFieldDiscretization::getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

void MEDCouplingFieldDiscretization::getCellIdsHavingGaussLocalization(mcIdType locId, std::vector<mcIdType>& cellIds) const
{
  throw INTERP_KERNEL::Exception("Invalid method for the corresponding field discretization : available only for GaussPoint discretization !");
}

void MEDCouplingFieldDiscretization::RenumberEntitiesFromO2NArr(double eps, const mcIdType *old2NewPtr, mcIdType newNbOfEntity, DataArrayDouble *arr, const std::string& msg)
{
  if(!arr)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretization::RenumberEntitiesFromO2NArr : input array is NULL !");
  mcIdType oldNbOfElems=arr->getNumberOfTuples();
  std::size_t nbOfComp=arr->getNumberOfComponents();
  mcIdType newNbOfTuples=newNbOfEntity;
  MCAuto<DataArrayDouble> arrCpy=arr->deepCopy();
  const double *ptSrc=arrCpy->getConstPointer();
  arr->reAlloc(newNbOfTuples);
  double *ptToFill=arr->getPointer();
  std::fill(ptToFill,ptToFill+nbOfComp*newNbOfTuples,std::numeric_limits<double>::max());
  INTERP_KERNEL::AutoPtr<double> tmp=new double[nbOfComp];
  for(mcIdType i=0;i<oldNbOfElems;i++)
    {
      mcIdType newNb=old2NewPtr[i];
      if(newNb>=0)//if newNb<0 the node is considered as out.
        {
          if(std::find_if(ptToFill+newNb*nbOfComp,ptToFill+(newNb+1)*nbOfComp,std::bind(std::not_equal_to<double>(),std::placeholders::_1,std::numeric_limits<double>::max()))
          ==ptToFill+(newNb+1)*nbOfComp)
            std::copy(ptSrc+i*nbOfComp,ptSrc+(i+1)*nbOfComp,ptToFill+newNb*nbOfComp);
          else
            {
              std::transform(ptSrc+i*nbOfComp,ptSrc+(i+1)*nbOfComp,ptToFill+newNb*nbOfComp,(double *)tmp,std::minus<double>());
              std::transform((double *)tmp,((double *)tmp)+nbOfComp,(double *)tmp,[](double c){return fabs(c);});
              //if(!std::equal(ptSrc+i*nbOfComp,ptSrc+(i+1)*nbOfComp,ptToFill+newNb*nbOfComp))
              if(*std::max_element((double *)tmp,((double *)tmp)+nbOfComp)>eps)
                {
                  std::ostringstream oss;
                  oss << msg << " " << i << " and " << std::find(old2NewPtr,old2NewPtr+i,newNb)-old2NewPtr
                      << " have been merged and " << msg << " field on them are different !";
                  throw INTERP_KERNEL::Exception(oss.str().c_str());
                }
            }
        }
    }
}

void MEDCouplingFieldDiscretization::RenumberEntitiesFromN2OArr(const mcIdType *new2OldPtr, mcIdType new2OldSz, DataArrayDouble *arr, const std::string& msg)
{
  std::size_t nbOfComp=arr->getNumberOfComponents();
  MCAuto<DataArrayDouble> arrCpy=arr->deepCopy();
  const double *ptSrc=arrCpy->getConstPointer();
  arr->reAlloc(new2OldSz);
  double *ptToFill=arr->getPointer();
  for(mcIdType i=0;i<new2OldSz;i++)
    {
      mcIdType oldNb=new2OldPtr[i];
      std::copy(ptSrc+oldNb*nbOfComp,ptSrc+(oldNb+1)*nbOfComp,ptToFill+i*nbOfComp);
    }
}

MEDCouplingFieldDiscretization::~MEDCouplingFieldDiscretization()
{
}

TypeOfField MEDCouplingFieldDiscretizationP0::getEnum() const
{
  return TYPE;
}

/*!
 * This method is simply called by MEDCouplingFieldDiscretization::deepCopy. It performs the deep copy of \a this.
 *
 * \sa MEDCouplingFieldDiscretization::deepCopy.
 */
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
  if(!other)
    {
      reason="other spatial discretization is NULL, and this spatial discretization (P0) is defined.";
      return false;
    }
  const MEDCouplingFieldDiscretizationP0 *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationP0 *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason="Spatial discrtization of this is ON_CELLS, which is not the case of other.";
  return ret;
}

mcIdType MEDCouplingFieldDiscretizationP0::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::getNumberOfTuples : NULL input mesh !");
  return mesh->getNumberOfCells();
}

/*!
 * This method returns the number of tuples regarding exclusively the input code \b without \b using \b a \b mesh \b in \b input.
 * The input code coherency is also checked regarding spatial discretization of \a this.
 * If an incoherency is detected, an exception will be thrown. If the input code is coherent, the number of tuples expected is returned.
 * The number of tuples expected is equal to those to have a valid field lying on \a this and having a mesh fitting perfectly the input code (geometric type distribution).
 */
mcIdType MEDCouplingFieldDiscretizationP0::getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const
{
  if(code.size()%3!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::getNumberOfTuplesExpectedRegardingCode : invalid input code !");
  mcIdType nbOfSplit=ToIdType(idsPerType.size());
  mcIdType nbOfTypes=ToIdType(code.size()/3);
  mcIdType ret=0;
  for(mcIdType i=0;i<nbOfTypes;i++)
    {
      mcIdType nbOfEltInChunk=code[3*i+1];
      if(nbOfEltInChunk<0)
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::getNumberOfTuplesExpectedRegardingCode : invalid input code ! presence of negative value in a type !");
      mcIdType pos=code[3*i+2];
      if(pos!=-1)
        {
          if(pos<0 || pos>=nbOfSplit)
            {
              std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationP0::getNumberOfTuplesExpectedRegardingCode : input code points to pos " << pos << " in typeid " << i << " ! Should be in [0," << nbOfSplit << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          const DataArrayIdType *ids(idsPerType[pos]);
          if(!ids || !ids->isAllocated() || ids->getNumberOfComponents()!=1 || ids->getNumberOfTuples()!=nbOfEltInChunk || ids->getMinValueInArray()<0)
            {
              std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationP0::getNumberOfTuplesExpectedRegardingCode : input pfl chunck at pos " << pos << " should have " << i << " tuples and one component and with ids all >=0 !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      ret+=nbOfEltInChunk;
    }
  return ret;
}

mcIdType MEDCouplingFieldDiscretizationP0::getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::getNumberOfMeshPlaces : NULL input mesh !");
  return mesh->getNumberOfCells();
}

DataArrayIdType *MEDCouplingFieldDiscretizationP0::getOffsetArr(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::getOffsetArr : NULL input mesh !");
  std::size_t nbOfTuples=mesh->getNumberOfCells();
  DataArrayIdType *ret=DataArrayIdType::New();
  ret->alloc(nbOfTuples+1,1);
  ret->iota(0);
  return ret;
}

void MEDCouplingFieldDiscretizationP0::renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArray *>& arrays,
                                                             const mcIdType *old2NewBg, bool check)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::renumberArraysForCell : NULL input mesh !");
  const mcIdType *array=old2NewBg;
  if(check)
    array=DataArrayIdType::CheckAndPreparePermutation(old2NewBg,old2NewBg+mesh->getNumberOfCells());
  for(std::vector<DataArray *>::const_iterator it=arrays.begin();it!=arrays.end();it++)
    {
      if(*it)
        (*it)->renumberInPlace(array);
    }
  if(check)
    free(const_cast<mcIdType *>(array));
}

DataArrayDouble *MEDCouplingFieldDiscretizationP0::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::getLocalizationOfDiscValues : NULL input mesh !");
  return mesh->computeCellCenterOfMass();
}

void MEDCouplingFieldDiscretizationP0::computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const mcIdType *tupleIdsBg, const mcIdType *tupleIdsEnd,
                                                                          DataArrayIdType *&cellRestriction, DataArrayIdType *&trueTupleRestriction) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::computeMeshRestrictionFromTupleIds : NULL input mesh !");
  MCAuto<DataArrayIdType> tmp=DataArrayIdType::New();
  tmp->alloc(std::distance(tupleIdsBg,tupleIdsEnd),1);
  std::copy(tupleIdsBg,tupleIdsEnd,tmp->getPointer());
  MCAuto<DataArrayIdType> tmp2(tmp->deepCopy());
  cellRestriction=tmp.retn();
  trueTupleRestriction=tmp2.retn();
}

void MEDCouplingFieldDiscretizationP0::reprQuickOverview(std::ostream& stream) const
{
  stream << "P0 spatial discretization.";
}

void MEDCouplingFieldDiscretizationP0::checkCompatibilityWithNature(NatureOfField nat) const
{
}

void MEDCouplingFieldDiscretizationP0::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const
{
  if(!mesh || !da)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::checkCoherencyBetween : NULL input mesh or DataArray !");
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
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::getMeasureField : mesh instance specified is NULL !");
  return mesh->getMeasureField(isAbs);
}

void MEDCouplingFieldDiscretizationP0::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::getValueOn : NULL input mesh !");
  mcIdType id=mesh->getCellContainingPoint(loc,_precision);
  if(id==-1)
    throw INTERP_KERNEL::Exception("Specified point is detected outside of mesh : unable to apply P0::getValueOn !");
  arr->getTuple(id,res);
}

void MEDCouplingFieldDiscretizationP0::getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k, double *res) const
{
  const MEDCouplingCMesh *meshC=dynamic_cast<const MEDCouplingCMesh *>(mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("P0::getValueOnPos is only accessible for structured meshes !");
  mcIdType id=meshC->getCellIdFromPos(i,j,k);
  arr->getTuple(id,res);
}

DataArrayDouble *MEDCouplingFieldDiscretizationP0::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::getValueOnMulti : NULL input mesh !");
  MCAuto<DataArrayIdType> eltsArr,eltsIndexArr;
  mesh->getCellsContainingPoints(loc,nbOfPoints,_precision,eltsArr,eltsIndexArr);
  const mcIdType *elts(eltsArr->begin()),*eltsIndex(eltsIndexArr->begin());
  int spaceDim=mesh->getSpaceDimension();
  std::size_t nbOfComponents=arr->getNumberOfComponents();
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(nbOfPoints,nbOfComponents);
  double *ptToFill=ret->getPointer();
  for(mcIdType i=0;i<nbOfPoints;i++,ptToFill+=nbOfComponents)
    if(eltsIndex[i+1]-eltsIndex[i]>=1)
      arr->getTuple(elts[eltsIndex[i]],ptToFill);
    else
      {
        std::ostringstream oss; oss << "Point #" << i << " with coordinates : (";
        std::copy(loc+i*spaceDim,loc+(i+1)*spaceDim,std::ostream_iterator<double>(oss,", "));
        oss << ") detected outside mesh : unable to apply P0::getValueOnMulti ! ";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  return ret.retn();
}

/*!
 * Nothing to do. It's not a bug.
 */
void MEDCouplingFieldDiscretizationP0::renumberValuesOnNodes(double , const mcIdType *, mcIdType newNbOfNodes, DataArrayDouble *) const
{
}

void MEDCouplingFieldDiscretizationP0::renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const mcIdType *old2New, mcIdType newSz, DataArrayDouble *arr) const
{
  RenumberEntitiesFromO2NArr(epsOnVals,old2New,newSz,arr,"Cell");
}

void MEDCouplingFieldDiscretizationP0::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const mcIdType *new2old, mcIdType newSz, DataArrayDouble *arr) const
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
DataArrayIdType *MEDCouplingFieldDiscretizationP0::computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const mcIdType *startCellIds, const mcIdType *endCellIds) const
{
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(std::distance(startCellIds,endCellIds),1);
  std::copy(startCellIds,endCellIds,ret->getPointer());
  return ret.retn();
}

/*!
 * This method returns a submesh of 'mesh' instance constituting cell ids contained in array defined as an interval [start;end).
 * @param di is an array returned that specifies entity ids (here cells ids) in mesh 'mesh' of entity in returned submesh.
 * Example : The first cell id of returned mesh has the (*di)[0] id in 'mesh'
 *
 * \sa MEDCouplingFieldDiscretizationP0::buildSubMeshDataRange
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationP0::buildSubMeshData(const MEDCouplingMesh *mesh, const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::buildSubMeshData : NULL input mesh !");
  MCAuto<MEDCouplingMesh> ret=mesh->buildPart(start,end);
  MCAuto<DataArrayIdType> diSafe=DataArrayIdType::New();
  diSafe->alloc(std::distance(start,end),1);
  std::copy(start,end,diSafe->getPointer());
  di=diSafe.retn();
  return ret.retn();
}

/*!
 * This method is strictly equivalent to MEDCouplingFieldDiscretizationP0::buildSubMeshData except that it is optimized for input defined as a range of cell ids.
 * 
 * \param [out] beginOut Valid only if \a di is NULL
 * \param [out] endOut Valid only if \a di is NULL
 * \param [out] stepOut Valid only if \a di is NULL
 * \param [out] di is an array returned that specifies entity ids (nodes, cells, Gauss points... ) in array if no output range is foundable.
 *
 * \sa MEDCouplingFieldDiscretizationP0::buildSubMeshData
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationP0::buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP0::buildSubMeshDataRange : NULL input mesh !");
  MCAuto<MEDCouplingMesh> ret=mesh->buildPartRange(beginCellIds,endCellIds,stepCellIds);
  di=0; beginOut=beginCellIds; endOut=endCellIds; stepOut=stepCellIds;
  return ret.retn();
}

MCAuto<MEDCouplingFieldDiscretization> MEDCouplingFieldDiscretizationP0::aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const
{
  return EasyAggregate<MEDCouplingFieldDiscretizationP0>(fds);
}

mcIdType MEDCouplingFieldDiscretizationOnNodes::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationNodes::getNumberOfTuples : NULL input mesh !");
  return mesh->getNumberOfNodes();
}

/*!
 * This method returns the number of tuples regarding exclusively the input code \b without \b using \b a \b mesh \b in \b input.
 * The input code coherency is also checked regarding spatial discretization of \a this.
 * If an incoherency is detected, an exception will be thrown. If the input code is coherent, the number of tuples expected is returned.
 * The number of tuples expected is equal to those to have a valid field lying on \a this and having a mesh fitting perfectly the input code (geometric type distribution).
 */
mcIdType MEDCouplingFieldDiscretizationOnNodes::getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const
{
  if(code.size()%3!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationOnNodes::getNumberOfTuplesExpectedRegardingCode : invalid input code !");
  mcIdType nbOfSplit=ToIdType(idsPerType.size());
  mcIdType nbOfTypes=ToIdType(code.size()/3);
  mcIdType ret=0;
  for(mcIdType i=0;i<nbOfTypes;i++)
    {
      mcIdType nbOfEltInChunk=code[3*i+1];
      if(nbOfEltInChunk<0)
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationOnNodes::getNumberOfTuplesExpectedRegardingCode : invalid input code ! presence of negative value in a type !");
      mcIdType pos=code[3*i+2];
      if(pos!=-1)
        {
          if(pos<0 || pos>=nbOfSplit)
            {
              std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationOnNodes::getNumberOfTuplesExpectedRegardingCode : input code points to pos " << pos << " in typeid " << i << " ! Should be in [0," << nbOfSplit << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          const DataArrayIdType *ids(idsPerType[pos]);
          if(!ids || !ids->isAllocated() || ids->getNumberOfComponents()!=1 || ids->getNumberOfTuples()!=nbOfEltInChunk || ids->getMinValueInArray()<0)
            {
              std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationOnNodes::getNumberOfTuplesExpectedRegardingCode : input pfl chunck at pos " << pos << " should have " << i << " tuples and one component and with ids all >=0 !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      ret+=nbOfEltInChunk;
    }
  return ret;
}

mcIdType MEDCouplingFieldDiscretizationOnNodes::getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationNodes::getNumberOfMeshPlaces : NULL input mesh !");
  return mesh->getNumberOfNodes();
}

/*!
 * Nothing to do here.
 */
void MEDCouplingFieldDiscretizationOnNodes::renumberArraysForCell(const MEDCouplingMesh *, const std::vector<DataArray *>& arrays,
                                                                  const mcIdType *old2NewBg, bool check)
{
}

DataArrayIdType *MEDCouplingFieldDiscretizationOnNodes::getOffsetArr(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationNodes::getOffsetArr : NULL input mesh !");
  mcIdType nbOfTuples=mesh->getNumberOfNodes();
  DataArrayIdType *ret=DataArrayIdType::New();
  ret->alloc(nbOfTuples+1,1);
  ret->iota(0);
  return ret;
}

DataArrayDouble *MEDCouplingFieldDiscretizationOnNodes::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationNodes::getLocalizationOfDiscValues : NULL input mesh !");
  return mesh->getCoordinatesAndOwner();
}

void MEDCouplingFieldDiscretizationOnNodes::computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const mcIdType *tupleIdsBg, const mcIdType *tupleIdsEnd,
                                                                               DataArrayIdType *&cellRestriction, DataArrayIdType *&trueTupleRestriction) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationOnNodes::computeMeshRestrictionFromTupleIds : NULL input mesh !");
  MCAuto<DataArrayIdType> ret1=mesh->getCellIdsFullyIncludedInNodeIds(tupleIdsBg,tupleIdsEnd);
  const MEDCouplingUMesh *meshc=dynamic_cast<const MEDCouplingUMesh *>(mesh);
  if(!meshc)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationOnNodes::computeMeshRestrictionFromTupleIds : trying to subpart field on nodes by node ids ! Your mesh has to be unstructured !");
  MCAuto<MEDCouplingUMesh> meshPart=static_cast<MEDCouplingUMesh *>(meshc->buildPartOfMySelf(ret1->begin(),ret1->end(),true));
  MCAuto<DataArrayIdType> ret2=meshPart->computeFetchedNodeIds();
  cellRestriction=ret1.retn();
  trueTupleRestriction=ret2.retn();
}

void MEDCouplingFieldDiscretizationOnNodes::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const
{
  if(!mesh || !da)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationNodes::checkCoherencyBetween : NULL input mesh or DataArray !");
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
MEDCouplingMesh *MEDCouplingFieldDiscretizationOnNodes::buildSubMeshData(const MEDCouplingMesh *mesh, const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationNodes::buildSubMeshData : NULL input mesh !");
  DataArrayIdType *diTmp=0;
  MCAuto<MEDCouplingMesh> ret=mesh->buildPartAndReduceNodes(start,end,diTmp);
  MCAuto<DataArrayIdType> diTmpSafe(diTmp);
  MCAuto<DataArrayIdType> di2=diTmpSafe->invertArrayO2N2N2O(ret->getNumberOfNodes());
  di=di2.retn();
  return ret.retn();
}

/*!
 * This method is strictly equivalent to MEDCouplingFieldDiscretizationNodes::buildSubMeshData except that it is optimized for input defined as a range of cell ids.
 * 
 * \param [out] beginOut Valid only if \a di is NULL
 * \param [out] endOut Valid only if \a di is NULL
 * \param [out] stepOut Valid only if \a di is NULL
 * \param [out] di is an array returned that specifies entity ids (nodes, cells, Gauss points... ) in array if no output range is foundable.
 *
 * \sa MEDCouplingFieldDiscretizationNodes::buildSubMeshData
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationOnNodes::buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationOnNodes::buildSubMeshDataRange : NULL input mesh !");
  DataArrayIdType *diTmp=0;
  MCAuto<MEDCouplingMesh> ret=mesh->buildPartRangeAndReduceNodes(beginCellIds,endCellIds,stepCellIds,beginOut,endOut,stepOut,diTmp);
  if(diTmp)
    {
      MCAuto<DataArrayIdType> diTmpSafe(diTmp);
      MCAuto<DataArrayIdType> di2=diTmpSafe->invertArrayO2N2N2O(ret->getNumberOfNodes());
      di=di2.retn();
    }
  return ret.retn();
}

/*!
 * This method returns a tuple ids selection from cell ids selection [start;end).
 * This method is called by MEDCouplingFieldDiscretizationOnNodes::buildSubMeshData to return parameter \b di.
 * Here for P1 only nodes fetched by submesh of mesh[startCellIds:endCellIds) is returned !
 *
 * \return a newly allocated array containing ids to select into the DataArrayDouble of the field.
 * 
 */
DataArrayIdType *MEDCouplingFieldDiscretizationOnNodes::computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const mcIdType *startCellIds, const mcIdType *endCellIds) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP1::computeTupleIdsToSelectFromCellIds : NULL input mesh !");
  const MCAuto<MEDCouplingUMesh> umesh=mesh->buildUnstructured();
  MCAuto<MEDCouplingUMesh> umesh2=static_cast<MEDCouplingUMesh *>(umesh->buildPartOfMySelf(startCellIds,endCellIds,true));
  return umesh2->computeFetchedNodeIds();
}

void MEDCouplingFieldDiscretizationOnNodes::renumberValuesOnNodes(double epsOnVals, const mcIdType *old2NewPtr, mcIdType newNbOfNodes, DataArrayDouble *arr) const
{
  RenumberEntitiesFromO2NArr(epsOnVals,old2NewPtr,newNbOfNodes,arr,"Node");
}

/*!
 * Nothing to do it's not a bug.
 */
void MEDCouplingFieldDiscretizationOnNodes::renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const mcIdType *old2New, mcIdType newSz, DataArrayDouble *arr) const
{
}

/*!
 * Nothing to do it's not a bug.
 */
void MEDCouplingFieldDiscretizationOnNodes::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const mcIdType *new2old, mcIdType newSz, DataArrayDouble *arr) const
{
}

void MEDCouplingFieldDiscretizationOnNodes::getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k, double *res) const
{
  const MEDCouplingCMesh *meshC=dynamic_cast<const MEDCouplingCMesh *>(mesh);
  if(!meshC)
    throw INTERP_KERNEL::Exception("OnNodes::getValueOnPos(i,j,k) is only accessible for structured meshes !");
  mcIdType id=meshC->getNodeIdFromPos(i,j,k);
  arr->getTuple(id,res);
}

TypeOfField MEDCouplingFieldDiscretizationP1::getEnum() const
{
  return TYPE;
}

/*!
 * This method is simply called by MEDCouplingFieldDiscretization::deepCopy. It performs the deep copy of \a this.
 *
 * \sa MEDCouplingFieldDiscretization::deepCopy.
 */
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
  if(!other)
    {
      reason="other spatial discretization is NULL, and this spatial discretization (P1) is defined.";
      return false;
    }
  const MEDCouplingFieldDiscretizationP1 *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationP1 *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason="Spatial discrtization of this is ON_NODES, which is not the case of other.";
  return ret;
}

void MEDCouplingFieldDiscretizationP1::checkCompatibilityWithNature(NatureOfField nat) const
{
  if(nat!=IntensiveMaximum)
    throw INTERP_KERNEL::Exception("Invalid nature for P1 field  : expected IntensiveMaximum !");
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationP1::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP1::getMeasureField : mesh instance specified is NULL !");
  return mesh->getMeasureFieldOnNode(isAbs);
}

void MEDCouplingFieldDiscretizationP1::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP1::getValueOn : NULL input mesh !");
  mcIdType id=mesh->getCellContainingPoint(loc,_precision);
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
void MEDCouplingFieldDiscretizationP1::getValueInCell(const MEDCouplingMesh *mesh, mcIdType cellId, const DataArrayDouble *arr, const double *loc, double *res) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP1::getValueInCell : NULL input mesh !");
  std::vector<mcIdType> conn;
  std::vector<double> coo;
  mesh->getNodeIdsOfCell(cellId,conn);
  for(std::vector<mcIdType>::const_iterator iter=conn.begin();iter!=conn.end();iter++)
    mesh->getCoordinatesOfNode(*iter,coo);
  int spaceDim=mesh->getSpaceDimension();
  std::size_t nbOfNodes=conn.size();
  std::vector<const double *> vec(nbOfNodes);
  for(std::size_t i=0;i<nbOfNodes;i++)
    vec[i]=&coo[i*spaceDim];
  INTERP_KERNEL::AutoPtr<double> tmp=new double[nbOfNodes];
  INTERP_KERNEL::NormalizedCellType ct(mesh->getTypeOfCell(cellId));
  INTERP_KERNEL::barycentric_coords(ct,vec,loc,tmp);
  std::size_t sz=arr->getNumberOfComponents();
  INTERP_KERNEL::AutoPtr<double> tmp2=new double[sz];
  std::fill(res,res+sz,0.);
  for(std::size_t i=0;i<nbOfNodes;i++)
    {
      arr->getTuple(conn[i],(double *)tmp2);
      std::transform((double *)tmp2,((double *)tmp2)+sz,(double *)tmp2,std::bind(std::multiplies<double>(),std::placeholders::_1,tmp[i]));
      std::transform(res,res+sz,(double *)tmp2,res,std::plus<double>());
    }
}

DataArrayDouble *MEDCouplingFieldDiscretizationP1::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationP1::getValueOnMulti : NULL input mesh !");
  MCAuto<DataArrayIdType> eltsArr,eltsIndexArr;
  mesh->getCellsContainingPoints(loc,nbOfPoints,_precision,eltsArr,eltsIndexArr);
  const mcIdType *elts(eltsArr->begin()),*eltsIndex(eltsIndexArr->begin());
  int spaceDim=mesh->getSpaceDimension();
  std::size_t nbOfComponents=arr->getNumberOfComponents();
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc(nbOfPoints,nbOfComponents);
  double *ptToFill=ret->getPointer();
  for(mcIdType i=0;i<nbOfPoints;i++)
    if(eltsIndex[i+1]-eltsIndex[i]>=1)
      getValueInCell(mesh,elts[eltsIndex[i]],arr,loc+i*spaceDim,ptToFill+i*nbOfComponents);
    else
      {
        std::ostringstream oss; oss << "Point #" << i << " with coordinates : (";
        std::copy(loc+i*spaceDim,loc+(i+1)*spaceDim,std::ostream_iterator<double>(oss,", "));
        oss << ") detected outside mesh : unable to apply P1::getValueOnMulti ! ";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  return ret.retn();
}

void MEDCouplingFieldDiscretizationP1::reprQuickOverview(std::ostream& stream) const
{
  stream << "P1 spatial discretization.";
}

MCAuto<MEDCouplingFieldDiscretization> MEDCouplingFieldDiscretizationP1::aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const
{
  return EasyAggregate<MEDCouplingFieldDiscretizationP1>(fds);
}

MEDCouplingFieldDiscretizationPerCell::MEDCouplingFieldDiscretizationPerCell():_discr_per_cell(0)
{
}

MEDCouplingFieldDiscretizationPerCell::~MEDCouplingFieldDiscretizationPerCell()
{
  if(_discr_per_cell)
    _discr_per_cell->decrRef();
}

/*!
 * This constructor deep copies MEDCoupling::DataArrayIdType instance from other (if any).
 */
MEDCouplingFieldDiscretizationPerCell::MEDCouplingFieldDiscretizationPerCell(const MEDCouplingFieldDiscretizationPerCell& other, const mcIdType *startCellIds, const mcIdType *endCellIds):_discr_per_cell(0)
{
  DataArrayIdType *arr=other._discr_per_cell;
  if(arr)
    {
      if(startCellIds==0 && endCellIds==0)
        _discr_per_cell=arr->deepCopy();
      else
        _discr_per_cell=arr->selectByTupleIdSafe(startCellIds,endCellIds);
    }
}

MEDCouplingFieldDiscretizationPerCell::MEDCouplingFieldDiscretizationPerCell(const MEDCouplingFieldDiscretizationPerCell& other, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds):_discr_per_cell(0)
{
  DataArrayIdType *arr=other._discr_per_cell;
  if(arr)
    {
      _discr_per_cell=arr->selectByTupleIdSafeSlice(beginCellIds,endCellIds,stepCellIds);
    }
}

MEDCouplingFieldDiscretizationPerCell::MEDCouplingFieldDiscretizationPerCell(DataArrayIdType *dpc):_discr_per_cell(dpc)
{
  if(_discr_per_cell)
    _discr_per_cell->incrRef();
}

void MEDCouplingFieldDiscretizationPerCell::updateTime() const
{
  if(_discr_per_cell)
    updateTimeWith(*_discr_per_cell);
}

std::size_t MEDCouplingFieldDiscretizationPerCell::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(MEDCouplingFieldDiscretization::getHeapMemorySizeWithoutChildren());
  return ret;
}

std::vector<const BigMemoryObject *> MEDCouplingFieldDiscretizationPerCell::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret(MEDCouplingFieldDiscretization::getDirectChildrenWithNull());
  ret.push_back(_discr_per_cell);
  return ret;
}

void MEDCouplingFieldDiscretizationPerCell::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell has no discretization per cell !");
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell::checkCoherencyBetween : NULL input mesh or DataArray !");
  mcIdType nbOfTuples(_discr_per_cell->getNumberOfTuples());
  if(nbOfTuples!=mesh->getNumberOfCells())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell has a discretization per cell but it's not matching the underlying mesh !");
}

bool MEDCouplingFieldDiscretizationPerCell::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  if(!other)
    {
      reason="other spatial discretization is NULL, and this spatial discretization (PerCell) is defined.";
      return false;
    }
  const MEDCouplingFieldDiscretizationPerCell *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationPerCell *>(other);
  if(!otherC)
    {
      reason="Spatial discretization of this is ON_GAUSS, which is not the case of other.";
      return false;
    }
  if(_discr_per_cell==0)
    return otherC->_discr_per_cell==0;
  if(otherC->_discr_per_cell==0)
    return false;
  bool ret=_discr_per_cell->isEqualIfNotWhy(*otherC->_discr_per_cell,reason);
  if(!ret)
    reason.insert(0,"Field discretization per cell DataArrayIdType given the discid per cell :");
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
 * virtually by this method.
 */
void MEDCouplingFieldDiscretizationPerCell::renumberCells(const mcIdType *old2NewBg, bool check)
{
  mcIdType nbCells=_discr_per_cell->getNumberOfTuples();
  const mcIdType *array=old2NewBg;
  if(check)
    array=DataArrayIdType::CheckAndPreparePermutation(old2NewBg,old2NewBg+nbCells);
  //
  DataArrayIdType *dpc=_discr_per_cell->renumber(array);
  _discr_per_cell->decrRef();
  _discr_per_cell=dpc;
  //
  if(check)
    free(const_cast<mcIdType *>(array));
}

void MEDCouplingFieldDiscretizationPerCell::buildDiscrPerCellIfNecessary(const MEDCouplingMesh *mesh)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell::buildDiscrPerCellIfNecessary : NULL input mesh !");
  if(!_discr_per_cell)
    {
      _discr_per_cell=DataArrayIdType::New();
      mcIdType nbTuples=mesh->getNumberOfCells();
      _discr_per_cell->alloc(nbTuples,1);
      mcIdType *ptr=_discr_per_cell->getPointer();
      std::fill(ptr,ptr+nbTuples,DFT_INVALID_LOCID_VALUE);
    }
}

void MEDCouplingFieldDiscretizationPerCell::checkNoOrphanCells() const
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell::checkNoOrphanCells : no discretization defined !");
  MCAuto<DataArrayIdType> test( _discr_per_cell->findIdsEqual(DFT_INVALID_LOCID_VALUE));
  if(test->getNumberOfTuples()!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell::checkNoOrphanCells : presence of orphan cells !");
}

/*!
 * This method is useful when 'this' describes a field discretization with several gauss discretization on a \b same cell type.
 * For example same NORM_TRI3 cells having 6 gauss points and others with 12 gauss points.
 * This method returns 2 arrays with same size : the return value and 'locIds' output parameter.
 * For a given i into [0,locIds.size) ret[i] represents the set of cell ids of i_th set an locIds[i] represents the set of discretisation of the set.
 * The return vector contains a set of newly created instance to deal with.
 * The returned vector represents a \b partition of cells ids with a gauss discretization set.
 * 
 * If no descretization is set in 'this' and exception will be thrown.
 */
std::vector<DataArrayIdType *> MEDCouplingFieldDiscretizationPerCell::splitIntoSingleGaussDicrPerCellType(std::vector<mcIdType>& locIds) const
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationPerCell::splitIntoSingleGaussDicrPerCellType : no descretization set !");
  return _discr_per_cell->partitionByDifferentValues(locIds);
}

const DataArrayIdType *MEDCouplingFieldDiscretizationPerCell::getArrayOfDiscIds() const
{
  return _discr_per_cell;
}

void MEDCouplingFieldDiscretizationPerCell::setArrayOfDiscIds(const DataArrayIdType *adids)
{
  if(adids!=_discr_per_cell)
    {
      if(_discr_per_cell)
        _discr_per_cell->decrRef();
      _discr_per_cell=const_cast<DataArrayIdType *>(adids);
      if(_discr_per_cell)
        _discr_per_cell->incrRef();
      declareAsNew();
    }
}

MEDCouplingFieldDiscretizationGauss::MEDCouplingFieldDiscretizationGauss()
{
}

MEDCouplingFieldDiscretizationGauss::MEDCouplingFieldDiscretizationGauss(const MEDCouplingFieldDiscretizationGauss& other, const mcIdType *startCellIds, const mcIdType *endCellIds):MEDCouplingFieldDiscretizationPerCell(other,startCellIds,endCellIds),_loc(other._loc)
{
}

MEDCouplingFieldDiscretizationGauss::MEDCouplingFieldDiscretizationGauss(const MEDCouplingFieldDiscretizationGauss& other, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds):MEDCouplingFieldDiscretizationPerCell(other,beginCellIds,endCellIds,stepCellIds),_loc(other._loc)
{
}

TypeOfField MEDCouplingFieldDiscretizationGauss::getEnum() const
{
  return TYPE;
}

bool MEDCouplingFieldDiscretizationGauss::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  if(!other)
    {
      reason="other spatial discretization is NULL, and this spatial discretization (Gauss) is defined.";
      return false;
    }
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

/*!
 * This method is simply called by MEDCouplingFieldDiscretization::deepCopy. It performs the deep copy of \a this.
 *
 * \sa MEDCouplingFieldDiscretization::deepCopy.
 */
MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationGauss::clone() const
{
  return new MEDCouplingFieldDiscretizationGauss(*this);
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationGauss::clonePart(const mcIdType *startCellIds, const mcIdType *endCellIds) const
{
  return new MEDCouplingFieldDiscretizationGauss(*this,startCellIds,endCellIds);
}

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationGauss::clonePartRange(mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds) const
{
  return new MEDCouplingFieldDiscretizationGauss(*this,beginCellIds,endCellIds,stepCellIds);
}

std::string MEDCouplingFieldDiscretizationGauss::getStringRepr() const
{
  std::ostringstream oss; oss << REPR << "." << std::endl;
  if(_discr_per_cell)
    {
      if(_discr_per_cell->isAllocated())
        {
          oss << "Discretization per cell : ";
          std::copy(_discr_per_cell->begin(),_discr_per_cell->end(),std::ostream_iterator<mcIdType>(oss,", "));
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

std::size_t MEDCouplingFieldDiscretizationGauss::getHeapMemorySizeWithoutChildren() const
{
  std::size_t ret(MEDCouplingFieldDiscretizationPerCell::getHeapMemorySizeWithoutChildren());
  ret+=_loc.capacity()*sizeof(MEDCouplingGaussLocalization);
  for(std::vector<MEDCouplingGaussLocalization>::const_iterator it=_loc.begin();it!=_loc.end();it++)
    ret+=(*it).getMemorySize();
  return ret;
}

const char *MEDCouplingFieldDiscretizationGauss::getRepr() const
{
  return REPR;
}

/*!
 * This method returns the number of tuples regarding exclusively the input code \b without \b using \b a \b mesh \b in \b input.
 * The input code coherency is also checked regarding spatial discretization of \a this.
 * If an incoherency is detected, an exception will be thrown. If the input code is coherent, the number of tuples expected is returned.
 * The number of tuples expected is equal to those to have a valid field lying on \a this and having a mesh fitting perfectly the input code (geometric type distribution).
 */
mcIdType MEDCouplingFieldDiscretizationGauss::getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const
{
  if(!_discr_per_cell || !_discr_per_cell->isAllocated() || _discr_per_cell->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getNumberOfTuplesExpectedRegardingCode");
  if(code.size()%3!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getNumberOfTuplesExpectedRegardingCode : invalid input code !");
  mcIdType nbOfSplit=ToIdType(idsPerType.size());
  mcIdType nbOfTypes=ToIdType(code.size()/3);
  mcIdType ret(0);
  for(mcIdType i=0;i<nbOfTypes;i++)
    {
      mcIdType nbOfEltInChunk=code[3*i+1];
      if(nbOfEltInChunk<0)
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getNumberOfTuplesExpectedRegardingCode : invalid input code ! presence of negative value in a type !");
      mcIdType pos=code[3*i+2];
      if(pos!=-1)
        {
          if(pos<0 || pos>=nbOfSplit)
            {
              std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::getNumberOfTuplesExpectedRegardingCode : input code points to pos " << pos << " in typeid " << i << " ! Should be in [0," << nbOfSplit << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          const DataArrayIdType *ids(idsPerType[pos]);
          if(!ids || !ids->isAllocated() || ids->getNumberOfComponents()!=1 || ids->getNumberOfTuples()!=nbOfEltInChunk || ids->getMinValueInArray()<0)
            {
              std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::getNumberOfTuplesExpectedRegardingCode : input pfl chunck at pos " << pos << " should have " << i << " tuples and one component and with ids all >=0 !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      ret+=nbOfEltInChunk;
    }
  if(ret!=_discr_per_cell->getNumberOfTuples())
    {
      std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::getNumberOfTuplesExpectedRegardingCode : input code points to " << ret << " cells whereas discretization percell array lgth is " <<  _discr_per_cell->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return getNumberOfTuples(0);//0 is not an error ! It is to be sure that input mesh is not used
}

mcIdType MEDCouplingFieldDiscretizationGauss::getNumberOfTuples(const MEDCouplingMesh *) const
{
  mcIdType ret=0;
  if (_discr_per_cell == 0)
    throw INTERP_KERNEL::Exception("Discretization is not initialized!");
  const mcIdType *dcPtr=_discr_per_cell->getConstPointer();
  mcIdType nbOfTuples=_discr_per_cell->getNumberOfTuples();
  mcIdType maxSz=ToIdType(_loc.size());
  for(const mcIdType *w=dcPtr;w!=dcPtr+nbOfTuples;w++)
    {
      if(*w>=0 && *w<maxSz)
        ret+=_loc[*w].getNumberOfGaussPt();
      else
        {
          std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::getNumberOfTuples : At cell #" << std::distance(dcPtr,w) << " localization id is " << *w << " should be in [0," << maxSz << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret;
}

mcIdType MEDCouplingFieldDiscretizationGauss::getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getNumberOfMeshPlaces : NULL input mesh !");
  return mesh->getNumberOfCells();
}

/*!
 * This method is redevelopped for performance reasons, but it is equivalent to a call to MEDCouplingFieldDiscretizationGauss::buildNbOfGaussPointPerCellField
 * and a call to DataArrayDouble::computeOffsetsFull on the returned array.
 */
DataArrayIdType *MEDCouplingFieldDiscretizationGauss::getOffsetArr(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getOffsetArr : NULL input mesh !");
  mcIdType nbOfTuples=mesh->getNumberOfCells();
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  ret->alloc(nbOfTuples+1,1);
  mcIdType *retPtr(ret->getPointer());
  const mcIdType *start(_discr_per_cell->begin());
  if(_discr_per_cell->getNumberOfTuples()!=nbOfTuples)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getOffsetArr : mismatch between the mesh and the discretization ids array length !");
  mcIdType maxPossible=ToIdType(_loc.size());
  retPtr[0]=0;
  for(mcIdType i=0;i<nbOfTuples;i++,start++)
    {
      if(*start>=0 && *start<maxPossible)
        retPtr[i+1]=retPtr[i]+_loc[*start].getNumberOfGaussPt();
      else
        {
          std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::getOffsetArr : At position #" << i << " the locid = " << *start << " whereas it should be in [0," << maxPossible << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  return ret.retn();
}

void MEDCouplingFieldDiscretizationGauss::renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArray *>& arrays,
                                                                const mcIdType *old2NewBg, bool check)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::renumberArraysForCell : NULL input mesh !");
  const mcIdType *array=old2NewBg;
  if(check)
    array=DataArrayIdType::CheckAndPreparePermutation(old2NewBg,old2NewBg+mesh->getNumberOfCells());
  mcIdType nbOfCells=_discr_per_cell->getNumberOfTuples();
  mcIdType nbOfTuples=getNumberOfTuples(0);
  const mcIdType *dcPtr=_discr_per_cell->getConstPointer();
  mcIdType *array2=new mcIdType[nbOfTuples];//stores the final conversion array old2New to give to arrays in renumberInPlace.
  mcIdType *array3=new mcIdType[nbOfCells];//store for each cell in present dcp array (already renumbered) the offset needed by each cell in new numbering.
  array3[0]=0;
  for(mcIdType i=1;i<nbOfCells;i++)
    array3[i]=array3[i-1]+_loc[dcPtr[i-1]].getNumberOfGaussPt();
  mcIdType j=0;
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      mcIdType nbOfGaussPt=_loc[dcPtr[array[i]]].getNumberOfGaussPt();
      for(mcIdType k=0;k<nbOfGaussPt;k++,j++)
        array2[j]=array3[array[i]]+k;
    }
  delete [] array3;
  for(std::vector<DataArray *>::const_iterator it=arrays.begin();it!=arrays.end();it++)
    if(*it)
      (*it)->renumberInPlace(array2);
  delete [] array2;
  if(check)
    free(const_cast<mcIdType*>(array));
}

DataArrayDouble *MEDCouplingFieldDiscretizationGauss::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getLocalizationOfDiscValues : NULL input mesh !");
  checkNoOrphanCells();
  MCAuto<MEDCouplingUMesh> umesh=mesh->buildUnstructured();//in general do nothing
  mcIdType nbOfTuples=getNumberOfTuples(mesh);
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  int spaceDim=mesh->getSpaceDimension();
  ret->alloc(nbOfTuples,spaceDim);
  std::vector< mcIdType > locIds;
  std::vector<DataArrayIdType *> parts=splitIntoSingleGaussDicrPerCellType(locIds);
  std::vector< MCAuto<DataArrayIdType> > parts2(parts.size());
  std::copy(parts.begin(),parts.end(),parts2.begin());
  MCAuto<DataArrayIdType> offsets=buildNbOfGaussPointPerCellField();
  offsets->computeOffsets();
  const mcIdType *ptrOffsets=offsets->getConstPointer();
  const double *coords=umesh->getCoords()->getConstPointer();
  const mcIdType *connI=umesh->getNodalConnectivityIndex()->getConstPointer();
  const mcIdType *conn=umesh->getNodalConnectivity()->getConstPointer();
  double *valsToFill=ret->getPointer();
  for(std::size_t i=0;i<parts2.size();i++)
    {
      INTERP_KERNEL::GaussCoords calculator;
      //
      const MEDCouplingGaussLocalization& cli(_loc[locIds[i]]);//curLocInfo
      INTERP_KERNEL::NormalizedCellType typ(cli.getType());
      const std::vector<double>& wg(cli.getWeights());
      calculator.addGaussInfo(typ,INTERP_KERNEL::CellModel::GetCellModel(typ).getDimension(),
                              &cli.getGaussCoords()[0],ToIdType(wg.size()),&cli.getRefCoords()[0],
          INTERP_KERNEL::CellModel::GetCellModel(typ).getNumberOfNodes());
      //
      for(const mcIdType *w=parts2[i]->begin();w!=parts2[i]->end();w++)
        calculator.calculateCoords(cli.getType(),coords,spaceDim,conn+connI[*w]+1,valsToFill+spaceDim*(ptrOffsets[*w]));
    }
  ret->copyStringInfoFrom(*umesh->getCoords());
  return ret.retn();
}

void MEDCouplingFieldDiscretizationGauss::computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const mcIdType *tupleIdsBg, const mcIdType *tupleIdsEnd,
                                                                             DataArrayIdType *&cellRestriction, DataArrayIdType *&trueTupleRestriction) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::computeMeshRestrictionFromTupleIds : NULL input mesh !");
  MCAuto<DataArrayIdType> tmp=DataArrayIdType::New(); tmp->alloc(std::distance(tupleIdsBg,tupleIdsEnd),1);
  std::copy(tupleIdsBg,tupleIdsEnd,tmp->getPointer());
  tmp->sort(true);
  tmp=tmp->buildUnique();
  MCAuto<DataArrayIdType> nbOfNodesPerCell=buildNbOfGaussPointPerCellField();
  nbOfNodesPerCell->computeOffsetsFull();
  nbOfNodesPerCell->findIdsRangesInListOfIds(tmp,cellRestriction,trueTupleRestriction);
}

/*!
 * Empty : not a bug
 */
void MEDCouplingFieldDiscretizationGauss::checkCompatibilityWithNature(NatureOfField nat) const
{
}

void MEDCouplingFieldDiscretizationGauss::getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const
{
  mcIdType val=-1;
  if(_discr_per_cell)
    val=_discr_per_cell->getNumberOfTuples();
  tinyInfo.push_back(val);
  tinyInfo.push_back(ToIdType(_loc.size()));
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

void MEDCouplingFieldDiscretizationGauss::getSerializationIntArray(DataArrayIdType *& arr) const
{
  arr=0;
  if(_discr_per_cell)
    arr=_discr_per_cell;
}

void MEDCouplingFieldDiscretizationGauss::resizeForUnserialization(const std::vector<mcIdType>& tinyInfo, DataArrayIdType *& arr)
{
  mcIdType val=tinyInfo[0];
  if(val>=0)
    {
      _discr_per_cell=DataArrayIdType::New();
      _discr_per_cell->alloc(val,1);
    }
  else
    _discr_per_cell=0;
  arr=_discr_per_cell;
  commonUnserialization(tinyInfo);
}

void MEDCouplingFieldDiscretizationGauss::checkForUnserialization(const std::vector<mcIdType>& tinyInfo, const DataArrayIdType *arr)
{
  static const char MSG[]="MEDCouplingFieldDiscretizationGauss::checkForUnserialization : expect to have one not null DataArrayIdType !";
  mcIdType val=tinyInfo[0];
  if(val>=0)
    {
      if(!arr)
        throw INTERP_KERNEL::Exception(MSG);
      arr->checkNbOfTuplesAndComp(val,1,MSG);
      _discr_per_cell=const_cast<DataArrayIdType *>(arr);
      _discr_per_cell->incrRef();
    }
  else
    _discr_per_cell=0;
  commonUnserialization(tinyInfo);
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

double MEDCouplingFieldDiscretizationGauss::getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, mcIdType cellId, mcIdType nodeIdInCell, int compoId) const
{
  mcIdType offset=getOffsetOfCell(cellId);
  return da->getIJ(offset+nodeIdInCell,compoId);
}

void MEDCouplingFieldDiscretizationGauss::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const
{
  if(!mesh || !da)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::checkCoherencyBetween : NULL input mesh or DataArray !");
  MEDCouplingFieldDiscretizationPerCell::checkCoherencyBetween(mesh,da);
  for(std::vector<MEDCouplingGaussLocalization>::const_iterator iter=_loc.begin();iter!=_loc.end();iter++)
    (*iter).checkConsistencyLight();
  mcIdType nbOfDesc=ToIdType(_loc.size());
  mcIdType nbOfCells=mesh->getNumberOfCells();
  const mcIdType *dc=_discr_per_cell->getConstPointer();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      if(dc[i]>=nbOfDesc)
        {
          std::ostringstream oss; oss << "Cell # " << i << " of mesh \"" << mesh->getName() << "\" has an undefined gauss location ! Should never happened !";
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
  mcIdType nbOfTuples(getNumberOfTuples(mesh));
  if(nbOfTuples!=da->getNumberOfTuples())
    {
      std::ostringstream oss; oss << "Invalid number of tuples in the array : expecting " << nbOfTuples << " having " << da->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationGauss::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getMeasureField : mesh instance specified is NULL !");
  MCAuto<MEDCouplingUMesh> umesh(mesh->buildUnstructured());
  const double *coordsOfMesh( umesh->getCoords()->begin() );
  auto spaceDim(mesh->getSpaceDimension());
  auto meshDim(mesh->getMeshDimension());
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_GAUSS_PT);
  ret->setMesh(mesh);
  ret->setDiscretization(const_cast<MEDCouplingFieldDiscretizationGauss *>(this));
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getMeasureField : no discr per cell array not defined ! spatial localization is incorrect !");
  _discr_per_cell->checkAllocated();
  if(_discr_per_cell->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getMeasureField : no discr per cell array defined but with nb of components different from 1 !");
  if(_discr_per_cell->getNumberOfTuples()!=mesh->getNumberOfCells())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::getMeasureField : no discr per cell array defined but mismatch between nb of cells of mesh and size of spatial disr array !");
  MCAuto<DataArrayIdType> offset=getOffsetArr(mesh);
  MCAuto<DataArrayDouble> arr=DataArrayDouble::New(); arr->alloc(getNumberOfTuples(mesh),1);
  ret->setArray(arr);
  double *arrPtr=arr->getPointer();
  const mcIdType *offsetPtr=offset->getConstPointer();
  mcIdType maxGaussLoc=ToIdType(_loc.size());
  std::vector<mcIdType> locIds;
  std::vector<DataArrayIdType *> ids=splitIntoSingleGaussDicrPerCellType(locIds);
  std::vector< MCAuto<DataArrayIdType> > ids2(ids.size()); std::copy(ids.begin(),ids.end(),ids2.begin());
  for(std::size_t i=0;i<locIds.size();i++)
    {
      const DataArrayIdType *curIds=ids[i];
      mcIdType locId=locIds[i];
      if(locId>=0 && locId<maxGaussLoc)
        {
          const MEDCouplingGaussLocalization& loc=_loc[locId];
          mcIdType nbOfGaussPt=loc.getNumberOfGaussPt();
          INTERP_KERNEL::AutoPtr<double> weights=new double[nbOfGaussPt];
          for(const mcIdType *cellId=curIds->begin();cellId!=curIds->end();cellId++)
          {
            std::vector<mcIdType> conn;
            umesh->getNodeIdsOfCell(*cellId,conn);
            std::vector<double> ptsInCell; ptsInCell.reserve(conn.size()*loc.getDimension());
            std::for_each( conn.cbegin(), conn.cend(), [spaceDim,coordsOfMesh,&ptsInCell](mcIdType c) { ptsInCell.insert(ptsInCell.end(),coordsOfMesh+c*spaceDim,coordsOfMesh+(c+1)*spaceDim); } );
            std::size_t nbPtsInCell(ptsInCell.size()/spaceDim);
            INTERP_KERNEL::DenseMatrix jacobian(spaceDim,meshDim);
            MCAuto<DataArrayDouble> shapeFunc = loc.getDerivativeOfShapeFunctionValues();
            for(mcIdType iGPt = 0 ; iGPt < nbOfGaussPt ; ++iGPt)
            {
              for(auto i = 0 ; i < spaceDim ; ++i)
                for(auto j = 0 ; j < meshDim ; ++j)
                {
                  double res = 0.0;
                  for( std::size_t k = 0 ; k < nbPtsInCell ; ++k )
                    res += ptsInCell[spaceDim*k+i] * shapeFunc->getIJ(iGPt,meshDim*k+j);
                  jacobian[ i ][ j ] = res;
                }
              arrPtr[offsetPtr[*cellId]+iGPt]=std::abs( jacobian.toJacobian() )*loc.getWeight(FromIdType<int>(iGPt));
            }
          }
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::getMeasureField : Presence of localization id " << locId << " in cell #" << curIds->getIJ(0,0) << " ! Must be in [0," << maxGaussLoc << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

void MEDCouplingFieldDiscretizationGauss::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGauss::getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k, double *res) const
{
  throw INTERP_KERNEL::Exception("getValueOnPos(i,j,k) : Not applicable for Gauss points !");
}

DataArrayDouble *MEDCouplingFieldDiscretizationGauss::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const
{
  throw INTERP_KERNEL::Exception("getValueOnMulti : Not implemented yet for gauss points !");
}

MEDCouplingMesh *MEDCouplingFieldDiscretizationGauss::buildSubMeshData(const MEDCouplingMesh *mesh, const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::buildSubMeshData : NULL input mesh !");
  MCAuto<DataArrayIdType> diSafe=computeTupleIdsToSelectFromCellIds(mesh,start,end);
  MCAuto<MEDCouplingMesh> ret=mesh->buildPart(start,end);
  di=diSafe.retn();
  return ret.retn();
}

/*!
 * This method is strictly equivalent to MEDCouplingFieldDiscretizationGauss::buildSubMeshData except that it is optimized for input defined as a range of cell ids.
 * 
 * \param [out] beginOut Valid only if \a di is NULL
 * \param [out] endOut Valid only if \a di is NULL
 * \param [out] stepOut Valid only if \a di is NULL
 * \param [out] di is an array returned that specifies entity ids (nodes, cells, Gauss points... ) in array if no output range is foundable.
 *
 * \sa MEDCouplingFieldDiscretizationGauss::buildSubMeshData
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationGauss::buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const
{
  if(stepCellIds!=1)//even for stepCellIds==-1 the output will not be a range
    return MEDCouplingFieldDiscretization::buildSubMeshDataRange(mesh,beginCellIds,endCellIds,stepCellIds,beginOut,endOut,stepOut,di);
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::buildSubMeshDataRange : NULL input mesh !");
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::buildSubMeshDataRange : no discretization array set !");
  di=0; beginOut=0; endOut=0; stepOut=stepCellIds;
  const char msg[]="MEDCouplingFieldDiscretizationGauss::buildSubMeshDataRange : cell #";
  mcIdType nbOfTuples=_discr_per_cell->getNumberOfTuples();
  const mcIdType *w=_discr_per_cell->begin();
  mcIdType nbMaxOfLocId=ToIdType(_loc.size());
  for(mcIdType i=0;i<nbOfTuples;i++,w++)
    {
      if(*w!=DFT_INVALID_LOCID_VALUE)
        {
          if(*w>=0 && *w<nbMaxOfLocId)
            {
              mcIdType delta=_loc[*w].getNumberOfGaussPt();
              if(i<beginCellIds)
                beginOut+=delta;
              endOut+=delta;
              if(i>=endCellIds)
                break;
            }
          else
            { std::ostringstream oss; oss << msg << i << " has invalid id (" << *w << ") ! Should be in [0," << nbMaxOfLocId << ") !"; throw INTERP_KERNEL::Exception(oss.str().c_str()); }
        }
      else
        { std::ostringstream oss; oss << msg << i << " is detected as orphan !"; throw INTERP_KERNEL::Exception(oss.str().c_str()); }
    }
  MCAuto<MEDCouplingMesh> ret=mesh->buildPartRange(beginCellIds,endCellIds,stepCellIds);
  return ret.retn();
}

/*!
 * This method returns a tuple ids selection from cell ids selection [start;end).
 * This method is called by MEDCouplingFieldDiscretizationGauss::buildSubMeshData to return parameter \b di.
 *
 * \return a newly allocated array containing ids to select into the DataArrayDouble of the field.
 * 
 */
DataArrayIdType *MEDCouplingFieldDiscretizationGauss::computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const mcIdType *startCellIds, const mcIdType *endCellIds) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::computeTupleIdsToSelectFromCellIds : null mesh !");
  MCAuto<DataArrayIdType> nbOfNodesPerCell=buildNbOfGaussPointPerCellField();//check of _discr_per_cell not NULL pointer
  mcIdType nbOfCells(mesh->getNumberOfCells());
  if(_discr_per_cell->getNumberOfTuples()!=nbOfCells)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::computeTupleIdsToSelectFromCellIds : mismatch of nb of tuples of cell ids array and number of cells !");
  nbOfNodesPerCell->computeOffsetsFull();
  MCAuto<DataArrayIdType> sel=DataArrayIdType::New(); sel->useArray(startCellIds,false,DeallocType::CPP_DEALLOC,ToIdType(std::distance(startCellIds,endCellIds)),1);
  return sel->buildExplicitArrByRanges(nbOfNodesPerCell);
}

/*!
 * No implementation needed !
 */
void MEDCouplingFieldDiscretizationGauss::renumberValuesOnNodes(double , const mcIdType *, mcIdType newNbOfNodes, DataArrayDouble *) const
{
}

void MEDCouplingFieldDiscretizationGauss::renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const mcIdType *old2New, mcIdType newSz, DataArrayDouble *arr) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGauss::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const mcIdType *new2old, mcIdType newSz, DataArrayDouble *arr) const
{
  throw INTERP_KERNEL::Exception("Number of cells has changed and becomes higher with some cells that have been split ! Unable to conserve the Gauss field !");
}

MCAuto<MEDCouplingFieldDiscretization> MEDCouplingFieldDiscretizationGauss::aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const
{
  if(fds.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::aggregate : input array is empty");
  std::vector<MEDCouplingGaussLocalization> loc;//store the localizations for the output GaussDiscretization object
  std::vector< MCAuto<DataArrayIdType> > discPerCells(fds.size());
  std::size_t i(0);
  for(auto it=fds.begin();it!=fds.end();++it,++i)
    {
      const MEDCouplingFieldDiscretizationGauss *itc(dynamic_cast<const MEDCouplingFieldDiscretizationGauss *>(*it));
      if(!itc)
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::aggregate : same field discretization expected for all input discretizations !");
      //
      std::vector<MEDCouplingGaussLocalization> loc2(itc->_loc);
      std::vector<mcIdType> newLocId(loc2.size());
      for(std::size_t j=0;j<loc2.size();++j)
        {
          std::size_t k(0);
          for(;k<loc.size();++k)
            {
              if(loc2[j].isEqual(loc[k],1e-10))
                {
                  newLocId[j]=ToIdType(k);
                  break;
                }
            }
          if(k==loc.size())// current loc2[j]
            {
              newLocId[j]=ToIdType(loc.size());
              loc.push_back(loc2[j]);
            }
        }
      const DataArrayIdType *dpc(itc->_discr_per_cell);
      if(!dpc)
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::aggregate : Presence of nullptr array of disc per cell !");
      MCAuto<DataArrayIdType> dpc2(dpc->deepCopy());
      dpc2->transformWithIndArr(newLocId.data(),newLocId.data()+newLocId.size());
      discPerCells[i]=dpc2;
    }
  MCAuto<DataArrayIdType> dpc3(DataArrayIdType::Aggregate(ToConstVect(discPerCells)));
  MCAuto<MEDCouplingFieldDiscretizationGauss> ret(new MEDCouplingFieldDiscretizationGauss(dpc3,loc));
  return DynamicCast<MEDCouplingFieldDiscretizationGauss,MEDCouplingFieldDiscretization>(ret);
}

void MEDCouplingFieldDiscretizationGauss::setGaussLocalizationOnType(const MEDCouplingMesh *mesh, INTERP_KERNEL::NormalizedCellType type, const std::vector<double>& refCoo,
                                                                     const std::vector<double>& gsCoo, const std::vector<double>& wg)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::setGaussLocalizationOnType : NULL input mesh !");
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
  if(ToIdType(cm.getDimension())!=mesh->getMeshDimension())
    {
      std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::setGaussLocalizationOnType : mismatch of dimensions ! MeshDim==" << mesh->getMeshDimension();
      oss << " whereas Type '" << cm.getRepr() << "' has dimension " << cm.getDimension() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  buildDiscrPerCellIfNecessary(mesh);
  mcIdType id=ToIdType(_loc.size());
  MEDCouplingGaussLocalization elt(type,refCoo,gsCoo,wg);
  _loc.push_back(elt);
  mcIdType *ptr=_discr_per_cell->getPointer();
  mcIdType nbCells=mesh->getNumberOfCells();
  for(mcIdType i=0;i<nbCells;i++)
    if(mesh->getTypeOfCell(i)==type)
      ptr[i]=id;
  zipGaussLocalizations();
}

void MEDCouplingFieldDiscretizationGauss::setGaussLocalizationOnCells(const MEDCouplingMesh *mesh, const mcIdType *begin, const mcIdType *end, const std::vector<double>& refCoo,
                                                                      const std::vector<double>& gsCoo, const std::vector<double>& wg)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::setGaussLocalizationOnCells : NULL input mesh !");
  buildDiscrPerCellIfNecessary(mesh);
  if(std::distance(begin,end)<1)
    throw INTERP_KERNEL::Exception("Size of [begin,end) must be equal or greater than 1 !");
  INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(*begin);
  MEDCouplingGaussLocalization elt(type,refCoo,gsCoo,wg);
  mcIdType id=ToIdType(_loc.size());
  mcIdType *ptr=_discr_per_cell->getPointer();
  for(const mcIdType *w=begin+1;w!=end;w++)
    {
      if(mesh->getTypeOfCell(*w)!=type)
        {
          std::ostringstream oss; oss << "The cell with id " << *w << " has been detected to be incompatible in the [begin,end) array specified !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  //
  for(const mcIdType *w2=begin;w2!=end;w2++)
    ptr[*w2]=id;
  //
  _loc.push_back(elt);
  zipGaussLocalizations();
}

void MEDCouplingFieldDiscretizationGauss::clearGaussLocalizations()
{
  if(_discr_per_cell)
    {
      _discr_per_cell->decrRef();
      _discr_per_cell=0;
    }
  _loc.clear();
}

void MEDCouplingFieldDiscretizationGauss::setGaussLocalization(mcIdType locId, const MEDCouplingGaussLocalization& loc)
{
  if(locId<0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::setGaussLocalization : localization id has to be >=0 !");
  mcIdType sz=ToIdType(_loc.size());
  MEDCouplingGaussLocalization gLoc(INTERP_KERNEL::NORM_ERROR);
  if(locId>=sz)
    _loc.resize(locId+1,gLoc);
  _loc[locId]=loc;
}

void MEDCouplingFieldDiscretizationGauss::resizeLocalizationVector(mcIdType newSz)
{
  if(newSz<0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::resizeLocalizationVector : new size has to be >=0 !");
  MEDCouplingGaussLocalization gLoc(INTERP_KERNEL::NORM_ERROR);
  _loc.resize(newSz,gLoc);
}

MEDCouplingGaussLocalization& MEDCouplingFieldDiscretizationGauss::getGaussLocalization(mcIdType locId)
{
  checkLocalizationId(locId);
  return _loc[locId];
}

mcIdType MEDCouplingFieldDiscretizationGauss::getNbOfGaussLocalization() const
{
  return ToIdType(_loc.size());
}

mcIdType MEDCouplingFieldDiscretizationGauss::getGaussLocalizationIdOfOneCell(mcIdType cellId) const
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("No Gauss localization still set !");
  mcIdType locId=_discr_per_cell->begin()[cellId];
  if(locId<0)
    throw INTERP_KERNEL::Exception("No Gauss localization set for the specified cell !");
  return locId;
}

mcIdType MEDCouplingFieldDiscretizationGauss::getGaussLocalizationIdOfOneType(INTERP_KERNEL::NormalizedCellType type) const
{
  std::set<mcIdType> ret=getGaussLocalizationIdsOfOneType(type);
  if(ret.empty())
    throw INTERP_KERNEL::Exception("No gauss discretization found for the specified type !");
  if(ret.size()>1)
    throw INTERP_KERNEL::Exception("Several gauss discretizations have been found for the specified type !");
  return *ret.begin();
}

std::set<mcIdType> MEDCouplingFieldDiscretizationGauss::getGaussLocalizationIdsOfOneType(INTERP_KERNEL::NormalizedCellType type) const
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("No Gauss localization still set !");
  std::set<mcIdType> ret;
  mcIdType id=0;
  for(std::vector<MEDCouplingGaussLocalization>::const_iterator iter=_loc.begin();iter!=_loc.end();iter++,id++)
    if((*iter).getType()==type)
      ret.insert(id);
  return ret;
}

void MEDCouplingFieldDiscretizationGauss::getCellIdsHavingGaussLocalization(mcIdType locId, std::vector<mcIdType>& cellIds) const
{
  if(locId<0 || locId>=ToIdType(_loc.size()))
    throw INTERP_KERNEL::Exception("Invalid locId given : must be in range [0:getNbOfGaussLocalization()) !");
  mcIdType nbOfTuples=_discr_per_cell->getNumberOfTuples();
  const mcIdType *ptr=_discr_per_cell->getConstPointer();
  for(mcIdType i=0;i<nbOfTuples;i++)
    if(ptr[i]==locId)
      cellIds.push_back(i);
}

const MEDCouplingGaussLocalization& MEDCouplingFieldDiscretizationGauss::getGaussLocalization(mcIdType locId) const
{
  checkLocalizationId(locId);
  return _loc[locId];
}

void MEDCouplingFieldDiscretizationGauss::checkLocalizationId(mcIdType locId) const
{
  if(locId<0 || locId>=ToIdType(_loc.size()))
    throw INTERP_KERNEL::Exception("Invalid locId given : must be in range [0:getNbOfGaussLocalization()) !");
}

mcIdType MEDCouplingFieldDiscretizationGauss::getOffsetOfCell(mcIdType cellId) const
{
  mcIdType ret=0;
  const mcIdType *start=_discr_per_cell->getConstPointer();
  for(const mcIdType *w=start;w!=start+cellId;w++)
    ret+=_loc[*w].getNumberOfGaussPt();
  return ret;
}

/*!
 * This method do the assumption that there is no orphan cell. If there is an exception is thrown.
 * This method makes the assumption too that '_discr_per_cell' is defined. If not an exception is thrown.
 * This method returns a newly created array with number of tuples equals to '_discr_per_cell->getNumberOfTuples' and number of components equal to 1.
 * The i_th tuple in returned array is the number of gauss point if the corresponding cell.
 */
DataArrayIdType *MEDCouplingFieldDiscretizationGauss::buildNbOfGaussPointPerCellField() const
{
  if(!_discr_per_cell)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGauss::buildNbOfGaussPointPerCellField : no discretization array set !");
  mcIdType nbOfTuples=_discr_per_cell->getNumberOfTuples();
  MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
  const mcIdType *w=_discr_per_cell->begin();
  ret->alloc(nbOfTuples,1);
  mcIdType *valsToFill=ret->getPointer();
  mcIdType nbMaxOfLocId=ToIdType(_loc.size());
  for(mcIdType i=0;i<nbOfTuples;i++,w++)
    if(*w!=DFT_INVALID_LOCID_VALUE)
      {
        if(*w>=0 && *w<nbMaxOfLocId)
          valsToFill[i]=_loc[*w].getNumberOfGaussPt();
        else
          {
            std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::buildNbOfGaussPointPerCellField : cell #" << i << " has invalid id (" << *w << ") ! Should be in [0," << nbMaxOfLocId << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    else
      {
        std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGauss::buildNbOfGaussPointPerCellField : cell #" << i << " is detected as orphan !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  return ret.retn();
}

void MEDCouplingFieldDiscretizationGauss::reprQuickOverview(std::ostream& stream) const
{
  stream << "Gauss points spatial discretization.";
}

/*!
 * This method makes the assumption that _discr_per_cell is set.
 * This method reduces as much as possible number size of _loc.
 * This method is useful when several set on same cells has been done and that some Gauss Localization are no more used.
 */
void MEDCouplingFieldDiscretizationGauss::zipGaussLocalizations()
{
  const mcIdType *start=_discr_per_cell->begin();
  mcIdType nbOfTuples=_discr_per_cell->getNumberOfTuples();
  INTERP_KERNEL::AutoPtr<mcIdType> tmp=new mcIdType[_loc.size()];
  std::fill((mcIdType *)tmp,(mcIdType *)tmp+_loc.size(),-2);
  for(const mcIdType *w=start;w!=start+nbOfTuples;w++)
    if(*w>=0)
      tmp[*w]=1;
  mcIdType fid=0;
  for(mcIdType i=0;i<ToIdType(_loc.size());i++)
    if(tmp[i]!=-2)
      tmp[i]=fid++;
  if(fid==ToIdType(_loc.size()))
    return;
  // zip needed
  mcIdType *start2=_discr_per_cell->getPointer();
  for(mcIdType *w2=start2;w2!=start2+nbOfTuples;w2++)
    if(*w2>=0)
      *w2=tmp[*w2];
  std::vector<MEDCouplingGaussLocalization> tmpLoc;
  for(mcIdType i=0;i<ToIdType(_loc.size());i++)
    if(tmp[i]!=-2)
      tmpLoc.push_back(_loc[i]);
  _loc=tmpLoc;
}

void MEDCouplingFieldDiscretizationGauss::commonUnserialization(const std::vector<mcIdType>& tinyInfo)
{
  mcIdType nbOfLoc=tinyInfo[1];
  _loc.clear();
  mcIdType dim=tinyInfo[2];
  mcIdType delta=-1;
  if(nbOfLoc>0)
    delta=(ToIdType(tinyInfo.size())-3)/nbOfLoc;
  for(mcIdType i=0;i<nbOfLoc;i++)
    {
      std::vector<mcIdType> tmp(tinyInfo.begin()+3+i*delta,tinyInfo.begin()+3+(i+1)*delta);
      MEDCouplingGaussLocalization elt=MEDCouplingGaussLocalization::BuildNewInstanceFromTinyInfo(dim,tmp);
      _loc.push_back(elt);
    }
}

MEDCouplingFieldDiscretizationGaussNE::MEDCouplingFieldDiscretizationGaussNE()
{
}

TypeOfField MEDCouplingFieldDiscretizationGaussNE::getEnum() const
{
  return TYPE;
}

/*!
 * This method is simply called by MEDCouplingFieldDiscretization::deepCopy. It performs the deep copy of \a this.
 *
 * \sa MEDCouplingFieldDiscretization::deepCopy.
 */
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
  if(!other)
    {
      reason="other spatial discretization is NULL, and this spatial discretization (GaussNE) is defined.";
      return false;
    }
  const MEDCouplingFieldDiscretizationGaussNE *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationGaussNE *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason="Spatial discrtization of this is ON_GAUSS_NE, which is not the case of other.";
  return ret;
}

/*!
 * This method returns the number of tuples regarding exclusively the input code \b without \b using \b a \b mesh \b in \b input.
 * The input code coherency is also checked regarding spatial discretization of \a this.
 * If an incoherency is detected, an exception will be thrown. If the input code is coherent, the number of tuples expected is returned.
 * The number of tuples expected is equal to those to have a valid field lying on \a this and having a mesh fitting perfectly the input code (geometric type distribution).
 */
mcIdType MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuplesExpectedRegardingCode(const std::vector<mcIdType>& code, const std::vector<const DataArrayIdType *>& idsPerType) const
{
  if(code.size()%3!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuplesExpectedRegardingCode : invalid input code !");
  mcIdType nbOfSplit=ToIdType(idsPerType.size());
  mcIdType nbOfTypes=ToIdType(code.size()/3);
  mcIdType ret(0);
  for(mcIdType i=0;i<nbOfTypes;i++)
    {
      const INTERP_KERNEL::CellModel& cm(INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)code[3*i]));
      if(cm.isDynamic())
        {
          std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuplesExpectedRegardingCode : At pos #" << i << " the geometric type " << cm.getRepr() << " is dynamic ! There are not managed by GAUSS_NE field discretization !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
      mcIdType nbOfEltInChunk=code[3*i+1];
      if(nbOfEltInChunk<0)
        throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuplesExpectedRegardingCode : invalid input code ! presence of negative value in a type !");
      mcIdType pos=code[3*i+2];
      if(pos!=-1)
        {
          if(pos<0 || pos>=nbOfSplit)
            {
              std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuplesExpectedRegardingCode : input code points to pos " << pos << " in typeid " << i << " ! Should be in [0," << nbOfSplit << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
          const DataArrayIdType *ids(idsPerType[pos]);
          if(!ids || !ids->isAllocated() || ids->getNumberOfComponents()!=1 || ids->getNumberOfTuples()!=nbOfEltInChunk || ids->getMinValueInArray()<0)
            {
              std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuplesExpectedRegardingCode : input pfl chunck at pos " << pos << " should have " << i << " tuples and one component and with ids all >=0 !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      ret+=nbOfEltInChunk*ToIdType(cm.getNumberOfNodes());
    }
  return ret;
}

mcIdType MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::getNumberOfTuples : NULL input mesh !");
  mcIdType ret=0;
  mcIdType nbOfCells=mesh->getNumberOfCells();
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      if(cm.isDynamic())
        throw INTERP_KERNEL::Exception("Not implemented yet Gauss node on elements for polygons and polyedrons !");
      ret+=cm.getNumberOfNodes();
    }
  return ret;
}

mcIdType MEDCouplingFieldDiscretizationGaussNE::getNumberOfMeshPlaces(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::getNumberOfMeshPlaces : NULL input mesh !");
  return mesh->getNumberOfCells();
}

DataArrayIdType *MEDCouplingFieldDiscretizationGaussNE::getOffsetArr(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::getOffsetArr : NULL input mesh !");
  mcIdType nbOfTuples=mesh->getNumberOfCells();
  DataArrayIdType *ret=DataArrayIdType::New();
  ret->alloc(nbOfTuples+1,1);
  mcIdType *retPtr=ret->getPointer();
  retPtr[0]=0;
  for(mcIdType i=0;i<nbOfTuples;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      if(cm.isDynamic())
        throw INTERP_KERNEL::Exception("Not implemented yet Gauss node on elements for polygons and polyedrons !");
      retPtr[i+1]=retPtr[i]+cm.getNumberOfNodes();
    }
  return ret;
}

void MEDCouplingFieldDiscretizationGaussNE::renumberArraysForCell(const MEDCouplingMesh *mesh, const std::vector<DataArray *>& arrays,
                                                                  const mcIdType *old2NewBg, bool check)
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::renumberArraysForCell : NULL input mesh !");
  const mcIdType *array=old2NewBg;
  if(check)
    array=DataArrayIdType::CheckAndPreparePermutation(old2NewBg,old2NewBg+mesh->getNumberOfCells());
  mcIdType nbOfCells=mesh->getNumberOfCells();
  mcIdType nbOfTuples=getNumberOfTuples(mesh);
  mcIdType *array2=new mcIdType[nbOfTuples];//stores the final conversion array old2New to give to arrays in renumberInPlace.
  mcIdType *array3=new mcIdType[nbOfCells];//store for each cell in after renumbering the offset needed by each cell in new numbering.
  array3[0]=0;
  for(mcIdType i=1;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(ToIdType(std::distance(array,std::find(array,array+nbOfCells,i-1))));
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      array3[i]=array3[i-1]+cm.getNumberOfNodes();
    }
  mcIdType j=0;
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      for(mcIdType k=0;k<ToIdType(cm.getNumberOfNodes());k++,j++)
        array2[j]=array3[array[i]]+k;
    }
  delete [] array3;
  for(std::vector<DataArray *>::const_iterator it=arrays.begin();it!=arrays.end();it++)
    if(*it)
      (*it)->renumberInPlace(array2);
  delete [] array2;
  if(check)
    free(const_cast<mcIdType *>(array));
}

DataArrayDouble *MEDCouplingFieldDiscretizationGaussNE::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::getLocalizationOfDiscValues : NULL input mesh !");
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  MCAuto<MEDCouplingUMesh> umesh=mesh->buildUnstructured();//in general do nothing
  mcIdType nbOfTuples=getNumberOfTuples(umesh);
  int spaceDim=mesh->getSpaceDimension();
  ret->alloc(nbOfTuples,spaceDim);
  const double *coords=umesh->getCoords()->begin();
  const mcIdType *connI=umesh->getNodalConnectivityIndex()->getConstPointer();
  const mcIdType *conn=umesh->getNodalConnectivity()->getConstPointer();
  mcIdType nbCells=umesh->getNumberOfCells();
  double *retPtr=ret->getPointer();
  for(mcIdType i=0;i<nbCells;i++,connI++)
    for(const mcIdType *w=conn+connI[0]+1;w!=conn+connI[1];w++)
      if(*w>=0)
        retPtr=std::copy(coords+(*w)*spaceDim,coords+((*w)+1)*spaceDim,retPtr);
  return ret.retn();
}

/*!
 * Reimplemented from MEDCouplingFieldDiscretization::integral for performance reason. The default implementation is valid too for GAUSS_NE spatial discretization.
 */
void MEDCouplingFieldDiscretizationGaussNE::integral(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, bool isWAbs, double *res) const
{
  if(!mesh || !arr)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::integral : input mesh or array is null !");
  std::size_t nbOfCompo=arr->getNumberOfComponents();
  std::fill(res,res+nbOfCompo,0.);
  //
  MCAuto<MEDCouplingFieldDouble> vol=mesh->getMeasureField(isWAbs);
  std::set<INTERP_KERNEL::NormalizedCellType> types=mesh->getAllGeoTypes();
  MCAuto<DataArrayIdType> nbOfNodesPerCell=mesh->computeNbOfNodesPerCell();
  nbOfNodesPerCell->computeOffsetsFull();
  const double *arrPtr=arr->begin(),*volPtr=vol->getArray()->begin();
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator it=types.begin();it!=types.end();it++)
    {
      std::size_t wArrSz=-1;
      const double *wArr=GetWeightArrayFromGeometricType(*it,wArrSz);
      INTERP_KERNEL::AutoPtr<double> wArr2=new double[wArrSz];
      double sum=std::accumulate(wArr,wArr+wArrSz,0.);
      std::transform(wArr,wArr+wArrSz,(double *)wArr2,std::bind(std::multiplies<double>(),std::placeholders::_1,1./sum));        
      MCAuto<DataArrayIdType> ids=mesh->giveCellsWithType(*it);
      MCAuto<DataArrayIdType> ids2=ids->buildExplicitArrByRanges(nbOfNodesPerCell);
      const mcIdType *ptIds2=ids2->begin(),*ptIds=ids->begin();
      mcIdType nbOfCellsWithCurGeoType=ids->getNumberOfTuples();
      for(mcIdType i=0;i<nbOfCellsWithCurGeoType;i++,ptIds++,ptIds2+=wArrSz)
        {
          for(std::size_t k=0;k<nbOfCompo;k++)
            {
              double tmp=0.;
              for(std::size_t j=0;j<wArrSz;j++)
                tmp+=arrPtr[nbOfCompo*ptIds2[j]+k]*wArr2[j];
              res[k]+=tmp*volPtr[*ptIds];
            }
        }
    }
}

const double *MEDCouplingFieldDiscretizationGaussNE::GetWeightArrayFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType, std::size_t& lgth)
{
  switch(geoType)
  {
    case INTERP_KERNEL::NORM_POINT1:
      lgth=sizeof(FGP_POINT1)/sizeof(double);
      return FGP_POINT1;
    case INTERP_KERNEL::NORM_SEG2:
      lgth=sizeof(FGP_SEG2)/sizeof(double);
      return FGP_SEG2;
    case INTERP_KERNEL::NORM_SEG3:
      lgth=sizeof(FGP_SEG3)/sizeof(double);
      return FGP_SEG3;
    case INTERP_KERNEL::NORM_SEG4:
      lgth=sizeof(FGP_SEG4)/sizeof(double);
      return FGP_SEG4;
    case INTERP_KERNEL::NORM_TRI3:
      lgth=sizeof(FGP_TRI3)/sizeof(double);
      return FGP_TRI3;
    case INTERP_KERNEL::NORM_TRI6:
      lgth=sizeof(FGP_TRI6)/sizeof(double);
      return FGP_TRI6;
    case INTERP_KERNEL::NORM_TRI7:
      lgth=sizeof(FGP_TRI7)/sizeof(double);
      return FGP_TRI7;
    case INTERP_KERNEL::NORM_QUAD4:
      lgth=sizeof(FGP_QUAD4)/sizeof(double);
      return FGP_QUAD4;
    case INTERP_KERNEL::NORM_QUAD8:
      lgth=sizeof(FGP_QUAD8)/sizeof(double);
      return FGP_QUAD8;
    case INTERP_KERNEL::NORM_QUAD9:
      lgth=sizeof(FGP_QUAD9)/sizeof(double);
      return FGP_QUAD9;
    case INTERP_KERNEL::NORM_TETRA4:
      lgth=sizeof(FGP_TETRA4)/sizeof(double);
      return FGP_TETRA4;
    case INTERP_KERNEL::NORM_TETRA10:
      lgth=sizeof(FGP_TETRA10)/sizeof(double);
      return FGP_TETRA10;
    case INTERP_KERNEL::NORM_PENTA6:
      lgth=sizeof(FGP_PENTA6)/sizeof(double);
      return FGP_PENTA6;
    case INTERP_KERNEL::NORM_PENTA15:
      lgth=sizeof(FGP_PENTA15)/sizeof(double);
      return FGP_PENTA15;
    case INTERP_KERNEL::NORM_PENTA18:
      lgth=sizeof(FGP_PENTA18)/sizeof(double);
      return FGP_PENTA18;
    case INTERP_KERNEL::NORM_HEXA8:
      lgth=sizeof(FGP_HEXA8)/sizeof(double);
      return FGP_HEXA8;
    case INTERP_KERNEL::NORM_HEXA20:
      lgth=sizeof(FGP_HEXA20)/sizeof(double);
      return FGP_HEXA20;
    case INTERP_KERNEL::NORM_HEXA27:
      lgth=sizeof(FGP_HEXA27)/sizeof(double);
      return FGP_HEXA27;
    case INTERP_KERNEL::NORM_PYRA5:
      lgth=sizeof(FGP_PYRA5)/sizeof(double);
      return FGP_PYRA5;
    case INTERP_KERNEL::NORM_PYRA13:
      lgth=sizeof(FGP_PYRA13)/sizeof(double);
      return FGP_PYRA13;
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::GetWeightArrayFromGeometricType : only SEG[2,3,4], TRI[3,6,7], QUAD[4,9], TETRA[4,10], PENTA[6,15,18], HEXA[8,20,27], PYRA[5,13] supported !");
  }
}

const double *MEDCouplingFieldDiscretizationGaussNE::GetRefCoordsFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType, std::size_t& lgth)
{
  switch(geoType)
  {
    case INTERP_KERNEL::NORM_POINT1:
      lgth=0;
      return 0;
    case INTERP_KERNEL::NORM_SEG2:
      lgth=sizeof(REF_SEG2)/sizeof(double);
      return REF_SEG2;
    case INTERP_KERNEL::NORM_SEG3:
      lgth=sizeof(REF_SEG3)/sizeof(double);
      return REF_SEG3;
    case INTERP_KERNEL::NORM_SEG4:
      lgth=sizeof(REF_SEG4)/sizeof(double);
      return REF_SEG4;
    case INTERP_KERNEL::NORM_TRI3:
      lgth=sizeof(REF_TRI3)/sizeof(double);
      return REF_TRI3;
    case INTERP_KERNEL::NORM_TRI6:
      lgth=sizeof(REF_TRI6)/sizeof(double);
      return REF_TRI6;
    case INTERP_KERNEL::NORM_TRI7:
      lgth=sizeof(REF_TRI7)/sizeof(double);
      return REF_TRI7;
    case INTERP_KERNEL::NORM_QUAD4:
      lgth=sizeof(REF_QUAD4)/sizeof(double);
      return REF_QUAD4;
    case INTERP_KERNEL::NORM_QUAD8:
      lgth=sizeof(REF_QUAD8)/sizeof(double);
      return REF_QUAD8;
    case INTERP_KERNEL::NORM_QUAD9:
      lgth=sizeof(REF_QUAD9)/sizeof(double);
      return REF_QUAD9;
    case INTERP_KERNEL::NORM_TETRA4:
      lgth=sizeof(REF_TETRA4)/sizeof(double);
      return REF_TETRA4;
    case INTERP_KERNEL::NORM_TETRA10:
      lgth=sizeof(REF_TETRA10)/sizeof(double);
      return REF_TETRA10;
    case INTERP_KERNEL::NORM_PENTA6:
      lgth=sizeof(REF_PENTA6)/sizeof(double);
      return REF_PENTA6;
    case INTERP_KERNEL::NORM_PENTA15:
      lgth=sizeof(REF_PENTA15)/sizeof(double);
      return REF_PENTA15;
    case INTERP_KERNEL::NORM_PENTA18:
      lgth=sizeof(REF_PENTA18)/sizeof(double);
      return REF_PENTA18;
    case INTERP_KERNEL::NORM_HEXA8:
      lgth=sizeof(REF_HEXA8)/sizeof(double);
      return REF_HEXA8;
    case INTERP_KERNEL::NORM_HEXA20:
      lgth=sizeof(REF_HEXA20)/sizeof(double);
      return REF_HEXA20;
    case INTERP_KERNEL::NORM_HEXA27:
      lgth=sizeof(REF_HEXA27)/sizeof(double);
      return REF_HEXA27;
    case INTERP_KERNEL::NORM_PYRA5:
      lgth=sizeof(REF_PYRA5)/sizeof(double);
      return REF_PYRA5;
    case INTERP_KERNEL::NORM_PYRA13:
      lgth=sizeof(REF_PYRA13)/sizeof(double);
      return REF_PYRA13;
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::GetRefCoordsFromGeometricType : only SEG[2,3,4], TRI[3,6,7], QUAD[4,8,9], TETRA[4,10], PENTA[6,15,18], HEXA[8,20,27], PYRA[5,13] supported !");
  }
}

const double *MEDCouplingFieldDiscretizationGaussNE::GetLocsFromGeometricType(INTERP_KERNEL::NormalizedCellType geoType, std::size_t& lgth)
{
  switch(geoType)
  {
    case INTERP_KERNEL::NORM_POINT1:
      {
        lgth=0;
        return 0;
      }
    case INTERP_KERNEL::NORM_SEG2:
      {
        lgth=sizeof(LOC_SEG2)/sizeof(double);
        return LOC_SEG2;
      }
    case INTERP_KERNEL::NORM_SEG3:
      {
        lgth=sizeof(LOC_SEG3)/sizeof(double);
        return LOC_SEG3;
      }
    case INTERP_KERNEL::NORM_SEG4:
      {
        lgth=sizeof(LOC_SEG4)/sizeof(double);
        return LOC_SEG4;
      }
    case INTERP_KERNEL::NORM_TRI3:
      {
        lgth=sizeof(LOC_TRI3)/sizeof(double);
        return LOC_TRI3;
      }
    case INTERP_KERNEL::NORM_TRI6:
      {
        lgth=sizeof(LOC_TRI6)/sizeof(double);
        return LOC_TRI6;
      }
    case INTERP_KERNEL::NORM_TRI7:
      {
        lgth=sizeof(LOC_TRI7)/sizeof(double);
        return LOC_TRI7;
      }
    case INTERP_KERNEL::NORM_QUAD4:
      {
        lgth=sizeof(LOC_QUAD4)/sizeof(double);
        return LOC_QUAD4;
      }
    case INTERP_KERNEL::NORM_QUAD8:
      {
        lgth=sizeof(LOC_QUAD8)/sizeof(double);
        return LOC_QUAD8;
      }
    case INTERP_KERNEL::NORM_QUAD9:
      {
        lgth=sizeof(LOC_QUAD9)/sizeof(double);
        return LOC_QUAD9;
      }
    case INTERP_KERNEL::NORM_TETRA4:
      {
        lgth=sizeof(LOC_TETRA4)/sizeof(double);
        return LOC_TETRA4;
      }
    case INTERP_KERNEL::NORM_TETRA10:
      {
        lgth=sizeof(LOC_TETRA10)/sizeof(double);
        return LOC_TETRA10;
      }
    case INTERP_KERNEL::NORM_PENTA6:
      {
        lgth=sizeof(LOC_PENTA6)/sizeof(double);
        return LOC_PENTA6;
      }
    case INTERP_KERNEL::NORM_PENTA15:
      {
        lgth=sizeof(LOC_PENTA15)/sizeof(double);
        return LOC_PENTA15;
      }
    case INTERP_KERNEL::NORM_PENTA18:
      {
        lgth=sizeof(LOC_PENTA18)/sizeof(double);
        return LOC_PENTA18;
      }
    case INTERP_KERNEL::NORM_HEXA8:
      {
        lgth=sizeof(LOC_HEXA8)/sizeof(double);
        return LOC_HEXA8;
      }
    case INTERP_KERNEL::NORM_HEXA20:
      {
        lgth=sizeof(LOC_HEXA20)/sizeof(double);
        return LOC_HEXA20;
      }
    case INTERP_KERNEL::NORM_HEXA27:
      {
        lgth=sizeof(LOC_HEXA27)/sizeof(double);
        return LOC_HEXA27;
      }
    case INTERP_KERNEL::NORM_PYRA5:
      {
        lgth=sizeof(LOC_PYRA5)/sizeof(double);
        return LOC_PYRA5;
      }
    case INTERP_KERNEL::NORM_PYRA13:
      {
        lgth=sizeof(LOC_PYRA13)/sizeof(double);
        return LOC_PYRA13;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::GetLocsFromGeometricType : only SEG[2,3,4], TRI[3,6,7], QUAD[4,8,9], TETRA[4,10], PENTA[6,15,18], HEXA[8,20,27], PYRA[5,13] supported !");
  }
}

void MEDCouplingFieldDiscretizationGaussNE::computeMeshRestrictionFromTupleIds(const MEDCouplingMesh *mesh, const mcIdType *tupleIdsBg, const mcIdType *tupleIdsEnd,
                                                                               DataArrayIdType *&cellRestriction, DataArrayIdType *&trueTupleRestriction) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::computeMeshRestrictionFromTupleIds : NULL input mesh !");
  MCAuto<DataArrayIdType> tmp=DataArrayIdType::New(); tmp->alloc(std::distance(tupleIdsBg,tupleIdsEnd),1);
  std::copy(tupleIdsBg,tupleIdsEnd,tmp->getPointer());
  tmp->sort(true);
  tmp=tmp->buildUnique();
  MCAuto<DataArrayIdType> nbOfNodesPerCell=mesh->computeNbOfNodesPerCell();
  nbOfNodesPerCell->computeOffsetsFull();
  nbOfNodesPerCell->findIdsRangesInListOfIds(tmp,cellRestriction,trueTupleRestriction);
}

void MEDCouplingFieldDiscretizationGaussNE::checkCompatibilityWithNature(NatureOfField nat) const
{
}

double MEDCouplingFieldDiscretizationGaussNE::getIJK(const MEDCouplingMesh *mesh, const DataArrayDouble *da, mcIdType cellId, mcIdType nodeIdInCell, int compoId) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::getIJK : NULL input mesh !");
  mcIdType offset=0;
  for(mcIdType i=0;i<cellId;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      offset+=cm.getNumberOfNodes();
    }
  return da->getIJ(offset+nodeIdInCell,compoId);
}

void MEDCouplingFieldDiscretizationGaussNE::checkCoherencyBetween(const MEDCouplingMesh *mesh, const DataArray *da) const
{
  mcIdType nbOfTuples(getNumberOfTuples(mesh));
  if(nbOfTuples!=da->getNumberOfTuples())
    {
      std::ostringstream oss; oss << "Invalid number of tuples in the array : expecting " << nbOfTuples << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationGaussNE::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::getMeasureField : mesh instance specified is NULL !");
  MCAuto<MEDCouplingFieldDouble> vol=mesh->getMeasureField(isAbs);
  const double *volPtr=vol->getArray()->begin();
  MCAuto<MEDCouplingFieldDouble> ret=MEDCouplingFieldDouble::New(ON_GAUSS_NE);
  ret->setMesh(mesh);
  //
  std::set<INTERP_KERNEL::NormalizedCellType> types=mesh->getAllGeoTypes();
  MCAuto<DataArrayIdType> nbOfNodesPerCell=mesh->computeNbOfNodesPerCell();
  mcIdType nbTuples=nbOfNodesPerCell->accumulate((std::size_t)0);
  nbOfNodesPerCell->computeOffsetsFull();
  MCAuto<DataArrayDouble> arr=DataArrayDouble::New(); arr->alloc(nbTuples,1);
  ret->setArray(arr);
  double *arrPtr=arr->getPointer();
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator it=types.begin();it!=types.end();it++)
    {
      std::size_t wArrSz=-1;
      const double *wArr=GetWeightArrayFromGeometricType(*it,wArrSz);
      INTERP_KERNEL::AutoPtr<double> wArr2=new double[wArrSz];
      double sum=std::accumulate(wArr,wArr+wArrSz,0.);
      std::transform(wArr,wArr+wArrSz,(double *)wArr2,std::bind(std::multiplies<double>(),std::placeholders::_1,1./sum));     
      MCAuto<DataArrayIdType> ids=mesh->giveCellsWithType(*it);
      MCAuto<DataArrayIdType> ids2=ids->buildExplicitArrByRanges(nbOfNodesPerCell);
      const mcIdType *ptIds2=ids2->begin(),*ptIds=ids->begin();
      mcIdType nbOfCellsWithCurGeoType=ids->getNumberOfTuples();
      for(mcIdType i=0;i<nbOfCellsWithCurGeoType;i++,ptIds++)
        for(std::size_t j=0;j<wArrSz;j++,ptIds2++)
          arrPtr[*ptIds2]=wArr2[j]*volPtr[*ptIds];
    }
  ret->synchronizeTimeWithSupport();
  return ret.retn();
}

void MEDCouplingFieldDiscretizationGaussNE::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGaussNE::getValueOnPos(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, mcIdType i, mcIdType j, mcIdType k, double *res) const
{
  throw INTERP_KERNEL::Exception("getValueOnPos(i,j,k) : Not applicable for Gauss points !");
}

DataArrayDouble *MEDCouplingFieldDiscretizationGaussNE::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const
{
  throw INTERP_KERNEL::Exception("getValueOnMulti : Not implemented for Gauss NE !");
}

MEDCouplingMesh *MEDCouplingFieldDiscretizationGaussNE::buildSubMeshData(const MEDCouplingMesh *mesh, const mcIdType *start, const mcIdType *end, DataArrayIdType *&di) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::buildSubMeshData : NULL input mesh !");
  MCAuto<DataArrayIdType> diSafe=computeTupleIdsToSelectFromCellIds(mesh,start,end);
  MCAuto<MEDCouplingMesh> ret=mesh->buildPart(start,end);
  di=diSafe.retn();
  return ret.retn();
}

/*!
 * This method is strictly equivalent to MEDCouplingFieldDiscretizationGauss::buildSubMeshData except that it is optimized for input defined as a range of cell ids.
 * 
 * \param [out] beginOut Valid only if \a di is NULL
 * \param [out] endOut Valid only if \a di is NULL
 * \param [out] stepOut Valid only if \a di is NULL
 * \param [out] di is an array returned that specifies entity ids (nodes, cells, Gauss points... ) in array if no output range is foundable.
 *
 * \sa MEDCouplingFieldDiscretizationGauss::buildSubMeshData
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationGaussNE::buildSubMeshDataRange(const MEDCouplingMesh *mesh, mcIdType beginCellIds, mcIdType endCellIds, mcIdType stepCellIds, mcIdType& beginOut, mcIdType& endOut, mcIdType& stepOut, DataArrayIdType *&di) const
{
  if(stepCellIds!=1)//even for stepCellIds==-1 the output will not be a range
    return MEDCouplingFieldDiscretization::buildSubMeshDataRange(mesh,beginCellIds,endCellIds,stepCellIds,beginOut,endOut,stepOut,di);
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::buildSubMeshDataRange : NULL input mesh !");
  mcIdType nbOfCells=mesh->getNumberOfCells();
  di=0; beginOut=0; endOut=0; stepOut=stepCellIds;
  const char msg[]="MEDCouplingFieldDiscretizationGaussNE::buildSubMeshDataRange : cell #";
  for(mcIdType i=0;i<nbOfCells;i++)
    {
      INTERP_KERNEL::NormalizedCellType type=mesh->getTypeOfCell(i);
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(type);
      if(cm.isDynamic())
        { std::ostringstream oss; oss << msg << i << " presence of dynamic cell (polygons and polyedrons) ! Not implemented !"; throw INTERP_KERNEL::Exception(oss.str().c_str()); }
      mcIdType delta=cm.getNumberOfNodes();
      if(i<beginCellIds)
        beginOut+=delta;
      endOut+=delta;
      if(i>=endCellIds)
        break;
    }
  MCAuto<MEDCouplingMesh> ret=mesh->buildPartRange(beginCellIds,endCellIds,stepCellIds);
  return ret.retn();
}


/*!
 * This method returns a tuple ids selection from cell ids selection [start;end).
 * This method is called by MEDCouplingFieldDiscretizationGaussNE::buildSubMeshData to return parameter \b di.
 *
 * \return a newly allocated array containing ids to select into the DataArrayDouble of the field.
 * 
 */
DataArrayIdType *MEDCouplingFieldDiscretizationGaussNE::computeTupleIdsToSelectFromCellIds(const MEDCouplingMesh *mesh, const mcIdType *startCellIds, const mcIdType *endCellIds) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationGaussNE::computeTupleIdsToSelectFromCellIds : null mesh !");
  MCAuto<DataArrayIdType> nbOfNodesPerCell=mesh->computeNbOfNodesPerCell();
  nbOfNodesPerCell->computeOffsetsFull();
  MCAuto<DataArrayIdType> sel=DataArrayIdType::New(); sel->useArray(startCellIds,false,DeallocType::CPP_DEALLOC,std::distance(startCellIds,endCellIds),1);
  return sel->buildExplicitArrByRanges(nbOfNodesPerCell);
}

/*!
 * No implementation needed !
 */
void MEDCouplingFieldDiscretizationGaussNE::renumberValuesOnNodes(double , const mcIdType *, mcIdType newNbOfNodes, DataArrayDouble *) const
{
}

void MEDCouplingFieldDiscretizationGaussNE::renumberValuesOnCells(double epsOnVals, const MEDCouplingMesh *mesh, const mcIdType *old2New, mcIdType newSz, DataArrayDouble *arr) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

MCAuto<MEDCouplingFieldDiscretization> MEDCouplingFieldDiscretizationGaussNE::aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const
{
  return EasyAggregate<MEDCouplingFieldDiscretizationGaussNE>(fds);
}

void MEDCouplingFieldDiscretizationGaussNE::renumberValuesOnCellsR(const MEDCouplingMesh *mesh, const mcIdType *new2old, mcIdType newSz, DataArrayDouble *arr) const
{
  throw INTERP_KERNEL::Exception("Not implemented yet !");
}

void MEDCouplingFieldDiscretizationGaussNE::reprQuickOverview(std::ostream& stream) const
{
  stream << "Gauss points on nodes per element spatial discretization.";
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

/*!
 * This method is simply called by MEDCouplingFieldDiscretization::deepCopy. It performs the deep copy of \a this.
 *
 * \sa MEDCouplingFieldDiscretization::deepCopy.
 */
MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationKriging::clone() const
{
  return new MEDCouplingFieldDiscretizationKriging;
}

std::string MEDCouplingFieldDiscretizationKriging::getStringRepr() const
{
  return std::string(REPR);
}

void MEDCouplingFieldDiscretizationKriging::checkCompatibilityWithNature(NatureOfField nat) const
{
  if(nat!=IntensiveMaximum)
    throw INTERP_KERNEL::Exception("Invalid nature for Kriging field : expected IntensiveMaximum !");
}

bool MEDCouplingFieldDiscretizationKriging::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  if(!other)
    {
      reason="other spatial discretization is NULL, and this spatial discretization (Kriginig) is defined.";
      return false;
    }
  const MEDCouplingFieldDiscretizationKriging *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationKriging *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason="Spatial discrtization of this is ON_NODES_KR, which is not the case of other.";
  return ret;
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationKriging::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::getMeasureField : mesh instance specified is NULL !");
  throw INTERP_KERNEL::Exception("getMeasureField on FieldDiscretizationKriging : not implemented yet !");
}

void MEDCouplingFieldDiscretizationKriging::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, double *res) const
{
  MCAuto<DataArrayDouble> res2=MEDCouplingFieldDiscretizationKriging::getValueOnMulti(arr,mesh,loc,1);
  std::copy(res2->begin(),res2->end(),res);
}

DataArrayDouble *MEDCouplingFieldDiscretizationKriging::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfTargetPoints) const
{
  if(!arr || !arr->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::getValueOnMulti : input array is null or not allocated !");
  mcIdType nbOfRows=getNumberOfMeshPlaces(mesh);
  if(arr->getNumberOfTuples()!=nbOfRows)
    {
      std::ostringstream oss; oss << "MEDCouplingFieldDiscretizationKriging::getValueOnMulti : input array does not have correct number of tuples ! Excepted " << nbOfRows << " having " << arr->getNumberOfTuples() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  mcIdType nbCols(-1);
  std::size_t nbCompo=arr->getNumberOfComponents();
  MCAuto<DataArrayDouble> m(computeEvaluationMatrixOnGivenPts(mesh,loc,nbOfTargetPoints,nbCols));
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbOfTargetPoints,nbCompo);
  INTERP_KERNEL::matrixProduct(m->begin(),nbOfTargetPoints,nbCols,arr->begin(),nbOfRows,ToIdType(nbCompo),ret->getPointer());
  return ret.retn();
}

void MEDCouplingFieldDiscretizationKriging::reprQuickOverview(std::ostream& stream) const
{
  stream << "Kriging spatial discretization.";
}

MCAuto<MEDCouplingFieldDiscretization> MEDCouplingFieldDiscretizationKriging::aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const
{
  return EasyAggregate<MEDCouplingFieldDiscretizationKriging>(fds);
}

/*!
 * Returns the matrix of size nbRows = \a nbOfTargetPoints and \a nbCols = \a nbCols. This matrix is useful if 
 * 
 * \return the new result matrix to be deallocated by the caller.
 */
DataArrayDouble *MEDCouplingFieldDiscretizationKriging::computeEvaluationMatrixOnGivenPts(const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfTargetPoints, mcIdType& nbCols) const
{
  mcIdType isDrift(-1),nbRows(-1);
  MCAuto<DataArrayDouble> matrixInv(computeInverseMatrix(mesh,isDrift,nbRows));
  //
  MCAuto<DataArrayDouble> coords=getLocalizationOfDiscValues(mesh);
  mcIdType nbOfPts(coords->getNumberOfTuples());
  std::size_t dimension(coords->getNumberOfComponents());
  MCAuto<DataArrayDouble> locArr=DataArrayDouble::New();
  locArr->useArray(loc,false,DeallocType::CPP_DEALLOC,nbOfTargetPoints,dimension);
  nbCols=nbOfPts;
  //
  MCAuto<DataArrayDouble> matrix2=coords->buildEuclidianDistanceDenseMatrixWith(locArr);
  operateOnDenseMatrix(mesh->getSpaceDimension(),nbOfTargetPoints*nbOfPts,matrix2->getPointer());
  //
  MCAuto<DataArrayDouble> matrix3=DataArrayDouble::New();
  matrix3->alloc(nbOfTargetPoints*nbRows,1);
  double *work=matrix3->getPointer();
  const double *workCst(matrix2->begin()),*workCst2(loc);
  for(mcIdType i=0;i<nbOfTargetPoints;i++,workCst+=nbOfPts,workCst2+=isDrift-1)
    {
      for(mcIdType j=0;j<nbOfPts;j++)
        work[i*nbRows+j]=workCst[j];
      work[i*nbRows+nbOfPts]=1.0;
      for(mcIdType j=0;j<isDrift-1;j++)
        work[i*nbRows+(nbOfPts+1+j)]=workCst2[j];
    }
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbOfTargetPoints,nbRows);
  INTERP_KERNEL::matrixProduct(matrix3->begin(),nbOfTargetPoints,nbRows,matrixInv->begin(),nbRows,nbRows,ret->getPointer());
  MCAuto<DataArrayDouble> ret2(DataArrayDouble::New());
  ret2->alloc(nbOfTargetPoints*nbOfPts,1);
  workCst=ret->begin(); work=ret2->getPointer();
  for(mcIdType i=0;i<nbOfTargetPoints;i++,workCst+=nbRows)
    work=std::copy(workCst,workCst+nbOfPts,work);
  return ret2.retn();
}

/*!
 * This method returns the square matrix of size \a matSz that is the inverse of the kriging matrix. The returned matrix can returned all the coeffs of kriging
 * when multiplied by the vector of values attached to each point.
 * 
 * \param [out] isDrift return if drift coefficients are present in the returned vector of coefficients. If different from 0 there is presence of drift coefficients.
 * \param [out] matSz the size of returned square matrix
 * \return the new result matrix to be deallocated by the caller.
 * \sa computeMatrix
 */
DataArrayDouble *MEDCouplingFieldDiscretizationKriging::computeInverseMatrix(const MEDCouplingMesh *mesh, mcIdType& isDrift, mcIdType& matSz) const
{
  MCAuto<DataArrayDouble> matrixWithDrift(computeMatrix(mesh,isDrift,matSz));
  MCAuto<DataArrayDouble> matrixInv(DataArrayDouble::New());
  matrixInv->alloc(matSz*matSz,1);
  INTERP_KERNEL::inverseMatrix(matrixWithDrift->getConstPointer(),matSz,matrixInv->getPointer());
  return matrixInv.retn();
}

/*!
 * This method computes the kriging matrix.
 * \return the new result matrix to be deallocated by the caller.
 * \sa computeInverseMatrix
 */
DataArrayDouble *MEDCouplingFieldDiscretizationKriging::computeMatrix(const MEDCouplingMesh *mesh, mcIdType& isDrift, mcIdType& matSz) const
{
  if(!mesh)
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::computeMatrix : NULL input mesh !");
    MCAuto<DataArrayDouble> coords(getLocalizationOfDiscValues(mesh));
    mcIdType nbOfPts(coords->getNumberOfTuples());
    MCAuto<DataArrayDouble> matrix(coords->buildEuclidianDistanceDenseMatrix());
    operateOnDenseMatrix(mesh->getSpaceDimension(),nbOfPts*nbOfPts,matrix->getPointer());
    // Drift
    MCAuto<DataArrayDouble> matrixWithDrift(performDrift(matrix,coords,isDrift));
    matSz=nbOfPts+isDrift;
    return matrixWithDrift.retn();
}

/*!
 * This method computes coefficients to apply to each representing points of \a mesh, that is to say the nodes of \a mesh given a field array \a arr whose
 * number of tuples should be equal to the number of representing points in \a mesh.
 * 
 * \param [in] mesh is the sources of nodes on which kriging will be done regarding the parameters and the value of \c this->getSpaceDimension()
 * \param [in] arr input field DataArrayDouble whose number of tuples must be equal to the number of nodes in \a mesh
 * \param [out] isDrift return if drift coefficients are present in the returned vector of coefficients. If different from 0 there is presence of drift coefficients.
 *              Whatever the value of \a isDrift the number of tuples of returned DataArrayDouble  will be equal to \c arr->getNumberOfTuples() + \a isDrift.
 * \return a newly allocated array containing coefficients including or not drift coefficient at the end depending the value of \a isDrift parameter.
 */
DataArrayDouble *MEDCouplingFieldDiscretizationKriging::computeVectorOfCoefficients(const MEDCouplingMesh *mesh, const DataArrayDouble *arr, mcIdType& isDrift) const
{
  mcIdType nbRows(-1);
  MCAuto<DataArrayDouble> matrixInv(computeInverseMatrix(mesh,isDrift,nbRows));
  MCAuto<DataArrayDouble> KnewiK(DataArrayDouble::New());
  KnewiK->alloc(nbRows*1,1);
  MCAuto<DataArrayDouble> arr2(PerformDriftOfVec(arr,isDrift));
  INTERP_KERNEL::matrixProduct(matrixInv->getConstPointer(),nbRows,nbRows,arr2->getConstPointer(),arr2->getNumberOfTuples(),1,KnewiK->getPointer());
  return KnewiK.retn();
}

/*!
 * Apply \f f(x) on each element x in \a matrixPtr. \a matrixPtr is expected to be a dense matrix represented by a chunck of memory of size at least equal to \a nbOfElems.
 *
 * \param [in] spaceDimension space dimension of the input mesh on which the Kriging has to be performed
 * \param [in] nbOfElems is the result of the product of nb of rows and the nb of columns of matrix \a matrixPtr
 * \param [in,out] matrixPtr is the dense matrix whose on each values the operation will be applied
 */
void MEDCouplingFieldDiscretizationKriging::operateOnDenseMatrix(int spaceDimension, mcIdType nbOfElems, double *matrixPtr) const
{
  switch(spaceDimension)
  {
    case 1:
      {
        OperateOnDenseMatrixH3(nbOfElems,matrixPtr);
        break;
      }
    case 2:
      {
        OperateOnDenseMatrixH2Ln(nbOfElems,matrixPtr);
        break;
      }
    case 3:
      {
        //nothing here : it is not a bug g(h)=h with spaceDim 3.
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::operateOnDenseMatrix : only dimension 1, 2 and 3 implemented !");
  }
}

void MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH3(mcIdType nbOfElems, double *matrixPtr)
{
  for(mcIdType i=0;i<nbOfElems;i++)
    {
      double val=matrixPtr[i];
      matrixPtr[i]=val*val*val;
    }
}

void MEDCouplingFieldDiscretizationKriging::OperateOnDenseMatrixH2Ln(mcIdType nbOfElems, double *matrixPtr)
{
  for(mcIdType i=0;i<nbOfElems;i++)
    {
      double val=matrixPtr[i];
      if(val!=0.)
        matrixPtr[i]=val*val*log(val);
    }
}

/*!
 * Performs a drift to the rectangular input matrix \a matr.
 * This method generate a dense matrix starting from an input dense matrix \a matr and input array \a arr.
 * \param [in] matr The rectangular dense matrix (with only one component). The number of rows of \a matr must be equal to the number of tuples of \a arr
 * \param [in] arr The array of coords to be appended in the input dense matrix \a matr. Typically arr is an array of coordinates.
 * \param [out] delta the delta of number of columns between returned dense matrix and input dense matrix \a matr. \a delta is equal to number of components of \a arr + 1.
 * \sa performDrift
 */
DataArrayDouble *MEDCouplingFieldDiscretizationKriging::PerformDriftRect(const DataArrayDouble *matr, const DataArrayDouble *arr, mcIdType& delta)
{
  if(!matr || !matr->isAllocated() || matr->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::PerformDriftRect : invalid input dense matrix ! Must be allocated not NULL and with exactly one component !");
  if(!arr || !arr->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::PerformDriftRect : invalid input array of coordiantes ! Must be allocated and not NULL !");
  std::size_t spaceDimension(arr->getNumberOfComponents());
  mcIdType nbOfPts(arr->getNumberOfTuples()),nbOfEltInMatrx(matr->getNumberOfTuples());
  delta=ToIdType(spaceDimension)+1;
  mcIdType nbOfCols(nbOfEltInMatrx/nbOfPts);
  if(nbOfEltInMatrx%nbOfPts!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::PerformDriftRect : size of input dense matrix and input arrays mismatch ! NbOfElems in matrix % nb of tuples in array must be equal to 0 !");
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(nbOfPts*(nbOfCols+delta));
  double *retPtr(ret->getPointer());
  const double *mPtr(matr->begin()),*aPtr(arr->begin());
  for(mcIdType i=0;i<nbOfPts;i++,aPtr+=spaceDimension,mPtr+=nbOfCols)
    {
      retPtr=std::copy(mPtr,mPtr+nbOfCols,retPtr);
      *retPtr++=1.;
      retPtr=std::copy(aPtr,aPtr+spaceDimension,retPtr);
    }
  return ret.retn();
}

/*!
 * \return a newly allocated array having \a isDrift more tuples than \a arr.
 * \sa computeVectorOfCoefficients
 */
DataArrayDouble *MEDCouplingFieldDiscretizationKriging::PerformDriftOfVec(const DataArrayDouble *arr, mcIdType isDrift)
{
  if(!arr || !arr->isAllocated() || arr->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::PerformDriftOfVec : input array must be not NULL allocated and with one component !");
  if(isDrift<0)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::PerformDriftOfVec : isDrift parameter must be >=0 !");
  MCAuto<DataArrayDouble> arr2(DataArrayDouble::New());
  arr2->alloc((arr->getNumberOfTuples()+isDrift)*1,1);
  double *work(std::copy(arr->begin(),arr->end(),arr2->getPointer()));
  std::fill(work,work+isDrift,0.);
  return arr2.retn();
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
 * \sa MEDCouplingFieldDiscretizationKriging::PerformDriftRect
 */
DataArrayDouble *MEDCouplingFieldDiscretizationKriging::performDrift(const DataArrayDouble *matr, const DataArrayDouble *arr, mcIdType& delta) const
{
  std::size_t spaceDimension(arr->getNumberOfComponents());
  delta=ToIdType(spaceDimension)+1;
  mcIdType szOfMatrix(arr->getNumberOfTuples());
  if(szOfMatrix*szOfMatrix!=matr->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationKriging::performDrift : invalid size");
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  ret->alloc((szOfMatrix+delta)*(szOfMatrix+delta),1);
  const double *srcWork=matr->getConstPointer();
  const double *srcWork2=arr->getConstPointer();
  double *destWork=ret->getPointer();
  for(mcIdType i=0;i<szOfMatrix;i++)
    {
      destWork=std::copy(srcWork,srcWork+szOfMatrix,destWork);
      srcWork+=szOfMatrix;
      *destWork++=1.;
      destWork=std::copy(srcWork2,srcWork2+spaceDimension,destWork);
      srcWork2+=spaceDimension;
    }
  std::fill(destWork,destWork+szOfMatrix,1.); destWork+=szOfMatrix;
  std::fill(destWork,destWork+spaceDimension+1,0.); destWork+=spaceDimension+1;
  MCAuto<DataArrayDouble> arrNoI=arr->toNoInterlace();
  srcWork2=arrNoI->getConstPointer();
  for(std::size_t i=0;i<spaceDimension;i++)
    {
      destWork=std::copy(srcWork2,srcWork2+szOfMatrix,destWork);
      srcWork2+=szOfMatrix;
      std::fill(destWork,destWork+spaceDimension+1,0.);
      destWork+=spaceDimension+1;
    }
  //
  return ret.retn();
}
