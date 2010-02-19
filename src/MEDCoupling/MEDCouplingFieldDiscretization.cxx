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
#include "MEDCouplingFieldDiscretization.hxx"
#include "MEDCouplingPointSet.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <limits>
#include <algorithm>

using namespace ParaMEDMEM;

const char MEDCouplingFieldDiscretizationP0::REPR[]="P0";

const TypeOfField MEDCouplingFieldDiscretizationP0::TYPE=ON_CELLS;

const char MEDCouplingFieldDiscretizationP1::REPR[]="P1";

const TypeOfField MEDCouplingFieldDiscretizationP1::TYPE=ON_NODES;

MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretization::New(TypeOfField type)
{
  switch(type)
    {
    case MEDCouplingFieldDiscretizationP0::TYPE:
      return new MEDCouplingFieldDiscretizationP0;
    case MEDCouplingFieldDiscretizationP1::TYPE:
      return new MEDCouplingFieldDiscretizationP1;
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
  throw INTERP_KERNEL::Exception("Representation does not match with any field discretization !");
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

bool MEDCouplingFieldDiscretizationP0::isEqual(const MEDCouplingFieldDiscretization *other) const
{
  const MEDCouplingFieldDiscretizationP0 *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationP0 *>(other);
  return otherC!=0;
}

int MEDCouplingFieldDiscretizationP0::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfCells();
}

DataArrayDouble *MEDCouplingFieldDiscretizationP0::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  return mesh->getBarycenterAndOwner();
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

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationP0::getWeightingField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  return mesh->getMeasureField(isAbs);
}

/*!
 * Nothing to do. It's not a bug.
 */
void MEDCouplingFieldDiscretizationP0::renumberValuesOnNodes(const DataArrayInt *old2New, DataArrayDouble *arr) const
{
}

/*!
 * This method returns a submesh of 'mesh' instance constituting cell ids contained in array defined as an interval [start;end).
 * @ param di is an array returned that specifies entity ids (here cells ids) in mesh 'mesh' of entity in returned submesh.
 * Example : The first cell id of returned mesh has the (*di)[0] id in 'mesh'
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationP0::buildSubMeshData(const int *start, const int *end, const MEDCouplingMesh *mesh, DataArrayInt *&di) const
{
  MEDCouplingPointSet* ret=((const MEDCouplingPointSet *) mesh)->buildPartOfMySelf(start,end,false);
  di=DataArrayInt::New();
  di->alloc(end-start,1);
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

bool MEDCouplingFieldDiscretizationP1::isEqual(const MEDCouplingFieldDiscretization *other) const
{
  const MEDCouplingFieldDiscretizationP1 *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationP1 *>(other);
  return otherC!=0;
}

int MEDCouplingFieldDiscretizationP1::getNumberOfTuples(const MEDCouplingMesh *mesh) const
{
  return mesh->getNumberOfNodes();
}

DataArrayDouble *MEDCouplingFieldDiscretizationP1::getLocalizationOfDiscValues(const MEDCouplingMesh *mesh) const
{
  return mesh->getCoordinatesAndOwner();
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
      message << " cells in mesh and " << da->getNumberOfTuples() << " tuples in field !";
      throw INTERP_KERNEL::Exception(message.str().c_str());
    }
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationP1::getWeightingField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  return mesh->getMeasureFieldOnNode(isAbs);
}

void MEDCouplingFieldDiscretizationP1::renumberValuesOnNodes(const DataArrayInt *old2New, DataArrayDouble *arr) const
{
  int oldNbOfElems=old2New->getNbOfElems();
  const int *old2NewPtr=old2New->getConstPointer();
  int nbOfComp=arr->getNumberOfComponents();
  int newNbOfTuples=(*std::max_element(old2NewPtr,old2NewPtr+oldNbOfElems))+1;
  DataArrayDouble *arrCpy=arr->deepCopy();
  const double *ptSrc=arrCpy->getConstPointer();
  arr->reAlloc(newNbOfTuples);
  double *ptToFill=arr->getPointer();
  std::fill(ptToFill,ptToFill+nbOfComp*newNbOfTuples,std::numeric_limits<double>::max());
  for(int i=0;i<oldNbOfElems;i++)
    {
      int newNb=old2NewPtr[i];
      if(std::find_if(ptToFill+newNb*nbOfComp,ptToFill+(newNb+1)*nbOfComp,std::bind2nd(std::not_equal_to<double>(),std::numeric_limits<double>::max()))
         ==ptToFill+(newNb+1)*nbOfComp)
        std::copy(ptSrc+i*nbOfComp,ptSrc+(i+1)*nbOfComp,ptToFill+newNb*nbOfComp);
      else
        {
          if(!std::equal(ptSrc+i*nbOfComp,ptSrc+(i+1)*nbOfComp,ptToFill+newNb*nbOfComp))
            {
              arrCpy->decrRef();
              std::ostringstream oss;
              oss << "Node " << i << " and " << std::find(old2NewPtr,old2NewPtr+i,newNb)-old2NewPtr
                  << " have been merged and nodal field on them are different !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
    }
  arrCpy->decrRef();
}

/*!
 * This method invert array 'di' that is a conversion map from Old to New node numbering to New to Old node numbering.
 */
DataArrayInt *MEDCouplingFieldDiscretizationP1::invertArrayO2N2N2O(const MEDCouplingMesh *mesh, const DataArrayInt *di)
{
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(mesh->getNumberOfNodes(),1);
  int nbOfOldNodes=di->getNumberOfTuples();
  const int *old2New=di->getConstPointer();
  int *pt=ret->getPointer();
  for(int i=0;i!=nbOfOldNodes;i++)
    if(old2New[i]!=-1)
      pt[old2New[i]]=i;
  return ret;
}

/*!
 * This method returns a submesh of 'mesh' instance constituting cell ids contained in array defined as an interval [start;end).
* @ param di is an array returned that specifies entity ids (here nodes ids) in mesh 'mesh' of entity in returned submesh.
 * Example : The first node id of returned mesh has the (*di)[0] id in 'mesh'
 */
MEDCouplingMesh *MEDCouplingFieldDiscretizationP1::buildSubMeshData(const int *start, const int *end, const MEDCouplingMesh *mesh, DataArrayInt *&di) const
{
  MEDCouplingPointSet* ret=((const MEDCouplingPointSet *) mesh)->buildPartOfMySelf(start,end,true);
  DataArrayInt *diInv=ret->zipCoordsTraducer();
  di=invertArrayO2N2N2O(ret,diInv);
  diInv->decrRef();
  return ret;
}
