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
// Author : Anthony Geay (EDF R&D)

#include "MEDCouplingFieldDiscretizationOnNodesFE.hxx"
#include "MEDCouplingNormalizedUnstructuredMesh.txx"
#include "InterpKernelDenseMatrix.hxx"
#include "InterpKernelRootsMultiDim.hxx"
#include "MEDCouplingUMesh.hxx"
#include "InterpolationHelper.txx"
#include "InterpKernelGaussCoords.hxx"

#include <sstream>

using namespace MEDCoupling;

const char MEDCouplingFieldDiscretizationOnNodesFE::REPR[]="FE";

std::string MEDCouplingFieldDiscretizationOnNodesFE::getStringRepr() const
{
  return std::string(REPR);
}

void MEDCouplingFieldDiscretizationOnNodesFE::reprQuickOverview(std::ostream& stream) const
{
  stream << "NodeFE spatial discretization.";
}

MCAuto<MEDCouplingFieldDiscretization> MEDCouplingFieldDiscretizationOnNodesFE::aggregate(std::vector<const MEDCouplingFieldDiscretization *>& fds) const
{
  return EasyAggregate<MEDCouplingFieldDiscretizationOnNodesFE>(fds);
}

bool MEDCouplingFieldDiscretizationOnNodesFE::isEqualIfNotWhy(const MEDCouplingFieldDiscretization *other, double eps, std::string& reason) const
{
  if(!other)
    {
      reason="other spatial discretization is NULL, and this spatial discretization (Node FE) is defined.";
      return false;
    }
  const MEDCouplingFieldDiscretizationOnNodesFE *otherC=dynamic_cast<const MEDCouplingFieldDiscretizationOnNodesFE *>(other);
  bool ret=otherC!=0;
  if(!ret)
    reason="Spatial discrtization of this is ON_NODES_FE, which is not the case of other.";
  return ret;
}

/*!
 * This method is simply called by MEDCouplingFieldDiscretization::deepCopy. It performs the deep copy of \a this.
 *
 * \sa MEDCouplingFieldDiscretization::deepCopy.
 */
MEDCouplingFieldDiscretization *MEDCouplingFieldDiscretizationOnNodesFE::clone() const
{
  return new MEDCouplingFieldDiscretizationOnNodesFE;
}

void MEDCouplingFieldDiscretizationOnNodesFE::checkCompatibilityWithNature(NatureOfField nat) const
{
  if(nat!=IntensiveMaximum)
    throw INTERP_KERNEL::Exception("Invalid nature for NodeFE field : expected IntensiveMaximum !");
}

MEDCouplingFieldDouble *MEDCouplingFieldDiscretizationOnNodesFE::getMeasureField(const MEDCouplingMesh *mesh, bool isAbs) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDCouplingFieldDiscretizationOnNodesFE::getMeasureField : mesh instance specified is NULL !");
  throw INTERP_KERNEL::Exception("getMeasureField on MEDCouplingFieldDiscretizationOnNodesFE : not implemented yet !");
}

class Functor
{
private:
  const MEDCouplingGaussLocalization* _gl;
  std::size_t _nb_pts_in_cell;
  const double *_pts_in_cell;
  const double *_point;
public:
  Functor(const MEDCouplingGaussLocalization& gl, std::size_t nbPtsInCell, const double *ptsInCell, const double point[3]):_gl(&gl),_nb_pts_in_cell(nbPtsInCell),
  _pts_in_cell(ptsInCell),_point(point) { }
  std::vector<double> operator()(const std::vector<double>& x)
  {
    MEDCouplingGaussLocalization gl(_gl->getType(),_gl->getRefCoords(),x,{1.0});
    MCAuto<DataArrayDouble> shapeFunc = gl.getShapeFunctionValues();
    const double *shapeFuncPtr( shapeFunc->begin() );
    std::vector<double> ret(3,0);
    for(std::size_t iPt = 0; iPt < _nb_pts_in_cell ; ++iPt)
    {
      for(short iDim = 0 ; iDim < 3 ; ++iDim)
        ret[iDim] += shapeFuncPtr[iPt] * _pts_in_cell[3*iPt + iDim];
    }
    ret[0] -= _point[0]; ret[1] -= _point[1]; ret[2] -= _point[2];
    return ret;
  }
};

bool IsInside3D(const MEDCouplingGaussLocalization& gl, const std::vector<double>& ptsInCell, const double locInReal[3], double locInRef[3])
{
  constexpr double EPS_IN_OUT = 1e-12;
  std::size_t nbPtsInCell(ptsInCell.size()/3);
  bool ret(false);
  const double *refCoo(gl.getRefCoords().data());
  INTERP_KERNEL::NormalizedCellType ct(gl.getType());
  Functor func(gl,nbPtsInCell,ptsInCell.data(),locInReal);

  auto myJacobian = [&gl,nbPtsInCell,ptsInCell](const std::vector<double>& x, const std::vector<double>&, INTERP_KERNEL::DenseMatrix& jacobian)
  {
    MEDCouplingGaussLocalization mygl(gl.getType(),gl.getRefCoords(),x,{1.0});
    MCAuto<DataArrayDouble> shapeFunc = mygl.getDerivativeOfShapeFunctionValues();
    for(std::size_t i = 0 ; i < 3 ; ++i)
      for(std::size_t j = 0 ; j < 3 ; ++j)
      {
        double res = 0.0;
        for( std::size_t k = 0 ; k < nbPtsInCell ; ++k )
          res += ptsInCell[k*3+i] * shapeFunc->getIJ(0,3*k+j);
        jacobian[ i ][ j ] = res;
      }
  };

  // loop on refcoords as initialization point for Newton algo. vini is the initialization vector of Newton.
  for(std::size_t attemptId = 0 ; attemptId < nbPtsInCell ; ++attemptId)
  {
    std::vector<double> vini(refCoo + attemptId*3, refCoo + (attemptId+1)*3);
    try
    {
      bool check(true);
      //INTERP_KERNEL::SolveWithNewton(vini,check,func);
      INTERP_KERNEL::SolveWithNewtonWithJacobian(vini,check,func,myJacobian);
      ret = (check==false);//looks strange but OK regarding newt (SolveWithNewton) at page 387 of numerical recipes for semantic of check parameter
    }
    catch( INTERP_KERNEL::Exception& ex )
      { ret = false; }// Something get wrong during Newton process
    if(ret)
    {//Newton has converged. Now check if it converged to a point inside cell
      if( ! INTERP_KERNEL::GaussInfo::IsInOrOutForReference(ct,vini.data(),EPS_IN_OUT) )
      {// converged but locInReal has been detected outside of cell
        ret = false;
      }
    }
    if(ret)
    {
      locInRef[0] = vini[0]; locInRef[1] = vini[1]; locInRef[2] = vini[2];
      return ret;
    }
  }
  std::fill(locInRef,locInRef+3,std::numeric_limits<double>::max());
  return false;
}

void MEDCouplingFieldDiscretizationOnNodesFE::getValueOn(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *ptsCoo, double *res) const
{
  MCAuto<DataArrayDouble> res2( this->getValueOnMulti(arr,mesh,ptsCoo,1) );
  std::copy(res2->begin(),res2->end(),res);
}

void MEDCouplingFieldDiscretizationOnNodesFE::GetRefCoordOfListOf3DPtsIn3D(const MEDCouplingUMesh *umesh, const double *ptsCoo, mcIdType nbOfPts,
  std::function<void(const MEDCouplingGaussLocalization&, const std::vector<mcIdType>&)> customFunc)
{
  const double *coordsOfMesh( umesh->getCoords()->begin() );
  MEDCouplingNormalizedUnstructuredMesh<3,3> mesh_wrapper(umesh);
  BBTreeStandAlone<3,mcIdType> tree( INTERP_KERNEL::BuildBBTree( mesh_wrapper ) );
  for(mcIdType iPt = 0 ; iPt < nbOfPts ; ++iPt)
  {
    std::vector<mcIdType> elems;
    tree.getElementsAroundPoint(ptsCoo+3*iPt,elems);
    bool found(false);
    for(auto cellId = elems.cbegin() ; cellId != elems.cend() && !found ; ++cellId)
    {
      INTERP_KERNEL::NormalizedCellType gt( umesh->getTypeOfCell(*cellId) );
      std::vector<mcIdType> conn;
      umesh->getNodeIdsOfCell(*cellId,conn);
      MCAuto<DataArrayDouble> refCoo( MEDCouplingGaussLocalization::GetDefaultReferenceCoordinatesOf(gt) );
      std::vector<double> refCooCpp(refCoo->begin(),refCoo->end());
      std::vector<double> gsCoo(ptsCoo + iPt*3,ptsCoo + (iPt+1)*3); 
      MEDCouplingGaussLocalization gl(gt,refCooCpp,{0,0,0},{1.});
      std::vector<double> ptsInCell; ptsInCell.reserve(conn.size()*gl.getDimension());
      std::for_each( conn.cbegin(), conn.cend(), [coordsOfMesh,&ptsInCell](mcIdType c) { ptsInCell.insert(ptsInCell.end(),coordsOfMesh+c*3,coordsOfMesh+(c+1)*3); } );
      std::vector<double> locInRef(3);
      if( IsInside3D(gl,ptsInCell,gsCoo.data(),locInRef.data()) )
      {
        gl.setGaussCoords(locInRef);
        customFunc(gl,conn);
        found = true;
      }
    }
    if(!found)
      THROW_IK_EXCEPTION("getValueOnMulti on MEDCouplingFieldDiscretizationOnNodesFE : fail to locate point #" << iPt << " X=" << ptsCoo[0] << " Y=" << ptsCoo[1] << " Z=" << ptsCoo[2] << " !");
  }
}

const MEDCouplingUMesh *MEDCouplingFieldDiscretizationOnNodesFE::checkConfig3D(const MEDCouplingMesh *mesh) const
{
  const MEDCouplingUMesh *umesh( dynamic_cast<const MEDCouplingUMesh *>(mesh) );
  if( !umesh )
    THROW_IK_EXCEPTION("getValueOn : not implemented yet for type != MEDCouplingUMesh !");
  if(umesh->getSpaceDimension() != 3 || umesh->getMeshDimension() != 3)
    THROW_IK_EXCEPTION("getValueOn : implemented only for meshes with spacedim == 3 and meshdim == 3 !");
  umesh->checkConsistency();
  return umesh;
}

DataArrayDouble *MEDCouplingFieldDiscretizationOnNodesFE::getValueOnMulti(const DataArrayDouble *arr, const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfTargetPoints) const
{
  if(!arr || !arr->isAllocated())
    throw INTERP_KERNEL::Exception("getValueOnMulti : input array is null or not allocated !");
  mcIdType nbOfRows=getNumberOfMeshPlaces(mesh);
  if(arr->getNumberOfTuples()!=nbOfRows)
  {
    THROW_IK_EXCEPTION( "getValueOnMulti : input array does not have correct number of tuples ! Excepted " << nbOfRows << " having " << arr->getNumberOfTuples() << " !")
  }
  const MEDCouplingUMesh *umesh = checkConfig3D(mesh);
  std::size_t nbCompo( arr->getNumberOfComponents() );
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbOfTargetPoints,nbCompo);
  double *res( ret->getPointer() );

  auto arrayFeeder = [arr,&res,nbCompo](const MEDCouplingGaussLocalization& gl, const std::vector<mcIdType>& conn)
  {
    MCAuto<DataArrayDouble> resVector( gl.getShapeFunctionValues() );
    {
      std::for_each(res,res+nbCompo,[](double& v) { v = 0.0; });
      for(std::size_t iComp = 0 ; iComp < nbCompo ; ++iComp)
        for(int iPt = 0 ; iPt < gl.getNumberOfPtsInRefCell(); ++iPt)
        {
          {
            res[iComp] += resVector->getIJ(0,iPt) * arr->getIJ(conn[iPt],iComp);
          }
        }
      res += nbCompo;
    }
  };

  GetRefCoordOfListOf3DPtsIn3D(umesh,loc,nbOfTargetPoints,arrayFeeder);
  return ret.retn();
}

/*!
 * Returns for each \a nbOfPoints point in \a loc its coordinate in reference element.
 */
MCAuto<DataArrayDouble> MEDCouplingFieldDiscretizationOnNodesFE::getCooInRefElement(const MEDCouplingMesh *mesh, const double *loc, mcIdType nbOfPoints) const
{
  const MEDCouplingUMesh *umesh = checkConfig3D(mesh);
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New());
  ret->alloc(nbOfPoints,3);
  double *retPtr(ret->getPointer() );
  
  auto arrayFeeder = [&retPtr](const MEDCouplingGaussLocalization& gl, const std::vector<mcIdType>& conn)
  {
    std::vector<double> resVector( gl.getGaussCoords() );
    {
      std::copy(resVector.begin(),resVector.end(),retPtr);
      retPtr += 3;
    }
  };

  GetRefCoordOfListOf3DPtsIn3D(umesh,loc,nbOfPoints,arrayFeeder);
  return ret;
}
