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

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingCMesh.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingSkyLineArray.hxx"
#include "CellModel.hxx"
#include "VolSurfUser.txx"
#include "InterpolationUtils.hxx"
#include "PointLocatorAlgos.txx"
#include "BBTree.txx"
#include "BBTreeDst.txx"
#include "SplitterTetra.hxx"
#include "DiameterCalculator.hxx"
#include "DirectedBoundingBox.hxx"
#include "InterpKernelMatrixTools.hxx"
#include "InterpKernelMeshQuality.hxx"
#include "InterpKernelCellSimplify.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelAutoPtr.hxx"
#include "InterpKernelGeo2DNode.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "MEDCouplingUMesh_internal.hxx"

#include <sstream>
#include <fstream>
#include <numeric>
#include <cstring>
#include <limits>
#include <list>

using namespace MEDCoupling;

/*!
 * This method checks that all arrays are set. If yes nothing done if no an exception is thrown.
 */
void MEDCouplingUMesh::checkFullyDefined() const
{
  if(!_nodal_connec_index || !_nodal_connec || !_coords)
    throw INTERP_KERNEL::Exception("Reverse nodal connectivity computation requires full connectivity and coordinates set in unstructured mesh.");
}

/*!
 * This method checks that all connectivity arrays are set. If yes nothing done if no an exception is thrown.
 */
void MEDCouplingUMesh::checkConnectivityFullyDefined() const
{
  if(!_nodal_connec_index || !_nodal_connec)
    throw INTERP_KERNEL::Exception("Reverse nodal connectivity computation requires full connectivity set in unstructured mesh.");
}

void MEDCouplingUMesh::reprConnectivityOfThisLL(std::ostringstream& stream) const
{
  if(_nodal_connec!=0 && _nodal_connec_index!=0)
    {
      int nbOfCells=getNumberOfCells();
      const int *c=_nodal_connec->getConstPointer();
      const int *ci=_nodal_connec_index->getConstPointer();
      for(int i=0;i<nbOfCells;i++)
        {
          const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)c[ci[i]]);
          stream << "Cell #" << i << " " << cm.getRepr() << " : ";
          std::copy(c+ci[i]+1,c+ci[i+1],std::ostream_iterator<int>(stream," "));
          stream << "\n";
        }
    }
  else
    stream << "Connectivity not defined !\n";
}


/*!
 * This method implements policy 0 of virtual method MEDCoupling::MEDCouplingUMesh::simplexize.
 */
DataArrayInt *MEDCouplingUMesh::simplexizePol0()
{
  checkConnectivityFullyDefined();
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplexizePol0 : this policy is only available for mesh with meshdim == 2 !");
  int nbOfCells=getNumberOfCells();
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  int nbOfCutCells=getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD4);
  ret->alloc(nbOfCells+nbOfCutCells,1);
  if(nbOfCutCells==0) { ret->iota(0); return ret.retn(); }
  int *retPt=ret->getPointer();
  MCAuto<DataArrayInt> newConn=DataArrayInt::New();
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New();
  newConnI->alloc(nbOfCells+nbOfCutCells+1,1);
  newConn->alloc(getNodalConnectivityArrayLen()+3*nbOfCutCells,1);
  int *pt=newConn->getPointer();
  int *ptI=newConnI->getPointer();
  ptI[0]=0;
  const int *oldc=_nodal_connec->begin();
  const int *ci=_nodal_connec_index->begin();
  for(int i=0;i<nbOfCells;i++,ci++)
    {
      if((INTERP_KERNEL::NormalizedCellType)oldc[ci[0]]==INTERP_KERNEL::NORM_QUAD4)
        {
          const int tmp[8]={(int)INTERP_KERNEL::NORM_TRI3,oldc[ci[0]+1],oldc[ci[0]+2],oldc[ci[0]+3],
            (int)INTERP_KERNEL::NORM_TRI3,oldc[ci[0]+1],oldc[ci[0]+3],oldc[ci[0]+4]};
          pt=std::copy(tmp,tmp+8,pt);
          ptI[1]=ptI[0]+4;
          ptI[2]=ptI[0]+8;
          *retPt++=i;
          *retPt++=i;
          ptI+=2;
        }
      else
        {
          pt=std::copy(oldc+ci[0],oldc+ci[1],pt);
          ptI[1]=ptI[0]+ci[1]-ci[0];
          ptI++;
          *retPt++=i;
        }
    }
  _nodal_connec->decrRef();
  _nodal_connec=newConn.retn();
  _nodal_connec_index->decrRef();
  _nodal_connec_index=newConnI.retn();
  computeTypes();
  updateTime();
  return ret.retn();
}

/*!
 * This method implements policy 1 of virtual method MEDCoupling::MEDCouplingUMesh::simplexize.
 */
DataArrayInt *MEDCouplingUMesh::simplexizePol1()
{
  checkConnectivityFullyDefined();
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplexizePol0 : this policy is only available for mesh with meshdim == 2 !");
  int nbOfCells=getNumberOfCells();
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  int nbOfCutCells=getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD4);
  ret->alloc(nbOfCells+nbOfCutCells,1);
  if(nbOfCutCells==0) { ret->iota(0); return ret.retn(); }
  int *retPt=ret->getPointer();
  MCAuto<DataArrayInt> newConn=DataArrayInt::New();
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New();
  newConnI->alloc(nbOfCells+nbOfCutCells+1,1);
  newConn->alloc(getNodalConnectivityArrayLen()+3*nbOfCutCells,1);
  int *pt=newConn->getPointer();
  int *ptI=newConnI->getPointer();
  ptI[0]=0;
  const int *oldc=_nodal_connec->begin();
  const int *ci=_nodal_connec_index->begin();
  for(int i=0;i<nbOfCells;i++,ci++)
    {
      if((INTERP_KERNEL::NormalizedCellType)oldc[ci[0]]==INTERP_KERNEL::NORM_QUAD4)
        {
          const int tmp[8]={(int)INTERP_KERNEL::NORM_TRI3,oldc[ci[0]+1],oldc[ci[0]+2],oldc[ci[0]+4],
            (int)INTERP_KERNEL::NORM_TRI3,oldc[ci[0]+2],oldc[ci[0]+3],oldc[ci[0]+4]};
          pt=std::copy(tmp,tmp+8,pt);
          ptI[1]=ptI[0]+4;
          ptI[2]=ptI[0]+8;
          *retPt++=i;
          *retPt++=i;
          ptI+=2;
        }
      else
        {
          pt=std::copy(oldc+ci[0],oldc+ci[1],pt);
          ptI[1]=ptI[0]+ci[1]-ci[0];
          ptI++;
          *retPt++=i;
        }
    }
  _nodal_connec->decrRef();
  _nodal_connec=newConn.retn();
  _nodal_connec_index->decrRef();
  _nodal_connec_index=newConnI.retn();
  computeTypes();
  updateTime();
  return ret.retn();
}

/*!
 * This method implements policy INTERP_KERNEL::PLANAR_FACE_5 of virtual method MEDCoupling::MEDCouplingUMesh::simplexize.
 */
DataArrayInt *MEDCouplingUMesh::simplexizePlanarFace5()
{
  checkConnectivityFullyDefined();
  if(getMeshDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplexizePlanarFace5 : this policy is only available for mesh with meshdim == 3 !");
  int nbOfCells=getNumberOfCells();
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  int nbOfCutCells=getNumberOfCellsWithType(INTERP_KERNEL::NORM_HEXA8);
  ret->alloc(nbOfCells+4*nbOfCutCells,1);
  if(nbOfCutCells==0) { ret->iota(0); return ret.retn(); }
  int *retPt=ret->getPointer();
  MCAuto<DataArrayInt> newConn=DataArrayInt::New();
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New();
  newConnI->alloc(nbOfCells+4*nbOfCutCells+1,1);
  newConn->alloc(getNodalConnectivityArrayLen()+16*nbOfCutCells,1);//21
  int *pt=newConn->getPointer();
  int *ptI=newConnI->getPointer();
  ptI[0]=0;
  const int *oldc=_nodal_connec->begin();
  const int *ci=_nodal_connec_index->begin();
  for(int i=0;i<nbOfCells;i++,ci++)
    {
      if((INTERP_KERNEL::NormalizedCellType)oldc[ci[0]]==INTERP_KERNEL::NORM_HEXA8)
        {
          for(int j=0;j<5;j++,pt+=5,ptI++)
            {
              pt[0]=(int)INTERP_KERNEL::NORM_TETRA4;
              pt[1]=oldc[ci[0]+INTERP_KERNEL::SPLIT_NODES_5_WO[4*j+0]+1]; pt[2]=oldc[ci[0]+INTERP_KERNEL::SPLIT_NODES_5_WO[4*j+1]+1]; pt[3]=oldc[ci[0]+INTERP_KERNEL::SPLIT_NODES_5_WO[4*j+2]+1]; pt[4]=oldc[ci[0]+INTERP_KERNEL::SPLIT_NODES_5_WO[4*j+3]+1];
              *retPt++=i;
              ptI[1]=ptI[0]+5;
            }
        }
      else
        {
          pt=std::copy(oldc+ci[0],oldc+ci[1],pt);
          ptI[1]=ptI[0]+ci[1]-ci[0];
          ptI++;
          *retPt++=i;
        }
    }
  _nodal_connec->decrRef();
  _nodal_connec=newConn.retn();
  _nodal_connec_index->decrRef();
  _nodal_connec_index=newConnI.retn();
  computeTypes();
  updateTime();
  return ret.retn();
}

/*!
 * This method implements policy INTERP_KERNEL::PLANAR_FACE_6 of virtual method MEDCoupling::MEDCouplingUMesh::simplexize.
 */
DataArrayInt *MEDCouplingUMesh::simplexizePlanarFace6()
{
  checkConnectivityFullyDefined();
  if(getMeshDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::simplexizePlanarFace6 : this policy is only available for mesh with meshdim == 3 !");
  int nbOfCells=getNumberOfCells();
  MCAuto<DataArrayInt> ret=DataArrayInt::New();
  int nbOfCutCells=getNumberOfCellsWithType(INTERP_KERNEL::NORM_HEXA8);
  ret->alloc(nbOfCells+5*nbOfCutCells,1);
  if(nbOfCutCells==0) { ret->iota(0); return ret.retn(); }
  int *retPt=ret->getPointer();
  MCAuto<DataArrayInt> newConn=DataArrayInt::New();
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New();
  newConnI->alloc(nbOfCells+5*nbOfCutCells+1,1);
  newConn->alloc(getNodalConnectivityArrayLen()+21*nbOfCutCells,1);
  int *pt=newConn->getPointer();
  int *ptI=newConnI->getPointer();
  ptI[0]=0;
  const int *oldc=_nodal_connec->begin();
  const int *ci=_nodal_connec_index->begin();
  for(int i=0;i<nbOfCells;i++,ci++)
    {
      if((INTERP_KERNEL::NormalizedCellType)oldc[ci[0]]==INTERP_KERNEL::NORM_HEXA8)
        {
          for(int j=0;j<6;j++,pt+=5,ptI++)
            {
              pt[0]=(int)INTERP_KERNEL::NORM_TETRA4;
              pt[1]=oldc[ci[0]+INTERP_KERNEL::SPLIT_NODES_6_WO[4*j+0]+1]; pt[2]=oldc[ci[0]+INTERP_KERNEL::SPLIT_NODES_6_WO[4*j+1]+1]; pt[3]=oldc[ci[0]+INTERP_KERNEL::SPLIT_NODES_6_WO[4*j+2]+1]; pt[4]=oldc[ci[0]+INTERP_KERNEL::SPLIT_NODES_6_WO[4*j+3]+1];
              *retPt++=i;
              ptI[1]=ptI[0]+5;
            }
        }
      else
        {
          pt=std::copy(oldc+ci[0],oldc+ci[1],pt);
          ptI[1]=ptI[0]+ci[1]-ci[0];
          ptI++;
          *retPt++=i;
        }
    }
  _nodal_connec->decrRef();
  _nodal_connec=newConn.retn();
  _nodal_connec_index->decrRef();
  _nodal_connec_index=newConnI.retn();
  computeTypes();
  updateTime();
  return ret.retn();
}

/*!
 * Tessellates \a this 2D mesh by dividing not straight edges of quadratic faces,
 * so that the number of cells remains the same. Quadratic faces are converted to
 * polygons. This method works only for 2D meshes in
 * 2D space. If no cells are quadratic (INTERP_KERNEL::NORM_QUAD8,
 * INTERP_KERNEL::NORM_TRI6, INTERP_KERNEL::NORM_QPOLYG ), \a this mesh remains unchanged.
 * \warning This method can lead to a huge amount of nodes if \a eps is very low.
 *  \param [in] eps - specifies the maximal angle (in radians) between 2 sub-edges of
 *         a polylinized edge constituting the input polygon.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If \a this->getMeshDimension() != 2.
 *  \throw If \a this->getSpaceDimension() != 2.
 */
void MEDCouplingUMesh::tessellate2DInternal(double eps)
{
  checkFullyDefined();
  if(getMeshDimension()!=2 || getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2DInternal works on umeshes with meshdim equal to 2 and spaceDim equal to 2 too!");
  double epsa=fabs(eps);
  if(epsa<std::numeric_limits<double>::min())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2DInternal : epsilon is null ! Please specify a higher epsilon. If too tiny it can lead to a huge amount of nodes and memory !");
  MCAuto<DataArrayInt> desc1(DataArrayInt::New()),descIndx1(DataArrayInt::New()),revDesc1(DataArrayInt::New()),revDescIndx1(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> mDesc(buildDescendingConnectivity2(desc1,descIndx1,revDesc1,revDescIndx1));
  revDesc1=0; revDescIndx1=0;
  mDesc->tessellate2D(eps);
  subDivide2DMesh(mDesc->_nodal_connec->begin(),mDesc->_nodal_connec_index->begin(),desc1->begin(),descIndx1->begin());
  setCoords(mDesc->getCoords());
}

/*!
 * Tessellates \a this 1D mesh in 2D space by dividing not straight quadratic edges.
 * \warning This method can lead to a huge amount of nodes if \a eps is very low.
 *  \param [in] eps - specifies the maximal angle (in radian) between 2 sub-edges of
 *         a sub-divided edge.
 *  \throw If the coordinates array is not set.
 *  \throw If the nodal connectivity of cells is not defined.
 *  \throw If \a this->getMeshDimension() != 1.
 *  \throw If \a this->getSpaceDimension() != 2.
 */
void MEDCouplingUMesh::tessellate2DCurveInternal(double eps)
{
  checkFullyDefined();
  if(getMeshDimension()!=1 || getSpaceDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2DCurveInternal works on umeshes with meshdim equal to 1 and spaceDim equal to 2 too!");
  double epsa=fabs(eps);
  if(epsa<std::numeric_limits<double>::min())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::tessellate2DCurveInternal : epsilon is null ! Please specify a higher epsilon. If too tiny it can lead to a huge amount of nodes and memory !");
  INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision arcPrec(1.e-10);  // RAII
  int nbCells=getNumberOfCells();
  int nbNodes=getNumberOfNodes();
  const int *conn=_nodal_connec->begin();
  const int *connI=_nodal_connec_index->begin();
  const double *coords=_coords->begin();
  std::vector<double> addCoo;
  std::vector<int> newConn;//no direct DataArrayInt because interface with Geometric2D
  MCAuto<DataArrayInt> newConnI(DataArrayInt::New());
  newConnI->alloc(nbCells+1,1);
  int *newConnIPtr=newConnI->getPointer();
  *newConnIPtr=0;
  int tmp1[3];
  INTERP_KERNEL::Node *tmp2[3];
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(int i=0;i<nbCells;i++,newConnIPtr++)
    {
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
      if(cm.isQuadratic())
        {//assert(connI[i+1]-connI[i]-1==3)
          tmp1[0]=conn[connI[i]+1+0]; tmp1[1]=conn[connI[i]+1+1]; tmp1[2]=conn[connI[i]+1+2];
          tmp2[0]=new INTERP_KERNEL::Node(coords[2*tmp1[0]],coords[2*tmp1[0]+1]);
          tmp2[1]=new INTERP_KERNEL::Node(coords[2*tmp1[1]],coords[2*tmp1[1]+1]);
          tmp2[2]=new INTERP_KERNEL::Node(coords[2*tmp1[2]],coords[2*tmp1[2]+1]);
          INTERP_KERNEL::EdgeArcCircle *eac=INTERP_KERNEL::EdgeArcCircle::BuildFromNodes(tmp2[0],tmp2[2],tmp2[1]);
          if(eac)
            {
              eac->tesselate(tmp1,nbNodes,epsa,newConn,addCoo);
              types.insert((INTERP_KERNEL::NormalizedCellType)newConn[newConnIPtr[0]]);
              delete eac;
              newConnIPtr[1]=(int)newConn.size();
            }
          else
            {
              types.insert(INTERP_KERNEL::NORM_SEG2);
              newConn.push_back(INTERP_KERNEL::NORM_SEG2);
              newConn.insert(newConn.end(),conn+connI[i]+1,conn+connI[i]+3);
              newConnIPtr[1]=newConnIPtr[0]+3;
            }
        }
      else
        {
          types.insert((INTERP_KERNEL::NormalizedCellType)conn[connI[i]]);
          newConn.insert(newConn.end(),conn+connI[i],conn+connI[i+1]);
          newConnIPtr[1]=newConnIPtr[0]+3;
        }
    }
  if(addCoo.empty() && newConn.size()==_nodal_connec->getNumberOfTuples())//nothing happens during tessellation : no update needed
    return ;
  _types=types;
  DataArrayInt::SetArrayIn(newConnI,_nodal_connec_index);
  MCAuto<DataArrayInt> newConnArr=DataArrayInt::New();
  newConnArr->alloc((int)newConn.size(),1);
  std::copy(newConn.begin(),newConn.end(),newConnArr->getPointer());
  DataArrayInt::SetArrayIn(newConnArr,_nodal_connec);
  MCAuto<DataArrayDouble> newCoords=DataArrayDouble::New();
  newCoords->alloc(nbNodes+((int)addCoo.size())/2,2);
  double *work=std::copy(_coords->begin(),_coords->end(),newCoords->getPointer());
  std::copy(addCoo.begin(),addCoo.end(),work);
  DataArrayDouble::SetArrayIn(newCoords,_coords);
  updateTime();
}


/*!
 * This private method is used to subdivide edges of a mesh with meshdim==2. If \a this has no a meshdim equal to 2 an exception will be thrown.
 * This method completly ignore coordinates.
 * \param nodeSubdived is the nodal connectivity of subdivision of edges
 * \param nodeIndxSubdived is the nodal connectivity index of subdivision of edges
 * \param desc is descending connectivity in format specified in MEDCouplingUMesh::buildDescendingConnectivity2
 * \param descIndex is descending connectivity index in format specified in MEDCouplingUMesh::buildDescendingConnectivity2
 */
void MEDCouplingUMesh::subDivide2DMesh(const int *nodeSubdived, const int *nodeIndxSubdived, const int *desc, const int *descIndex)
{
  checkFullyDefined();
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::subDivide2DMesh : works only on umesh with meshdim==2 !");
  int nbOfCells=getNumberOfCells();
  int *connI=_nodal_connec_index->getPointer();
  int newConnLgth=0;
  for(int i=0;i<nbOfCells;i++,connI++)
    {
      int offset=descIndex[i];
      int nbOfEdges=descIndex[i+1]-offset;
      //
      bool ddirect=desc[offset+nbOfEdges-1]>0;
      int eedgeId=std::abs(desc[offset+nbOfEdges-1])-1;
      int ref=ddirect?nodeSubdived[nodeIndxSubdived[eedgeId+1]-1]:nodeSubdived[nodeIndxSubdived[eedgeId]+1];
      for(int j=0;j<nbOfEdges;j++)
        {
          bool direct=desc[offset+j]>0;
          int edgeId=std::abs(desc[offset+j])-1;
          if(!INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)nodeSubdived[nodeIndxSubdived[edgeId]]).isQuadratic())
            {
              int id1=nodeSubdived[nodeIndxSubdived[edgeId]+1];
              int id2=nodeSubdived[nodeIndxSubdived[edgeId+1]-1];
              int ref2=direct?id1:id2;
              if(ref==ref2)
                {
                  int nbOfSubNodes=nodeIndxSubdived[edgeId+1]-nodeIndxSubdived[edgeId]-1;
                  newConnLgth+=nbOfSubNodes-1;
                  ref=direct?id2:id1;
                }
              else
                {
                  std::ostringstream oss; oss << "MEDCouplingUMesh::subDivide2DMesh : On polygon #" << i << " edgeid #" << j << " subedges mismatch : end subedge k!=start subedge k+1 !";
                  throw INTERP_KERNEL::Exception(oss.str());
                }
            }
          else
            {
              throw INTERP_KERNEL::Exception("MEDCouplingUMesh::subDivide2DMesh : this method only subdivides into linear edges !");
            }
        }
      newConnLgth++;//+1 is for cell type
      connI[1]=newConnLgth;
    }
  //
  MCAuto<DataArrayInt> newConn=DataArrayInt::New();
  newConn->alloc(newConnLgth,1);
  int *work=newConn->getPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      *work++=INTERP_KERNEL::NORM_POLYGON;
      int offset=descIndex[i];
      int nbOfEdges=descIndex[i+1]-offset;
      for(int j=0;j<nbOfEdges;j++)
        {
          bool direct=desc[offset+j]>0;
          int edgeId=std::abs(desc[offset+j])-1;
          if(direct)
            work=std::copy(nodeSubdived+nodeIndxSubdived[edgeId]+1,nodeSubdived+nodeIndxSubdived[edgeId+1]-1,work);
          else
            {
              int nbOfSubNodes=nodeIndxSubdived[edgeId+1]-nodeIndxSubdived[edgeId]-1;
              std::reverse_iterator<const int *> it(nodeSubdived+nodeIndxSubdived[edgeId+1]);
              work=std::copy(it,it+nbOfSubNodes-1,work);
            }
        }
    }
  DataArrayInt::SetArrayIn(newConn,_nodal_connec);
  _types.clear();
  if(nbOfCells>0)
    _types.insert(INTERP_KERNEL::NORM_POLYGON);
}

/*!
 * Keeps from \a this only cells which constituing point id are in the ids specified by [ \a begin,\a end ).
 * The resulting cell ids are stored at the end of the 'cellIdsKept' parameter.
 * Parameter \a fullyIn specifies if a cell that has part of its nodes in ids array is kept or not.
 * If \a fullyIn is true only cells whose ids are \b fully contained in [ \a begin,\a end ) tab will be kept.
 *
 * \param [in] begin input start of array of node ids.
 * \param [in] end input end of array of node ids.
 * \param [in] fullyIn input that specifies if all node ids must be in [ \a begin,\a end ) array to consider cell to be in.
 * \param [in,out] cellIdsKeptArr array where all candidate cell ids are put at the end.
 */
void MEDCouplingUMesh::fillCellIdsToKeepFromNodeIds(const int *begin, const int *end, bool fullyIn, DataArrayInt *&cellIdsKeptArr) const
{
  MCAuto<DataArrayInt> cellIdsKept=DataArrayInt::New(); cellIdsKept->alloc(0,1);
  checkConnectivityFullyDefined();
  int tmp=-1;
  int sz=getNodalConnectivity()->getMaxValue(tmp); sz=std::max(sz,0)+1;
  std::vector<bool> fastFinder(sz,false);
  for(const int *work=begin;work!=end;work++)
    if(*work>=0 && *work<sz)
      fastFinder[*work]=true;
  int nbOfCells=getNumberOfCells();
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *connIndex=getNodalConnectivityIndex()->getConstPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      int ref=0,nbOfHit=0;
      for(const int *work2=conn+connIndex[i]+1;work2!=conn+connIndex[i+1];work2++)
        if(*work2>=0)
          {
            ref++;
            if(fastFinder[*work2])
              nbOfHit++;
          }
      if((ref==nbOfHit && fullyIn) || (nbOfHit!=0 && !fullyIn))
        cellIdsKept->pushBackSilent(i);
    }
  cellIdsKeptArr=cellIdsKept.retn();
}

/*!
 * This method works on a 3D curve linear mesh that is to say (meshDim==1 and spaceDim==3).
 * If it is not the case an exception will be thrown.
 * This method is non const because the coordinate of \a this can be appended with some new points issued from
 * intersection of plane defined by ('origin','vec').
 * This method has one in/out parameter : 'cut3DCurve'.
 * Param 'cut3DCurve' is expected to be of size 'this->getNumberOfCells()'. For each i in [0,'this->getNumberOfCells()')
 * if cut3DCurve[i]==-2, it means that for cell #i in \a this nothing has been detected previously.
 * if cut3DCurve[i]==-1, it means that cell#i has been already detected to be fully part of plane defined by ('origin','vec').
 * This method will throw an exception if \a this contains a non linear segment.
 */
void MEDCouplingUMesh::split3DCurveWithPlane(const double *origin, const double *vec, double eps, std::vector<int>& cut3DCurve)
{
  checkFullyDefined();
  if(getMeshDimension()!=1 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split3DCurveWithPlane works on umeshes with meshdim equal to 1 and spaceDim equal to 3 !");
  int ncells=getNumberOfCells();
  int nnodes=getNumberOfNodes();
  double vec2[3],vec3[3],vec4[3];
  double normm=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  if(normm<1e-6)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split3DCurveWithPlane : parameter 'vec' should have a norm2 greater than 1e-6 !");
  vec2[0]=vec[0]/normm; vec2[1]=vec[1]/normm; vec2[2]=vec[2]/normm;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connI=_nodal_connec_index->getConstPointer();
  const double *coo=_coords->getConstPointer();
  std::vector<double> addCoo;
  for(int i=0;i<ncells;i++)
    {
      if(conn[connI[i]]==(int)INTERP_KERNEL::NORM_SEG2)
        {
          if(cut3DCurve[i]==-2)
            {
              int st=conn[connI[i]+1],endd=conn[connI[i]+2];
              vec3[0]=coo[3*endd]-coo[3*st]; vec3[1]=coo[3*endd+1]-coo[3*st+1]; vec3[2]=coo[3*endd+2]-coo[3*st+2];
              double normm2=sqrt(vec3[0]*vec3[0]+vec3[1]*vec3[1]+vec3[2]*vec3[2]);
              double colin=std::abs((vec3[0]*vec2[0]+vec3[1]*vec2[1]+vec3[2]*vec2[2])/normm2);
              if(colin>eps)//if colin<=eps -> current SEG2 is colinear to the input plane
                {
                  const double *st2=coo+3*st;
                  vec4[0]=st2[0]-origin[0]; vec4[1]=st2[1]-origin[1]; vec4[2]=st2[2]-origin[2];
                  double pos=-(vec4[0]*vec2[0]+vec4[1]*vec2[1]+vec4[2]*vec2[2])/((vec3[0]*vec2[0]+vec3[1]*vec2[1]+vec3[2]*vec2[2]));
                  if(pos>eps && pos<1-eps)
                    {
                      int nNode=((int)addCoo.size())/3;
                      vec4[0]=st2[0]+pos*vec3[0]; vec4[1]=st2[1]+pos*vec3[1]; vec4[2]=st2[2]+pos*vec3[2];
                      addCoo.insert(addCoo.end(),vec4,vec4+3);
                      cut3DCurve[i]=nnodes+nNode;
                    }
                }
            }
        }
      else
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::split3DCurveWithPlane : this method is only available for linear cell (NORM_SEG2) !");
    }
  if(!addCoo.empty())
    {
      int newNbOfNodes=nnodes+((int)addCoo.size())/3;
      MCAuto<DataArrayDouble> coo2=DataArrayDouble::New();
      coo2->alloc(newNbOfNodes,3);
      double *tmp=coo2->getPointer();
      tmp=std::copy(_coords->begin(),_coords->end(),tmp);
      std::copy(addCoo.begin(),addCoo.end(),tmp);
      DataArrayDouble::SetArrayIn(coo2,_coords);
    }
}

/*!
 * This method incarnates the policy 0 for MEDCouplingUMesh::buildExtrudedMesh method.
 * \param mesh1D is the input 1D mesh used for translation computation.
 * \return newCoords new coords filled by this method.
 */
DataArrayDouble *MEDCouplingUMesh::fillExtCoordsUsingTranslation(const MEDCouplingUMesh *mesh1D, bool isQuad) const
{
  int oldNbOfNodes=getNumberOfNodes();
  int nbOf1DCells=mesh1D->getNumberOfCells();
  int spaceDim=getSpaceDimension();
  DataArrayDouble *ret=DataArrayDouble::New();
  std::vector<bool> isQuads;
  int nbOfLevsInVec=isQuad?2*nbOf1DCells+1:nbOf1DCells+1;
  ret->alloc(oldNbOfNodes*nbOfLevsInVec,spaceDim);
  double *retPtr=ret->getPointer();
  const double *coords=getCoords()->getConstPointer();
  double *work=std::copy(coords,coords+spaceDim*oldNbOfNodes,retPtr);
  std::vector<int> v;
  std::vector<double> c;
  double vec[3];
  v.reserve(3);
  c.reserve(6);
  for(int i=0;i<nbOf1DCells;i++)
    {
      v.resize(0);
      mesh1D->getNodeIdsOfCell(i,v);
      c.resize(0);
      mesh1D->getCoordinatesOfNode(v[isQuad?2:1],c);
      mesh1D->getCoordinatesOfNode(v[0],c);
      std::transform(c.begin(),c.begin()+spaceDim,c.begin()+spaceDim,vec,std::minus<double>());
      for(int j=0;j<oldNbOfNodes;j++)
        work=std::transform(vec,vec+spaceDim,retPtr+spaceDim*(i*oldNbOfNodes+j),work,std::plus<double>());
      if(isQuad)
        {
          c.resize(0);
          mesh1D->getCoordinatesOfNode(v[1],c);
          mesh1D->getCoordinatesOfNode(v[0],c);
          std::transform(c.begin(),c.begin()+spaceDim,c.begin()+spaceDim,vec,std::minus<double>());
          for(int j=0;j<oldNbOfNodes;j++)
            work=std::transform(vec,vec+spaceDim,retPtr+spaceDim*(i*oldNbOfNodes+j),work,std::plus<double>());
        }
    }
  ret->copyStringInfoFrom(*getCoords());
  return ret;
}

/*!
 * This method incarnates the policy 1 for MEDCouplingUMesh::buildExtrudedMesh method.
 * \param mesh1D is the input 1D mesh used for translation and automatic rotation computation.
 * \return newCoords new coords filled by this method.
 */
DataArrayDouble *MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation(const MEDCouplingUMesh *mesh1D, bool isQuad) const
{
  if(mesh1D->getSpaceDimension()==2)
    return fillExtCoordsUsingTranslAndAutoRotation2D(mesh1D,isQuad);
  if(mesh1D->getSpaceDimension()==3)
    return fillExtCoordsUsingTranslAndAutoRotation3D(mesh1D,isQuad);
  throw INTERP_KERNEL::Exception("Not implemented rotation and translation alg. for spacedim other than 2 and 3 !");
}

/*!
 * This method incarnates the policy 1 for MEDCouplingUMesh::buildExtrudedMesh method.
 * \param mesh1D is the input 1D mesh used for translation and automatic rotation computation.
 * \return newCoords new coords filled by this method.
 */
DataArrayDouble *MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation2D(const MEDCouplingUMesh *mesh1D, bool isQuad) const
{
  if(isQuad)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation2D : not implemented for quadratic cells !");
  int oldNbOfNodes=getNumberOfNodes();
  int nbOf1DCells=mesh1D->getNumberOfCells();
  if(nbOf1DCells<2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation2D : impossible to detect any angle of rotation ! Change extrusion policy 1->0 !");
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  int nbOfLevsInVec=nbOf1DCells+1;
  ret->alloc(oldNbOfNodes*nbOfLevsInVec,2);
  double *retPtr=ret->getPointer();
  retPtr=std::copy(getCoords()->getConstPointer(),getCoords()->getConstPointer()+getCoords()->getNbOfElems(),retPtr);
  MCAuto<MEDCouplingUMesh> tmp=MEDCouplingUMesh::New();
  MCAuto<DataArrayDouble> tmp2=getCoords()->deepCopy();
  tmp->setCoords(tmp2);
  const double *coo1D=mesh1D->getCoords()->getConstPointer();
  const int *conn1D=mesh1D->getNodalConnectivity()->getConstPointer();
  const int *connI1D=mesh1D->getNodalConnectivityIndex()->getConstPointer();
  for(int i=1;i<nbOfLevsInVec;i++)
    {
      const double *begin=coo1D+2*conn1D[connI1D[i-1]+1];
      const double *end=coo1D+2*conn1D[connI1D[i-1]+2];
      const double *third=i+1<nbOfLevsInVec?coo1D+2*conn1D[connI1D[i]+2]:coo1D+2*conn1D[connI1D[i-2]+1];
      const double vec[2]={end[0]-begin[0],end[1]-begin[1]};
      tmp->translate(vec);
      double tmp3[2],radius,alpha,alpha0;
      const double *p0=i+1<nbOfLevsInVec?begin:third;
      const double *p1=i+1<nbOfLevsInVec?end:begin;
      const double *p2=i+1<nbOfLevsInVec?third:end;
      INTERP_KERNEL::EdgeArcCircle::GetArcOfCirclePassingThru(p0,p1,p2,tmp3,radius,alpha,alpha0);
      double cosangle=i+1<nbOfLevsInVec?(p0[0]-tmp3[0])*(p1[0]-tmp3[0])+(p0[1]-tmp3[1])*(p1[1]-tmp3[1]):(p2[0]-tmp3[0])*(p1[0]-tmp3[0])+(p2[1]-tmp3[1])*(p1[1]-tmp3[1]);
      double angle=acos(cosangle/(radius*radius));
      tmp->rotate(end,0,angle);
      retPtr=std::copy(tmp2->getConstPointer(),tmp2->getConstPointer()+tmp2->getNbOfElems(),retPtr);
    }
  return ret.retn();
}

/*!
 * This method incarnates the policy 1 for MEDCouplingUMesh::buildExtrudedMesh method.
 * \param mesh1D is the input 1D mesh used for translation and automatic rotation computation.
 * \return newCoords new coords filled by this method.
 */
DataArrayDouble *MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation3D(const MEDCouplingUMesh *mesh1D, bool isQuad) const
{
  if(isQuad)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation3D : not implemented for quadratic cells !");
  int oldNbOfNodes=getNumberOfNodes();
  int nbOf1DCells=mesh1D->getNumberOfCells();
  if(nbOf1DCells<2)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::fillExtCoordsUsingTranslAndAutoRotation3D : impossible to detect any angle of rotation ! Change extrusion policy 1->0 !");
  MCAuto<DataArrayDouble> ret=DataArrayDouble::New();
  int nbOfLevsInVec=nbOf1DCells+1;
  ret->alloc(oldNbOfNodes*nbOfLevsInVec,3);
  double *retPtr=ret->getPointer();
  retPtr=std::copy(getCoords()->getConstPointer(),getCoords()->getConstPointer()+getCoords()->getNbOfElems(),retPtr);
  MCAuto<MEDCouplingUMesh> tmp=MEDCouplingUMesh::New();
  MCAuto<DataArrayDouble> tmp2=getCoords()->deepCopy();
  tmp->setCoords(tmp2);
  const double *coo1D=mesh1D->getCoords()->getConstPointer();
  const int *conn1D=mesh1D->getNodalConnectivity()->getConstPointer();
  const int *connI1D=mesh1D->getNodalConnectivityIndex()->getConstPointer();
  for(int i=1;i<nbOfLevsInVec;i++)
    {
      const double *begin=coo1D+3*conn1D[connI1D[i-1]+1];
      const double *end=coo1D+3*conn1D[connI1D[i-1]+2];
      const double *third=i+1<nbOfLevsInVec?coo1D+3*conn1D[connI1D[i]+2]:coo1D+3*conn1D[connI1D[i-2]+1];
      const double vec[3]={end[0]-begin[0],end[1]-begin[1],end[2]-begin[2]};
      tmp->translate(vec);
      double tmp3[2],radius,alpha,alpha0;
      const double *p0=i+1<nbOfLevsInVec?begin:third;
      const double *p1=i+1<nbOfLevsInVec?end:begin;
      const double *p2=i+1<nbOfLevsInVec?third:end;
      double vecPlane[3]={
        (p1[1]-p0[1])*(p2[2]-p1[2])-(p1[2]-p0[2])*(p2[1]-p1[1]),
        (p1[2]-p0[2])*(p2[0]-p1[0])-(p1[0]-p0[0])*(p2[2]-p1[2]),
        (p1[0]-p0[0])*(p2[1]-p1[1])-(p1[1]-p0[1])*(p2[0]-p1[0]),
      };
      double norm=sqrt(vecPlane[0]*vecPlane[0]+vecPlane[1]*vecPlane[1]+vecPlane[2]*vecPlane[2]);
      if(norm>1.e-7)
        {
          vecPlane[0]/=norm; vecPlane[1]/=norm; vecPlane[2]/=norm;
          double norm2=sqrt(vecPlane[0]*vecPlane[0]+vecPlane[1]*vecPlane[1]);
          double vec2[2]={vecPlane[1]/norm2,-vecPlane[0]/norm2};
          double s2=norm2;
          double c2=cos(asin(s2));
          double m[3][3]={
            {vec2[0]*vec2[0]*(1-c2)+c2, vec2[0]*vec2[1]*(1-c2), vec2[1]*s2},
            {vec2[0]*vec2[1]*(1-c2), vec2[1]*vec2[1]*(1-c2)+c2, -vec2[0]*s2},
            {-vec2[1]*s2, vec2[0]*s2, c2}
          };
          double p0r[3]={m[0][0]*p0[0]+m[0][1]*p0[1]+m[0][2]*p0[2], m[1][0]*p0[0]+m[1][1]*p0[1]+m[1][2]*p0[2], m[2][0]*p0[0]+m[2][1]*p0[1]+m[2][2]*p0[2]};
          double p1r[3]={m[0][0]*p1[0]+m[0][1]*p1[1]+m[0][2]*p1[2], m[1][0]*p1[0]+m[1][1]*p1[1]+m[1][2]*p1[2], m[2][0]*p1[0]+m[2][1]*p1[1]+m[2][2]*p1[2]};
          double p2r[3]={m[0][0]*p2[0]+m[0][1]*p2[1]+m[0][2]*p2[2], m[1][0]*p2[0]+m[1][1]*p2[1]+m[1][2]*p2[2], m[2][0]*p2[0]+m[2][1]*p2[1]+m[2][2]*p2[2]};
          INTERP_KERNEL::EdgeArcCircle::GetArcOfCirclePassingThru(p0r,p1r,p2r,tmp3,radius,alpha,alpha0);
          double cosangle=i+1<nbOfLevsInVec?(p0r[0]-tmp3[0])*(p1r[0]-tmp3[0])+(p0r[1]-tmp3[1])*(p1r[1]-tmp3[1]):(p2r[0]-tmp3[0])*(p1r[0]-tmp3[0])+(p2r[1]-tmp3[1])*(p1r[1]-tmp3[1]);
          double angle=acos(cosangle/(radius*radius));
          tmp->rotate(end,vecPlane,angle);
        }
      retPtr=std::copy(tmp2->getConstPointer(),tmp2->getConstPointer()+tmp2->getNbOfElems(),retPtr);
    }
  return ret.retn();
}

/*!
 * This method is private because not easy to use for end user. This method is const contrary to
 * MEDCouplingUMesh::buildExtrudedMesh method because this->_coords are expected to contain
 * the coords sorted slice by slice.
 * \param isQuad specifies presence of quadratic cells.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildExtrudedMeshFromThisLowLev(int nbOfNodesOf1Lev, bool isQuad) const
{
  int nbOf1DCells(getNumberOfNodes()/nbOfNodesOf1Lev-1);
  int nbOf2DCells(getNumberOfCells());
  int nbOf3DCells(nbOf2DCells*nbOf1DCells);
  MEDCouplingUMesh *ret(MEDCouplingUMesh::New("Extruded",getMeshDimension()+1));
  const int *conn(_nodal_connec->begin()),*connI(_nodal_connec_index->begin());
  MCAuto<DataArrayInt> newConn(DataArrayInt::New()),newConnI(DataArrayInt::New());
  newConnI->alloc(nbOf3DCells+1,1);
  int *newConnIPtr(newConnI->getPointer());
  *newConnIPtr++=0;
  std::vector<int> newc;
  for(int j=0;j<nbOf2DCells;j++)
    {
      AppendExtrudedCell(conn+connI[j],conn+connI[j+1],nbOfNodesOf1Lev,isQuad,newc);
      *newConnIPtr++=(int)newc.size();
    }
  newConn->alloc((int)(newc.size())*nbOf1DCells,1);
  int *newConnPtr(newConn->getPointer());
  int deltaPerLev(isQuad?2*nbOfNodesOf1Lev:nbOfNodesOf1Lev);
  newConnIPtr=newConnI->getPointer();
  for(int iz=0;iz<nbOf1DCells;iz++)
    {
      if(iz!=0)
        std::transform(newConnIPtr+1,newConnIPtr+1+nbOf2DCells,newConnIPtr+1+iz*nbOf2DCells,std::bind2nd(std::plus<int>(),newConnIPtr[iz*nbOf2DCells]));
      const int *posOfTypeOfCell(newConnIPtr);
      for(std::vector<int>::const_iterator iter=newc.begin();iter!=newc.end();iter++,newConnPtr++)
        {
          int icell((int)(iter-newc.begin()));//std::distance unfortunately cannot been called here in C++98
          if(icell!=*posOfTypeOfCell)
            {
              if(*iter!=-1)
                *newConnPtr=(*iter)+iz*deltaPerLev;
              else
                *newConnPtr=-1;
            }
          else
            {
              *newConnPtr=*iter;
              posOfTypeOfCell++;
            }
        }
    }
  ret->setConnectivity(newConn,newConnI,true);
  ret->setCoords(getCoords());
  return ret;
}


/*!
 * This method find in candidate pool defined by 'candidates' the cells equal following the polycy 'compType'.
 * If any true is returned and the results will be put at the end of 'result' output parameter. If not false is returned
 * and result remains unchanged.
 * The semantic of 'compType' is specified in MEDCouplingPointSet::zipConnectivityTraducer method.
 * If in 'candidates' pool -1 value is considered as an empty value.
 * WARNING this method returns only ONE set of result !
 */
bool MEDCouplingUMesh::AreCellsEqualInPool(const std::vector<int>& candidates, int compType, const int *conn, const int *connI, DataArrayInt *result)
{
  if(candidates.size()<1)
    return false;
  bool ret=false;
  std::vector<int>::const_iterator iter=candidates.begin();
  int start=(*iter++);
  for(;iter!=candidates.end();iter++)
    {
      int status=AreCellsEqual(conn,connI,start,*iter,compType);
      if(status!=0)
        {
          if(!ret)
            {
              result->pushBackSilent(start);
              ret=true;
            }
          if(status==1)
            result->pushBackSilent(*iter);
          else
            result->pushBackSilent(status==2?(*iter+1):-(*iter+1));
        }
    }
  return ret;
}

/*!
 * This is the low algorithm of MEDCouplingUMesh::buildPartOfMySelf.
 * Keeps from \a this only cells which constituing point id are in the ids specified by [ \a begin,\a end ).
 * The return newly allocated mesh will share the same coordinates as \a this.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelfKeepCoords(const int *begin, const int *end) const
{
  checkConnectivityFullyDefined();
  int ncell=getNumberOfCells();
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New();
  ret->_mesh_dim=_mesh_dim;
  ret->setCoords(_coords);
  std::size_t nbOfElemsRet=std::distance(begin,end);
  int *connIndexRet=(int *)malloc((nbOfElemsRet+1)*sizeof(int));
  connIndexRet[0]=0;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int newNbring=0;
  for(const int *work=begin;work!=end;work++,newNbring++)
    {
      if(*work>=0 && *work<ncell)
        connIndexRet[newNbring+1]=connIndexRet[newNbring]+connIndex[*work+1]-connIndex[*work];
      else
        {
          free(connIndexRet);
          std::ostringstream oss; oss << "MEDCouplingUMesh::buildPartOfMySelfKeepCoords : On pos #" << std::distance(begin,work) << " input cell id =" << *work << " should be in [0," << ncell << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  int *connRet=(int *)malloc(connIndexRet[nbOfElemsRet]*sizeof(int));
  int *connRetWork=connRet;
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(const int *work=begin;work!=end;work++)
    {
      types.insert((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*work]]);
      connRetWork=std::copy(conn+connIndex[*work],conn+connIndex[*work+1],connRetWork);
    }
  MCAuto<DataArrayInt> connRetArr=DataArrayInt::New();
  connRetArr->useArray(connRet,true,C_DEALLOC,connIndexRet[nbOfElemsRet],1);
  MCAuto<DataArrayInt> connIndexRetArr=DataArrayInt::New();
  connIndexRetArr->useArray(connIndexRet,true,C_DEALLOC,(int)nbOfElemsRet+1,1);
  ret->setConnectivity(connRetArr,connIndexRetArr,false);
  ret->_types=types;
  ret->copyTinyInfoFrom(this);
  return ret.retn();
}

/*!
 * This is the low algorithm of MEDCouplingUMesh::buildPartOfMySelfSlice.
 * CellIds are given using range specified by a start an end and step.
 */
MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelfKeepCoordsSlice(int start, int end, int step) const
{
  checkFullyDefined();
  int ncell=getNumberOfCells();
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New();
  ret->_mesh_dim=_mesh_dim;
  ret->setCoords(_coords);
  int newNbOfCells=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"MEDCouplingUMesh::buildPartOfMySelfKeepCoordsSlice : ");
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New(); newConnI->alloc(newNbOfCells+1,1);
  int *newConnIPtr=newConnI->getPointer(); *newConnIPtr=0;
  int work=start;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  for(int i=0;i<newNbOfCells;i++,newConnIPtr++,work+=step)
    {
      if(work>=0 && work<ncell)
        {
          newConnIPtr[1]=newConnIPtr[0]+connIndex[work+1]-connIndex[work];
        }
      else
        {
          std::ostringstream oss; oss << "MEDCouplingUMesh::buildPartOfMySelfKeepCoordsSlice : On pos #" << i << " input cell id =" << work << " should be in [0," << ncell << ") !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
    }
  MCAuto<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(newConnIPtr[0],1);
  int *newConnPtr=newConn->getPointer();
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  work=start;
  for(int i=0;i<newNbOfCells;i++,newConnIPtr++,work+=step)
    {
      types.insert((INTERP_KERNEL::NormalizedCellType)conn[connIndex[work]]);
      newConnPtr=std::copy(conn+connIndex[work],conn+connIndex[work+1],newConnPtr);
    }
  ret->setConnectivity(newConn,newConnI,false);
  ret->_types=types;
  ret->copyTinyInfoFrom(this);
  return ret.retn();
}


int MEDCouplingFastNbrer(int id, unsigned nb, const INTERP_KERNEL::CellModel& cm, bool compute, const int *conn1, const int *conn2)
{
  return id;
}

int MEDCouplingOrientationSensitiveNbrer(int id, unsigned nb, const INTERP_KERNEL::CellModel& cm, bool compute, const int *conn1, const int *conn2)
{
  if(!compute)
    return id+1;
  else
    {
      if(cm.getOrientationStatus(nb,conn1,conn2))
        return id+1;
      else
        return -(id+1);
    }
}


/*!
 * Implementes \a conversionType 0 for meshes with meshDim = 1, of MEDCouplingUMesh::convertLinearCellsToQuadratic method.
 * \return a newly created DataArrayInt instance that the caller should deal with containing cell ids of converted cells.
 * \sa MEDCouplingUMesh::convertLinearCellsToQuadratic.
 */
DataArrayInt *MEDCouplingUMesh::convertLinearCellsToQuadratic1D0(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const
{
  MCAuto<DataArrayDouble> bary=computeCellCenterOfMass();
  MCAuto<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(0,1);
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New(); newConnI->alloc(1,1); newConnI->setIJ(0,0,0);
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(0,1);
  int nbOfCells=getNumberOfCells();
  int nbOfNodes=getNumberOfNodes();
  const int *cPtr=_nodal_connec->begin();
  const int *icPtr=_nodal_connec_index->begin();
  int lastVal=0,offset=nbOfNodes;
  for(int i=0;i<nbOfCells;i++,icPtr++)
    {
      INTERP_KERNEL::NormalizedCellType type=(INTERP_KERNEL::NormalizedCellType)cPtr[*icPtr];
      if(type==INTERP_KERNEL::NORM_SEG2)
        {
          types.insert(INTERP_KERNEL::NORM_SEG3);
          newConn->pushBackSilent((int)INTERP_KERNEL::NORM_SEG3);
          newConn->pushBackValsSilent(cPtr+icPtr[0]+1,cPtr+icPtr[0]+3);
          newConn->pushBackSilent(offset++);
          lastVal+=4;
          newConnI->pushBackSilent(lastVal);
          ret->pushBackSilent(i);
        }
      else
        {
          types.insert(type);
          lastVal+=(icPtr[1]-icPtr[0]);
          newConnI->pushBackSilent(lastVal);
          newConn->pushBackValsSilent(cPtr+icPtr[0],cPtr+icPtr[1]);
        }
    }
  MCAuto<DataArrayDouble> tmp=bary->selectByTupleIdSafe(ret->begin(),ret->end());
  coords=DataArrayDouble::Aggregate(getCoords(),tmp); conn=newConn.retn(); connI=newConnI.retn();
  return ret.retn();
}

DataArrayInt *MEDCouplingUMesh::convertLinearCellsToQuadratic2DAnd3D0(const MEDCouplingUMesh *m1D, const DataArrayInt *desc, const DataArrayInt *descI, DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const
{
  MCAuto<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(0,1);
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New(); newConnI->alloc(1,1); newConnI->setIJ(0,0,0);
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(0,1);
  //
  const int *descPtr(desc->begin()),*descIPtr(descI->begin());
  DataArrayInt *conn1D=0,*conn1DI=0;
  std::set<INTERP_KERNEL::NormalizedCellType> types1D;
  DataArrayDouble *coordsTmp=0;
  MCAuto<DataArrayInt> ret1D=m1D->convertLinearCellsToQuadratic1D0(conn1D,conn1DI,coordsTmp,types1D); ret1D=0;
  MCAuto<DataArrayDouble> coordsTmpSafe(coordsTmp);
  MCAuto<DataArrayInt> conn1DSafe(conn1D),conn1DISafe(conn1DI);
  const int *c1DPtr=conn1D->begin();
  const int *c1DIPtr=conn1DI->begin();
  int nbOfCells=getNumberOfCells();
  const int *cPtr=_nodal_connec->begin();
  const int *icPtr=_nodal_connec_index->begin();
  int lastVal=0;
  for(int i=0;i<nbOfCells;i++,icPtr++,descIPtr++)
    {
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)cPtr[*icPtr];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
      if(!cm.isQuadratic())
        {
          INTERP_KERNEL::NormalizedCellType typ2=cm.getQuadraticType();
          types.insert(typ2); newConn->pushBackSilent(typ2);
          newConn->pushBackValsSilent(cPtr+icPtr[0]+1,cPtr+icPtr[1]);
          for(const int *d=descPtr+descIPtr[0];d!=descPtr+descIPtr[1];d++)
            newConn->pushBackSilent(c1DPtr[c1DIPtr[*d]+3]);
          lastVal+=(icPtr[1]-icPtr[0])+(descIPtr[1]-descIPtr[0]);
          newConnI->pushBackSilent(lastVal);
          ret->pushBackSilent(i);
        }
      else
        {
          types.insert(typ);
          lastVal+=(icPtr[1]-icPtr[0]);
          newConnI->pushBackSilent(lastVal);
          newConn->pushBackValsSilent(cPtr+icPtr[0],cPtr+icPtr[1]);
        }
    }
  conn=newConn.retn(); connI=newConnI.retn(); coords=coordsTmpSafe.retn();
  return ret.retn();
}

/*!
 * Implementes \a conversionType 0 for meshes with meshDim = 2, of MEDCouplingUMesh::convertLinearCellsToQuadratic method.
 * \return a newly created DataArrayInt instance that the caller should deal with containing cell ids of converted cells.
 * \sa MEDCouplingUMesh::convertLinearCellsToQuadratic.
 */
DataArrayInt *MEDCouplingUMesh::convertLinearCellsToQuadratic2D0(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const
{
  MCAuto<DataArrayInt> desc(DataArrayInt::New()),descI(DataArrayInt::New()),tmp2(DataArrayInt::New()),tmp3(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> m1D=buildDescendingConnectivity(desc,descI,tmp2,tmp3); tmp2=0; tmp3=0;
  return convertLinearCellsToQuadratic2DAnd3D0(m1D,desc,descI,conn,connI,coords,types);
}

DataArrayInt *MEDCouplingUMesh::convertLinearCellsToQuadratic2D1(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const
{
  MCAuto<DataArrayInt> desc(DataArrayInt::New()),descI(DataArrayInt::New()),tmp2(DataArrayInt::New()),tmp3(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> m1D=buildDescendingConnectivity(desc,descI,tmp2,tmp3); tmp2=0; tmp3=0;
  //
  MCAuto<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(0,1);
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New(); newConnI->alloc(1,1); newConnI->setIJ(0,0,0);
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(0,1);
  //
  MCAuto<DataArrayDouble> bary=computeCellCenterOfMass();
  const int *descPtr(desc->begin()),*descIPtr(descI->begin());
  DataArrayInt *conn1D=0,*conn1DI=0;
  std::set<INTERP_KERNEL::NormalizedCellType> types1D;
  DataArrayDouble *coordsTmp=0;
  MCAuto<DataArrayInt> ret1D=m1D->convertLinearCellsToQuadratic1D0(conn1D,conn1DI,coordsTmp,types1D); ret1D=0;
  MCAuto<DataArrayDouble> coordsTmpSafe(coordsTmp);
  MCAuto<DataArrayInt> conn1DSafe(conn1D),conn1DISafe(conn1DI);
  const int *c1DPtr=conn1D->begin();
  const int *c1DIPtr=conn1DI->begin();
  int nbOfCells=getNumberOfCells();
  const int *cPtr=_nodal_connec->begin();
  const int *icPtr=_nodal_connec_index->begin();
  int lastVal=0,offset=coordsTmpSafe->getNumberOfTuples();
  for(int i=0;i<nbOfCells;i++,icPtr++,descIPtr++)
    {
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)cPtr[*icPtr];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
      if(!cm.isQuadratic())
        {
          INTERP_KERNEL::NormalizedCellType typ2=cm.getQuadraticType2();
          types.insert(typ2); newConn->pushBackSilent(typ2);
          newConn->pushBackValsSilent(cPtr+icPtr[0]+1,cPtr+icPtr[1]);
          for(const int *d=descPtr+descIPtr[0];d!=descPtr+descIPtr[1];d++)
            newConn->pushBackSilent(c1DPtr[c1DIPtr[*d]+3]);
          newConn->pushBackSilent(offset+ret->getNumberOfTuples());
          lastVal+=(icPtr[1]-icPtr[0])+(descIPtr[1]-descIPtr[0])+1;
          newConnI->pushBackSilent(lastVal);
          ret->pushBackSilent(i);
        }
      else
        {
          types.insert(typ);
          lastVal+=(icPtr[1]-icPtr[0]);
          newConnI->pushBackSilent(lastVal);
          newConn->pushBackValsSilent(cPtr+icPtr[0],cPtr+icPtr[1]);
        }
    }
  MCAuto<DataArrayDouble> tmp=bary->selectByTupleIdSafe(ret->begin(),ret->end());
  coords=DataArrayDouble::Aggregate(coordsTmpSafe,tmp); conn=newConn.retn(); connI=newConnI.retn();
  return ret.retn();
}

/*!
 * Implementes \a conversionType 0 for meshes with meshDim = 3, of MEDCouplingUMesh::convertLinearCellsToQuadratic method.
 * \return a newly created DataArrayInt instance that the caller should deal with containing cell ids of converted cells.
 * \sa MEDCouplingUMesh::convertLinearCellsToQuadratic.
 */
DataArrayInt *MEDCouplingUMesh::convertLinearCellsToQuadratic3D0(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const
{
  MCAuto<DataArrayInt> desc(DataArrayInt::New()),descI(DataArrayInt::New()),tmp2(DataArrayInt::New()),tmp3(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> m1D=explode3DMeshTo1D(desc,descI,tmp2,tmp3); tmp2=0; tmp3=0;
  return convertLinearCellsToQuadratic2DAnd3D0(m1D,desc,descI,conn,connI,coords,types);
}

DataArrayInt *MEDCouplingUMesh::convertLinearCellsToQuadratic3D1(DataArrayInt *&conn, DataArrayInt *&connI, DataArrayDouble *& coords, std::set<INTERP_KERNEL::NormalizedCellType>& types) const
{
  MCAuto<DataArrayInt> desc2(DataArrayInt::New()),desc2I(DataArrayInt::New()),tmp2(DataArrayInt::New()),tmp3(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> m2D=buildDescendingConnectivityGen<MinusOneSonsGeneratorBiQuadratic>(desc2,desc2I,tmp2,tmp3,MEDCouplingFastNbrer); tmp2=0; tmp3=0;
  MCAuto<DataArrayInt> desc1(DataArrayInt::New()),desc1I(DataArrayInt::New()),tmp4(DataArrayInt::New()),tmp5(DataArrayInt::New());
  MCAuto<MEDCouplingUMesh> m1D=explode3DMeshTo1D(desc1,desc1I,tmp4,tmp5); tmp4=0; tmp5=0;
  //
  MCAuto<DataArrayInt> newConn=DataArrayInt::New(); newConn->alloc(0,1);
  MCAuto<DataArrayInt> newConnI=DataArrayInt::New(); newConnI->alloc(1,1); newConnI->setIJ(0,0,0);
  MCAuto<DataArrayInt> ret=DataArrayInt::New(),ret2=DataArrayInt::New(); ret->alloc(0,1); ret2->alloc(0,1);
  //
  MCAuto<DataArrayDouble> bary=computeCellCenterOfMass();
  const int *descPtr(desc1->begin()),*descIPtr(desc1I->begin()),*desc2Ptr(desc2->begin()),*desc2IPtr(desc2I->begin());
  DataArrayInt *conn1D=0,*conn1DI=0,*conn2D=0,*conn2DI=0;
  std::set<INTERP_KERNEL::NormalizedCellType> types1D,types2D;
  DataArrayDouble *coordsTmp=0,*coordsTmp2=0;
  MCAuto<DataArrayInt> ret1D=m1D->convertLinearCellsToQuadratic1D0(conn1D,conn1DI,coordsTmp,types1D); ret1D=DataArrayInt::New(); ret1D->alloc(0,1);
  MCAuto<DataArrayInt> conn1DSafe(conn1D),conn1DISafe(conn1DI);
  MCAuto<DataArrayDouble> coordsTmpSafe(coordsTmp);
  MCAuto<DataArrayInt> ret2D=m2D->convertLinearCellsToQuadratic2D1(conn2D,conn2DI,coordsTmp2,types2D); ret2D=DataArrayInt::New(); ret2D->alloc(0,1);
  MCAuto<DataArrayDouble> coordsTmp2Safe(coordsTmp2);
  MCAuto<DataArrayInt> conn2DSafe(conn2D),conn2DISafe(conn2DI);
  const int *c1DPtr=conn1D->begin(),*c1DIPtr=conn1DI->begin(),*c2DPtr=conn2D->begin(),*c2DIPtr=conn2DI->begin();
  int nbOfCells=getNumberOfCells();
  const int *cPtr=_nodal_connec->begin();
  const int *icPtr=_nodal_connec_index->begin();
  int lastVal=0,offset=coordsTmpSafe->getNumberOfTuples();
  for(int i=0;i<nbOfCells;i++,icPtr++,descIPtr++,desc2IPtr++)
    {
      INTERP_KERNEL::NormalizedCellType typ=(INTERP_KERNEL::NormalizedCellType)cPtr[*icPtr];
      const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(typ);
      if(!cm.isQuadratic())
        {
          INTERP_KERNEL::NormalizedCellType typ2=cm.getQuadraticType2();
          if(typ2==INTERP_KERNEL::NORM_ERROR)
            {
              std::ostringstream oss; oss << "MEDCouplingUMesh::convertLinearCellsToQuadratic3D1 : On cell #" << i << " the linear cell type does not support advanced quadratization !";
              throw INTERP_KERNEL::Exception(oss.str());
            }
          types.insert(typ2); newConn->pushBackSilent(typ2);
          newConn->pushBackValsSilent(cPtr+icPtr[0]+1,cPtr+icPtr[1]);
          for(const int *d=descPtr+descIPtr[0];d!=descPtr+descIPtr[1];d++)
            newConn->pushBackSilent(c1DPtr[c1DIPtr[*d]+3]);
          for(const int *d=desc2Ptr+desc2IPtr[0];d!=desc2Ptr+desc2IPtr[1];d++)
            {
              int nodeId2=c2DPtr[c2DIPtr[(*d)+1]-1];
              int tmpPos=newConn->getNumberOfTuples();
              newConn->pushBackSilent(nodeId2);
              ret2D->pushBackSilent(nodeId2); ret1D->pushBackSilent(tmpPos);
            }
          newConn->pushBackSilent(offset+ret->getNumberOfTuples());
          lastVal+=(icPtr[1]-icPtr[0])+(descIPtr[1]-descIPtr[0])+(desc2IPtr[1]-desc2IPtr[0])+1;
          newConnI->pushBackSilent(lastVal);
          ret->pushBackSilent(i);
        }
      else
        {
          types.insert(typ);
          lastVal+=(icPtr[1]-icPtr[0]);
          newConnI->pushBackSilent(lastVal);
          newConn->pushBackValsSilent(cPtr+icPtr[0],cPtr+icPtr[1]);
        }
    }
  MCAuto<DataArrayInt> diffRet2D=ret2D->getDifferentValues();
  MCAuto<DataArrayInt> o2nRet2D=diffRet2D->invertArrayN2O2O2N(coordsTmp2Safe->getNumberOfTuples());
  coordsTmp2Safe=coordsTmp2Safe->selectByTupleId(diffRet2D->begin(),diffRet2D->end());
  MCAuto<DataArrayDouble> tmp=bary->selectByTupleIdSafe(ret->begin(),ret->end());
  std::vector<const DataArrayDouble *> v(3); v[0]=coordsTmpSafe; v[1]=coordsTmp2Safe; v[2]=tmp;
  int *c=newConn->getPointer();
  const int *cI(newConnI->begin());
  for(const int *elt=ret1D->begin();elt!=ret1D->end();elt++)
    c[*elt]=o2nRet2D->getIJ(c[*elt],0)+offset;
  offset=coordsTmp2Safe->getNumberOfTuples();
  for(const int *elt=ret->begin();elt!=ret->end();elt++)
    c[cI[(*elt)+1]-1]+=offset;
  coords=DataArrayDouble::Aggregate(v); conn=newConn.retn(); connI=newConnI.retn();
  return ret.retn();
}

DataArrayInt *MEDCouplingUMesh::buildUnionOf2DMeshLinear(const MEDCouplingUMesh *skin, const DataArrayInt *n2o) const
{
  int nbOfNodesExpected(skin->getNumberOfNodes());
  const int *n2oPtr(n2o->begin());
  MCAuto<DataArrayInt> revNodal(DataArrayInt::New()),revNodalI(DataArrayInt::New());
  skin->getReverseNodalConnectivity(revNodal,revNodalI);
  const int *revNodalPtr(revNodal->begin()),*revNodalIPtr(revNodalI->begin());
  const int *nodalPtr(skin->getNodalConnectivity()->begin());
  const int *nodalIPtr(skin->getNodalConnectivityIndex()->begin());
  MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(nbOfNodesExpected+1,1);
  int *work(ret->getPointer());  *work++=INTERP_KERNEL::NORM_POLYGON;
  if(nbOfNodesExpected<1)
    return ret.retn();
  int prevCell(0),prevNode(nodalPtr[nodalIPtr[0]+1]);
  *work++=n2oPtr[prevNode];
  for(int i=1;i<nbOfNodesExpected;i++)
    {
      if(nodalIPtr[prevCell+1]-nodalIPtr[prevCell]==3)
        {
          std::set<int> conn(nodalPtr+nodalIPtr[prevCell]+1,nodalPtr+nodalIPtr[prevCell]+3);
          conn.erase(prevNode);
          if(conn.size()==1)
            {
              int curNode(*(conn.begin()));
              *work++=n2oPtr[curNode];
              std::set<int> shar(revNodalPtr+revNodalIPtr[curNode],revNodalPtr+revNodalIPtr[curNode+1]);
              shar.erase(prevCell);
              if(shar.size()==1)
                {
                  prevCell=*(shar.begin());
                  prevNode=curNode;
                }
              else
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMeshLinear : presence of unexpected 2 !");
            }
          else
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMeshLinear : presence of unexpected 1 !");
        }
      else
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMeshLinear : presence of unexpected cell !");
    }
  return ret.retn();
}

DataArrayInt *MEDCouplingUMesh::buildUnionOf2DMeshQuadratic(const MEDCouplingUMesh *skin, const DataArrayInt *n2o) const
{
  int nbOfNodesExpected(skin->getNumberOfNodes());
  int nbOfTurn(nbOfNodesExpected/2);
  const int *n2oPtr(n2o->begin());
  MCAuto<DataArrayInt> revNodal(DataArrayInt::New()),revNodalI(DataArrayInt::New());
  skin->getReverseNodalConnectivity(revNodal,revNodalI);
  const int *revNodalPtr(revNodal->begin()),*revNodalIPtr(revNodalI->begin());
  const int *nodalPtr(skin->getNodalConnectivity()->begin());
  const int *nodalIPtr(skin->getNodalConnectivityIndex()->begin());
  MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(nbOfNodesExpected+1,1);
  int *work(ret->getPointer());  *work++=INTERP_KERNEL::NORM_QPOLYG;
  if(nbOfNodesExpected<1)
    return ret.retn();
  int prevCell(0),prevNode(nodalPtr[nodalIPtr[0]+1]);
  *work=n2oPtr[prevNode]; work[nbOfTurn]=n2oPtr[nodalPtr[nodalIPtr[0]+3]]; work++;
  for(int i=1;i<nbOfTurn;i++)
    {
      if(nodalIPtr[prevCell+1]-nodalIPtr[prevCell]==4)
        {
          std::set<int> conn(nodalPtr+nodalIPtr[prevCell]+1,nodalPtr+nodalIPtr[prevCell]+3);
          conn.erase(prevNode);
          if(conn.size()==1)
            {
              int curNode(*(conn.begin()));
              *work=n2oPtr[curNode];
              std::set<int> shar(revNodalPtr+revNodalIPtr[curNode],revNodalPtr+revNodalIPtr[curNode+1]);
              shar.erase(prevCell);
              if(shar.size()==1)
                {
                  int curCell(*(shar.begin()));
                  work[nbOfTurn]=n2oPtr[nodalPtr[nodalIPtr[curCell]+3]];
                  prevCell=curCell;
                  prevNode=curNode;
                  work++;
                }
              else
                throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMeshQuadratic : presence of unexpected 2 !");
            }
          else
            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMeshQuadratic : presence of unexpected 1 !");
        }
      else
        throw INTERP_KERNEL::Exception("MEDCouplingUMesh::buildUnionOf2DMeshQuadratic : presence of unexpected cell !");
    }
  return ret.retn();
}

MEDCouplingUMesh *MEDCouplingUMesh::MergeUMeshesLL(const std::vector<const MEDCouplingUMesh *>& a)
{
  if(a.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::MergeUMeshes : input array must be NON EMPTY !");
  std::vector<const MEDCouplingUMesh *>::const_iterator it=a.begin();
  int meshDim=(*it)->getMeshDimension();
  int nbOfCells=(*it)->getNumberOfCells();
  int meshLgth=(*it++)->getNodalConnectivityArrayLen();
  for(;it!=a.end();it++)
    {
      if(meshDim!=(*it)->getMeshDimension())
        throw INTERP_KERNEL::Exception("Mesh dimensions mismatches, MergeUMeshes impossible !");
      nbOfCells+=(*it)->getNumberOfCells();
      meshLgth+=(*it)->getNodalConnectivityArrayLen();
    }
  std::vector<const MEDCouplingPointSet *> aps(a.size());
  std::copy(a.begin(),a.end(),aps.begin());
  MCAuto<DataArrayDouble> pts=MergeNodesArray(aps);
  MCAuto<MEDCouplingUMesh> ret=MEDCouplingUMesh::New("merge",meshDim);
  ret->setCoords(pts);
  MCAuto<DataArrayInt> c=DataArrayInt::New();
  c->alloc(meshLgth,1);
  int *cPtr=c->getPointer();
  MCAuto<DataArrayInt> cI=DataArrayInt::New();
  cI->alloc(nbOfCells+1,1);
  int *cIPtr=cI->getPointer();
  *cIPtr++=0;
  int offset=0;
  int offset2=0;
  for(it=a.begin();it!=a.end();it++)
    {
      int curNbOfCell=(*it)->getNumberOfCells();
      const int *curCI=(*it)->_nodal_connec_index->begin();
      const int *curC=(*it)->_nodal_connec->begin();
      cIPtr=std::transform(curCI+1,curCI+curNbOfCell+1,cIPtr,std::bind2nd(std::plus<int>(),offset));
      for(int j=0;j<curNbOfCell;j++)
        {
          const int *src=curC+curCI[j];
          *cPtr++=*src++;
          for(;src!=curC+curCI[j+1];src++,cPtr++)
            {
              if(*src!=-1)
                *cPtr=*src+offset2;
              else
                *cPtr=-1;
            }
        }
      offset+=curCI[curNbOfCell];
      offset2+=(*it)->getNumberOfNodes();
    }
  //
  ret->setConnectivity(c,cI,true);
  return ret.retn();
}


/*!
 * \param [in] pt the start pointer (included) of the coordinates of the point
 * \param [in] cellIdsBg the start pointer (included) of cellIds
 * \param [in] cellIdsEnd the end pointer (excluded) of cellIds
 * \param [in] nc nodal connectivity
 * \param [in] ncI nodal connectivity index
 * \param [in,out] ret0 the min distance between \a this and the external input point
 * \param [out] cellId that corresponds to minimal distance. If the closer node is not linked to any cell in \a this -1 is returned.
 * \sa MEDCouplingUMesh::distanceToPoint, MEDCouplingUMesh::distanceToPoints
 */
void MEDCouplingUMesh::DistanceToPoint3DSurfAlg(const double *pt, const int *cellIdsBg, const int *cellIdsEnd, const double *coords, const int *nc, const int *ncI, double& ret0, int& cellId)
{
  cellId=-1;
  ret0=std::numeric_limits<double>::max();
  for(const int *zeCell=cellIdsBg;zeCell!=cellIdsEnd;zeCell++)
    {
      switch((INTERP_KERNEL::NormalizedCellType)nc[ncI[*zeCell]])
      {
        case INTERP_KERNEL::NORM_TRI3:
          {
            double tmp=INTERP_KERNEL::DistanceFromPtToTriInSpaceDim3(pt,coords+3*nc[ncI[*zeCell]+1],coords+3*nc[ncI[*zeCell]+2],coords+3*nc[ncI[*zeCell]+3]);
            if(tmp<ret0)
              { ret0=tmp; cellId=*zeCell; }
            break;
          }
        case INTERP_KERNEL::NORM_QUAD4:
        case INTERP_KERNEL::NORM_POLYGON:
          {
            double tmp=INTERP_KERNEL::DistanceFromPtToPolygonInSpaceDim3(pt,nc+ncI[*zeCell]+1,nc+ncI[*zeCell+1],coords);
            if(tmp<ret0)
              { ret0=tmp; cellId=*zeCell; }
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::distanceToPoint3DSurfAlg : not managed cell type ! Supporting TRI3, QUAD4 and POLYGON !");
      }
    }
}

/*!
 * \param [in] pt the start pointer (included) of the coordinates of the point
 * \param [in] cellIdsBg the start pointer (included) of cellIds
 * \param [in] cellIdsEnd the end pointer (excluded) of cellIds
 * \param [in] nc nodal connectivity
 * \param [in] ncI nodal connectivity index
 * \param [in,out] ret0 the min distance between \a this and the external input point
 * \param [out] cellId that corresponds to minimal distance. If the closer node is not linked to any cell in \a this -1 is returned.
 * \sa MEDCouplingUMesh::distanceToPoint, MEDCouplingUMesh::distanceToPoints
 */
void MEDCouplingUMesh::DistanceToPoint2DCurveAlg(const double *pt, const int *cellIdsBg, const int *cellIdsEnd, const double *coords, const int *nc, const int *ncI, double& ret0, int& cellId)
{
  cellId=-1;
  ret0=std::numeric_limits<double>::max();
  for(const int *zeCell=cellIdsBg;zeCell!=cellIdsEnd;zeCell++)
    {
      switch((INTERP_KERNEL::NormalizedCellType)nc[ncI[*zeCell]])
      {
        case INTERP_KERNEL::NORM_SEG2:
          {
            std::size_t uselessEntry=0;
            double tmp=INTERP_KERNEL::SquareDistanceFromPtToSegInSpaceDim2(pt,coords+2*nc[ncI[*zeCell]+1],coords+2*nc[ncI[*zeCell]+2],uselessEntry);
            tmp=sqrt(tmp);
            if(tmp<ret0)
              { ret0=tmp; cellId=*zeCell; }
            break;
          }
        default:
          throw INTERP_KERNEL::Exception("MEDCouplingUMesh::distanceToPoint2DCurveAlg : not managed cell type ! Supporting SEG2 !");
      }
    }
}
DataArrayInt *MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeedAlg(std::vector<bool>& fetched, const int *seedBg, const int *seedEnd, const DataArrayInt *arrIn, const DataArrayInt *arrIndxIn, int nbOfDepthPeeling, int& nbOfDepthPeelingPerformed)
{
  nbOfDepthPeelingPerformed=0;
  if(!seedBg || !seedEnd || !arrIn || !arrIndxIn)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeedAlg : some input pointer is NULL !");
  int nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
  std::vector<bool> fetched2(nbOfTuples,false);
  int i=0;
  for(const int *seedElt=seedBg;seedElt!=seedEnd;seedElt++,i++)
    {
      if(*seedElt>=0 && *seedElt<nbOfTuples)
        { fetched[*seedElt]=true; fetched2[*seedElt]=true; }
      else
        { std::ostringstream oss; oss << "MEDCouplingUMesh::ComputeSpreadZoneGraduallyFromSeedAlg : At pos #" << i << " of seeds value is " << *seedElt << "! Should be in [0," << nbOfTuples << ") !"; throw INTERP_KERNEL::Exception(oss.str()); }
    }
  const int *arrInPtr=arrIn->begin();
  const int *arrIndxPtr=arrIndxIn->begin();
  int targetNbOfDepthPeeling=nbOfDepthPeeling!=-1?nbOfDepthPeeling:std::numeric_limits<int>::max();
  std::vector<int> idsToFetch1(seedBg,seedEnd);
  std::vector<int> idsToFetch2;
  std::vector<int> *idsToFetch=&idsToFetch1;
  std::vector<int> *idsToFetchOther=&idsToFetch2;
  while(!idsToFetch->empty() && nbOfDepthPeelingPerformed<targetNbOfDepthPeeling)
    {
      for(std::vector<int>::const_iterator it=idsToFetch->begin();it!=idsToFetch->end();it++)
        for(const int *it2=arrInPtr+arrIndxPtr[*it];it2!=arrInPtr+arrIndxPtr[*it+1];it2++)
          if(!fetched[*it2])
            { fetched[*it2]=true; fetched2[*it2]=true; idsToFetchOther->push_back(*it2); }
      std::swap(idsToFetch,idsToFetchOther);
      idsToFetchOther->clear();
      nbOfDepthPeelingPerformed++;
    }
  int lgth=(int)std::count(fetched2.begin(),fetched2.end(),true);
  i=0;
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(lgth,1);
  int *retPtr=ret->getPointer();
  for(std::vector<bool>::const_iterator it=fetched2.begin();it!=fetched2.end();it++,i++)
    if(*it)
      *retPtr++=i;
  return ret.retn();
}

/*!
 * This method put in zip format into parameter 'zipFrmt' in full interlace mode.
 * This format is often asked by INTERP_KERNEL algorithms to avoid many indirections into coordinates array.
 */
void MEDCouplingUMesh::FillInCompact3DMode(int spaceDim, int nbOfNodesInCell, const int *conn, const double *coo, double *zipFrmt)
{
  double *w=zipFrmt;
  if(spaceDim==3)
    for(int i=0;i<nbOfNodesInCell;i++)
      w=std::copy(coo+3*conn[i],coo+3*conn[i]+3,w);
  else if(spaceDim==2)
    {
      for(int i=0;i<nbOfNodesInCell;i++)
        {
          w=std::copy(coo+2*conn[i],coo+2*conn[i]+2,w);
          *w++=0.;
        }
    }
  else
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::FillInCompact3DMode : Invalid spaceDim specified : must be 2 or 3 !");
}

/*!
 * This method takes in input a cell defined by its MEDcouplingUMesh connectivity [ \a connBg , \a connEnd ) and returns its extruded cell by inserting the result at the end of ret.
 * \param nbOfNodesPerLev in parameter that specifies the number of nodes of one slice of global dataset
 * \param isQuad specifies the policy of connectivity.
 * @ret in/out parameter in which the result will be append
 */
void MEDCouplingUMesh::AppendExtrudedCell(const int *connBg, const int *connEnd, int nbOfNodesPerLev, bool isQuad, std::vector<int>& ret)
{
  INTERP_KERNEL::NormalizedCellType flatType=(INTERP_KERNEL::NormalizedCellType)connBg[0];
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(flatType);
  ret.push_back(cm.getExtrudedType());
  int deltaz=isQuad?2*nbOfNodesPerLev:nbOfNodesPerLev;
  switch(flatType)
  {
    case INTERP_KERNEL::NORM_POINT1:
      {
        ret.push_back(connBg[1]);
        ret.push_back(connBg[1]+nbOfNodesPerLev);
        break;
      }
    case INTERP_KERNEL::NORM_SEG2:
      {
        int conn[4]={connBg[1],connBg[2],connBg[2]+deltaz,connBg[1]+deltaz};
        ret.insert(ret.end(),conn,conn+4);
        break;
      }
    case INTERP_KERNEL::NORM_SEG3:
      {
        int conn[8]={connBg[1],connBg[3],connBg[3]+deltaz,connBg[1]+deltaz,connBg[2],connBg[3]+nbOfNodesPerLev,connBg[2]+deltaz,connBg[1]+nbOfNodesPerLev};
        ret.insert(ret.end(),conn,conn+8);
        break;
      }
    case INTERP_KERNEL::NORM_QUAD4:
      {
        int conn[8]={connBg[1],connBg[2],connBg[3],connBg[4],connBg[1]+deltaz,connBg[2]+deltaz,connBg[3]+deltaz,connBg[4]+deltaz};
        ret.insert(ret.end(),conn,conn+8);
        break;
      }
    case INTERP_KERNEL::NORM_TRI3:
      {
        int conn[6]={connBg[1],connBg[2],connBg[3],connBg[1]+deltaz,connBg[2]+deltaz,connBg[3]+deltaz};
        ret.insert(ret.end(),conn,conn+6);
        break;
      }
    case INTERP_KERNEL::NORM_TRI6:
      {
        int conn[15]={connBg[1],connBg[2],connBg[3],connBg[1]+deltaz,connBg[2]+deltaz,connBg[3]+deltaz,connBg[4],connBg[5],connBg[6],connBg[4]+deltaz,connBg[5]+deltaz,connBg[6]+deltaz,
          connBg[1]+nbOfNodesPerLev,connBg[2]+nbOfNodesPerLev,connBg[3]+nbOfNodesPerLev};
        ret.insert(ret.end(),conn,conn+15);
        break;
      }
    case INTERP_KERNEL::NORM_QUAD8:
      {
        int conn[20]={
          connBg[1],connBg[2],connBg[3],connBg[4],connBg[1]+deltaz,connBg[2]+deltaz,connBg[3]+deltaz,connBg[4]+deltaz,
          connBg[5],connBg[6],connBg[7],connBg[8],connBg[5]+deltaz,connBg[6]+deltaz,connBg[7]+deltaz,connBg[8]+deltaz,
          connBg[1]+nbOfNodesPerLev,connBg[2]+nbOfNodesPerLev,connBg[3]+nbOfNodesPerLev,connBg[4]+nbOfNodesPerLev
        };
        ret.insert(ret.end(),conn,conn+20);
        break;
      }
    case INTERP_KERNEL::NORM_POLYGON:
      {
        std::back_insert_iterator< std::vector<int> > ii(ret);
        std::copy(connBg+1,connEnd,ii);
        *ii++=-1;
        std::reverse_iterator<const int *> rConnBg(connEnd);
        std::reverse_iterator<const int *> rConnEnd(connBg+1);
        std::transform(rConnBg,rConnEnd,ii,std::bind2nd(std::plus<int>(),deltaz));
        std::size_t nbOfRadFaces=std::distance(connBg+1,connEnd);
        for(std::size_t i=0;i<nbOfRadFaces;i++)
          {
            *ii++=-1;
            int conn[4]={connBg[(i+1)%nbOfRadFaces+1],connBg[i+1],connBg[i+1]+deltaz,connBg[(i+1)%nbOfRadFaces+1]+deltaz};
            std::copy(conn,conn+4,ii);
          }
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("A flat type has been detected that has not its extruded representation !");
  }
}


/*!
 * This method is part of the Slice3D algorithm. It is the first step of assembly process, ones coordinates have been computed (by MEDCouplingUMesh::split3DCurveWithPlane method).
 * This method allows to compute given the status of 3D curve cells and the descending connectivity 3DSurf->3DCurve to deduce the intersection of each 3D surf cells
 * with a plane. The result will be put in 'cut3DSuf' out parameter.
 * \param [in] cut3DCurve  input paramter that gives for each 3DCurve cell if it owns fully to the plane or partially.
 * \param [out] nodesOnPlane, returns all the nodes that are on the plane.
 * \param [in] nodal3DSurf is the nodal connectivity of 3D surf mesh.
 * \param [in] nodalIndx3DSurf is the nodal connectivity index of 3D surf mesh.
 * \param [in] nodal3DCurve is the nodal connectivity of 3D curve mesh.
 * \param [in] nodal3DIndxCurve is the nodal connectivity index of 3D curve mesh.
 * \param [in] desc is the descending connectivity 3DSurf->3DCurve
 * \param [in] descIndx is the descending connectivity index 3DSurf->3DCurve
 * \param [out] cut3DSuf input/output param.
 */
void MEDCouplingUMesh::AssemblyForSplitFrom3DCurve(const std::vector<int>& cut3DCurve, std::vector<int>& nodesOnPlane, const int *nodal3DSurf, const int *nodalIndx3DSurf,
                                                   const int *nodal3DCurve, const int *nodalIndx3DCurve,
                                                   const int *desc, const int *descIndx,
                                                   std::vector< std::pair<int,int> >& cut3DSurf)
{
  std::set<int> nodesOnP(nodesOnPlane.begin(),nodesOnPlane.end());
  int nbOf3DSurfCell=(int)cut3DSurf.size();
  for(int i=0;i<nbOf3DSurfCell;i++)
    {
      std::vector<int> res;
      int offset=descIndx[i];
      int nbOfSeg=descIndx[i+1]-offset;
      for(int j=0;j<nbOfSeg;j++)
        {
          int edgeId=desc[offset+j];
          int status=cut3DCurve[edgeId];
          if(status!=-2)
            {
              if(status>-1)
                res.push_back(status);
              else
                {
                  res.push_back(nodal3DCurve[nodalIndx3DCurve[edgeId]+1]);
                  res.push_back(nodal3DCurve[nodalIndx3DCurve[edgeId]+2]);
                }
            }
        }
      switch(res.size())
      {
        case 2:
          {
            cut3DSurf[i].first=res[0]; cut3DSurf[i].second=res[1];
            break;
          }
        case 1:
        case 0:
          {
            std::set<int> s1(nodal3DSurf+nodalIndx3DSurf[i]+1,nodal3DSurf+nodalIndx3DSurf[i+1]);
            std::set_intersection(nodesOnP.begin(),nodesOnP.end(),s1.begin(),s1.end(),std::back_insert_iterator< std::vector<int> >(res));
            if(res.size()==2)
              {
                cut3DSurf[i].first=res[0]; cut3DSurf[i].second=res[1];
              }
            else
              {
                cut3DSurf[i].first=-1; cut3DSurf[i].second=-1;
              }
            break;
          }
        default:
          {// case when plane is on a multi colinear edge of a polyhedron
            if((int)res.size()==2*nbOfSeg)
              {
                cut3DSurf[i].first=-2; cut3DSurf[i].second=i;
              }
            else
              throw INTERP_KERNEL::Exception("MEDCouplingUMesh::AssemblyPointsFrom3DCurve : unexpected situation !");
          }
      }
    }
}


/*!
 * \a this is expected to be a mesh with spaceDim==3 and meshDim==3. If not an exception will be thrown.
 * This method is part of the Slice3D algorithm. It is the second step of assembly process, ones coordinates have been computed (by MEDCouplingUMesh::split3DCurveWithPlane method).
 * This method allows to compute given the result of 3D surf cells with plane and the descending connectivity 3D->3DSurf to deduce the intersection of each 3D cells
 * with a plane. The result will be put in 'nodalRes' 'nodalResIndx' and 'cellIds' out parameters.
 * \param cut3DSurf  input paramter that gives for each 3DSurf its intersection with plane (result of MEDCouplingUMesh::AssemblyForSplitFrom3DCurve).
 * \param desc is the descending connectivity 3D->3DSurf
 * \param descIndx is the descending connectivity index 3D->3DSurf
 */
void MEDCouplingUMesh::assemblyForSplitFrom3DSurf(const std::vector< std::pair<int,int> >& cut3DSurf,
                                                  const int *desc, const int *descIndx,
                                                  DataArrayInt *nodalRes, DataArrayInt *nodalResIndx, DataArrayInt *cellIds) const
{
  checkFullyDefined();
  if(getMeshDimension()!=3 || getSpaceDimension()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingUMesh::assemblyForSplitFrom3DSurf works on umeshes with meshdim equal to 3 and spaceDim equal to 3 too!");
  const int *nodal3D(_nodal_connec->begin()),*nodalIndx3D(_nodal_connec_index->begin());
  int nbOfCells(getNumberOfCells());
  for(int i=0;i<nbOfCells;i++)
    {
      std::map<int, std::set<int> > m;
      int offset=descIndx[i];
      int nbOfFaces=descIndx[i+1]-offset;
      int start=-1;
      int end=-1;
      for(int j=0;j<nbOfFaces;j++)
        {
          const std::pair<int,int>& p=cut3DSurf[desc[offset+j]];
          if(p.first!=-1 && p.second!=-1)
            {
              if(p.first!=-2)
                {
                  start=p.first; end=p.second;
                  m[p.first].insert(p.second);
                  m[p.second].insert(p.first);
                }
              else
                {
                  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel((INTERP_KERNEL::NormalizedCellType)nodal3D[nodalIndx3D[i]]);
                  int sz=nodalIndx3D[i+1]-nodalIndx3D[i]-1;
                  INTERP_KERNEL::AutoPtr<int> tmp=new int[sz];
                  INTERP_KERNEL::NormalizedCellType cmsId;
                  unsigned nbOfNodesSon=cm.fillSonCellNodalConnectivity2(j,nodal3D+nodalIndx3D[i]+1,sz,tmp,cmsId);
                  start=tmp[0]; end=tmp[nbOfNodesSon-1];
                  for(unsigned k=0;k<nbOfNodesSon;k++)
                    {
                      m[tmp[k]].insert(tmp[(k+1)%nbOfNodesSon]);
                      m[tmp[(k+1)%nbOfNodesSon]].insert(tmp[k]);
                    }
                }
            }
        }
      if(m.empty())
        continue;
      std::vector<int> conn(1,(int)INTERP_KERNEL::NORM_POLYGON);
      int prev=end;
      while(end!=start)
        {
          std::map<int, std::set<int> >::const_iterator it=m.find(start);
          const std::set<int>& s=(*it).second;
          std::set<int> s2; s2.insert(prev);
          std::set<int> s3;
          std::set_difference(s.begin(),s.end(),s2.begin(),s2.end(),inserter(s3,s3.begin()));
          if(s3.size()==1)
            {
              int val=*s3.begin();
              conn.push_back(start);
              prev=start;
              start=val;
            }
          else
            start=end;
        }
      conn.push_back(end);
      if(conn.size()>3)
        {
          nodalRes->insertAtTheEnd(conn.begin(),conn.end());
          nodalResIndx->pushBackSilent(nodalRes->getNumberOfTuples());
          cellIds->pushBackSilent(i);
        }
    }
}


void MEDCouplingUMesh::ComputeAllTypesInternal(std::set<INTERP_KERNEL::NormalizedCellType>& types, const DataArrayInt *nodalConnec, const DataArrayInt *nodalConnecIndex)
{
  if(nodalConnec && nodalConnecIndex)
    {
      types.clear();
      const int *conn(nodalConnec->begin()),*connIndex(nodalConnecIndex->begin());
      int nbOfElem(nodalConnecIndex->getNbOfElems()-1);
      if(nbOfElem>0)
        for(const int *pt=connIndex;pt!=connIndex+nbOfElem;pt++)
          types.insert((INTERP_KERNEL::NormalizedCellType)conn[*pt]);
    }
}
