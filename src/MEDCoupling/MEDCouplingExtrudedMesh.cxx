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
#include "MEDCouplingExtrudedMesh.hxx"
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingMemArray.hxx"
#include "CellModel.hxx"

#include <limits>
#include <algorithm>
#include <functional>
#include <iterator>
#include <cmath>
#include <list>
#include <set>

using namespace ParaMEDMEM;

/*!
 * Build an extruded mesh instance from 3D and 2D unstructured mesh lying on the \b same \b coords.
 * @param mesh3D 3D unstructured mesh.
 * @param mesh2D 2D unstructured mesh lying on the same coordinates than mesh3D. \b Warning mesh2D is \b not \b const
 * because the mesh is aggregated and potentially modified by rotate or translate method.
 * @param cell2DId Id of cell in mesh2D mesh where the computation of 1D mesh will be done.
 */
MEDCouplingExtrudedMesh *MEDCouplingExtrudedMesh::New(const MEDCouplingUMesh *mesh3D, MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception)
{
  return new MEDCouplingExtrudedMesh(mesh3D,mesh2D,cell2DId);
}

MEDCouplingMeshType MEDCouplingExtrudedMesh::getType() const
{
  return EXTRUDED;
}

MEDCouplingExtrudedMesh::MEDCouplingExtrudedMesh(const MEDCouplingUMesh *mesh3D, MEDCouplingUMesh *mesh2D, int cell2DId) throw(INTERP_KERNEL::Exception)
try:_mesh2D(mesh2D),_mesh1D(MEDCouplingUMesh::New()),_mesh3D_ids(0),_cell_2D_id(cell2DId)
{
  if(_mesh2D!=0)
    _mesh2D->incrRef();
  computeExtrusion(mesh3D);
}
catch(INTERP_KERNEL::Exception&)
  {
  }

bool MEDCouplingExtrudedMesh::isStructured() const
{
  return false;
}

int MEDCouplingExtrudedMesh::getNumberOfCells() const
{
  return _mesh2D->getNumberOfCells()*_mesh1D->getNumberOfCells();
}

int MEDCouplingExtrudedMesh::getNumberOfNodes() const
{
  return _mesh2D->getNumberOfNodes();
}

int MEDCouplingExtrudedMesh::getSpaceDimension() const
{
  return 3;
}

int MEDCouplingExtrudedMesh::getMeshDimension() const
{
  return 3;
}

INTERP_KERNEL::NormalizedCellType MEDCouplingExtrudedMesh::getTypeOfCell(int cellId) const
{
  int nbOfCells2D=_mesh2D->getNumberOfCells();
  int locId=cellId%nbOfCells2D;
  INTERP_KERNEL::NormalizedCellType tmp=_mesh2D->getTypeOfCell(locId);
  return INTERP_KERNEL::CellModel::getCellModel(tmp).getExtrudedType();
}

void MEDCouplingExtrudedMesh::getNodeIdsOfCell(int cellId, std::vector<int>& conn) const
{
  int nbOfCells2D=_mesh2D->getNumberOfCells();
  int nbOfNodes2D=_mesh2D->getNumberOfNodes();
  int locId=cellId%nbOfCells2D;
  int lev=cellId/nbOfCells2D;
  std::vector<int> tmp,tmp2;
  _mesh2D->getNodeIdsOfCell(locId,tmp);
  tmp2=tmp;
  std::transform(tmp.begin(),tmp.end(),tmp.begin(),std::bind2nd(std::plus<int>(),nbOfNodes2D*lev));
  std::transform(tmp2.begin(),tmp2.end(),tmp2.begin(),std::bind2nd(std::plus<int>(),nbOfNodes2D*(lev+1)));
  conn.insert(conn.end(),tmp.begin(),tmp.end());
  conn.insert(conn.end(),tmp2.begin(),tmp2.end());
}

void MEDCouplingExtrudedMesh::getCoordinatesOfNode(int nodeId, std::vector<double>& coo) const
{
  int nbOfNodes2D=_mesh2D->getNumberOfNodes();
  int locId=nodeId%nbOfNodes2D;
  int lev=nodeId/nbOfNodes2D;
  std::vector<double> tmp,tmp2;
  _mesh2D->getCoordinatesOfNode(locId,tmp);
  tmp2=tmp;
  int spaceDim=_mesh1D->getSpaceDimension();
  const double *z=_mesh1D->getCoords()->getConstPointer();
  std::transform(tmp.begin(),tmp.end(),z+lev*spaceDim,tmp.begin(),std::plus<double>());
  std::transform(tmp2.begin(),tmp2.end(),z+(lev+1)*spaceDim,tmp2.begin(),std::plus<double>());
  coo.insert(coo.end(),tmp.begin(),tmp.end());
  coo.insert(coo.end(),tmp2.begin(),tmp2.end());
}

void MEDCouplingExtrudedMesh::checkCoherency() const throw (INTERP_KERNEL::Exception)
{
}

void MEDCouplingExtrudedMesh::getBoundingBox(double *bbox) const
{
  double bbox2D[6];
  _mesh2D->getBoundingBox(bbox2D);
  const double *nodes1D=_mesh1D->getCoords()->getConstPointer();
  int nbOfNodes1D=_mesh1D->getNumberOfNodes();
  double bbox1DMin[3],bbox1DMax[3],tmp[3];
  std::fill(bbox1DMin,bbox1DMin+3,std::numeric_limits<double>::max());
  std::fill(bbox1DMax,bbox1DMax+3,-(std::numeric_limits<double>::max()));
  for(int i=0;i<nbOfNodes1D;i++)
    {
      std::transform(nodes1D+3*i,nodes1D+3*(i+1),bbox1DMin,bbox1DMin,static_cast<const double& (*)(const double&, const double&)>(std::min<double>));
      std::transform(nodes1D+3*i,nodes1D+3*(i+1),bbox1DMax,bbox1DMax,static_cast<const double& (*)(const double&, const double&)>(std::max<double>));
    }
  std::transform(bbox1DMax,bbox1DMax+3,bbox1DMin,tmp,std::minus<double>());
  int id=std::max_element(tmp,tmp+3)-tmp;
  bbox[0]=bbox1DMin[0]; bbox[1]=bbox1DMax[0];
  bbox[2]=bbox1DMin[1]; bbox[3]=bbox1DMax[1];
  bbox[4]=bbox1DMin[2]; bbox[5]=bbox1DMax[2];
  bbox[2*id+1]+=tmp[id];
}

void MEDCouplingExtrudedMesh::updateTime()
{
  if(_mesh2D)
    {
      updateTimeWith(*_mesh2D);
    }
  if(_mesh1D)
    {
      updateTimeWith(*_mesh1D);
    }
}

MEDCouplingFieldDouble *MEDCouplingExtrudedMesh::getMeasureField(bool) const
{
  //not implemented yet
  return 0;
}

MEDCouplingFieldDouble *MEDCouplingExtrudedMesh::getMeasureFieldOnNode(bool) const
{
  //not implemented yet
  return 0;
}

MEDCouplingFieldDouble *MEDCouplingExtrudedMesh::buildOrthogonalField() const
{
  //not implemented yet
  throw INTERP_KERNEL::Exception("MEDCouplingExtrudedMesh::buildOrthogonalField not implemented yet !");
}

int MEDCouplingExtrudedMesh::getCellContainingPoint(const double *pos, double eps) const
{
  //not implemented yet
  return -1;
}

MEDCouplingExtrudedMesh::~MEDCouplingExtrudedMesh()
{
  if(_mesh2D)
    _mesh2D->decrRef();
  if(_mesh1D)
    _mesh1D->decrRef();
  if(_mesh3D_ids)
    _mesh3D_ids->decrRef();
}

void MEDCouplingExtrudedMesh::computeExtrusion(const MEDCouplingUMesh *mesh3D) throw(INTERP_KERNEL::Exception)
{
  const char errMsg1[]="2D mesh is empty unable to compute extrusion !";
  const char errMsg2[]="Coords between 2D and 3D meshes are not the same ! Try MEDCouplingPointSet::tryToShareSameCoords method";
  const char errMsg3[]="No chance to find extrusion pattern in mesh3D,mesh2D couple because nbCells3D%nbCells2D!=0 !";
  if(_mesh2D==0 || mesh3D==0)
    throw INTERP_KERNEL::Exception(errMsg1);
  if(_mesh2D->getCoords()!=mesh3D->getCoords())
    throw INTERP_KERNEL::Exception(errMsg2);
  if(mesh3D->getNumberOfCells()%_mesh2D->getNumberOfCells()!=0)
    throw INTERP_KERNEL::Exception(errMsg3);
  if(!_mesh3D_ids)
    _mesh3D_ids=DataArrayInt::New();
  if(!_mesh1D)
    _mesh1D=MEDCouplingUMesh::New();
  computeExtrusionAlg(mesh3D);
}

void MEDCouplingExtrudedMesh::build1DExtrusion(int idIn3DDesc, int newId, int nbOf1DLev, MEDCouplingUMesh *subMesh,
                                               const int *desc3D, const int *descIndx3D,
                                               const int *revDesc3D, const int *revDescIndx3D,
                                               bool computeMesh1D) throw(INTERP_KERNEL::Exception)
{
  int nbOf2DCells=_mesh2D->getNumberOfCells();
  int start=revDescIndx3D[idIn3DDesc];
  int end=revDescIndx3D[idIn3DDesc+1];
  if(end-start!=1)
    {
      std::ostringstream ost; ost << "Invalid bases 2D mesh specified : 2D cell # " <<  idIn3DDesc;
      ost << " shared by more than 1 3D cell !!!";
      throw INTERP_KERNEL::Exception(ost.str().c_str());
    }
  int current3DCell=revDesc3D[start];
  int current2DCell=idIn3DDesc;
  int *mesh3DIDs=_mesh3D_ids->getPointer();
  mesh3DIDs[newId]=current3DCell;
  const int *conn2D=subMesh->getNodalConnectivity()->getConstPointer();
  const int *conn2DIndx=subMesh->getNodalConnectivityIndex()->getConstPointer();
  for(int i=1;i<nbOf1DLev;i++)
    {
      std::vector<int> conn(conn2D+conn2DIndx[current2DCell]+1,conn2D+conn2DIndx[current2DCell+1]);
      std::sort(conn.begin(),conn.end());
      if(computeMesh1D)
        computeBaryCenterOfFace(conn,i-1);
      current2DCell=findOppositeFaceOf(current2DCell,current3DCell,conn,
                                       desc3D,descIndx3D,conn2D,conn2DIndx);
      start=revDescIndx3D[current2DCell];
      end=revDescIndx3D[current2DCell+1];
      if(end-start!=2)
        {
          std::ostringstream ost; ost << "Expecting to have 2 3D cells attached to 2D cell " << current2DCell << "!";
          ost << " : Impossible or call tryToShareSameCoords method !";
          throw INTERP_KERNEL::Exception(ost.str().c_str());
        }
      if(revDesc3D[start]!=current3DCell)
        current3DCell=revDesc3D[start];
      else
        current3DCell=revDesc3D[start+1];
      mesh3DIDs[i*nbOf2DCells+newId]=current3DCell;
    }
  if(computeMesh1D)
    {
      std::vector<int> conn(conn2D+conn2DIndx[current2DCell]+1,conn2D+conn2DIndx[current2DCell+1]);
      std::sort(conn.begin(),conn.end());
      computeBaryCenterOfFace(conn,nbOf1DLev-1);
      current2DCell=findOppositeFaceOf(current2DCell,current3DCell,conn,
                                       desc3D,descIndx3D,conn2D,conn2DIndx);
      conn.clear();
      conn.insert(conn.end(),conn2D+conn2DIndx[current2DCell]+1,conn2D+conn2DIndx[current2DCell+1]);
      std::sort(conn.begin(),conn.end());
      computeBaryCenterOfFace(conn,nbOf1DLev);
    }
}

int MEDCouplingExtrudedMesh::findOppositeFaceOf(int current2DCell, int current3DCell, const std::vector<int>& connSorted,
                                                const int *desc3D, const int *descIndx3D,
                                                const int *conn2D, const int *conn2DIndx) throw(INTERP_KERNEL::Exception)
{
  int start=descIndx3D[current3DCell];
  int end=descIndx3D[current3DCell+1];
  bool found=false;
  for(const int *candidate2D=desc3D+start;candidate2D!=desc3D+end && !found;candidate2D++)
    {
      if(*candidate2D!=current2DCell)
        {
          std::vector<int> conn2(conn2D+conn2DIndx[*candidate2D]+1,conn2D+conn2DIndx[*candidate2D+1]);
          std::sort(conn2.begin(),conn2.end());
          std::list<int> intersect;
          std::set_intersection(connSorted.begin(),connSorted.end(),conn2.begin(),conn2.end(),
                                std::insert_iterator< std::list<int> >(intersect,intersect.begin()));
          if(intersect.empty())
            return *candidate2D;
        }
    }
  std::ostringstream ost; ost << "Impossible to find an opposite 2D face of face # " <<  current2DCell;
  ost << " in 3D cell # " << current3DCell << " : Impossible or call tryToShareSameCoords method !";
  throw INTERP_KERNEL::Exception(ost.str().c_str());
}

void MEDCouplingExtrudedMesh::computeBaryCenterOfFace(const std::vector<int>& nodalConnec, int lev1DId)
{
  double *zoneToUpdate=_mesh1D->getCoords()->getPointer()+lev1DId*3;
  std::fill(zoneToUpdate,zoneToUpdate+3,0.);
  const double *coords=_mesh2D->getCoords()->getConstPointer();
  for(std::vector<int>::const_iterator iter=nodalConnec.begin();iter!=nodalConnec.end();iter++)
    std::transform(zoneToUpdate,zoneToUpdate+3,coords+3*(*iter),zoneToUpdate,std::plus<double>());
  std::transform(zoneToUpdate,zoneToUpdate+3,zoneToUpdate,std::bind2nd(std::multiplies<double>(),(double)(1./nodalConnec.size())));
}

int MEDCouplingExtrudedMesh::findCorrespCellByNodalConn(const std::vector<int>& nodalConnec, const int *revNodalPtr, const int *revNodalIndxPtr) throw(INTERP_KERNEL::Exception)
{
  std::vector<int>::const_iterator iter=nodalConnec.begin();
  std::set<int> s1(revNodalPtr+revNodalIndxPtr[*iter],revNodalPtr+revNodalIndxPtr[*iter+1]);
  iter++;
  for(;iter!=nodalConnec.end();iter++)
    {
      std::set<int> s2(revNodalPtr+revNodalIndxPtr[*iter],revNodalPtr+revNodalIndxPtr[*iter+1]);
      std::set<int> s3;
      std::set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(),std::insert_iterator< std::set<int> >(s3,s3.end()));
      s1=s3;
    }
  if(s1.size()==1)
    return *(s1.begin());
  std::ostringstream ostr;
  ostr << "Cell with nodal connec : ";
  std::copy(nodalConnec.begin(),nodalConnec.end(),std::ostream_iterator<int>(ostr," "));
  ostr << " is not part of mesh";
  throw INTERP_KERNEL::Exception(ostr.str().c_str());
}

/*!
 * This method is callable on 1Dmeshes (meshDim==1 && spaceDim==3) returned by MEDCouplingExtrudedMesh::getMesh1D typically.
 * These 1Dmeshes (meshDim==1 && spaceDim==3) have a special semantic because these meshes do not specify a static location but a translation along a path.
 * This method checks that 'm1' and 'm2' are compatible, if not an exception is thrown. In case these meshes ('m1' and 'm2') are compatible 2 corresponding meshes
 * are created ('m1r' and 'm2r') that can be used for interpolation.
 * @param m1 input mesh with meshDim==1 and spaceDim==3
 * @param m2 input mesh with meshDim==1 and spaceDim==3
 * @param eps tolerance acceptable to determine compatibility
 * @param m1r output mesh with ref count equal to 1 with meshDim==1 and spaceDim==1
 * @param m2r output mesh with ref count equal to 1 with meshDim==1 and spaceDim==1
 * @param v is the output normalized vector of the common direction of 'm1' and 'm2'  
 * @throw in case that m1 and m2 are not compatible each other.
 */
void MEDCouplingExtrudedMesh::project1DMeshes(const MEDCouplingUMesh *m1, const MEDCouplingUMesh *m2, double eps,
                                              MEDCouplingUMesh *&m1r, MEDCouplingUMesh *&m2r, double *v) throw(INTERP_KERNEL::Exception)
{
  m1r=m1->clone(true);
  m2r=m2->clone(true);
  m1r->changeSpaceDimension(1);
  m2r->changeSpaceDimension(1);
}

void MEDCouplingExtrudedMesh::rotate(const double *center, const double *vector, double angle)
{
  _mesh2D->rotate(center,vector,angle);
  _mesh1D->rotate(center,vector,angle);
}

void MEDCouplingExtrudedMesh::translate(const double *vector)
{
  _mesh2D->translate(vector);
  _mesh1D->translate(vector);
}

MEDCouplingMesh *MEDCouplingExtrudedMesh::mergeMyselfWith(const MEDCouplingMesh *other) const
{
  // not implemented yet !
  return 0;
}

DataArrayDouble *MEDCouplingExtrudedMesh::getCoordinatesAndOwner() const
{
  DataArrayDouble *arr2D=_mesh2D->getCoords();
  DataArrayDouble *arr1D=_mesh1D->getCoords();
  DataArrayDouble *ret=DataArrayDouble::New();
  ret->alloc(getNumberOfNodes(),3);
  int nbOf1DLev=_mesh1D->getNumberOfNodes();
  int nbOf2DNodes=_mesh2D->getNumberOfNodes();
  const double *ptSrc=arr2D->getConstPointer();
  double *pt=ret->getPointer();
  std::copy(ptSrc,ptSrc+3*nbOf2DNodes,pt);
  for(int i=1;i<nbOf1DLev;i++)
    {
      std::copy(ptSrc,ptSrc+3*nbOf2DNodes,pt+3*i*nbOf2DNodes);
      double vec[3];
      std::copy(arr1D->getConstPointer()+3*i,arr1D->getConstPointer()+3*(i+1),vec);
      std::transform(arr1D->getConstPointer()+3*(i-1),arr1D->getConstPointer()+3*i,vec,vec,std::minus<double>());
      for(int j=0;j<nbOf2DNodes;j++)
        std::transform(vec,vec+3,pt+3*(i*nbOf2DNodes+j),pt+3*(i*nbOf2DNodes+j),std::plus<double>());
    }
  return ret;
}

DataArrayDouble *MEDCouplingExtrudedMesh::getBarycenterAndOwner() const
{
  //not yet implemented
  return 0;
}

void MEDCouplingExtrudedMesh::computeExtrusionAlg(const MEDCouplingUMesh *mesh3D) throw(INTERP_KERNEL::Exception)
{
  _mesh3D_ids->alloc(mesh3D->getNumberOfCells(),1);
  int nbOf1DLev=mesh3D->getNumberOfCells()/_mesh2D->getNumberOfCells();
  _mesh1D->setMeshDimension(1);
  _mesh1D->allocateCells(nbOf1DLev);
  int tmpConn[2];
  for(int i=0;i<nbOf1DLev;i++)
    {
      tmpConn[0]=i;
      tmpConn[1]=i+1;
      _mesh1D->insertNextCell(INTERP_KERNEL::NORM_SEG2,2,tmpConn);
    }
  _mesh1D->finishInsertingCells();
  DataArrayDouble *myCoords=DataArrayDouble::New();
  myCoords->alloc(nbOf1DLev+1,3);
  _mesh1D->setCoords(myCoords);
  myCoords->decrRef();
  DataArrayInt *desc,*descIndx,*revDesc,*revDescIndx;
  desc=DataArrayInt::New(); descIndx=DataArrayInt::New(); revDesc=DataArrayInt::New(); revDescIndx=DataArrayInt::New();
  MEDCouplingUMesh *subMesh=mesh3D->buildDescendingConnectivity(desc,descIndx,revDesc,revDescIndx);
  DataArrayInt *revNodal2D,*revNodalIndx2D;
  revNodal2D=DataArrayInt::New(); revNodalIndx2D=DataArrayInt::New();
  subMesh->getReverseNodalConnectivity(revNodal2D,revNodalIndx2D);
  const int *nodal2D=_mesh2D->getNodalConnectivity()->getConstPointer();
  const int *nodal2DIndx=_mesh2D->getNodalConnectivityIndex()->getConstPointer();
  const int *revNodal2DPtr=revNodal2D->getConstPointer();
  const int *revNodalIndx2DPtr=revNodalIndx2D->getConstPointer();
  const int *descP=desc->getConstPointer();
  const int *descIndxP=descIndx->getConstPointer();
  const int *revDescP=revDesc->getConstPointer();
  const int *revDescIndxP=revDescIndx->getConstPointer();
  //
  int nbOf2DCells=_mesh2D->getNumberOfCells();
  for(int i=0;i<nbOf2DCells;i++)
    {
      int idInSubMesh;
       std::vector<int> nodalConnec(nodal2D+nodal2DIndx[i]+1,nodal2D+nodal2DIndx[i+1]);
       try
        {
          idInSubMesh=findCorrespCellByNodalConn(nodalConnec,revNodal2DPtr,revNodalIndx2DPtr);
        }
       catch(INTERP_KERNEL::Exception& e)
         {
           std::ostringstream ostr; ostr << "mesh2D cell # " << i << " is not part of any cell of 3D mesh !\n";
           ostr << e.what();
           throw INTERP_KERNEL::Exception(ostr.str().c_str());
         }
       build1DExtrusion(idInSubMesh,i,nbOf1DLev,subMesh,descP,descIndxP,revDescP,revDescIndxP,i==_cell_2D_id);
    }
  //
  revNodal2D->decrRef();
  revNodalIndx2D->decrRef();
  subMesh->decrRef();
  desc->decrRef();
  descIndx->decrRef();
  revDesc->decrRef();
  revDescIndx->decrRef();
}
