// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingStructuredMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingUMesh.hxx"

#include <numeric>

using namespace ParaMEDMEM;

MEDCouplingStructuredMesh::MEDCouplingStructuredMesh()
{
}

MEDCouplingStructuredMesh::MEDCouplingStructuredMesh(const MEDCouplingStructuredMesh& other, bool deepCopy):MEDCouplingMesh(other)
{
}

MEDCouplingStructuredMesh::~MEDCouplingStructuredMesh()
{
}

std::size_t MEDCouplingStructuredMesh::getHeapMemorySize() const
{
  return MEDCouplingMesh::getHeapMemorySize();
}

void MEDCouplingStructuredMesh::copyTinyStringsFrom(const MEDCouplingMesh *other) throw(INTERP_KERNEL::Exception)
{
  MEDCouplingMesh::copyTinyStringsFrom(other);
}

bool MEDCouplingStructuredMesh::isEqualIfNotWhy(const MEDCouplingMesh *other, double prec, std::string& reason) const throw(INTERP_KERNEL::Exception)
{
  return MEDCouplingMesh::isEqualIfNotWhy(other,prec,reason);
}

INTERP_KERNEL::NormalizedCellType MEDCouplingStructuredMesh::getTypeOfCell(int cellId) const
{
  switch(getMeshDimension())
    {
    case 3:
      return INTERP_KERNEL::NORM_HEXA8;
    case 2:
      return INTERP_KERNEL::NORM_QUAD4;
    case 1:
      return INTERP_KERNEL::NORM_SEG2;
    default:
      throw INTERP_KERNEL::Exception("Unexpected dimension for MEDCouplingCurveLinearMesh::getTypeOfCell !");
    }
}

std::set<INTERP_KERNEL::NormalizedCellType> MEDCouplingStructuredMesh::getAllGeoTypes() const
{
  std::set<INTERP_KERNEL::NormalizedCellType> ret2;
  ret2.insert(getTypeOfCell(0));
  return ret2;
}

int MEDCouplingStructuredMesh::getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType type) const
{
  int ret=getNumberOfCells();
  if(type==getTypeOfCell(0))
    return ret;
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getTypeOfCell(0));
  std::ostringstream oss; oss << "MEDCouplingStructuredMesh::getNumberOfCellsWithType : no specified type ! Type available is " << cm.getRepr() << " !";
  throw INTERP_KERNEL::Exception(oss.str().c_str());
}

DataArrayInt *MEDCouplingStructuredMesh::giveCellsWithType(INTERP_KERNEL::NormalizedCellType type) const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  if(getTypeOfCell(0)==type)
    {
      ret->alloc(getNumberOfCells(),1);
      ret->iota(0);
    }
  else
    ret->alloc(0,1);
  return ret.retn();
}

DataArrayInt *MEDCouplingStructuredMesh::computeNbOfNodesPerCell() const throw(INTERP_KERNEL::Exception)
{
  int nbCells=getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbCells,1);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getTypeOfCell(0));
  ret->fillWithValue((int)cm.getNumberOfNodes());
  return ret.retn();
}

DataArrayInt *MEDCouplingStructuredMesh::computeNbOfFacesPerCell() const throw(INTERP_KERNEL::Exception)
{
  int nbCells=getNumberOfCells();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New();
  ret->alloc(nbCells,1);
  const INTERP_KERNEL::CellModel& cm=INTERP_KERNEL::CellModel::GetCellModel(getTypeOfCell(0));
  ret->fillWithValue((int)cm.getNumberOfSons());
  return ret.retn();
}

void MEDCouplingStructuredMesh::getNodeIdsOfCell(int cellId, std::vector<int>& conn) const
{
  int meshDim=getMeshDimension();
  int tmpCell[3],tmpNode[3];
  getSplitCellValues(tmpCell);
  getSplitNodeValues(tmpNode);
  int tmp2[3];
  GetPosFromId(cellId,meshDim,tmpCell,tmp2);
  switch(meshDim)
    {
    case 1:
      conn.push_back(tmp2[0]); conn.push_back(tmp2[0]+1);
      break;
    case 2:
      conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]); conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+1);
      conn.push_back((tmp2[1]+1)*(tmpCell[1]+1)+tmp2[0]+1); conn.push_back((tmp2[1]+1)*(tmpCell[1]+1)+tmp2[0]);
      break;
    case 3:
      conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+tmp2[2]*tmpNode[2]); conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+1+tmp2[2]*tmpNode[2]);
      conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+1+tmp2[2]*tmpNode[2]); conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+tmp2[2]*tmpNode[2]);
      conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+(tmp2[2]+1)*tmpNode[2]); conn.push_back(tmp2[1]*tmpCell[1]+tmp2[0]+1+(tmp2[2]+1)*tmpNode[2]);
      conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+1+(tmp2[2]+1)*tmpNode[2]); conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+(tmp2[2]+1)*tmpNode[2]);
      break;
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::getNodeIdsOfCell : big problem spacedim must be in 1,2 or 3 !");
    };
}

/*!
 * See MEDCouplingUMesh::getDistributionOfTypes for more information
 */
std::vector<int> MEDCouplingStructuredMesh::getDistributionOfTypes() const throw(INTERP_KERNEL::Exception)
{
  //only one type of cell
  std::vector<int> ret(3);
  ret[0]=getTypeOfCell(0);
  ret[1]=getNumberOfCells();
  ret[2]=0; //ret[3*k+2]==0 because it has no sense here
  return ret;
}

/*!
 * See MEDCouplingUMesh::checkTypeConsistencyAndContig for more information
 */
DataArrayInt *MEDCouplingStructuredMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  if(code.empty())
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkTypeConsistencyAndContig : code is empty, should not !");
  std::size_t sz=code.size();
  if(sz!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkTypeConsistencyAndContig : code should be of size 3 exactly !");

  int nbCells=getNumberOfCellsWithType((INTERP_KERNEL::NormalizedCellType)code[0]);
  if(code[2]==-1)
    {
      if(code[1]==nbCells)
        return 0;
      else
        throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkTypeConsistencyAndContig : number of cells mismatch !");
    }
  else
    {
      if(code[2]<-1) 
        throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkTypeConsistencyAndContig : code[2]<-1 mismatch !");
      if(code[2]>=(int)idsPerType.size()) 
        throw INTERP_KERNEL::Exception("MEDCouplingCurveLinearMesh::checkTypeConsistencyAndContig : code[2]>size idsPerType !");
      return idsPerType[code[2]]->deepCpy();
    }
}

/*!
 * See MEDCouplingUMesh::splitProfilePerType for more information
 */
void MEDCouplingStructuredMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  int nbCells=getNumberOfCells();
  code.resize(3);
  code[0]=(int)getTypeOfCell(0);
  code[1]=nbCells;
  code[2]=0;
  idsInPflPerType.push_back(profile->deepCpy());
  idsPerType.push_back(profile->deepCpy());
}

/*!
 * Creates a new unstructured mesh (MEDCouplingUMesh) from \a this structured one.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 * delete this array using decrRef() as it is no more needed. 
 *  \throw If \a this->getMeshDimension() is not among [1,2,3].
 */
MEDCouplingUMesh *MEDCouplingStructuredMesh::buildUnstructured() const throw(INTERP_KERNEL::Exception)
{
  int meshDim=getMeshDimension();
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New(getName(),meshDim);
  DataArrayDouble *coords=getCoordinatesAndOwner();
  ret->setCoords(coords);
  coords->decrRef();
  switch(meshDim)
    {
    case 1:
      fill1DUnstructuredMesh(ret);
      break;
    case 2:
      fill2DUnstructuredMesh(ret);
      break;
    case 3:
      fill3DUnstructuredMesh(ret);
      break;
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::buildUnstructured : big problem spacedim must be in 1,2 or 3 !");
    };
  return ret;
}

/*!
 * Creates a new MEDCouplingUMesh containing a part of cells of \a this mesh.
 * The cells to include to the
 * result mesh are specified by an array of cell ids.
 *  \param [in] start - an array of cell ids to include to the result mesh.
 *  \param [in] end - specifies the end of the array \a start, so that
 *              the last value of \a start is \a end[ -1 ].
 *  \return MEDCouplingMesh * - a new instance of MEDCouplingUMesh. The caller is to
 *         delete this mesh using decrRef() as it is no more needed. 
 */
MEDCouplingMesh *MEDCouplingStructuredMesh::buildPart(const int *start, const int *end) const
{
  MEDCouplingUMesh *um=buildUnstructured();
  MEDCouplingMesh *ret=um->buildPart(start,end);
  um->decrRef();
  return ret;
}

MEDCouplingMesh *MEDCouplingStructuredMesh::buildPartAndReduceNodes(const int *start, const int *end, DataArrayInt*& arr) const
{
  MEDCouplingUMesh *um=buildUnstructured();
  MEDCouplingMesh *ret=um->buildPartAndReduceNodes(start,end,arr);
  um->decrRef();
  return ret;
}

DataArrayInt *MEDCouplingStructuredMesh::simplexize(int policy) throw(INTERP_KERNEL::Exception)
{
  throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::simplexize : not available for Cartesian mesh !");
}

/*!
 * Returns a new MEDCouplingFieldDouble holding normal vectors to cells of \a this
 * 2D mesh. The computed vectors have 3 components and are normalized.
 *  \return MEDCouplingFieldDouble * - a new instance of MEDCouplingFieldDouble on
 *          cells and one time. The caller is to delete this field using decrRef() as
 *          it is no more needed.
 *  \throw If \a this->getMeshDimension() != 2.
 */
MEDCouplingFieldDouble *MEDCouplingStructuredMesh::buildOrthogonalField() const
{
  if(getMeshDimension()!=2)
    throw INTERP_KERNEL::Exception("Expected a MEDCouplingStructuredMesh with meshDim == 2 !");
  MEDCouplingFieldDouble *ret=MEDCouplingFieldDouble::New(ON_CELLS,NO_TIME);
  DataArrayDouble *array=DataArrayDouble::New();
  int nbOfCells=getNumberOfCells();
  array->alloc(nbOfCells,3);
  double *vals=array->getPointer();
  for(int i=0;i<nbOfCells;i++)
    { vals[3*i]=0.; vals[3*i+1]=0.; vals[3*i+2]=1.; }
  ret->setArray(array);
  array->decrRef();
  ret->setMesh(this);
  return ret;
}

void MEDCouplingStructuredMesh::fill1DUnstructuredMesh(MEDCouplingUMesh *m) const
{
  int nbOfCells=-1;
  getNodeGridStructure(&nbOfCells);
  nbOfCells--;
  DataArrayInt *connI=DataArrayInt::New();
  connI->alloc(nbOfCells+1,1);
  int *ci=connI->getPointer();
  DataArrayInt *conn=DataArrayInt::New();
  conn->alloc(3*nbOfCells,1);
  ci[0]=0;
  int *cp=conn->getPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      cp[3*i]=(int)INTERP_KERNEL::NORM_SEG2;
      cp[3*i+1]=i;
      cp[3*i+2]=i+1;
      ci[i+1]=3*(i+1);
    }
  m->setConnectivity(conn,connI,true);
  conn->decrRef();
  connI->decrRef();
}

void MEDCouplingStructuredMesh::fill2DUnstructuredMesh(MEDCouplingUMesh *m) const
{
  int ns[2];
  getNodeGridStructure(ns);
  int n1=ns[0]-1;
  int n2=ns[1]-1;
  DataArrayInt *connI=DataArrayInt::New();
  connI->alloc(n1*n2+1,1);
  int *ci=connI->getPointer();
  DataArrayInt *conn=DataArrayInt::New();
  conn->alloc(5*n1*n2,1);
  ci[0]=0;
  int *cp=conn->getPointer();
  int pos=0;
  for(int j=0;j<n2;j++)
    for(int i=0;i<n1;i++,pos++)
      {
        cp[5*pos]=(int)INTERP_KERNEL::NORM_QUAD4;
        cp[5*pos+1]=i+1+j*(n1+1);
        cp[5*pos+2]=i+j*(n1+1);
        cp[5*pos+3]=i+(j+1)*(n1+1);
        cp[5*pos+4]=i+1+(j+1)*(n1+1);
        ci[pos+1]=5*(pos+1);
    }
  m->setConnectivity(conn,connI,true);
  conn->decrRef();
  connI->decrRef();
}

void MEDCouplingStructuredMesh::fill3DUnstructuredMesh(MEDCouplingUMesh *m) const
{
  int ns[3];
  getNodeGridStructure(ns);
  int n1=ns[0]-1;
  int n2=ns[1]-1;
  int n3=ns[2]-1;
  DataArrayInt *connI=DataArrayInt::New();
  connI->alloc(n1*n2*n3+1,1);
  int *ci=connI->getPointer();
  DataArrayInt *conn=DataArrayInt::New();
  conn->alloc(9*n1*n2*n3,1);
  ci[0]=0;
  int *cp=conn->getPointer();
  int pos=0;
  for(int k=0;k<n3;k++)
    for(int j=0;j<n2;j++)
      for(int i=0;i<n1;i++,pos++)
        {
          cp[9*pos]=(int)INTERP_KERNEL::NORM_HEXA8;
          int tmp=(n1+1)*(n2+1);
          cp[9*pos+1]=i+1+j*(n1+1)+k*tmp;
          cp[9*pos+2]=i+j*(n1+1)+k*tmp;
          cp[9*pos+3]=i+(j+1)*(n1+1)+k*tmp;
          cp[9*pos+4]=i+1+(j+1)*(n1+1)+k*tmp;
          cp[9*pos+5]=i+1+j*(n1+1)+(k+1)*tmp;
          cp[9*pos+6]=i+j*(n1+1)+(k+1)*tmp;
          cp[9*pos+7]=i+(j+1)*(n1+1)+(k+1)*tmp;
          cp[9*pos+8]=i+1+(j+1)*(n1+1)+(k+1)*tmp;
          ci[pos+1]=9*(pos+1);
        }
  m->setConnectivity(conn,connI,true);
  conn->decrRef();
  connI->decrRef();
}

/*!
 * Returns a cell id by its (i,j,k) index. The cell is located between the i-th and
 * ( i + 1 )-th nodes along X axis etc.
 *  \param [in] i - a index of node coordinates array along X axis.
 *  \param [in] j - a index of node coordinates array along Y axis.
 *  \param [in] k - a index of node coordinates array along Z axis.
 *  \return int - a cell id in \a this mesh.
 */
int MEDCouplingStructuredMesh::getCellIdFromPos(int i, int j, int k) const
{
  int tmp[3]={i,j,k};
  int tmp2[3];
  int meshDim=getMeshDimension();
  getSplitCellValues(tmp2);
  std::transform(tmp,tmp+meshDim,tmp2,tmp,std::multiplies<int>());
  return std::accumulate(tmp,tmp+meshDim,0);
}

/*!
 * Returns a node id by its (i,j,k) index.
 *  \param [in] i - a index of node coordinates array along X axis.
 *  \param [in] j - a index of node coordinates array along Y axis.
 *  \param [in] k - a index of node coordinates array along Z axis.
 *  \return int - a node id in \a this mesh.
 */
int MEDCouplingStructuredMesh::getNodeIdFromPos(int i, int j, int k) const
{
  int tmp[3]={i,j,k};
  int tmp2[3];
  int meshDim=getMeshDimension();
  getSplitNodeValues(tmp2);
  std::transform(tmp,tmp+meshDim,tmp2,tmp,std::multiplies<int>());
  return std::accumulate(tmp,tmp+meshDim,0);
}

void MEDCouplingStructuredMesh::GetPosFromId(int nodeId, int meshDim, const int *split, int *res)
{
  int work=nodeId;
  for(int i=meshDim-1;i>=0;i--)
    {
      int pos=work/split[i];
      work=work%split[i];
      res[i]=pos;
    }
}
