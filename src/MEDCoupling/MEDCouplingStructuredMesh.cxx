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
#include "MEDCoupling1GTUMesh.hxx"
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
  return GetGeoTypeGivenMeshDimension(getMeshDimension());
}

INTERP_KERNEL::NormalizedCellType MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension(int meshDim) throw(INTERP_KERNEL::Exception)
{
  switch(meshDim)
    {
    case 3:
      return INTERP_KERNEL::NORM_HEXA8;
    case 2:
      return INTERP_KERNEL::NORM_QUAD4;
    case 1:
      return INTERP_KERNEL::NORM_SEG2;
    default:
      throw INTERP_KERNEL::Exception("Unexpected dimension for MEDCouplingStructuredMesh::GetGeoTypeGivenMeshDimension !");
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

/*!
 * This method computes effective number of nodes per cell. That is to say nodes appearing several times in nodal connectivity of a cell,
 * will be counted only once here whereas it will be counted several times in MEDCouplingMesh::computeNbOfNodesPerCell method.
 * Here for structured mesh it returns exactly as MEDCouplingStructuredMesh::computeNbOfNodesPerCell does.
 *
 * \return DataArrayInt * - new object to be deallocated by the caller.
 */
DataArrayInt *MEDCouplingStructuredMesh::computeEffectiveNbOfNodesPerCell() const throw(INTERP_KERNEL::Exception)
{
  return computeNbOfNodesPerCell();
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
      conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]); conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+1);
      conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+1); conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]);
      break;
    case 3:
      conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+tmp2[2]*tmpNode[2]); conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+1+tmp2[2]*tmpNode[2]);
      conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+1+tmp2[2]*tmpNode[2]); conn.push_back((tmp2[1]+1)*tmpNode[1]+tmp2[0]+tmp2[2]*tmpNode[2]);
      conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+(tmp2[2]+1)*tmpNode[2]); conn.push_back(tmp2[1]*tmpNode[1]+tmp2[0]+1+(tmp2[2]+1)*tmpNode[2]);
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
  ret[2]=-1; //ret[3*k+2]==-1 because it has no sense here
  return ret;
}

/*!
 * This method tries to minimize at most the number of deep copy.
 * So if \a idsPerType is not empty it can be returned directly (without copy, but with ref count incremented) in return.
 * 
 * See MEDCouplingUMesh::checkTypeConsistencyAndContig for more information
 */
DataArrayInt *MEDCouplingStructuredMesh::checkTypeConsistencyAndContig(const std::vector<int>& code, const std::vector<const DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  int nbOfCells=getNumberOfCells();
  if(code.size()!=3)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : invalid input code should be exactly of size 3 !");
  if(code[0]!=(int)getTypeOfCell(0))
    {
      std::ostringstream oss; oss << "MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : Mismatch of geometric type ! Asking for " << code[0] << " whereas the geometric type is \a this is " << getTypeOfCell(0) << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(code[2]==-1)
    {
      if(code[1]==nbOfCells)
        return 0;
      else
        {
          std::ostringstream oss; oss << "MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : mismatch between the number of cells in this (" << nbOfCells << ") and the number of non profile (" << code[1] << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(code[2]!=0)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : single geo type mesh ! 0 or -1 is expected at pos #2 of input code !");
  if(idsPerType.size()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : input code points to DataArrayInt #0 whereas the size of idsPerType is not equal to 1 !");
  const DataArrayInt *pfl=idsPerType[0];
  if(!pfl)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : the input code points to a NULL DataArrayInt at rank 0 !");
  if(pfl->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::checkTypeConsistencyAndContig : input profile should have exactly one component !");
  pfl->checkAllIdsInRange(0,nbOfCells);
  pfl->incrRef();
  return const_cast<DataArrayInt *>(pfl);
}

/*!
 * This method is the opposite of MEDCouplingUMesh::checkTypeConsistencyAndContig method. Given a list of cells in \a profile it returns a list of sub-profiles sorted by geo type.
 * The result is put in the array \a idsPerType. In the returned parameter \a code, foreach i \a code[3*i+2] refers (if different from -1) to a location into the \a idsPerType.
 * This method has 1 input \a profile and 3 outputs \a code \a idsInPflPerType and \a idsPerType.
 * 
 * \param [out] code is a vector of size 3*n where n is the number of different geometric type in \a this \b reduced to the profile \a profile. \a code has exactly the same semantic than in MEDCouplingUMesh::checkTypeConsistencyAndContig method.
 * \param [out] idsInPflPerType is a vector of size of different geometric type in the subpart defined by \a profile of \a this ( equal to \a code.size()/3). For each i,
 *              \a idsInPflPerType[i] stores the tuple ids in \a profile that correspond to the geometric type code[3*i+0]
 * \param [out] idsPerType is a vector of size of different sub profiles needed to be defined to represent the profile \a profile for a given geometric type.
 *              This vector can be empty in case of all geometric type cells are fully covered in ascending in the given input \a profile.
 * 
 * \warning for performance reasons no deep copy will be performed, if \a profile can been used as this in output parameters \a idsInPflPerType and \a idsPerType.
 *
 * \throw if \a profile has not exactly one component. It throws too, if \a profile contains some values not in [0,getNumberOfCells()) or if \a this is not fully defined
 *
 *  \b Example1: <br>
 *          - Before \a this has 3 cells \a profile contains [0,1,2]
 *          - After \a code contains [NORM_...,nbCells,-1], \a idsInPflPerType [[0,1,2]] and \a idsPerType is empty <br>
 * 
 *  \b Example2: <br>
 *          - Before \a this has 3 cells \a profile contains [1,2]
 *          - After \a code contains [NORM_...,nbCells,0], \a idsInPflPerType [[0,1]] and \a idsPerType is [[1,2]] <br>

 */
void MEDCouplingStructuredMesh::splitProfilePerType(const DataArrayInt *profile, std::vector<int>& code, std::vector<DataArrayInt *>& idsInPflPerType, std::vector<DataArrayInt *>& idsPerType) const throw(INTERP_KERNEL::Exception)
{
  if(!profile || !profile->isAllocated())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::splitProfilePerType : input profile is NULL or not allocated !");
  if(profile->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::splitProfilePerType : input profile should have exactly one component !");
  int nbTuples=profile->getNumberOfTuples();
  int nbOfCells=getNumberOfCells();
  code.resize(3); idsInPflPerType.resize(1);
  code[0]=(int)getTypeOfCell(0); code[1]=nbOfCells;
  idsInPflPerType.resize(1);
  if(profile->isIdentity() && nbTuples==nbOfCells)
    {
      code[2]=-1;
      idsInPflPerType[0]=0;
      idsPerType.clear();
      return ;
    }
  code[1]=profile->getNumberOfTuples();
  code[2]=0;
  profile->checkAllIdsInRange(0,nbOfCells);
  idsPerType.resize(1);
  idsPerType[0]=profile->deepCpy();
  idsInPflPerType[0]=DataArrayInt::Range(0,nbTuples,1);
}

/*!
 * Creates a new unstructured mesh (MEDCoupling1SGTUMesh) from \a this structured one.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 * delete this array using decrRef() as it is no more needed. 
 *  \throw If \a this->getMeshDimension() is not among [1,2,3].
 */
MEDCoupling1SGTUMesh *MEDCouplingStructuredMesh::build1SGTUnstructured() const throw(INTERP_KERNEL::Exception)
{
  int meshDim=getMeshDimension(); 
  if(meshDim<0 || meshDim>3)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::build1SGTUnstructured : meshdim must be in [1,2,3] !");
  MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> coords(getCoordinatesAndOwner());
  int ns[3];
  getNodeGridStructure(ns);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(Build1GTNodalConnectivity(ns,ns+meshDim));
  MEDCouplingAutoRefCountObjectPtr<MEDCoupling1SGTUMesh> ret(MEDCoupling1SGTUMesh::New(getName().c_str(),GetGeoTypeGivenMeshDimension(meshDim)));
  ret->setNodalConnectivity(conn); ret->setCoords(coords);
  return ret.retn();
}

/*!
 * Creates a new unstructured mesh (MEDCouplingUMesh) from \a this structured one.
 *  \return MEDCouplingUMesh * - a new instance of MEDCouplingUMesh. The caller is to
 * delete this array using decrRef() as it is no more needed. 
 *  \throw If \a this->getMeshDimension() is not among [1,2,3].
 */
MEDCouplingUMesh *MEDCouplingStructuredMesh::buildUnstructured() const throw(INTERP_KERNEL::Exception)
{
  MEDCouplingAutoRefCountObjectPtr<MEDCoupling1SGTUMesh> ret0(build1SGTUnstructured());
  return ret0->buildUnstructured();
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
  std::vector<int> cgs(getCellGridStructure());
  std::vector< std::pair<int,int> > cellPartFormat,nodePartFormat;
  if(IsPartStructured(start,end,cgs,cellPartFormat))
    {
      MEDCouplingAutoRefCountObjectPtr<MEDCouplingStructuredMesh> ret(buildStructuredSubPart(cellPartFormat));
      nodePartFormat=cellPartFormat;
      for(std::vector< std::pair<int,int> >::iterator it=nodePartFormat.begin();it!=nodePartFormat.end();it++)
        (*it).second++;
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp1(BuildExplicitIdsFrom(getNodeGridStructure(),nodePartFormat));
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp2(DataArrayInt::New()); tmp2->alloc(getNumberOfNodes(),1);
      tmp2->fillWithValue(-1);
      MEDCouplingAutoRefCountObjectPtr<DataArrayInt> tmp3(DataArrayInt::New()); tmp3->alloc(tmp1->getNumberOfTuples(),1); tmp3->iota(0);
      tmp2->setPartOfValues3(tmp3,tmp1->begin(),tmp1->end(),0,1,1);
      arr=tmp2.retn();
      return ret.retn();
    }
  else
    {
      MEDCouplingUMesh *um=buildUnstructured();
      MEDCouplingMesh *ret=um->buildPartAndReduceNodes(start,end,arr);
      um->decrRef();
      return ret;
    }
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

/*!
 * \return DataArrayInt * - newly allocated instance of nodal connectivity compatible for MEDCoupling1SGTMesh instance
 */
DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivity(const int *nodeStBg, const int *nodeStEnd) throw(INTERP_KERNEL::Exception)
{
  std::size_t dim=std::distance(nodeStBg,nodeStEnd);
  switch(dim)
    {
    case 1:
      return Build1GTNodalConnectivity1D(nodeStBg);
    case 2:
      return Build1GTNodalConnectivity2D(nodeStBg);
    case 3:
      return Build1GTNodalConnectivity3D(nodeStBg);
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::Build1GTNodalConnectivity : only dimension in [1,2,3] supported !");
    }
}

DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivity1D(const int *nodeStBg) throw(INTERP_KERNEL::Exception)
{
  int nbOfCells(*nodeStBg-1);
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(DataArrayInt::New());
  conn->alloc(2*nbOfCells,1);
  int *cp=conn->getPointer();
  for(int i=0;i<nbOfCells;i++)
    {
      cp[2*i+0]=i;
      cp[2*i+1]=i+1;
    }
  return conn.retn();
}

DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivity2D(const int *nodeStBg) throw(INTERP_KERNEL::Exception)
{
  int n1=nodeStBg[0]-1;
  int n2=nodeStBg[1]-1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(DataArrayInt::New());
  conn->alloc(4*n1*n2,1);
  int *cp=conn->getPointer();
  int pos=0;
  for(int j=0;j<n2;j++)
    for(int i=0;i<n1;i++,pos++)
      {
        cp[4*pos+0]=i+1+j*(n1+1);
        cp[4*pos+1]=i+j*(n1+1);
        cp[4*pos+2]=i+(j+1)*(n1+1);
        cp[4*pos+3]=i+1+(j+1)*(n1+1);
    }
  return conn.retn();
}

DataArrayInt *MEDCouplingStructuredMesh::Build1GTNodalConnectivity3D(const int *nodeStBg) throw(INTERP_KERNEL::Exception)
{
  int n1=nodeStBg[0]-1;
  int n2=nodeStBg[1]-1;
  int n3=nodeStBg[2]-1;
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> conn(DataArrayInt::New());
  conn->alloc(8*n1*n2*n3,1);
  int *cp=conn->getPointer();
  int pos=0;
  for(int k=0;k<n3;k++)
    for(int j=0;j<n2;j++)
      for(int i=0;i<n1;i++,pos++)
        {
          int tmp=(n1+1)*(n2+1);
          cp[8*pos+0]=i+1+j*(n1+1)+k*tmp;
          cp[8*pos+1]=i+j*(n1+1)+k*tmp;
          cp[8*pos+2]=i+(j+1)*(n1+1)+k*tmp;
          cp[8*pos+3]=i+1+(j+1)*(n1+1)+k*tmp;
          cp[8*pos+4]=i+1+j*(n1+1)+(k+1)*tmp;
          cp[8*pos+5]=i+j*(n1+1)+(k+1)*tmp;
          cp[8*pos+6]=i+(j+1)*(n1+1)+(k+1)*tmp;
          cp[8*pos+7]=i+1+(j+1)*(n1+1)+(k+1)*tmp;
        }
  return conn.retn();
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

std::vector<int> MEDCouplingStructuredMesh::getCellGridStructure() const throw(INTERP_KERNEL::Exception)
{
  std::vector<int> ret(getNodeGridStructure());
  std::transform(ret.begin(),ret.end(),ret.begin(),std::bind2nd(std::plus<int>(),-1));
  return ret;
}

/*!
 * This method states if given part ids [ \a startIds, \a stopIds) and a structure \a st returns if it can be considered as a structured dataset.
 * If true is returned \a partCompactFormat will contain the information to build the corresponding part.
 *
 * \sa MEDCouplingStructuredMesh::BuildExplicitIdsFrom
 */
bool MEDCouplingStructuredMesh::IsPartStructured(const int *startIds, const int *stopIds, const std::vector<int>& st, std::vector< std::pair<int,int> >& partCompactFormat) throw(INTERP_KERNEL::Exception)
{
  int dim((int)st.size());
  partCompactFormat.resize(dim);
  if(dim<1 || dim>3)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::isPartStructured : input structure must be of dimension in [1,2,3] !");
  std::vector<int> tmp2(dim),tmp(dim),tmp3(dim),tmp4(dim); tmp2[0]=1;
  for(int i=1;i<dim;i++)
    tmp2[i]=tmp2[i-1]*st[i-1];
  std::size_t sz(std::distance(startIds,stopIds));
  if(sz==0)
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::IsPartStructured : empty input !");
  GetPosFromId(*startIds,dim,&tmp2[0],&tmp[0]);
  partCompactFormat.resize(dim);
  for(int i=0;i<dim;i++)
    partCompactFormat[i].first=tmp[i];
  if(tmp[dim-1]<0 || tmp[dim-1]>=st[dim-1])
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::IsPartStructured : first id in input is not in valid range !");
  if(sz==1)
    {
      for(int i=0;i<dim;i++)
        partCompactFormat[i].second=tmp[i]+1;
      return true;
    }
  GetPosFromId(startIds[sz-1],dim,&tmp2[0],&tmp3[0]);
  int szExp(1);
  for(int i=0;i<dim;i++)
    {
      if(tmp3[i]<0 || tmp3[i]>=st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::IsPartStructured : last id in input is not in valid range !");
      partCompactFormat[i].second=tmp3[i]+1;
      tmp4[i]=partCompactFormat[i].second-partCompactFormat[i].first;
      if(tmp4[i]<=0)
        return false;
      szExp*=tmp4[i];
    }
  if(szExp!=(int)sz)
    return false;
  const int *w(startIds);
  switch(dim)
    {
    case 3:
      {
        for(int i=0;i<tmp4[2];i++)
          {
            int a=tmp2[2]*(partCompactFormat[2].first+i);
            for(int j=0;j<tmp4[1];j++)
              {
                int b=tmp2[1]*(partCompactFormat[1].first+j);
                for(int k=0;k<tmp4[0];k++,w++)
                  {
                    if(partCompactFormat[0].first+k+b+a!=*w)
                      return false;
                  }
              }
          }
        return true;
      }
    case 2:
      {
        for(int j=0;j<tmp4[1];j++)
          {
            int b=tmp2[1]*(partCompactFormat[1].first+j);
            for(int k=0;k<tmp4[0];k++,w++)
              {
                if(partCompactFormat[0].first+k+b!=*w)
                  return false;
              }
          }
        return true;
      }
    case 1:
      {
        for(int k=0;k<tmp4[0];k++,w++)
          {
            if(partCompactFormat[0].first+k!=*w)
              return false;
          }
        return true;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::IsPartStructured : internal error !");
    }
}

/*!
 * This method builds the explicit entity array from the structure in \a st and the range in \a partCompactFormat.
 *If the range contains invalid values regarding sructure an exception will be thrown.
 *
 * \return DataArrayInt * - a new object.
 * \sa MEDCouplingStructuredMesh::IsPartStructured
 */
DataArrayInt *MEDCouplingStructuredMesh::BuildExplicitIdsFrom(const std::vector<int>& st, const std::vector< std::pair<int,int> >& partCompactFormat) throw(INTERP_KERNEL::Exception)
{
  if(st.size()!=partCompactFormat.size())
    throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : input arrays must have the same size !");
  int nbOfItems(1);
  std::vector<int> dims(st.size());
  for(std::size_t i=0;i<st.size();i++)
    {
      if(partCompactFormat[i].first<0 || partCompactFormat[i].first>st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : invalid input range 1 !");
      if(partCompactFormat[i].second<0 || partCompactFormat[i].second>st[i])
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : invalid input range 2 !");
      if(partCompactFormat[i].second<=partCompactFormat[i].first)
        throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : invalid input range 3 !");
      dims[i]=partCompactFormat[i].second-partCompactFormat[i].first;
      nbOfItems*=dims[i];
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New());
  ret->alloc(nbOfItems,1);
  int *pt(ret->getPointer());
  switch(st.size())
    {
    case 3:
      {
        for(int i=0;i<dims[2];i++)
          {
            int a=(partCompactFormat[2].first+i)*st[0]*st[1];
            for(int j=0;j<dims[1];j++)
              {
                int b=(partCompactFormat[1].first+j)*st[0];
                for(int k=0;k<dims[0];k++,pt++)
                  *pt=partCompactFormat[0].first+k+b+a;
              }
          }
        break;
      }
    case 2:
      {
        for(int j=0;j<dims[1];j++)
          {
            int b=(partCompactFormat[1].first+j)*st[0];
            for(int k=0;k<dims[0];k++,pt++)
              *pt=partCompactFormat[0].first+k+b;
          }
        break;
      }
    case 1:
      {
        for(int k=0;k<dims[0];k++,pt++)
          *pt=partCompactFormat[0].first+k;
        break;
      }
    default:
      throw INTERP_KERNEL::Exception("MEDCouplingStructuredMesh::BuildExplicitIdsFrom : Dimension supported are 1,2 or 3 !");
    }
  return ret.retn();
}
