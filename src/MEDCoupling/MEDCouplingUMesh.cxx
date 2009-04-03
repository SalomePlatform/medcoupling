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
#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "CellModel.hxx"
#include "VolSurfFormulae.hxx"

#include <sstream>

using namespace ParaMEDMEM;

const char MEDCouplingUMesh::PART_OF_NAME[]="PartOf_";

MEDCouplingUMesh *MEDCouplingUMesh::New()
{
  return new MEDCouplingUMesh;
}

MEDCouplingUMesh *MEDCouplingUMesh::clone(bool recDeepCpy) const
{
  return new MEDCouplingUMesh(*this,recDeepCpy);
}

void MEDCouplingUMesh::updateTime()
{
  MEDCouplingPointSet::updateTime();
  if(_nodal_connec)
    {
      updateTimeWith(*_nodal_connec);
    }
  if(_nodal_connec_index)
    {
      updateTimeWith(*_nodal_connec_index);
    }
}

MEDCouplingUMesh::MEDCouplingUMesh():_iterator(-1),_mesh_dim(-1),
                                     _nodal_connec(0),_nodal_connec_index(0)
{
}

void MEDCouplingUMesh::checkCoherency() const throw(INTERP_KERNEL::Exception)
{
  for(std::set<INTERP_KERNEL::NormalizedCellType>::const_iterator iter=_types.begin();iter!=_types.end();iter++)
    {
      if(INTERP_KERNEL::CellModel::getCellModel(*iter).getDimension()!=_mesh_dim)
        {
          std::ostringstream message;
          message << "Mesh invalid because dimension is " << _mesh_dim << " and there is presence of cell(s) with type " << (*iter);
          throw INTERP_KERNEL::Exception(message.str().c_str());
        }
    }
}

void MEDCouplingUMesh::setMeshDimension(unsigned meshDim)
{
  _mesh_dim=meshDim;
  declareAsNew();
}

void MEDCouplingUMesh::allocateCells(int nbOfCells)
{
  if(_nodal_connec_index)
    {
      _nodal_connec_index->decrRef();
    }
  if(_nodal_connec)
    {
      _nodal_connec->decrRef();
    }

  _nodal_connec_index=DataArrayInt::New();
  _nodal_connec_index->alloc(nbOfCells+1,1);
  int *pt=_nodal_connec_index->getPointer();
  pt[0]=0;
  _nodal_connec=DataArrayInt::New();
  _nodal_connec->alloc(2*nbOfCells,1);
  _iterator=0;
  _types.clear();
  declareAsNew();
}

void MEDCouplingUMesh::insertNextCell(INTERP_KERNEL::NormalizedCellType type, int size, const int *nodalConnOfCell)
{
  int *pt=_nodal_connec_index->getPointer();
  int idx=pt[_iterator];

  _nodal_connec->writeOnPlace(idx,type,nodalConnOfCell,size);
  _types.insert(type);
  pt[++_iterator]=idx+size+1;
}

void MEDCouplingUMesh::finishInsertingCells()
{
  const int *pt=_nodal_connec_index->getConstPointer();
  int idx=pt[_iterator];

  _nodal_connec->reAlloc(idx);
  _nodal_connec_index->reAlloc(_iterator+1);
  _iterator=-1;
}

bool MEDCouplingUMesh::isEqual(const MEDCouplingMesh *other, double prec) const
{
  checkFullyDefined();
  const MEDCouplingUMesh *otherC=dynamic_cast<const MEDCouplingUMesh *>(other);
  if(!otherC)
    return false;
  otherC->checkFullyDefined();
  if(!MEDCouplingMesh::isEqual(other,prec))
    return false;
  if(_mesh_dim!=otherC->_mesh_dim)
    return false;
  if(_types!=otherC->_types)
    return false;
  if(!areCoordsEqual(*otherC,prec))
    return false;
  return true;
}

/*!
 * \b WARNING this method do the assumption that connectivity lies on the coordinates set.
 * For speed reasons no check of this will be done.
 */
void MEDCouplingUMesh::getReverseNodalConnectivity(DataArrayInt *revNodal, DataArrayInt *revNodalIndx) const
{
  checkFullyDefined();
  int nbOfNodes=getNumberOfNodes();
  int *revNodalIndxPtr=new int[nbOfNodes+1];
  revNodalIndx->useArray(revNodalIndxPtr,true,CPP_DEALLOC,nbOfNodes+1,1);
  std::fill(revNodalIndxPtr,revNodalIndxPtr+nbOfNodes+1,0);
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int nbOfCells=getNumberOfCells();
  int nbOfEltsInRevNodal=0;
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      const int *strtNdlConnOfCurCell=conn+connIndex[eltId]+1;
      const int *endNdlConnOfCurCell=conn+connIndex[eltId+1];
      for(const int *iter=strtNdlConnOfCurCell;iter!=endNdlConnOfCurCell;iter++)
        if(*iter>=0)//for polyhedrons
          {
            nbOfEltsInRevNodal++;
            revNodalIndxPtr[(*iter)+1]++;
          }
    }
  std::transform(revNodalIndxPtr+1,revNodalIndxPtr+nbOfNodes+1,revNodalIndxPtr,revNodalIndxPtr+1,std::plus<int>());
  int *revNodalPtr=new int[nbOfEltsInRevNodal];
  revNodal->useArray(revNodalPtr,true,CPP_DEALLOC,nbOfEltsInRevNodal,1);
  std::fill(revNodalPtr,revNodalPtr+nbOfEltsInRevNodal,-1);
  for(int eltId=0;eltId<nbOfCells;eltId++)
    {
      const int *strtNdlConnOfCurCell=conn+connIndex[eltId]+1;
      const int *endNdlConnOfCurCell=conn+connIndex[eltId+1];
      for(const int *iter=strtNdlConnOfCurCell;iter!=endNdlConnOfCurCell;iter++)
        if(*iter>=0)//for polyhedrons
          *std::find_if(revNodalPtr+revNodalIndxPtr[*iter],revNodalPtr+revNodalIndxPtr[*iter+1],std::bind2nd(std::equal_to<int>(),-1))=eltId;
    }
}

void MEDCouplingUMesh::zipCoords()
{
  checkFullyDefined();
  DataArrayInt *traducer=zipCoordsTraducer();
  traducer->decrRef();
}

struct MEDCouplingAccVisit
{
  MEDCouplingAccVisit():_new_nb_of_nodes(0) { }
  int operator()(int val) { if(val!=-1) return _new_nb_of_nodes++; else return -1; }
  int _new_nb_of_nodes;
};

/*!
 * Array returned is the correspondance old to new.
 * The maximum value stored in returned array is the number of nodes of 'this' minus 1 after call of this method.
 * The size of returned array is the number of nodes of the old (previous to the call of this method) number of nodes.
 * -1 values in returned array means that the corresponding old node is no more used.
 */
DataArrayInt *MEDCouplingUMesh::zipCoordsTraducer()
{
  DataArrayInt *ret=DataArrayInt::New();
  int nbOfNodes=getNumberOfNodes();
  int spaceDim=getSpaceDimension();
  int *traducer=new int[nbOfNodes];
  std::fill(traducer,traducer+nbOfNodes,-1);
  ret->useArray(traducer,true,CPP_DEALLOC,nbOfNodes,1);
  int nbOfCells=getNumberOfCells();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int *conn=_nodal_connec->getPointer();
  for(int i=0;i<nbOfCells;i++)
    for(int j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        traducer[conn[j]]=1;
  int newNbOfNodes=std::count(traducer,traducer+nbOfNodes,1);
  std::transform(traducer,traducer+nbOfNodes,traducer,MEDCouplingAccVisit());
  for(int i=0;i<nbOfCells;i++)
    for(int j=connIndex[i]+1;j<connIndex[i+1];j++)
      if(conn[j]>=0)
        conn[j]=traducer[conn[j]];
  DataArrayDouble *newCoords=DataArrayDouble::New();
  double *newCoordsPtr=new double[newNbOfNodes*spaceDim];
  const double *oldCoordsPtr=_coords->getConstPointer();
  newCoords->useArray(newCoordsPtr,true,CPP_DEALLOC,newNbOfNodes,spaceDim);
  int *work=std::find_if(traducer,traducer+nbOfNodes,std::bind2nd(std::not_equal_to<int>(),-1));
  for(;work!=traducer+nbOfNodes;work=std::find_if(work,traducer+nbOfNodes,std::bind2nd(std::not_equal_to<int>(),-1)))
    {
      newCoordsPtr=std::copy(oldCoordsPtr+spaceDim*(work-traducer),oldCoordsPtr+spaceDim*(work-traducer+1),newCoordsPtr);
      work++;
    }
  setCoords(newCoords);
  newCoords->decrRef();
  return ret;
}

MEDCouplingPointSet *MEDCouplingUMesh::buildPartOfMySelf(const int *start, const int *end, bool keepCoords) const
{
  MEDCouplingUMesh *ret=buildPartOfMySelfKeepCoords(start,end);
  if(!keepCoords)
    ret->zipCoords();
  return ret;
}

/*!
 * Given a boundary box 'bbox' returns elements 'elems' contained in this 'bbox'.
 * Warning 'elems' is incremented during the call so if elems is not empty before call returned elements will be
 * added in 'elems' parameter.
 */
void MEDCouplingUMesh::giveElemsInBoundingBox(const double *bbox, double eps, std::vector<int>& elems)
{
  int dim=getSpaceDimension();
  double* elem_bb=new double[2*dim];
  const int* conn      = getNodalConnectivity()->getConstPointer();
  const int* conn_index= getNodalConnectivityIndex()->getConstPointer();
  const double* coords = getCoords()->getConstPointer();
  int nbOfCells=getNumberOfCells();
  for ( int ielem=0; ielem<nbOfCells;ielem++ )
    {
      for (int i=0; i<dim; i++)
        {
          elem_bb[i*2]=std::numeric_limits<double>::max();
          elem_bb[i*2+1]=-std::numeric_limits<double>::max();
        }

      for (int inode=conn_index[ielem]+1; inode<conn_index[ielem+1]; inode++)//+1 due to offset of cell type.
        {
          int node= conn[inode];
     
          for (int idim=0; idim<dim; idim++)
            {
              if ( coords[node*dim+idim] < elem_bb[idim*2] )
                {
                  elem_bb[idim*2] = coords[node*dim+idim] ;
                }
              if ( coords[node*dim+idim] > elem_bb[idim*2+1] )
                {
                  elem_bb[idim*2+1] = coords[node*dim+idim] ;
                }
            }
        }
      if (intersectsBoundingBox(elem_bb, bbox, dim, eps))
        {
          elems.push_back(ielem);
        }
    }
  delete [] elem_bb;
}

INTERP_KERNEL::NormalizedCellType MEDCouplingUMesh::getTypeOfCell(int cellId) const
{
  const int *ptI=_nodal_connec_index->getConstPointer();
  const int *pt=_nodal_connec->getConstPointer();
  return (INTERP_KERNEL::NormalizedCellType) pt[ptI[cellId]];
}

int MEDCouplingUMesh::getNumberOfNodesInCell(int cellId) const
{
  int *ptI=_nodal_connec_index->getPointer();
  return ptI[cellId+1]-ptI[cellId]-1;
}

void MEDCouplingUMesh::setConnectivity(DataArrayInt *conn, DataArrayInt *connIndex, bool isComputingTypes)
{
  DataArrayInt::setArrayIn(conn,_nodal_connec);
  DataArrayInt::setArrayIn(connIndex,_nodal_connec_index);
  if(isComputingTypes)
    computeTypes();
}

MEDCouplingUMesh::MEDCouplingUMesh(const MEDCouplingUMesh& other, bool deepCpy):MEDCouplingPointSet(other,deepCpy),_iterator(-1),_mesh_dim(other._mesh_dim),
                                                                                _nodal_connec(0),_nodal_connec_index(0),
                                                                                _types(other._types)
{
  if(other._nodal_connec)
    _nodal_connec=other._nodal_connec->performCpy(deepCpy);
  if(other._nodal_connec_index)
    _nodal_connec_index=other._nodal_connec_index->performCpy(deepCpy);
}

MEDCouplingUMesh::~MEDCouplingUMesh()
{
  if(_nodal_connec)
    _nodal_connec->decrRef();
  if(_nodal_connec_index)
    _nodal_connec_index->decrRef();
}

void MEDCouplingUMesh::computeTypes()
{
  if(_nodal_connec && _nodal_connec_index)
    {
      _types.clear();
      const int *conn=_nodal_connec->getConstPointer();
      const int *connIndex=_nodal_connec_index->getConstPointer();
      int nbOfElem=_nodal_connec_index->getNbOfElems()-1;
      for(const int *pt=connIndex;pt!=connIndex+nbOfElem;pt++)
        _types.insert((INTERP_KERNEL::NormalizedCellType)conn[*pt]);
    }
}

/*!
 * This method checks that all arrays are set. If yes nothing done if no an exception is thrown.
 */
void MEDCouplingUMesh::checkFullyDefined() const throw(INTERP_KERNEL::Exception)
{
  if(!_nodal_connec_index || !_nodal_connec || !_coords)
    throw INTERP_KERNEL::Exception("Reverse nodal connectivity computation requires full connectivity and coordinates set in unstructured mesh.");
}

int MEDCouplingUMesh::getNumberOfCells() const
{ 
  if(_nodal_connec_index)
    if(_iterator==-1)
      return _nodal_connec_index->getNumberOfTuples()-1;
    else
      return _iterator;
  else
    throw INTERP_KERNEL::Exception("Unable to get number of cells because no connectivity specified !");
}

int MEDCouplingUMesh::getMeshLength() const
{
  return _nodal_connec->getNbOfElems();
}

void MEDCouplingUMesh::getTinySerializationInformation(std::vector<int>& tinyInfo) const
{
  tinyInfo.resize(5);
  tinyInfo[0] = getSpaceDimension();
  tinyInfo[1] = getMeshDimension();
  tinyInfo[2] = getNumberOfNodes();
  tinyInfo[3] = getNumberOfCells();
  tinyInfo[4] = getMeshLength();
}

/*!
 * @param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 */
void MEDCouplingUMesh::resizeForSerialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2)
{
  a1->alloc(tinyInfo[4]+tinyInfo[3]+1,1);
  a2->alloc(tinyInfo[2],tinyInfo[0]);
}

void MEDCouplingUMesh::serialize(DataArrayInt *&a1, DataArrayDouble *&a2)
{
  a1=DataArrayInt::New();
  a1->alloc(getMeshLength()+getNumberOfCells()+1,1);
  int *ptA1=a1->getPointer();
  const int *conn=getNodalConnectivity()->getConstPointer();
  const int *index=getNodalConnectivityIndex()->getConstPointer();
  ptA1=std::copy(index,index+getNumberOfCells()+1,ptA1);
  std::copy(conn,conn+getMeshLength(),ptA1);
  a2=getCoords();
  a2->incrRef();
}

/*!
 * @param tinyInfo must be equal to the result given by getTinySerializationInformation method.
 */
MEDCouplingPointSet *MEDCouplingUMesh::buildObjectFromUnserialization(const std::vector<int>& tinyInfo, DataArrayInt *a1, DataArrayDouble *a2)
{
  MEDCouplingUMesh* meshing = MEDCouplingUMesh::New() ;
  // Coordinates
  meshing->setCoords(a2);
  // Connectivity
  const int *recvBuffer=a1->getConstPointer();
  DataArrayInt* myConnecIndex=DataArrayInt::New();
  myConnecIndex->alloc(tinyInfo[3]+1,1);
  std::copy(recvBuffer,recvBuffer+tinyInfo[3]+1,myConnecIndex->getPointer());
  DataArrayInt* myConnec=DataArrayInt::New();
  myConnec->alloc(tinyInfo[4],1);
  std::copy(recvBuffer+tinyInfo[3]+1,recvBuffer+tinyInfo[3]+1+tinyInfo[4],myConnec->getPointer());
  meshing->setConnectivity(myConnec, myConnecIndex) ;
  myConnec->decrRef();
  myConnecIndex->decrRef();
  //
  meshing->setMeshDimension(tinyInfo[1]);
  return meshing;
}

MEDCouplingUMesh *MEDCouplingUMesh::buildPartOfMySelfKeepCoords(const int *start, const int *end) const
{
  checkFullyDefined();
  MEDCouplingUMesh *ret=MEDCouplingUMesh::New();
  std::ostringstream stream; stream << PART_OF_NAME << getName();
  ret->setName(stream.str().c_str());
  ret->_mesh_dim=_mesh_dim;
  ret->setCoords(_coords);
  int nbOfElemsRet=end-start;
  int *connIndexRet=new int[nbOfElemsRet+1];
  connIndexRet[0]=0;
  const int *conn=_nodal_connec->getConstPointer();
  const int *connIndex=_nodal_connec_index->getConstPointer();
  int newNbring=0;
  for(const int *work=start;work!=end;work++,newNbring++)
    connIndexRet[newNbring+1]=connIndexRet[newNbring]+connIndex[*work+1]-connIndex[*work];
  int *connRet=new int[connIndexRet[nbOfElemsRet]];
  int *connRetWork=connRet;
  std::set<INTERP_KERNEL::NormalizedCellType> types;
  for(const int *work=start;work!=end;work++)
    {
      types.insert((INTERP_KERNEL::NormalizedCellType)conn[connIndex[*work]]);
      connRetWork=std::copy(conn+connIndex[*work],conn+connIndex[*work+1],connRetWork);
    }
  DataArrayInt *connRetArr=DataArrayInt::New();
  connRetArr->useArray(connRet,true,CPP_DEALLOC,connIndexRet[nbOfElemsRet],1);
  DataArrayInt *connIndexRetArr=DataArrayInt::New();
  connIndexRetArr->useArray(connIndexRet,true,CPP_DEALLOC,nbOfElemsRet+1,1);
  ret->setConnectivity(connRetArr,connIndexRetArr,false);
  ret->_types=types;
  connRetArr->decrRef();
  connIndexRetArr->decrRef();
  return ret;
}

/*!
 * brief returns the volumes of the cells underlying the field \a field
 *
 * For 2D geometries, the returned field contains the areas.
 * For 3D geometries, the returned field contains the volumes.
 *
 * param field field on which cells the volumes are required
 * return field containing the volumes, area or length depending the meshdimension.
 */
MEDCouplingFieldDouble *MEDCouplingUMesh::getMeasureField() const
{
  int ipt, type ;
  int nbelem       = getNumberOfCells() ;
  int dim_mesh     = getMeshDimension();
  int dim_space    = getSpaceDimension() ;
  const double *coords = getCoords()->getConstPointer() ;
  const int *connec = getNodalConnectivity()->getConstPointer() ;
  const int *connec_index = getNodalConnectivityIndex()->getConstPointer() ;


  MEDCouplingFieldDouble* field = MEDCouplingFieldDouble::New(ON_CELLS);
  DataArrayDouble* array = DataArrayDouble::New() ;
  array->alloc(nbelem, 1) ;
  double *area_vol = array->getPointer() ;

  switch (dim_mesh)
    {
    case 2: // getting the areas
      for ( int iel=0 ; iel<nbelem ; iel++ )
        {
          ipt = connec_index[iel] ;
          type = connec[ipt] ;

          switch ( type )
            {
            case INTERP_KERNEL::NORM_TRI3 :
            case INTERP_KERNEL::NORM_TRI6 :
              {
                int N1 = connec[ipt+1];
                int N2 = connec[ipt+2];
                int N3 = connec[ipt+3];

                area_vol[iel]=INTERP_KERNEL::calculateAreaForTria(coords+(dim_space*N1),
                                                                  coords+(dim_space*N2),
                                                                  coords+(dim_space*N3),
                                                                  dim_space);
              }
              break ;

            case INTERP_KERNEL::NORM_QUAD4 :
            case INTERP_KERNEL::NORM_QUAD8 :
              {
                int N1 = connec[ipt+1];
                int N2 = connec[ipt+2];
                int N3 = connec[ipt+3];
                int N4 = connec[ipt+4];

                area_vol[iel]=INTERP_KERNEL::calculateAreaForQuad(coords+dim_space*N1,
                                                                  coords+dim_space*N2,
                                                                  coords+dim_space*N3,
                                                                  coords+dim_space*N4,
                                                                  dim_space) ;
              }
              break ;

            case INTERP_KERNEL::NORM_POLYGON :
              {
                // We must remember that the first item is the type. That's
                // why we substract 1 to get the number of nodes of this polygon
                int size = connec_index[iel+1] - connec_index[iel] - 1 ;

                const double **pts = new const double *[size] ;

                for ( int inod=0 ; inod<size ; inod++ )
                  {
                    // Remember the first item is the type
                    pts[inod] = coords+dim_space*connec[ipt+inod+1] ;
                  }

                area_vol[iel]=INTERP_KERNEL::calculateAreaForPolyg(pts,size,dim_space);
                delete [] pts;
              }
              break ;

            default :
              throw INTERP_KERNEL::Exception("Bad Support to get Areas on it !");

            } // End of switch

        } // End of the loop over the cells
      break;
    case 3: // getting the volumes
      for ( int iel=0 ; iel<nbelem ; iel++ )
        {
          ipt = connec_index[iel] ;
          type = connec[ipt] ;

          switch ( type )
            {
            case INTERP_KERNEL::NORM_TETRA4 :
            case INTERP_KERNEL::NORM_TETRA10 :
              {
                int N1 = connec[ipt+1];
                int N2 = connec[ipt+2];
                int N3 = connec[ipt+3];
                int N4 = connec[ipt+4];

                area_vol[iel]=INTERP_KERNEL::calculateVolumeForTetra(coords+dim_space*N1,
                                                                     coords+dim_space*N2,
                                                                     coords+dim_space*N3,
                                                                     coords+dim_space*N4) ;
              }
              break ;

            case INTERP_KERNEL::NORM_PYRA5 :
            case INTERP_KERNEL::NORM_PYRA13 :
              {
                int N1 = connec[ipt+1];
                int N2 = connec[ipt+2];
                int N3 = connec[ipt+3];
                int N4 = connec[ipt+4];
                int N5 = connec[ipt+5];

                area_vol[iel]=INTERP_KERNEL::calculateVolumeForPyra(coords+dim_space*N1,
                                                                    coords+dim_space*N2,
                                                                    coords+dim_space*N3,
                                                                    coords+dim_space*N4,
                                                                    coords+dim_space*N5) ;
              }
              break ;

            case INTERP_KERNEL::NORM_PENTA6 :
            case INTERP_KERNEL::NORM_PENTA15 :
              {
                int N1 = connec[ipt+1];
                int N2 = connec[ipt+2];
                int N3 = connec[ipt+3];
                int N4 = connec[ipt+4];
                int N5 = connec[ipt+5];
                int N6 = connec[ipt+6];

                area_vol[iel]=INTERP_KERNEL::calculateVolumeForPenta(coords+dim_space*N1,
                                                                     coords+dim_space*N2,
                                                                     coords+dim_space*N3,
                                                                     coords+dim_space*N4,
                                                                     coords+dim_space*N5,
                                                                     coords+dim_space*N6) ;
              }
              break ;

            case INTERP_KERNEL::NORM_HEXA8 :
            case INTERP_KERNEL::NORM_HEXA20 :
              {
                int N1 = connec[ipt+1];
                int N2 = connec[ipt+2];
                int N3 = connec[ipt+3];
                int N4 = connec[ipt+4];
                int N5 = connec[ipt+5];
                int N6 = connec[ipt+6];
                int N7 = connec[ipt+7];
                int N8 = connec[ipt+8];

                area_vol[iel]=INTERP_KERNEL::calculateVolumeForHexa(coords+dim_space*N1,
                                                                    coords+dim_space*N2,
                                                                    coords+dim_space*N3,
                                                                    coords+dim_space*N4,
                                                                    coords+dim_space*N5,
                                                                    coords+dim_space*N6,
                                                                    coords+dim_space*N7,
                                                                    coords+dim_space*N8) ;
              }
              break ;

            case INTERP_KERNEL::NORM_POLYHED :
              {
                throw INTERP_KERNEL::Exception("Not yet implemented !");
              }
              break ;

            default:
              throw INTERP_KERNEL::Exception("Bad Support to get Volume on it !");
            }
        }
      break;
    default:
      throw INTERP_KERNEL::Exception("interpolation is not available for this dimension");
    }
  
  field->setArray(array) ;
  array->decrRef();
  field->setMesh(const_cast<MEDCouplingUMesh *>(this));  
  return field ;
}
