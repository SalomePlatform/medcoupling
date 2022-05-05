// Copyright (C) 2021-2022  CEA/DEN, EDF R&D
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

#include "MeshFormatWriter.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileField.hxx"
#include "MEDFileData.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "libmesh5.hxx"
#include "MEDMESHConverterUtilities.hxx"
#include <cstring>
#include <algorithm>
#include <map>
#include <cstdlib>
#include <fstream>


namespace MEDCoupling {
  
MeshFormatWriter::MeshFormatWriter()
{}
MeshFormatWriter::MeshFormatWriter(const std::string& meshFileName,
                                   const std::vector<std::string>& fieldFileNames):_meshFileName(meshFileName),
                                   _fieldFileNames(fieldFileNames)
{}
MeshFormatWriter::~MeshFormatWriter()
{}
void MeshFormatWriter::setMeshFileName(const std::string& meshFileName)
{
    _meshFileName = meshFileName;
}
void MeshFormatWriter::setFieldFileNames(const std::vector<std::string>& fieldFileNames)
{
    _fieldFileNames = fieldFileNames;
}
void MeshFormatWriter::setMEDFileDS(MEDCoupling::MEDFileData* mfd)
{

    if(!mfd)
    {
        addMessage( MeshFormat::Comment(" MEDFileData is nullptr! ") << _meshFileName, /*fatal=*/true );
        return;
    }
    if ( !mfd->getNumberOfMeshes())
    {
        addMessage( MeshFormat::Comment("No Mesh in MEDFileData! ") << _meshFileName, /*fatal=*/true );
        return;
    }
    if ( mfd->getNumberOfMeshes() > 1)
    {
        addMessage( MeshFormat::Comment("More than One Mesh in File! ") << _meshFileName, /*fatal=*/true );
        return;
    }

    MEDCoupling::MEDFileMeshes* meshes  = mfd->getMeshes();
    _mesh  = meshes->getMeshAtPos(0);

    _mesh->incrRef();
    MEDCoupling::MEDFileFields* fields = mfd->getFields();

    for (int i = 0; i<fields->getNumberOfFields(); i++ )
    {
        MEDCoupling::MEDFileAnyTypeFieldMultiTS* field = fields->getFieldAtPos(i);
        MEDCoupling::MEDFileFieldMultiTS * f = dynamic_cast<MEDCoupling::MEDFileFieldMultiTS *>(field);
        _fields.push_back(f);
    }


}

void MeshFormatWriter::write()
{

    MeshFormat::Localizer loc;

    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingMesh > mesh = _mesh->getMeshAtLevel( 1 );
    MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh = mesh->buildUnstructured();
    _dim = umesh->getSpaceDimension();

    if (_dim != 2 && _dim != 3)
    {
        addMessage( MeshFormat::Comment("Only 3D or 2D mesh allowed! ") << _meshFileName, /*fatal=*/true );
        return ;
    }


    _version = sizeof(double) < 8 ? 1 : 2;
    _writer = MeshFormat::MeshFormatParser();
    _myCurrentOpenFile = _meshFileName;
    _myCurrentFileId = _writer.GmfOpenMesh( _meshFileName.c_str(), GmfWrite, _version, _dim );
    if ( !_myCurrentFileId )
    {
        if ( MeshFormat::isMeshExtensionCorrect( _meshFileName ))
        {
            addMessage( MeshFormat::Comment("Can't open for writing ") << _meshFileName, /*fatal=*/true );
            return;
        }

        else
        {
            addMessage( MeshFormat::Comment("Not '.mesh' or '.meshb' extension of file ") << _meshFileName, /*fatal=*/true );
            return;
        }

    }


    perform();
    _writer.GmfCloseMesh(_myCurrentFileId);
    _myCurrentFileId = -1;
    _myCurrentOpenFile = "";
    if (_fields.size()) performFields();


}

MeshFormat::Status MeshFormatWriter::perform()
{


    MEDCoupling::MCAuto< MEDCoupling::MEDCouplingMesh > mesh1;
    MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh1;
    MEDCoupling::MCAuto< MEDCoupling::MEDCouplingMesh > mesh2;
    MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh2;
    MEDCoupling::MCAuto< MEDCoupling::MEDCouplingMesh > mesh3;
    MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh3;

    std::vector<int> dims = _mesh->getNonEmptyLevelsExt();
    int dim = _mesh->getMeshDimension();
    bool threeDElements = false;
    bool twoDElements = false;
    bool OneDElements = false;
    if (dims.size() != 0)
    {
        bool maxLevelDimElments = ( std::find(dims.begin(), dims.end(), 0) != dims.end() );
        bool nextToMaxLevelDimElments = ( std::find(dims.begin(), dims.end(), -1) != dims.end() );
        bool nextToNextToMaxLevelDimElments = (std::find(dims.begin(), dims.end(), -2) != dims.end() );
        threeDElements = (dim == 3) ? maxLevelDimElments : false ;
        twoDElements =  (dim == 3) ? nextToMaxLevelDimElments : maxLevelDimElments ;
        OneDElements = (dim == 3) ? nextToNextToMaxLevelDimElments : nextToMaxLevelDimElments;
    }

    MEDCoupling::mcIdType nbEdgesNSEG2 = 0;
    MEDCoupling::mcIdType nbEdgesNSEG3 = 0;
    MEDCoupling::mcIdType nbTRI3 = 0;
    MEDCoupling::mcIdType nbTRI6 = 0;
    MEDCoupling::mcIdType nbQUAD4 = 0;
    MEDCoupling::mcIdType nbQUAD8 = 0;
    MEDCoupling::mcIdType nbQUAD9 = 0;
    MEDCoupling::mcIdType nbTETRA4 = 0;
    MEDCoupling::mcIdType nbTETRA10 = 0;
    MEDCoupling::mcIdType nbPYRA5 = 0;
    MEDCoupling::mcIdType nbHEXA8 = 0;
    MEDCoupling::mcIdType nbHEXA20 = 0;
    MEDCoupling::mcIdType nbHEXA27 = 0;
    MEDCoupling::mcIdType nbPENTA6 = 0;

    if (OneDElements)
    {
        mesh1 = _mesh->getMeshAtLevel( 1-dim );
        umesh1 = mesh1->buildUnstructured();
        nbEdgesNSEG2 = umesh1->getNumberOfCellsWithType(INTERP_KERNEL::NORM_SEG2);

        nbEdgesNSEG3 = umesh1->getNumberOfCellsWithType(INTERP_KERNEL::NORM_SEG3);

    }
    if (twoDElements)
    {
        mesh2 = _mesh->getMeshAtLevel( 2-dim );
        umesh2 = mesh2->buildUnstructured();
        nbTRI3 = umesh2->getNumberOfCellsWithType(INTERP_KERNEL::NORM_TRI3);

        nbTRI6 = umesh2->getNumberOfCellsWithType(INTERP_KERNEL::NORM_TRI6);

        nbQUAD4 = umesh2->getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD4);

        nbQUAD8 = umesh2->getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD8);

        nbQUAD9 = umesh2->getNumberOfCellsWithType(INTERP_KERNEL::NORM_QUAD9);

    }
    if (threeDElements)
    {

        mesh3 = _mesh->getMeshAtLevel( 3-dim );
        umesh3 = mesh3->buildUnstructured();
        nbTETRA4 = umesh3->getNumberOfCellsWithType(INTERP_KERNEL::NORM_TETRA4);

        nbTETRA10 = umesh3->getNumberOfCellsWithType(INTERP_KERNEL::NORM_TETRA10);

        nbPYRA5 = umesh3->getNumberOfCellsWithType(INTERP_KERNEL::NORM_PYRA5);

        nbHEXA8 = umesh3->getNumberOfCellsWithType(INTERP_KERNEL::NORM_HEXA8);

        nbHEXA20 = umesh3->getNumberOfCellsWithType(INTERP_KERNEL::NORM_HEXA20);

        nbHEXA27 = umesh3->getNumberOfCellsWithType(INTERP_KERNEL::NORM_HEXA27);

        nbPENTA6 = umesh3->getNumberOfCellsWithType(INTERP_KERNEL::NORM_PENTA6);

    }



    MEDCoupling::MCAuto< MEDCoupling::MEDCouplingMesh > mesh0 = _mesh->getMeshAtLevel(1);
    MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh0 = mesh0->buildUnstructured();

    // nodes
    getNodes(umesh0);

    // edges SEG2
    if ( nbEdgesNSEG2 > 0)
    {
        getNSEG2( nbEdgesNSEG2, umesh1);
    }
    //~// nodes of quadratic edges SEG3
    if ( nbEdgesNSEG3 > 0)
    {
        getNSEG3( nbEdgesNSEG3, umesh1);
    }

    // triangles TRI3
    if ( nbTRI3 > 0)
    {
        getTRI3( nbTRI3, umesh2);
    }

    // nodes of quadratic triangles
    // triangles TRI6
    if ( nbTRI6 > 0)
    {
        getTRI6( nbTRI6, umesh2);

    }
    //~// quadrangles QUAD4
    if ( nbQUAD4 > 0)
    {
        getQUAD4(nbQUAD4, umesh2);
    }

    //~// quadrangles quadratic QUAD8
    if ( nbQUAD8 > 0)
    {
        getQUAD8( nbQUAD8, umesh2);
    }

    //~// quadrangles quadratic QUAD9
    if ( nbQUAD9 > 0)
    {
        getQUAD9( nbQUAD9, umesh2);
    }

    // terahedra TETRA4
    if ( nbTETRA4 > 0)
    {
        getTETRA4( nbTETRA4, umesh3);
    }
    //~terahedra TETRA10
    if ( nbTETRA10 > 0)
    {
        getTETRA10( nbTETRA10, umesh3);
    }

    //~// pyramids 5
    if ( nbPYRA5 > 0)
    {
        getPYRA5(nbPYRA5, umesh3);
    }
    //~// hexahedra 8
    if ( nbHEXA8 > 0)
    {
        getHEXA8( nbHEXA8, umesh3);
    }

    //~// hexahedra 20
    if ( nbHEXA20 > 0)
    {
        getHEXA20( nbHEXA20, umesh3);
    }

    //~// hexahedra 27
    if ( nbHEXA27 > 0)
    {
        getHEXA27 (nbHEXA27, umesh3);
    }

    // prism
    if ( nbPENTA6 > 0)
    {
        getPENTA6(nbPENTA6, umesh3);
    }

    linkFamilyToCells();
    writeCells();


    return MeshFormat::DRS_OK;
}


MeshFormat::Status MeshFormatWriter::performFields()
{

    MeshFormat::Status status = MeshFormat::Status::DRS_OK;

    if (_fields.size() != _fieldFileNames.size() )
    {
        addMessage( MeshFormat::Comment(" Number of fields and number of input *.sol files must be equal ") << _meshFileName, /*fatal=*/true );
        status = MeshFormat::Status::DRS_FAIL;
        return status;
    }


    int  dim = _mesh->getMeshDimension();  // dim mesh  field lying to
    std::vector<std::string>::const_iterator fieldFileIt = _fieldFileNames.begin();
    int iField = 0;
    std::vector<int> levs {0} ;
    for (; fieldFileIt !=_fieldFileNames.end();  ++fieldFileIt)
    {
        // Open files
        _myCurrentOpenFile = *fieldFileIt;

        MEDCoupling::MEDFileFieldMultiTS* f = _fields[iField];

        if(!f)
            continue;// why???
        if ( f->getMeshName() == _mesh->getName() )
        {

            std::vector< std::vector<MEDCoupling::TypeOfField> > fTypes = f->getTypesOfFieldAvailable();
            std::vector< std::pair<int,int> >  iters = f->getIterations();
            const std::vector<std::string>& compInfo = f->getInfo();
            std::pair<int,int> it = iters[0];

            //~// Open File for writing
            _myCurrentFileId = _writer.GmfOpenMesh( fieldFileIt->c_str(), GmfWrite, _version, _dim );

            if ( fTypes[0].size() == 1 && fTypes[0][0] == MEDCoupling::ON_NODES )
            {
                setFieldOnNodes(f, it.first,  it.second, compInfo.size());
            }
            else
            {
                setFieldOnCells( f, it.first,  it.second, levs );
            }
            //~// Close File
            _writer.GmfCloseMesh( _myCurrentFileId );
        }


        iField++;
    }

    return status;
}

MeshFormat::Status MeshFormatWriter::setFieldOnNodes(MEDCoupling::MEDFileFieldMultiTS * f, int iteration, int order, size_t compSize)
{


    std::vector<INTERP_KERNEL::NormalizedCellType> types;
    std::vector< std::vector<MEDCoupling::TypeOfField> > typesF;
    std::vector< std::vector<std::string> > pfls, locs;
    std::vector< std::vector< std::pair<mcIdType,mcIdType> > > valsVec;
    valsVec = f->getFieldSplitedByType( iteration, order, _mesh->getName().c_str(),
                                        types, typesF, pfls, locs);
    // believe that there can be only one type in a nodal field,
    // so do not perform a loop on types
    const MEDCoupling::DataArrayDouble* valsArray = f->getUndergroundDataArray(iteration, order);
    int typTab[] = { getGmfSolKwd((int)compSize, _dim) };
    _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfSolAtVertices, (int)valsVec[0][0].second, 1, typTab);
    double* valTab0 = new double[compSize];
    double* valTab;
    for ( size_t i = valsVec[0][0].first; i < (std::size_t)valsVec[0][0].second; ++i )
    {

        for ( size_t j = 0; j < compSize; ++j )
            valTab0[j] = valsArray->getIJ( i, j );
            
    if (compSize == 9 || compSize == 4){ // full matrix ==>uper triangular matrix
      extractSymetricTensor(valTab0, valTab);
      _writer.GmfSetLin( _myCurrentFileId, MeshFormat::GmfSolAtVertices, valTab);
      delete [] valTab;
      
    }
    else  
        _writer.GmfSetLin( _myCurrentFileId, MeshFormat::GmfSolAtVertices, valTab0);

    }
    delete [] valTab0;

    return MeshFormat::Status::DRS_OK;

}

MeshFormat::Status MeshFormatWriter::setFieldOnCells(MEDCoupling::MEDFileFieldMultiTS * f, int iteration, int order, std::vector<int> levs )
{

    int  dim = _mesh->getMeshDimension();  // dim mesh  field lying to
    int absDim = f->getNonEmptyLevels(iteration, order,  f->getMeshName(), levs);

    MEDCoupling::MEDCouplingFieldDouble**  cellToNodeFldb  = new MEDCoupling::MEDCouplingFieldDouble* [(int)levs.size()] ;
    MEDCoupling::MEDCouplingFieldDouble**  fldb  = new MEDCoupling::MEDCouplingFieldDouble* [(int)levs.size()] ;

    for (size_t k = 0; k<levs.size(); k++) fldb[k] = f->field( iteration, order,_mesh );

    // turn it node discretization
    for (size_t l = 0; l < levs.size(); l++) cellToNodeFldb[l] = fldb[l]->cellToNodeDiscretization() ;

    for(size_t j =0; j < levs.size(); j++ )
    {

        const mcIdType pointsNumber = cellToNodeFldb[j]->getNumberOfTuples();
        const mcIdType nbComp = (int) cellToNodeFldb[j]->getNumberOfComponents() ;

        MEDCoupling::DataArrayDouble* timeStamp = cellToNodeFldb[j]->getArray();
        double* values = timeStamp->getPointer();

        int typ = getGmfSolKwd((int)nbComp, _dim) ;
        if(typ == -1)
        {
            addMessage( MeshFormat::Comment(" error with Number of Component   ") << nbComp, /*fatal=*/true );
            return MeshFormat::Status::DRS_FAIL;
        }

        int typTab[] = {typ};
        _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfSolAtVertices, pointsNumber, 1, typTab);


        double *valTab;
        for (int i = 0; i < pointsNumber ; i++ )
        {

            double valTab0[10]; //max 3x3 matrix +1 for safety ;)

            std::copy(values, values+nbComp, valTab0);

            if (nbComp == 9 || nbComp == 4) // full matrix ==>uper triangular matrix
            {
        extractSymetricTensor(valTab0, valTab);
                _writer.GmfSetLin( _myCurrentFileId, MeshFormat::GmfSolAtVertices, valTab);
                delete [] valTab;
            }
            else //sym mat, scalar or vec
            {
                valTab = new double[nbComp];
                std::copy(valTab0, valTab0+nbComp, valTab);
                _writer.GmfSetLin( _myCurrentFileId, MeshFormat::GmfSolAtVertices, valTab);
                delete [] valTab;
            }
            values+=nbComp;

        }


    }

    for(size_t i = 0; i < levs.size(); i++ )  fldb[i]->decrRef();
    for(size_t i = 0; i < levs.size(); i++ )  cellToNodeFldb[i]->decrRef();
    delete [] cellToNodeFldb;
    delete [] fldb;

    return MeshFormat::Status::DRS_OK;
}
 /*\
 |*| extract the upper triangular matrix  of fullTensor
 |*| if _dim == 2 fill symTensor with values at index 0, 1 & 3 of fullTensor
 |*| |x0 x1|   
 |*| |x2 x3|  
 |*| if _dim == 3 fill symTensor with values at index 0, 1, 2, 4, 5 & 8 of fullTensor
 |*| |x0 x1 x2|
 |*| |x3 x4 x5|
 |*| |x6 x7 x8|
 \*/
void MeshFormatWriter::extractSymetricTensor(double fullTensor[], double*& symTensor)
{
     symTensor = new double[_dim*(_dim+1)/2];
  for (int ii =0; ii<_dim; ii++)
    for (int jj =ii; jj<_dim; jj++)
    {
      int kk = _dim*(_dim-1)/2- (_dim-ii)*(_dim-ii-1)/2+jj;
      symTensor[kk] = fullTensor[ii+jj*_dim];
    }  
}
int MeshFormatWriter::getGmfSolKwd(const int nbComp, const int dim)
{
    if (nbComp== 1) return GmfSca;
    else if( dim == nbComp) return GmfVec;
    else if (dim*(dim+1)/2 == nbComp || dim*dim == nbComp ) return GmfSymMat;
    //~else if (dim*dim == nbComp) return GmfMat; // Not valid in mg-adapt if not sym
    else  return -1;
}
bool MeshFormatWriter::checkFileName()
{
    bool ret = true;
    return ret;
}
bool MeshFormatWriter::checkFieldFileName()
{
    bool ret = true;
    return ret;

}

std::string MeshFormatWriter::getMeshFileName() const
{
    return _meshFileName;
}


std::vector<std::string> MeshFormatWriter::getFieldFileNames() const
{
    return _fieldFileNames;
}

MeshFormat::Status MeshFormatWriter::addMessage(const std::string& msg,
        const bool         isFatal/*=false*/)
{
    if ( isFatal )
        _myErrorMessages.clear(); // warnings are useless if a fatal error encounters

    _myErrorMessages.push_back( msg );

    //~MESSAGE(msg);
#ifdef _DEBUG_
    std::cout << msg << std::endl;
#endif
    return ( _myStatus = isFatal ? MeshFormat::DRS_FAIL : MeshFormat::DRS_WARN_SKIP_ELEM );
}


void MeshFormatWriter::forward_shift(std::vector<MEDCoupling::mcIdType> &conn)
{
    std::vector<MEDCoupling::mcIdType>::iterator it = conn.begin();
    for (; it != conn.end(); ++it) *it = *it+1;
}


void MeshFormatWriter::getNodes(MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh0)
{
    MEDCoupling::mcIdType nbNodes = 0;
    nbNodes = umesh0->getNumberOfNodes();


    // nodes
    MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> coordArray = umesh0->getCoordinatesAndOwner();
    double* coordPrt = coordArray->getPointer();
    _writer.GmfSetKwd( _myCurrentFileId, MeshFormat::GmfVertices, nbNodes );
    double xyz[3];
    int j = (int)nbNodes;

    int idNode = 0;

    while ( j >0 )
    {

        std::copy(coordPrt, coordPrt+_dim, xyz);

        MeshFormatNode e(xyz[0], xyz[1], xyz[2], idNode);
        _idNodeToNode.insert(std::pair <int, MeshFormatNode> (idNode, e));

        coordPrt+= _dim;
        j--;
        idNode++;
    }
    linkFamilyToNodes();
    std::map <int, MeshFormatNode>::iterator itNode = _idNodeToNode.begin();
    for (; itNode!= _idNodeToNode.end(); ++itNode)
        _dim == 3?  _writer.GmfSetLin( _myCurrentFileId, MeshFormat::GmfVertices, itNode->second.xyz[0],
                     itNode->second.xyz[1], itNode->second.xyz[2], std::abs(itNode->second._famId) ) :
                      _writer.GmfSetLin( _myCurrentFileId, MeshFormat::GmfVertices, itNode->second.xyz[0],
                     itNode->second.xyz[1], std::abs(itNode->second._famId) );
}


void MeshFormatWriter::getNSEG2(MEDCoupling::mcIdType nbEdgesNSEG2, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh1)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh1->giveCellsWithType(INTERP_KERNEL::NORM_SEG2);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh1->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);

        MeshFormatCell e(INTERP_KERNEL::NORM_SEG2, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_SEG2, idCellToCell) );
}


void MeshFormatWriter::getNSEG3( MEDCoupling::mcIdType nbEdgesNSEG3, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh1)
{


    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh1->giveCellsWithType(INTERP_KERNEL::NORM_SEG3);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {

        std::vector<MEDCoupling::mcIdType> conn;
        umesh1->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_SEG3, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_SEG3, idCellToCell) );
}


void MeshFormatWriter::getTRI3( MEDCoupling::mcIdType nbTRI3, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh2)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh2->giveCellsWithType(INTERP_KERNEL::NORM_TRI3);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh2->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_TRI3, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_TRI3, idCellToCell) );
}


void MeshFormatWriter::getTRI6( MEDCoupling::mcIdType nbTRI6, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh2)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh2->giveCellsWithType(INTERP_KERNEL::NORM_TRI6);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh2->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_TRI6, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }

    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_TRI6, idCellToCell) );
}

void MeshFormatWriter::getQUAD4( MEDCoupling::mcIdType nbQUAD4, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh2)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh2->giveCellsWithType(INTERP_KERNEL::NORM_QUAD4);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh2->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_QUAD4, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_QUAD4, idCellToCell) );
}

void MeshFormatWriter::getQUAD8(MEDCoupling::mcIdType nbQUAD8, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh2)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh2->giveCellsWithType(INTERP_KERNEL::NORM_QUAD8);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh2->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_QUAD8, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_QUAD8, idCellToCell) );
}

void MeshFormatWriter::getQUAD9(MEDCoupling::mcIdType nbQUAD9, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh2)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh2->giveCellsWithType(INTERP_KERNEL::NORM_QUAD9);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh2->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_QUAD9, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }

    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_QUAD9, idCellToCell) );
}

void MeshFormatWriter::getTETRA4(MEDCoupling::mcIdType nbTETRA4, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh3)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh3->giveCellsWithType(INTERP_KERNEL::NORM_TETRA4);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh3->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_TETRA4, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }

    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_TETRA4, idCellToCell) );
}

void MeshFormatWriter::getTETRA10(MEDCoupling::mcIdType nbTETRA10, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh3)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh3->giveCellsWithType(INTERP_KERNEL::NORM_TETRA10);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh3->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_TETRA10, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_TETRA10, idCellToCell) );
}

void MeshFormatWriter::getPYRA5(MEDCoupling::mcIdType nbPYRA5, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh3)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh3->giveCellsWithType(INTERP_KERNEL::NORM_PYRA5);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh3->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_PYRA5, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_PYRA5, idCellToCell) );
}

void MeshFormatWriter::getHEXA8(MEDCoupling::mcIdType nbHEXA8, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh3)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh3->giveCellsWithType(INTERP_KERNEL::NORM_HEXA8);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh3->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_HEXA8, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_HEXA8, idCellToCell) );
}
void MeshFormatWriter::getHEXA20(MEDCoupling::mcIdType nbHEXA20, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh3)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh3->giveCellsWithType(INTERP_KERNEL::NORM_HEXA20);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh3->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_HEXA20, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_HEXA20, idCellToCell) );
}
void MeshFormatWriter::getHEXA27(MEDCoupling::mcIdType nbHEXA27, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh3)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh3->giveCellsWithType(INTERP_KERNEL::NORM_HEXA27);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh3->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_HEXA27, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_HEXA27, idCellToCell) );
}

void MeshFormatWriter::getPENTA6(MEDCoupling::mcIdType nbPENTA6, MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh3)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> elementId = umesh3->giveCellsWithType(INTERP_KERNEL::NORM_PENTA6);
    std::map<int, MeshFormatCell> idCellToCell;
    for ( const mcIdType *it=elementId->begin(); it!=elementId->end(); it++ )
    {
        std::vector<MEDCoupling::mcIdType> conn;
        umesh3->getNodeIdsOfCell(*it,  conn) ;
        forward_shift(conn);
        MeshFormatCell e(INTERP_KERNEL::NORM_PENTA6, (int)*it);
        e.setConn(conn);
        idCellToCell.insert(std::pair <int, MeshFormatCell> (*it, e));
    }
    _typeToIdCellToCell.insert(std::pair <INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >(INTERP_KERNEL::NORM_PENTA6, idCellToCell) );
}


void MeshFormatWriter::linkFamilyToNodes()
{

    std::map<std::string,mcIdType> famInfos =  _mesh->getFamilyInfo();
    std::map<std::string,mcIdType>::const_iterator famIt =  famInfos.begin();
    for (; famIt != famInfos.end(); ++famIt)
    {
        if(!famIt->second) continue; //FAMILLE_ZERO
        MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> nodeIds =  _mesh->getNodeFamilyArr(famIt->first);
        const MEDCoupling::mcIdType * nodeIdsIt = nodeIds->begin(), * famIDEnd = nodeIds->end();
        for(; nodeIdsIt< famIDEnd; ++nodeIdsIt) {

            std::map <int, MeshFormatNode>::iterator itNode = _idNodeToNode.find((int)*nodeIdsIt);
            if (itNode == _idNodeToNode.end()) continue;
            else itNode->second._famId =(int) famIt->second;


        }
    }
}



void MeshFormatWriter::linkFamilyToCells()
{

    std::vector<int> levs =  _mesh->getNonEmptyLevels();
    for (size_t iDim = 0; iDim < levs.size(); iDim++ )
    {
        int meshDimRelToMax = levs[iDim];
        MEDCoupling::MCAuto< MEDCoupling::MEDCouplingMesh > mesh = _mesh->getMeshAtLevel( meshDimRelToMax);
        MEDCoupling::MCAuto< MEDCoupling::MEDCouplingUMesh > umesh0 = mesh->buildUnstructured();
        const MEDCoupling::DataArrayIdType * famIds = _mesh->getFamilyFieldAtLevel(meshDimRelToMax);
        const MEDCoupling::mcIdType * famID = famIds->begin(), *famIDEnd = famIds->end();
        for (; famID < famIDEnd; ++famID)
        {
            if (!(*famID)) continue; // "FAMILLE_ZERO"
            std::string famName = _mesh->getFamilyNameGivenId(*famID);
            MEDCoupling::MCAuto<MEDCoupling::DataArrayIdType> cellIds =  _mesh->getFamilyArr( meshDimRelToMax, famName);
            const MEDCoupling::mcIdType * cellIdsIt = cellIds->begin(), *cellIDEnd = cellIds->end();
            for(; cellIdsIt< cellIDEnd; ++cellIdsIt)
            {
                INTERP_KERNEL::NormalizedCellType type = umesh0->getTypeOfCell(*cellIdsIt); //TODO
                std::map<INTERP_KERNEL::NormalizedCellType, std::map <int, MeshFormatCell> >::iterator itCellMap = _typeToIdCellToCell.find(type);
                if (itCellMap == _typeToIdCellToCell.end()) continue;
                else
                {
                    std::map <int, MeshFormatCell>::iterator itCell = itCellMap->second.find((int)*cellIdsIt);
                    if (itCell == itCellMap->second.end()) continue;
                    else itCell->second._famId = (int)*famID;
                }

            }



        }
    }
}
void MeshFormatWriter::writeCells()
{

    std::map < INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> >::iterator typeCellMapIt = _typeToIdCellToCell.begin();
    for (; typeCellMapIt!= _typeToIdCellToCell.end(); ++typeCellMapIt)
    {
        std::map<int, MeshFormatCell>::iterator cellMapIt = typeCellMapIt->second.begin();
        switch (typeCellMapIt->first)
        {
        case INTERP_KERNEL::NORM_SEG2 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfEdges, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfEdges, conn[0], conn[1], std::abs(cellMapIt->second._famId) );
            }
            break;
        }
        case INTERP_KERNEL::NORM_SEG3 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfEdges, (int)typeCellMapIt->second.size());
            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfEdges, conn[0], conn[1], std::abs(cellMapIt->second._famId) );
            }
            cellMapIt = typeCellMapIt->second.begin();

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtEdges, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtEdges, cellMapIt->first+1, 1, conn[2] );
            }
            break;
        }
        case INTERP_KERNEL::NORM_TRI3 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfTriangles, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfTriangles, conn[0], conn[1], conn[2], std::abs(cellMapIt->second._famId) );
            }
            break;
        }
        case INTERP_KERNEL::NORM_TRI6 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfTriangles, (int)typeCellMapIt->second.size());
            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfTriangles, conn[0], conn[1], conn[2], std::abs(cellMapIt->second._famId) );
            }

            cellMapIt = typeCellMapIt->second.begin();
            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtTriangles, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtTriangles, cellMapIt->first+1, 3, conn[3], conn[4], conn[5] );
            }
            break;
        }
        case INTERP_KERNEL::NORM_QUAD4 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfQuadrilaterals, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfQuadrilaterals, conn[0], conn[1], conn[2], conn[3], std::abs(cellMapIt->second._famId) );
            }
            break;
        }
        case INTERP_KERNEL::NORM_QUAD8 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfQuadrilaterals, (int)typeCellMapIt->second.size());
            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfQuadrilaterals, conn[0], conn[1], conn[2], conn[3], std::abs(cellMapIt->second._famId) );
            }
            cellMapIt = typeCellMapIt->second.begin();
            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtQuadrilaterals, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtQuadrilaterals, cellMapIt->first+1, 4, conn[4], conn[5],
                                  conn[6], conn[7] );
            }
            break;
        }
        case INTERP_KERNEL::NORM_QUAD9 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfQuadrilaterals, (int)typeCellMapIt->second.size());
            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfQuadrilaterals, conn[0], conn[1], conn[2], conn[3], std::abs(cellMapIt->second._famId) );
            }
            cellMapIt = typeCellMapIt->second.begin();
            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtQuadrilaterals, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtQuadrilaterals,cellMapIt->first+1, 5, conn[4], conn[5],
                                  conn[6], conn[7], conn[8] );
            }
            break;
        }
        case INTERP_KERNEL::NORM_TETRA4 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfTetrahedra, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfTetrahedra, conn[0], conn[2], conn[1], conn[3],  std::abs(cellMapIt->second._famId) );
            }
            break;
        }
        case INTERP_KERNEL::NORM_TETRA10 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtTetrahedra, (int)typeCellMapIt->second.size());
            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfTetrahedra, conn[0], conn[2], conn[1], conn[3], std::abs(cellMapIt->second._famId) );
            }
            cellMapIt = typeCellMapIt->second.begin();
            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfTetrahedra, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtTetrahedra, cellMapIt->first+1, 6, conn[6], conn[5],
                                  conn[4], conn[7], conn[8], conn[9] );
            }
            break;
        }
        case INTERP_KERNEL::NORM_PYRA5 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfPyramids, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfPyramids, conn[3], conn[2], conn[1], conn[0], conn[4], std::abs(cellMapIt->second._famId) );

            }
            break;
        }
        case INTERP_KERNEL::NORM_HEXA8 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfHexahedra, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfHexahedra, conn[0], conn[3], conn[2], conn[1], conn[4], conn[7], conn[6], conn[5], std::abs(cellMapIt->second._famId) );

            }
            break;
        }
        case INTERP_KERNEL::NORM_HEXA20 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfHexahedra, (int)typeCellMapIt->second.size());
            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfHexahedra, conn[0], conn[3], conn[2], conn[1],
                                  conn[4], conn[7], conn[6], conn[5], std::abs(cellMapIt->second._famId) );
            }
            cellMapIt = typeCellMapIt->second.begin();
            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtHexahedra, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtHexahedra,  cellMapIt->first+1, 12, conn[11], conn[10], conn[9],
                                  conn[8], conn[15], conn[14], conn[13], conn[12], conn[16], conn[19], conn[18], conn[17] );
            }
            break;
        }
        case INTERP_KERNEL::NORM_HEXA27 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfHexahedra, (int)typeCellMapIt->second.size());
            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtHexahedra, conn[0], conn[3], conn[2], conn[1],
                                  conn[4], conn[7], conn[6], conn[5], std::abs(cellMapIt->second._famId) );
            }
            cellMapIt = typeCellMapIt->second.begin();
            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtHexahedra, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtHexahedra, cellMapIt->first+1, 19, conn[11], conn[10], conn[9],
                                  conn[8], conn[15], conn[14], conn[13], conn[12], conn[16], conn[19], conn[18], conn[17],
                                  conn[20], conn[24], conn[23], conn[22], conn[21], conn[25], conn[26], std::abs(cellMapIt->second._famId) );
            }
            break;
        }
        case INTERP_KERNEL::NORM_PENTA6 :
        {

            _writer.GmfSetKwd(_myCurrentFileId, MeshFormat::GmfPrisms, (int)typeCellMapIt->second.size());

            for (; cellMapIt != typeCellMapIt->second.end(); ++cellMapIt)
            {
                std::vector<MEDCoupling::mcIdType> conn = cellMapIt->second.conn;
                _writer.GmfSetLin(_myCurrentFileId, MeshFormat::GmfPrisms, conn[0], conn[2], conn[1], conn[3], conn[5], conn[4], std::abs(cellMapIt->second._famId) );
            }
            break;
        }
        }
    }
}
}
