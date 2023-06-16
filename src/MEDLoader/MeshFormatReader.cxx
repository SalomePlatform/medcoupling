// Copyright (C) 2021-2023  CEA/DEN, EDF R&D
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

#include "MeshFormatReader.hxx"

#include "MEDFileMesh.hxx"
#include "MEDFileField.hxx"
#include "MEDFileData.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "libmesh5.hxx"
#include "MEDMESHConverterUtilities.hxx"
#include <cstring>
#include <fstream>

namespace MEDCoupling {
MeshFormatReader::MeshFormatReader():_myMeshName("MESH")
{}  
MeshFormatReader::MeshFormatReader(const std::string& meshFileName,
                                   const std::vector<std::string>& fieldFileName):_myFile(meshFileName),
                                   _myFieldFileNames(fieldFileName),
                                   _myCurrentFileId(0),
                                   _myMeshName("MESH")
{}

MeshFormatReader::~MeshFormatReader()
{}

MEDCoupling::MCAuto<MEDCoupling::MEDFileData> MeshFormatReader::loadInMedFileDS()
{
    _myStatus = perform();
    if(_myStatus != MeshFormat::DRS_OK) return 0;

    if ( !_uMesh->getName().c_str() || strlen( _uMesh->getName().c_str() ) == 0 )
        _uMesh->setName( _myMeshName );
    if ( _myFieldFileNames.size() ) performFields();


    MEDCoupling::MCAuto< MEDCoupling::MEDFileMeshes > meshes = MEDCoupling::MEDFileMeshes::New();
    _myMed = MEDCoupling::MEDFileData::New();
    meshes->pushMesh( _uMesh );
    _myMed->setMeshes( meshes );

    if ( _fields ) _myMed->setFields( _fields );

    return _myMed.retn();
}

//================================================================================
/*!
 * \brief Read a GMF file
 */
//================================================================================

MeshFormat::Status MeshFormatReader::perform()
{
    MeshFormat::Localizer loc;

    MeshFormat::Status status = MeshFormat::DRS_OK;


    // open the file
     _reader = MeshFormat::MeshFormatParser();
    _myCurrentOpenFile = _myFile;
    _myCurrentFileId = _reader.GmfOpenMesh( _myFile.c_str(), GmfRead, &_version, &_dim );
    if ( !_myCurrentFileId )
    {
        if ( MeshFormat::isMeshExtensionCorrect( _myFile ))
        {

            return addMessage( MeshFormat::Comment("Can't open for reading ") << _myFile,
                               /*fatal=*/true );
        }
        else
            return addMessage( MeshFormat::Comment("Not '.mesh' or '.meshb' extension of file ") << _myFile,
                               /*fatal=*/true );
    }


    // Read nodes


    MEDCoupling::DataArrayDouble* coordArray = MEDCoupling::DataArrayDouble::New();
    setNodes(coordArray);


    int nbEdges = _reader.GmfStatKwd(_myCurrentFileId, MeshFormat::GmfEdges);
    int nbTria = _reader.GmfStatKwd(_myCurrentFileId, MeshFormat::GmfTriangles);
    int nbQuad = _reader.GmfStatKwd(_myCurrentFileId, MeshFormat::GmfQuadrilaterals);
    int nbTet = _reader.GmfStatKwd( _myCurrentFileId, MeshFormat::GmfTetrahedra );
    int nbPyr = _reader.GmfStatKwd(_myCurrentFileId, MeshFormat::GmfPyramids);
    int nbHex = _reader.GmfStatKwd(_myCurrentFileId, MeshFormat::GmfHexahedra);
    int nbPrism = _reader.GmfStatKwd(_myCurrentFileId, MeshFormat::GmfPrisms);

    _dim1NbEl = nbEdges;
    _dim2NbEl = nbTria + nbQuad;
    _dim3NbEl = nbTet + nbPyr + nbHex + nbPrism;
    bool okdim1 = (nbEdges > 0), okdim2= (_dim2NbEl > 0), okdim3= (_dim3NbEl > 0);

    MEDCoupling::MEDCouplingUMesh* dimMesh1;
    MEDCoupling::MEDCouplingUMesh* dimMesh2;
    MEDCoupling::MEDCouplingUMesh* dimMesh3;

    // dim 1
    if (okdim1 ) {
        dimMesh1 = MEDCoupling::MEDCouplingUMesh::New();
        dimMesh1->setCoords( coordArray );
        dimMesh1->allocateCells( nbEdges);
        dimMesh1->setMeshDimension( 1 );
    }
    // dim 2
    if (okdim2) {
        dimMesh2 = MEDCoupling::MEDCouplingUMesh::New();
        dimMesh2->setCoords( coordArray );
        dimMesh2->allocateCells(_dim2NbEl);
        dimMesh2->setMeshDimension( 2 );
    }
    // dim 3
    if (okdim3) {
        dimMesh3 = MEDCoupling::MEDCouplingUMesh::New();
        dimMesh3->setCoords( coordArray );
        dimMesh3->allocateCells(_dim3NbEl);
        dimMesh3->setMeshDimension( 3 );
    }

    // Read elements

    int iN[28]; // 28 - nb nodes in HEX27 (+ 1 for safety :)

    /* Read edges */

    if ( nbEdges)
    {
        setEdges(dimMesh1, nbEdges);
        dimMesh1->decrRef();
    }

    /* Read triangles */
    if (nbTria)
    {
        setTriangles(dimMesh2, nbTria);
    }

    /* Read quadrangles */
    if (nbQuad)
    {
        setQuadrangles( dimMesh2, nbQuad);
    }

    if (okdim2 ) {
        dimMesh2->finishInsertingCells();
        _uMesh->setMeshAtLevel( 2 - _dim, dimMesh2 );
        dimMesh2->sortCellsInMEDFileFrmt();
        dimMesh2->decrRef();
    }
    /* Read terahedra */
    if (  nbTet )
    {
        setTetrahedras( dimMesh3, nbTet);
    }

    /* Read pyramids */
    if ( nbPyr )
    {
        setPyramids( dimMesh3, nbPyr);
    }

    /* Read hexahedra */
    if ( nbHex )
    {
        setHexahedras(dimMesh3, nbHex);
    }

    /* Read prism */
    //~const int prismIDShift = myMesh->GetMeshInfo().NbElements();
    if ( nbPrism )
    {
        setPrisms(dimMesh3, nbPrism);
    }

    if (okdim3 ) {
        dimMesh3->finishInsertingCells();
        _uMesh->setMeshAtLevel( 3 - _dim, dimMesh3 );
        dimMesh3->decrRef();
    }

    buildFamilies();
    coordArray->decrRef();
    _reader.GmfCloseMesh(_myCurrentFileId);
    _myCurrentFileId = -1;
    _myCurrentOpenFile ="";
    return status;
}

MeshFormat::Status MeshFormatReader::performFields()
{

    MeshFormat::Status status = MeshFormat::DRS_OK;
    _fields = MEDCoupling::MEDFileFields::New();


    const MeshFormat::GmfKwdCod meshFormatSol[1] = { MeshFormat::GmfSolAtVertices};  // ONLY VERTEX SOL FOR NOW

    int dim, version;
    std::vector<std::string>::const_iterator fieldFileIt = _myFieldFileNames.begin();

    for (; fieldFileIt !=_myFieldFileNames.end();  ++fieldFileIt)
    {
        _myCurrentOpenFile = *fieldFileIt;
        _myCurrentFileId = _reader.GmfOpenMesh( fieldFileIt->c_str(), GmfRead, &version, &dim );
        if ( !_myCurrentFileId )
        {
            if ( MeshFormat::isMeshExtensionCorrect( *fieldFileIt ))
            {

                return addMessage( MeshFormat::Comment("Can't open for reading ") << *fieldFileIt,
                                   /*fatal=*/true );
            }
            else
                return addMessage( MeshFormat::Comment("Not '.sol' or '.solb' extension of file ") << _myFile,
                                   /*fatal=*/true );
        }

        if (version != _version)
        {

            addMessage( MeshFormat::Comment("Warning sol file version is different than mesh file version ") << *fieldFileIt,
                               /*fatal=*/false );
        }
        if (dim != _dim)
        {

            return addMessage( MeshFormat::Comment("Error sol file must have same dimension as mesh file ") << *fieldFileIt,
                               /*fatal=*/true );
        }


        MeshFormat::GmfKwdCod kwd = meshFormatSol[0];
        int NmbSol, NmbTypes, NmbReals, TypesTab[ GmfMaxTyp ];
        NmbSol = _reader.GmfStatKwd( _myCurrentFileId, kwd, &NmbTypes, &NmbReals, TypesTab );
        if(NmbSol)
        {

            for(int i=0; i<NmbTypes; i++)
            {
                _reader.GmfGotoKwd(_myCurrentFileId, kwd);
                switch(TypesTab[i])
                {
                case GmfSca:
                {

                    setFields(kwd, NmbSol, 1);
                    break;
                }
                case GmfVec:
                {

                    setFields(kwd, NmbSol, dim);
                    break;
                }
                case GmfSymMat :
                {

                    setFields(kwd, NmbSol, dim*(dim+1)/2 );
                    break;
                }
                case GmfMat:
                {

                    setFields(kwd, NmbSol, dim*dim);
                    break;
                }
                }

            }



        }

  
    
        _reader.GmfCloseMesh(_myCurrentFileId);
        _myCurrentFileId = -1;
        _myCurrentOpenFile ="";
    }
    return status;
}

void MeshFormatReader::setFields( MeshFormat::GmfKwdCod kwd, int nmbSol, int nbComp)
{

    MEDCoupling::MCAuto<MEDCoupling::DataArrayDouble> fieldValues =  MEDCoupling::DataArrayDouble::New();
    fieldValues->alloc(nmbSol, nbComp);
    double* values = fieldValues->getPointer();

    int ref;
    double *val = new double[nbComp];

    bool isOnAll = (_uMesh->getNumberOfNodes()== nmbSol);


    for(int i = 1; i<= nmbSol; i++)
    {
        callParserGetLin(kwd, val, nbComp, &ref);

        std::copy(val, val+nbComp, values);

        values+=nbComp;
    }
    values=fieldValues->getPointer(); // to delete

    delete []val;

    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingFieldDouble> timeStamp;
    MEDCoupling::MCAuto<MEDCoupling::MEDFileFieldMultiTS> tsField =  MEDCoupling::MEDFileFieldMultiTS::New();
    MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> dimMesh ;
    int  dimRel;
    MEDCoupling::TypeOfField typeOfField;
    setTypeOfFieldAndDimRel(kwd,  &typeOfField,  &dimRel );

    timeStamp = MEDCoupling::MEDCouplingFieldDouble::New(typeOfField);
    dimMesh = _uMesh->getMeshAtLevel(dimRel);


    timeStamp->setMesh( dimMesh );
    std::string name = "Field_on_Vertex";
    timeStamp->setName(name);
    timeStamp->setArray(fieldValues);
    // set an order
    const int nbTS = tsField->getNumberOfTS();
    if ( nbTS > 0 )
        timeStamp->setOrder( nbTS );

    // add the time-stamp
    timeStamp->checkConsistencyLight();

    if (isOnAll)
    {
        tsField->appendFieldNoProfileSBT( timeStamp );
    }


    _fields->pushField( tsField);

}

void MeshFormatReader::setMeshName(const std::string& theMeshName)
{
    _myMeshName = theMeshName;
}

std::string MeshFormatReader::getMeshName() const
{
    return _myMeshName;
}


void MeshFormatReader::setFile(const std::string& theFileName)
{
    _myFile = theFileName;
}

void MeshFormatReader::setFieldFileNames(const std::vector<std::string>& theFieldFileNames)
{
    _myFieldFileNames = theFieldFileNames;
}
std::vector<std::string> MeshFormatReader::getFieldFileNames() const
{
    return _myFieldFileNames;
}

MeshFormat::Status MeshFormatReader::addMessage(const std::string& msg,
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

void MeshFormatReader::backward_shift(mcIdType* tab, int size)
{
    for(int i = 0; i<size; i++) tab[i] = tab[i]-1;
}


MeshFormat::Status MeshFormatReader::setNodes( MEDCoupling::DataArrayDouble* coordArray)
{

    MeshFormat::Status status = MeshFormat::DRS_OK;
    int nbNodes = _reader.GmfStatKwd(_myCurrentFileId, MeshFormat::GmfVertices);
    if ( nbNodes < 1 )
        return addMessage( "No nodes in the mesh", /*fatal=*/true );

    _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfVertices);

    int ref;
    _uMesh = MEDCoupling::MEDFileUMesh::New();
    coordArray->alloc( nbNodes, _dim );
    double* coordPrt = coordArray->getPointer();
    std::vector<double> nCoords (_dim+2, 0.0) ;

    double* coordPointer  = &nCoords[0];

    if ( _version != GmfFloat )
    {
        double x, y, z;

        for ( int i = 1; i <= nbNodes; ++i )
        {
            _dim == 2 ? _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfVertices, &x, &y, &ref) : \
            _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfVertices, &x, &y, &z, &ref);
            nCoords[0] = x;
            nCoords[1] = y;
            nCoords[2] = z;
            std::copy(coordPointer, coordPointer+_dim, coordPrt);
            MeshFormatElement e(MeshFormat::GmfVertices, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 1);
            coordPrt += _dim;

        }

    }
    else
    {
        float x, y, z;

        for ( int i = 1; i <= nbNodes; ++i )
        {
            _dim == 2 ? _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfVertices, &x, &y, &ref) :
            _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfVertices, &x, &y, &z, &ref);
            nCoords[0] = x;
            nCoords[1] = y;
            nCoords[2] = z;
            std::copy(coordPointer, coordPointer+_dim, coordPrt);
            MeshFormatElement e(MeshFormat::GmfVertices, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 1);
            coordPrt += _dim;
        }
    }
    _uMesh->setCoords( coordArray );

    return status;

}


void MeshFormatReader::setEdges( MEDCoupling::MEDCouplingUMesh* dimMesh1, int nbEdges)
{

    int iN[28]; // 28 - nb nodes in HEX27 (+ 1 for safety :)
    int ref;
    // read extra vertices for quadratic edges
    std::vector<int> quadNodesAtEdges( nbEdges + 1, -1 );
    if ( int nbQuadEdges = _reader.GmfStatKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtEdges))
    {
        _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtEdges);
        for ( int i = 1; i <= nbQuadEdges; ++i )
        {
            _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtEdges, &iN[0], &iN[1], &iN[2]);
            if ( iN[1] >= 1 )
                quadNodesAtEdges[ iN[0] ] = iN[2];
        }
    }
    // create edges
    _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfEdges);
    for ( int i = 1; i <= nbEdges; ++i )
    {
        _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfEdges, &iN[0], &iN[1], &ref);
        const int midN = quadNodesAtEdges[ i ];
        if ( midN > 0 )
        {

            mcIdType nodalConnPerCell[3] = {iN[0], iN[1], midN};
            backward_shift(nodalConnPerCell, 3);
            dimMesh1->insertNextCell(INTERP_KERNEL::NORM_SEG3, 3,nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfEdges, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 1 -_dim);
        }
        else
        {

            mcIdType nodalConnPerCell[2] = {iN[0], iN[1]};
            backward_shift(nodalConnPerCell, 2);
            dimMesh1->insertNextCell(INTERP_KERNEL::NORM_SEG2, 2,nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfEdges, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 1 -_dim);
        }
    }

    dimMesh1->finishInsertingCells();
    dimMesh1->sortCellsInMEDFileFrmt();
    _uMesh->setMeshAtLevel( 1 - _dim, dimMesh1 );
}

void MeshFormatReader::setTriangles(MEDCoupling::MEDCouplingUMesh* dimMesh2, int nbTria)
{
    int iN[28]; // 28 - nb nodes in HEX27 (+ 1 for safety :)
    int ref;
    // read extra vertices for quadratic triangles
    std::vector< std::vector<int> > quadNodesAtTriangles( nbTria + 1 );
    if ( int nbQuadTria = _reader.GmfStatKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtTriangles ))
    {
        _reader.GmfGotoKwd( _myCurrentFileId, MeshFormat::GmfExtraVerticesAtTriangles );
        for ( int i = 1; i <= nbQuadTria; ++i )
        {
            _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtTriangles,
                             &iN[0], &iN[1], &iN[2], &iN[3], &iN[4],
                             &iN[5]); // iN[5] - preview TRIA7
            if ( iN[0] <= nbTria )
            {
                std::vector<int>& nodes = quadNodesAtTriangles[ iN[0] ];
                nodes.insert( nodes.end(), & iN[2], & iN[5+1] );
                nodes.resize( iN[1] );
            }
        }
    }

    // create triangles

    _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfTriangles);

    for ( int i = 1; i <= nbTria; ++i )
    {
        _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfTriangles, &iN[0], &iN[1], &iN[2], &ref);
        std::vector<int>& midN = quadNodesAtTriangles[ i ];
        if ( midN.size() >= 3 )
        {

            mcIdType nodalConnPerCell[6] = {iN[0], iN[1], iN[2], midN[0], midN[1], midN[2]};
            backward_shift(nodalConnPerCell, 6);
            dimMesh2->insertNextCell(INTERP_KERNEL::NORM_TRI6, 6, nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfTriangles, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 2 -_dim);
        }
        else
        {
            mcIdType nodalConnPerCell[3] = {iN[0], iN[1], iN[2]};
            backward_shift(nodalConnPerCell, 3);
            dimMesh2->insertNextCell(INTERP_KERNEL::NORM_TRI3,3,nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfTriangles, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 2 -_dim);
        }
        if ( !midN.empty() ) MeshFormat::FreeVector( midN );
    }


}

void MeshFormatReader::setQuadrangles( MEDCoupling::MEDCouplingUMesh* dimMesh2, int nbQuad)
{
    int iN[28]; // 28 - nb nodes in HEX27 (+ 1 for safety :)
    int ref;
    // read extra vertices for quadratic quadrangles
    std::vector< std::vector<int> > quadNodesAtQuadrilaterals( nbQuad + 1 );
    if ( int nbQuadQuad = _reader.GmfStatKwd( _myCurrentFileId, MeshFormat::GmfExtraVerticesAtQuadrilaterals ))
    {
        _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtQuadrilaterals);
        for ( int i = 1; i <= nbQuadQuad; ++i )
        {
            _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtQuadrilaterals,
                             &iN[0], &iN[1], &iN[2], &iN[3], &iN[4], &iN[5], &iN[6]);
            if ( iN[0] <= nbQuad )
            {
                std::vector<int>& nodes = quadNodesAtQuadrilaterals[ iN[0] ];
                nodes.insert( nodes.end(), & iN[2], & iN[6+1] );
                nodes.resize( iN[1] );
            }
        }
    }
    // create quadrangles
    _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfQuadrilaterals);
    for ( int i = 1; i <= nbQuad; ++i )
    {
        _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfQuadrilaterals, &iN[0], &iN[1], &iN[2], &iN[3], &ref);
        std::vector<int>& midN = quadNodesAtQuadrilaterals[ i ];
        if ( midN.size() == 8-4 ) // QUAD8
        {
            mcIdType nodalConnPerCell[8] = {iN[0], iN[1], iN[2], iN[3],
                                            midN[0], midN[1], midN[2], midN[3]
                                           };
            backward_shift(nodalConnPerCell, 8);
            dimMesh2->insertNextCell(INTERP_KERNEL::NORM_QUAD8,8, nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfQuadrilaterals, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 2 -_dim);
        }
        else if ( midN.size() > 8-4 ) // QUAD9
        {
            mcIdType nodalConnPerCell[9] = {iN[0], iN[1], iN[2], iN[3],
                                            midN[0], midN[1], midN[2], midN[3], midN[4]
                                           };
            backward_shift(nodalConnPerCell, 9);
            dimMesh2->insertNextCell(INTERP_KERNEL::NORM_QUAD9,9, nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfQuadrilaterals, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 2 -_dim);
        }
        else // QUAD4
        {
            mcIdType nodalConnPerCell[4] = {iN[0], iN[1], iN[2], iN[3]};
            backward_shift(nodalConnPerCell, 4);
            dimMesh2->insertNextCell(INTERP_KERNEL::NORM_QUAD4, 4, nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfQuadrilaterals, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 2 -_dim);
        }
        if ( !midN.empty() ) MeshFormat::FreeVector( midN );
    }

}

void MeshFormatReader::setTetrahedras( MEDCoupling::MEDCouplingUMesh* dimMesh3, int nbTet)
{
    int iN[28]; // 28 - nb nodes in HEX27 (+ 1 for safety :)
    int ref;
    // read extra vertices for quadratic tetrahedra
    std::vector< std::vector<int> > quadNodesAtTetrahedra( nbTet + 1 );
    if ( int nbQuadTetra = _reader.GmfStatKwd( _myCurrentFileId, MeshFormat::GmfExtraVerticesAtTetrahedra ))
    {
        _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtTetrahedra);
        for ( int i = 1; i <= nbQuadTetra; ++i )
        {
            _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtTetrahedra,
                             &iN[0], &iN[1], &iN[2], &iN[3], &iN[4], &iN[5], &iN[6], &iN[7]);
            if ( iN[0] <= nbTet )
            {
                std::vector<int>& nodes = quadNodesAtTetrahedra[ iN[0] ];
                nodes.insert( nodes.end(), & iN[2], & iN[7+1] );
                nodes.resize( iN[1] );
            }
        }
    }
    // create tetrahedra
    _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfTetrahedra);
    for ( int i = 1; i <= nbTet; ++i )
    {
        _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfTetrahedra, &iN[0], &iN[1], &iN[2], &iN[3], &ref);
        std::vector<int>& midN = quadNodesAtTetrahedra[ i ];
        if ( midN.size() >= 10-4 ) // TETRA10
        {
            mcIdType nodalConnPerCell[10] = {iN[0], iN[2], iN[1], iN[3],
                                             midN[2], midN[1], midN[0], midN[3], midN[5], midN[4]
                                            };
            backward_shift(nodalConnPerCell, 10);
            dimMesh3->insertNextCell(INTERP_KERNEL::NORM_TETRA10, 10,nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfTetrahedra, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 3 -_dim);
        }
        else // TETRA4
        {
            mcIdType nodalConnPerCell[4] = {iN[0], iN[2], iN[1], iN[3]};
            backward_shift(nodalConnPerCell, 4);
            dimMesh3->insertNextCell(INTERP_KERNEL::NORM_TETRA4,4, nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfTetrahedra, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 3 -_dim);
        }
        if ( !midN.empty() ) MeshFormat::FreeVector( midN );
    }

}


void MeshFormatReader::setPyramids( MEDCoupling::MEDCouplingUMesh* dimMesh3, int nbPyr)
{
    int iN[28]; // 28 - nb nodes in HEX27 (+ 1 for safety :)
    int ref;
    _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfPyramids);
    for ( int i = 1; i <= nbPyr; ++i )
    {
        _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfPyramids, &iN[0], &iN[1], &iN[2], &iN[3], &iN[4], &ref);
        mcIdType nodalConnPerCell[5] = {iN[3], iN[2], iN[1], iN[0], iN[4]};
        backward_shift(nodalConnPerCell, 5);
        dimMesh3->insertNextCell(INTERP_KERNEL::NORM_PYRA5, 5,nodalConnPerCell);
        MeshFormatElement e(MeshFormat::GmfPyramids, i-1);
        _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 3 -_dim);
    }

}


void MeshFormatReader::setHexahedras( MEDCoupling::MEDCouplingUMesh* dimMesh3, int nbHex)
{
    int iN[28]; // 28 - nb nodes in HEX27 (+ 1 for safety :)
    int ref;
    // read extra vertices for quadratic hexahedra
    std::vector< std::vector<int> > quadNodesAtHexahedra( nbHex + 1 );
    if ( int nbQuadHexa = _reader.GmfStatKwd( _myCurrentFileId, MeshFormat::GmfExtraVerticesAtHexahedra ))
    {
        _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtHexahedra);
        for ( int i = 1; i <= nbQuadHexa; ++i )
        {
            _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfExtraVerticesAtHexahedra, &iN[0], &iN[1], // Hexa Id, Nb extra vertices
                             &iN[2], &iN[3], &iN[4], &iN[5],
                             &iN[6], &iN[7], &iN[8], &iN[9],
                             &iN[10], &iN[11], &iN[12], &iN[13], // HEXA20
                             &iN[14],
                             &iN[15], &iN[16], &iN[17], &iN[18],
                             &iN[19],
                             &iN[20]);                          // HEXA27
            if ( iN[0] <= nbHex )
            {
                std::vector<int>& nodes = quadNodesAtHexahedra[ iN[0] ];
                nodes.insert( nodes.end(), & iN[2], & iN[20+1] );
                nodes.resize( iN[1] );
            }
        }
    }
    // create hexhedra


    _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfHexahedra);
    for ( int i = 1; i <= nbHex; ++i )
    {
        _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfHexahedra, &iN[0], &iN[1], &iN[2], &iN[3],
                         &iN[4], &iN[5], &iN[6], &iN[7], &ref);
        std::vector<int>& midN = quadNodesAtHexahedra[ i ];
        if ( midN.size() == 20-8 ) // HEXA20
        {
            mcIdType nodalConnPerCell[20] = {iN[0], iN[3], iN[2], iN[1],
                                             iN[4], iN[7], iN[6], iN[5],
                                             midN[3], midN[2], midN[1], midN[0],
                                             midN[7], midN[6], midN[5], midN[4],
                                             midN[8], midN[11], midN[10], midN[9]
                                            };
            backward_shift(nodalConnPerCell, 20);
            dimMesh3->insertNextCell(INTERP_KERNEL::NORM_HEXA20,20, nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfHexahedra, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 3 -_dim);


        }
        else if ( midN.size() >= 27-8 ) // HEXA27
        {
            mcIdType nodalConnPerCell[27] = {iN[0], iN[3], iN[2], iN[1],
                                             iN[4], iN[7], iN[6], iN[5],
                                             midN[3], midN[2], midN[1], midN[0],
                                             midN[7], midN[6], midN[5], midN[4],
                                             midN[8], midN[11], midN[10], midN[9],
                                             midN[12],
                                             midN[16], midN[15], midN[14], midN[13],
                                             midN[17],
                                             midN[18]
                                            };
            backward_shift(nodalConnPerCell, 27);
            dimMesh3->insertNextCell(INTERP_KERNEL::NORM_HEXA27,27, nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfHexahedra, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 3 -_dim);

        }
        else // HEXA8
        {
            mcIdType nodalConnPerCell[8] = {iN[0], iN[3], iN[2], iN[1],
                                            iN[4], iN[7], iN[6], iN[5]
                                           };
            backward_shift(nodalConnPerCell, 8);
            dimMesh3->insertNextCell(INTERP_KERNEL::NORM_HEXA8,8, nodalConnPerCell);
            MeshFormatElement e(MeshFormat::GmfHexahedra, i-1);
            _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 3 -_dim);

        }
        if ( !midN.empty() ) MeshFormat::FreeVector( midN );
    }


}



void MeshFormatReader::setPrisms(MEDCoupling::MEDCouplingUMesh* dimMesh3, int nbPrism)
{
    int iN[28]; // 28 - nb nodes in HEX27 (+ 1 for safety :)
    int ref;
    _reader.GmfGotoKwd(_myCurrentFileId, MeshFormat::GmfPrisms);
    for ( int i = 1; i <= nbPrism; ++i )
    {
        _reader.GmfGetLin(_myCurrentFileId, MeshFormat::GmfPrisms, &iN[0], &iN[1], &iN[2], &iN[3], &iN[4], &iN[5], &ref);
        mcIdType nodalConnPerCell[8] = {iN[0], iN[2], iN[1], iN[3], iN[5], iN[4]};
        backward_shift(nodalConnPerCell, 8);
        dimMesh3->insertNextCell(INTERP_KERNEL::NORM_PENTA6, 6,nodalConnPerCell);
        MeshFormatElement e(MeshFormat::GmfPrisms, i-1);
        _fams.insert(std::pair <int, MeshFormatElement> (ref, e), 3 -_dim);

    }


}

void MeshFormatReader::callParserGetLin( MeshFormat::GmfKwdCod kwd,  double* val, int valSize, int* ref)
{
    switch(valSize)
    {
    case 1:
    {
        _reader.GmfGetLin(_myCurrentFileId, kwd, &val[0], ref);
        break;
    }
    case 2:
    {
        _reader.GmfGetLin(_myCurrentFileId, kwd, &val[0], &val[1], ref);
        break;
    }
    case 3:
    {
        _reader.GmfGetLin(_myCurrentFileId, kwd, &val[0], &val[1], &val[2], ref);
        break;
    }
    case 4:
    {
        _reader.GmfGetLin(_myCurrentFileId, kwd, &val[0], &val[1], &val[2], &val[3], ref);
        break;
    }
    case 6:
    {
        _reader.GmfGetLin(_myCurrentFileId, kwd, &val[0], &val[1], &val[2], &val[3], &val[4], &val[5], ref);
        break;
    }
    case 9:
    {
        _reader.GmfGetLin(_myCurrentFileId, kwd, &val[0], &val[1], &val[2], &val[3], &val[4], &val[5], &val[6], &val[7], &val[8], ref);
        break;
    }
    }
}

void MeshFormatReader::setTypeOfFieldAndDimRel(MeshFormat::GmfKwdCod kwd, MEDCoupling::TypeOfField* typeOfField, int* dimRel )
{
    switch (kwd)
    {
    case MeshFormat::GmfSolAtVertices :
    {
        *typeOfField = MEDCoupling::ON_NODES;
        *dimRel = 1 ;
        break;
    }
    case MeshFormat::GmfSolAtEdges :
    {
        *typeOfField = MEDCoupling::ON_CELLS;
        *dimRel = 1 - _dim ;
        break;
    }
    case MeshFormat::GmfSolAtTriangles :
    case MeshFormat::GmfSolAtQuadrilaterals :
    {
        *typeOfField = MEDCoupling::ON_CELLS;
        *dimRel = 2 - _dim;
        break;
    }
    case MeshFormat::GmfSolAtTetrahedra :
    case MeshFormat::GmfSolAtPrisms :
    case MeshFormat::GmfSolAtHexahedra :
    {
        *typeOfField = MEDCoupling::ON_CELLS;
        *dimRel = 3 - _dim;
        break;
    }
    }
}


INTERP_KERNEL::NormalizedCellType MeshFormatReader::toMedType(MeshFormat::GmfKwdCod kwd)
{
    INTERP_KERNEL::NormalizedCellType type;
    //~switch (kwd)
    //~{
    //~case MeshFormat::GmfEdges :
    //~{

    //~type = INTERP_KERNEL::NORM_SEG2;
    //~break;
    //~}
    //~case MeshFormat::GmfTriangles :
    //~{
    //~type = INTERP_KERNEL::NORM_TRI3;
    //~break;
    //~}
    //~case MeshFormat::GmfQuadrilaterals :
    //~{
    //~type = INTERP_KERNEL::NORM_QUAD;
    //~break;
    //~}
    //~case MeshFormat::GmfSolAtTetrahedra :
    //~case MeshFormat::GmfSolAtPrisms :
    //~case MeshFormat::GmfSolAtHexahedra :
    //~{
    //~*typeOfField = MEDCoupling::ON_CELLS;
    //~*dimRel = 3 - _dim;
    //~break;
    //~}
    //~}
    return type;
}


void MeshFormatReader::buildFamilies()
{
    buildNodesFamilies();
    buildCellsFamilies();
}

void MeshFormatReader::buildCellsFamilies()
{
    std::vector<int> levs =  _uMesh->getNonEmptyLevels();
    for (size_t iDim = 0; iDim<levs.size(); iDim++ )
    {
        int dimRelMax = levs[iDim];
        std::map <int, std::vector<MeshFormatElement>* > famDim = _fams.getMapAtLevel(dimRelMax);
        std::map <int, std::vector<MeshFormatElement>* >::const_iterator _meshFormatFamsIt = famDim.begin();
        std::vector< const MEDCoupling::DataArrayIdType* > fams;
        MEDCoupling::DataArrayIdType* cellIds = MEDCoupling::DataArrayIdType::New();
        cellIds->alloc(_uMesh->getSizeAtLevel(dimRelMax), 1);
        cellIds->fillWithZero();

        for(; _meshFormatFamsIt!= famDim.end(); ++_meshFormatFamsIt)
        {
            const int famId = _meshFormatFamsIt->first;
            std::string famName ="FromMeshGemsFormatAttributFamily_"+std::to_string(famId);
            std::vector <MeshFormatElement>* cellsInFam = _meshFormatFamsIt->second;
            if (!famId) continue;
            std::vector <MeshFormatElement>::iterator cellsInFamIt = cellsInFam->begin();

            _uMesh->addFamily(famName, famId);
            for ( ; cellsInFamIt !=cellsInFam->end(); ++cellsInFamIt)
            {
                cellIds->setIJ(cellsInFamIt->_id, 0, famId);
            }

        }
        _uMesh->setFamilyFieldArr(dimRelMax, cellIds->deepCopy());
        cellIds->decrRef();


    }  
}

void MeshFormatReader::buildNodesFamilies()
{
  std::vector<int> levs =  _uMesh->getNonEmptyLevels();
     int dimRelMax = 1;
  std::map <int, std::vector<MeshFormatElement>* > famDim = _fams.getMapAtLevel(dimRelMax);
  std::map <int, std::vector<MeshFormatElement>* >::const_iterator _meshFormatFamsIt = famDim.begin();
  std::vector< const MEDCoupling::DataArrayIdType* > fams;
  MEDCoupling::DataArrayIdType* cellIds = MEDCoupling::DataArrayIdType::New();
  cellIds->alloc(_uMesh->getSizeAtLevel(dimRelMax), 1);
  cellIds->fillWithZero();

  for(; _meshFormatFamsIt!= famDim.end(); ++_meshFormatFamsIt)
  {
    const int famId = _meshFormatFamsIt->first;
    if (!famId) continue; 
    bool thisIsACellFamily = false;
    
    for (size_t iDim = 0; iDim<levs.size(); iDim++ )
    {
      int dimMesh = levs[iDim];
      std::map <int, std::vector<MeshFormatElement>* > famDimAtLevel = _fams.getMapAtLevel(dimMesh);  
      std::map <int, std::vector<MeshFormatElement>* >::iterator famDimAtLevelId = famDimAtLevel.find(famId);  
      if (famDimAtLevelId != famDimAtLevel.end())
      {
        thisIsACellFamily = true;
        break;
      }
      
    }
    
    if (thisIsACellFamily) continue;
    std::string famName ="FromMeshGemsFormatAttributFamily_"+std::to_string(famId);
    std::vector <MeshFormatElement>* cellsInFam = _meshFormatFamsIt->second;
    std::vector <MeshFormatElement>::iterator cellsInFamIt = cellsInFam->begin();

    _uMesh->addFamily(famName, famId);
    for ( ; cellsInFamIt !=cellsInFam->end(); ++cellsInFamIt)
    {
      cellIds->setIJ(cellsInFamIt->_id, 0, famId);
    }

  }
  _uMesh->setFamilyFieldArr(dimRelMax, cellIds->deepCopy());
  cellIds->decrRef();
}
}
