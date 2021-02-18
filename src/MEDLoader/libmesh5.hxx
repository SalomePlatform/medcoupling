// Copyright (C) 2021  CEA/DEN, EDF R&D
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

#ifndef MESHFORMATPARSER_HXX
#define MESHFORMATPARSER_HXX
/*----------------------------------------------------------*/
/*                                                                                                                      */
/*                                              LIBMESH V 5.46                                          */
/*                                                                                                                      */
/*----------------------------------------------------------*/
/*                                                                                                                      */
/*      Description:            handle .meshb file format I/O           */
/*      Author:                         Loic MARECHAL                                           */
/*      Creation date:          feb 16 2007                                                     */
/*      Last modification:      dec 12 2020                                                     */
/*                                                                                                                      */
/*----------------------------------------------------------*/


/*----------------------------------------------------------*/
/* Defines                                                                                                      */
/*----------------------------------------------------------*/


#define GmfStrSiz 1024
#define GmfMaxTyp 1000
#define GmfMaxKwd 81
#define GmfMshVer 1
#define GmfRead 1
#define GmfWrite 2
#define GmfSca 1
#define GmfVec 2
#define GmfSymMat 3
#define GmfMat 4
#define GmfFloat 1
#define GmfDouble 2


/*----------------------------------------------------------*/
/* Defines                                                                                                      */
/*----------------------------------------------------------*/

#define Asc 1
#define Bin 2
#define MshFil 4
#define SolFil 8
#define MaxMsh 100
#define InfKwd 1
#define RegKwd 2
#define SolKwd 3
#define WrdSiz 4
#define BufSiz 10000


// see MeshGems/Docs/meshgems_formats_description.pdf
extern const char* GmfKwdFmt[ GmfMaxKwd + 1 ][4];
/*----------------------------------------------------------*/
/* Structures                                                                                           */
/*----------------------------------------------------------*/
namespace  MeshFormat {
typedef struct
{
    int typ, SolSiz, NmbWrd, NmbLin, NmbTyp, TypTab[ GmfMaxTyp ];
    long pos;
    char fmt[ GmfMaxTyp*9 ];
} KwdSct;

typedef struct
{
    int dim, ver, mod, typ, cod, pos;
    long NexKwdPos, siz;
    KwdSct KwdTab[ GmfMaxKwd + 1 ];
    FILE *hdl;
    int *IntBuf;
    float *FltBuf;
    unsigned char *buf;
    char FilNam[ GmfStrSiz ];
    double DblBuf[1000/8];
    unsigned char blk[ BufSiz + 1000 ];
} GmfMshSct;





// see MeshGems/Docs/meshgems_formats_description.pdf
enum GmfKwdCod
{
    GmfReserved1, \
    GmfVersionFormatted, \
    GmfReserved2, \
    GmfDimension, \
    GmfVertices, \
    GmfEdges, \
    GmfTriangles, \
    GmfQuadrilaterals, \
    GmfTetrahedra, \
    GmfPrisms, \
    GmfHexahedra, \
    GmfIterationsAll, \
    GmfTimesAll, \
    GmfCorners, \
    GmfRidges, \
    GmfRequiredVertices, \
    GmfRequiredEdges, \
    GmfRequiredTriangles, \
    GmfRequiredQuadrilaterals, \
    GmfTangentAtEdgeVertices, \
    GmfNormalAtVertices, \
    GmfNormalAtTriangleVertices, \
    GmfNormalAtQuadrilateralVertices, \
    GmfAngleOfCornerBound, \
    GmfTrianglesP2, \
    GmfEdgesP2, \
    GmfSolAtPyramids, \
    GmfQuadrilateralsQ2, \
    GmfISolAtPyramids, \
    GmfSubDomainFromGeom, \
    GmfTetrahedraP2, \
    GmfFault_NearTri, \
    GmfFault_Inter, \
    GmfHexahedraQ2, \
    GmfExtraVerticesAtEdges, \
    GmfExtraVerticesAtTriangles, \
    GmfExtraVerticesAtQuadrilaterals, \
    GmfExtraVerticesAtTetrahedra, \
    GmfExtraVerticesAtPrisms, \
    GmfExtraVerticesAtHexahedra, \
    GmfVerticesOnGeometricVertices, \
    GmfVerticesOnGeometricEdges, \
    GmfVerticesOnGeometricTriangles, \
    GmfVerticesOnGeometricQuadrilaterals, \
    GmfEdgesOnGeometricEdges, \
    GmfFault_FreeEdge, \
    GmfPolyhedra, \
    GmfPolygons, \
    GmfFault_Overlap, \
    GmfPyramids, \
    GmfBoundingBox, \
    GmfBody, \
    GmfPrivateTable, \
    GmfFault_BadShape, \
    GmfEnd, \
    GmfTrianglesOnGeometricTriangles, \
    GmfTrianglesOnGeometricQuadrilaterals, \
    GmfQuadrilateralsOnGeometricTriangles, \
    GmfQuadrilateralsOnGeometricQuadrilaterals, \
    GmfTangents, \
    GmfNormals, \
    GmfTangentAtVertices, \
    GmfSolAtVertices, \
    GmfSolAtEdges, \
    GmfSolAtTriangles, \
    GmfSolAtQuadrilaterals, \
    GmfSolAtTetrahedra, \
    GmfSolAtPrisms, \
    GmfSolAtHexahedra, \
    GmfDSolAtVertices, \
    GmfISolAtVertices, \
    GmfISolAtEdges, \
    GmfISolAtTriangles, \
    GmfISolAtQuadrilaterals, \
    GmfISolAtTetrahedra, \
    GmfISolAtPrisms, \
    GmfISolAtHexahedra, \
    GmfIterations, \
    GmfTime, \
    GmfFault_SmallTri, \
    GmfCoarseHexahedra, \
    GmfFault_MultipleEdge
};



class MeshFormatParser {
    /*----------------------------------------------------------*/
    /* External procedures                                                                          */
    /*----------------------------------------------------------*/
public :
    MeshFormatParser();
    int GmfOpenMesh(const char *, int, ...);
    int GmfCloseMesh(int);
    int GmfStatKwd(int, int, ...);
    int GmfGotoKwd(int, int);
    int GmfSetKwd(int, int, ...);
    void GmfGetLin(int, int, ...);
    void GmfSetLin(int, int, ...);
private :


    /*----------------------------------------------------------*/
    /*  private procedures methods                                                       */
    /*----------------------------------------------------------*/

    void ScaWrd(GmfMshSct *, unsigned char *);
    void ScaDblWrd(GmfMshSct *, unsigned char *);
    void ScaBlk(GmfMshSct *, unsigned char *, int);
    long GetPos(GmfMshSct *);
    void RecWrd(GmfMshSct *, unsigned char *);
    void RecDblWrd(GmfMshSct *, unsigned char *);
    void RecBlk(GmfMshSct *, unsigned char *, int);
    void SetPos(GmfMshSct *, long);
    int ScaKwdTab(GmfMshSct *);
    void ExpFmt(GmfMshSct *, int);
    void ScaKwdHdr(GmfMshSct *, int);

    void GmfCpyLin(int, int, int);


    int GmfIniFlg;
    GmfMshSct *GmfMshTab[ MaxMsh + 1 ];





};

}
#endif // MESHFORMATPARSER_HXX
