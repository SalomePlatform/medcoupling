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

#ifndef __INTERPKERNELGAUSSCOORDS_HXX__
#define __INTERPKERNELGAUSSCOORDS_HXX__

#include "INTERPKERNELDefines.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <vector>

namespace INTERP_KERNEL 
{
  typedef std::vector<double> DataVector;
  typedef std::vector<int>    IndexVector;

  //Class to store Gauss Points information
  class GaussInfo 
  {
  public:
    INTERPKERNEL_EXPORT GaussInfo( NormalizedCellType theGeometry,
                                   const DataVector& theGaussCoord,
                                   int theNbGauss,
                                   const DataVector& theReferenceCoord,
                                   int theNbRef
                                   );
    INTERPKERNEL_EXPORT ~GaussInfo();

    INTERPKERNEL_EXPORT NormalizedCellType getCellType() const;    

    INTERPKERNEL_EXPORT int getGaussCoordDim() const;
    INTERPKERNEL_EXPORT int getReferenceCoordDim() const;
    INTERPKERNEL_EXPORT DataVector getGaussCoords() const { return _my_gauss_coord; }
    INTERPKERNEL_EXPORT DataVector getRefCoords() const { return _my_reference_coord; }
    INTERPKERNEL_EXPORT NormalizedCellType getGeoType() const { return _my_geometry; }

    INTERPKERNEL_EXPORT int getNbGauss() const;
    INTERPKERNEL_EXPORT int getNbRef() const;

    INTERPKERNEL_EXPORT GaussInfo convertToLinear() const;

    INTERPKERNEL_EXPORT const double* getFunctionValues( const int theGaussId ) const;

    INTERPKERNEL_EXPORT void initLocalInfo();
    
    INTERPKERNEL_EXPORT static std::vector<double> NormalizeCoordinatesIfNecessary(NormalizedCellType ct, int inputDim, const std::vector<double>& inputArray);

  public:
    static const double SEG2_REF[2];
    static const double SEG3_REF[3];
    static const double TRIA3A_REF[6];
    static const double TRIA3B_REF[6];
    static const double TRIA6A_REF[12];
    static const double TRIA6B_REF[12];
    static const double TRIA7A_REF[14];
    static const double QUAD4A_REF[8];
    static const double QUAD4B_REF[8];
    static const double QUAD8A_REF[16];
    static const double QUAD8B_REF[16];
    static const double QUAD9A_REF[18];
    static const double TETRA4A_REF[12];
    static const double TETRA4B_REF[12];
    static const double TETRA10A_REF[30];
    static const double TETRA10B_REF[30];
    static const double PYRA5A_REF[15];
    static const double PYRA5B_REF[15];
    static const double PYRA13A_REF[39];
    static const double PYRA13B_REF[39];
    static const double PENTA6A_REF[18];
    static const double PENTA6B_REF[18];
    static const double PENTA15A_REF[45];
    static const double PENTA15B_REF[45];
    static const double HEXA8A_REF[24];
    static const double HEXA8B_REF[24];
    static const double HEXA20A_REF[60];
    static const double HEXA20B_REF[60];
    static const double HEXA27A_REF[81];
  protected:
    static bool IsSatisfy(const std::vector<double>& ref1, const std::vector<double>& ref2);
    bool isSatisfy();
    
    void point1Init();
    
    //1D
    void seg2Init();
    void seg3Init();

    //2D
    void tria3aInit();
    void tria3bInit();
    void tria6aInit();
    void tria6bInit();
    void tria7aInit();

    void quad4aInit();
    static void Quad4aInit(GaussInfo& obj) { obj.quad4aInit(); }
    void quad4bInit();
    static void Quad4bInit(GaussInfo& obj) { obj.quad4bInit(); }
    void quad4cInit();
    static void Quad4cInit(GaussInfo& obj) { obj.quad4cInit(); }
    void quad4DegSeg2Init();
    static void Quad4DegSeg2Init(GaussInfo& obj) { obj.quad4DegSeg2Init(); }
    void quad8aInit();
    void quad8bInit();
    void quad9aInit();

    //3D
    void tetra4aInit();
    void tetra4bInit();
    void tetra10aInit();
    void tetra10bInit();

    void pyra5aInit();
    void pyra5bInit();
    void pyra13aInit();
    void pyra13bInit();

    void penta6aInit();
    static void Penta6aInit(GaussInfo& obj) { obj.penta6aInit(); }
    void penta6bInit();
    static void Penta6bInit(GaussInfo& obj) { obj.penta6bInit(); }
    void penta6DegTria3aInit();
    static void Penta6DegTria3aInit(GaussInfo& obj) { obj.penta6DegTria3aInit(); }
    void penta6DegTria3bInit();
    static void Penta6DegTria3bInit(GaussInfo& obj) { obj.penta6DegTria3bInit(); }
    
    void penta15aInit();
    static void Penta15aInit(GaussInfo& obj) { obj.penta15aInit(); }
    void penta15bInit();
    static void Penta15bInit(GaussInfo& obj) { obj.penta15bInit(); }

    void hexa8aInit();
    static void Hexa8aInit(GaussInfo& obj) { obj.hexa8aInit(); }
    void hexa8bInit();
    static void Hexa8bInit(GaussInfo& obj) { obj.hexa8bInit(); }
    void hexa8DegQuad4aInit();
    static void Hexa8DegQuad4aInit(GaussInfo& obj) { obj.hexa8DegQuad4aInit(); }
    void hexa8DegQuad4bInit();
    static void Hexa8DegQuad4bInit(GaussInfo& obj) { obj.hexa8DegQuad4bInit(); }
    void hexa8DegQuad4cInit();
    static void Hexa8DegQuad4cInit(GaussInfo& obj) { obj.hexa8DegQuad4cInit(); }
    void hexa20aInit();
    void hexa20bInit();
    void hexa27aInit();

  private:
    //INFORMATION from MEDMEM
    NormalizedCellType _my_geometry;               //Cell type

    int                _my_nb_gauss;                //Nb of the gauss points for element
    DataVector         _my_gauss_coord;             //Gauss coordinates

    int                _my_nb_ref;                  //Nb of the nodes for element:
                                                 //NORM_SEG2 - 2
                                                 //NORM_SEG3 - 3
                                                 //NORM_TRI3 - 3
                                                 //.............

    DataVector         _my_reference_coord;         //Reference coordinates

    //LOCAL INFORMATION
    DataVector         _my_local_reference_coord;    //Vector to store reference coordinates
    int                _my_local_ref_dim;            //Dimension of the local reference coordinates:
                                                 // (x)       - 1D case
                                                 // (x, y)    - 2D case
                                                 // (x, y, z) - 3D case
    int                _my_local_nb_ref;             //Nb of the local reference coordinates

    DataVector         _my_function_value;          //Shape Function values
  };


  //Class for calculation of the coordinates of the gauss points 
  class GaussCoords 
  {
  public:

    INTERPKERNEL_EXPORT GaussCoords();
    INTERPKERNEL_EXPORT ~GaussCoords();

    INTERPKERNEL_EXPORT void addGaussInfo( NormalizedCellType theGeometry,
                                           int coordDim,
                                           const double* theGaussCoord,
                                           int theNbGauss,
                                           const double* theReferenceCoord,
                                           int theNbRef);

    INTERPKERNEL_EXPORT double* calculateCoords( NormalizedCellType theGeometry, 
                                                 const double* theNodeCoords, 
                                                 const int theSpaceDim,
                                                 const int* theIndex);

    INTERPKERNEL_EXPORT void calculateCoords( NormalizedCellType theGeometry, 
                                              const double* theNodeCoords, 
                                              const int theSpaceDim,
                                              const int* theIndex,
                                              double *result);
  private:
    const GaussInfo *getInfoGivenCellType(NormalizedCellType cellType);
    void calculateCoordsAlg(const GaussInfo *info, const double* theNodeCoords, const int theSpaceDim, const int *theIndex,
                            double *result);
  private:
    typedef std::vector<GaussInfo*> GaussInfoVector;
    GaussInfoVector _my_gauss_info;
  };
}
#endif //INTERPKERNELGAUSSCOORDS
