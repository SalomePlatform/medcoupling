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

    INTERPKERNEL_EXPORT int getNbGauss() const;
    INTERPKERNEL_EXPORT int getNbRef() const;

    INTERPKERNEL_EXPORT const double* getFunctionValues( const int theGaussId ) const;

    INTERPKERNEL_EXPORT void initLocalInfo();

  protected:

    bool isSatisfy();

    //1D
    void seg2Init();
    void seg3Init();

    //2D
    void tria3aInit();
    void tria3bInit();
    void tria6aInit();
    void tria6bInit();

    void quad4aInit();
    void quad4bInit();
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
    void penta6bInit();
    void penta15aInit();
    void penta15bInit();

    void hexa8aInit();
    void hexa8bInit();
    void hexa20aInit();
    void hexa20bInit();


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
