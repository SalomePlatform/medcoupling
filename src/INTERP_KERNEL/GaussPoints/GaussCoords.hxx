//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#ifndef __INTERPKERNELGAUSS_HXX__
#define __INTERPKERNELGAUSS_HXX__

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
    GaussInfo( NormalizedCellType theGeometry,
               const DataVector& theGaussCoord,
               int theNbGauss,
               const DataVector& theReferenceCoord,
               int theNbRef
               );
    ~GaussInfo();

    NormalizedCellType GetCellType() const;    

    int GetGaussCoordDim() const;
    int GetReferenceCoordDim() const;

    int GetNbGauss() const;
    int GetNbRef() const;

    const double* GetFunctionValues( const int theGaussId ) const;

    void InitLocalInfo() throw (INTERP_KERNEL::Exception);

  protected:

    bool isSatisfy();

    //1D
    void Seg2Init();
    void Seg3Init();

    //2D
    void Tria3aInit();
    void Tria3bInit();
    void Tria6aInit();
    void Tria6bInit();

    void Quad4aInit();
    void Quad4bInit();
    void Quad8aInit();
    void Quad8bInit();

    //3D
    void Tetra4aInit();
    void Tetra4bInit();
    void Tetra10aInit();
    void Tetra10bInit();

    void Pyra5aInit();
    void Pyra5bInit();
    void Pyra13aInit();
    void Pyra13bInit();

    void Penta6aInit();
    void Penta6bInit();
    void Penta15aInit();
    void Penta15bInit();

    void Hexa8aInit();
    void Hexa8bInit();
    void Hexa20aInit();
    void Hexa20bInit();


  private:
    //INFORMATION from MEDMEM
    NormalizedCellType myGeometry;               //Cell type

    int                myNbGauss;                //Nb of the gauss points for element
    DataVector         myGaussCoord;             //Gauss coordinates

    int                myNbRef;                  //Nb of the nodes for element:
                                                 //NORM_SEG2 - 2
                                                 //NORM_SEG3 - 3
                                                 //NORM_TRI3 - 3
                                                 //.............

    DataVector         myReferenceCoord;         //Reference coordinates

    //LOCAL INFORMATION
    DataVector         myLocalReferenceCoord;    //Vector to store reference coordinates
    int                myLocalRefDim;            //Dimension of the local reference coordinates:
                                                 // (x)       - 1D case
                                                 // (x, y)    - 2D case
                                                 // (x, y, z) - 3D case
    int                myLocalNbRef;             //Nb of the local reference coordinates

    DataVector         myFunctionValue;          //Shape Function values
  };


  //Class for calculation of the coordinates of the gauss points 
  class GaussCoords 
  {
  public:

    GaussCoords();
    ~GaussCoords();

    void addGaussInfo( NormalizedCellType theGeometry,
                       int coordDim,
                       const double* theGaussCoord,
                       int theNbGauss,
                       const double* theReferenceCoord,
                       int theNbRef) throw (INTERP_KERNEL::Exception);

    double* CalculateCoords( NormalizedCellType theGeometry, 
                             const double* theNodeCoords, 
                             const int theDimSpace,
                             const int* theIndex
                             ) throw(INTERP_KERNEL::Exception);
  private:
    typedef std::vector<GaussInfo*> GaussInfoVector;
    GaussInfoVector myGaussInfo;
  };
}
#endif //INTERPKERNELGAUSS
