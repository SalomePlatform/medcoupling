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

//Local includes
#include "InterpKernelGaussCoords.hxx"
#include "CellModel.hxx"

//STL includes
#include <math.h>
#include <algorithm>
#include <sstream>

using namespace INTERP_KERNEL;

const double GaussInfo::SEG2_REF[2]={-1.0, 1.0};

const double GaussInfo::SEG3_REF[3]={-1.0, 1.0, 0.0};

const double GaussInfo::TRIA3A_REF[6]={-1.0, 1.0, -1.0, -1.0, 1.0, -1.0};

const double GaussInfo::TRIA3B_REF[6]={0.0, 0.0, 1.0, 0.0, 0.0, 1.0};

const double GaussInfo::TRIA6A_REF[12]={-1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 0.0};

const double GaussInfo::TRIA6B_REF[12]={0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5};

const double GaussInfo::TRIA7A_REF[14]={0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.3333333333333333, 0.3333333333333333};

const double GaussInfo::QUAD4A_REF[8]={-1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0};

const double GaussInfo::QUAD4B_REF[8]={-1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0};

const double GaussInfo::QUAD8A_REF[16]={-1.0, 1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0, 1.0};

const double GaussInfo::QUAD8B_REF[16]={-1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0};

const double GaussInfo::QUAD9A_REF[18]={-1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0};

const double GaussInfo::TETRA4A_REF[12]={0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0};

const double GaussInfo::TETRA4B_REF[12]={0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};

const double GaussInfo::TETRA10A_REF[30]={0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5, 0.0, 0.0};

const double GaussInfo::TETRA10B_REF[30]={0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.5};

const double GaussInfo::PYRA5A_REF[15]={1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0};

const double GaussInfo::PYRA5B_REF[15]={1.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

const double GaussInfo::PYRA13A_REF[39]={1.0, 0.0, 0.0, 0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.5, 0.0, -0.5, 0.5, 0.0, -0.5, -0.5, 0.0, 0.5, -0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5, -0.5, 0.0, 0.5, 0.0, -0.5, 0.5};

const double GaussInfo::PYRA13B_REF[39]={1.0, 0.0, 0.0, 0.0, -1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.5, -0.5, 0.0, -0.5, -0.5, 0.0, -0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, -0.5, 0.5, -0.5, 0.0, 0.5, 0.0, 0.5, 0.5};

const double GaussInfo::PENTA6A_REF[18]={-1.0, 1.0, 0.0, -1.0, -0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0};

const double GaussInfo::PENTA6B_REF[18]={-1.0, 1.0, 0.0, -1.0, 0.0, 0.0, -1.0, -0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0};

const double GaussInfo::PENTA15A_REF[45]={-1.0, 1.0, 0.0, -1.0, -0.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 0.0, 0.0, -1.0, 0.5, 0.5, -1.0, 0.0, 0.5, -1.0, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.5, 0.5, 1.0, 0.0, 0.5, 1.0, 0.5, 0.0};

const double GaussInfo::PENTA15B_REF[45]={-1.0, 1.0, 0.0, -1.0, 0.0, 0.0, -1.0, -0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, -1.0, 0.5, 0.0, -1.0, 0.0, 0.5, -1.0, 0.5, 0.5, 1.0, 0.5, 0.0, 1.0, 0.0, 0.5, 1.0, 0.5, 0.5, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0};

const double GaussInfo::HEXA8A_REF[24]={-1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0};

const double GaussInfo::HEXA8B_REF[24]={-1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0};

const double GaussInfo::HEXA20A_REF[60]={-1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, -1.0, -1.0, -1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, 1.0, 0.0, -1.0, -1.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, -1.0, 0.0, 1.0, 1.0, 0.0, -1.0, 1.0, 0.0, 0.0, -1.0, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0, 1.0, -1.0, 0.0, 1.0};

const double GaussInfo::HEXA20B_REF[60]={-1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 0.0, -1.0, 0.0, 1.0, -1.0, 1.0, 0.0, -1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, -1.0, 1.0, -1.0, -1.0, 0.0, -1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, -1.0, 0.0};

const double GaussInfo::HEXA27A_REF[81]={-1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 1.0, -1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 0.0, -1.0, 0.0, 1.0, -1.0, 1.0, 0.0, -1.0, 0.0, -1.0, -1.0, -1.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 0.0, -1.0, 1.0, -1.0, -1.0, 0.0, -1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, -1.0, 0.0, 0.0, 0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0};

//Define common part of the code in the MACRO
//---------------------------------------------------------------
#define LOCAL_COORD_MACRO_BEGIN                                         \
  _my_local_reference_coord.resize( _my_local_ref_dim*_my_local_nb_ref );           \
  for( int refId = 0; refId < _my_local_nb_ref; refId++ )                   \
    {                                                                   \
      double* coords = &_my_local_reference_coord[ refId*_my_local_ref_dim ];   \
      switch(refId)                                                     \
        {

//---------------------------------------------------------------
#define LOCAL_COORD_MACRO_END                   \
  }                                             \
}

//---------------------------------------------------------------
#define SHAPE_FUN_MACRO_BEGIN                                           \
  for( int gaussId     = 0 ; gaussId < _my_nb_gauss ; gaussId++ )          \
    {                                                                   \
      double* funValue =  &_my_function_value[ gaussId * _my_nb_ref ];        \
      const double* gc = &_my_gauss_coord[ gaussId * getGaussCoordDim() ];

//---------------------------------------------------------------
#define SHAPE_FUN_MACRO_END                     \
  }

#define CHECK_MACRO                                                        \
  if( ! aSatify )                                                          \
    {                                                                      \
      std::ostringstream stream;                                           \
      stream << "Error in the gauss localization for the cell with type "; \
      stream << cellModel.getRepr();                                       \
      stream << " !!!";                                                    \
      throw INTERP_KERNEL::Exception(stream.str().c_str());                \
    }


//---------------------------------------------------------------
static bool IsEqual(double theLeft, double theRight) 
{
  static double EPS = 1.0E-3;
  if(fabs(theLeft) + fabs(theRight) > EPS)
    return fabs(theLeft-theRight)/(fabs(theLeft)+fabs(theRight)) < EPS;
  return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////
//                                GAUSS INFO CLASS                                            //
////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * Constructor of the GaussInfo
 */
GaussInfo::GaussInfo( NormalizedCellType theGeometry,
                      const DataVector& theGaussCoord,
                      int theNbGauss,
                      const DataVector& theReferenceCoord,
                      int theNbRef ) :
  _my_geometry(theGeometry),
  _my_nb_gauss(theNbGauss),
  _my_gauss_coord(theGaussCoord),
  _my_nb_ref(theNbRef),
  _my_reference_coord(theReferenceCoord)
{

  //Allocate shape function values
  _my_function_value.resize( _my_nb_gauss * _my_nb_ref );
}

/*!
 * Destructor
 */
GaussInfo::~GaussInfo()
{
}

/*!
 * Return dimension of the gauss coordinates
 */
int GaussInfo::getGaussCoordDim() const 
{
  if( _my_nb_gauss ) 
    {
      return _my_gauss_coord.size()/_my_nb_gauss;
    }
  else 
    {
      return 0;
    }
}

/*!
 * Return dimension of the reference coordinates
 */
int GaussInfo::getReferenceCoordDim() const 
{
  if( _my_nb_ref ) 
    {
      return _my_reference_coord.size()/_my_nb_ref;
    }
  else 
    {
      return 0;
    }
}

/*!
 * Return type of the cell.
 */
NormalizedCellType GaussInfo::getCellType() const 
{
  return _my_geometry;
}

/*!
 * Return Nb of the gauss points.
 */
int GaussInfo::getNbGauss() const 
{
  return _my_nb_gauss;
}

/*!
 * Return Nb of the reference coordinates.
 */
int GaussInfo::getNbRef() const 
{
  return _my_nb_ref;
}

GaussInfo GaussInfo::convertToLinear() const
{
  switch(_my_geometry)
    {
    case NORM_SEG3:
      {
        std::vector<double> a(SEG3_REF,SEG3_REF+3);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(SEG2_REF,SEG2_REF+2);
            return GaussInfo(NORM_SEG2,_my_gauss_coord,getNbGauss(),c,2);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for SEG3 !");
      }
    case NORM_TRI6:
      {
        std::vector<double> a(TRIA6A_REF,TRIA6A_REF+12),b(TRIA6B_REF,TRIA6B_REF+12);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(TRIA3A_REF,TRIA3A_REF+6);
            return GaussInfo(NORM_TRI3,_my_gauss_coord,getNbGauss(),c,3);
          }
        if(IsSatisfy(b,_my_reference_coord))
          {
            std::vector<double> c(TRIA3B_REF,TRIA3B_REF+6);
            return GaussInfo(NORM_TRI3,_my_gauss_coord,getNbGauss(),c,3);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for TRI6 !");
      }
    case NORM_TRI7:
      {
        std::vector<double> a(TRIA7A_REF,TRIA7A_REF+14);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(TRIA3B_REF,TRIA3B_REF+6);
            return GaussInfo(NORM_TRI3,_my_gauss_coord,getNbGauss(),c,3);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for TRI7 !");
      }
    case NORM_QUAD8:
      {
        std::vector<double> a(QUAD8A_REF,QUAD8A_REF+16),b(QUAD8B_REF,QUAD8B_REF+16);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(QUAD4A_REF,QUAD4A_REF+8);
            return GaussInfo(NORM_QUAD4,_my_gauss_coord,getNbGauss(),c,4);
          }
        if(IsSatisfy(b,_my_reference_coord))
          {
            std::vector<double> c(QUAD4B_REF,QUAD4B_REF+8);
            return GaussInfo(NORM_QUAD4,_my_gauss_coord,getNbGauss(),c,4);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for QUAD8 !");
      }
    case NORM_QUAD9:
      {
        std::vector<double> a(QUAD9A_REF,QUAD9A_REF+18);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(QUAD4B_REF,QUAD4B_REF+8);
            return GaussInfo(NORM_QUAD4,_my_gauss_coord,getNbGauss(),c,4);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for QUAD9 !");
      }
    case NORM_TETRA10:
      {
        std::vector<double> a(TETRA10A_REF,TETRA10A_REF+30),b(TETRA10B_REF,TETRA10B_REF+30);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(TETRA4A_REF,TETRA4A_REF+12);
            return GaussInfo(NORM_TETRA4,_my_gauss_coord,getNbGauss(),c,4);
          }
        if(IsSatisfy(b,_my_reference_coord))
          {
            std::vector<double> c(TETRA4B_REF,TETRA4B_REF+12);
            return GaussInfo(NORM_TETRA4,_my_gauss_coord,getNbGauss(),c,4);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for TETRA10 !");
      }
    case NORM_PYRA13:
      {
        std::vector<double> a(PYRA13A_REF,PYRA13A_REF+39),b(PYRA13B_REF,PYRA13B_REF+39);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(PYRA5A_REF,PYRA5A_REF+15);
            return GaussInfo(NORM_PYRA5,_my_gauss_coord,getNbGauss(),c,5);
          }
        if(IsSatisfy(b,_my_reference_coord))
          {
            std::vector<double> c(PYRA5B_REF,PYRA5B_REF+15);
            return GaussInfo(NORM_PYRA5,_my_gauss_coord,getNbGauss(),c,5);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for PYRA13 !");
      }
    case NORM_PENTA15:
      {
        std::vector<double> a(PENTA15A_REF,PENTA15A_REF+45),b(PENTA15B_REF,PENTA15B_REF+45);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(PENTA6A_REF,PENTA6A_REF+18);
            return GaussInfo(NORM_PENTA6,_my_gauss_coord,getNbGauss(),c,6);
          }
        if(IsSatisfy(b,_my_reference_coord))
          {
            std::vector<double> c(PENTA6B_REF,PENTA6B_REF+18);
            return GaussInfo(NORM_PENTA6,_my_gauss_coord,getNbGauss(),c,6);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for PENTA15 !");
      }
    case NORM_HEXA20:
      {
        std::vector<double> a(HEXA20A_REF,HEXA20A_REF+60),b(HEXA20B_REF,HEXA20B_REF+60);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(HEXA8A_REF,HEXA8A_REF+24);
            return GaussInfo(NORM_HEXA8,_my_gauss_coord,getNbGauss(),c,8);
          }
        if(IsSatisfy(b,_my_reference_coord))
          {
            std::vector<double> c(HEXA8B_REF,HEXA8B_REF+24);
            return GaussInfo(NORM_HEXA8,_my_gauss_coord,getNbGauss(),c,8);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for HEXA20 !");
      }
    case NORM_HEXA27:
      {
        std::vector<double> a(HEXA27A_REF,HEXA27A_REF+81);
        if(IsSatisfy(a,_my_reference_coord))
          {
            std::vector<double> c(HEXA8B_REF,HEXA8B_REF+24);
            return GaussInfo(NORM_HEXA8,_my_gauss_coord,getNbGauss(),c,8);
          }
        throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not recognized pattern for HEXA27 !");
      }
    default:
      throw INTERP_KERNEL::Exception("GaussInfo::convertToLinear : not implemented yet for other types than TETRA10, HEXA20, HEXA27, TRI3, QUAD8, QUAD8, PYRA13, PENTA15 !");
    }
}


bool GaussInfo::IsSatisfy(const std::vector<double>& ref1, const std::vector<double>& ref2)
{
  std::size_t sz(ref1.size());
  if(sz!=ref2.size())
    return false;
  for(std::size_t i=0;i<sz;i++)
    if(!IsEqual(ref1[i],ref2[i]))
      return false;
  return true;
}

/*!
 * Check coordinates
 */
bool GaussInfo::isSatisfy() 
{

  bool anIsSatisfy = ((_my_local_nb_ref == _my_nb_ref) && (_my_local_ref_dim == getReferenceCoordDim()));
  //Check coordinates
  if(anIsSatisfy)
    {
      for( int refId = 0; refId < _my_local_nb_ref; refId++ ) 
        {
          double* refCoord = &_my_reference_coord[ refId*_my_local_ref_dim ];
          double* localRefCoord = &_my_local_reference_coord[ refId*_my_local_ref_dim ];
          bool anIsEqual = false;
          for( int dimId = 0; dimId < _my_local_ref_dim; dimId++ ) 
            {
              anIsEqual = IsEqual( localRefCoord[dimId], refCoord[dimId]);
              if(!anIsEqual ) 
                {
                  return false;
                }
            }
        }
    }
  return anIsSatisfy;
}

std::vector<double> GaussInfo::NormalizeCoordinatesIfNecessary(NormalizedCellType ct, int inputDim, const std::vector<double>& inputArray)
{
  std::size_t sz(inputArray.size()),dim((std::size_t)inputDim);
  if(dim==0)
    throw INTERP_KERNEL::Exception("GaussInfo::NormalizeCoordinatesIfNecessary : invalid dimension ! Must be !=0 !");
  if(sz%dim!=0)
    throw INTERP_KERNEL::Exception("GaussInfo::NormalizeCoordinatesIfNecessary : invalid input array ! Inconsistent with the given dimension !");
  const CellModel& cm(CellModel::GetCellModel(ct));
  std::size_t baseDim((std::size_t)cm.getDimension());
  if(baseDim==dim)
    return inputArray;
  std::size_t nbOfItems(sz/dim);
  std::vector<double> ret(nbOfItems*baseDim);
  if(baseDim>dim)
    {
      for(std::size_t i=0;i<nbOfItems;i++)
        {
          std::size_t j=0;
          for(;j<dim;j++)
            ret[i*baseDim+j]=inputArray[i*dim+j];
          for(;j<baseDim;j++)
            ret[i*baseDim+j]=0.;
        }
    }
  else
    {
      for(std::size_t i=0;i<nbOfItems;i++)
        {
          std::size_t j=0;
          for(;j<baseDim;j++)
            ret[i*baseDim+j]=inputArray[i*dim+j];
        }
    }
  return ret;
}

typedef void (*MapToShapeFunction)(GaussInfo& obj);

/*!
 * Initialize the internal vectors
 */
void GaussInfo::initLocalInfo()
{
  bool aSatify = false;
  const CellModel& cellModel(CellModel::GetCellModel(_my_geometry));
  switch( _my_geometry ) 
    {
    case NORM_POINT1:
      _my_local_ref_dim = 0;
      _my_local_nb_ref  = 1;
      point1Init();
      aSatify = isSatisfy();
      CHECK_MACRO;
      break;

    case NORM_SEG2:
      _my_local_ref_dim = 1;
      _my_local_nb_ref  = 2;
      seg2Init();
      aSatify = isSatisfy();
      CHECK_MACRO;
      break;

    case NORM_SEG3:
      _my_local_ref_dim = 1;
      _my_local_nb_ref  = 3;
      seg3Init();
      aSatify = isSatisfy();
      CHECK_MACRO;
      break;

    case NORM_TRI3:
      _my_local_ref_dim = 2;
      _my_local_nb_ref  = 3;
      tria3aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          tria3bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_TRI6:
      _my_local_ref_dim = 2;
      _my_local_nb_ref  = 6;
      tria6aInit();
      aSatify = isSatisfy();
      if(!aSatify)
        {
          tria6bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_TRI7:
      _my_local_ref_dim = 2;
      _my_local_nb_ref  = 7;
      tria7aInit();
      aSatify = isSatisfy();
      CHECK_MACRO;
      break;

    case NORM_QUAD4:
      {
        _my_local_ref_dim = 2;
        _my_local_nb_ref  = 4;
        MapToShapeFunction QUAD4PTR[]={Quad4aInit,Quad4bInit,Quad4cInit,Quad4DegSeg2Init};
        std::size_t NB_OF_QUAD4PTR(sizeof(QUAD4PTR)/sizeof(MapToShapeFunction));
        for(std::size_t i=0;i<NB_OF_QUAD4PTR && !aSatify;i++)
          {
            (QUAD4PTR[i])(*this);
            aSatify = isSatisfy();
          }
        CHECK_MACRO;
        break;
      }
      break;

    case NORM_QUAD8:
      _my_local_ref_dim = 2;
      _my_local_nb_ref  = 8;
      quad8aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          quad8bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_QUAD9:
      _my_local_ref_dim = 2;
      _my_local_nb_ref  = 9;
      quad9aInit();
      aSatify = isSatisfy();
      CHECK_MACRO;
      break;

    case NORM_TETRA4:
      _my_local_ref_dim = 3;
      _my_local_nb_ref  = 4;
      tetra4aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          tetra4bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_TETRA10:
      _my_local_ref_dim = 3;
      _my_local_nb_ref  = 10;
      tetra10aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          tetra10bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_PYRA5:
      _my_local_ref_dim = 3;
      _my_local_nb_ref  = 5;
      pyra5aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          pyra5bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_PYRA13:
      _my_local_ref_dim = 3;
      _my_local_nb_ref  = 13;
      pyra13aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          pyra13bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_PENTA6:
      {
        _my_local_ref_dim = 3;
        _my_local_nb_ref  = 6;
        MapToShapeFunction PENTA6PTR[]={Penta6aInit,Penta6bInit,Penta6DegTria3aInit,Penta6DegTria3bInit};
        std::size_t NB_OF_PENTA6PTR(sizeof(PENTA6PTR)/sizeof(MapToShapeFunction));
        for(std::size_t i=0;i<NB_OF_PENTA6PTR && !aSatify;i++)
          {
            (PENTA6PTR[i])(*this);
            aSatify = isSatisfy();
          }
        CHECK_MACRO;
        break;
      }


      _my_local_ref_dim = 3;
      _my_local_nb_ref  = 6;
      penta6aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          penta6bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_PENTA15:
      {
        _my_local_ref_dim = 3;
        _my_local_nb_ref  = 15;
        MapToShapeFunction PENTA15PTR[]={Penta15aInit,Penta15bInit};
        std::size_t NB_OF_PENTA15PTR(sizeof(PENTA15PTR)/sizeof(MapToShapeFunction));
        for(std::size_t i=0;i<NB_OF_PENTA15PTR && !aSatify;i++)
          {
            (PENTA15PTR[i])(*this);
            aSatify = isSatisfy();
          }
        CHECK_MACRO;
        break;
      }

    case NORM_HEXA8:
      {
        _my_local_ref_dim = 3;
        _my_local_nb_ref  = 8;
        MapToShapeFunction HEXA8PTR[]={Hexa8aInit,Hexa8bInit,Hexa8DegQuad4aInit,Hexa8DegQuad4bInit,Hexa8DegQuad4cInit};
        std::size_t NB_OF_HEXA8PTR(sizeof(HEXA8PTR)/sizeof(MapToShapeFunction));
        for(std::size_t i=0;i<NB_OF_HEXA8PTR && !aSatify;i++)
          {
            (HEXA8PTR[i])(*this);
            aSatify = isSatisfy();
          }
        CHECK_MACRO;
        break;
      }

    case NORM_HEXA20:
      _my_local_ref_dim = 3;
      _my_local_nb_ref  = 20;
      hexa20aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          hexa20bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_HEXA27:
      _my_local_ref_dim = 3;
      _my_local_nb_ref  = 27;
      hexa27aInit();
      aSatify = isSatisfy();
      CHECK_MACRO
      break;

    default:
      throw INTERP_KERNEL::Exception("Not managed cell type !");
      break;
    }
}

/**
 * Return shape function value by node id
 */
const double* GaussInfo::getFunctionValues( const int theGaussId ) const 
{
  return &_my_function_value[ _my_nb_ref*theGaussId ];
}

void GaussInfo::point1Init()
{
  double *funValue(&_my_function_value[0]);
  funValue[0] = 1. ;
}

/*!
 * Init Segment 2 Reference coordinates ans Shape function.
 */
void GaussInfo::seg2Init() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = SEG2_REF[0];
    break;
  case 1:
    coords[0] = SEG2_REF[1];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.5*(1.0 - gc[0]);
  funValue[1] = 0.5*(1.0 + gc[0]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Segment 3 Reference coordinates ans Shape function.
 */
void GaussInfo::seg3Init() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = SEG3_REF[0];
    break;
  case 1:
    coords[0] = SEG3_REF[1];
    break;
  case 2:
    coords[0] = SEG3_REF[2];
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = -0.5*(1.0 - gc[0])*gc[0];
  funValue[1] = 0.5*(1.0 + gc[0])*gc[0];
  funValue[2] = (1.0 + gc[0])*(1.0 - gc[0]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Triangle Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::tria3aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = TRIA3A_REF[0];
    coords[1] = TRIA3A_REF[1];
    break;
  case 1:
    coords[0] = TRIA3A_REF[2];
    coords[1] = TRIA3A_REF[3];
    break;
  case 2:
    coords[0] = TRIA3A_REF[4];
    coords[1] = TRIA3A_REF[5];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.5*(1.0 + gc[1]);
  funValue[1] = -0.5*(gc[0] + gc[1]);
  funValue[2] = 0.5*(1.0 + gc[0]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Triangle Reference coordinates ans Shape function.
 * Case B.
 */
void GaussInfo::tria3bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = TRIA3B_REF[0];
    coords[1] = TRIA3B_REF[1];
    break;
  case 1:
    coords[0] = TRIA3B_REF[2];
    coords[1] = TRIA3B_REF[3];
    break;
  case 2:
    coords[0] = TRIA3B_REF[4];
    coords[1] = TRIA3B_REF[5];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 1.0 - gc[0] - gc[1];
  funValue[1] = gc[0];
  funValue[2] = gc[1];
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadratic Triangle Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::tria6aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = TRIA6A_REF[0];
    coords[1] = TRIA6A_REF[1];
    break;
  case 1:
    coords[0] = TRIA6A_REF[2];
    coords[1] = TRIA6A_REF[3];
    break;
  case 2:
    coords[0] = TRIA6A_REF[4];
    coords[1] = TRIA6A_REF[5];
    break;
  case 3:
    coords[0] = TRIA6A_REF[6];
    coords[1] = TRIA6A_REF[7];
    break;
  case 4:
    coords[0] = TRIA6A_REF[8];
    coords[1] = TRIA6A_REF[9];
    break;
  case 5:
    coords[0] = TRIA6A_REF[10];
    coords[1] = TRIA6A_REF[11];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.5*(1.0 + gc[1])*gc[1];
  funValue[1] = 0.5*(gc[0] + gc[1])*(gc[0] + gc[1] + 1);
  funValue[2] = 0.5*(1.0 + gc[0])*gc[0];
  funValue[3] = -1.0*(1.0 + gc[1])*(gc[0] + gc[1]);
  funValue[4] = -1.0*(1.0 + gc[0])*(gc[0] + gc[1]);
  funValue[5] = (1.0 + gc[1])*(1.0 + gc[1]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadratic Triangle Reference coordinates ans Shape function.
 * Case B.
 */
void GaussInfo::tria6bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = TRIA6B_REF[0];
    coords[1] = TRIA6B_REF[1];
    break;
  case 1:
    coords[0] = TRIA6B_REF[2];
    coords[1] = TRIA6B_REF[3];
    break;
  case 2:
    coords[0] = TRIA6B_REF[4];
    coords[1] = TRIA6B_REF[5];
    break;
  case 3:
    coords[0] = TRIA6B_REF[6];
    coords[1] = TRIA6B_REF[7];
    break;
  case 4:
    coords[0] = TRIA6B_REF[8];
    coords[1] = TRIA6B_REF[9];
    break;
  case 5:
    coords[0] = TRIA6B_REF[10];
    coords[1] = TRIA6B_REF[11];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = (1.0 - gc[0] - gc[1])*(1.0 - 2.0*gc[0] - 2.0*gc[1]);
  funValue[1] = gc[0]*(2.0*gc[0] - 1.0);
  funValue[2] = gc[1]*(2.0*gc[1] - 1.0);
  funValue[3] = 4.0*gc[0]*(1.0 - gc[0] - gc[1]);
  funValue[4] = 4.0*gc[0]*gc[1];
  funValue[5] = 4.0*gc[1]*(1.0 - gc[0] - gc[1]);
  SHAPE_FUN_MACRO_END;
}

void GaussInfo::tria7aInit()
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = TRIA7A_REF[0];
    coords[1] = TRIA7A_REF[1];
    break;
  case 1:
    coords[0] = TRIA7A_REF[2];
    coords[1] = TRIA7A_REF[3];
    break;
  case 2:
    coords[0] = TRIA7A_REF[4];
    coords[1] = TRIA7A_REF[5];
    break;
  case 3:
    coords[0] = TRIA7A_REF[6];
    coords[1] = TRIA7A_REF[7];
    break;
  case 4:
    coords[0] = TRIA7A_REF[8];
    coords[1] = TRIA7A_REF[9];
    break;
  case 5:
    coords[0] = TRIA7A_REF[10];
    coords[1] = TRIA7A_REF[11];
    break;
  case 6:
    coords[0] = TRIA7A_REF[12];
    coords[1] = TRIA7A_REF[13];
    break;
  LOCAL_COORD_MACRO_END;

  SHAPE_FUN_MACRO_BEGIN;
  funValue[0]=1-3*(gc[0]+gc[1])+2*(gc[0]*gc[0]+gc[1]*gc[1])+7*gc[0]*gc[1]-3*gc[0]*gc[1]*(gc[0]+gc[1]);
  funValue[1]=gc[0]*(-1+2*gc[0]+3*gc[1]-3*gc[1]*(gc[0]+gc[1]));
  funValue[2]=gc[1]*(-1.+3.*gc[0]+2.*gc[1]-3.*gc[0]*(gc[0]+gc[1]));
  funValue[3]=4*gc[0]*(1-gc[0]-4*gc[1]+3*gc[1]*(gc[0]+gc[1]));
  funValue[4]=4*gc[0]*gc[1]*(-2+3*(gc[0]+gc[1]));
  funValue[5]=4*gc[1]*(1-4*gc[0]-gc[1]+3*gc[0]*(gc[0]+gc[1]));
  funValue[6]=27*gc[0]*gc[1]*(1-gc[0]-gc[1]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadrangle Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::quad4aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = QUAD4A_REF[0];
    coords[1] = QUAD4A_REF[1];
    break;
  case 1:
    coords[0] = QUAD4A_REF[2];
    coords[1] = QUAD4A_REF[3];
    break;
  case 2:
    coords[0] = QUAD4A_REF[4];
    coords[1] = QUAD4A_REF[5];
    break;
  case 3:
    coords[0] = QUAD4A_REF[6];
    coords[1] = QUAD4A_REF[7];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.25*(1.0 + gc[1])*(1.0 - gc[0]);
  funValue[1] = 0.25*(1.0 - gc[1])*(1.0 - gc[0]);
  funValue[2] = 0.25*(1.0 - gc[1])*(1.0 + gc[0]);
  funValue[3] = 0.25*(1.0 + gc[0])*(1.0 + gc[1]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadrangle Reference coordinates ans Shape function.
 * Case B.
 */
void GaussInfo::quad4bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = QUAD4B_REF[0];
    coords[1] = QUAD4B_REF[1];
    break;
  case 1:
    coords[0] = QUAD4B_REF[2];
    coords[1] = QUAD4B_REF[3];
    break;
  case 2:
    coords[0] = QUAD4B_REF[4];
    coords[1] = QUAD4B_REF[5];
    break;
  case 3:
    coords[0] = QUAD4B_REF[6];
    coords[1] = QUAD4B_REF[7];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.25*(1.0 - gc[0])*(1.0 - gc[1]);
  funValue[1] = 0.25*(1.0 + gc[0])*(1.0 - gc[1]);
  funValue[2] = 0.25*(1.0 + gc[0])*(1.0 + gc[1]);
  funValue[3] = 0.25*(1.0 - gc[0])*(1.0 + gc[1]);
  SHAPE_FUN_MACRO_END;
}

void GaussInfo::quad4cInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] =  1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] = -1.0;
   break;

   LOCAL_COORD_MACRO_END;

   SHAPE_FUN_MACRO_BEGIN;
   funValue[0] = 0.25*(1.0 - gc[0])*(1.0 - gc[1]);
   funValue[1] = 0.25*(1.0 - gc[0])*(1.0 + gc[1]);
   funValue[2] = 0.25*(1.0 + gc[0])*(1.0 + gc[1]);
   funValue[3] = 0.25*(1.0 + gc[0])*(1.0 - gc[1]);
   SHAPE_FUN_MACRO_END;
}

/*!
 * This shapefunc map is same as degenerated seg2Init
 */
void GaussInfo::quad4DegSeg2Init()
{
  LOCAL_COORD_MACRO_BEGIN;
 case  0:
   coords[0] = -1.0;
   coords[1] =  0.0;
   break;
 case  1:
   coords[0] =  1.0;
   coords[1] =  0.0;
   break;
 case  2:
   coords[0] =  0.0;
   coords[1] =  0.0;
   break;
 case  3:
   coords[0] =  0.0;
   coords[1] =  0.0;
   break;
   LOCAL_COORD_MACRO_END;

   SHAPE_FUN_MACRO_BEGIN;
   funValue[0] = 0.5*(1.0 - gc[0]);
   funValue[1] = 0.5*(1.0 + gc[0]);
   funValue[2] = 0.;
   funValue[3] = 0.;
   SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadratic Quadrangle Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::quad8aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = QUAD8A_REF[0];
    coords[1] = QUAD8A_REF[1];
    break;
  case 1:
    coords[0] = QUAD8A_REF[2];
    coords[1] = QUAD8A_REF[3];
    break;
  case 2:
    coords[0] = QUAD8A_REF[4];
    coords[1] = QUAD8A_REF[5];
    break;
  case 3:
    coords[0] = QUAD8A_REF[6];
    coords[1] = QUAD8A_REF[7];
    break;
  case 4:
    coords[0] = QUAD8A_REF[8];
    coords[1] = QUAD8A_REF[9];
    break;
  case 5:
    coords[0] = QUAD8A_REF[10];
    coords[1] = QUAD8A_REF[11];
    break;
  case 6:
    coords[0] = QUAD8A_REF[12];
    coords[1] = QUAD8A_REF[13];
    break;
  case 7:
    coords[0] = QUAD8A_REF[14];
    coords[1] = QUAD8A_REF[15];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.25*(1.0 + gc[1])*(1.0 - gc[0])*(gc[1] - gc[0] - 1.0);
  funValue[1] = 0.25*(1.0 - gc[1])*(1.0 - gc[0])*(-gc[1] - gc[0] - 1.0);
  funValue[2] = 0.25*(1.0 - gc[1])*(1.0 + gc[0])*(-gc[1] + gc[0] - 1.0);
  funValue[3] = 0.25*(1.0 + gc[1])*(1.0 + gc[0])*(gc[1] + gc[0] - 1.0);
  funValue[4] = 0.5*(1.0 - gc[0])*(1.0 - gc[1])*(1.0 + gc[1]);
  funValue[5] = 0.5*(1.0 - gc[1])*(1.0 - gc[0])*(1.0 + gc[0]);
  funValue[6] = 0.5*(1.0 + gc[0])*(1.0 - gc[1])*(1.0 + gc[1]);
  funValue[7] = 0.5*(1.0 + gc[1])*(1.0 - gc[0])*(1.0 + gc[0]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadratic Quadrangle Reference coordinates ans Shape function.
 * Case B.
 */
void GaussInfo::quad8bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = QUAD8B_REF[0];
    coords[1] = QUAD8B_REF[1];
    break;
  case 1:
    coords[0] = QUAD8B_REF[2];
    coords[1] = QUAD8B_REF[3];
    break;
  case 2:
    coords[0] = QUAD8B_REF[4];
    coords[1] = QUAD8B_REF[5];
    break;
  case 3:
    coords[0] = QUAD8B_REF[6];
    coords[1] = QUAD8B_REF[7];
    break;
  case 4:
    coords[0] = QUAD8B_REF[8];
    coords[1] = QUAD8B_REF[9];
    break;
  case 5:
    coords[0] = QUAD8B_REF[10];
    coords[1] = QUAD8B_REF[11];
    break;
  case 6:
    coords[0] = QUAD8B_REF[12];
    coords[1] = QUAD8B_REF[13];
    break;
  case 7:
    coords[0] = QUAD8B_REF[14];
    coords[1] = QUAD8B_REF[15];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.25*(1.0 - gc[0])*(1.0 - gc[1])*(-1.0 - gc[0] - gc[1]);
  funValue[1] = 0.25*(1.0 + gc[0])*(1.0 - gc[1])*(-1.0 + gc[0] - gc[1]);
  funValue[2] = 0.25*(1.0 + gc[0])*(1.0 + gc[1])*(-1.0 + gc[0] + gc[1]);
  funValue[3] = 0.25*(1.0 - gc[0])*(1.0 + gc[1])*(-1.0 - gc[0] + gc[1]);
  funValue[4] = 0.5*(1.0 - gc[0]*gc[0])*(1.0 - gc[1]);
  funValue[5] = 0.5*(1.0 - gc[1]*gc[1])*(1.0 + gc[0]);
  funValue[6] = 0.5*(1.0 - gc[0]*gc[0])*(1.0 + gc[1]);
  funValue[7] = 0.5*(1.0 - gc[1]*gc[1])*(1.0 - gc[0]);
  SHAPE_FUN_MACRO_END;
}

void GaussInfo::quad9aInit()
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = QUAD9A_REF[0];
    coords[1] = QUAD9A_REF[1];
    break;
  case 1:
    coords[0] = QUAD9A_REF[2];
    coords[1] = QUAD9A_REF[3];
    break;
  case 2:
    coords[0] = QUAD9A_REF[4];
    coords[1] = QUAD9A_REF[5];
    break;
  case 3:
    coords[0] = QUAD9A_REF[6];
    coords[1] = QUAD9A_REF[7];
    break;
  case 4:
    coords[0] = QUAD9A_REF[8];
    coords[1] = QUAD9A_REF[9];
    break;
  case 5:
    coords[0] = QUAD9A_REF[10];
    coords[1] = QUAD9A_REF[11];
    break;
  case 6:
    coords[0] = QUAD9A_REF[12];
    coords[1] = QUAD9A_REF[13];
    break;
  case 7:
    coords[0] = QUAD9A_REF[14];
    coords[1] = QUAD9A_REF[15];
    break;
  case 8:
    coords[0] = QUAD9A_REF[16];
    coords[1] = QUAD9A_REF[17];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.25*gc[0]*gc[1]*(gc[0]-1.)*(gc[1]-1.);
  funValue[1] = 0.25*gc[0]*gc[1]*(gc[0]+1.)*(gc[1]-1.);
  funValue[2] = 0.25*gc[0]*gc[1]*(gc[0]+1.)*(gc[1]+1.);
  funValue[3] = 0.25*gc[0]*gc[1]*(gc[0]-1.)*(gc[1]+1.);
  funValue[4] = 0.5*(1.-gc[0]*gc[0])*gc[1]*(gc[1]-1.);
  funValue[5] = 0.5*gc[0]*(gc[0]+1.)*(1.-gc[1]*gc[1]);
  funValue[6] = 0.5*(1.-gc[0]*gc[0])*gc[1]*(gc[1]+1.);
  funValue[7] = 0.5*gc[0]*(gc[0]-1.)*(1.-gc[1]*gc[1]);
  funValue[8] = (1.-gc[0]*gc[0])*(1.-gc[1]*gc[1]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Tetrahedron Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::tetra4aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = TETRA4A_REF[0];
    coords[1] = TETRA4A_REF[1];
    coords[2] = TETRA4A_REF[2];
    break;
  case 1:
    coords[0] = TETRA4A_REF[3];
    coords[1] = TETRA4A_REF[4];
    coords[2] = TETRA4A_REF[5];
    break;
  case 2:
    coords[0] = TETRA4A_REF[6];
    coords[1] = TETRA4A_REF[7];
    coords[2] = TETRA4A_REF[8];
    break;
  case 3:
    coords[0] = TETRA4A_REF[9];
    coords[1] = TETRA4A_REF[10];
    coords[2] = TETRA4A_REF[11];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = gc[1];
  funValue[1] = gc[2];
  funValue[2] = 1.0 - gc[0] - gc[1] - gc[2];
  funValue[3] = gc[0];
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Tetrahedron Reference coordinates ans Shape function.
 * Case B.
 */
void GaussInfo::tetra4bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = TETRA4B_REF[0];
    coords[1] = TETRA4B_REF[1];
    coords[2] = TETRA4B_REF[2];
    break;
  case 1:
    coords[0] = TETRA4B_REF[3];
    coords[1] = TETRA4B_REF[4];
    coords[2] = TETRA4B_REF[5];
    break;
  case 2:
    coords[0] = TETRA4B_REF[6];
    coords[1] = TETRA4B_REF[7];
    coords[2] = TETRA4B_REF[8];
    break;
  case 3:
    coords[0] = TETRA4B_REF[9];
    coords[1] = TETRA4B_REF[10];
    coords[2] = TETRA4B_REF[11];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = gc[1];
  funValue[2] = gc[2];
  funValue[1] = 1.0 - gc[0] - gc[1] - gc[2];
  funValue[3] = gc[0];
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadratic Tetrahedron Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::tetra10aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = TETRA10A_REF[0];
    coords[1] = TETRA10A_REF[1];
    coords[2] = TETRA10A_REF[2];
    break;
  case 1:
    coords[0] = TETRA10A_REF[3];
    coords[1] = TETRA10A_REF[4];
    coords[2] = TETRA10A_REF[5];
    break;
  case 2:
    coords[0] = TETRA10A_REF[6];
    coords[1] = TETRA10A_REF[7];
    coords[2] = TETRA10A_REF[8];
    break;
  case 3:
    coords[0] = TETRA10A_REF[9];
    coords[1] = TETRA10A_REF[10];
    coords[2] = TETRA10A_REF[11];
    break;
  case 4:
    coords[0] = TETRA10A_REF[12];
    coords[1] = TETRA10A_REF[13];
    coords[2] = TETRA10A_REF[14];
    break;
  case 5:
    coords[0] = TETRA10A_REF[15];
    coords[1] = TETRA10A_REF[16];
    coords[2] = TETRA10A_REF[17];
    break;
  case 6:
    coords[0] = TETRA10A_REF[18];
    coords[1] = TETRA10A_REF[19];
    coords[2] = TETRA10A_REF[20];
    break;
  case 7:
    coords[0] = TETRA10A_REF[21];
    coords[1] = TETRA10A_REF[22];
    coords[2] = TETRA10A_REF[23];
    break;
  case 8:
    coords[0] = TETRA10A_REF[24];
    coords[1] = TETRA10A_REF[25];
    coords[2] = TETRA10A_REF[26];
    break;
  case 9:
    coords[0] = TETRA10A_REF[27];
    coords[1] = TETRA10A_REF[28];
    coords[2] = TETRA10A_REF[29];
    break;
  LOCAL_COORD_MACRO_END;

  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = gc[1]*(2.0*gc[1] - 1.0);
  funValue[1] = gc[2]*(2.0*gc[2] - 1.0);
  funValue[2] = (1.0 - gc[0] - gc[1] - gc[2])*(1.0 - 2.0*gc[0] - 2.0*gc[1] - 2.0*gc[2]);
  funValue[3] = gc[0]*(2.0*gc[0] - 1.0);
  funValue[4] = 4.0*gc[1]*gc[2];
  funValue[5] = 4.0*gc[2]*(1.0 - gc[0] - gc[1] - gc[2]);
  funValue[6] = 4.0*gc[1]*(1.0 - gc[0] - gc[1] - gc[2]);
  funValue[7] = 4.0*gc[0]*gc[1];
  funValue[8] = 4.0*gc[0]*gc[2];
  funValue[9] = 4.0*gc[0]*(1.0 - gc[0] - gc[1] - gc[2]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadratic Tetrahedron Reference coordinates ans Shape function.
 * Case B.
 */
void GaussInfo::tetra10bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = TETRA10B_REF[0];
    coords[1] = TETRA10B_REF[1];
    coords[2] = TETRA10B_REF[2];
    break;
  case 1:
    coords[0] = TETRA10B_REF[3];
    coords[1] = TETRA10B_REF[4];
    coords[2] = TETRA10B_REF[5];
    break;
  case 2:
    coords[0] = TETRA10B_REF[6];
    coords[1] = TETRA10B_REF[7];
    coords[2] = TETRA10B_REF[8];
    break;
  case 3:
    coords[0] = TETRA10B_REF[9];
    coords[1] = TETRA10B_REF[10];
    coords[2] = TETRA10B_REF[11];
    break;
  case 4:
    coords[0] = TETRA10B_REF[12];
    coords[1] = TETRA10B_REF[13];
    coords[2] = TETRA10B_REF[14];
    break;
  case 5:
    coords[0] = TETRA10B_REF[15];
    coords[1] = TETRA10B_REF[16];
    coords[2] = TETRA10B_REF[17];
    break;
  case 6:
    coords[0] = TETRA10B_REF[18];
    coords[1] = TETRA10B_REF[19];
    coords[2] = TETRA10B_REF[20];
    break;
  case 7:
    coords[0] = TETRA10B_REF[21];
    coords[1] = TETRA10B_REF[22];
    coords[2] = TETRA10B_REF[23];
    break;
  case 8:
    coords[0] = TETRA10B_REF[24];
    coords[1] = TETRA10B_REF[25];
    coords[2] = TETRA10B_REF[26];
    break;
  case 9:
    coords[0] = TETRA10B_REF[27];
    coords[1] = TETRA10B_REF[28];
    coords[2] = TETRA10B_REF[29];
    break;
  LOCAL_COORD_MACRO_END;
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = gc[1]*(2.0*gc[1] - 1.0);
  funValue[2] = gc[2]*(2.0*gc[2] - 1.0);
  funValue[1] = (1.0 - gc[0] - gc[1] - gc[2])*(1.0 - 2.0*gc[0] - 2.0*gc[1] - 2.0*gc[2]);
  funValue[3] = gc[0]*(2.0*gc[0] - 1.0);
  funValue[6] = 4.0*gc[1]*gc[2];
  funValue[5] = 4.0*gc[2]*(1.0 - gc[0] - gc[1] - gc[2]);
  funValue[4] = 4.0*gc[1]*(1.0 - gc[0] - gc[1] - gc[2]);
  funValue[7] = 4.0*gc[0]*gc[1];
  funValue[9] = 4.0*gc[0]*gc[2];
  funValue[8] = 4.0*gc[0]*(1.0 - gc[0] - gc[1] - gc[2]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Pyramid Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::pyra5aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = PYRA5A_REF[0];
    coords[1] = PYRA5A_REF[1];
    coords[2] = PYRA5A_REF[2];
    break;
  case 1:
    coords[0] = PYRA5A_REF[3];
    coords[1] = PYRA5A_REF[4];
    coords[2] = PYRA5A_REF[5];
    break;
  case 2:
    coords[0] = PYRA5A_REF[6];
    coords[1] = PYRA5A_REF[7];
    coords[2] = PYRA5A_REF[8];
    break;
  case 3:
    coords[0] = PYRA5A_REF[9];
    coords[1] = PYRA5A_REF[10];
    coords[2] = PYRA5A_REF[11];
    break;
  case 4:
    coords[0] = PYRA5A_REF[12];
    coords[1] = PYRA5A_REF[13];
    coords[2] = PYRA5A_REF[14];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.25*(-gc[0] + gc[1] - 1.0)*(-gc[0] - gc[1] - 1.0)*(1.0 - gc[2]);
  funValue[1] = 0.25*(-gc[0] - gc[1] - 1.0)*(+gc[0] - gc[1] - 1.0)*(1.0 - gc[2]);
  funValue[2] = 0.25*(+gc[0] + gc[1] - 1.0)*(+gc[0] - gc[1] - 1.0)*(1.0 - gc[2]);
  funValue[3] = 0.25*(+gc[0] + gc[1] - 1.0)*(-gc[0] + gc[1] - 1.0)*(1.0 - gc[2]);
  funValue[4] = gc[2];
  SHAPE_FUN_MACRO_END;
}
/*!
 * Init Pyramid Reference coordinates ans Shape function.
 * Case B.
 */
void GaussInfo::pyra5bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = PYRA5B_REF[0];
    coords[1] = PYRA5B_REF[1];
    coords[2] = PYRA5B_REF[2];
    break;
  case 1:
    coords[0] = PYRA5B_REF[3];
    coords[1] = PYRA5B_REF[4];
    coords[2] = PYRA5B_REF[5];
    break;
  case 2:
    coords[0] = PYRA5B_REF[6];
    coords[1] = PYRA5B_REF[7];
    coords[2] = PYRA5B_REF[8];
    break;
  case 3:
    coords[0] = PYRA5B_REF[9];
    coords[1] = PYRA5B_REF[10];
    coords[2] = PYRA5B_REF[11];
    break;
  case 4:
    coords[0] = PYRA5B_REF[12];
    coords[1] = PYRA5B_REF[13];
    coords[2] = PYRA5B_REF[14];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.25*(-gc[0] + gc[1] - 1.0)*(-gc[0] - gc[1] - 1.0)*(1.0 - gc[2]);
  funValue[3] = 0.25*(-gc[0] - gc[1] - 1.0)*(+gc[0] - gc[1] - 1.0)*(1.0 - gc[2]);
  funValue[2] = 0.25*(+gc[0] + gc[1] - 1.0)*(+gc[0] - gc[1] - 1.0)*(1.0 - gc[2]);
  funValue[1] = 0.25*(+gc[0] + gc[1] - 1.0)*(-gc[0] + gc[1] - 1.0)*(1.0 - gc[2]);
  funValue[4] = gc[2];
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadratic Pyramid Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::pyra13aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = PYRA13A_REF[0];
    coords[1] = PYRA13A_REF[1];
    coords[2] = PYRA13A_REF[2];
    break;
  case 1:
    coords[0] = PYRA13A_REF[3];
    coords[1] = PYRA13A_REF[4];
    coords[2] = PYRA13A_REF[5];
    break;
  case 2:
    coords[0] = PYRA13A_REF[6];
    coords[1] = PYRA13A_REF[7];
    coords[2] = PYRA13A_REF[8];
    break;
  case 3:
    coords[0] = PYRA13A_REF[9];
    coords[1] = PYRA13A_REF[10];
    coords[2] = PYRA13A_REF[11];
    break;
  case 4:
    coords[0] = PYRA13A_REF[12];
    coords[1] = PYRA13A_REF[13];
    coords[2] = PYRA13A_REF[14];
    break;
  case 5:
    coords[0] = PYRA13A_REF[15];
    coords[1] = PYRA13A_REF[16];
    coords[2] = PYRA13A_REF[17];
    break;
  case 6:
    coords[0] = PYRA13A_REF[18];
    coords[1] = PYRA13A_REF[19];
    coords[2] = PYRA13A_REF[20];
    break;
  case 7:
    coords[0] = PYRA13A_REF[21];
    coords[1] = PYRA13A_REF[22];
    coords[2] = PYRA13A_REF[23];
    break;
  case 8:
    coords[0] = PYRA13A_REF[24];
    coords[1] = PYRA13A_REF[25];
    coords[2] = PYRA13A_REF[26];
    break;
  case 9:
    coords[0] = PYRA13A_REF[27];
    coords[1] = PYRA13A_REF[28];
    coords[2] = PYRA13A_REF[29];
    break;
  case 10:
    coords[0] = PYRA13A_REF[30];
    coords[1] = PYRA13A_REF[31];
    coords[2] = PYRA13A_REF[32];
    break;
  case 11:
    coords[0] = PYRA13A_REF[33];
    coords[1] = PYRA13A_REF[34];
    coords[2] = PYRA13A_REF[35];
    break;
  case 12:
    coords[0] = PYRA13A_REF[36];
    coords[1] = PYRA13A_REF[37];
    coords[2] = PYRA13A_REF[38];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.5*(-gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] - gc[1] + gc[2] - 1.0)*
    (gc[0] - 0.5)/(1.0 - gc[2]);
  funValue[1] = 0.5*(-gc[0] - gc[1] + gc[2] - 1.0)*(+gc[0] - gc[1] + gc[2] - 1.0)*
    (gc[1] - 0.5)/(1.0 - gc[2]);
  funValue[2] = 0.5*(+gc[0] - gc[1] + gc[2] - 1.0)*(+gc[0] + gc[1] + gc[2] - 1.0)*
    (-gc[0] - 0.5)/(1.0 - gc[2]);
  funValue[3] = 0.5*(+gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] + gc[1] + gc[2] - 1.0)*
    (-gc[1] - 0.5)/(1.0 - gc[2]);
  
  funValue[4] = 2.0*gc[2]*(gc[2] - 0.5);
  
  funValue[5] = 0.5*(-gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] - gc[1] + gc[2] - 1.0)*
    (gc[0] - gc[1] + gc[2] - 1.0)/(1.0 - gc[2]);
  funValue[6] = 0.5*(-gc[0] - gc[1] + gc[2] - 1.0)*(gc[0] - gc[1] + gc[2] - 1.0)*
    (gc[0] + gc[1] + gc[2] - 1.0)/(1.0 - gc[2]);
  funValue[7] = 0.5*(gc[0] - gc[1] + gc[2] - 1.0)*(gc[0] + gc[1] + gc[2] - 1.0)*
    (-gc[0] + gc[1] + gc[2] - 1.0)/(1.0 - gc[2]);
  funValue[8] = 0.5*(gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] + gc[1] + gc[2] - 1.0)*
    (-gc[0] - gc[1] + gc[2] - 1.0)/(1.0 - gc[2]);
  
  funValue[9] = 0.5*gc[2]*(-gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] - gc[1] + gc[2] - 1.0)/
    (1.0 - gc[2]);
  funValue[10] = 0.5*gc[2]*(-gc[0] - gc[1] + gc[2] - 1.0)*(gc[0] - gc[1] + gc[2] - 1.0)/
    (1.0 - gc[2]);
  funValue[11] = 0.5*gc[2]*(gc[0] - gc[1] + gc[2] - 1.0)*(gc[0] + gc[1] + gc[2] - 1.0)/
    (1.0 - gc[2]);
  funValue[12] = 0.5*gc[2]*(gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] + gc[1] + gc[2] - 1.0)/
    (1.0 - gc[2]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Quadratic Pyramid Reference coordinates ans Shape function.
 * Case B.
 */
void GaussInfo::pyra13bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = PYRA13B_REF[0];
    coords[1] = PYRA13B_REF[1];
    coords[2] = PYRA13B_REF[2];
    break;
  case 1:
    coords[0] = PYRA13B_REF[3];
    coords[1] = PYRA13B_REF[4];
    coords[2] = PYRA13B_REF[5];
    break;
  case 2:
    coords[0] = PYRA13B_REF[6];
    coords[1] = PYRA13B_REF[7];
    coords[2] = PYRA13B_REF[8];
    break;
  case 3:
    coords[0] = PYRA13B_REF[9];
    coords[1] = PYRA13B_REF[10];
    coords[2] = PYRA13B_REF[11];
    break;
  case 4:
    coords[0] = PYRA13B_REF[12];
    coords[1] = PYRA13B_REF[13];
    coords[2] = PYRA13B_REF[14];
    break;
  case 5:
    coords[0] = PYRA13B_REF[15];
    coords[1] = PYRA13B_REF[16];
    coords[2] = PYRA13B_REF[17];
    break;
  case 6:
    coords[0] = PYRA13B_REF[18];
    coords[1] = PYRA13B_REF[19];
    coords[2] = PYRA13B_REF[20];
    break;
  case 7:
    coords[0] = PYRA13B_REF[21];
    coords[1] = PYRA13B_REF[22];
    coords[2] = PYRA13B_REF[23];
    break;
  case 8:
    coords[0] = PYRA13B_REF[24];
    coords[1] = PYRA13B_REF[25];
    coords[2] = PYRA13B_REF[26];
    break;
  case 9:
    coords[0] = PYRA13B_REF[27];
    coords[1] = PYRA13B_REF[28];
    coords[2] = PYRA13B_REF[29];
    break;
  case 10:
    coords[0] = PYRA13B_REF[30];
    coords[1] = PYRA13B_REF[31];
    coords[2] = PYRA13B_REF[32];
    break;
  case 11:
    coords[0] = PYRA13B_REF[33];
    coords[1] = PYRA13B_REF[34];
    coords[2] = PYRA13B_REF[35];
    break;
  case 12:
    coords[0] = PYRA13B_REF[36];
    coords[1] = PYRA13B_REF[37];
    coords[2] = PYRA13B_REF[38];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] =0.5*(-gc[0]+gc[1]+gc[2]-1.0)*(-gc[0]-gc[1]+gc[2]-1.0)*(gc[0]-0.5)/(1.0-gc[2]);
  funValue[1] =0.5*(+gc[0]+gc[1]+gc[2]-1.0)*(-gc[0]+gc[1]+gc[2]-1.0)*(-gc[1]-0.5)/(1.0-gc[2]);
  funValue[2] =0.5*(+gc[0]-gc[1]+gc[2]-1.0)*(+gc[0]+gc[1]+gc[2]-1.0)*(-gc[0]-0.5)/(1.0-gc[2]);
  funValue[3] =0.5*(-gc[0]-gc[1]+gc[2]-1.0)*(+gc[0]-gc[1]+gc[2]-1.0)*(gc[1]-0.5)/(1.0-gc[2]);
  
  funValue[4] =2.*gc[2]*(gc[2]-0.5);
  
  funValue[5] =-0.5*(gc[0]+gc[1]+gc[2]-1.0)*(-gc[0]+gc[1]+gc[2]-1.0)*(-gc[0]-gc[1]+gc[2]-1.0)/(1.0-gc[2]);
  funValue[6] =-0.5*(gc[0]-gc[1]+gc[2]-1.0)*(gc[0]+gc[1]+gc[2]-1.0)*(-gc[0]+gc[1]+gc[2]-1.0)/(1.0-gc[2]);
  funValue[7] =-0.5*(-gc[0]-gc[1]+gc[2]-1.0)*(gc[0]-gc[1]+gc[2]-1.0)*(gc[0]+gc[1]+gc[2]-1.0)/(1.0-gc[2]);
  funValue[8] =-0.5*(-gc[0]+gc[1]+gc[2]-1.0)*(-gc[0]-gc[1]+gc[2]-1.0)*(gc[0]-gc[1]+gc[2]-1.0)/(1.0-gc[2]);
  
  funValue[9] =gc[2]*(-gc[0]+gc[1]+gc[2]-1.0)*(-gc[0]-gc[1]+gc[2]-1.0)/(1.0-gc[2]);
  funValue[10]=gc[2]*(gc[0]+gc[1]+gc[2]-1.0)*(-gc[0]+gc[1]+gc[2]-1.0)/(1.0-gc[2]);
  funValue[11]=gc[2]*(gc[0]-gc[1]+gc[2]-1.0)*(gc[0]+gc[1]+gc[2]-1.0)/(1.0-gc[2]);
  funValue[12]=gc[2]*(-gc[0]-gc[1]+gc[2]-1.0)*(gc[0]-gc[1]+gc[2]-1.0)/(1.0-gc[2]);
  
  SHAPE_FUN_MACRO_END;
}


/*!
 * Init Pentahedron Reference coordinates and Shape function.
 * Case A.
 */
void GaussInfo::penta6aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = PENTA6A_REF[0];
    coords[1] = PENTA6A_REF[1];
    coords[2] = PENTA6A_REF[2];
    break;
  case 1:
    coords[0] = PENTA6A_REF[3];
    coords[1] = PENTA6A_REF[4];
    coords[2] = PENTA6A_REF[5];
    break;
  case 2:
    coords[0] = PENTA6A_REF[6];
    coords[1] = PENTA6A_REF[7];
    coords[2] = PENTA6A_REF[8];
    break;
  case 3:
    coords[0] = PENTA6A_REF[9];
    coords[1] = PENTA6A_REF[10];
    coords[2] = PENTA6A_REF[11];
    break;
  case 4:
    coords[0] = PENTA6A_REF[12];
    coords[1] = PENTA6A_REF[13];
    coords[2] = PENTA6A_REF[14];
    break;
  case 5:
    coords[0] = PENTA6A_REF[15];
    coords[1] = PENTA6A_REF[16];
    coords[2] = PENTA6A_REF[17];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.5*gc[1]*(1.0 - gc[0]);
  funValue[1] = 0.5*gc[2]*(1.0 - gc[0]);
  funValue[2] = 0.5*(1.0 - gc[1] - gc[2])*(1.0 - gc[0]);
  
  funValue[3] = 0.5*gc[1]*(gc[0] + 1.0);
  funValue[4] = 0.5*gc[2]*(gc[0] + 1.0);
  funValue[5] = 0.5*(1.0 - gc[1] - gc[2])*(1.0 + gc[0]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Pentahedron Reference coordinates and Shape function.
 * Case B.
 */
void GaussInfo::penta6bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = PENTA6B_REF[0];
    coords[1] = PENTA6B_REF[1];
    coords[2] = PENTA6B_REF[2];
    break;
  case 1:
    coords[0] = PENTA6B_REF[3];
    coords[1] = PENTA6B_REF[4];
    coords[2] = PENTA6B_REF[5];
    break;
  case 2:
    coords[0] = PENTA6B_REF[6];
    coords[1] = PENTA6B_REF[7];
    coords[2] = PENTA6B_REF[8];
    break;
  case 3:
    coords[0] = PENTA6B_REF[9];
    coords[1] = PENTA6B_REF[10];
    coords[2] = PENTA6B_REF[11];
    break;
  case 4:
    coords[0] = PENTA6B_REF[12];
    coords[1] = PENTA6B_REF[13];
    coords[2] = PENTA6B_REF[14];
    break;
  case 5:
    coords[0] = PENTA6B_REF[15];
    coords[1] = PENTA6B_REF[16];
    coords[2] = PENTA6B_REF[17];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.5*gc[1]*(1.0 - gc[0]);
  funValue[2] = 0.5*gc[2]*(1.0 - gc[0]);
  funValue[1] = 0.5*(1.0 - gc[1] - gc[2])*(1.0 - gc[0]);
  funValue[3] = 0.5*gc[1]*(gc[0] + 1.0);
  funValue[5] = 0.5*gc[2]*(gc[0] + 1.0);
  funValue[4] = 0.5*(1.0 - gc[1] - gc[2])*(1.0 + gc[0]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * This shapefunc map is same as degenerated tria3aInit
 */
void GaussInfo::penta6DegTria3aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  5:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
   LOCAL_COORD_MACRO_END;

   SHAPE_FUN_MACRO_BEGIN;
   funValue[0] = 0.5*(1.0 + gc[1]);
   funValue[1] = -0.5*(gc[0] + gc[1]);
   funValue[2] = 0.5*(1.0 + gc[0]);
   funValue[3] = 0.;
   funValue[4] = 0.;
   funValue[5] = 0.;
   SHAPE_FUN_MACRO_END;
}

/*!
 * This shapefunc map is same as degenerated tria3bInit
 */
void GaussInfo::penta6DegTria3bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
 case  0:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  5:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
   LOCAL_COORD_MACRO_END;

   SHAPE_FUN_MACRO_BEGIN;
   funValue[0] = 1.0 - gc[0] - gc[1];
   funValue[1] = gc[0];
   funValue[2] = gc[1];
   funValue[3] = 0.;
   funValue[4] = 0.;
   funValue[5] = 0.;
   SHAPE_FUN_MACRO_END;
}

/*!
 * Init Pentahedron Reference coordinates and Shape function.
 * Case A.
 */
void GaussInfo::penta15aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = PENTA15A_REF[0];
    coords[1] = PENTA15A_REF[1];
    coords[2] = PENTA15A_REF[2];
    break;
  case 1:
    coords[0] = PENTA15A_REF[3];
    coords[1] = PENTA15A_REF[4];
    coords[2] = PENTA15A_REF[5];
    break;
  case 2:
    coords[0] = PENTA15A_REF[6];
    coords[1] = PENTA15A_REF[7];
    coords[2] = PENTA15A_REF[8];
    break;
  case 3:
    coords[0] = PENTA15A_REF[9];
    coords[1] = PENTA15A_REF[10];
    coords[2] = PENTA15A_REF[11];
    break;
  case 4:
    coords[0] = PENTA15A_REF[12];
    coords[1] = PENTA15A_REF[13];
    coords[2] = PENTA15A_REF[14];
    break;
  case 5:
    coords[0] = PENTA15A_REF[15];
    coords[1] = PENTA15A_REF[16];
    coords[2] = PENTA15A_REF[17];
    break;
  case 6:
    coords[0] = PENTA15A_REF[18];
    coords[1] = PENTA15A_REF[19];
    coords[2] = PENTA15A_REF[20];
    break;
  case 7:
    coords[0] = PENTA15A_REF[21];
    coords[1] = PENTA15A_REF[22];
    coords[2] = PENTA15A_REF[23];
    break;
  case 8:
    coords[0] = PENTA15A_REF[24];
    coords[1] = PENTA15A_REF[25];
    coords[2] = PENTA15A_REF[26];
    break;
  case 9:
    coords[0] = PENTA15A_REF[27];
    coords[1] = PENTA15A_REF[28];
    coords[2] = PENTA15A_REF[29];
    break;
  case 10:
    coords[0] = PENTA15A_REF[30];
    coords[1] = PENTA15A_REF[31];
    coords[2] = PENTA15A_REF[32];
    break;
  case 11:
    coords[0] = PENTA15A_REF[33];
    coords[1] = PENTA15A_REF[34];
    coords[2] = PENTA15A_REF[35];
    break;
  case 12:
    coords[0] = PENTA15A_REF[36];
    coords[1] = PENTA15A_REF[37];
    coords[2] = PENTA15A_REF[38];
    break;
  case 13:
    coords[0] = PENTA15A_REF[39];
    coords[1] = PENTA15A_REF[40];
    coords[2] = PENTA15A_REF[41];
    break;
  case 14:
    coords[0] = PENTA15A_REF[42];
    coords[1] = PENTA15A_REF[43];
    coords[2] = PENTA15A_REF[44];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.5*gc[1]*(1.0 - gc[0])*(2.0*gc[1] - 2.0 - gc[0]);
  funValue[1] = 0.5*gc[2]*(1.0 - gc[0])*(2.0*gc[2] - 2.0 - gc[0]);
  funValue[2] = 0.5*(gc[0] - 1.0)*(1.0 - gc[1] - gc[2])*(gc[0] + 2.0*gc[1] + 2.0*gc[2]);
  
  funValue[3] = 0.5*gc[1]*(1.0 + gc[0])*(2.0*gc[1] - 2.0 + gc[0]);
  funValue[4] = 0.5*gc[2]*(1.0 + gc[0])*(2.0*gc[2] - 2.0 + gc[0]);
  funValue[5] = 0.5*(-gc[0] - 1.0)*(1.0 - gc[1] - gc[2])*(-gc[0] + 2.0*gc[1] + 2.0*gc[2]);
  
  funValue[6] = 2.0*gc[1]*gc[2]*(1.0 - gc[0]);
  funValue[7] = 2.0*gc[2]*(1.0 - gc[1] - gc[2])*(1.0 - gc[0]);
  funValue[8] = 2.0*gc[1]*(1.0 - gc[1] - gc[2])*(1.0 - gc[0]);
  
  funValue[9] = gc[1]*(1.0 - gc[0]*gc[0]);
  funValue[10] = gc[2]*(1.0 - gc[0]*gc[0]);
  funValue[11] = (1.0 - gc[1] - gc[2])*(1.0 - gc[0]*gc[0]);
  
  funValue[12] = 2.0*gc[1]*gc[2]*(1.0 + gc[0]);
  funValue[13] = 2.0*gc[2]*(1.0 - gc[1] - gc[2])*(1.0 + gc[0]);
  funValue[14] = 2.0*gc[1]*(1.0 - gc[1] - gc[2])*(1.0 + gc[0]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Qaudratic Pentahedron Reference coordinates and Shape function.
 * Case B.
 */
void GaussInfo::penta15bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = PENTA15B_REF[0];
    coords[1] = PENTA15B_REF[1];
    coords[2] = PENTA15B_REF[2];
    break;
  case 1:
    coords[0] = PENTA15B_REF[3];
    coords[1] = PENTA15B_REF[4];
    coords[2] = PENTA15B_REF[5];
    break;
  case 2:
    coords[0] = PENTA15B_REF[6];
    coords[1] = PENTA15B_REF[7];
    coords[2] = PENTA15B_REF[8];
    break;
  case 3:
    coords[0] = PENTA15B_REF[9];
    coords[1] = PENTA15B_REF[10];
    coords[2] = PENTA15B_REF[11];
    break;
  case 4:
    coords[0] = PENTA15B_REF[12];
    coords[1] = PENTA15B_REF[13];
    coords[2] = PENTA15B_REF[14];
    break;
  case 5:
    coords[0] = PENTA15B_REF[15];
    coords[1] = PENTA15B_REF[16];
    coords[2] = PENTA15B_REF[17];
    break;
  case 6:
    coords[0] = PENTA15B_REF[18];
    coords[1] = PENTA15B_REF[19];
    coords[2] = PENTA15B_REF[20];
    break;
  case 7:
    coords[0] = PENTA15B_REF[21];
    coords[1] = PENTA15B_REF[22];
    coords[2] = PENTA15B_REF[23];
    break;
  case 8:
    coords[0] = PENTA15B_REF[24];
    coords[1] = PENTA15B_REF[25];
    coords[2] = PENTA15B_REF[26];
    break;
  case 9:
    coords[0] = PENTA15B_REF[27];
    coords[1] = PENTA15B_REF[28];
    coords[2] = PENTA15B_REF[29];
    break;
  case 10:
    coords[0] = PENTA15B_REF[30];
    coords[1] = PENTA15B_REF[31];
    coords[2] = PENTA15B_REF[32];
    break;
  case 11:
    coords[0] = PENTA15B_REF[33];
    coords[1] = PENTA15B_REF[34];
    coords[2] = PENTA15B_REF[35];
    break;
  case 12:
    coords[0] = PENTA15B_REF[36];
    coords[1] = PENTA15B_REF[37];
    coords[2] = PENTA15B_REF[38];
    break;
  case 13:
    coords[0] = PENTA15B_REF[39];
    coords[1] = PENTA15B_REF[40];
    coords[2] = PENTA15B_REF[41];
    break;
  case 14:
    coords[0] = PENTA15B_REF[42];
    coords[1] = PENTA15B_REF[43];
    coords[2] = PENTA15B_REF[44];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.5*gc[1]*(1.0 - gc[0])*(2.0*gc[1] - 2.0 - gc[0]);
  funValue[2] = 0.5*gc[2]*(1.0 - gc[0])*(2.0*gc[2] - 2.0 - gc[0]);
  funValue[1] = 0.5*(gc[0] - 1.0)*(1.0 - gc[1] - gc[2])*(gc[0] + 2.0*gc[1] + 2.0*gc[2]);
  
  funValue[3] = 0.5*gc[1]*(1.0 + gc[0])*(2.0*gc[1] - 2.0 + gc[0]);
  funValue[5] = 0.5*gc[2]*(1.0 + gc[0])*(2.0*gc[2] - 2.0 + gc[0]);
  funValue[4] = 0.5*(-gc[0] - 1.0)*(1.0 - gc[1] - gc[2])*(-gc[0] + 2.0*gc[1] + 2.0*gc[2]);
  
  funValue[8] = 2.0*gc[1]*gc[2]*(1.0 - gc[0]);
  funValue[7] = 2.0*gc[2]*(1.0 - gc[1] - gc[2])*(1.0 - gc[0]);
  funValue[6] = 2.0*gc[1]*(1.0 - gc[1] - gc[2])*(1.0 - gc[0]);
  
  funValue[12] = gc[1]*(1.0 - gc[0]*gc[0]);
  funValue[14] = gc[2]*(1.0 - gc[0]*gc[0]);
  funValue[13] = (1.0 - gc[1] - gc[2])*(1.0 - gc[0]*gc[0]);
  
  funValue[11] = 2.0*gc[1]*gc[2]*(1.0 + gc[0]);
  funValue[10] = 2.0*gc[2]*(1.0 - gc[1] - gc[2])*(1.0 + gc[0]);
  funValue[9]  = 2.0*gc[1]*(1.0 - gc[1] - gc[2])*(1.0 + gc[0]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Hehahedron Reference coordinates and Shape function.
 * Case A.
 */
void GaussInfo::hexa8aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = HEXA8A_REF[0];
    coords[1] = HEXA8A_REF[1];
    coords[2] = HEXA8A_REF[2];
    break;
  case 1:
    coords[0] = HEXA8A_REF[3];
    coords[1] = HEXA8A_REF[4];
    coords[2] = HEXA8A_REF[5];
    break;
  case 2:
    coords[0] = HEXA8A_REF[6];
    coords[1] = HEXA8A_REF[7];
    coords[2] = HEXA8A_REF[8];
    break;
  case 3:
    coords[0] = HEXA8A_REF[9];
    coords[1] = HEXA8A_REF[10];
    coords[2] = HEXA8A_REF[11];
    break;
  case 4:
    coords[0] = HEXA8A_REF[12];
    coords[1] = HEXA8A_REF[13];
    coords[2] = HEXA8A_REF[14];
    break;
  case 5:
    coords[0] = HEXA8A_REF[15];
    coords[1] = HEXA8A_REF[16];
    coords[2] = HEXA8A_REF[17];
    break;
  case 6:
    coords[0] = HEXA8A_REF[18];
    coords[1] = HEXA8A_REF[19];
    coords[2] = HEXA8A_REF[20];
    break;
  case 7:
    coords[0] = HEXA8A_REF[21];
    coords[1] = HEXA8A_REF[22];
    coords[2] = HEXA8A_REF[23];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.125*(1.0 - gc[0])*(1.0 - gc[1])*(1.0 - gc[2]);
  funValue[1] = 0.125*(1.0 + gc[0])*(1.0 - gc[1])*(1.0 - gc[2]);
  funValue[2] = 0.125*(1.0 + gc[0])*(1.0 + gc[1])*(1.0 - gc[2]);
  funValue[3] = 0.125*(1.0 - gc[0])*(1.0 + gc[1])*(1.0 - gc[2]);
  
  funValue[4] = 0.125*(1.0 - gc[0])*(1.0 - gc[1])*(1.0 + gc[2]);
  funValue[5] = 0.125*(1.0 + gc[0])*(1.0 - gc[1])*(1.0 + gc[2]);
  funValue[6] = 0.125*(1.0 + gc[0])*(1.0 + gc[1])*(1.0 + gc[2]);
  funValue[7] = 0.125*(1.0 - gc[0])*(1.0 + gc[1])*(1.0 + gc[2]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Hehahedron Reference coordinates and Shape function.
 * Case B.
 */
void GaussInfo::hexa8bInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = HEXA8B_REF[0];
    coords[1] = HEXA8B_REF[1];
    coords[2] = HEXA8B_REF[2];
    break;
  case 1:
    coords[0] = HEXA8B_REF[3];
    coords[1] = HEXA8B_REF[4];
    coords[2] = HEXA8B_REF[5];
    break;
  case 2:
    coords[0] = HEXA8B_REF[6];
    coords[1] = HEXA8B_REF[7];
    coords[2] = HEXA8B_REF[8];
    break;
  case 3:
    coords[0] = HEXA8B_REF[9];
    coords[1] = HEXA8B_REF[10];
    coords[2] = HEXA8B_REF[11];
    break;
  case 4:
    coords[0] = HEXA8B_REF[12];
    coords[1] = HEXA8B_REF[13];
    coords[2] = HEXA8B_REF[14];
    break;
  case 5:
    coords[0] = HEXA8B_REF[15];
    coords[1] = HEXA8B_REF[16];
    coords[2] = HEXA8B_REF[17];
    break;
  case 6:
    coords[0] = HEXA8B_REF[18];
    coords[1] = HEXA8B_REF[19];
    coords[2] = HEXA8B_REF[20];
    break;
  case 7:
    coords[0] = HEXA8B_REF[21];
    coords[1] = HEXA8B_REF[22];
    coords[2] = HEXA8B_REF[23];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.125*(1.0 - gc[0])*(1.0 - gc[1])*(1.0 - gc[2]);
  funValue[3] = 0.125*(1.0 + gc[0])*(1.0 - gc[1])*(1.0 - gc[2]);
  funValue[2] = 0.125*(1.0 + gc[0])*(1.0 + gc[1])*(1.0 - gc[2]);
  funValue[1] = 0.125*(1.0 - gc[0])*(1.0 + gc[1])*(1.0 - gc[2]);
  
  funValue[4] = 0.125*(1.0 - gc[0])*(1.0 - gc[1])*(1.0 + gc[2]);
  funValue[7] = 0.125*(1.0 + gc[0])*(1.0 - gc[1])*(1.0 + gc[2]);
  funValue[6] = 0.125*(1.0 + gc[0])*(1.0 + gc[1])*(1.0 + gc[2]);
  funValue[5] = 0.125*(1.0 - gc[0])*(1.0 + gc[1])*(1.0 + gc[2]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * This shapefunc map is same as degenerated quad4bInit
 */
void GaussInfo::hexa8DegQuad4aInit()
{
  LOCAL_COORD_MACRO_BEGIN;
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  5:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  6:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  7:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.25*(1.0 + gc[1])*(1.0 - gc[0]);
  funValue[1] = 0.25*(1.0 - gc[1])*(1.0 - gc[0]);
  funValue[2] = 0.25*(1.0 - gc[1])*(1.0 + gc[0]);
  funValue[3] = 0.25*(1.0 + gc[0])*(1.0 + gc[1]);
  funValue[4] = 0.;
  funValue[5] = 0.;
  funValue[6] = 0.;
  funValue[7] = 0.;
  SHAPE_FUN_MACRO_END;
}

/*!
 * This shapefunc map is same as degenerated quad4bInit
 */
void GaussInfo::hexa8DegQuad4bInit()
{
  LOCAL_COORD_MACRO_BEGIN;
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  5:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  6:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  7:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
   LOCAL_COORD_MACRO_END;

   SHAPE_FUN_MACRO_BEGIN;
   funValue[0] = 0.25*(1.0 - gc[0])*(1.0 - gc[1]);
   funValue[1] = 0.25*(1.0 + gc[0])*(1.0 - gc[1]);
   funValue[2] = 0.25*(1.0 + gc[0])*(1.0 + gc[1]);
   funValue[3] = 0.25*(1.0 - gc[0])*(1.0 + gc[1]);
   funValue[4] = 0.;
   funValue[5] = 0.;
   funValue[6] = 0.;
   funValue[7] = 0.;
   SHAPE_FUN_MACRO_END;
}

/*!
 * This shapefunc map is same as degenerated quad4cInit
 */
void GaussInfo::hexa8DegQuad4cInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  5:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  6:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  7:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;

   LOCAL_COORD_MACRO_END;

   SHAPE_FUN_MACRO_BEGIN;
   funValue[0] = 0.25*(1.0 - gc[0])*(1.0 - gc[1]);
   funValue[1] = 0.25*(1.0 - gc[0])*(1.0 + gc[1]);
   funValue[2] = 0.25*(1.0 + gc[0])*(1.0 + gc[1]);
   funValue[3] = 0.25*(1.0 + gc[0])*(1.0 - gc[1]);
   funValue[4] = 0. ;
   funValue[5] = 0. ;
   funValue[6] = 0. ;
   funValue[7] = 0. ;
   SHAPE_FUN_MACRO_END;
}

/*!
 * Init Qaudratic Hehahedron Reference coordinates and Shape function.
 * Case A.
 */
void GaussInfo::hexa20aInit() 
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = HEXA20A_REF[0];
    coords[1] = HEXA20A_REF[1];
    coords[2] = HEXA20A_REF[2];
    break;
  case 1:
    coords[0] = HEXA20A_REF[3];
    coords[1] = HEXA20A_REF[4];
    coords[2] = HEXA20A_REF[5];
    break;
  case 2:
    coords[0] = HEXA20A_REF[6];
    coords[1] = HEXA20A_REF[7];
    coords[2] = HEXA20A_REF[8];
    break;
  case 3:
    coords[0] = HEXA20A_REF[9];
    coords[1] = HEXA20A_REF[10];
    coords[2] = HEXA20A_REF[11];
    break;
  case 4:
    coords[0] = HEXA20A_REF[12];
    coords[1] = HEXA20A_REF[13];
    coords[2] = HEXA20A_REF[14];
    break;
  case 5:
    coords[0] = HEXA20A_REF[15];
    coords[1] = HEXA20A_REF[16];
    coords[2] = HEXA20A_REF[17];
    break;
  case 6:
    coords[0] = HEXA20A_REF[18];
    coords[1] = HEXA20A_REF[19];
    coords[2] = HEXA20A_REF[20];
    break;
  case 7:
    coords[0] = HEXA20A_REF[21];
    coords[1] = HEXA20A_REF[22];
    coords[2] = HEXA20A_REF[23];
    break;
  case 8:
    coords[0] = HEXA20A_REF[24];
    coords[1] = HEXA20A_REF[25];
    coords[2] = HEXA20A_REF[26];
    break;
  case 9:
    coords[0] = HEXA20A_REF[27];
    coords[1] = HEXA20A_REF[28];
    coords[2] = HEXA20A_REF[29];
    break;
  case 10:
    coords[0] = HEXA20A_REF[30];
    coords[1] = HEXA20A_REF[31];
    coords[2] = HEXA20A_REF[32];
    break;
  case 11:
    coords[0] = HEXA20A_REF[33];
    coords[1] = HEXA20A_REF[34];
    coords[2] = HEXA20A_REF[35];
    break;
  case 12:
    coords[0] = HEXA20A_REF[36];
    coords[1] = HEXA20A_REF[37];
    coords[2] = HEXA20A_REF[38];
    break;
  case 13:
    coords[0] = HEXA20A_REF[39];
    coords[1] = HEXA20A_REF[40];
    coords[2] = HEXA20A_REF[41];
    break;
  case 14:
    coords[0] = HEXA20A_REF[42];
    coords[1] = HEXA20A_REF[43];
    coords[2] = HEXA20A_REF[44];
    break;
  case 15:
    coords[0] = HEXA20A_REF[45];
    coords[1] = HEXA20A_REF[46];
    coords[2] = HEXA20A_REF[47];
    break;
  case 16:
    coords[0] = HEXA20A_REF[48];
    coords[1] = HEXA20A_REF[49];
    coords[2] = HEXA20A_REF[50];
    break;
  case 17:
    coords[0] = HEXA20A_REF[51];
    coords[1] = HEXA20A_REF[52];
    coords[2] = HEXA20A_REF[53];
    break;
  case 18:
    coords[0] = HEXA20A_REF[54];
    coords[1] = HEXA20A_REF[55];
    coords[2] = HEXA20A_REF[56];
    break;
  case 19:
    coords[0] = HEXA20A_REF[57];
    coords[1] = HEXA20A_REF[58];
    coords[2] = HEXA20A_REF[59];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  funValue[0] = 0.125*(1.0 - gc[0])*(1.0 - gc[1])*(1.0 - gc[2])*
    (-2.0 - gc[0] - gc[1] - gc[2]);
  funValue[1] = 0.125*(1.0 + gc[0])*(1.0 - gc[1])*(1.0 - gc[2])*
    (-2.0 + gc[0] - gc[1] - gc[2]);
  funValue[2] = 0.125*(1.0 + gc[0])*(1.0 + gc[1])*(1.0 - gc[2])*
    (-2.0 + gc[0] + gc[1] - gc[2]);
  funValue[3] = 0.125*(1.0 - gc[0])*(1.0 + gc[1])*(1.0 - gc[2])*
    (-2.0 - gc[0] + gc[1] - gc[2]);
  funValue[4] = 0.125*(1.0 - gc[0])*(1.0 - gc[1])*(1.0 + gc[2])*
    (-2.0 - gc[0] - gc[1] + gc[2]);
  funValue[5] = 0.125*(1.0 + gc[0])*(1.0 - gc[1])*(1.0 + gc[2])*
    (-2.0 + gc[0] - gc[1] + gc[2]);
  funValue[6] = 0.125*(1.0 + gc[0])*(1.0 + gc[1])*(1.0 + gc[2])*
    (-2.0 + gc[0] + gc[1] + gc[2]);
  funValue[7] = 0.125*(1.0 - gc[0])*(1.0 + gc[1])*(1.0 + gc[2])*
    (-2.0 - gc[0] + gc[1] + gc[2]);
  
  funValue[8] = 0.25*(1.0 - gc[0]*gc[0])*(1.0 - gc[1])*(1.0 - gc[2]);
  funValue[9] = 0.25*(1.0 - gc[1]*gc[1])*(1.0 + gc[0])*(1.0 - gc[2]);
  funValue[10] = 0.25*(1.0 - gc[0]*gc[0])*(1.0 + gc[1])*(1.0 - gc[2]);
  funValue[11] = 0.25*(1.0 - gc[1]*gc[1])*(1.0 - gc[0])*(1.0 - gc[2]);
  funValue[12] = 0.25*(1.0 - gc[2]*gc[2])*(1.0 - gc[0])*(1.0 - gc[1]);
  funValue[13] = 0.25*(1.0 - gc[2]*gc[2])*(1.0 + gc[0])*(1.0 - gc[1]);
  funValue[14] = 0.25*(1.0 - gc[2]*gc[2])*(1.0 + gc[0])*(1.0 + gc[1]);
  funValue[15] = 0.25*(1.0 - gc[2]*gc[2])*(1.0 - gc[0])*(1.0 + gc[1]);
  funValue[16] = 0.25*(1.0 - gc[0]*gc[0])*(1.0 - gc[1])*(1.0 + gc[2]);
  funValue[17] = 0.25*(1.0 - gc[1]*gc[1])*(1.0 + gc[0])*(1.0 + gc[2]);
  funValue[18] = 0.25*(1.0 - gc[0]*gc[0])*(1.0 + gc[1])*(1.0 + gc[2]);
  funValue[19] = 0.25*(1.0 - gc[1]*gc[1])*(1.0 - gc[0])*(1.0 + gc[2]);
  SHAPE_FUN_MACRO_END;
}

/*!
 * Init Qaudratic Hehahedron Reference coordinates and Shape function.
 * Case B.
 */
void GaussInfo::hexa20bInit()
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = HEXA20B_REF[0];
    coords[1] = HEXA20B_REF[1];
    coords[2] = HEXA20B_REF[2];
    break;
  case 1:
    coords[0] = HEXA20B_REF[3];
    coords[1] = HEXA20B_REF[4];
    coords[2] = HEXA20B_REF[5];
    break;
  case 2:
    coords[0] = HEXA20B_REF[6];
    coords[1] = HEXA20B_REF[7];
    coords[2] = HEXA20B_REF[8];
    break;
  case 3:
    coords[0] = HEXA20B_REF[9];
    coords[1] = HEXA20B_REF[10];
    coords[2] = HEXA20B_REF[11];
    break;
  case 4:
    coords[0] = HEXA20B_REF[12];
    coords[1] = HEXA20B_REF[13];
    coords[2] = HEXA20B_REF[14];
    break;
  case 5:
    coords[0] = HEXA20B_REF[15];
    coords[1] = HEXA20B_REF[16];
    coords[2] = HEXA20B_REF[17];
    break;
  case 6:
    coords[0] = HEXA20B_REF[18];
    coords[1] = HEXA20B_REF[19];
    coords[2] = HEXA20B_REF[20];
    break;
  case 7:
    coords[0] = HEXA20B_REF[21];
    coords[1] = HEXA20B_REF[22];
    coords[2] = HEXA20B_REF[23];
    break;
  case 8:
    coords[0] = HEXA20B_REF[24];
    coords[1] = HEXA20B_REF[25];
    coords[2] = HEXA20B_REF[26];
    break;
  case 9:
    coords[0] = HEXA20B_REF[27];
    coords[1] = HEXA20B_REF[28];
    coords[2] = HEXA20B_REF[29];
    break;
  case 10:
    coords[0] = HEXA20B_REF[30];
    coords[1] = HEXA20B_REF[31];
    coords[2] = HEXA20B_REF[32];
    break;
  case 11:
    coords[0] = HEXA20B_REF[33];
    coords[1] = HEXA20B_REF[34];
    coords[2] = HEXA20B_REF[35];
    break;
  case 12:
    coords[0] = HEXA20B_REF[36];
    coords[1] = HEXA20B_REF[37];
    coords[2] = HEXA20B_REF[38];
    break;
  case 13:
    coords[0] = HEXA20B_REF[39];
    coords[1] = HEXA20B_REF[40];
    coords[2] = HEXA20B_REF[41];
    break;
  case 14:
    coords[0] = HEXA20B_REF[42];
    coords[1] = HEXA20B_REF[43];
    coords[2] = HEXA20B_REF[44];
    break;
  case 15:
    coords[0] = HEXA20B_REF[45];
    coords[1] = HEXA20B_REF[46];
    coords[2] = HEXA20B_REF[47];
    break;
  case 16:
    coords[0] = HEXA20B_REF[48];
    coords[1] = HEXA20B_REF[49];
    coords[2] = HEXA20B_REF[50];
    break;
  case 17:
    coords[0] = HEXA20B_REF[51];
    coords[1] = HEXA20B_REF[52];
    coords[2] = HEXA20B_REF[53];
    break;
  case 18:
    coords[0] = HEXA20B_REF[54];
    coords[1] = HEXA20B_REF[55];
    coords[2] = HEXA20B_REF[56];
    break;
  case 19:
    coords[0] = HEXA20B_REF[57];
    coords[1] = HEXA20B_REF[58];
    coords[2] = HEXA20B_REF[59];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  
  funValue[0] = 0.125*(1.0 - gc[0])*(1.0 - gc[1])*(1.0 - gc[2])*
    (-2.0 - gc[0] - gc[1] - gc[2]);
  funValue[3] = 0.125*(1.0 + gc[0])*(1.0 - gc[1])*(1.0 - gc[2])*
    (-2.0 + gc[0] - gc[1] - gc[2]);
  funValue[2] = 0.125*(1.0 + gc[0])*(1.0 + gc[1])*(1.0 - gc[2])*
    (-2.0 + gc[0] + gc[1] - gc[2]);
  funValue[1] = 0.125*(1.0 - gc[0])*(1.0 + gc[1])*(1.0 - gc[2])*
    (-2.0 - gc[0] + gc[1] - gc[2]);
  funValue[4] = 0.125*(1.0 - gc[0])*(1.0 - gc[1])*(1.0 + gc[2])*
    (-2.0 - gc[0] - gc[1] + gc[2]);
  funValue[7] = 0.125*(1.0 + gc[0])*(1.0 - gc[1])*(1.0 + gc[2])*
    (-2.0 + gc[0] - gc[1] + gc[2]);
  funValue[6] = 0.125*(1.0 + gc[0])*(1.0 + gc[1])*(1.0 + gc[2])*
    (-2.0 + gc[0] + gc[1] + gc[2]);
  funValue[5] = 0.125*(1.0 - gc[0])*(1.0 + gc[1])*(1.0 + gc[2])*
    (-2.0 - gc[0] + gc[1] + gc[2]);
  
  funValue[11] = 0.25*(1.0 - gc[0]*gc[0])*(1.0 - gc[1])*(1.0 - gc[2]);
  funValue[10] = 0.25*(1.0 - gc[1]*gc[1])*(1.0 + gc[0])*(1.0 - gc[2]);
  funValue[9] = 0.25*(1.0 - gc[0]*gc[0])*(1.0 + gc[1])*(1.0 - gc[2]);
  funValue[8] = 0.25*(1.0 - gc[1]*gc[1])*(1.0 - gc[0])*(1.0 - gc[2]);
  funValue[16] = 0.25*(1.0 - gc[2]*gc[2])*(1.0 - gc[0])*(1.0 - gc[1]);
  funValue[19] = 0.25*(1.0 - gc[2]*gc[2])*(1.0 + gc[0])*(1.0 - gc[1]);
  funValue[18] = 0.25*(1.0 - gc[2]*gc[2])*(1.0 + gc[0])*(1.0 + gc[1]);
  funValue[17] = 0.25*(1.0 - gc[2]*gc[2])*(1.0 - gc[0])*(1.0 + gc[1]);
  funValue[15] = 0.25*(1.0 - gc[0]*gc[0])*(1.0 - gc[1])*(1.0 + gc[2]);
  funValue[14] = 0.25*(1.0 - gc[1]*gc[1])*(1.0 + gc[0])*(1.0 + gc[2]);
  funValue[13] = 0.25*(1.0 - gc[0]*gc[0])*(1.0 + gc[1])*(1.0 + gc[2]);
  funValue[12] = 0.25*(1.0 - gc[1]*gc[1])*(1.0 - gc[0])*(1.0 + gc[2]);
  SHAPE_FUN_MACRO_END;
}

void GaussInfo::hexa27aInit()
{
  LOCAL_COORD_MACRO_BEGIN;
  case 0:
    coords[0] = HEXA27A_REF[0];
    coords[1] = HEXA27A_REF[1];
    coords[2] = HEXA27A_REF[2];
    break;
  case 1:
    coords[0] = HEXA27A_REF[3];
    coords[1] = HEXA27A_REF[4];
    coords[2] = HEXA27A_REF[5];
    break;
  case 2:
    coords[0] = HEXA27A_REF[6];
    coords[1] = HEXA27A_REF[7];
    coords[2] = HEXA27A_REF[8];
    break;
  case 3:
    coords[0] = HEXA27A_REF[9];
    coords[1] = HEXA27A_REF[10];
    coords[2] = HEXA27A_REF[11];
    break;
  case 4:
    coords[0] = HEXA27A_REF[12];
    coords[1] = HEXA27A_REF[13];
    coords[2] = HEXA27A_REF[14];
    break;
  case 5:
    coords[0] = HEXA27A_REF[15];
    coords[1] = HEXA27A_REF[16];
    coords[2] = HEXA27A_REF[17];
    break;
  case 6:
    coords[0] = HEXA27A_REF[18];
    coords[1] = HEXA27A_REF[19];
    coords[2] = HEXA27A_REF[20];
    break;
  case 7:
    coords[0] = HEXA27A_REF[21];
    coords[1] = HEXA27A_REF[22];
    coords[2] = HEXA27A_REF[23];
    break;
  case 8:
    coords[0] = HEXA27A_REF[24];
    coords[1] = HEXA27A_REF[25];
    coords[2] = HEXA27A_REF[26];
    break;
  case 9:
    coords[0] = HEXA27A_REF[27];
    coords[1] = HEXA27A_REF[28];
    coords[2] = HEXA27A_REF[29];
    break;
  case 10:
    coords[0] = HEXA27A_REF[30];
    coords[1] = HEXA27A_REF[31];
    coords[2] = HEXA27A_REF[32];
    break;
  case 11:
    coords[0] = HEXA27A_REF[33];
    coords[1] = HEXA27A_REF[34];
    coords[2] = HEXA27A_REF[35];
    break;
  case 12:
    coords[0] = HEXA27A_REF[36];
    coords[1] = HEXA27A_REF[37];
    coords[2] = HEXA27A_REF[38];
    break;
  case 13:
    coords[0] = HEXA27A_REF[39];
    coords[1] = HEXA27A_REF[40];
    coords[2] = HEXA27A_REF[41];
    break;
  case 14:
    coords[0] = HEXA27A_REF[42];
    coords[1] = HEXA27A_REF[43];
    coords[2] = HEXA27A_REF[44];
    break;
  case 15:
    coords[0] = HEXA27A_REF[45];
    coords[1] = HEXA27A_REF[46];
    coords[2] = HEXA27A_REF[47];
    break;
  case 16:
    coords[0] = HEXA27A_REF[48];
    coords[1] = HEXA27A_REF[49];
    coords[2] = HEXA27A_REF[50];
    break;
  case 17:
    coords[0] = HEXA27A_REF[51];
    coords[1] = HEXA27A_REF[52];
    coords[2] = HEXA27A_REF[53];
    break;
  case 18:
    coords[0] = HEXA27A_REF[54];
    coords[1] = HEXA27A_REF[55];
    coords[2] = HEXA27A_REF[56];
    break;
  case 19:
    coords[0] = HEXA27A_REF[57];
    coords[1] = HEXA27A_REF[58];
    coords[2] = HEXA27A_REF[59];
    break;
  case 20:
    coords[0] = HEXA27A_REF[60];
    coords[1] = HEXA27A_REF[61];
    coords[2] = HEXA27A_REF[62];
    break;
  case 21:
    coords[0] = HEXA27A_REF[63];
    coords[1] = HEXA27A_REF[64];
    coords[2] = HEXA27A_REF[65];
    break;
  case 22:
    coords[0] = HEXA27A_REF[66];
    coords[1] = HEXA27A_REF[67];
    coords[2] = HEXA27A_REF[68];
    break;
  case 23:
    coords[0] = HEXA27A_REF[69];
    coords[1] = HEXA27A_REF[70];
    coords[2] = HEXA27A_REF[71];
    break;
  case 24:
    coords[0] = HEXA27A_REF[72];
    coords[1] = HEXA27A_REF[73];
    coords[2] = HEXA27A_REF[74];
    break;
  case 25:
    coords[0] = HEXA27A_REF[75];
    coords[1] = HEXA27A_REF[76];
    coords[2] = HEXA27A_REF[77];
    break;
  case 26:
    coords[0] = HEXA27A_REF[78];
    coords[1] = HEXA27A_REF[79];
    coords[2] = HEXA27A_REF[80];
    break;
  LOCAL_COORD_MACRO_END;
  
  SHAPE_FUN_MACRO_BEGIN;
  
  funValue[0] =0.125*gc[0]*(gc[0]-1.)*gc[1]*(gc[1]-1.)*gc[2]*(gc[2]-1.);
  funValue[1] =0.125*gc[0]*(gc[0]-1.)*gc[1]*(gc[1]+1.)*gc[2]*(gc[2]-1.);
  funValue[2] =0.125*gc[0]*(gc[0]+1.)*gc[1]*(gc[1]+1.)*gc[2]*(gc[2]-1.);
  funValue[3] =0.125*gc[0]*(gc[0]+1.)*gc[1]*(gc[1]-1.)*gc[2]*(gc[2]-1.);
  funValue[4] =0.125*gc[0]*(gc[0]-1.)*gc[1]*(gc[1]-1.)*gc[2]*(gc[2]+1.);
  funValue[5] =0.125*gc[0]*(gc[0]-1.)*gc[1]*(gc[1]+1.)*gc[2]*(gc[2]+1.);
  funValue[6] =0.125*gc[0]*(gc[0]+1.)*gc[1]*(gc[1]+1.)*gc[2]*(gc[2]+1.);
  funValue[7] =0.125*gc[0]*(gc[0]+1.)*gc[1]*(gc[1]-1.)*gc[2]*(gc[2]+1.);
  funValue[8] =0.25*gc[0]*(gc[0]-1.)*(1.-gc[1]*gc[1])*gc[2]*(gc[2]-1.);
  funValue[9] =0.25*(1.-gc[0]*gc[0])*gc[1]*(gc[1]+1.)*gc[2]*(gc[2]-1.);
  funValue[10]=0.25*gc[0]*(gc[0]+1.)*(1.-gc[1]*gc[1])*gc[2]*(gc[2]-1.);
  funValue[11]=0.25*(1.-gc[0]*gc[0])*gc[1]*(gc[1]-1.)*gc[2]*(gc[2]-1.);
  funValue[12]=0.25*gc[0]*(gc[0]-1.)*(1.-gc[1]*gc[1])*gc[2]*(gc[2]+1.);
  funValue[13]=0.25*(1.-gc[0]*gc[0])*gc[1]*(gc[1]+1.)*gc[2]*(gc[2]+1.);
  funValue[14]=0.25*gc[0]*(gc[0]+1.)*(1.-gc[1]*gc[1])*gc[2]*(gc[2]+1.);
  funValue[15]=0.25*(1.-gc[0]*gc[0])*gc[1]*(gc[1]-1.)*gc[2]*(gc[2]+1.);
  funValue[16]=0.25*gc[0]*(gc[0]-1.)*gc[1]*(gc[1]-1.)*(1.-gc[2]*gc[2]);
  funValue[17]=0.25*gc[0]*(gc[0]-1.)*gc[1]*(gc[1]+1.)*(1.-gc[2]*gc[2]);
  funValue[18]=0.25*gc[0]*(gc[0]+1.)*gc[1]*(gc[1]+1.)*(1.-gc[2]*gc[2]);
  funValue[19]=0.25*gc[0]*(gc[0]+1.)*gc[1]*(gc[1]-1.)*(1.-gc[2]*gc[2]);
  funValue[20]=0.5*(1.-gc[0]*gc[0])*(1.-gc[1]*gc[1])*gc[2]*(gc[2]-1.);
  funValue[21]=0.5*gc[0]*(gc[0]-1.)*(1.-gc[1]*gc[1])*(1.-gc[2]*gc[2]);
  funValue[22]=0.5*(1.-gc[0]*gc[0])*gc[1]*(gc[1]+1.)*(1.-gc[2]*gc[2]);
  funValue[23]=0.5*gc[0]*(gc[0]+1.)*(1.-gc[1]*gc[1])*(1.-gc[2]*gc[2]);
  funValue[24]=0.5*(1.-gc[0]*gc[0])*gc[1]*(gc[1]-1.)*(1.-gc[2]*gc[2]);
  funValue[25]=0.5*(1.-gc[0]*gc[0])*(1.-gc[1]*gc[1])*gc[2]*(gc[2]+1.);
  funValue[26]=(1.-gc[0]*gc[0])*(1.-gc[1]*gc[1])*(1.-gc[2]*gc[2]);
  
  SHAPE_FUN_MACRO_END;
}

////////////////////////////////////////////////////////////////////////////////////////////////
//                                GAUSS COORD CLASS                                           //
////////////////////////////////////////////////////////////////////////////////////////////////
/*!
 * Constructor
 */
GaussCoords::GaussCoords()
{
}

/*!
 * Destructor
 */
GaussCoords::~GaussCoords()
{
  GaussInfoVector::iterator it = _my_gauss_info.begin();
  for( ; it != _my_gauss_info.end(); it++ ) 
    {
      if((*it) != NULL)
        delete (*it);
    }
}

/*!
 * Add Gauss localization info 
 */
void GaussCoords::addGaussInfo( NormalizedCellType theGeometry,
                                int coordDim,
                                const double* theGaussCoord,
                                int theNbGauss,
                                const double* theReferenceCoord,
                                int theNbRef)
{
  GaussInfoVector::iterator it = _my_gauss_info.begin();
  for( ; it != _my_gauss_info.end(); it++ ) 
    {
      if( (*it)->getCellType() == theGeometry ) 
        {
          break;
        }
    }

  DataVector aGaussCoord;
  for(int i = 0 ; i < theNbGauss*coordDim; i++ )
    aGaussCoord.push_back(theGaussCoord[i]);

  DataVector aReferenceCoord;
  for(int i = 0 ; i < theNbRef*coordDim; i++ )
    aReferenceCoord.push_back(theReferenceCoord[i]);


  GaussInfo* info = new GaussInfo( theGeometry, aGaussCoord, theNbGauss, aReferenceCoord, theNbRef);
  info->initLocalInfo();

  //If info with cell type doesn't exist add it
  if( it == _my_gauss_info.end() ) 
    {
      _my_gauss_info.push_back(info);

      // If information exists, update it
    }
  else 
    {
      int index = std::distance(_my_gauss_info.begin(),it);
      delete (*it);
      _my_gauss_info[index] = info;
    }
}


/*!
 * Calculate gauss points coordinates
 */
double* GaussCoords::calculateCoords( NormalizedCellType theGeometry, 
                                      const double *theNodeCoords, 
                                      const int theSpaceDim,
                                      const int *theIndex)
{
  const GaussInfo *info = getInfoGivenCellType(theGeometry);
  int nbCoords = theSpaceDim * info->getNbGauss();
  double *aCoords = new double[nbCoords];
  calculateCoordsAlg(info,theNodeCoords,theSpaceDim,theIndex,aCoords);
  return aCoords;
}


void GaussCoords::calculateCoords( NormalizedCellType theGeometry, const double *theNodeCoords, const int theSpaceDim, const int *theIndex, double *result)
{
  const GaussInfo *info = getInfoGivenCellType(theGeometry);
  calculateCoordsAlg(info,theNodeCoords,theSpaceDim,theIndex,result);
}

void GaussCoords::calculateCoordsAlg(const GaussInfo *info, const double* theNodeCoords, const int theSpaceDim, const int *theIndex, double *result)
{
  int aConn = info->getNbRef();

  int nbCoords = theSpaceDim * info->getNbGauss();
  std::fill(result,result+nbCoords,0.);

  for( int gaussId = 0; gaussId < info->getNbGauss(); gaussId++ ) 
    {
      double *coord=result+gaussId*theSpaceDim;
      const double *function=info->getFunctionValues(gaussId);
      for ( int connId = 0; connId < aConn ; connId++ ) 
        {
          const double* nodeCoord = theNodeCoords + (theIndex[connId]*theSpaceDim);
          for( int dimId = 0; dimId < theSpaceDim; dimId++ )
            coord[dimId] += nodeCoord[dimId]*function[connId];
        }
    }
}

const GaussInfo *GaussCoords::getInfoGivenCellType(NormalizedCellType cellType)
{
  GaussInfoVector::const_iterator it = _my_gauss_info.begin();
  //Try to find gauss localization info
  for( ; it != _my_gauss_info.end() ; it++ ) 
    if( (*it)->getCellType()==cellType) 
      return (*it);
  throw INTERP_KERNEL::Exception("Can't find gauss localization information !");
}
