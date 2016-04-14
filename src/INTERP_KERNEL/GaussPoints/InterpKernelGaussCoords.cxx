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
  case  0:
    coords[0] = -1.0;
    break;
  case  1:
    coords[0] =  1.0;
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
  case  0:
    coords[0] = -1.0;
    break;
  case  1:
    coords[0] =  1.0;
    break;
  case  2:
    coords[0] =  0.0;
    break;
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
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] = -1.0;
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
 case  0:
   coords[0] =  0.0;
   coords[1] =  0.0;
   break;
 case  1:
   coords[0] =  1.0;
   coords[1] =  0.0;
   break;
 case  2:
   coords[0] =  0.0;
   coords[1] =  1.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] = -1.0;
   break;
 case  3:
   coords[0] = -1.0;
   coords[1] =  1.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] = -1.0;
   break;
 case  5:
   coords[0] =  0.0;
   coords[1] =  0.0;
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
 case  0:
   coords[0] =  0.0;
   coords[1] =  0.0;
   break;

 case  1:
   coords[0] =  1.0;
   coords[1] =  0.0;
   break;

 case  2:
   coords[0] =  0.0;
   coords[1] =  1.0;
   break;

 case  3:
   coords[0] =  0.5;
   coords[1] =  0.0;
   break;

 case  4:
   coords[0] =  0.5;
   coords[1] =  0.5;
   break;

 case  5:
   coords[0] =  0.0;
   coords[1] =  0.5;
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
   coords[0] = 0.0;
   coords[1] = 0.0;
   break;
 case 1:
   coords[0] = 1.0;
   coords[1] = 0.0;
   break;
 case 2:
   coords[0] = 0.0;
   coords[1] = 1.0;
   break;
 case 3:
   coords[0] = 0.5;
   coords[1] = 0.0;
   break;
 case 4:
   coords[0] = 0.5;
   coords[1] = 0.5;
   break;
 case 5:
   coords[0] = 0.0;
   coords[1] = 0.5;
   break;
 case 6:
   coords[0] = 0.3333333333333333;
   coords[1] = 0.3333333333333333;
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
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] = -1.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  1.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   break;
 case  1:
   coords[0] =  1.0;
   coords[1] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   break;
 case  3:
   coords[0] = -1.0;
   coords[1] =  1.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] = -1.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  1.0;
   break;
 case  4:
   coords[0] = -1.0;
   coords[1] =  0.0;
   break;
 case  5:
   coords[0] =  0.0;
   coords[1] = -1.0;
   break;
 case  6:
   coords[0] =  1.0;
   coords[1] =  0.0;
   break;
 case  7:
   coords[0] =  0.0;
   coords[1] =  1.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   break;
 case  1:
   coords[0] =  1.0;
   coords[1] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   break;
 case  3:
   coords[0] = -1.0;
   coords[1] =  1.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] = -1.0;
   break;
 case  5:
   coords[0] =  1.0;
   coords[1] =  0.0;
   break;
 case  6:
   coords[0] =  0.0;
   coords[1] =  1.0;
   break;
 case  7:
   coords[0] = -1.0;
   coords[1] =  0.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   break;
 case  1:
   coords[0] =  1.0;
   coords[1] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   break;
 case  3:
   coords[0] = -1.0;
   coords[1] =  1.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] = -1.0;
   break;
 case  5:
   coords[0] =  1.0;
   coords[1] =  0.0;
   break;
 case  6:
   coords[0] =  0.0;
   coords[1] =  1.0;
   break;
 case  7:
   coords[0] = -1.0;
   coords[1] =  0.0;
   break;
 case  8:
   coords[0] =  0.0;
   coords[1] =  0.0;
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
 case  0:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case  2:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
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
 case  0:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case  1:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
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
 case  0:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case  2:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.5;
   coords[2] =  0.5;
   break;
 case  5:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case  6:
   coords[0] =  0.0;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  7:
   coords[0] =  0.5;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  8:
   coords[0] =  0.5;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case  9:
   coords[0] =  0.5;
   coords[1] =  0.0;
   coords[2] =  0.0;
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
 case  0:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case  1:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  6:
   coords[0] =  0.0;
   coords[1] =  0.5;
   coords[2] =  0.5;
   break;
 case  5:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  7:
   coords[0] =  0.5;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  9:
   coords[0] =  0.5;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case  8:
   coords[0] =  0.5;
   coords[1] =  0.0;
   coords[2] =  0.0;
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
 case  0:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
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
 case  0:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
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
 case  0:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;

 case  5:
   coords[0] =  0.5;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  6:
   coords[0] = -0.5;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  7:
   coords[0] = -0.5;
   coords[1] = -0.5;
   coords[2] =  0.0;
   break;
 case  8:
   coords[0] =  0.5;
   coords[1] = -0.5;
   coords[2] =  0.0;
   break;
 case  9:
   coords[0] =  0.5;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case 10:
   coords[0] =  0.0;
   coords[1] =  0.5;
   coords[2] =  0.5;
   break;
 case 11:
   coords[0] = -0.5;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case 12:
   coords[0] =  0.0;
   coords[1] = -0.5;
   coords[2] =  0.5;
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
 case  0:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case  5:
   coords[0] =  0.5;
   coords[1] = -0.5;
   coords[2] =  0.0;
   break;
 case  6:
   coords[0] = -0.5;
   coords[1] = -0.5;
   coords[2] =  0.0;
   break;
 case  7:
   coords[0] = -0.5;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  8:
   coords[0] =  0.5;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  9:
   coords[0] =  0.5;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case 10:
   coords[0] =  0.0;
   coords[1] = -0.5;
   coords[2] =  0.5;
   break;
 case 11:
   coords[0] = -0.5;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case 12:
   coords[0] =  0.0;
   coords[1] =  0.5;
   coords[2] =  0.5;
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
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] = -0.0;
   coords[2] =  1.0;
   break;
 case  2:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case  5:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] = -1.0;
   coords[1] = -0.0;
   coords[2] =  1.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  5:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case  4:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] = -0.0;
   coords[2] =  1.0;
   break;
 case  2:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  4:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case  5:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;

 case  6:
   coords[0] = -1.0;
   coords[1] =  0.5;
   coords[2] =  0.5;
   break;
 case  7:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case  8:
   coords[0] = -1.0;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  9:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case 10:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case 11:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case 12:
   coords[0] =  1.0;
   coords[1] =  0.5;
   coords[2] =  0.5;
   break;
 case 13:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case 14:
   coords[0] =  1.0;
   coords[1] =  0.5;
   coords[2] =  0.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  2:
   coords[0] = -1.0;
   coords[1] = -0.0;
   coords[2] =  1.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case  5:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case  4:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;

 case  8:
   coords[0] = -1.0;
   coords[1] =  0.5;
   coords[2] =  0.5;
   break;
 case  7:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case  6:
   coords[0] = -1.0;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case 12:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case 14:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case 13:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case 11:
   coords[0] =  1.0;
   coords[1] =  0.5;
   coords[2] =  0.5;
   break;
 case 10:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case  9:
   coords[0] =  1.0;
   coords[1] =  0.5;
   coords[2] =  0.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case  1:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case  3:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case  4:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case  5:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case  6:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;
 case  7:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case  4:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case  7:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case  6:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;
 case  5:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case  1:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case  3:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case  4:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case  5:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case  6:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;
 case  7:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;

 case  8:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case  9:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] = -1.0;
   break;
 case 10:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case 11:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] = -1.0;
   break;
 case 12:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case 13:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case 14:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case 15:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case 16:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case 17:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case 18:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;
 case 19:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
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
 case  0:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case  3:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case  2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case  1:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case  4:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case  7:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case  6:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;
 case  5:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;

 case 11:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case 10:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] = -1.0;
   break;
 case  9:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case  8:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] = -1.0;
   break;
 case 16:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case 19:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case 18:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case 17:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case 15:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case 14:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case 13:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;
 case 12:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
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
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case 1:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case 2:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case 3:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case 4:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case 5:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;
 case 6:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;
 case 7:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case 8:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] = -1.0;
   break;
 case 9:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] = -1.0;
   break;
 case 10:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] = -1.0;
   break;
 case 11:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] = -1.0;
   break;
 case 12:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case 13:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  1.0;
   break;
 case 14:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case 15:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] =  1.0;
   break;
 case 16:
   coords[0] = -1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case 17:
   coords[0] = -1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case 18:
   coords[0] =  1.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case 19:
   coords[0] =  1.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case 20:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] = -1.0;
   break;
 case 21:
   coords[0] = -1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case 22:
   coords[0] =  0.0;
   coords[1] =  1.0;
   coords[2] =  0.0;
   break;
 case 23:
   coords[0] =  1.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
   break;
 case 24:
   coords[0] =  0.0;
   coords[1] = -1.0;
   coords[2] =  0.0;
   break;
 case 25:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  1.0;
   break;
 case 26:
   coords[0] =  0.0;
   coords[1] =  0.0;
   coords[2] =  0.0;
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
