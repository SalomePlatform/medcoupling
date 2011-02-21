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

//Local includes
#include "GaussCoords.hxx"
#include "CellModel.hxx"
using namespace INTERP_KERNEL;

//STL includes
#include <math.h>
#include <algorithm>
#include <sstream>

//#define MYDEBUG

#ifdef MYDEBUG
#include <iostream>
#endif

//Define common part of the code in the MACRO
//---------------------------------------------------------------
#define LOCAL_COORD_MACRO_BEGIN                                         \
  myLocalReferenceCoord.resize( myLocalRefDim*myLocalNbRef );           \
  for( int refId = 0; refId < myLocalNbRef; refId++ )                   \
    {                                                                   \
      double* coords = &myLocalReferenceCoord[ refId*myLocalRefDim ];   \
      switch(refId)                                                     \
        {

//---------------------------------------------------------------
#define LOCAL_COORD_MACRO_END                   \
  }                                             \
}

//---------------------------------------------------------------
#define SHAPE_FUN_MACRO_BEGIN                                           \
  for( int gaussId     = 0 ; gaussId < myNbGauss ; gaussId++ )          \
    {                                                                   \
      double* funValue =  &myFunctionValue[ gaussId * myNbRef ];        \
      const double* gc = &myGaussCoord[ gaussId * GetGaussCoordDim() ];

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
  myGeometry(theGeometry),
  myGaussCoord(theGaussCoord),
  myNbGauss(theNbGauss),
  myReferenceCoord(theReferenceCoord),
  myNbRef(theNbRef) 
{

  //Allocate shape function values
  myFunctionValue.resize( myNbGauss * myNbRef );
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
int GaussInfo::GetGaussCoordDim() const 
{
  if( myNbGauss ) 
    {
      return myGaussCoord.size()/myNbGauss;
    }
  else 
    {
      return 0;
    }
}

/*!
 * Return dimension of the reference coordinates
 */
int GaussInfo::GetReferenceCoordDim() const 
{
  if( myNbRef ) 
    {
      return myReferenceCoord.size()/myNbRef;
    } else 
    {
      return 0;
    }
}

/*!
 * Return type of the cell.
 */
NormalizedCellType GaussInfo::GetCellType() const 
{
  return myGeometry;
}

/*!
 * Return Nb of the gauss points.
 */
int GaussInfo::GetNbGauss() const 
{
  return myNbGauss;
}

/*!
 * Return Nb of the reference coordinates.
 */
int GaussInfo::GetNbRef() const 
{
  return myNbRef;
}

/*!
 * Check coordinates
 */
bool GaussInfo::isSatisfy() 
{

  bool anIsSatisfy = ((myLocalNbRef == myNbRef) && (myLocalRefDim == GetReferenceCoordDim()));
  //Check coordinates
  if(anIsSatisfy)
    {
      for( int refId = 0; refId < myLocalNbRef; refId++ ) 
        {
          double* refCoord = &myReferenceCoord[ refId*myLocalRefDim ];
          double* localRefCoord = &myLocalReferenceCoord[ refId*myLocalRefDim ];
          bool anIsEqual = false;
          for( int dimId = 0; dimId < myLocalRefDim; dimId++ ) 
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

/*!
 * Initialize the internal vectors
 */
void GaussInfo::InitLocalInfo() throw (INTERP_KERNEL::Exception) 
{
  bool aSatify = false;
  const CellModel& cellModel=CellModel::getCellModel(myGeometry);
  switch( myGeometry ) 
    {
    case NORM_SEG2:
      myLocalRefDim = 1;
      myLocalNbRef  = 2;
      Seg2Init();
      aSatify = isSatisfy();
      CHECK_MACRO;
      break;

    case NORM_SEG3:
      myLocalRefDim = 1;
      myLocalNbRef  = 3;
      Seg3Init();
      aSatify = isSatisfy();
      CHECK_MACRO;
      break;

    case NORM_TRI3:
      myLocalRefDim = 2;
      myLocalNbRef  = 3;
      Tria3aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Tria3bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_TRI6:
      myLocalRefDim = 2;
      myLocalNbRef  = 6;
      Tria6aInit();
      aSatify = isSatisfy();
      if(!aSatify)
        {
          Tria6bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_QUAD4:
      myLocalRefDim = 2;
      myLocalNbRef  = 4;
      Quad4aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Quad4bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_QUAD8:
      myLocalRefDim = 2;
      myLocalNbRef  = 8;
      Quad8aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Quad8bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_TETRA4:
      myLocalRefDim = 3;
      myLocalNbRef  = 4;
      Tetra4aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Tetra4bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_TETRA10:
      myLocalRefDim = 3;
      myLocalNbRef  = 10;
      Tetra10aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Tetra10bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_PYRA5:
      myLocalRefDim = 3;
      myLocalNbRef  = 5;
      Pyra5aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Pyra5bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_PYRA13:
      myLocalRefDim = 3;
      myLocalNbRef  = 13;
      Pyra13aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Pyra13bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_PENTA6:
      myLocalRefDim = 3;
      myLocalNbRef  = 6;
      Penta6aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Penta6bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_PENTA15:
      myLocalRefDim = 3;
      myLocalNbRef  = 15;
      Penta15aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Penta15bInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_HEXA8:
      myLocalRefDim = 3;
      myLocalNbRef  = 8;
      Hexa8aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Hexa8aInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    case NORM_HEXA20:
      myLocalRefDim = 3;
      myLocalNbRef  = 20;
      Hexa20aInit();
      aSatify = isSatisfy();

      if(!aSatify)
        {
          Hexa20aInit();
          aSatify = isSatisfy();
          CHECK_MACRO;
        }
      break;

    default:
      throw INTERP_KERNEL::Exception("Bad Cell type !!!");
      break;
    }
}

/**
 * Return shape function value by node id
 */
const double* GaussInfo::GetFunctionValues( const int theGaussId ) const 
{
  return &myFunctionValue[ myNbRef*theGaussId ];
}

/*!
 * Init Segment 2 Reference coordinates ans Shape function.
 */
void GaussInfo::Seg2Init() 
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
void GaussInfo::Seg3Init() 
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
   funValue[0] = 0.5*(1.0 - gc[0])*gc[0];
   funValue[1] = 0.5*(1.0 + gc[0])*gc[0];
   funValue[2] = (1.0 + gc[0])*(1.0 - gc[0]);
   SHAPE_FUN_MACRO_END;
}

/*!
 * Init Triangle Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::Tria3aInit() 
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
void GaussInfo::Tria3bInit() 
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
void GaussInfo::Tria6aInit() 
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
void GaussInfo::Tria6bInit() 
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

/*!
 * Init Quadrangle Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::Quad4aInit() 
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
void GaussInfo::Quad4bInit() 
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


/*!
 * Init Quadratic Quadrangle Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::Quad8aInit() 
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
void GaussInfo::Quad8bInit() 
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

/*!
 * Init Tetrahedron Reference coordinates ans Shape function.
 * Case A.
 */
void GaussInfo::Tetra4aInit() 
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
void GaussInfo::Tetra4bInit() 
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
void GaussInfo::Tetra10aInit() 
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
void GaussInfo::Tetra10bInit() 
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
void GaussInfo::Pyra5aInit() 
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
void GaussInfo::Pyra5bInit() 
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
void GaussInfo::Pyra13aInit() 
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
void GaussInfo::Pyra13bInit() 
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
 case  8:
   coords[0] =  0.5;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  7:
   coords[0] = -0.5;
   coords[1] =  0.5;
   coords[2] =  0.0;
   break;
 case  6:
   coords[0] = -0.5;
   coords[1] = -0.5;
   coords[2] =  0.0;
   break;
 case  5:
   coords[0] =  0.5;
   coords[1] = -0.5;
   coords[2] =  0.0;
   break;
 case  9:
   coords[0] =  0.5;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case 12:
   coords[0] =  0.0;
   coords[1] =  0.5;
   coords[2] =  0.5;
   break;
 case 11:
   coords[0] = -0.5;
   coords[1] =  0.0;
   coords[2] =  0.5;
   break;
 case 10:
   coords[0] =  0.0;
   coords[1] = -0.5;
   coords[2] =  0.5;
   break;
   LOCAL_COORD_MACRO_END;

   SHAPE_FUN_MACRO_BEGIN;
   funValue[0] = 0.5*(-gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] - gc[1] + gc[2] - 1.0)*
     (gc[0] - 0.5)/(1.0 - gc[2]);
   funValue[3] = 0.5*(-gc[0] - gc[1] + gc[2] - 1.0)*(+gc[0] - gc[1] + gc[2] - 1.0)*
     (gc[1] - 0.5)/(1.0 - gc[2]);
   funValue[2] = 0.5*(+gc[0] - gc[1] + gc[2] - 1.0)*(+gc[0] + gc[1] + gc[2] - 1.0)*
     (-gc[0] - 0.5)/(1.0 - gc[2]);
   funValue[1] = 0.5*(+gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] + gc[1] + gc[2] - 1.0)*
     (-gc[1] - 0.5)/(1.0 - gc[2]);

   funValue[4] = 2.0*gc[2]*(gc[2] - 0.5);

   funValue[8] = 0.5*(-gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] - gc[1] + gc[2] - 1.0)*
     (gc[0] - gc[1] + gc[2] - 1.0)/(1.0 - gc[2]);
   funValue[7] = 0.5*(-gc[0] - gc[1] + gc[2] - 1.0)*(gc[0] - gc[1] + gc[2] - 1.0)*
     (gc[0] + gc[1] + gc[2] - 1.0)/(1.0 - gc[2]);
   funValue[6] = 0.5*(gc[0] - gc[1] + gc[2] - 1.0)*(gc[0] + gc[1] + gc[2] - 1.0)*
     (-gc[0] + gc[1] + gc[2] - 1.0)/(1.0 - gc[2]);
   funValue[5] = 0.5*(gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] + gc[1] + gc[2] - 1.0)*
     (-gc[0] - gc[1] + gc[2] - 1.0)/(1.0 - gc[2]);

   funValue[9] = 0.5*gc[2]*(-gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] - gc[1] + gc[2] - 1.0)/
     (1.0 - gc[2]);
   funValue[12] = 0.5*gc[2]*(-gc[0] - gc[1] + gc[2] - 1.0)*(gc[0] - gc[1] + gc[2] - 1.0)/
     (1.0 - gc[2]);
   funValue[11] = 0.5*gc[2]*(gc[0] - gc[1] + gc[2] - 1.0)*(gc[0] + gc[1] + gc[2] - 1.0)/
     (1.0 - gc[2]);
   funValue[10] = 0.5*gc[2]*(gc[0] + gc[1] + gc[2] - 1.0)*(-gc[0] + gc[1] + gc[2] - 1.0)/
     (1.0 - gc[2]);
   SHAPE_FUN_MACRO_END;
}


/*!
 * Init Pentahedron Reference coordinates and Shape function.
 * Case A.
 */
void GaussInfo::Penta6aInit() 
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
void GaussInfo::Penta6bInit() 
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
 * Init Pentahedron Reference coordinates and Shape function.
 * Case A.
 */
void GaussInfo::Penta15aInit() 
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
void GaussInfo::Penta15bInit() 
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
void GaussInfo::Hexa8aInit() 
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
void GaussInfo::Hexa8bInit() 
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
 * Init Qaudratic Hehahedron Reference coordinates and Shape function.
 * Case A.
 */
void GaussInfo::Hexa20aInit() 
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
void GaussInfo::Hexa20bInit()
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
  GaussInfoVector::iterator it = myGaussInfo.begin();
  for( ; it != myGaussInfo.end(); it++ ) 
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
                                int theNbRef) throw (INTERP_KERNEL::Exception) 
{
  GaussInfoVector::iterator it = myGaussInfo.begin();
  for( ; it != myGaussInfo.end(); it++ ) 
    {
      if( (*it)->GetCellType() == theGeometry ) 
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
  info->InitLocalInfo();

  //If info with cell type doesn't exist add it
  if( it == myGaussInfo.end() ) 
    {
      myGaussInfo.push_back(info);

      // If information exists, update it
    } else 
    {
      int index = std::distance(myGaussInfo.begin(),it);
      delete (*it);
      myGaussInfo[index] = info;
    }
}


/*!
 * Calculate gauss points coordinates
 */
double* GaussCoords::CalculateCoords( NormalizedCellType theGeometry, 
                                      const double* theNodeCoords, 
                                      const int theSpaceDim,
                                      const int* theIndex) 
  throw (INTERP_KERNEL::Exception) 
{
  GaussInfoVector::const_iterator it = myGaussInfo.begin();
  GaussInfoVector::const_iterator it_end = myGaussInfo.end();

  //Try to find gauss localization info
  for( ; it != it_end ; it++ ) 
    {
      if( (*it)->GetCellType() == theGeometry ) 
        {
          break;
        }
    }
  if(it == it_end) 
    {
      throw INTERP_KERNEL::Exception("Can't find gauss localization information !!!");
    }
  const GaussInfo* info = (*it);
  int aConn = info->GetNbRef();

  int nbCoords = theSpaceDim * info->GetNbGauss();
  double* aCoords = new double [nbCoords];
  for( int i = 0; i < nbCoords; i++ )
    aCoords[i] = 0.0;

#ifdef MYDEBUG
  std::cout<<"#######################################################"<<std::endl;
  std::cout<<"Cell type : "<<info->GetCellType()<<std::endl;
#endif

  for( int gaussId = 0; gaussId < info->GetNbGauss(); gaussId++ ) 
    {
#ifdef MYDEBUG
      std::cout<<"Gauss ID = "<<gaussId<<std::endl;
#endif

      double* coord = &aCoords[ gaussId*theSpaceDim ];
      const double* function = info->GetFunctionValues(gaussId);
      for ( int connId = 0; connId < aConn ; connId++ ) 
        {
          const double* nodeCoord = theNodeCoords + (theIndex[connId]*theSpaceDim);

#ifdef MYDEBUG
          std::cout<<"Function Value["<<connId<<"] = "<<function[connId]<<std::endl;

          std::cout<<"Node #"<<connId <<" ";
          for( int dimId = 0; dimId < theSpaceDim; dimId++ ) 
            {
              switch(dimId)
                {
                case 0:
                  std::cout<<"( "<<nodeCoord[dimId];
                  break;
                case 1:
                  std::cout<<", "<<nodeCoord[dimId];
                  break;
                case 2:
                  std::cout<<", "<<nodeCoord[dimId]<<" )";
                  break;
                }
            }
          std::cout<<std::endl;
#endif

          for( int dimId = 0; dimId < theSpaceDim; dimId++ ) 
            {
              coord[dimId] += nodeCoord[dimId]*function[connId];
            }
        }
    }
  return aCoords;
}
