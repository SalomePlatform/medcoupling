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
// File      : SauvMedConvertor.cxx
// Created   : Tue Aug 16 14:43:20 2011
// Author    : Edward AGAPOV (eap)
//

#include "SauvMedConvertor.hxx"

#include "CellModel.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileField.hxx"
#include "MEDFileData.hxx"
#include "MEDCouplingFieldDouble.hxx"

#include <iostream>
#include <cassert>
#include <cmath>
#include <queue>
#include <limits>

#include <cstdlib>
#include <cstring>
#include <fcntl.h>

#ifdef WIN32
#include <io.h>
#else
#include <unistd.h>
#endif

#ifdef HAS_XDR
#include <rpc/types.h>
#include <rpc/xdr.h>
#endif

using namespace SauvUtilities;
using namespace MEDCoupling;

namespace
{
  // for ASCII file reader
  const int GIBI_MaxOutputLen = 150; // max length of a line in the sauve file
  const int GIBI_BufferSize   = 16184; // for buffered reading

  using namespace INTERP_KERNEL;

  const size_t MaxMedCellType = NORM_ERROR;
  const size_t NbGibiCellTypes = 47;
  const TCellType GibiTypeToMed[NbGibiCellTypes] =
    {
      /*1 */ NORM_POINT1 ,/*2 */ NORM_SEG2   ,/*3 */ NORM_SEG3   ,/*4 */ NORM_TRI3   ,/*5 */ NORM_ERROR  ,
      /*6 */ NORM_TRI6   ,/*7 */ NORM_ERROR  ,/*8 */ NORM_QUAD4  ,/*9 */ NORM_ERROR  ,/*10*/ NORM_QUAD8  ,
      /*11*/ NORM_ERROR  ,/*12*/ NORM_ERROR  ,/*13*/ NORM_ERROR  ,/*14*/ NORM_HEXA8  ,/*15*/ NORM_HEXA20 ,
      /*16*/ NORM_PENTA6 ,/*17*/ NORM_PENTA15,/*18*/ NORM_ERROR  ,/*19*/ NORM_ERROR  ,/*20*/ NORM_ERROR  ,
      /*21*/ NORM_ERROR  ,/*22*/ NORM_ERROR  ,/*23*/ NORM_TETRA4 ,/*24*/ NORM_TETRA10,/*25*/ NORM_PYRA5  ,
      /*26*/ NORM_PYRA13 ,/*27*/ NORM_ERROR  ,/*28*/ NORM_ERROR  ,/*29*/ NORM_ERROR  ,/*30*/ NORM_ERROR  ,
      /*31*/ NORM_ERROR  ,/*32*/ NORM_ERROR  ,/*33*/ NORM_ERROR  ,/*34*/ NORM_ERROR  ,/*35*/ NORM_ERROR  ,
      /*36*/ NORM_ERROR  ,/*37*/ NORM_ERROR  ,/*38*/ NORM_ERROR  ,/*39*/ NORM_ERROR  ,/*40*/ NORM_ERROR  ,
      /*41*/ NORM_ERROR  ,/*42*/ NORM_ERROR  ,/*43*/ NORM_ERROR  ,/*44*/ NORM_ERROR  ,/*45*/ NORM_ERROR  ,
      /*46*/ NORM_ERROR  ,/*47*/ NORM_ERROR
    };

  //================================================================================
  /*!
   * \brief Return dimension of a group
   */
  //================================================================================

  unsigned getDim( const Group* grp )
  {
    return SauvUtilities::getDimension( grp->_groups.empty() ? grp->_cellType : grp->_groups[0]->_cellType );
  }

  //================================================================================
  /*!
   * \brief Converts connectivity of quadratic elements
   */
  //================================================================================

  inline void ConvertQuadratic( const INTERP_KERNEL::NormalizedCellType type,
                                const Cell &                            aCell )
  {
    if ( const int * conn = getGibi2MedQuadraticInterlace( type ))
      {
        Cell* ma = const_cast<Cell*>(&aCell);
        std::vector< Node* > new_nodes( ma->_nodes.size() );
        for (std:: size_t i = 0; i < new_nodes.size(); ++i )
          new_nodes[ i ] = ma->_nodes[ conn[ i ]];
        ma->_nodes.swap( new_nodes );
      }
  }

  //================================================================================
  /*!
   * \brief Returns a vector of pairs of node indices to inverse a med volume element
   */
  //================================================================================

  void getReverseVector (const INTERP_KERNEL::NormalizedCellType type,
                         std::vector<std::pair<int,int> > &                swapVec )
  {
    swapVec.clear();

    switch ( type )
      {
      case NORM_TETRA4:
        swapVec.resize(1);
        swapVec[0] = std::make_pair( 1, 2 );
        break;
      case NORM_PYRA5:
        swapVec.resize(1);
        swapVec[0] = std::make_pair( 1, 3 );
        break;
      case NORM_PENTA6:
        swapVec.resize(2);
        swapVec[0] = std::make_pair( 1, 2 );
        swapVec[1] = std::make_pair( 4, 5 );
        break;
      case NORM_HEXA8:
        swapVec.resize(2);
        swapVec[0] = std::make_pair( 1, 3 );
        swapVec[1] = std::make_pair( 5, 7 );
        break;
      case NORM_TETRA10:
        swapVec.resize(3);
        swapVec[0] = std::make_pair( 1, 2 );
        swapVec[1] = std::make_pair( 4, 6 );
        swapVec[2] = std::make_pair( 8, 9 );
        break;
      case NORM_PYRA13:
        swapVec.resize(4);
        swapVec[0] = std::make_pair( 1, 3 );
        swapVec[1] = std::make_pair( 5, 8 );
        swapVec[2] = std::make_pair( 6, 7 );
        swapVec[3] = std::make_pair( 10, 12 );
        break;
      case NORM_PENTA15:
        swapVec.resize(4);
        swapVec[0] = std::make_pair( 1, 2 );
        swapVec[1] = std::make_pair( 4, 5 );
        swapVec[2] = std::make_pair( 6, 8 );
        swapVec[3] = std::make_pair( 9, 11 );
        break;
      case NORM_HEXA20:
        swapVec.resize(7);
        swapVec[0] = std::make_pair( 1, 3 );
        swapVec[1] = std::make_pair( 5, 7 );
        swapVec[2] = std::make_pair( 8, 11 );
        swapVec[3] = std::make_pair( 9, 10 );
        swapVec[4] = std::make_pair( 12, 15 );
        swapVec[5] = std::make_pair( 13, 14 );
        swapVec[6] = std::make_pair( 17, 19 );
        break;
        //   case NORM_SEG3: no need to reverse edges
        //     swapVec.resize(1);
        //     swapVec[0] = std::make_pair( 1, 2 );
        //     break;
      case NORM_TRI6:
        swapVec.resize(2);
        swapVec[0] = std::make_pair( 1, 2 );
        swapVec[1] = std::make_pair( 3, 5 );
        break;
      case NORM_QUAD8:
        swapVec.resize(3);
        swapVec[0] = std::make_pair( 1, 3 );
        swapVec[1] = std::make_pair( 4, 7 );
        swapVec[2] = std::make_pair( 5, 6 );
        break;
      default:;
      }
  }

  //================================================================================
  /*!
   * \brief Inverses element orientation using vector of indices to swap
   */
  //================================================================================

  inline void reverse(const Cell & aCell, const std::vector<std::pair<int,int> > & swapVec )
  {
    Cell* ma = const_cast<Cell*>(&aCell);
    for ( unsigned i = 0; i < swapVec.size(); ++i )
      std::swap( ma->_nodes[ swapVec[i].first ],
                 ma->_nodes[ swapVec[i].second ]);
    if ( swapVec.empty() )
      ma->_reverse = true;
    else
      ma->_reverse = false;
  }
  //================================================================================
  /*!
   * \brief Comparator of cells by number used for ordering cells within a med group
   */
  struct TCellByIDCompare
  {
    bool operator () (const Cell* i1, const Cell* i2) const
    {
      return i1->_number < i2->_number;
    }
  };
  typedef std::map< const Cell*, unsigned, TCellByIDCompare > TCellToOrderMap;

  //================================================================================
  /*!
   * \brief Fill Group::_relocTable if necessary
   */
  //================================================================================

  void setRelocationTable( Group* grp, TCellToOrderMap& cell2order )
  {
    if ( !grp->_isProfile ) return;

    // check if relocation table is necessary
    bool isRelocated = false;
    unsigned newOrder = 0;
    TCellToOrderMap::iterator c2oIt = cell2order.begin(), c2oEnd = cell2order.end();
    for ( ; !isRelocated && c2oIt != c2oEnd; ++c2oIt, ++newOrder )
      isRelocated = ( c2oIt->second != newOrder );

    if ( isRelocated )
      {
        grp->_relocTable.resize( cell2order.size() );
        for ( newOrder = 0, c2oIt = cell2order.begin(); c2oIt != c2oEnd; ++c2oIt, ++newOrder )
          grp->_relocTable[ c2oIt->second ] = newOrder;
      }
  }
}

namespace // define default GAUSS points
{
  typedef std::vector<double> TDoubleVector;
  typedef double*             TCoordSlice;
  typedef int                 TInt;
  //---------------------------------------------------------------
  //! Shape function definitions
  //---------------------------------------------------------------
  struct TShapeFun
  {
    TInt myDim;
    TInt myNbRef;
    TDoubleVector myRefCoord;

    TShapeFun(TInt theDim = 0, TInt theNbRef = 0)
      :myDim(theDim),myNbRef(theNbRef),myRefCoord(theNbRef*theDim) {}

    TInt GetNbRef() const { return myNbRef; }

    TCoordSlice GetCoord(TInt theRefId) { return &myRefCoord[0] + theRefId * myDim; }
  };
  //---------------------------------------------------------------
  /*!
   * \brief Description of family of integration points
   */
  //---------------------------------------------------------------
  struct TGaussDef
  {
    int           myType;      //!< element geometry (EGeometrieElement or med_geometrie_element)
    TDoubleVector myRefCoords; //!< description of reference points
    TDoubleVector myCoords;    //!< coordinates of Gauss points
    TDoubleVector myWeights;   //!< weights, len(weights)==<nb of gauss points>

    /*!
     * \brief Creates definition of gauss points family
     *  \param geomType - element geometry (EGeometrieElement or med_geometrie_element)
     *  \param nbPoints - nb gauss point
     *  \param variant - [1-3] to choose the variant of definition
     * 
     * Throws in case of invalid parameters
     * variant == 1 refers to "Fonctions de forme et points d'integration 
     *              des elements finis" v7.4 by J. PELLET, X. DESROCHES, 15/09/05
     * variant == 2 refers to the same doc v6.4 by J.P. LEFEBVRE, X. DESROCHES, 03/07/03
     * variant == 3 refers to the same doc v6.4, second variant for 2D elements
     */
    TGaussDef(const int geomType, const int nbPoints, const int variant=1);

    int dim() const { return SauvUtilities::getDimension( NormalizedCellType( myType )); }
    int nbPoints() const { return myWeights.capacity(); }

  private:
    void add(const double x, const double weight);
    void add(const double x, const double y, const double weight);
    void add(const double x, const double y, const double z, const double weight);
    void setRefCoords(const TShapeFun& aShapeFun) { myRefCoords = aShapeFun.myRefCoord; }
  };
  struct TSeg2a: TShapeFun {
    TSeg2a();
  };
  struct TSeg3a: TShapeFun {
    TSeg3a();
  };
  struct TTria3a: TShapeFun {
    TTria3a();
  };
  struct TTria6a: TShapeFun {
    TTria6a();
  };
  struct TTria3b: TShapeFun {
    TTria3b();
  };
  struct TTria6b: TShapeFun {
    TTria6b();
  };
  struct TQuad4a: TShapeFun {
    TQuad4a();
  };
  struct TQuad8a: TShapeFun {
    TQuad8a();
  };
  struct TQuad4b: TShapeFun {
    TQuad4b();
  };
  struct TQuad8b: TShapeFun {
    TQuad8b();
  };
  struct TTetra4a: TShapeFun {
    TTetra4a();
  };
  struct TTetra10a: TShapeFun {
    TTetra10a();
  };
  struct TTetra4b: TShapeFun {
    TTetra4b();
  };
  struct TTetra10b: TShapeFun {
    TTetra10b();
  };
  struct THexa8a: TShapeFun {
    THexa8a();
  };
  struct THexa20a: TShapeFun {
    THexa20a(TInt theDim = 3, TInt theNbRef = 20);
  };
  struct THexa27a: THexa20a {
    THexa27a();
  };
  struct THexa8b: TShapeFun {
    THexa8b();
  };
  struct THexa20b: TShapeFun {
    THexa20b(TInt theDim = 3, TInt theNbRef = 20);
  };
  struct TPenta6a: TShapeFun {
    TPenta6a();
  };
  struct TPenta6b: TShapeFun {
    TPenta6b();
  };
  struct TPenta15a: TShapeFun {
    TPenta15a();
  };
  struct TPenta15b: TShapeFun {
    TPenta15b();
  };
  struct TPyra5a: TShapeFun {
    TPyra5a();
  };
  struct TPyra5b: TShapeFun {
    TPyra5b();
  };
  struct TPyra13a: TShapeFun {
    TPyra13a();
  };
  struct TPyra13b: TShapeFun {
    TPyra13b();
  };

  void TGaussDef::add(const double x, const double weight)
  {
    if ( dim() != 1 )
      THROW_IK_EXCEPTION("TGaussDef: dim() != 1");
    if ( myWeights.capacity() == myWeights.size() )
      THROW_IK_EXCEPTION("TGaussDef: Extra gauss point");
    myCoords.push_back( x );
    myWeights.push_back( weight );
  }
  void TGaussDef::add(const double x, const double y, const double weight)
  {
    if ( dim() != 2 )
      THROW_IK_EXCEPTION("TGaussDef: dim() != 2");
    if ( myWeights.capacity() == myWeights.size() )
      THROW_IK_EXCEPTION("TGaussDef: Extra gauss point");
    myCoords.push_back( x );
    myCoords.push_back( y );
    myWeights.push_back( weight );
  }
  void TGaussDef::add(const double x, const double y, const double z, const double weight)
  {
    if ( dim() != 3 )
      THROW_IK_EXCEPTION("TGaussDef: dim() != 3");
    if ( myWeights.capacity() == myWeights.size() )
      THROW_IK_EXCEPTION("TGaussDef: Extra gauss point");
    myCoords.push_back( x );
    myCoords.push_back( y );
    myCoords.push_back( z );
    myWeights.push_back( weight );
  }

  //---------------------------------------------------------------
  TSeg2a::TSeg2a():TShapeFun(1,2)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0; break;
      case  1: aCoord[0] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TSeg3a::TSeg3a():TShapeFun(1,3)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0; break;
      case  1: aCoord[0] =  1.0; break;
      case  2: aCoord[0] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TTria3a::TTria3a():
    TShapeFun(2,3)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] =  1.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] = -1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TTria6a::TTria6a():TShapeFun(2,6)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] =  1.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] = -1.0; break;

      case  3: aCoord[0] = -1.0;  aCoord[1] =  1.0; break;
      case  4: aCoord[0] =  0.0;  aCoord[1] = -1.0; break;
      case  5: aCoord[0] =  0.0;  aCoord[1] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TTria3b::TTria3b():
    TShapeFun(2,3)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] =  0.0;  aCoord[1] =  0.0; break;
      case  1: aCoord[0] =  1.0;  aCoord[1] =  0.0; break;
      case  2: aCoord[0] =  0.0;  aCoord[1] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TTria6b::TTria6b():
    TShapeFun(2,6)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] =  0.0;  aCoord[1] =  0.0; break;
      case  1: aCoord[0] =  1.0;  aCoord[1] =  0.0; break;
      case  2: aCoord[0] =  0.0;  aCoord[1] =  1.0; break;

      case  3: aCoord[0] =  0.5;  aCoord[1] =  0.0; break;
      case  4: aCoord[0] =  0.5;  aCoord[1] =  0.5; break;
      case  5: aCoord[0] =  0.0;  aCoord[1] =  0.5; break;
      }
    }
  }
  //---------------------------------------------------------------
  TQuad4a::TQuad4a():
    TShapeFun(2,4)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] =  1.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] = -1.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TQuad8a::TQuad8a():
    TShapeFun(2,8)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] =  1.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] = -1.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  1.0; break;

      case  4: aCoord[0] = -1.0;  aCoord[1] =  0.0; break;
      case  5: aCoord[0] =  0.0;  aCoord[1] = -1.0; break;
      case  6: aCoord[0] =  1.0;  aCoord[1] =  0.0; break;
      case  7: aCoord[0] =  0.0;  aCoord[1] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TQuad4b::TQuad4b():
    TShapeFun(2,4)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] = -1.0; break;
      case  1: aCoord[0] =  1.0;  aCoord[1] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] =  1.0; break;
      case  3: aCoord[0] = -1.0;  aCoord[1] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TQuad8b::TQuad8b():
    TShapeFun(2,8)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] = -1.0; break;
      case  1: aCoord[0] =  1.0;  aCoord[1] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] =  1.0; break;
      case  3: aCoord[0] = -1.0;  aCoord[1] =  1.0; break;

      case  4: aCoord[0] =  0.0;  aCoord[1] = -1.0; break;
      case  5: aCoord[0] =  1.0;  aCoord[1] =  0.0; break;
      case  6: aCoord[0] =  0.0;  aCoord[1] =  1.0; break;
      case  7: aCoord[0] = -1.0;  aCoord[1] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TTetra4a::TTetra4a():
    TShapeFun(3,4)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  1: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case  2: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TTetra10a::TTetra10a():
    TShapeFun(3,10)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  1: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case  2: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;

      case  4: aCoord[0] =  0.0;  aCoord[1] =  0.5;  aCoord[2] =  0.5; break;
      case  5: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case  6: aCoord[0] =  0.0;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;

      case  7: aCoord[0] =  0.5;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      case  8: aCoord[0] =  0.5;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case  9: aCoord[0] =  0.5;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TTetra4b::TTetra4b():
    TShapeFun(3,4)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  2: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case  1: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TTetra10b::TTetra10b():
    TShapeFun(3,10)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  2: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case  1: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;

      case  6: aCoord[0] =  0.0;  aCoord[1] =  0.5;  aCoord[2] =  0.5; break;
      case  5: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case  4: aCoord[0] =  0.0;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;

      case  7: aCoord[0] =  0.5;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      case  9: aCoord[0] =  0.5;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case  8: aCoord[0] =  0.5;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  THexa8a::THexa8a():
    TShapeFun(3,8)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case  1: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case  3: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case  4: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case  5: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case  6: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;
      case  7: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  THexa20a::THexa20a(TInt theDim, TInt theNbRef):
    TShapeFun(theDim,theNbRef)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case  1: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case  3: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case  4: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case  5: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case  6: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;
      case  7: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;

      case  8: aCoord[0] =  0.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case  9: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] = -1.0; break;
      case 10: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case 11: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] = -1.0; break;
      case 12: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] =  0.0; break;
      case 13: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] =  0.0; break;
      case 14: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case 15: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case 16: aCoord[0] =  0.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case 17: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case 18: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;
      case 19: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  THexa27a::THexa27a():
    THexa20a(3,27)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case 20: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] = -1.0; break;
      case 21: aCoord[0] =  0.0;  aCoord[1] = -1.0;  aCoord[2] =  0.0; break;
      case 22: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case 23: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case 24: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case 25: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case 26: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  THexa8b::THexa8b():
    TShapeFun(3,8)
  {
    TInt aNbRef = GetNbRef();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case  4: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case  7: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case  6: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;
      case  5: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  THexa20b::THexa20b(TInt theDim, TInt theNbRef):
    TShapeFun(theDim,theNbRef)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case  2: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case  4: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case  7: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case  6: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;
      case  5: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;

      case 11: aCoord[0] =  0.0;  aCoord[1] = -1.0;  aCoord[2] = -1.0; break;
      case 10: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] = -1.0; break;
      case  9: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] = -1.0; break;
      case  8: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] = -1.0; break;
      case 16: aCoord[0] = -1.0;  aCoord[1] = -1.0;  aCoord[2] =  0.0; break;
      case 19: aCoord[0] =  1.0;  aCoord[1] = -1.0;  aCoord[2] =  0.0; break;
      case 18: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case 17: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case 15: aCoord[0] =  0.0;  aCoord[1] = -1.0;  aCoord[2] =  1.0; break;
      case 14: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case 13: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  1.0; break;
      case 12: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TPenta6a::TPenta6a():
    TShapeFun(3,6)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] = -0.0;  aCoord[2] =  1.0; break;
      case  2: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  4: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case  5: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TPenta6b::TPenta6b():
    TShapeFun(3,6)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  2: aCoord[0] = -1.0;  aCoord[1] = -0.0;  aCoord[2] =  1.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  5: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case  4: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TPenta15a::TPenta15a():
    TShapeFun(3,15)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] = -0.0;  aCoord[2] =  1.0; break;
      case  2: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  4: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case  5: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;

      case  6: aCoord[0] = -1.0;  aCoord[1] =  0.5;  aCoord[2] =  0.5; break;
      case  7: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case  8: aCoord[0] = -1.0;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      case  9: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case 10: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case 11: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case 12: aCoord[0] =  1.0;  aCoord[1] =  0.5;  aCoord[2] =  0.5; break;
      case 13: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case 14: aCoord[0] =  1.0;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TPenta15b::TPenta15b():
    TShapeFun(3,15)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] = -1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  2: aCoord[0] = -1.0;  aCoord[1] = -0.0;  aCoord[2] =  1.0; break;
      case  1: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  1.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  5: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case  4: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;

      case  8: aCoord[0] = -1.0;  aCoord[1] =  0.5;  aCoord[2] =  0.5; break;
      case  7: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case  6: aCoord[0] = -1.0;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      case 12: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case 14: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      case 13: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case 11: aCoord[0] =  1.0;  aCoord[1] =  0.5;  aCoord[2] =  0.5; break;
      case 10: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case  9: aCoord[0] =  1.0;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TPyra5a::TPyra5a():
    TShapeFun(3,5)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  1: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  2: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  0.0;  aCoord[1] = -1.0;  aCoord[2] =  0.0; break;
      case  4: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TPyra5b::TPyra5b():
    TShapeFun(3,5)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){        
      case  0: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  2: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  1: aCoord[0] =  0.0;  aCoord[1] = -1.0;  aCoord[2] =  0.0; break;
      case  4: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;
      }
    }
  }
  //---------------------------------------------------------------
  TPyra13a::TPyra13a():
    TShapeFun(3,13)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  1: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  2: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  0.0;  aCoord[1] = -1.0;  aCoord[2] =  0.0; break;
      case  4: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;

      case  5: aCoord[0] =  0.5;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      case  6: aCoord[0] = -0.5;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      case  7: aCoord[0] = -0.5;  aCoord[1] = -0.5;  aCoord[2] =  0.0; break;
      case  8: aCoord[0] =  0.5;  aCoord[1] = -0.5;  aCoord[2] =  0.0; break;
      case  9: aCoord[0] =  0.5;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case 10: aCoord[0] =  0.0;  aCoord[1] =  0.5;  aCoord[2] =  0.5; break;
      case 11: aCoord[0] = -0.5;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case 12: aCoord[0] =  0.0;  aCoord[1] = -0.5;  aCoord[2] =  0.5; break;
      }
    }
  }
  //---------------------------------------------------------------
  TPyra13b::TPyra13b():
    TShapeFun(3,13)
  {
    TInt aNbRef = myRefCoord.size();
    for(TInt aRefId = 0; aRefId < aNbRef; aRefId++){
      TCoordSlice aCoord = GetCoord(aRefId);
      switch(aRefId){
      case  0: aCoord[0] =  1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  3: aCoord[0] =  0.0;  aCoord[1] =  1.0;  aCoord[2] =  0.0; break;
      case  2: aCoord[0] = -1.0;  aCoord[1] =  0.0;  aCoord[2] =  0.0; break;
      case  1: aCoord[0] =  0.0;  aCoord[1] = -1.0;  aCoord[2] =  0.0; break;
      case  4: aCoord[0] =  0.0;  aCoord[1] =  0.0;  aCoord[2] =  1.0; break;

      case  8: aCoord[0] =  0.5;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      case  7: aCoord[0] = -0.5;  aCoord[1] =  0.5;  aCoord[2] =  0.0; break;
      case  6: aCoord[0] = -0.5;  aCoord[1] = -0.5;  aCoord[2] =  0.0; break;
      case  5: aCoord[0] =  0.5;  aCoord[1] = -0.5;  aCoord[2] =  0.0; break;
      case  9: aCoord[0] =  0.5;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case 12: aCoord[0] =  0.0;  aCoord[1] =  0.5;  aCoord[2] =  0.5; break;
      case 11: aCoord[0] = -0.5;  aCoord[1] =  0.0;  aCoord[2] =  0.5; break;
      case 10: aCoord[0] =  0.0;  aCoord[1] = -0.5;  aCoord[2] =  0.5; break;
      }
    }
  }
  /*!
   * \brief Fill definition of gauss points family
   */

  TGaussDef::TGaussDef(const int geom, const int nbGauss, const int variant)
  {
    myType = geom;
    myCoords .reserve( nbGauss * dim() );
    myWeights.reserve( nbGauss );

    switch ( geom ) {

    case NORM_SEG2:
    case NORM_SEG3:
      if (geom == NORM_SEG2) setRefCoords( TSeg2a() );
      else                   setRefCoords( TSeg3a() );
      switch ( nbGauss ) {
      case 1: {
        add( 0.0, 2.0 ); break;
      }
      case 2: {
        const double a = 0.577350269189626;
        add(  a,  1.0 );
        add( -a,  1.0 ); break;
      }
      case 3: {
        const double a = 0.774596669241;
        const double P1 = 1./1.8;
        const double P2 = 1./1.125;
        add( -a,  P1 );
        add(  0,  P2 ); 
        add(  a,  P1 ); break;
      }
      case 4: {
        const double a  = 0.339981043584856, b  = 0.861136311594053;
        const double P1 = 0.652145154862546, P2 = 0.347854845137454 ;
        add(  a,  P1 );
        add( -a,  P1 );
        add(  b,  P2 ); 
        add( -b,  P2 ); break;
      }
      default:
        THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for SEG"<<nbGauss);
      }
      break;

    case NORM_TRI3:
    case NORM_TRI6:
      if ( variant == 1 ) {
        if (geom == NORM_TRI3) setRefCoords( TTria3b() );
        else                   setRefCoords( TTria6b() );
        switch ( nbGauss ) {
        case 1: { // FPG1
          add( 1/3., 1/3., 1/2. ); break;
        }
        case 3: { // FPG3
          // what about COT3 ???
          add( 1/6., 1/6., 1/6. );
          add( 2/3., 1/6., 1/6. );
          add( 1/6., 2/3., 1/6. ); break;
        }
        case 4: { // FPG4
          add( 1/5., 1/5.,  25/(24*4.) );
          add( 3/5., 1/5.,  25/(24*4.) );
          add( 1/5., 3/5.,  25/(24*4.) );
          add( 1/3., 1/3., -27/(24*4.) ); break;
        }
        case 6: { // FPG6
          const double P1 = 0.11169079483905, P2 = 0.0549758718227661;
          const double a  = 0.445948490915965, b = 0.091576213509771;
          add(     b,     b, P2 ); 
          add( 1-2*b,     b, P2 );
          add(     b, 1-2*b, P2 );
          add(     a, 1-2*a, P1 );
          add(     a,     a, P1 ); 
          add( 1-2*a,     a, P1 ); break;
        }
        case 7: { // FPG7
          const double A  = 0.470142064105115;
          const double B  = 0.101286507323456;
          const double P1 = 0.066197076394253;
          const double P2 = 0.062969590272413;
          add(  1/3.,  1/3., 9/80. ); 
          add(     A,     A, P1 ); 
          add( 1-2*A,     A, P1 );
          add(     A, 1-2*A, P1 );
          add(     B,     B, P2 ); 
          add( 1-2*B,     B, P2 );
          add(     B, 1-2*B, P2 ); break;
        }
        case 12: { // FPG12
          const double A  = 0.063089014491502;
          const double B  = 0.249286745170910;
          const double C  = 0.310352451033785;
          const double D  = 0.053145049844816;
          const double P1 = 0.025422453185103;
          const double P2 = 0.058393137863189;
          const double P3 = 0.041425537809187;
          add(     A,     A, P1 ); 
          add( 1-2*A,     A, P1 );
          add(     A, 1-2*A, P1 );
          add(     B,     B, P2 ); 
          add( 1-2*B,     B, P2 );
          add(     B, 1-2*B, P2 );
          add(     C,     D, P3 );
          add(     D,     C, P3 );
          add( 1-C-D,     C, P3 );
          add( 1-C-D,     D, P3 );
          add(     C, 1-C-D, P3 );
          add(     D, 1-C-D, P3 ); break;
        }
        default:
          THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for TRIA, variant 1: "
                             <<nbGauss);
        }
      }
      else if ( variant == 2 ) {
        if (geom == NORM_TRI3) setRefCoords( TTria3a() );
        else                   setRefCoords( TTria6a() );
        switch ( nbGauss ) {
        case 1: {
          add( -1/3., -1/3., 2. ); break;
        }
        case 3: {
          add( -2/3.,  1/3., 2/3. );
          add( -2/3., -2/3., 2/3. );
          add(  1/3., -2/3., 2/3. ); break;
        }
        case 6: {
          const double P1 = 0.11169079483905, P2 = 0.0549758718227661;
          const double A  = 0.445948490915965, B = 0.091576213509771;
          add( 2*B-1, 1-4*B, 4*P2 ); 
          add( 2*B-1, 2*B-1, 4*P2 );
          add( 1-4*B, 2*B-1, 4*P2 );
          add( 1-4*A, 2*A-1, 4*P1 );
          add( 2*A-1, 1-4*A, 4*P1 ); 
          add( 2*A-1, 2*A-1, 4*P1 ); break;
        }
        default:
          THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for TRIA, variant 2: "
                             <<nbGauss);
        }
      }
      else if ( variant == 3 ) {
        if (geom == NORM_TRI3) setRefCoords( TTria3b() );
        else                   setRefCoords( TTria6b() );
        switch ( nbGauss ) {
        case 4: {
          add( 1/3., 1/3., -27/96 );
          add( 0.2 , 0.2 ,  25/96 );
          add( 0.6 , 0.2 ,  25/96 );
          add( 0.2 , 0.6 ,  25/96 ); break;
        }
        default:
          THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for TRIA, variant 3: "
                             <<nbGauss);
        }
      }
      break;

    case NORM_QUAD4:
    case NORM_QUAD8:
      if ( variant == 1 ) {
        if (geom == NORM_QUAD4) setRefCoords( TQuad4b() );
        else                    setRefCoords( TQuad8b() );
        switch ( nbGauss ) {
        case 1: { // FPG1
          add(  0,  0,  4 ); break;
        }
        case 4: { // FPG4
          const double a = 1/sqrt(3.);
          add( -a, -a,  1 );
          add(  a, -a,  1 );
          add(  a,  a,  1 );
          add( -a,  a,  1 ); break;
        }
        case 5: { // out from the 3 specs
          const double a = 0.774596669241483;
          add( -a, -a,  0.5 );
          add(  a, -a,  0.5 );
          add(  a,  a,  0.5 );
          add( -a,  a,  0.5 );
          add(  0,  0,  2.0 ); break;
        }
        case 9: { // FPG9
          const double a = 0.774596669241483;
          add( -a, -a,  25/81. );
          add(  a, -a,  25/81. );
          add(  a,  a,  25/81. );
          add( -a,  a,  25/81. );
          add( 0., -a,  40/81. );
          add(  a, 0.,  40/81. );
          add( 0.,  a,  40/81. );
          add( -a, 0.,  40/81. );
          add( 0., 0.,  64/81. ); break;
        }
        default:
          THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for QUAD, variant 1: "
                             <<nbGauss);
        }
      }
      else if ( variant == 2 ) {
        if (geom == NORM_QUAD4) setRefCoords( TQuad4a() );
        else                    setRefCoords( TQuad8a() );
        switch ( nbGauss ) {
        case 4: {
          const double a = 1/sqrt(3.);
          add( -a,  a,  1 );
          add( -a, -a,  1 );
          add(  a, -a,  1 );
          add(  a,  a,  1 ); break;
        }
        case 9: {
          const double a = 0.774596669241483;
          add( -a,  a,  25/81. );
          add( -a, -a,  25/81. );
          add(  a, -a,  25/81. );
          add(  a,  a,  25/81. );
          add( -a, 0.,  40/81. );
          add( 0., -a,  40/81. );
          add(  a, 0.,  40/81. );
          add( 0.,  a,  40/81. );
          add( 0., 0.,  64/81. ); break;
        }
        default:
          THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for QUAD, variant 1: "
                             <<nbGauss);
        }
      }
      else if ( variant == 3 ) {
        if (geom == NORM_QUAD4) setRefCoords( TQuad4b() );
        else                    setRefCoords( TQuad8b() );
        switch ( nbGauss ) {
        case 4: {
          const double a = 3/sqrt(3.);
          add( -a, -a,  1 );
          add( -a,  a,  1 );
          add(  a, -a,  1 );
          add(  a,  a,  1 ); break;
        }
        case 9: {
          const double a = sqrt(3/5.), c1 = 5/9., c2 = 8/9.;
          const double c12 = c1*c2, c22 = c2*c2, c1c2 = c1*c2;
          add( -a, -a,  c12  );
          add( -a, 0.,  c1c2 );
          add( -a,  a,  c12  );
          add( 0., -a,  c1c2 );
          add( 0., 0.,  c22  );
          add( 0.,  a,  c1c2 );
          add(  a, -a,  c12  );
          add(  a, 0.,  c1c2 );
          add(  a,  a,  c12  ); break;
        }
        default:
          THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for QUAD, variant 3: "
                             <<nbGauss);
        }
      }
      break;

    case NORM_TETRA4:
    case NORM_TETRA10:
      if (geom == NORM_TETRA4) setRefCoords( TTetra4a() );
      else                 setRefCoords( TTetra10a() );
      switch ( nbGauss ) {
      case 4: { // FPG4
        const double a = (5 - sqrt(5.))/20., b = (5 + 3*sqrt(5.))/20.;
        add(  a,  a,  a,  1/24. );
        add(  a,  a,  b,  1/24. );
        add(  a,  b,  a,  1/24. );
        add(  b,  a,  a,  1/24. ); break;
      }
      case 5: { // FPG5
        const double a = 0.25, b = 1/6., c = 0.5;
        add(  a,  a,  a, -2/15. );
        add(  b,  b,  b,  3/40. );
        add(  b,  b,  c,  3/40. );
        add(  b,  c,  b,  3/40. );
        add(  c,  b,  b,  3/40. ); break;
      }
      case 15: { // FPG15
        const double a = 0.25;
        const double b1 = (7 + sqrt(15.))/34., c1 = (13 + 3*sqrt(15.))/34., d = (5 - sqrt(15.))/20.;
        const double b2 = (7 - sqrt(15.))/34., c2 = (13 - 3*sqrt(15.))/34., e = (5 + sqrt(15.))/20.;
        const double P1 = (2665 - 14*sqrt(15.))/226800.;
        const double P2 = (2665 + 14*sqrt(15.))/226800.;
        add(  a,  a,  a,  8/405.);//_____
        add( b1, b1, b1,  P1    );
        add( b1, b1, c1,  P1    );
        add( b1, c1, b1,  P1    );
        add( c1, b1, b1,  P1    );//_____
        add( b2, b2, b2,  P2    );
        add( b2, b2, c2,  P2    );
        add( b2, c2, b2,  P2    );
        add( c2, b2, b2,  P2    );//_____
        add(  d,  d,  e,  5/567.);
        add(  d,  e,  d,  5/567.);
        add(  e,  d,  d,  5/567.);
        add(  d,  e,  e,  5/567.);
        add(  e,  d,  e,  5/567.);
        add(  e,  e,  d,  5/567.);
        break;
      }
      default:
        THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for TETRA: "<<nbGauss);
      }
      break;

    case NORM_PYRA5:
    case NORM_PYRA13:
      if (geom == NORM_PYRA5) setRefCoords( TPyra5a() );
      else                setRefCoords( TPyra13a() );
      switch ( nbGauss ) {
      case 5: { // FPG5
        const double h1 = 0.1531754163448146;
        const double h2 = 0.6372983346207416;
        add(  .5,  0.,  h1,  2/15. );
        add(  0.,  .5,  h1,  2/15. );
        add( -.5,  0.,  h1,  2/15. );
        add(  0., -.5,  h1,  2/15. );
        add(  0.,  0.,  h2,  2/15. ); break;
      }
      case 6: { // FPG6
        const double p1 = 0.1024890634400000 ;
        const double p2 = 0.1100000000000000 ;
        const double p3 = 0.1467104129066667 ;
        const double a  = 0.5702963741068025 ;
        const double h1 = 0.1666666666666666 ;
        const double h2 = 0.08063183038464675;
        const double h3 = 0.6098484849057127 ;
        add(  a, 0.,  h1,  p1 );
        add( 0.,  a,  h1,  p1 );
        add( -a, 0.,  h1,  p1 );
        add( 0., -a,  h1,  p1 );
        add( 0., 0.,  h2,  p2 );
        add( 0., 0.,  h3,  p3 ); break;
      }
      case 27: { // FPG27
        const double a1  = 0.788073483; 
        const double b6  = 0.499369002; 
        const double b1  = 0.848418011; 
        const double c8  = 0.478508449; 
        const double c1  = 0.652816472; 
        const double d12 = 0.032303742; 
        const double d1  = 1.106412899;
        double z = 1/2., fz = b1/2*(1 - z);
        add(  0.,  0.,   z,  a1 ); // 1
        add(  fz,  fz,   z,  b6 ); // 2
        add( -fz,  fz,   z,  b6 ); // 3
        add( -fz, -fz,   z,  b6 ); // 4
        add(  fz, -fz,   z,  b6 ); // 5
        z = (1 - b1)/2.;
        add(  0.,  0.,   z,  b6 ); // 6
        z = (1 + b1)/2.;
        add(  0.,  0.,   z,  b6 ); // 7
        z = (1 - c1)/2.; fz = c1*(1 - z);
        add(  fz,  0.,   z,  c8 ); // 8
        add(  0.,  fz,   z,  c8 ); // 9
        add( -fz,  0.,   z,  c8 ); // 10
        add(  0., -fz,   z,  c8 ); // 11
        z = (1 + c1)/2.; fz = c1*(1 - z);
        add(  fz,  0.,   z,  c8 ); // 12
        add(  0.,  fz,   z,  c8 ); // 13
        add( -fz,  0.,   z,  c8 ); // 14
        add(  0., -fz,   z,  c8 ); // 15
        z = (1 - d1)/2., fz = d1/2*(1 - z);
        add(  fz,  fz,   z,  d12); // 16
        add( -fz,  fz,   z,  d12); // 17
        add( -fz, -fz,   z,  d12); // 18
        add(  fz, -fz,   z,  d12); // 19
        z = 1/2.; fz = d1*(1 - z);
        add(  fz,  0.,   z,  d12); // 20
        add(  0.,  fz,   z,  d12); // 21
        add( -fz,  0.,   z,  d12); // 22
        add(  0., -fz,   z,  d12); // 23
        z = (1 + d1)/2., fz = d1/2*(1 - z);
        add(  fz,  fz,   z,  d12); // 24
        add( -fz,  fz,   z,  d12); // 25
        add( -fz, -fz,   z,  d12); // 26
        add(  fz, -fz,   z,  d12); // 27
        break;
      }
      default:
        THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for PYRA: "<<nbGauss);
      }
      break;
    case NORM_PENTA6:
    case NORM_PENTA15:
      if (geom == NORM_PENTA6) setRefCoords( TPenta6a() );
      else                 setRefCoords( TPenta15a() );
      switch ( nbGauss ) {
      case 6: { // FPG6
        const double a = sqrt(3.)/3.;
        add( -a, .5, .5,  1/6. );
        add( -a, 0., .5,  1/6. );
        add( -a, .5, 0.,  1/6. );
        add(  a, .5, .5,  1/6. );
        add(  a, 0., .5,  1/6. );
        add(  a, .5, 0.,  1/6. ); break;
      }
      case 8: { // FPG8
        const double a = 0.577350269189626;
        add( -a, 1/3., 1/3., -27/96. );
        add( -a,  0.6,  0.2,  25/96. );
        add( -a,  0.2,  0.6,  25/96. );
        add( -a,  0.2,  0.2,  25/96. );
        add( +a, 1/3., 1/3., -27/96. );
        add( +a,  0.6,  0.2,  25/96. );
        add( +a,  0.2,  0.6,  25/96. );
        add( +a,  0.2,  0.2,  25/96. ); break;
      }
      case 21: { // FPG21
        const double d = sqrt(3/5.), c1 = 5/9., c2 = 8/9.; // d <=> alfa
        const double a = (6 + sqrt(15.))/21.;
        const double b = (6 - sqrt(15.))/21.;
        const double P1 = (155 + sqrt(15.))/2400.;
        const double P2 = (155 - sqrt(15.))/2400.;  //___
        add( -d,  1/3.,  1/3., c1*9/80. );//___
        add( -d,     a,     a, c1*P1    );
        add( -d, 1-2*a,     a, c1*P1    );
        add( -d,     a, 1-2*a, c1*P1    );//___
        add( -d,     b,     b, c1*P2    );
        add( -d, 1-2*b,     b, c1*P2    );
        add( -d,     b, 1-2*b, c1*P2    );//___
        add( 0.,  1/3.,  1/3., c2*9/80. );//___
        add( 0.,     a,     a, c2*P1    );
        add( 0., 1-2*a,     a, c2*P1    );
        add( 0.,     a, 1-2*a, c2*P1    );//___
        add( 0.,     b,     b, c2*P2    );
        add( 0., 1-2*b,     b, c2*P2    );
        add( 0.,     b, 1-2*b, c2*P2    );//___
        add(  d,  1/3.,  1/3., c1*9/80. );//___
        add(  d,     a,     a, c1*P1    );
        add(  d, 1-2*a,     a, c1*P1    );
        add(  d,     a, 1-2*a, c1*P1    );//___
        add(  d,     b,     b, c1*P2    );
        add(  d, 1-2*b,     b, c1*P2    );
        add(  d,     b, 1-2*b, c1*P2    );//___
        break;
      }
      default:
        THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for PENTA: " <<nbGauss);
      }
      break;

    case NORM_HEXA8:
    case NORM_HEXA20:
      if (geom == NORM_HEXA8) setRefCoords( THexa8a() );
      else                    setRefCoords( THexa20a() );
      switch ( nbGauss ) {
      case 8: { // FPG8
        const double a = sqrt(3.)/3.;
        add( -a, -a, -a,  1. );
        add( -a, -a,  a,  1. );
        add( -a,  a, -a,  1. );
        add( -a,  a,  a,  1. );
        add(  a, -a, -a,  1. );
        add(  a, -a,  a,  1. );
        add(  a,  a, -a,  1. );
        add(  a,  a,  a,  1. ); break;
      }
      case 27: { // FPG27
        const double a = sqrt(3/5.), c1 = 5/9., c2 = 8/9.;
        const double c12 = c1*c1, c13 = c1*c1*c1;
        const double c22 = c2*c2, c23 = c2*c2*c2;
        add( -a, -a, -a,   c13  ); // 1
        add( -a, -a, 0., c12*c2 ); // 2
        add( -a, -a,  a,   c13  ); // 3
        add( -a, 0., -a, c12*c2 ); // 4
        add( -a, 0., 0., c1*c22 ); // 5
        add( -a, 0.,  a, c12*c2 ); // 6
        add( -a,  a, -a,   c13  ); // 7
        add( -a,  a, 0., c12*c2 ); // 8
        add( -a,  a,  a,   c13  ); // 9
        add( 0., -a, -a, c12*c2 ); // 10
        add( 0., -a, 0., c1*c22 ); // 11
        add( 0., -a,  a, c12*c2 ); // 12
        add( 0., 0., -a, c1*c22 ); // 13
        add( 0., 0., 0.,   c23  ); // 14
        add( 0., 0.,  a, c1*c22 ); // 15
        add( 0.,  a, -a, c12*c2 ); // 16
        add( 0.,  a, 0., c1*c22 ); // 17
        add( 0.,  a,  a, c12*c2 ); // 18
        add(  a, -a, -a,   c13  ); // 19
        add(  a, -a, 0., c12*c2 ); // 20
        add(  a, -a,  a,   c13  ); // 21
        add(  a, 0., -a, c12*c2 ); // 22
        add(  a, 0., 0., c1*c22 ); // 23
        add(  a, 0.,  a, c12*c2 ); // 24
        add(  a,  a, -a,   c13  ); // 25
        add(  a,  a, 0., c12*c2 ); // 26
        add(  a,  a,  a,   c13  ); // 27
        break;
      }
      default:
        THROW_IK_EXCEPTION("TGaussDef: Invalid nb of gauss points for PENTA: " <<nbGauss);
      }
      break;

    default:
      THROW_IK_EXCEPTION("TGaussDef: unexpected EGeometrieElement: "<< geom);
    }

    if ( myWeights.capacity() != myWeights.size() )
      THROW_IK_EXCEPTION("TGaussDef: Not all gauss points defined");
  }
}
  
//================================================================================
/*!
 * \brief Return dimension for the given cell type
 */
//================================================================================

unsigned SauvUtilities::getDimension( INTERP_KERNEL::NormalizedCellType type )
{
  return type == NORM_ERROR ? -1 : INTERP_KERNEL::CellModel::GetCellModel( type ).getDimension();
}

//================================================================================
/*!
 * \brief Returns interlace array to transform a quadratic GIBI element to a MED one.
 *        i-th array item gives node index in GIBI connectivity for i-th MED node
 */
//================================================================================

const int * SauvUtilities::getGibi2MedQuadraticInterlace( INTERP_KERNEL::NormalizedCellType type )
{
  static std::vector<const int*> conn;
  static const int hexa20 [] = {0,6,4,2, 12,18,16,14, 7,5,3,1, 19,17,15,13, 8,11,10,9};
  static const int penta15[] = {0,2,4, 9,11,13, 1,3,5, 10,12,14, 6,8,7};
  static const int pyra13 [] = {0,2,4,6, 12, 1,3,5,7, 8,9,10,11};
  static const int tetra10[] = {0,2,4, 9, 1,3,5, 6,7,8};
  static const int quad8  [] = {0,2,4,6, 1,3,5,7};
  static const int tria6  [] = {0,2,4, 1,3,5};
  static const int seg3   [] = {0,2,1};
  if ( conn.empty() )
    {
      conn.resize( MaxMedCellType + 1, 0 );
      conn[ NORM_HEXA20 ] = hexa20;
      conn[ NORM_PENTA15] = penta15;
      conn[ NORM_PYRA13 ] = pyra13;
      conn[ NORM_TETRA10] = tetra10;
      conn[ NORM_SEG3   ] = seg3;
      conn[ NORM_TRI6   ] = tria6;
      conn[ NORM_QUAD8  ] = quad8;
    }
  return conn[ type ];
}

//================================================================================
/*!
 * \brief Avoid coping sortedNodeIDs
 */
//================================================================================

Cell::Cell(const Cell& ma)
  : _nodes(ma._nodes), _reverse(ma._reverse), _sortedNodeIDs(0), _number(ma._number)
{
  if ( ma._sortedNodeIDs )
    {
      _sortedNodeIDs = new int[ _nodes.size() ];
      std::copy( ma._sortedNodeIDs, ma._sortedNodeIDs + _nodes.size(), _sortedNodeIDs );
    }
}

//================================================================================
/*!
 * \brief Rerturn the i-th link of face
 */
//================================================================================

SauvUtilities::Link Cell::link(int i) const
{
  int i2 = ( i + 1 ) % _nodes.size();
  if ( _reverse )
    return std::make_pair( _nodes[i2]->_number, _nodes[i]->_number );
  else
    return std::make_pair( _nodes[i]->_number, _nodes[i2]->_number );
}

//================================================================================
/*!
 * \brief creates if needed and return _sortedNodeIDs
 */
//================================================================================

const TID* Cell::getSortedNodes() const
{
  if ( !_sortedNodeIDs )
    {
      size_t l=_nodes.size();
      _sortedNodeIDs = new int[ l ];

      for (size_t i=0; i!=l; ++i)
        _sortedNodeIDs[i]=_nodes[i]->_number;
      std::sort( _sortedNodeIDs, _sortedNodeIDs + l );
    }
  return _sortedNodeIDs;
}

//================================================================================
/*!
 * \brief Compare sorted ids of cell nodes
 */
//================================================================================

bool Cell::operator< (const Cell& ma) const
{
  if ( _nodes.size() == 1 )
    return _nodes[0] < ma._nodes[0];

  const int* v1 = getSortedNodes();
  const int* v2 = ma.getSortedNodes();
  for ( const int* vEnd = v1 + _nodes.size(); v1 < vEnd; ++v1, ++v2 )
    if(*v1 != *v2)
      return *v1 < *v2;
  return false;
}

//================================================================================
/*!
 * \brief Dump a Cell
 */
//================================================================================

std::ostream& SauvUtilities::operator<< (std::ostream& os, const SauvUtilities::Cell& ma)
{
  os << "cell " << ma._number << " (" << ma._nodes.size() << " nodes) : < " << ma._nodes[0]->_number;
  for( size_t i=1; i!=ma._nodes.size(); ++i)
    os << ", " << ma._nodes[i]->_number;
#ifdef _DEBUG_
  os << " > sortedNodes: ";
  if ( ma._sortedNodeIDs ) {
    os << "< ";
    for( size_t i=0; i!=ma._nodes.size(); ++i)
      os << ( i ? ", " : "" ) << ma._sortedNodeIDs[i];
    os << " >";
  }
  else {
    os << "NULL";
  }
#endif
  return os;
}

//================================================================================
/*!
 * \brief Return nb of elements in the group
 */
//================================================================================

int Group::size() const
{
  int sizze = 0;
  if ( !_relocTable.empty() )
    sizze =  _relocTable.size();
  else if ( _medGroup )
    sizze = _medGroup->getNumberOfTuples();
  else if ( !_cells.empty() )
    sizze = _cells.size();
  else
    for ( size_t i = 0; i < _groups.size(); ++i )
      sizze += _groups[i]->size();
  return sizze;
}

//================================================================================
/*!
 * \brief Conver gibi element type to med one
 */
//================================================================================

INTERP_KERNEL::NormalizedCellType SauvUtilities::gibi2medGeom( size_t gibiType )
{
  if ( gibiType < 1 || gibiType > NbGibiCellTypes )
    return NORM_ERROR;

  return GibiTypeToMed[ gibiType - 1 ];
}

//================================================================================
/*!
 * \brief Conver med element type to gibi one
 */
//================================================================================

int SauvUtilities::med2gibiGeom( INTERP_KERNEL::NormalizedCellType medGeomType )
{
  for ( unsigned int i = 0; i < NbGibiCellTypes; i++ )
    if ( GibiTypeToMed[ i ] == medGeomType )
      return i + 1;

  return -1;
}

//================================================================================
/*!
 * \brief Remember the file name
 */
//================================================================================

FileReader::FileReader(const char* fileName):_fileName(fileName),_iRead(0),_nbToRead(0)
{
}

//================================================================================
/*!
 * \brief Constructor of ASCII sauve file reader
 */
//================================================================================

ASCIIReader::ASCIIReader(const char* fileName)
  :FileReader(fileName),
   _file(-1)
{
}

//================================================================================
/*!
 * \brief Return true
 */
//================================================================================

bool ASCIIReader::isASCII() const
{
  return true;
}

//================================================================================
/*!
 * \brief Try to open an ASCII file
 */
//================================================================================

bool ASCIIReader::open()
{
#ifdef WIN32
  _file = ::_open (_fileName.c_str(), _O_RDONLY|_O_BINARY);
#else
  _file = ::open (_fileName.c_str(), O_RDONLY);
#endif
  if (_file >= 0)
    {
      _start  = new char [GIBI_BufferSize]; // working buffer beginning
      //_tmpBuf = new char [GIBI_MaxOutputLen];
      _ptr    = _start;
      _eptr   = _start;
      _lineNb = 0;
    }
  else
    {
      //THROW_IK_EXCEPTION("Can't open file "<<_fileName << " fd: " << _file);
    }
  return (_file >= 0);
}

//================================================================================
/*!
 * \brief Close the file
 */
//================================================================================

ASCIIReader::~ASCIIReader()
{
  if (_file >= 0)
    {
      ::close (_file);
      if (_start != 0L)
        {
          delete [] _start;
          //delete [] _tmpBuf;
          _start = 0;
        }
      _file = -1;
    }
}

//================================================================================
/*!
 * \brief Return a next line of the file
 */
//================================================================================

bool ASCIIReader::getNextLine (char* & line, bool raiseOEF /*= true*/ )
{
  if ( getLine( line )) return true;
  if ( raiseOEF )
    THROW_IK_EXCEPTION("Unexpected EOF on ln "<<_lineNb);
  return false;
}

//================================================================================
/*!
 * \brief Read a next line of the file if necessary
 */
//================================================================================

bool ASCIIReader::getLine(char* & line)
{
  bool aResult = true;
  // Check the state of the buffer;
  // if there is too little left, read the next portion of data
  int nBytesRest = _eptr - _ptr;
  if (nBytesRest < GIBI_MaxOutputLen)
    {
      if (nBytesRest > 0)
        {
          // move the remaining portion to the buffer beginning
          for ( int i = 0; i < nBytesRest; ++i )
            _start[i] = _ptr[i];
          //memcpy (_tmpBuf, _ptr, nBytesRest);
          //memcpy (_start, _tmpBuf, nBytesRest);
        }
      else
        {
          nBytesRest = 0;
        }
      _ptr = _start;
      const int nBytesRead = ::read (_file,
                                     &_start [nBytesRest],
                                     GIBI_BufferSize - nBytesRest);
      nBytesRest += nBytesRead;
      _eptr = &_start [nBytesRest];
    }
  // Check the buffer for the end-of-line
  char * ptr = _ptr;
  while (true)
    {
      // Check for end-of-the-buffer, the ultimate criterion for termination
      if (ptr >= _eptr)
        {
          if (nBytesRest <= 0)
            aResult = false;
          else
            _eptr[-1] = '\0';
          break;
        }
      // seek the line-feed character
      if (ptr[0] == '\n')
        {
          if (ptr[-1] == '\r')
            ptr[-1] = '\0';
          ptr[0] = '\0';
          ++ptr;
          break;
        }
      ++ptr;
    }
  // Output the result
  line = _ptr;
  _ptr = ptr;
  _lineNb++;

  return aResult;
}

//================================================================================
/*!
 * \brief Prepare for iterating over given nb of values
 *  \param nbToRead - nb of fields to read
 *  \param nbPosInLine - nb of fields in one line
 *  \param width - field width
 *  \param shift - shift from line beginning to the field start
 */
//================================================================================

void ASCIIReader::init( int nbToRead, int nbPosInLine, int width, int shift /*= 0*/ )
{
  _nbToRead    = nbToRead;
  _nbPosInLine = nbPosInLine;
  _width       = width;
  _shift       = shift;
  _iPos = _iRead = 0;
  if ( _nbToRead )
    {
      getNextLine( _curPos );
      _curPos = _curPos + _shift;
    }
  else
    {
      _curPos = 0;
    }
  _curLocale.clear();
}

//================================================================================
/*!
 * \brief Prepare for iterating over given nb of string values
 */
//================================================================================

void ASCIIReader::initNameReading(int nbValues, int width /*= 8*/)
{
  init( nbValues, 72 / ( width + 1 ), width, 1 );
}

//================================================================================
/*!
 * \brief Prepare for iterating over given nb of integer values
 */
//================================================================================

void ASCIIReader::initIntReading(int nbValues)
{
  init( nbValues, 10, 8 );
}

//================================================================================
/*!
 * \brief Prepare for iterating over given nb of real values
 */
//================================================================================

void ASCIIReader::initDoubleReading(int nbValues)
{
  init( nbValues, 3, 22 );

  // Correction 2 of getDouble(): set "C" numeric locale to read numbers
  // with dot decimal point separator, as it is in SAUVE files
  _curLocale = setlocale(LC_NUMERIC, "C");
}

//================================================================================
/*!
 * \brief Return true if not all values have been read
 */
//================================================================================

bool ASCIIReader::more() const
{
  bool result = false;
  if ( _iRead < _nbToRead)
    {
      if ( _curPos ) result = true;
    }
  return result;
}

//================================================================================
/*!
 * \brief Go to the nex value
 */
//================================================================================

void ASCIIReader::next()
{
  if ( !more() )
    THROW_IK_EXCEPTION("SauvUtilities::ASCIIReader::next(): no more() values to read");
  ++_iRead;
  ++_iPos;
  if ( _iRead < _nbToRead )
    {
      if ( _iPos >= _nbPosInLine )
        {
          getNextLine( _curPos );
          _curPos = _curPos + _shift;
          _iPos = 0;
        }
      else
        {
          _curPos = _curPos + _width + _shift;
        }
    }
  else
    {
      _curPos = 0;
      if ( !_curLocale.empty() )
        {
          setlocale(LC_NUMERIC, _curLocale.c_str());
          _curLocale.clear();
        }
    }
}

//================================================================================
/*!
 * \brief Return the current integer value
 */
//================================================================================

int ASCIIReader::getInt() const
{
  // fix for two glued ints (issue 0021009):
  // Line nb    |   File contents
  // ------------------------------------------------------------------------------------
  // 53619905   |       1       2       6       8
  // 53619906   |                                                                SCALAIRE
  // 53619907   |    -63312600499       1       0       0       0      -2       0       2
  //   where -63312600499 is actualy -633 and 12600499
  char hold=_curPos[_width];
  _curPos[_width] = '\0';
  int result = atoi( _curPos );
  _curPos[_width] = hold;
  return result;
  //return atoi(str());
}

//================================================================================
/*!
 * \brief Return the current float value
 */
//================================================================================

float ASCIIReader::getFloat() const
{
  return getDouble();
}

//================================================================================
/*!
 * \brief Return the current double value
 */
//================================================================================

double ASCIIReader::getDouble() const
{
  //std::string aStr (_curPos);

  // Correction: add missing 'E' specifier
  // int aPosStart = aStr.find_first_not_of(" \t");
  // if (aPosStart < (int)aStr.length()) {
  //   int aPosSign = aStr.find_first_of("+-", aPosStart + 1); // pass the leading symbol, as it can be a sign
  //   if (aPosSign < (int)aStr.length()) {
  //     if (aStr[aPosSign - 1] != 'e' && aStr[aPosSign - 1] != 'E')
  //       aStr.insert(aPosSign, "E", 1);
  //   }
  // }

  // Different Correction (more optimal)
  // Sample:
  //  0.00000000000000E+00 -2.37822406690632E+01  6.03062748797469E+01
  //  7.70000000000000-100  7.70000000000000+100  7.70000000000000+100
  //0123456789012345678901234567890123456789012345678901234567890123456789
  const size_t posE = 18;
  std::string aStr (_curPos);
  if ( aStr.find('E') == std::string::npos && aStr.find('e') == std::string::npos )
    {
      if ( aStr.size() < posE+1 )
        THROW_IK_EXCEPTION("No more doubles (line #" << lineNb() << ")");
      aStr.insert( posE, "E", 1 );
      return atof(aStr.c_str());
    }
  return atof( _curPos );
}

//================================================================================
/*!
 * \brief Return the current string value
 */
//================================================================================

std::string ASCIIReader::getName() const
{
  int len = _width;
  while (( _curPos[len-1] == ' ' || _curPos[len-1] == 0) && len > 0 )
    len--;
  return std::string( _curPos, len );
}

//================================================================================
/*!
 * \brief Constructor of a binary sauve file reader
 */
//================================================================================

XDRReader::XDRReader(const char* fileName) :FileReader(fileName), _xdrs_file(NULL)
{
}

//================================================================================
/*!
 * \brief Close the XDR sauve file
 */
//================================================================================

XDRReader::~XDRReader()
{
#ifdef HAS_XDR  
  if ( _xdrs_file )
    {
      xdr_destroy((XDR*)_xdrs);
      free((XDR*)_xdrs);
      ::fclose(_xdrs_file);
      _xdrs_file = NULL;
    }
#endif
}

//================================================================================
/*!
 * \brief Return false
 */
//================================================================================

bool XDRReader::isASCII() const
{
  return false;
}

//================================================================================
/*!
 * \brief Try to open an XRD file
 */
//================================================================================

bool XDRReader::open()
{
  bool xdr_ok = false;
#ifdef HAS_XDR
#ifdef WIN32
  if ((_xdrs_file = ::fopen(_fileName.c_str(), "rb")))
#else 
    if ((_xdrs_file = ::fopen(_fileName.c_str(), "r")))
#endif
      {
        _xdrs = (XDR *)malloc(sizeof(XDR));
        xdrstdio_create((XDR*)_xdrs, _xdrs_file, XDR_DECODE);

        const int maxsize = 10;
        char icha[maxsize+1];
        char* icha2 = icha;
        if (( xdr_ok = xdr_string((XDR*)_xdrs, &icha2, maxsize)))
          {
            icha[maxsize] = '\0';
            xdr_ok = (strcmp(icha, "CASTEM XDR") == 0);
          }
        if ( !xdr_ok )
          {
            xdr_destroy((XDR*)_xdrs);
            free((XDR*)_xdrs);
            fclose(_xdrs_file);
            _xdrs_file = NULL;
          }
      }
#endif
  return xdr_ok;
}

//================================================================================
/*!
 * \brief A stub
 */
//================================================================================

bool XDRReader::getNextLine (char* &, bool )
{
  return true;
}

//================================================================================
/*!
 * \brief Prepare for iterating over given nb of values
 *  \param nbToRead - nb of fields to read
 *  \param width - field width
 */
//================================================================================

void XDRReader::init( int nbToRead, int width/*=0*/ )
{
  if(_iRead < _nbToRead)
    {
      std::cout << "_iRead, _nbToRead : " << _iRead << " " << _nbToRead << std::endl;
      std::cout << "Unfinished iteration before new one !" << std::endl;
      THROW_IK_EXCEPTION("SauvUtilities::XDRReader::init(): Unfinished iteration before new one !");
    }
  _iRead    = 0;
  _nbToRead = nbToRead;
  _width    = width;
}

//================================================================================
/*!
 * \brief Prepare for iterating over given nb of string values
 */
//================================================================================

void XDRReader::initNameReading(int nbValues, int width)
{
  init( nbValues, width );
  _xdr_kind = _xdr_kind_char;
  if(nbValues*width)
    {
      unsigned int nels = nbValues*width;
      _xdr_cvals = (char*)malloc((nels+1)*sizeof(char));
#ifdef HAS_XDR
      xdr_string((XDR*)_xdrs, &_xdr_cvals, nels);
#endif
      _xdr_cvals[nels] = '\0';
    }
}

//================================================================================
/*!
 * \brief Prepare for iterating over given nb of integer values
 */
//================================================================================

void XDRReader::initIntReading(int nbValues)
{
  init( nbValues );
  _xdr_kind = _xdr_kind_int;
  if(nbValues)
    {
#ifdef HAS_XDR
      unsigned int nels = nbValues;
      unsigned int actual_nels;
      _xdr_ivals = (int*)malloc(nels*sizeof(int));
      xdr_array((XDR*)_xdrs, (char **)&_xdr_ivals, &actual_nels, nels, sizeof(int), (xdrproc_t)xdr_int);
#endif
    }
}

//================================================================================
/*!
 * \brief Prepare for iterating over given nb of real values
 */
//================================================================================

void XDRReader::initDoubleReading(int nbValues)
{
  init( nbValues );
  _xdr_kind = _xdr_kind_double;
  if(nbValues)
    {
#ifdef HAS_XDR
      unsigned int nels = nbValues;
      unsigned int actual_nels;
      _xdr_dvals = (double*)malloc(nels*sizeof(double));
      xdr_array((XDR*)_xdrs, (char **)&_xdr_dvals, &actual_nels, nels, sizeof(double), (xdrproc_t)xdr_double);
#endif
    }
}

//================================================================================
/*!
 * \brief Return true if not all values have been read
 */
//================================================================================

bool XDRReader::more() const
{
  return _iRead < _nbToRead;
}

//================================================================================
/*!
 * \brief Go to the nex value
 */
//================================================================================

void XDRReader::next()
{
  if ( !more() )
    THROW_IK_EXCEPTION("SauvUtilities::XDRReader::next(): no more() values to read");

  ++_iRead;
  if ( _iRead < _nbToRead )
    {
    }
  else
    {
      if(_xdr_kind == _xdr_kind_char) free(_xdr_cvals);
      if(_xdr_kind == _xdr_kind_int) free(_xdr_ivals);
      if(_xdr_kind == _xdr_kind_double) free(_xdr_dvals);
      _xdr_kind = _xdr_kind_null;
    }
}

//================================================================================
/*!
 * \brief Return the current integer value
 */
//================================================================================

int XDRReader::getInt() const
{
  if(_iRead < _nbToRead)
    {
      return _xdr_ivals[_iRead];
    }
  else
    {
      int result = 0;
#ifdef HAS_XDR
      xdr_int((XDR*)_xdrs, &result);
#endif
      return result;
    }
}

//================================================================================
/*!
 * \brief Return the current float value
 */
//================================================================================

float  XDRReader::getFloat() const
{
  float result = 0;
#ifdef HAS_XDR
  xdr_float((XDR*)_xdrs, &result);
#endif
  return result;
}

//================================================================================
/*!
 * \brief Return the current double value
 */
//================================================================================

double XDRReader::getDouble() const
{
  if(_iRead < _nbToRead)
    {
      return _xdr_dvals[_iRead];
    }
  else
    {
      double result = 0;
#ifdef HAS_XDR
      xdr_double((XDR*)_xdrs, &result);
#endif
      return result;
    }
}

//================================================================================
/*!
 * \brief Return the current string value
 */
//================================================================================

std::string XDRReader::getName() const
{
  int len = _width;
  char* s = _xdr_cvals + _iRead*_width;
  while (( s[len-1] == ' ' || s[len-1] == 0) && len > 0 )
    len--;
  return std::string( s, len );
}

//================================================================================
/*!
 * \brief Throw an exception if not all needed data is present
 */
//================================================================================

void IntermediateMED::checkDataAvailability() const
{
  if ( _spaceDim == 0 )
    THROW_IK_EXCEPTION("Wrong file format"); // it is the first record in the sauve file

  if ( _groups.empty() )
    THROW_IK_EXCEPTION("No elements have been read");

  if ( _points.empty() || _nbNodes == 0 )
    THROW_IK_EXCEPTION("Nodes of elements are not filled");

  if ( _coords.empty() )
    THROW_IK_EXCEPTION("Node coordinates are missing");

  if ( _coords.size() < _nbNodes * _spaceDim )
    THROW_IK_EXCEPTION("Nodes and coordinates mismatch");
}

//================================================================================
/*!
 * \brief Safely adds a new Group
 */
//================================================================================

Group* IntermediateMED::addNewGroup(std::vector<SauvUtilities::Group*>* groupsToFix)
{
  if ( _groups.size() == _groups.capacity() ) // re-allocation would occure
    {
      std::vector<Group> newGroups( _groups.size() );
      newGroups.push_back( Group() );

      for ( size_t i = 0; i < _groups.size(); ++i )
        {
          // avoid copying _cells
          std::vector<const Cell*> cells;
          cells.swap( _groups[i]._cells );
          newGroups[i] = _groups[i];
          newGroups[i]._cells.swap( cells );

          // correct pointers to sub-groups
          for ( size_t j = 0; j < _groups[i]._groups.size(); ++j )
            {
              int iG = _groups[i]._groups[j] - &_groups[0];
              newGroups[i]._groups[j] = & newGroups[ iG ];
            }
        }

      // fix given groups
      if ( groupsToFix )
        for ( size_t i = 0; i < groupsToFix->size(); ++i )
          if ( (*groupsToFix)[i] )
            {
              int iG = (*groupsToFix)[i] - &_groups[0];
              (*groupsToFix)[i] = & newGroups[ iG ];
            }

      // fix field supports
      for ( int isNode = 0; isNode < 2; ++isNode )
        {
          std::vector<DoubleField* >& fields = isNode ? _nodeFields : _cellFields;
          for ( size_t i = 0; i < fields.size(); ++i )
            {
              if ( !fields[i] ) continue;
              for ( size_t j = 0; j < fields[i]->_sub.size(); ++j )
                if ( fields[i]->_sub[j]._support )
                  {
                    int iG = fields[i]->_sub[j]._support - &_groups[0];
                    fields[i]->_sub[j]._support = & newGroups[ iG ];
                  }
              if ( fields[i]->_group )
                {
                  int iG = fields[i]->_group - &_groups[0];
                  fields[i]->_group = & newGroups[ iG ];
                }
            }
        }

      _groups.swap( newGroups );
    }
  else
    {
      _groups.push_back( Group() );
    }
  return &_groups.back();
}

//================================================================================
/*!
 * \brief Makes MEDCoupling::MEDFileData from self
 */
//================================================================================

MEDCoupling::MEDFileData* IntermediateMED::convertInMEDFileDS()
{
  MCAuto< MEDFileUMesh >  mesh   = makeMEDFileMesh();
  MCAuto< MEDFileFields > fields = makeMEDFileFields(mesh);

  MCAuto< MEDFileMeshes > meshes = MEDFileMeshes::New();
  MCAuto< MEDFileData >  medData = MEDFileData::New();
  meshes->pushMesh( mesh );
  medData->setMeshes( meshes );
  if ( fields ) medData->setFields( fields );

  return medData.retn();
}

//================================================================================
/*!
 * \brief Creates MEDCoupling::MEDFileUMesh from its data
 */
//================================================================================

MEDCoupling::MEDFileUMesh* IntermediateMED::makeMEDFileMesh()
{
  // check if all needed piles are present
  checkDataAvailability();

  decreaseHierarchicalDepthOfSubgroups();

  // set long names (before orienting!)
  setGroupLongNames();

  // fix element orientation
  if ( _spaceDim == 2 || _spaceDim == 1 )
    orientElements2D();
  else if ( _spaceDim == 3 )
    orientElements3D();

  // process groups
  eraseUselessGroups();
  //detectMixDimGroups();

  // assign IDs
  _points.numberNodes();
  numberElements();

  // make the med mesh

  MEDFileUMesh* mesh = MEDFileUMesh::New();

  DataArrayDouble *coords = getCoords();
  setConnectivity( mesh, coords );
  setGroups( mesh );

  coords->decrRef();

  if ( !mesh->getName().c_str() || strlen( mesh->getName().c_str() ) == 0 )
    mesh->setName( "MESH" );

  return mesh;
}

//================================================================================
/*!
 * \brief Set long names to groups
 */
//================================================================================

void IntermediateMED::setGroupLongNames()
{
  if ( _listGIBItoMED_mail.empty() )
    return;

  // IMP 0023285: only keep the meshes named in the table MED_MAIL
  // clear all group names
  for ( size_t i = 0; i < _groups.size(); ++i )
    if ( !_groups[i]._isProfile )
      _groups[i]._name.clear();


  // IMP 0020434: mapping GIBI names to MED names
  // set med names to objects (mesh, fields, support, group or other)

  std::set<int> treatedGroups;

  std::list<nameGIBItoMED>::iterator itGIBItoMED = _listGIBItoMED_mail.begin();
  for (; itGIBItoMED != _listGIBItoMED_mail.end(); itGIBItoMED++)
    {
      if ( (int)_groups.size() < itGIBItoMED->gibi_id ) continue;

      SauvUtilities::Group & grp = _groups[itGIBItoMED->gibi_id - 1];

      // if there are several names for grp then the 1st name is the name
      // of grp and the rest ones are names of groups referring grp (issue 0021311)
      const bool isRefName = !treatedGroups.insert( itGIBItoMED->gibi_id ).second;
      if ( !isRefName )
        {
          grp._name = _mapStrings[ itGIBItoMED->med_id ];
        }
      else if ( !grp._refNames.empty() && grp._refNames.back().empty() )
        {
          for ( unsigned i = 0; i < grp._refNames.size(); ++i )
            if ( grp._refNames[i].empty() )
              grp._refNames[i] = _mapStrings[ (*itGIBItoMED).med_id ];
        }
      else
        {
          grp._refNames.push_back( _mapStrings[ (*itGIBItoMED).med_id ]);
        }
    }

  // IMP 0023285: only keep the meshes named in the table MED_MAIL
  // remove all cells belonging to non-named groups only

  // use Cell::_reverse to mark cells to keep
  for ( size_t i = 0; i < _groups.size(); ++i )
    {
      SauvUtilities::Group & grp = _groups[i];
      if ( grp._isProfile || !grp._name.empty() )
        {
          for ( size_t iC = 0; iC < grp._cells.size(); ++iC )
            grp._cells[iC]->_reverse = true;

          for (size_t j = 0; j < grp._groups.size(); ++j )
            for ( size_t iC = 0; iC < grp._groups[j]->_cells.size(); ++iC )
              grp._groups[j]->_cells[iC]->_reverse = true;
        }
    }
  // remove non-marked cells (with _reverse == false)
  CellsByDimIterator cellsIt( *this );
  while ( cellsIt.nextType() )
    {
      std::set<Cell> & cells = _cellsByType[ cellsIt.type() ];
      std::set<Cell>::iterator cIt = cells.begin();
      while ( cIt != cells.end() )
        if ( cIt->_reverse )
          {
            cIt->_reverse = false;
            ++cIt;
          }
        else
          {
            cells.erase( cIt++ );
          }
    }
}

//================================================================================
/*!
 * \brief Set long names to fields
 */
//================================================================================

void IntermediateMED::setFieldLongNames(std::set< std::string >& usedNames)
{
  std::list<nameGIBItoMED>::iterator itGIBItoMED = _listGIBItoMED_cham.begin();
  for (; itGIBItoMED != _listGIBItoMED_cham.end(); itGIBItoMED++)
    {
      if (itGIBItoMED->gibi_pile == PILE_FIELD)
        {
          _cellFields[itGIBItoMED->gibi_id - 1]->_name = _mapStrings[itGIBItoMED->med_id];
        }
      else if (itGIBItoMED->gibi_pile == PILE_NODES_FIELD)
        {
          _nodeFields[itGIBItoMED->gibi_id - 1]->_name = _mapStrings[itGIBItoMED->med_id];
        }
    } // iterate on _listGIBItoMED_cham

  for (itGIBItoMED =_listGIBItoMED_comp.begin(); itGIBItoMED != _listGIBItoMED_comp.end(); itGIBItoMED++)
    {
      std::string medName  = _mapStrings[itGIBItoMED->med_id];
      std::string gibiName = _mapStrings[itGIBItoMED->gibi_id];

      bool name_found = false;
      for ( int isNodal = 0; isNodal < 2 && !name_found; ++isNodal )
        {
          std::vector<DoubleField* > & fields = isNodal ? _nodeFields : _cellFields;
          for ( size_t ifi = 0; ifi < fields.size() && !name_found; ifi++)
            {
              if (medName.find( fields[ifi]->_name + "." ) == 0 )
                {
                  std::vector<DoubleField::_Sub_data>& aSubDs = fields[ifi]->_sub;
                  int nbSub = aSubDs.size();
                  for (int isu = 0; isu < nbSub; isu++)
                    for (int ico = 0; ico < aSubDs[isu].nbComponents(); ico++)
                      {
                        if (aSubDs[isu].compName(ico) == gibiName)
                          {
                            std::string medNameCompo = medName.substr( fields[ifi]->_name.size() + 1 );
                            fields[ifi]->_sub[isu].compName(ico) = medNameCompo;
                          }
                      }
                }
            }
        }
    } // iterate on _listGIBItoMED_comp

  for ( size_t i = 0; i < _nodeFields.size() ; i++)
    usedNames.insert( _nodeFields[i]->_name );
  for ( size_t i = 0; i < _cellFields.size() ; i++)
    usedNames.insert( _cellFields[i]->_name );
}

//================================================================================
/*!
 * \brief Decrease hierarchical depth of subgroups
 */
//================================================================================

void IntermediateMED::decreaseHierarchicalDepthOfSubgroups()
{
  for (size_t i=0; i!=_groups.size(); ++i)
    {
      Group& grp = _groups[i];
      for (size_t j = 0; j < grp._groups.size(); ++j )
        {
          Group & sub_grp = *grp._groups[j];
          if ( !sub_grp._groups.empty() )
            {
              // replace j with its 1st subgroup
              grp._groups[j] = sub_grp._groups[0];
              // push back the rest subs
              grp._groups.insert( grp._groups.end(), ++sub_grp._groups.begin(), sub_grp._groups.end() );
            }
        }
      // remove empty sub-_groups
      std::vector< Group* > newSubGroups;
      newSubGroups.reserve( grp._groups.size() );
      for (size_t j = 0; j < grp._groups.size(); ++j )
        if ( !grp._groups[j]->empty() )
          newSubGroups.push_back( grp._groups[j] );
      if ( newSubGroups.size() < grp._groups.size() )
        grp._groups.swap( newSubGroups );
    }
}

//================================================================================
/*!
 * \brief Erase _groups that won't be converted
 */
//================================================================================

void IntermediateMED::eraseUselessGroups()
{
  // propagate _isProfile=true to sub-groups of composite groups
  // for (size_t int i=0; i!=_groups.size(); ++i)
  // {
  //   Group* grp = _groups[i];
  //   if ( grp->_isProfile && !grp->_groups.empty() )
  //     for (size_t j = 0; j < grp->_groups.size(); ++j )
  //       grp->_groups[j]->_isProfile=true;
  // }
  std::set<Group*> groups2convert;
  // keep not named sub-groups of field supports
  for (size_t i=0; i!=_groups.size(); ++i)
    {
      Group& grp = _groups[i];
      if ( grp._isProfile && !grp._groups.empty() )
        groups2convert.insert( grp._groups.begin(), grp._groups.end() );
    }

  // keep named groups and their subgroups
  for (size_t i=0; i!=_groups.size(); ++i)
    {
      Group& grp = _groups[i];
      if ( !grp._name.empty() && !grp.empty() )
        {
          groups2convert.insert( &grp );
          groups2convert.insert( grp._groups.begin(), grp._groups.end() );
        }
    }
  // erase groups that are not in groups2convert and not _isProfile
  for (size_t i=0; i!=_groups.size(); ++i)
    {
      Group* grp = &_groups[i];
      if ( !grp->_isProfile && !groups2convert.count( grp ) )
        {
          grp->_cells.clear();
          grp->_groups.clear();
        }
    }
}

//================================================================================
/*!
 * \brief Detect _groups of mixed dimension
 */
//================================================================================

void IntermediateMED::detectMixDimGroups()
{
  //hasMixedCells = false;
  for ( size_t i=0; i < _groups.size(); ++i )
    {
      Group& grp = _groups[i];
      if ( grp._groups.size() < 2 )
        continue;

      // check if sub-groups have different dimension
      unsigned dim1 = getDim( &grp );
      for ( size_t j = 1; j  < grp._groups.size(); ++j )
        {
          unsigned dim2 = getDim( grp._groups[j] );
          if ( dim1 != dim2 )
            {
              grp._cells.clear();
              grp._groups.clear();
              if ( !grp._name.empty() )
                std::cout << "Erase a group with elements of different dim |" << grp._name << "|"<< std::endl;
              break;
            }
        }
    }
}

//================================================================================
/*!
 * \brief Fix connectivity of elements in 2D space
 */
//================================================================================

void IntermediateMED::orientElements2D()
{
  std::set<Cell>::const_iterator elemIt, elemEnd;
  std::vector< std::pair<int,int> > swapVec;

  // ------------------------------------
  // fix connectivity of quadratic edges
  // ------------------------------------
  std::set<Cell>& quadEdges = _cellsByType[ INTERP_KERNEL::NORM_SEG3 ];
  if ( !quadEdges.empty() )
    {
      elemIt = quadEdges.begin(), elemEnd = quadEdges.end();
      for ( ; elemIt != elemEnd; ++elemIt )
        ConvertQuadratic( INTERP_KERNEL::NORM_SEG3, *elemIt );
    }

  CellsByDimIterator faceIt( *this, 2 );
  while ( const std::set<Cell > * faces = faceIt.nextType() )
    {
      TCellType cellType = faceIt.type();
      bool isQuadratic = getGibi2MedQuadraticInterlace( cellType );

      getReverseVector( cellType, swapVec );

      // ------------------------------------
      // fix connectivity of quadratic faces
      // ------------------------------------
      if ( isQuadratic )
        for ( elemIt = faces->begin(), elemEnd = faces->end(); elemIt != elemEnd; elemIt++ )
          ConvertQuadratic( cellType, *elemIt );

      // --------------------------
      // orient faces clockwise
      // --------------------------
      // COMMENTED for issue 0022612 note 17739
      // int iQuad = isQuadratic ? 2 : 1;
      // for ( elemIt = faces->begin(), elemEnd = faces->end(); elemIt != elemEnd; elemIt++ )
      //   {
      //     // look for index of the most left node
      //     int iLeft = 0, iNode, nbNodes = elemIt->_nodes.size() / iQuad;
      //     double x, minX = nodeCoords( elemIt->_nodes[0] )[0];
      //     for ( iNode = 1; iNode < nbNodes; ++iNode )
      //       if (( x = nodeCoords( elemIt->_nodes[ iNode ])[ 0 ]) < minX )
      //         minX = x, iLeft = iNode;

      //     // indeces of the nodes neighboring the most left one
      //     int iPrev = ( iLeft - 1 < 0 ) ? nbNodes - 1 : iLeft - 1;
      //     int iNext = ( iLeft + 1 == nbNodes ) ? 0 : iLeft + 1;
      //     // find components of prev-left and left-next vectors
      //     double xP = nodeCoords( elemIt->_nodes[ iPrev ])[ 0 ];
      //     double yP = nodeCoords( elemIt->_nodes[ iPrev ])[ 1 ];
      //     double xN = nodeCoords( elemIt->_nodes[ iNext ])[ 0 ];
      //     double yN = nodeCoords( elemIt->_nodes[ iNext ])[ 1 ];
      //     double xL = nodeCoords( elemIt->_nodes[ iLeft ])[ 0 ];
      //     double yL = nodeCoords( elemIt->_nodes[ iLeft ])[ 1 ];
      //     double xPL = xL - xP, yPL = yL - yP; // components of prev-left vector
      //     double xLN = xN - xL, yLN = yN - yL; // components of left-next vector
      //     // normalise y of the vectors
      //     double modPL = sqrt ( xPL * xPL + yPL * yPL );
      //     double modLN = sqrt ( xLN * xLN + yLN * yLN );
      //     if ( modLN > std::numeric_limits<double>::min() &&
      //          modPL > std::numeric_limits<double>::min() )
      //       {
      //         yPL /= modPL;
      //         yLN /= modLN;
      //         // summary direction of neighboring links must be positive
      //         bool clockwise = ( yPL + yLN > 0 );
      //         if ( !clockwise )
      //           reverse( *elemIt, swapVec );
      //       }
      //   }
    }
}

//================================================================================
/*!
 * \brief Fix connectivity of elements in 3D space
 */
//================================================================================

void IntermediateMED::orientElements3D()
{
  // set _reverse flags of faces
  // COMMENTED for issue 0022612 note 17739
  //orientFaces3D();

  // -----------------
  // fix connectivity
  // -----------------

  std::set<Cell>::const_iterator elemIt, elemEnd;
  std::vector< std::pair<int,int> > swapVec;

  for ( int dim = 1; dim <= 3; ++dim )
    {
      CellsByDimIterator cellsIt( *this, dim );
      while ( const std::set<Cell > * elems = cellsIt.nextType() )
        {
          TCellType cellType = cellsIt.type();
          bool isQuadratic = getGibi2MedQuadraticInterlace( cellType );
          getReverseVector( cellType, swapVec );

          elemIt = elems->begin(), elemEnd = elems->end();
          for ( ; elemIt != elemEnd; elemIt++ )
            {
              // GIBI connectivity -> MED one
              if( isQuadratic )
                ConvertQuadratic( cellType, *elemIt );

              // reverse faces
              if ( elemIt->_reverse )
                reverse ( *elemIt, swapVec );
            }
        }
    }

  // COMMENTED for issue 0022612 note 17739
  //orientVolumes();
}

//================================================================================
/*!
 * \brief Orient equally (by setting _reverse flag) all connected faces in 3D space
 */
//================================================================================

void IntermediateMED::orientFaces3D()
{
  // fill map of links and their faces
  std::set<const Cell*> faces;
  std::map<const Cell*, Group*> fgm;
  std::map<Link, std::list<const Cell*> > linkFacesMap;
  std::map<Link, std::list<const Cell*> >::iterator lfIt, lfIt2;

  for (size_t i=0; i!=_groups.size(); ++i)
    {
      Group& grp = _groups[i];
      if ( !grp._cells.empty() && getDimension( grp._cellType ) == 2 )
        for ( size_t j = 0; j < grp._cells.size(); ++j )
          if ( faces.insert( grp._cells[j] ).second )
            {
              for ( size_t k = 0; k < grp._cells[j]->_nodes.size(); ++k )
                linkFacesMap[ grp._cells[j]->link( k ) ].push_back( grp._cells[j] );
              fgm.insert( std::make_pair( grp._cells[j], &grp ));
            }
    }
  // dump linkFacesMap
  //     for ( lfIt = linkFacesMap.begin(); lfIt!=linkFacesMap.end(); lfIt++) {
  //       cout<< "LINK: " << lfIt->first.first << "-" << lfIt->first.second << std::endl;
  //       std::list<const Cell*> & fList = lfIt->second;
  //       std::list<const Cell*>::iterator fIt = fList.begin();
  //       for ( ; fIt != fList.end(); fIt++ )
  //         cout << "\t" << **fIt << fgm[*fIt]->nom << std::endl;
  //     }

  // Each oriented link must appear in one face only, else a face is reversed.

  std::queue<const Cell*> faceQueue; /* the queue contains well oriented faces
                                        whose neighbors orientation is to be checked */
  bool manifold = true;
  while ( !linkFacesMap.empty() )
    {
      if ( faceQueue.empty() )
        {
          assert( !linkFacesMap.begin()->second.empty() );
          faceQueue.push( linkFacesMap.begin()->second.front() );
        }
      while ( !faceQueue.empty() )
        {
          const Cell* face = faceQueue.front();
          faceQueue.pop();

          // loop on links of <face>
          for ( int i = 0; i < (int)face->_nodes.size(); ++i )
            {
              Link link = face->link( i );
              // find the neighbor faces
              lfIt = linkFacesMap.find( link );
              int nbFaceByLink = 0;
              std::list< const Cell* > ml;
              if ( lfIt != linkFacesMap.end() )
                {
                  std::list<const Cell*> & fList = lfIt->second;
                  std::list<const Cell*>::iterator fIt = fList.begin();
                  assert( fIt != fList.end() );
                  for ( ; fIt != fList.end(); fIt++, nbFaceByLink++ )
                    {
                      ml.push_back( *fIt );
                      if ( *fIt != face ) // wrongly oriented neighbor face
                        {
                          const Cell* badFace = *fIt;
                          // reverse and remove badFace from linkFacesMap
                          for ( int j = 0; j < (int)badFace->_nodes.size(); ++j )
                            {
                              Link badlink = badFace->link( j );
                              if ( badlink == link ) continue;
                              lfIt2 = linkFacesMap.find( badlink );
                              if ( lfIt2 != linkFacesMap.end() )
                                {
                                  std::list<const Cell*> & ff = lfIt2->second;
                                  std::list<const Cell*>::iterator lfIt3 = find( ff.begin(), ff.end(), badFace );
                                  // check if badFace has been found,
                                  // else we can't erase it
                                  // case of degenerated face in edge
                                  if (lfIt3 != ff.end())
                                    {
                                      ff.erase( lfIt3 );
                                      if ( ff.empty() )
                                        linkFacesMap.erase( lfIt2 );
                                    }
                                }
                            }
                          badFace->_reverse = true; // reverse
                          //INFOS_MED( "REVERSE " << *badFace );
                          faceQueue.push( badFace );
                        }
                    }
                  linkFacesMap.erase( lfIt );
                }
              // add good neighbors to the queue
              Link revLink( link.second, link.first );
              lfIt = linkFacesMap.find( revLink );
              if ( lfIt != linkFacesMap.end() )
                {
                  std::list<const Cell*> & fList = lfIt->second;
                  std::list<const Cell*>::iterator fIt = fList.begin();
                  for ( ; fIt != fList.end(); fIt++, nbFaceByLink++ )
                    {
                      ml.push_back( *fIt );
                      if ( *fIt != face )
                        faceQueue.push( *fIt );
                    }
                  linkFacesMap.erase( lfIt );
                }
              if ( nbFaceByLink > 2 )
                {
                  if ( manifold )
                    {
                      std::list<const Cell*>::iterator ii = ml.begin();
                      std::cout << nbFaceByLink << " faces by 1 link:" << std::endl;
                      for( ; ii!= ml.end(); ii++ )
                        std::cout << "in sub-mesh <" << fgm[ *ii ]->_name << "> " << **ii << std::endl;
                    }
                  manifold = false;
                }
            } // loop on links of the being checked face
        } // loop on the face queue
    } // while ( !linkFacesMap.empty() )

  if ( !manifold )
    std::cout << " -> Non manifold mesh, faces orientation may be incorrect" << std::endl;
}

//================================================================================
/*!
 * \brief Orient volumes according to MED conventions:
 * normal of a bottom (first) face should be outside
 */
//================================================================================

void IntermediateMED::orientVolumes()
{
  std::set<Cell>::const_iterator elemIt, elemEnd;
  std::vector< std::pair<int,int> > swapVec;

  CellsByDimIterator cellsIt( *this, 3 );
  while ( const std::set<Cell > * elems = cellsIt.nextType() )
    {
      TCellType cellType = cellsIt.type();
      elemIt = elems->begin(), elemEnd = elems->end();
      int nbBottomNodes = 0;
      switch ( cellType )
        {
        case NORM_TETRA4:
        case NORM_TETRA10:
        case NORM_PENTA6:
        case NORM_PENTA15:
          nbBottomNodes = 3; break;
        case NORM_PYRA5:
        case NORM_PYRA13:
        case NORM_HEXA8:
        case NORM_HEXA20:
          nbBottomNodes = 4; break;
        default: continue;
        }
      getReverseVector( cellType, swapVec );

      for ( ; elemIt != elemEnd; elemIt++ )
        {
          // find a normal to the bottom face
          const double* n[4];
          n[0] = nodeCoords( elemIt->_nodes[0]); // 3 bottom nodes
          n[1] = nodeCoords( elemIt->_nodes[1]);
          n[2] = nodeCoords( elemIt->_nodes[2]);
          n[3] = nodeCoords( elemIt->_nodes[nbBottomNodes]); // a top node
          double vec01[3]; // vector n[0]-n[1]
          vec01[0] = n[1][0] - n[0][0];
          vec01[1] = n[1][1] - n[0][1];
          vec01[2] = n[1][2] - n[0][2];
          double vec02 [3]; // vector n[0]-n[2]
          vec02[0] = n[2][0] - n[0][0];
          vec02[1] = n[2][1] - n[0][1];
          vec02[2] = n[2][2] - n[0][2];
          double normal [3]; // vec01 ^ vec02
          normal[0] = vec01[1] * vec02[2] - vec01[2] * vec02[1];
          normal[1] = vec01[2] * vec02[0] - vec01[0] * vec02[2];
          normal[2] = vec01[0] * vec02[1] - vec01[1] * vec02[0];
          // check if the 102 angle is convex
          if ( nbBottomNodes > 3 )
            {
              const double* n3 = nodeCoords( elemIt->_nodes[nbBottomNodes-1] );// last bottom node
              double vec03 [3];  // vector n[0]-n3
              vec03[0] = n3[0] - n[0][0];
              vec03[1] = n3[1] - n[0][1];
              vec03[2] = n3[2] - n[0][2];
              if ( fabs( normal[0]+normal[1]+normal[2] ) <= std::numeric_limits<double>::max() ) // vec01 || vec02
                {
                  normal[0] = vec01[1] * vec03[2] - vec01[2] * vec03[1]; // vec01 ^ vec03
                  normal[1] = vec01[2] * vec03[0] - vec01[0] * vec03[2];
                  normal[2] = vec01[0] * vec03[1] - vec01[1] * vec03[0];
                }
              else
                {
                  double vec [3]; // normal ^ vec01
                  vec[0] = normal[1] * vec01[2] - normal[2] * vec01[1];
                  vec[1] = normal[2] * vec01[0] - normal[0] * vec01[2];
                  vec[2] = normal[0] * vec01[1] - normal[1] * vec01[0];
                  double dot2 = vec[0]*vec03[0] + vec[1]*vec03[1] + vec[2]*vec03[2]; // vec*vec03
                  if ( dot2 < 0 ) // concave -> reverse normal
                    {
                      normal[0] *= -1;
                      normal[1] *= -1;
                      normal[2] *= -1;
                    }
                }
            }
          // direction from top to bottom
          double tbDir[3];
          tbDir[0] = n[0][0] - n[3][0];
          tbDir[1] = n[0][1] - n[3][1];
          tbDir[2] = n[0][2] - n[3][2];

          // compare 2 directions: normal and top-bottom
          double dot = normal[0]*tbDir[0] + normal[1]*tbDir[1] + normal[2]*tbDir[2];
          if ( dot < 0. ) // need reverse
            reverse( *elemIt, swapVec );

        } // loop on volumes of one geometry
    } // loop on 3D geometry types

}

//================================================================================
/*!
 * \brief Assign new IDs to nodes by skipping not used nodes and return their number
 */
//================================================================================

int NodeContainer::numberNodes()
{
  int id = 1;
  for ( size_t i = 0; i < _nodes.size(); ++i )
    for ( size_t j = 0; j < _nodes[i].size(); ++j )
      if ( _nodes[i][j].isUsed() )
        _nodes[i][j]._number = id++;
  return id-1;
}


//================================================================================
/*!
 * \brief Assign new IDs to elements
 */
//================================================================================

void IntermediateMED::numberElements()
{
  std::set<Cell>::const_iterator elemIt, elemEnd;

  // numbering _cells of type NORM_POINT1 by node number
  {
    const std::set<Cell>& points = _cellsByType[ INTERP_KERNEL::NORM_POINT1 ];
    elemIt = points.begin(), elemEnd = points.end();
    for ( ; elemIt != elemEnd; ++elemIt )
      elemIt->_number = elemIt->_nodes[0]->_number;
  }

  // numbering 1D-3D _cells
  for ( int dim = 1; dim <= 3; ++dim )
    {
      // check if re-numeration is needed (to try to keep elem oreder as in sauve file )
      bool ok = true, renumEntity = false;
      CellsByDimIterator cellsIt( *this, dim );
      int prevNbElems = 0;
      while ( const std::set<Cell> * typeCells = cellsIt.nextType() )
        {
          TID minNumber = std::numeric_limits<TID>::max(), maxNumber = 0;
          for ( elemIt = typeCells->begin(), elemEnd = typeCells->end(); elemIt!=elemEnd; ++elemIt)
            {
              if ( elemIt->_number < minNumber ) minNumber = elemIt->_number;
              if ( elemIt->_number > maxNumber ) maxNumber = elemIt->_number;
            }
          TID typeSize = typeCells->size();
          if ( typeSize != maxNumber - minNumber + 1 )
            ok = false;
          if ( prevNbElems+1 != (int)minNumber )
            ok = false;
          if ( prevNbElems != 0 && minNumber == 1 )
            renumEntity = true;

          prevNbElems += typeSize;
        }

      if ( ok && renumEntity ) // each geom type was numerated separately
        {
          cellsIt.init( dim );
          prevNbElems = cellsIt.nextType()->size(); // no need to renumber the first type
          while ( const std::set<Cell> * typeCells = cellsIt.nextType() )
            {
              for ( elemIt = typeCells->begin(), elemEnd = typeCells->end(); elemIt!=elemEnd; ++elemIt)
                elemIt->_number += prevNbElems;
              prevNbElems += typeCells->size();
            }
        }
      if ( !ok )
        {
          int cellID=1;
          cellsIt.init( dim );
          while ( const std::set<Cell> * typeCells = cellsIt.nextType() )
            for ( elemIt = typeCells->begin(), elemEnd = typeCells->end(); elemIt!=elemEnd; ++elemIt)
              elemIt->_number = cellID++;
        }
    }
}

//================================================================================
/*!
 * \brief Creates coord array
 */
//================================================================================

MEDCoupling::DataArrayDouble * IntermediateMED::getCoords()
{
  DataArrayDouble* coordArray = DataArrayDouble::New();
  coordArray->alloc( _nbNodes, _spaceDim );
  double * coordPrt = coordArray->getPointer();
  for ( int i = 0, nb = _points.size(); i < nb; ++i )
    {
      Node* n = getNode( i+1 );
      if ( n->isUsed() )
        {
          const double* nCoords = nodeCoords( n );
          std::copy( nCoords, nCoords+_spaceDim, coordPrt );
          coordPrt += _spaceDim;
        }
    }
  return coordArray;
}

//================================================================================
/*!
 * \brief Sets connectivity of elements to the mesh
 *  \param mesh - mesh to fill in
 *  \param coords - coordinates that must be shared by all meshes of different dim
 */
//================================================================================

void IntermediateMED::setConnectivity( MEDCoupling::MEDFileUMesh*    mesh,
                                       MEDCoupling::DataArrayDouble* coords )
{
  int meshDim = 0;

  mesh->setCoords( coords );

  std::set<Cell>::const_iterator elemIt, elemEnd;
  for ( int dim = 3; dim > 0; --dim )
    {
      CellsByDimIterator dimCells( *this, dim );

      int nbOfCells = 0;
      while ( const std::set<Cell > * cells = dimCells.nextType() )
        nbOfCells += cells->size();
      if ( nbOfCells == 0 )
        continue;

      if ( !meshDim ) meshDim = dim;

      MEDCouplingUMesh* dimMesh = MEDCouplingUMesh::New();
      dimMesh->setCoords( coords );
      dimMesh->setMeshDimension( dim );
      dimMesh->allocateCells( nbOfCells );

      int prevNbCells = 0;
      dimCells.init( dim );
      while ( const std::set<Cell > * cells = dimCells.nextType() )
        {
          // fill connectivity array to take into account order of elements in the sauv file
          const int nbCellNodes = cells->begin()->_nodes.size();
          std::vector< TID > connectivity( cells->size() * nbCellNodes );
          int * nodalConnOfCell;
          for ( elemIt = cells->begin(), elemEnd = cells->end(); elemIt != elemEnd; ++elemIt )
            {
              const Cell& cell = *elemIt;
              const int index = cell._number - 1 - prevNbCells;
              nodalConnOfCell = &connectivity[ index * nbCellNodes ];
              if ( cell._reverse )
                for ( int i = nbCellNodes-1; i >= 0; --i )
                  *nodalConnOfCell++ = cell._nodes[i]->_number - 1;
              else
                for ( int i = 0; i < nbCellNodes; ++i )
                  *nodalConnOfCell++ = cell._nodes[i]->_number - 1;
            }
          prevNbCells += cells->size();

          // fill dimMesh
          TCellType cellType = dimCells.type();
          nodalConnOfCell = &connectivity[0];
          for ( size_t i = 0; i < cells->size(); ++i, nodalConnOfCell += nbCellNodes )
            dimMesh->insertNextCell( cellType, nbCellNodes, nodalConnOfCell );
        }
      dimMesh->finishInsertingCells();
      mesh->setMeshAtLevel( dim - meshDim, dimMesh );
      dimMesh->decrRef();
    }
}

//================================================================================
/*!
 * \brief Fill in the mesh with groups
 *  \param mesh - mesh to fill in
 */
//================================================================================

void IntermediateMED::setGroups( MEDCoupling::MEDFileUMesh* mesh )
{
  bool isMeshNameSet = false;
  const int meshDim = mesh->getMeshDimension();
  for ( int dim = 0; dim <= meshDim; ++dim )
    {
      const int meshDimRelToMaxExt = ( dim == 0 ? 1 : dim - meshDim );

      std::vector<const DataArrayInt *> medGroups;
      std::vector<MCAuto<DataArrayInt> > refGroups;
      for ( size_t i = 0; i < _groups.size(); ++i )
        {
          Group& grp = _groups[i];
          if ( (int)getDim( &grp ) != dim &&
               grp._groups.empty() ) // to allow groups on diff dims
            continue;
          // convert only named groups or field supports
          if ( grp.empty() || (grp._name.empty() && !grp._isProfile ))
            continue;
          //if ( grp._medGroup ) continue; // already converted

          // sort cells by ID and remember their initial order in the group
          TCellToOrderMap cell2order;
          unsigned orderInGroup = 0;
          std::vector< Group* > groupVec;
          if ( grp._groups.empty() ) groupVec.push_back( & grp );
          else                       groupVec = grp._groups;
          for ( size_t iG = 0; iG < groupVec.size(); ++iG )
            {
              Group* aG = groupVec[ iG ];
              if ( (int)getDim( aG ) != dim )
                continue;
              for ( size_t iC = 0; iC < aG->_cells.size(); ++iC )
                cell2order.insert( cell2order.end(), std::make_pair( aG->_cells[iC], orderInGroup++ ));
            }
          if ( cell2order.empty() )
            continue;
          bool isSelfIntersect = ( orderInGroup != cell2order.size() );
          if ( isSelfIntersect ) // self intersecting group
            {
              std::ostringstream msg;
              msg << "Self intersecting sub-mesh: id = " << i+1
                  << ", name = |" << grp._name << "|" << std::endl
                  << " nb unique elements = " << cell2order.size() << std::endl
                  << " total nb elements  = " << orderInGroup;
              if ( grp._isProfile )
                {
                  THROW_IK_EXCEPTION( msg.str() );
                }
              else
                {
                  std::cout << msg.str() << std::endl;
                }
            }
          // create a med group
          grp._medGroup = DataArrayInt::New();
          grp._medGroup->setName( grp._name.c_str() );
          grp._medGroup->alloc( cell2order.size(), /*nbOfCompo=*/1 );
          int * idsPtr = grp._medGroup->getPointer();
          TCellToOrderMap::iterator cell2orderIt, cell2orderEnd = cell2order.end();
          for ( cell2orderIt = cell2order.begin(); cell2orderIt != cell2orderEnd; ++cell2orderIt )
            *idsPtr++ = (*cell2orderIt).first->_number - 1;

          // try to set the mesh name
          if ( !isMeshNameSet &&
               dim == meshDim &&
               !grp._name.empty() &&
               grp.size() == mesh->getSizeAtLevel( meshDimRelToMaxExt ))
            {
              mesh->setName( grp._name.c_str() );
              isMeshNameSet = true;
            }
          if ( !grp._name.empty() )
            {
              medGroups.push_back( grp._medGroup );
            }
          // set relocation table
          setRelocationTable( &grp, cell2order );

          // Issue 0021311. Use case: a gibi group has references (recorded in pile 1)
          // and several names (pile 27) refer (pile 10) to this group.
          // We create a copy of this group per each named reference
          std::set<std::string> uniqueNames;
          uniqueNames.insert( grp._name );
          for ( unsigned iRef = 0 ; iRef < grp._refNames.size(); ++iRef )
            if ( !grp._refNames[ iRef ].empty() &&
                 uniqueNames.insert( grp._refNames[ iRef ]).second ) // for name uniqueness (23155)
              {
                refGroups.push_back( grp._medGroup->deepCopy() );
                refGroups.back()->setName( grp._refNames[ iRef ].c_str() );
                medGroups.push_back( refGroups.back() );
              }
        }
      mesh->setGroupsAtLevel( meshDimRelToMaxExt, medGroups );
    }
}

//================================================================================
/*!
 * \brief Return true if the group is on all elements and return its relative dimension
 */
//================================================================================

bool IntermediateMED::isOnAll( const Group* grp, int & dimRel ) const
{
  int dim = getDim( grp );

  int nbElems = 0;
  if ( dim == 0 )
    {
      nbElems = _nbNodes;
      dimRel  = 0;
    }
  else
    {
      CellsByDimIterator dimCells( *this, dim );
      while ( const std::set<Cell > * cells = dimCells.nextType() )
        nbElems += cells->size();

      int meshDim = 3;
      for ( ; meshDim > 0; --meshDim )
        {
          dimCells.init( meshDim );
          if ( dimCells.nextType() )
            break;
        }
      dimRel = dim - meshDim;
    }

  bool onAll = ( nbElems == grp->size() );
  return onAll;
}

//================================================================================
/*!
 * \brief Makes fields from own data
 */
//================================================================================

MEDCoupling::MEDFileFields * IntermediateMED::makeMEDFileFields(MEDCoupling::MEDFileUMesh* mesh)
{
  if ( _nodeFields.empty() && _cellFields.empty() ) return 0;

  // set long names
  std::set< std::string > usedFieldNames;
  setFieldLongNames(usedFieldNames);

  MEDFileFields* fields = MEDFileFields::New();

  for ( size_t i = 0; i < _nodeFields.size(); ++i )
    setFields( _nodeFields[i], fields, mesh, i+1, usedFieldNames );

  for ( size_t i = 0; i < _cellFields.size(); ++i )
    setFields( _cellFields[i], fields, mesh, i+1, usedFieldNames );

  return fields;
}

//================================================================================
/*!
 * \brief Make med fields from a SauvUtilities::DoubleField
 */
//================================================================================

void IntermediateMED::setFields( SauvUtilities::DoubleField* fld,
                                 MEDCoupling::MEDFileFields*  medFields,
                                 MEDCoupling::MEDFileUMesh*   mesh,
                                 const TID                   castemID,
                                 std::set< std::string >&    usedFieldNames)
{
  bool sameNbGauss = true;
  if ( !fld || !fld->isMedCompatible( sameNbGauss )) return;

  if ( !sameNbGauss )
    fld->splitSubWithDiffNbGauss();

  // if ( !fld->hasCommonSupport() ):
  //     each sub makes MEDFileFieldMultiTS
  // else:
  //     unite several subs into a MEDCouplingFieldDouble

  const bool uniteSubs = fld->hasCommonSupport() && sameNbGauss;
  if ( !uniteSubs )
    std::cout << "Castem field #" << castemID << " <" << fld->_name
              << "> is incompatible with MED format, so we split it into several fields:" << std::endl;

  for ( size_t iSub = 0; iSub < fld->_sub.size(); )
    {
      // set field name
      if ( !uniteSubs || fld->_name.empty() )
        makeFieldNewName( usedFieldNames, fld );

      // allocate values
      DataArrayDouble * values = DataArrayDouble::New();
      values->alloc( fld->getNbTuples(iSub), fld->_sub[iSub].nbComponents() );

      // set values
      double * valPtr = values->getPointer();
      if ( uniteSubs )
        {
          int nbElems = fld->_group->size();
          for ( int elemShift = 0; elemShift < nbElems && iSub < fld->_sub.size(); )
            elemShift += fld->setValues( valPtr, iSub++, elemShift );
          setTS( fld, values, medFields, mesh );
        }
      else
        {
          fld->setValues( valPtr, iSub );
          setTS( fld, values, medFields, mesh, iSub++ );

          std::cout << fld->_name << " with compoments";
          for ( size_t i = 0; i < (size_t)fld->_sub[iSub-1].nbComponents(); ++i )
            std::cout << " " << fld->_sub[iSub-1]._comp_names[ i ];
          std::cout << std::endl;
        }
    }
}

//================================================================================
/*!
 * \brief Store value array of a field into med fields
 */
//================================================================================

void IntermediateMED::setTS( SauvUtilities::DoubleField*  fld,
                             MEDCoupling::DataArrayDouble* values,
                             MEDCoupling::MEDFileFields*   medFields,
                             MEDCoupling::MEDFileUMesh*    mesh,
                             const int                    iSub)
{
  // treat a field support
  const Group* support = fld->getSupport( iSub );
  int dimRel;
  const bool onAll = isOnAll( support, dimRel );
  if ( !onAll && support->_name.empty() )
    {
      const_cast<Group*>(support)->_name += "PFL_" + fld->_name;
      support->_medGroup->setName( support->_name.c_str() );
    }

  // make and fill a time-stamp

  MEDCouplingFieldDouble * timeStamp = MEDCouplingFieldDouble::New( fld->getMedType( iSub ),
                                                                    fld->getMedTimeDisc() );
  timeStamp->setName( fld->_name.c_str() );
  timeStamp->setDescription( fld->_description.c_str() );
  // set the mesh
  if ( onAll )
    {
      MCAuto
        < MEDCouplingUMesh > dimMesh = mesh->getMeshAtLevel( dimRel );
      timeStamp->setMesh( dimMesh );
    }
  else if ( timeStamp->getTypeOfField() == MEDCoupling::ON_NODES )
    {
      DataArrayDouble * coo = mesh->getCoords();
      MCAuto
        <DataArrayDouble> subCoo = coo->selectByTupleId(support->_medGroup->begin(),
                                                        support->_medGroup->end());
      MCAuto< MEDCouplingUMesh > nodeSubMesh =
        MEDCouplingUMesh::Build0DMeshFromCoords( subCoo );
      timeStamp->setMesh( nodeSubMesh );
    }
  else
    {
      MCAuto
        < MEDCouplingUMesh > dimMesh = mesh->getMeshAtLevel( dimRel );
      MCAuto
        <MEDCouplingMesh> subMesh = dimMesh->buildPart(support->_medGroup->begin(),
                                                       support->_medGroup->end());
      timeStamp->setMesh( subMesh);
    }
  // set values
  for ( size_t i = 0; i < (size_t)fld->_sub[iSub].nbComponents(); ++i )
    values->setInfoOnComponent( i, fld->_sub[iSub]._comp_names[ i ].c_str() );
  timeStamp->setArray( values );
  values->decrRef();
  // set gauss points
  if ( timeStamp->getTypeOfField() == MEDCoupling::ON_GAUSS_PT )
    {
      TGaussDef gaussDef( fld->_sub[iSub]._support->_cellType,
                          fld->_sub[iSub].nbGauss() );
      timeStamp->setGaussLocalizationOnType( fld->_sub[iSub]._support->_cellType,
                                             gaussDef.myRefCoords,
                                             gaussDef.myCoords,
                                             gaussDef.myWeights );
    }
  // get a field to add the time-stamp
  bool isNewMedField = false;
  if ( !fld->_curMedField || fld->_name != fld->_curMedField->getName() )
    {
      fld->_curMedField = MEDFileFieldMultiTS::New();
      isNewMedField = true;
    }

  // set an order
  const int nbTS = fld->_curMedField->getNumberOfTS();
  if ( nbTS > 0 )
    timeStamp->setOrder( nbTS );

  // add the time-stamp
  timeStamp->checkConsistencyLight();
  if ( onAll )
    fld->_curMedField->appendFieldNoProfileSBT( timeStamp );
  else
    fld->_curMedField->appendFieldProfile( timeStamp, mesh, dimRel, support->_medGroup );
  timeStamp->decrRef();

  if ( isNewMedField ) // timeStamp must be added before this
    {
      medFields->pushField( fld->_curMedField );
    }
}

//================================================================================
/*!
 * \brief Make a new unique name for a field
 */
//================================================================================

void IntermediateMED::makeFieldNewName(std::set< std::string >&    usedNames,
                                       SauvUtilities::DoubleField* fld )
{
  std::string base = fld->_name;
  if ( base.empty() )
    {
      base = "F_";
    }
  else
    {
      std::string::size_type pos = base.rfind('_');
      if ( pos != std::string::npos )
        base = base.substr( 0, pos+1 );
      else
        base += '_';
    }

  int i = 1;
  do
    {
      fld->_name = base + SauvUtilities::toString( i++ );
    }
  while( !usedNames.insert( fld->_name ).second );
}

//================================================================================
/*!
 * \brief Split sub-components with different nb of gauss points into several sub-components
 *  \param [in,out] fld - a field to split if necessary
 */
//================================================================================

void DoubleField::splitSubWithDiffNbGauss()
{
  for ( size_t iSub = 0; iSub < _sub.size(); ++iSub )
    {
      if ( _sub[iSub].isSameNbGauss() ) continue;

      _sub.insert( _sub.begin() + iSub + 1, 1, _Sub_data() );
      _Sub_data & subToSplit = _sub[iSub];
      _Sub_data & subNew     = _sub[iSub+1];
      size_t iDiff = 1;
      while ( subToSplit._nb_gauss[ 0 ] == subToSplit._nb_gauss[ iDiff ] )
        ++iDiff;
      subNew._support = subToSplit._support;
      subNew._comp_names.assign( subToSplit._comp_names.begin() + iDiff,
                                 subToSplit._comp_names.end() );
      subNew._nb_gauss.assign  ( subToSplit._nb_gauss.begin() + iDiff,
                                 subToSplit._nb_gauss.end() );
      subToSplit._comp_names.resize( iDiff );
      subToSplit._nb_gauss.resize  ( iDiff );
    }
}

//================================================================================
/*!
 * \brief Return a vector ready to fill in
 */
//================================================================================

std::vector< double >& DoubleField::addComponent( int nb_values )
{
  _comp_values.push_back( std::vector< double >() );
  std::vector< double >& res = _comp_values.back();
  res.resize( nb_values );
  return res;
}

DoubleField::~DoubleField()
{
  if(_curMedField)
    _curMedField->decrRef();
}

//================================================================================
/*!
 * \brief Returns a supporting group
 */
//================================================================================

const Group* DoubleField::getSupport( const int iSub ) const
{
  return _group ? _group : _sub[iSub]._support;
}

//================================================================================
/*!
 * \brief Return true if each sub-component is a time stamp
 */
//================================================================================

bool DoubleField::isMultiTimeStamps() const
{
  if ( _sub.size() < 2 )
    return false;
  bool sameSupports = true;
  Group* grpp1 = _sub[0]._support;// grpp NOT grp because XDR under Windows defines grp...
  for ( size_t i = 1; i < _sub.size() && sameSupports; ++i )
    sameSupports = ( grpp1 == _sub[i]._support );

  return sameSupports;
}

//================================================================================
/*!
 * \brief True if the field can be converted into the med field
 */
//================================================================================

bool DoubleField::isMedCompatible(bool& sameNbGauss) const
{
  for ( size_t iSub = 0; iSub < _sub.size(); ++iSub )
    {
      if ( !getSupport(iSub) || !getSupport(iSub)->_medGroup )
        THROW_IK_EXCEPTION("SauvReader INTERNAL ERROR: NULL field support");

      sameNbGauss = true;
      if ( !_sub[iSub].isSameNbGauss() )
        {
          std::cout << "Field <" << _name << "> : different nb of gauss points in components" << std::endl;
          sameNbGauss = false;
          //return false;
        }
    }
  return true;
}

//================================================================================
/*!
 * \brief return true if all sub-components has same components and same nbGauss
 */
//================================================================================

bool DoubleField::hasSameComponentsBySupport() const
{
  std::vector< _Sub_data >::const_iterator sub_data = _sub.begin();
  const _Sub_data& first_sub_data = *sub_data;
  for ( ++sub_data ; sub_data != _sub.end(); ++sub_data )
    {
      if ( first_sub_data._comp_names != sub_data->_comp_names )
        return false; // diff names of components

      if ( first_sub_data._nb_gauss != sub_data->_nb_gauss &&
           first_sub_data._support->_cellType == sub_data->_support->_cellType)
        return false; // diff nb of gauss points on same cell type
    }
  return true;
}

//================================================================================
/*!
 * \brief Return type of MEDCouplingFieldDouble
 */
//================================================================================

MEDCoupling::TypeOfField DoubleField::getMedType( const int iSub ) const
{
  using namespace INTERP_KERNEL;

  const Group* grp = hasCommonSupport() ? _group : _sub[iSub]._support;
  if ( _sub[iSub].nbGauss() > 1 )
    {
      const CellModel& cm = CellModel::GetCellModel( _sub[iSub]._support->_cellType );
      return (int) cm.getNumberOfNodes() == _sub[iSub].nbGauss() ? ON_GAUSS_NE : ON_GAUSS_PT;
    }
  else
    {
      return getDim( grp ) == 0 ? ON_NODES : ON_CELLS;
    }
}

//================================================================================
/*!
 * \brief Return TypeOfTimeDiscretization
 */
//================================================================================

MEDCoupling::TypeOfTimeDiscretization DoubleField::getMedTimeDisc() const
{
  return ONE_TIME;
  // NO_TIME = 4,
  // ONE_TIME = 5,
  // LINEAR_TIME = 6,
  // CONST_ON_TIME_INTERVAL = 7
}

//================================================================================
/*!
 * \brief Return nb tuples to be used to allocate DataArrayDouble
 */
//================================================================================

int DoubleField::getNbTuples( const int iSub ) const
{
  int nb = 0;
  if ( hasCommonSupport() && !_group->_groups.empty() )
    for ( size_t i = 0; i < _group->_groups.size(); ++i )
      nb += _sub[i].nbGauss() * _sub[i]._support->size();
  else
    nb = _sub[iSub].nbGauss() * getSupport(iSub)->size();
  return nb;
}

//================================================================================
/*!
 * \brief Store values of a sub-component and return nb of elements in the iSub
 */
//================================================================================

int DoubleField::setValues( double * valPtr, const int iSub, const int elemShift ) const
{
  // find values for iSub
  int iComp = 0;
  for ( int iS = 0; iS < iSub; ++iS )
    iComp += _sub[iS].nbComponents();
  const std::vector< double > * compValues = &_comp_values[ iComp ];

  // Set values

  const std::vector< unsigned >& relocTable = getSupport( iSub )->_relocTable;

  const int nbElems      = _sub[iSub]._support->size();
  const int nbGauss      = _sub[iSub].nbGauss();
  const int nbComponents = _sub[iSub].nbComponents();
  const int nbValsByElem = nbComponents * nbGauss;

  // check nb values
  int nbVals = 0;
  for ( iComp = 0; iComp < nbComponents; ++iComp )
    nbVals += compValues[iComp].size();
  const bool isConstField = ( nbVals == nbComponents ); // one value per component (issue 22321)
  if ( !isConstField && nbVals != nbElems * nbValsByElem )
    THROW_IK_EXCEPTION("SauvMedConvertor.cxx: support size mismatches field size");

  // compute nb values in previous subs
  int valsShift = 0;
  for ( int iS = iSub-1, shift = elemShift; shift > 0; --iS)
    {
      int nbE = _sub[iS]._support->size();
      shift -= nbE;
      valsShift += nbE * _sub[iS].nbComponents() * _sub[iS].nbGauss();
    }

  if ( isConstField )
    for ( int iE = 0; iE < nbElems; ++iE )
      {
        int iMed = valsShift + nbValsByElem * ( relocTable.empty() ? iE : relocTable[iE+elemShift]-elemShift );
        for ( iComp = 0; iComp < nbComponents; ++iComp )
          valPtr[ iMed + iComp ] = compValues[iComp][ 0 ];
      }
  else
    for ( int iE = 0; iE < nbElems; ++iE )
      {
        int iMed = valsShift + nbValsByElem * ( relocTable.empty() ? iE : relocTable[iE+elemShift]-elemShift );
        for ( iComp = 0; iComp < nbComponents; ++iComp )
          for ( int iG = 0; iG < nbGauss; ++iG )
            valPtr[ iMed + iG * nbComponents + iComp ] = compValues[iComp][ iE * nbGauss + iG ];
      }
  return nbElems;
}

//================================================================================
/*!
 * \brief Destructor of IntermediateMED
 */
//================================================================================

IntermediateMED::~IntermediateMED()
{
  for ( size_t i = 0; i < _nodeFields.size(); ++i )
    if ( _nodeFields[i] )
      delete _nodeFields[i];
  _nodeFields.clear();

  for ( size_t i = 0; i < _cellFields.size(); ++i )
    if ( _cellFields[i] )
      delete _cellFields[i];
  _cellFields.clear();

  for ( size_t i = 0; i < _groups.size(); ++i )
    if ( _groups[i]._medGroup )
      _groups[i]._medGroup->decrRef();
}

//================================================================================
/*!
 * \brief CellsByDimIterator constructor
 */
CellsByDimIterator::CellsByDimIterator( const IntermediateMED & medi, int dimm)
{
  myImed = & medi;
  init( dimm );
}
/*!
 * \brief Initialize iteration on cells of given dimention
 */
void CellsByDimIterator::init(const int  dimm)
{
  myCurType = -1;
  myTypeEnd = INTERP_KERNEL::NORM_HEXA20 + 1;
  myDim = dimm;
}
/*!
 * \brief return next set of Cell's of required dimension
 */
const std::set< Cell > * CellsByDimIterator::nextType()
{
  while ( ++myCurType < myTypeEnd )
    if ( !myImed->_cellsByType[myCurType].empty() && ( myDim < 0 || dim(false) == myDim ))
      return & myImed->_cellsByType[myCurType];
  return 0;
}
/*!
 * \brief return dimension of cells returned by the last or further next()
 */
int CellsByDimIterator::dim(const bool last) const
{
  int typp = myCurType;
  if ( !last )
    while ( typp < myTypeEnd && myImed->_cellsByType[typp].empty() )
      ++typp;
  return typp < myTypeEnd ? getDimension( TCellType( typp )) : 4;
}
// END CellsByDimIterator ========================================================

