// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#ifdef WNT
#include <io.h>
#endif

#ifndef WNT
#define HAS_XDR
#endif

#ifdef HAS_XDR
#include <rpc/xdr.h>
#endif

using namespace SauvUtilities;
using namespace ParaMEDMEM;
using namespace std;

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
                                const Cell &                     aCell )
  {
    if ( const int * conn = getGibi2MedQuadraticInterlace( type ))
      {
        Cell* ma = (Cell*) & aCell;
        //cout << "###### BEFORE ConvertQuadratic() " << *ma << endl;
        vector< Node* > new_nodes( ma->_nodes.size() );
        for ( size_t i = 0; i < new_nodes.size(); ++i )
          new_nodes[ i ] = ma->_nodes[ conn[ i ]];
        ma->_nodes.swap( new_nodes );
        //cout << "###### AFTER ConvertQuadratic() " << *ma << endl;
      }
  }

  //================================================================================
  /*!
   * \brief Returns a vector of pairs of node indices to inverse a med volume element
   */
  //================================================================================

  void getReverseVector (const INTERP_KERNEL::NormalizedCellType type,
                         vector<pair<int,int> > &                swapVec )
  {
    swapVec.clear();

    switch ( type )
      {
      case NORM_TETRA4:
        swapVec.resize(1);
        swapVec[0] = make_pair( 1, 2 );
        break;
      case NORM_PYRA5:
        swapVec.resize(1);
        swapVec[0] = make_pair( 1, 3 );
        break;
      case NORM_PENTA6:
        swapVec.resize(2);
        swapVec[0] = make_pair( 1, 2 );
        swapVec[1] = make_pair( 4, 5 );
        break;
      case NORM_HEXA8:
        swapVec.resize(2);
        swapVec[0] = make_pair( 1, 3 );
        swapVec[1] = make_pair( 5, 7 );
        break;
      case NORM_TETRA10:
        swapVec.resize(3);
        swapVec[0] = make_pair( 1, 2 );
        swapVec[1] = make_pair( 4, 6 );
        swapVec[2] = make_pair( 8, 9 );
        break;
      case NORM_PYRA13:
        swapVec.resize(4);
        swapVec[0] = make_pair( 1, 3 );
        swapVec[1] = make_pair( 5, 8 );
        swapVec[2] = make_pair( 6, 7 );
        swapVec[3] = make_pair( 10, 12 );
        break;
      case NORM_PENTA15:
        swapVec.resize(4);
        swapVec[0] = make_pair( 1, 2 );
        swapVec[1] = make_pair( 4, 5 );
        swapVec[2] = make_pair( 6, 8 );
        swapVec[3] = make_pair( 9, 11 );
        break;
      case NORM_HEXA20:
        swapVec.resize(7);
        swapVec[0] = make_pair( 1, 3 );
        swapVec[1] = make_pair( 5, 7 );
        swapVec[2] = make_pair( 8, 11 );
        swapVec[3] = make_pair( 9, 10 );
        swapVec[4] = make_pair( 12, 15 );
        swapVec[5] = make_pair( 13, 14 );
        swapVec[6] = make_pair( 17, 19 );
        break;
        //   case NORM_SEG3: no need to reverse edges
        //     swapVec.resize(1);
        //     swapVec[0] = make_pair( 1, 2 );
        //     break;
      case NORM_TRI6:
        swapVec.resize(2);
        swapVec[0] = make_pair( 1, 2 );
        swapVec[1] = make_pair( 3, 5 );
        break;
      case NORM_QUAD8:
        swapVec.resize(3);
        swapVec[0] = make_pair( 1, 3 );
        swapVec[1] = make_pair( 4, 7 );
        swapVec[2] = make_pair( 5, 6 );
        break;
      default:;
      }
  }

  //================================================================================
  /*!
   * \brief Inverses element orientation using vector of indices to swap
   */
  //================================================================================

  inline void reverse(const Cell & aCell, const vector<pair<int,int> > & swapVec )
  {
    Cell* ma = (Cell*) & aCell;
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
   * \brief Comparator of cells by number used for ordering cells thinin a med group
   */
  struct TCellByIDCompare
  {
    bool operator () (const Cell* i1, const Cell* i2)
    {
      return i1->_number < i2->_number;
    }
  };
  typedef map< const Cell*, unsigned, TCellByIDCompare > TCellToOrderMap;

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
 * \brief Returns interlace array to transform a quadratic GIBI element to a MED one
 */
//================================================================================

const int * SauvUtilities::getGibi2MedQuadraticInterlace( INTERP_KERNEL::NormalizedCellType type )
{
  static vector<const int*> conn;
  static const int hexa20 [] = {0,6,4,2, 12,18,16,14, 7,5,3,1, 19,17,15,13, 8,11,10,9};
  static const int penta15[] = {0,2,4, 9,11,13, 1,3,5, 10,12,14, 6,7,3};
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
    return make_pair( _nodes[i2]->_number, _nodes[i]->_number );
  else
    return make_pair( _nodes[i]->_number, _nodes[i2]->_number );
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
  int size = 0;
  if ( !_relocTable.empty() )
    size =  _relocTable.size();
  else if ( _medGroup )
    size = _medGroup->getNumberOfTuples();
  else if ( !_cells.empty() )
    size = _cells.size();
  else
    for ( size_t i = 0; i < _groups.size(); ++i )
      size += _groups[i]->size();
  return size;
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
#ifdef WNT
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
  if ( _curPos[posE] != 'E' && _curPos[posE] != 'e' )
    {
      std::string aStr (_curPos);
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

string ASCIIReader::getName() const
{
  int len = _width;
  while (( _curPos[len-1] == ' ' || _curPos[len-1] == 0) && len > 0 )
    len--;
  return string( _curPos, len );
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
  if ((_xdrs_file = ::fopen(_fileName.c_str(), "r")))
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
      cout << "_iRead, _nbToRead : " << _iRead << " " << _nbToRead << endl;
      cout << "Unfinished iteration before new one !" << endl;
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
  return string( s, len );
}

//================================================================================
/*!
 * \brief Throw an exception if not all needed data is present
 */
//================================================================================

void IntermediateMED::checkDataAvailability() const throw(INTERP_KERNEL::Exception)
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
 * \brief Makes ParaMEDMEM::MEDFileData from self
 */
//================================================================================

ParaMEDMEM::MEDFileData* IntermediateMED::convertInMEDFileDS()
{
  MEDCouplingAutoRefCountObjectPtr< MEDFileUMesh >  mesh   = makeMEDFileMesh();
  MEDCouplingAutoRefCountObjectPtr< MEDFileFields > fields = makeMEDFileFields(mesh);

  MEDCouplingAutoRefCountObjectPtr< MEDFileMeshes > meshes = MEDFileMeshes::New();
  MEDCouplingAutoRefCountObjectPtr< MEDFileData >  medData = MEDFileData::New();
  meshes->pushMesh( mesh );
  medData->setMeshes( meshes );
  if ( fields ) medData->setFields( fields );

  medData->incrRef();
  return medData;
}

//================================================================================
/*!
 * \brief Creates ParaMEDMEM::MEDFileUMesh from its data
 */
//================================================================================

ParaMEDMEM::MEDFileUMesh* IntermediateMED::makeMEDFileMesh()
{
  // check if all needed piles are present
  checkDataAvailability();

  // set long names
  setGroupLongNames();

  // fix element orientation
  if ( _spaceDim == 2 )
    orientElements2D();
  else if ( _spaceDim == 3 )
    orientElements3D();

  // process groups
  decreaseHierarchicalDepthOfSubgroups();
  eraseUselessGroups();
  detectMixDimGroups();

  // assign IDs
  _points.numberNodes();
  numberElements();

  // make the med mesh

  MEDFileUMesh* mesh = MEDFileUMesh::New();

  DataArrayDouble *coords = getCoords();
  setConnectivity( mesh, coords );
  setGroups( mesh );

  coords->decrRef();

  if ( !mesh->getName() || strlen( mesh->getName() ) == 0 )
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
  // IMP 0020434: mapping GIBI names to MED names
  // set med names to objects (mesh, fields, support, group or other)

  set<int> treatedGroups;

  list<nameGIBItoMED>::iterator itGIBItoMED = _listGIBItoMED_mail.begin();
  for (; itGIBItoMED != _listGIBItoMED_mail.end(); itGIBItoMED++)
    {
      if ( (int)_groups.size() < itGIBItoMED->gibi_id ) continue;

      SauvUtilities::Group & grp = _groups[itGIBItoMED->gibi_id - 1];

      // if there are several names for grp then the 1st name is the name
      // of grp and the rest ones are names of groups referring grp (issue 0021311)
      const bool isRefName = !treatedGroups.insert( itGIBItoMED->gibi_id ).second;
      if ( !isRefName )
        grp._name = _mapStrings[ itGIBItoMED->med_id ];
      else
        for ( unsigned i = 0; i < grp._refNames.size(); ++i )
          if ( grp._refNames[i].empty() )
            grp._refNames[i] = _mapStrings[ (*itGIBItoMED).med_id ];
    }
}

//================================================================================
/*!
 * \brief Set long names to fields
 */
//================================================================================

void IntermediateMED::setFieldLongNames(set< string >& usedNames)
{
  list<nameGIBItoMED>::iterator itGIBItoMED = _listGIBItoMED_cham.begin();
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
      string medName  = _mapStrings[itGIBItoMED->med_id];
      string gibiName = _mapStrings[itGIBItoMED->gibi_id];

      bool name_found = false;
      for ( int isNodal = 0; isNodal < 2 && !name_found; ++isNodal )
        {
          vector<DoubleField* > & fields = isNodal ? _nodeFields : _cellFields;
          for ( size_t ifi = 0; ifi < fields.size() && !name_found; ifi++)
            {
              if (medName.find( fields[ifi]->_name + "." ) == 0 )
                {
                  vector<DoubleField::_Sub_data>& aSubDs = fields[ifi]->_sub;
                  int nbSub = aSubDs.size();
                  for (int isu = 0; isu < nbSub; isu++)
                    for (int ico = 0; ico < aSubDs[isu].nbComponents(); ico++)
                      {
                        if (aSubDs[isu].compName(ico) == gibiName)
                          {
                            string medNameCompo = medName.substr( fields[ifi]->_name.size() + 1 );
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
    vector< Group* > newSubGroups;
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
      unsigned dim2 = getDim( &grp );
      if ( dim1 != dim2 )
      {
        grp._cells.clear();
        grp._groups.clear();
        if ( !grp._name.empty() )
          cout << "Erase a group with elements of different dim |" << grp._name << "|"<< endl;
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
  set<Cell>::const_iterator elemIt, elemEnd;
  vector< pair<int,int> > swapVec;

  // ------------------------------------
  // fix connectivity of quadratic edges
  // ------------------------------------
  set<Cell>& quadEdges = _cellsByType[ INTERP_KERNEL::NORM_SEG3 ];
  if ( !quadEdges.empty() )
    {
      elemIt = quadEdges.begin(), elemEnd = quadEdges.end();
      for ( ; elemIt != elemEnd; ++elemIt )
        ConvertQuadratic( INTERP_KERNEL::NORM_SEG3, *elemIt );
    }

  CellsByDimIterator faceIt( *this, 2 );
  while ( const set<Cell > * faces = faceIt.nextType() )
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
      int iQuad = isQuadratic ? 2 : 1;
      for ( elemIt = faces->begin(), elemEnd = faces->end(); elemIt != elemEnd; elemIt++ )
        {
          // look for index of the most left node
          int iLeft = 0, iNode, nbNodes = elemIt->_nodes.size() / iQuad;
          double x, minX = nodeCoords( elemIt->_nodes[0] )[0];
          for ( iNode = 1; iNode < nbNodes; ++iNode )
            if (( x = nodeCoords( elemIt->_nodes[ iNode ])[ 0 ]) < minX )
              minX = x, iLeft = iNode;

          // indeces of the nodes neighboring the most left one
          int iPrev = ( iLeft - 1 < 0 ) ? nbNodes - 1 : iLeft - 1;
          int iNext = ( iLeft + 1 == nbNodes ) ? 0 : iLeft + 1;
          // find components of prev-left and left-next vectors
          double xP = nodeCoords( elemIt->_nodes[ iPrev ])[ 0 ];
          double yP = nodeCoords( elemIt->_nodes[ iPrev ])[ 1 ];
          double xN = nodeCoords( elemIt->_nodes[ iNext ])[ 0 ];
          double yN = nodeCoords( elemIt->_nodes[ iNext ])[ 1 ];
          double xL = nodeCoords( elemIt->_nodes[ iLeft ])[ 0 ];
          double yL = nodeCoords( elemIt->_nodes[ iLeft ])[ 1 ];
          double xPL = xL - xP, yPL = yL - yP; // components of prev-left vector
          double xLN = xN - xL, yLN = yN - yL; // components of left-next vector
          // normalise y of the vectors
          double modPL = sqrt ( xPL * xPL + yPL * yPL );
          double modLN = sqrt ( xLN * xLN + yLN * yLN );
          if ( modLN > std::numeric_limits<double>::min() &&
               modPL > std::numeric_limits<double>::min() )
            {
              yPL /= modPL;
              yLN /= modLN;
              // summary direction of neighboring links must be positive
              bool clockwise = ( yPL + yLN > 0 );
              if ( !clockwise )
                reverse( *elemIt, swapVec );
            }
        }
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
  orientFaces3D();

  // -----------------
  // fix connectivity
  // -----------------

  set<Cell>::const_iterator elemIt, elemEnd;
  vector< pair<int,int> > swapVec;

  for ( int dim = 1; dim <= 3; ++dim )
  {
    CellsByDimIterator cellsIt( *this, dim );
    while ( const set<Cell > * elems = cellsIt.nextType() )
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

  orientVolumes();
}

//================================================================================
/*!
 * \brief Orient equally (by setting _reverse flag) all connected faces in 3D space
 */
//================================================================================

void IntermediateMED::orientFaces3D()
{
  // fill map of links and their faces
  set<const Cell*> faces;
  map<const Cell*, Group*> fgm;
  map<Link, list<const Cell*> > linkFacesMap;
  map<Link, list<const Cell*> >::iterator lfIt, lfIt2;

  for (size_t i=0; i!=_groups.size(); ++i)
    {
      Group& grp = _groups[i];
      if ( !grp._cells.empty() && getDimension( grp._cellType ) == 2 )
        for ( size_t j = 0; j < grp._cells.size(); ++j )
          if ( faces.insert( grp._cells[j] ).second )
            {
              for ( size_t k = 0; k < grp._cells[j]->_nodes.size(); ++k )
                linkFacesMap[ grp._cells[j]->link( k ) ].push_back( grp._cells[j] );
              fgm.insert( make_pair( grp._cells[j], &grp ));
            }
    }
  // dump linkFacesMap
  //     for ( lfIt = linkFacesMap.begin(); lfIt!=linkFacesMap.end(); lfIt++) {
  //       cout<< "LINK: " << lfIt->first.first << "-" << lfIt->first.second << endl;
  //       list<const Cell*> & fList = lfIt->second;
  //       list<const Cell*>::iterator fIt = fList.begin();
  //       for ( ; fIt != fList.end(); fIt++ )
  //         cout << "\t" << **fIt << fgm[*fIt]->nom << endl;
  //     }

  // Each oriented link must appear in one face only, else a face is reversed.

  queue<const Cell*> faceQueue; /* the queue contains well oriented faces
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
              list< const Cell* > ml;
              if ( lfIt != linkFacesMap.end() )
                {
                  list<const Cell*> & fList = lfIt->second;
                  list<const Cell*>::iterator fIt = fList.begin();
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
                                  list<const Cell*> & ff = lfIt2->second;
                                  ff.erase( find( ff.begin(), ff.end(), badFace ));
                                  if ( ff.empty() )
                                    linkFacesMap.erase( lfIt2 );
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
                  list<const Cell*> & fList = lfIt->second;
                  list<const Cell*>::iterator fIt = fList.begin();
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
                      list<const Cell*>::iterator i = ml.begin();
                      cout << nbFaceByLink << " faces by 1 link:";
                      for( ; i!= ml.end(); i++ )
                        cout << "in sub-mesh " << fgm[ *i ]->_name << endl << **i;
                    }
                  manifold = false;
                }
            } // loop on links of the being checked face
        } // loop on the face queue
    } // while ( !linkFacesMap.empty() )

  if ( !manifold )
    cout << " -> Non manifold mesh, faces orientation may be incorrect" << endl;
}

//================================================================================
/*!
 * \brief Orient volumes according to MED conventions:
 * normal of a bottom (first) face should be outside
 */
//================================================================================

void IntermediateMED::orientVolumes()
{
  set<Cell>::const_iterator elemIt, elemEnd;
  vector< pair<int,int> > swapVec;

  CellsByDimIterator cellsIt( *this, 3 );
  while ( const set<Cell > * elems = cellsIt.nextType() )
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
              if ( fabs( normal[0]+normal[1]+normal[2] ) <= numeric_limits<double>::max() ) // vec01 || vec02
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
  set<Cell>::const_iterator elemIt, elemEnd;

  // numbering _cells of type NORM_POINT1 by node number
  {
    const set<Cell>& points = _cellsByType[ INTERP_KERNEL::NORM_POINT1 ];
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
      while ( const set<Cell> * typeCells = cellsIt.nextType() )
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
          if ( prevNbElems != 0 ) {
            if ( minNumber == 1 )
              renumEntity = true;
            else if ( prevNbElems+1 != (int)minNumber )
              ok = false;
          }
          prevNbElems += typeSize;
        }

      if ( ok && renumEntity ) // each geom type was numerated separately
        {
          cellsIt.init( dim );
          prevNbElems = cellsIt.nextType()->size(); // no need to renumber the first type
          while ( const set<Cell> * typeCells = cellsIt.nextType() )
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
          while ( const set<Cell> * typeCells = cellsIt.nextType() )
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

ParaMEDMEM::DataArrayDouble * IntermediateMED::getCoords()
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

void IntermediateMED::setConnectivity( ParaMEDMEM::MEDFileUMesh*    mesh,
                                       ParaMEDMEM::DataArrayDouble* coords )
{
  int meshDim = 0;

  mesh->setCoords( coords );

  set<Cell>::const_iterator elemIt, elemEnd;
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
        vector< TID > connectivity( cells->size() * nbCellNodes );
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

void IntermediateMED::setGroups( ParaMEDMEM::MEDFileUMesh* mesh )
{
  const int meshDim = mesh->getMeshDimension();
  for ( int dim = 0; dim <= meshDim; ++dim )
    {
      const int meshDimRelToMaxExt = ( dim == 0 ? 1 : dim - meshDim );

      vector<const DataArrayInt *> medGroups;
      vector<MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > refGroups;
      for ( size_t i = 0; i < _groups.size(); ++i )
        {
          Group& grp = _groups[i];
          if ( (int)getDim( &grp ) != dim )
            continue;
          // convert only named groups or field supports
          if ( grp.empty() || (grp._name.empty() && !grp._isProfile ))
            continue;
          //if ( grp._medGroup ) continue; // already converted

          // sort cells by ID and remember their initial order in the group
          TCellToOrderMap cell2order;
          unsigned orderInGroup = 0;
          vector< Group* > groupVec;
          if ( grp._groups.empty() ) groupVec.push_back( & grp );
          else                       groupVec = grp._groups;
          for ( size_t iG = 0; iG < groupVec.size(); ++iG )
            {
              Group* aG = groupVec[ iG ];
              for ( size_t iC = 0; iC < aG->_cells.size(); ++iC )
                cell2order.insert( cell2order.end(), make_pair( aG->_cells[iC], orderInGroup++ ));
            }
          bool isSelfIntersect = ( orderInGroup != cell2order.size() );
          if ( isSelfIntersect ) // self intersecting group
            {
              ostringstream msg;
              msg << "Self intersecting sub-mesh: id = " << i+1
                  << ", name = |" << grp._name << "|" << endl
                  << " nb unique elements = " << cell2order.size() << endl
                  << " total nb elements  = " << orderInGroup;
              if ( grp._isProfile )
                {
                  THROW_IK_EXCEPTION( msg.str() );
                }
              else
                {
                  cout << msg << endl;
                }
            }
          // create a med group
          grp._medGroup = DataArrayInt::New();
          grp._medGroup->setName( grp._name.c_str() );
          grp._medGroup->alloc( orderInGroup, /*nbOfCompo=*/1 );
          int * idsPrt = grp._medGroup->getPointer();
          TCellToOrderMap::iterator cell2orderIt, cell2orderEnd = cell2order.end();
          for ( cell2orderIt = cell2order.begin(); cell2orderIt != cell2orderEnd; ++cell2orderIt )
            *idsPrt++ = (*cell2orderIt).first->_number - 1;

          // try to set the mesh name
          if ( dim == meshDim &&
               !grp._name.empty() &&
               grp.size() == mesh->getSizeAtLevel( meshDimRelToMaxExt ))
            {
              mesh->setName( grp._name.c_str() );
            }
          else if ( !grp._name.empty() )
            {
              medGroups.push_back( grp._medGroup );
            }
          // set relocation table
          setRelocationTable( &grp, cell2order );

          // Issue 0021311. Use case: a gibi group has references (recorded in pile 1)
          // and several names (pile 27) refer (pile 10) to this group.
          // We create a copy of this group per each named reference
          for ( unsigned iRef = 0 ; iRef < grp._refNames.size(); ++iRef )
            if ( !grp._refNames[ iRef ].empty() )
              {
                refGroups.push_back( grp._medGroup->deepCpy() );
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
  CellsByDimIterator dimCells( *this, dim );
  while ( const set<Cell > * cells = dimCells.nextType() )
    nbElems += cells->size();

  const bool onAll = ( nbElems == grp->size() );

  if ( dim == 0 )
    dimRel = 0;
  else
    {
      int meshDim = 3;
      for ( ; meshDim > 0; --meshDim )
        {
          dimCells.init( meshDim );
          if ( dimCells.nextType() )
            break;
        }
      dimRel = dim - meshDim;
    }
  return onAll;
}

//================================================================================
/*!
 * \brief Makes fields from own data
 */
//================================================================================

ParaMEDMEM::MEDFileFields * IntermediateMED::makeMEDFileFields(ParaMEDMEM::MEDFileUMesh* mesh)
{
  if ( _nodeFields.empty() && _cellFields.empty() ) return 0;

  // set long names
  set< string > usedFieldNames;
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
                                 ParaMEDMEM::MEDFileFields*  medFields,
                                 ParaMEDMEM::MEDFileUMesh*   mesh,
                                 const TID                   castemID,
                                 set< string >&              usedFieldNames)
{
  if ( !fld || !fld->isMedCompatible() ) return;

  // if ( !fld->hasCommonSupport() ):
  //     each sub makes MEDFileFieldMultiTS
  // else:
  //     unite several subs into a MEDCouplingFieldDouble

  const bool uniteSubs = fld->hasCommonSupport();
  if ( !uniteSubs )
    cout << "Castem field #" << castemID << " " << fld->_name
         << " is incompatible with MED format, so we split it into several fields" << endl;

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
          for ( int elemShift = 0; elemShift < nbElems; )
            elemShift += fld->setValues( valPtr, iSub++, elemShift );
          setTS( fld, values, medFields, mesh );
        }
      else
        {
          fld->setValues( valPtr, iSub++ );
          setTS( fld, values, medFields, mesh, iSub );
        }
    }
}

//================================================================================
/*!
 * \brief Store value array of a field into med fields
 */
//================================================================================

void IntermediateMED::setTS( SauvUtilities::DoubleField*  fld,
                             ParaMEDMEM::DataArrayDouble* values,
                             ParaMEDMEM::MEDFileFields*   medFields,
                             ParaMEDMEM::MEDFileUMesh*    mesh,
                             const int                    iSub)
{
  // analyze a field support
  const Group* support = fld->getSupport();
  int dimRel;
  const bool onAll = isOnAll( support, dimRel );
  if ( !onAll && support->_name.empty() )
    {
      const_cast<Group*>(support)->_name += "PFL_" + fld->_name;
      support->_medGroup->setName( support->_name.c_str() );
    }

  // make a time-stamp
  MEDCouplingFieldDouble * timeStamp = MEDCouplingFieldDouble::New( fld->getMedType(),
                                                                    fld->getMedTimeDisc() );
  timeStamp->setName( fld->_name.c_str() );
  timeStamp->setDescription( fld->_description.c_str() );
  MEDCouplingAutoRefCountObjectPtr< MEDCouplingUMesh > dimMesh = mesh->getMeshAtLevel( dimRel );
  timeStamp->setMesh( dimMesh );
  for ( size_t i = 0; i < (size_t)fld->_sub[iSub].nbComponents(); ++i )
    values->setInfoOnComponent( i, fld->_sub[iSub]._comp_names[ i ].c_str() );
  timeStamp->setArray( values );
  values->decrRef();

  // get a field to add the time-stamp
  bool isNewMedField = false;
  if ( !fld->_curMedField || fld->_name != fld->_curMedField->getName() )
    {
      fld->_curMedField = MEDFileFieldMultiTS::New();
      isNewMedField = true;
    }

  // set an order
  timeStamp->setOrder( fld->_curMedField->getNumberOfTS() );

  // add the time-stamp
  if ( onAll )
    fld->_curMedField->appendFieldNoProfileSBT( timeStamp );
  else
    fld->_curMedField->appendFieldProfile( timeStamp, mesh, dimRel, support->_medGroup );
  timeStamp->decrRef();

  if ( isNewMedField ) // timeStamp must be added before this
    {
      medFields->pushField( fld->_curMedField );
      fld->_curMedField->decrRef();
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
  string base = fld->_name;
  if ( base.empty() )
    {
      base = "F_";
    }
  else
    {
      string::size_type pos = base.rfind('_');
      if ( pos != string::npos )
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
  Group* grp1 = _sub[0]._support;
  for ( size_t i = 1; i < _sub.size() && sameSupports; ++i )
    sameSupports = ( grp1 == _sub[i]._support );

  return sameSupports;
}

//================================================================================
/*!
 * \brief True if the field can be converted into the med field
 */
//================================================================================

bool DoubleField::isMedCompatible() const
{
  for ( size_t iSub = 0; iSub < _sub.size(); ++iSub )
    {
      if ( !getSupport(iSub) || !getSupport(iSub)->_medGroup )
        THROW_IK_EXCEPTION("SauvReader INTERNAL ERROR: NULL field support");

      if ( !_sub[iSub].isValidNbGauss() )
        {
          cout << "Skip field <" << _name << "> : different nb of gauss points in components" <<endl;
          return false;
        }
    }
  // check if there are no gauss or nbGauss() == nbCellNodes,
  // else we lack info on gauss point localization
  // TODO?
  return true;
}

//================================================================================
/*!
 * \brief return true if all sub-components has same components and same nbGauss
 */
//================================================================================

bool DoubleField::hasSameComponentsBySupport() const
{
  vector< _Sub_data >::const_iterator sub_data = _sub.begin();
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

ParaMEDMEM::TypeOfField DoubleField::getMedType( const int iSub ) const
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

ParaMEDMEM::TypeOfTimeDiscretization DoubleField::getMedTimeDisc() const
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
  const vector< double > * compValues = &_comp_values[ iComp ];

  const vector< unsigned >& relocTable = getSupport( iSub )->_relocTable;

  // Set values

  const int nbElems      = _sub[iSub]._support->size();
  const int nbGauss      = _sub[iSub].nbGauss();
  const int nbComponents = _sub[iSub].nbComponents();
  const int nbValsByElem = nbComponents * nbGauss;
  // check nb values
  int nbVals = 0;
  for ( iComp = 0; iComp < nbComponents; ++iComp )
    nbVals += compValues[iComp].size();
  if ( nbVals != nbElems * nbValsByElem )
    THROW_IK_EXCEPTION("SauvMedConvertor.cxx: support size mismatches field size");
  // compute nb values in previous subs
  int valsShift = 0;
  for ( int iS = iSub-1, shift = elemShift; shift > 0; )
  {
    int nbE = _sub[iS]._support->size();
    shift -= nbE;
    valsShift += nbE * _sub[iS].nbComponents() * _sub[iS].nbGauss();
  }
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
CellsByDimIterator::CellsByDimIterator( const IntermediateMED & medi, int dim)
{
  myImed = & medi;
  init( dim );
}
/*!
 * \brief Initialize iteration on cells of given dimention
 */
void CellsByDimIterator::init(const int  dim)
{
  myCurType = -1;
  myTypeEnd = INTERP_KERNEL::NORM_HEXA20 + 1;
  myDim = dim;
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
  int type = myCurType;
  if ( !last )
    while ( type < myTypeEnd && myImed->_cellsByType[type].empty() )
      ++type;
  return type < myTypeEnd ? getDimension( TCellType( type )) : 4;
}
// END CellsByDimIterator ========================================================

