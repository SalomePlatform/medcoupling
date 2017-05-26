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
// File      : SauvMedConvertor.hxx
// Created   : Tue Aug 16 14:14:02 2011
// Author    : Edward AGAPOV (eap)
//

#ifndef __SauvMedConvertor_HXX__
#define __SauvMedConvertor_HXX__

#include "InterpKernelException.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "SauvUtilities.hxx"
#include "MCType.hxx"

#include <vector>
#include <set>
#include <map>
#include <list>
#include <algorithm>

namespace MEDCoupling
{
  class DataArrayDouble;
  class DataArrayInt;
  class MEDFileData;
  class MEDFileFields;
  class MEDFileFieldMultiTS;
  class MEDFileUMesh;
}

namespace SauvUtilities
{
  struct IntermediateMED;

  // ==============================================================================
  typedef int                TID;  // an ID countered from 1
  typedef std::pair<TID,TID> Link; // a pair of node numbers

  typedef INTERP_KERNEL::NormalizedCellType TCellType;

  // ==============================================================================
  struct Node
  {
    TID    _number;
    size_t _coordID;

    Node():_number(0){}
    bool isUsed() const { return _number != 0; }
  };

  // ==============================================================================
  struct Cell
  {
    std::vector< Node* > _nodes;
    mutable bool         _reverse; // to reverse orienation of a face only
    mutable TID*         _sortedNodeIDs; // for comparison
    mutable TID          _number;

    Cell(size_t nnNodes=0) : _nodes(nnNodes),_reverse(false),_sortedNodeIDs(0),_number(0) {}
    Cell(const Cell& ma);
    void init() const { if ( _sortedNodeIDs ) delete [] _sortedNodeIDs; _sortedNodeIDs = 0; }
    ~Cell() { init(); }

    const TID* getSortedNodes() const; // creates if needed and return _sortedNodeIDs
    bool operator < (const Cell& ma) const;
    Link link(int i) const;

  private:
    Cell& operator=(const Cell& ma);
  };
  std::ostream& operator << (std::ostream& os, const Cell& ma);

  // ==============================================================================
  struct Group
  {
    TCellType                _cellType;
    std::string              _name;
    std::vector<const Cell*> _cells;
    std::vector< Group* >    _groups;    // des sous-groupes composant le Group
    bool                     _isProfile; // is a field support or not
    std::vector<std::string> _refNames;  /* names of groups referring this one;
                                            _refNames is resized according to nb of references
                                            while reading a group (pile 1) and it is filled with
                                            names while reading long names (pile 27); each named
                                            reference is converted into a copy of the medGroup
                                            (issue 0021311)
                                         */
    MEDCoupling::DataArrayInt* _medGroup;   // result of conversion
    std::vector< unsigned >   _relocTable; // for _cells[i] gives its index in _medGroup

    bool empty() const { return _cells.empty() && _groups.empty(); }
    int  size()  const;
    Group():_cellType(INTERP_KERNEL::NORM_ERROR), _isProfile(false), _medGroup(NULL) {}
  };

  // ==============================================================================
  struct DoubleField
  {
    // a field contains several subcomponents each referring to its own support and
    // having several named components
    // ----------------------------------------------------------------------------
    struct _Sub_data // a subcomponent
    // --------------------------------
    {
      Group*                   _support;    // support
      std::vector<std::string> _comp_names; // component names
      std::vector<int>         _nb_gauss;   // nb values per element in a component

      void setData( int nb_comp, Group* supp )
      { _support = supp; _comp_names.resize(nb_comp); _nb_gauss.resize(nb_comp,1); }
      int  nbComponents() const { return _comp_names.size(); }
      std::string & compName( int i_comp ) { return _comp_names[ i_comp ]; }
      bool isSameNbGauss() const { return *std::max_element( _nb_gauss.begin(), _nb_gauss.end() ) ==
          *std::min_element( _nb_gauss.begin(), _nb_gauss.end() ); }
      int  nbGauss() const { return _nb_gauss[0] ? _nb_gauss[0] : 1; }
      bool hasGauss() const { return nbGauss() > 1; }
    };
    // ----------------------------------------------------------------------------
    TID                      _idInFile;
    std::string              _name;
    std::string              _description; // field description
    std::vector< _Sub_data > _sub;
    Group*                   _group; /* if _group == NULL then each subcomponent makes a
                                        separate med field, else all subcomponents
                                        are converted into timestamps of one med field.
                                        The latter is possible only if nb of components in all subs
                                        is the same and supports of subcomponents do not overlap
                                     */
    std::vector< std::vector< double > > _comp_values;
    MEDCoupling::MEDFileFieldMultiTS*     _curMedField;

    DoubleField( int nb_sub, int total_nb_comp )
      : _sub(nb_sub), _group(NULL), _curMedField(NULL) { _comp_values.reserve( total_nb_comp ); }
    ~DoubleField();
    std::vector< double >& addComponent( int nb_values ); // return a vector ready to fill in
    bool hasCommonSupport() const { return _group; } // true if there is one support for all subs
    bool hasSameComponentsBySupport() const;

    bool isMultiTimeStamps() const;
    bool isMedCompatible(bool& sameNbGauss) const;
    MEDCoupling::TypeOfField getMedType( const int iSub=0 ) const;
    MEDCoupling::TypeOfTimeDiscretization getMedTimeDisc() const;
    int getNbTuples( const int iSub=0 ) const;
    int getNbValuesPerElement( const int iSub=0 ) const;
    int getNbGauss( const int iSub=0 ) const;
    const Group* getSupport( const int iSub=0 ) const;
    int setValues( double * valPtr, const int iSub, const int elemShift=0 ) const;
    void splitSubWithDiffNbGauss();

    //virtual void dump(std::ostream&) const;
    //virtual ~DoubleField() {}
  };
  // ==============================================================================
  /*!
   * \if developper
   * Iterator on set of Cell's of given dimension
   * \endif
   */
  class CellsByDimIterator
  {
  public:
    CellsByDimIterator( const IntermediateMED & medi, int dim=-1); // dim=-1 - for all dimensions
    void init(const int  dim=-1);

    //!< return next set of Cell's of required dimension
    const std::set<Cell > * nextType();
    //!< return dimension of Cell's returned by the last or further next()
    int dim(const bool last=true) const;
    //!< return type of Cell's returned by the last next()
    TCellType type() const { return TCellType( myCurType ); }

  private:
    const IntermediateMED* myImed;
    int myCurType, myTypeEnd;
    int myDim;
  };

  // ==============================================================================
  /*!
   * \if developper
   * Container of Node's. Prevents re-allocation at addition of Node's
   * \endif
   */
  class NodeContainer
  {
    std::vector< std::vector< Node > > _nodes;
  public:
    Node* getNode( const TID nID )
    {
      const size_t chunkSize = 1000;
      const size_t chunkID = (nID-1) / chunkSize;
      const size_t pos     = (nID-1) % chunkSize;
      if ( _nodes.size() < chunkID+1 )
      {
        std::vector< std::vector< Node > > newNodes(chunkID+1);
        for ( size_t i = 0; i < _nodes.size(); ++i )
          newNodes[i].swap( _nodes[i] );
        for ( size_t i = _nodes.size(); i < newNodes.size(); ++i )
          newNodes[i].resize( chunkSize );
        _nodes.swap( newNodes );
      }
      return & _nodes[chunkID][pos];
    }
    bool empty() const { return _nodes.empty(); }
    size_t size() const { return empty() ? 0 : _nodes.size() * _nodes[0].size(); }
    int numberNodes();
  };

  // ==============================================================================
  /*!
   * \if developper
   * Intermediate structure used to store data read from the Sauve format file.
   * The structure provides functions that transform the stored data to the MED format
   *
   * The elements inserted in maillage are ordered in order to avoid duplicated elements.
   * \endif
   */
  struct IntermediateMED
  {
    unsigned                   _spaceDim;
    unsigned                   _nbNodes;
    NodeContainer              _points;
    std::vector<double>        _coords;
    std::vector<Group>         _groups;
    std::vector<DoubleField* > _nodeFields;
    std::vector<DoubleField* > _cellFields;

    // IMP 0020434: mapping GIBI names to MED names
    std::list<nameGIBItoMED>  _listGIBItoMED_mail; // to read from table "MED_MAIL" of PILE_TABLES
    std::list<nameGIBItoMED>  _listGIBItoMED_cham; // to read from table "MED_CHAM" of PILE_TABLES
    std::list<nameGIBItoMED>  _listGIBItoMED_comp; // to read from table "MED_COMP" of PILE_TABLES
    std::map<int,std::string> _mapStrings;         // to read from PILE_STRINGS

    IntermediateMED(): _spaceDim(0), _nbNodes(0) {}
    ~IntermediateMED();

    Node* getNode( TID nID ) { return _points.getNode( nID ); }
    int getNbCellsOfType( TCellType type ) const { return _cellsByType[type].size(); }
    const Cell* insert(TCellType type, const Cell& ma) { return &( *_cellsByType[type].insert( ma ).first ); }
    Group* addNewGroup(std::vector<SauvUtilities::Group*>* groupsToFix=0);
    MEDCoupling::MEDFileData* convertInMEDFileDS();

  private:

    MEDCoupling::MEDFileUMesh* makeMEDFileMesh();
    MEDCoupling::DataArrayDouble * getCoords();
    void setConnectivity( MEDCoupling::MEDFileUMesh* mesh, MEDCoupling::DataArrayDouble* coords );
    void setGroups( MEDCoupling::MEDFileUMesh* mesh );
    MEDCoupling::MEDFileFields * makeMEDFileFields(MEDCoupling::MEDFileUMesh* mesh);
    void setFields( SauvUtilities::DoubleField*    fld,
                    MEDCoupling::MEDFileFields*     medFields,
                    MEDCoupling::MEDFileUMesh*      mesh,
                    const TID                      castemID,
                    std::set< std::string >&       usedNames);
    void setTS( SauvUtilities::DoubleField*  fld,
                MEDCoupling::DataArrayDouble* values,
                MEDCoupling::MEDFileFields*   medFields,
                MEDCoupling::MEDFileUMesh*    mesh,
                const int                    iSub=0);
    void checkDataAvailability() const;
    void setGroupLongNames();
    void setFieldLongNames(std::set< std::string >& usedNames);
    void makeFieldNewName(std::set< std::string >&    usedNames,
                          SauvUtilities::DoubleField* fld );
    void decreaseHierarchicalDepthOfSubgroups();
    void eraseUselessGroups();
    void detectMixDimGroups();
    void orientElements2D();
    void orientElements3D();
    void orientFaces3D();
    void orientVolumes();
    void numberElements();
    bool isOnAll( const Group* grp, int & dimRel ) const;
    const double* nodeCoords( const Node* n ) { return &_coords[ (n->_coordID-1) * _spaceDim ]; }

    // IntermediateMED()
    // { myNodesNumerated = myMaillesNumerated = myGroupsTreated = false; currentTypeMailles = 0; }
    // ~IntermediateMED();

    //bool myNodesNumerated, myMaillesNumerated;

    // mailles groupped by geom type; use insert() for filling in and
    // _CellsByDimIterator for exploring it
    //std::set<_Cell> maillage;
    std::set< Cell >  _cellsByType[ INTERP_KERNEL::NORM_HEXA20 + 1 ];
    friend class CellsByDimIterator;
  };

// ==============================================================================
  /*!
   * \brief ASCII sauve file reader
   */
  class ASCIIReader : public FileReader
  {
  public:
    ASCIIReader(const char* fileName);
    virtual ~ASCIIReader();
    virtual bool isASCII() const;
    virtual bool open();
    virtual bool getNextLine (char* & line, bool raiseOEF = true );
    virtual void initNameReading(int nbValues, int width = 8);
    virtual void initIntReading(int nbValues);
    virtual void initDoubleReading(int nbValues);
    virtual bool more() const;
    virtual void next();
    virtual int    getInt() const;
    virtual float  getFloat() const;
    virtual double getDouble() const;
    virtual std::string getName() const;
    int lineNb() const { return _lineNb; }

  private:

    bool getLine(char* & line);
    void init( int nbToRead, int nbPosInLine, int width, int shift = 0 );

    // getting a line from the file
    int   _file;
    char* _start; // working buffer beginning
    char* _ptr;
    char* _eptr;
    int   _lineNb;

    // line parsing
    int _iPos, _nbPosInLine, _width, _shift;
    char* _curPos;
  };
// ==============================================================================
  /*!
   * \brief XDR (binary) sauve file reader
   */
  class XDRReader : public FileReader
  {
  public:
    XDRReader(const char* fileName);
    virtual ~XDRReader();
    virtual bool isASCII() const;
    virtual bool open();
    virtual bool getNextLine (char* & line, bool raiseOEF = true );
    virtual void initNameReading(int nbValues, int width = 8);
    virtual void initIntReading(int nbValues);
    virtual void initDoubleReading(int nbValues);
    virtual bool more() const;
    virtual void next();
    virtual int    getInt() const;
    virtual float  getFloat() const;
    virtual double getDouble() const;
    virtual std::string getName() const;

  private:

    void init( int nbToRead, int width = 0 );

    FILE* _xdrs_file;
    void* _xdrs;
    int* _xdr_ivals;
    double* _xdr_dvals;
    char* _xdr_cvals;
    int _width;
    int _xdr_kind;
    enum
      {
        _xdr_kind_null,
        _xdr_kind_char,
        _xdr_kind_int,
        _xdr_kind_double
      };
  };
}

#endif
