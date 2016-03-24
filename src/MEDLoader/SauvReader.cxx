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
// File      : SauvReader.cxx
// Created   : Tue Aug 16 13:57:42 2011
// Author    : Edward AGAPOV (eap)
//

#include "SauvReader.hxx"

#include "SauvMedConvertor.hxx"
#include "MCAuto.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "MEDCouplingRefCountObject.hxx"

#include <cstring>
#include <sstream>
#include <iostream>

using namespace MEDCoupling;
using namespace SauvUtilities;
using namespace std;

#define GIBI_EQUAL(var_str, stat_str) (strncmp (var_str, stat_str, strlen(stat_str)) == 0)

//================================================================================
/*!
 * \brief Creates a reader of a given sauve file
 */
//================================================================================

SauvReader* SauvReader::New(const std::string& fileName)
{
  if ( fileName.empty() ) THROW_IK_EXCEPTION("Invalid file name");

  MEDCoupling::MCAuto< SauvUtilities::FileReader> parser;

  // try to open as XRD
  parser = new XDRReader( fileName.c_str() );
  if ( parser->open() )
    {
      SauvReader* reader = new SauvReader;
      reader->_fileReader = parser.retn();
      return reader;
    }

  // try to open as ASCII
  parser = new ASCIIReader( fileName.c_str() );
  if ( parser->open() )
    {
      SauvReader* reader = new SauvReader;
      reader->_fileReader = parser.retn();
      return reader;
    }

  THROW_IK_EXCEPTION("Unable to open file |"<< fileName << "|");
}
//================================================================================
/*!
 * \brief Destructor
 */
//================================================================================

SauvReader::~SauvReader()
{
  _fileReader->decrRef();
}

std::size_t SauvReader::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> SauvReader::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

//================================================================================
/*!
 * \brief Return current line of ASCII file to report an error
 */
//================================================================================

std::string SauvReader::lineNb() const
{
  if ( isASCII() )
    return string(" (line #") + SauvUtilities::toString
      ( static_cast<SauvUtilities::ASCIIReader*>( _fileReader )->lineNb() ) + ")";

  return "";
}

//================================================================================
/*!
 * \brief Reads contents of the sauve file and convert it to MEDFileData
 */
//================================================================================

MEDCoupling::MEDFileData * SauvReader::loadInMEDFileDS()
{
  SauvUtilities::IntermediateMED iMed; // intermadiate DS
  _iMed = &iMed;

  char* line; // a current line
  const char* enregistrement_type=" ENREGISTREMENT DE TYPE";

  while ( getNextLine(line, /*raiseOEF=*/false)) // external loop looking for "ENREGISTREMENT DE TYPE"
    {
      if ( isASCII() && !GIBI_EQUAL( line, enregistrement_type ))
        continue; // "ENREGISTREMENT DE TYPE" not found -> read the next line

      // read the number of a record
      int recordNumber;
      if ( isASCII() )
        recordNumber = atoi( line + strlen(enregistrement_type) + 1 );
      else
        recordNumber = getInt();

      // read the record
      if ( recordNumber == 2 )
        readRecord2();
      else if (recordNumber == 4 )
        readRecord4();
      else if (recordNumber == 7 )
        readRecord7();
      else if (recordNumber == 5 )
        break; // stop reading
      else
        if ( !isASCII() )
          THROW_IK_EXCEPTION("XDR : ENREGISTREMENT DE TYPE " << recordNumber << " not implemented!!!");
    }

  MEDCoupling::MEDFileData* medFileData = iMed.convertInMEDFileDS();

  return medFileData;
}

//================================================================================
/*!
 * \brief Reads "ENREGISTREMENT DE TYPE 4"
 */
//================================================================================

void SauvReader::readRecord4()
{
  if ( !isASCII() )
    {
      getInt(); // skip NIVEAU
      getInt(); // skip ERREUR
      _iMed->_spaceDim = getInt();
      getFloat(); // skip DENSITE
    }
  else
    {
      char* line;
      getNextLine(line);
      const char* s = " NIVEAU  15 NIVEAU ERREUR   0 DIMENSION";
      _iMed->_spaceDim = atoi( line + strlen( s ) + 1 );
      if ( !GIBI_EQUAL( line, " NIVEAU" ))
        THROW_IK_EXCEPTION( "Could not read space dimension" << lineNb() );
      }
  if ( _iMed->_spaceDim < 1 )
    THROW_IK_EXCEPTION( "Invalid space dimension:" << _iMed->_spaceDim );
}

//================================================================================
/*!
 * \brief Reads "ENREGISTREMENT DE TYPE 7"
 */
//================================================================================

void SauvReader::readRecord7()
{
  if ( !isASCII() )
    {
      getInt(); // skip NOMBRE INFO CASTEM2000
      getInt(); // skip IFOUR
      getInt(); // skip NIFOUR
      getInt(); // skip IFOMOD
      getInt(); // skip IECHO
      getInt(); // skip IIMPI
      getInt(); // skip IOSPI
      getInt(); // skip ISOTYP
      getInt(); // skip NSDPGE
    }
  else
    {
      // skip 3 lines:
      // NOMBRE INFO CASTEM2000   8
      // IFOUR   2 NIFOUR   0 IFOMOD   2 IECHO   1 IIMPI   0 IOSPI   0 ISOTYP   1
      // NSDPGE     0
      char* line;
      getNextLine(line);
      getNextLine(line);
      getNextLine(line);
    }
}

//================================================================================
/*!
 * \brief Reads the pile number, nb of objects and nb named of objects
 */
//================================================================================

int SauvReader::readPileNumber(int& nbNamedObjects, int& nbObjects)
{
  // FORMAT(' PILE NUMERO',I4,'NBRE ObjectS NOMMES',I8,'NBRE ObjectS',I8)
  int pileNumber;
  if ( !isASCII() )
    {
      initIntReading(3);
      pileNumber     = getInt(); next();
      nbNamedObjects = getInt(); next();
      nbObjects      = getInt(); next();
    }
  else
    {
      char* line;
      getNextLine(line);
      const char *s1 = " PILE NUMERO", *s2 = "NBRE ObjectS NOMMES", *s3 = "NBRE ObjectS";
      if ( ! GIBI_EQUAL( line, s1 ) )
        THROW_IK_EXCEPTION("Could not read the pile number " << lineNb() );
      line           = line + strlen(s1);
      pileNumber     = atoi( line );
      line           = line + 4 + strlen(s2);
      nbNamedObjects = atoi( line );
      line           = line + 8 + strlen(s3);
      nbObjects      = atoi( line );
    }
  if ( nbNamedObjects<0 )
    THROW_IK_EXCEPTION("Invalid nb of named objects: " << nbNamedObjects  << lineNb() );
  if ( nbObjects<0)
    THROW_IK_EXCEPTION("Invalid nb of objects: " << nbObjects  << lineNb() );
  // It appears to be a valid case
  // if ( nbObjects<nbNamedObjects)
  //   THROW_IK_EXCEPTION("In PILE " << pileNumber <<
  //                      " nb of objects is less than nb of named objects"  << lineNb() );
  return pileNumber;
}

//================================================================================
/*!
 * \brief Reads "ENREGISTREMENT DE TYPE 2"
 */
//================================================================================

void SauvReader::readRecord2()
{
  if ( _iMed->_spaceDim == 0 )
    THROW_IK_EXCEPTION("Missing ENREGISTREMENT DE TYPE   4");

  // read a pile number
  int pileNumber, nbNamedObjects, nbObjects;
  pileNumber = readPileNumber(nbNamedObjects, nbObjects);

  if ( !_encounteredPiles.insert( pileNumber ).second && // piles may repeat
       isASCII())
    return;

  // read object names and their indices
  vector<string> objectNames(nbNamedObjects);
  for ( initNameReading( nbNamedObjects ); more(); next() )
    objectNames[ index() ] = getName();

  vector<int> nameIndices(nbNamedObjects);
  for ( initIntReading( nbNamedObjects ); more(); next() )
    nameIndices[ index() ] = getInt();

  switch ( pileNumber )
    {
    case PILE_SOUS_MAILLAGE:
      read_PILE_SOUS_MAILLAGE(nbObjects, objectNames, nameIndices);
      break;
    case PILE_NODES_FIELD:
      read_PILE_NODES_FIELD(nbObjects, objectNames, nameIndices);
      break;
    case PILE_TABLES:
      read_PILE_TABLES(nbObjects, objectNames, nameIndices);
      break;
    case PILE_LREEL:
      read_PILE_LREEL(nbObjects, objectNames, nameIndices);
      break;
    case PILE_LOGIQUES:
      read_PILE_LOGIQUES(nbObjects, objectNames, nameIndices);
      break;
    case PILE_FLOATS:
      read_PILE_FLOATS(nbObjects, objectNames, nameIndices);
      break;
    case PILE_INTEGERS:
      read_PILE_INTEGERS(nbObjects, objectNames, nameIndices);
      break;
    case PILE_STRINGS:
      read_PILE_STRINGS(nbObjects, objectNames, nameIndices);
      break;
    case PILE_LMOTS:
      read_PILE_LMOTS(nbObjects, objectNames, nameIndices);
      break;
    case PILE_NOEUDS:
      read_PILE_NOEUDS(nbObjects, objectNames, nameIndices);
      break;
    case PILE_COORDONNEES:
      read_PILE_COORDONNEES(nbObjects, objectNames, nameIndices);
      break;
    case PILE_MODL:
      read_PILE_MODL(nbObjects, objectNames, nameIndices);
      break;
    case PILE_FIELD:
      read_PILE_FIELD(nbObjects, objectNames, nameIndices);
      break;
    default:
      if ( !isASCII() )
        THROW_IK_EXCEPTION("XDR : reading PILE " << pileNumber << " not implemented !!!");
    }
}

//================================================================================
/*!
 * \brief Reads "PILE NUMERO   1": gibi sub-meshes that are converted into med groups
 */
//================================================================================

void SauvReader::read_PILE_SOUS_MAILLAGE(const int                 nbObjects,
                                         std::vector<std::string>& objectNames,
                                         std::vector<int>&         nameIndices)
{
  _iMed->_groups.reserve(nbObjects*2); // fields may add some groups

  char* line;
  map<int,int> strangeGroupType;
  int i;

  for (int object=0; object!=nbObjects; ++object) // loop on sub-groups
    {
      initIntReading( 5 );
      int castemCellType = getIntNext();
      int nbSubGroups    = getIntNext();
      int nbReferences   = getIntNext();
      int nbNodesPerElem = getIntNext();
      int nbElements     = getIntNext();

      _iMed->_groups.push_back(Group());
      SauvUtilities::Group & group = _iMed->_groups.back();

      // Issue 0021311. Allocate places for names of referring groups
      // that will be possibly filled after reading long names from
      // PILE_TABLES and PILE_STRINGS
      group._refNames.resize( nbReferences );

      // castemCellType=0 corresponds to a sub-mesh composed of other sub-meshes
      if (castemCellType==0 && nbSubGroups>0)
        {
          group._groups.resize( nbSubGroups );
          for ( initIntReading( nbSubGroups ); more(); next() )
            group._groups[ index() ] = & _iMed->_groups[ getInt() - 1 ];
          //std::sort( group._groups.begin(), group._groups.end() ); // for _groups comparison in getFieldSupport()
        }
      // skip references
      if ( isASCII() )
        for ( i = 0; i < nbReferences; i += 10 ) // FORMAT(10I8)
          getNextLine(line);
      else
        for (initIntReading(nbReferences); more(); next());

      // skip colors
      if ( isASCII() )
        for ( i = 0; i < nbElements; i += 10 )
          getNextLine(line);
      else
        for (initIntReading(nbElements); more(); next());

      // not a composite group
      if (castemCellType>0 && nbSubGroups==0)
        {
          group._cellType = SauvUtilities::gibi2medGeom(castemCellType);

          initIntReading( nbElements * nbNodesPerElem );
          if ( group._cellType == INTERP_KERNEL::NORM_ERROR ) // look for group end
            {
              for ( ; more();  next());
              strangeGroupType.insert( make_pair( object, castemCellType ));
            }
          else
            {
              // if ( group._cellType == MED_POINT1 ) group._cellType = NORM_ERROR; // issue 21199

              // read connectivity of elements of a group
              SauvUtilities::Cell ma( nbNodesPerElem );
              SauvUtilities::Node* pNode;
              group._cells.resize( nbElements );
              for ( i = 0; i < nbElements; ++i )
                {
                  ma.init();
                  for ( int n = 0; n < nbNodesPerElem; ++n )
                    {
                      int nodeID = getIntNext();
                      pNode = _iMed->getNode( nodeID );
                      ma._nodes[n] = pNode;
                      _iMed->_nbNodes += ( !pNode->isUsed() );
                      pNode->_number = nodeID;
                    }
                  ma._number = _iMed->getNbCellsOfType( group._cellType ) + 1;
                  group._cells[i] = _iMed->insert( group._cellType, ma );
                }
            }
        }
    } // loop on groups

  // set group names
  for (i=0; i!=(int)objectNames.size(); ++i)
    {
      int grpID = nameIndices[i];
      SauvUtilities::Group & grp = _iMed->_groups[ grpID-1 ];
      if ( !grp._name.empty() ) // a group has several names
        { // create a group with subgroup grp and named grp.name
          SauvUtilities::Group* newGroup = _iMed->addNewGroup();
          newGroup->_groups.push_back( &_iMed->_groups[ grpID-1 ]);
          newGroup->_name = grp._name;
        }
      grp._name=objectNames[i];
#ifdef _DEBUG
      map<int,int>::iterator it = strangeGroupType.find( grpID - 1 );
      if ( it != strangeGroupType.end() )
        cout << "Skip " << grp._name << " of not supported CASTEM type: " << it->second << endl;
#endif
    }
} // read_PILE_SOUS_MAILLAGE()

//================================================================================
/*!
 * \brief Skip "PILE NUMERO  18" of XDR file
 */
//================================================================================

void SauvReader::read_PILE_LREEL (const int nbObjects, std::vector<std::string>&, std::vector<int>&)
{
  if ( isXRD() )
    {
      for (int object=0; object!=nbObjects; ++object) // pour chaque Group
        {
          initIntReading(1);
          int nb_vals = getIntNext();
          initDoubleReading(nb_vals);
          for(int i=0; i<nb_vals; i++) next();
        }
    }
}

//================================================================================
/*!
 * \brief Skip "PILE NUMERO  24" of XDR file
 */
//================================================================================

void SauvReader::read_PILE_LOGIQUES (const int, std::vector<std::string>&, std::vector<int>&)
{
  if ( isXRD() )
    {
      initIntReading(1);
      int nb_vals = getIntNext();
      initIntReading(nb_vals);
      for(int i=0; i<nb_vals; i++) next();
    }
}

//================================================================================
/*!
 * \brief Skip "PILE NUMERO  25" of XDR file
 */
//================================================================================

void SauvReader::read_PILE_FLOATS (const int, std::vector<std::string>&, std::vector<int>&)
{
  if ( isXRD() )
    {
      initIntReading(1);
      int nb_vals = getIntNext();
      initDoubleReading(nb_vals);
      for(int i=0; i<nb_vals; i++) next();
    }
}

//================================================================================
/*!
 * \brief Skip "PILE NUMERO  26" of XDR file
 */
//================================================================================

void SauvReader::read_PILE_INTEGERS (const int, std::vector<std::string>&, std::vector<int>&)
{
  if ( isXRD() )
    {
      initIntReading(1);
      int nb_vals = getIntNext();
      initIntReading(nb_vals);
      for(int i=0; i<nb_vals; i++) next();
    }
}

//================================================================================
/*!
 * \brief Skip "PILE NUMERO  29" of XDR file
 */
//================================================================================

void SauvReader::read_PILE_LMOTS (const int nbObjects, std::vector<std::string>&, std::vector<int>&)
{
  if ( isXRD() )
    {
      for (int object=0; object!=nbObjects; ++object) // pour chaque Group
        {
          initIntReading(2);
          int len = getIntNext();
          int nb_vals = getIntNext();
          int nb_char = len*nb_vals;
          int nb_char_tmp = 0;
          int fixed_length = 71;
          while (nb_char_tmp < nb_char)
            {
              int remain_len = nb_char - nb_char_tmp;
              int width;
              if ( remain_len > fixed_length )
                {
                  width = fixed_length;
                }
              else
                {
                  width = remain_len;
                }
              initNameReading(1, width);
              next();
              nb_char_tmp += width;
            }
        }
    }
}

//================================================================================
/*!
 * \brief Skip "PILE NUMERO  38" of XDR file
 */
//================================================================================

void SauvReader::read_PILE_MODL (const int nbObjects, std::vector<std::string>&, std::vector<int>&)
{
  if ( isXRD() )
    {
      for (int object=0; object!=nbObjects; ++object) // pour chaque Group
        {
          // see wrmodl.eso
          initIntReading(10);
          int n1  = getIntNext();
          int nm2 = getIntNext();
          int nm3 = getIntNext();
          int nm4 = getIntNext();
          int nm5 = getIntNext();
          int n45 = getIntNext();
          /*int nm6 =*/ getIntNext();
          /*int nm7 =*/ getIntNext();
          next();
          next();
          int nm1 = n1 * n45;
          int nm9 = n1 * 16;
          for (initIntReading(nm1); more(); next());
          for (initIntReading(nm9); more(); next());
          for (initNameReading(nm5, 8); more(); next());
          for (initNameReading(nm2, 8); more(); next());
          for (initNameReading(nm3, 8); more(); next());
          for (initIntReading(nm4); more(); next());
        }
    }
} // Fin case pile 38

//================================================================================
/*!
 * \brief Read "PILE NUMERO  32": links to node coordinates
 */
//================================================================================

void SauvReader::read_PILE_NOEUDS (const int nbObjects, std::vector<std::string>&, std::vector<int>&)
{
  initIntReading(1);
  int nb_indices = getIntNext();

  if (nb_indices != nbObjects)
    THROW_IK_EXCEPTION("Error of reading PILE NUMERO  " << PILE_NOEUDS << lineNb() );

  for ( initIntReading( nbObjects ); more(); next() )
    {
      int coordID = getInt();
      _iMed->getNode( index()+1 )->_coordID = coordID;
    }
}

//================================================================================
/*!
 * \brief Read "PILE NUMERO  33": node coordinates
 */
//================================================================================

void SauvReader::read_PILE_COORDONNEES (const int nbObjects, std::vector<std::string>&, std::vector<int>&)
{
  initIntReading(1);
  int nbReals = getIntNext();

  if ( nbReals < (int)(_iMed->_nbNodes*(_iMed->_spaceDim+1)) )
    THROW_IK_EXCEPTION("Error of reading PILE NUMERO  " << PILE_COORDONNEES << lineNb() );

  // there are coordinates + density for each node
  _iMed->_coords.resize( nbReals - nbReals/(_iMed->_spaceDim+1));
  double* coordPtr = &_iMed->_coords[0];

  initDoubleReading( nbReals );
  while ( more() )
    {
      for (unsigned j = 0; j < _iMed->_spaceDim; ++j, next())
        *coordPtr++ = getDouble();
      // skip density
      getDouble();
      next();
    }
}

//================================================================================
/*!
 * \brief Find or create a Group equal to a given field support
 */
//================================================================================

void SauvReader::setFieldSupport(const vector<SauvUtilities::Group*>& supports,
                                 SauvUtilities::DoubleField*          field)
{
  SauvUtilities::Group* group = NULL;
  set<SauvUtilities::Group*> sup_set( supports.begin(), supports.end() );
  if ( sup_set.size() == 1 ) // one or equal supports
    {
      group = supports[0];
    }
  else
    {
      // check if sub-components are on cells of different types
      map<int,int> nbGaussByCellType;
      for ( size_t i = 0; i < supports.size(); ++i )
        {
          map<int,int>::iterator ct2ng = nbGaussByCellType.find( supports[i]->_cellType );
          if ( ct2ng == nbGaussByCellType.end() )
            nbGaussByCellType[ supports[i]->_cellType ] = field->_sub[i].nbGauss();
          else if ( ct2ng->second != field->_sub[i].nbGauss() )
            return;
        }
      bool isSameCellType = ( nbGaussByCellType.size() == 1 );
      // try to find an existing composite group with the same sub-groups
      if ( isSameCellType )
        for ( size_t i = 0; i < _iMed->_groups.size() && !group; ++i )
          {
            Group & grp = _iMed->_groups[i];
            if (sup_set.size() == grp._groups.size())
              {
                bool sameOrder = true;
                for ( size_t j = 0; j < supports.size() && sameOrder; ++j )
                  sameOrder = ( supports[j] == grp._groups[ j % grp._groups.size() ]);
                if ( sameOrder )
                  group = & _iMed->_groups[i];
              }
          }
      if ( !group ) // no such a group, add a new one
        {
          vector<SauvUtilities::Group*> newGroups( supports.begin(),
                                                   supports.begin() + sup_set.size() );
          // check if supports includes newGroups in the same order
          bool sameOrder = true;
          for ( size_t j = newGroups.size(); j < supports.size() && sameOrder; ++j )
            sameOrder = ( supports[j] == newGroups[ j % newGroups.size() ]);
          if ( sameOrder )
            {
              group = _iMed->addNewGroup( & newGroups );
              group->_groups.swap( newGroups );
            }
        }
      // sort field sub-components and supports by cell type
      if ( group && !isSameCellType )
        {
          // sort groups
          vector<SauvUtilities::Group*>& groups = group->_groups;
          bool isModified = false, isSwapped = true;
          while ( isSwapped )
            {
              isSwapped = false;
              for ( size_t i = 1; i < groups.size(); ++i )
                {
                  int nbN1 = groups[i-1]->empty() ? 0 : groups[i-1]->_cells[0]->_nodes.size();
                  int nbN2 = groups[i  ]->empty() ? 0 : groups[i  ]->_cells[0]->_nodes.size();
                  if ( nbN1 > nbN2 )
                    {
                      isSwapped = isModified = true;
                      std::swap( groups[i], groups[i-1] );
                    }
                }
            }
          // relocate sub-components according to a new order of groups
          if ( isModified )
            {
              vector< DoubleField::_Sub_data > newSub   ( field->_sub.size() );
              vector< vector< double > >       newValues( field->_comp_values.size() );
              size_t iFromSub = 0, iNewSub = 0, iNewComp = 0;
              for ( ; iFromSub < field->_sub.size(); iFromSub += groups.size() )
                {
                  size_t iFromComp = iNewComp;
                  for ( size_t iG = 0; iG < groups.size(); ++iG )
                    {
                      size_t iComp = iFromComp;
                      for ( size_t iSub = iFromSub; iSub < field->_sub.size(); ++iSub )
                        if ( field->_sub[ iSub ]._support == groups[ iG ] )
                          {
                            newSub[ iNewSub++ ] = field->_sub[ iSub ];
                            int iC = 0, nbC = field->_sub[ iSub ].nbComponents();
                            for ( ; iC < nbC; ++iC )
                              newValues[ iNewComp++ ].swap( field->_comp_values[ iComp++ ]);
                            break;
                          }
                        else
                          {
                            iComp += field->_sub[ iSub ].nbComponents();
                          }
                    }
                }
              field->_sub.swap( newSub );
              field->_comp_values.swap( newValues );
            }
        }
    }
  if ( group )
    group->_isProfile = true;

  field->_group = group;
}

//================================================================================
/*!
 * \brief Set field names
 */
//================================================================================

void SauvReader::setFieldNames(const vector<SauvUtilities::DoubleField* >& fields,
                               const vector<string>&                       objets_nommes,
                               const vector<int>&                          indices_objets_nommes)
{
  unsigned i;
  for ( i = 0; i < indices_objets_nommes.size(); ++i )
    {
      int fieldIndex = indices_objets_nommes[ i ];
      if ( fields[ fieldIndex - 1 ] )
        fields[ fieldIndex - 1 ]->_name = objets_nommes[ i ];
    }
}

//================================================================================
/*!
 * \brief Read "PILE NUMERO   2": NODE FIELDS
 */
//================================================================================

void SauvReader::read_PILE_NODES_FIELD (const int                 nbObjects,
                                        std::vector<std::string>& objectNames,
                                        std::vector<int>&         nameIndices)
{
  _iMed->_nodeFields.resize( nbObjects, (SauvUtilities::DoubleField*) 0 );
  for (int object=0; object!=nbObjects; ++object) // loop on fields
    {
      // EXAMPLE ( with no values )

      // (1)       4       7       2       1
      // (2)     -88       0       3     -89       0       1     -90       0       2     -91
      // (2)       0       1
      // (3) FX   FY   FZ   FZ   FX   FY   FLX
      // (4)       0       0       0       0       0       0       0
      // (5)           cree  par  muc pri
      // (6)
      // (7)       2

      // (1): nb subcomponents, nb components(total), IFOUR, nb attributes
      int nb_sub, total_nb_comp, nb_attr;
      int i_sub, i_comp;
      initIntReading( 4 );
      nb_sub        = getIntNext();
      total_nb_comp = getIntNext();
      next(); // ignore IFOUR
      nb_attr       = getIntNext();
      if ( nb_sub < 0 || total_nb_comp < 0 || nb_attr < 0 )
        THROW_IK_EXCEPTION("Error of field reading " << lineNb());

      // (2) loop on subcomponents of a field, for each read
      // (a) support, (b) number of values and (c) number of components
      vector<Group*> supports( nb_sub );
      vector<int> nb_values  ( nb_sub );
      vector<int> nb_comps   ( nb_sub );
      int total_nb_values = 0;
      initIntReading( nb_sub * 3 );
      for ( i_sub = 0; i_sub < nb_sub; ++i_sub )
        {
          int supId = -getIntNext(); // (a) reference to support
          if ( supId < 1 || supId > (int)_iMed->_groups.size() )
            THROW_IK_EXCEPTION("Wrong mesh reference: "<< supId << lineNb() );
          supports[ i_sub ] = &_iMed->_groups[ supId-1 ]; // (a) reference to support

          nb_values[ i_sub ] = getIntNext();    // (b) nb points
          total_nb_values += nb_values[ i_sub ];
          if ( nb_values[ i_sub ] < 0 )
            THROW_IK_EXCEPTION(" Wrong nb of points: " << nb_values[ i_sub ]  << lineNb() );
          nb_comps[ i_sub ] = getInt(); next();     // (c) nb of components in i_sub
        }

      // create a field if there are values
      SauvUtilities::DoubleField* fdouble = 0;
      if ( total_nb_values > 0 )
        fdouble = new DoubleField( nb_sub, total_nb_comp );
      _iMed->_nodeFields[ object ] = fdouble;

      // (3) component names
      initNameReading( total_nb_comp, 4 );
      for ( i_sub = 0; i_sub < nb_sub; ++i_sub )
        {
          // store support id and nb components of a sub
          if ( fdouble )
            fdouble->_sub[ i_sub ].setData( nb_comps[ i_sub ], supports[ i_sub ] );
          for ( i_comp = 0; i_comp < nb_comps[ i_sub ]; ++i_comp, next() )
            {
              // store component name
              string compName = getName();
              if ( fdouble )
                fdouble->_sub[ i_sub ].compName( i_comp ) = compName;
            }
        }
      // (4) nb harmonics ( ignored )
      for ( initIntReading( total_nb_comp ); more(); next() );
      // (5) TYPE ( ignored )
      for (initNameReading(1, /*length=*/71); more(); next());
      // (6) TITRE ( ignored )
      for (initNameReading(1, /*length=*/71); more(); next());
      // (7) attributes ( ignored )
      for ( initIntReading( nb_attr ); more(); next() );

      for ( i_sub = 0; i_sub < nb_sub; ++i_sub )
        {
          // loop on components: read values
          initDoubleReading( nb_values[ i_sub ] * nb_comps[ i_sub ] );
          for ( i_comp = 0; i_comp < nb_comps[ i_sub ]; ++i_comp )
            {
              if ( fdouble )
                {
                  vector<double>& vals = fdouble->addComponent( nb_values[ i_sub ] );
                  for ( int i = 0; more() && i < nb_values[ i_sub ]; next(), ++i )
                    vals[ i ] = getDouble();
                }
              else
                {
                  for ( int i = 0; i < nb_values[ i_sub ]; next(), ++i );
                }
            }
        } // loop on subcomponents of a field

      // set a supporting group including all subs supports but only
      // if all subs have the same components
      if ( fdouble && fdouble->hasSameComponentsBySupport() )
        setFieldSupport( supports, fdouble );
      else
        for ( i_sub = 0; i_sub < nb_sub; ++i_sub )
          fdouble->_sub[ i_sub ]._support->_isProfile = true;

    } // end loop on field objects

  // set field names
  setFieldNames( _iMed->_nodeFields, objectNames, nameIndices );

}  // read_PILE_NODES_FIELD()

//================================================================================
/*!
 * \brief Read "PILE NUMERO  39": FIELDS
 */
//================================================================================

void SauvReader::read_PILE_FIELD (const int                 nbObjects,
                                  std::vector<std::string>& objectNames,
                                  std::vector<int>&         nameIndices)
{
  // REAL EXAMPLE

  // (1)        1       2       6      16
  // (2)                                                         CARACTERISTIQUES
  // (3)      -15  317773       4       0       0       0      -2       0       3
  // (4)             317581
  // (5)  0
  // (6)   317767  317761  317755  317815
  // (7)  YOUN     NU       H        SIGY
  // (8)  REAL*8            REAL*8            REAL*8            REAL*8
  // (9)        1       1       0       0
  // (10)  2.00000000000000E+05
  // (11)       1       1       0       0
  // (12)  3.30000000000000E-01
  // (13)       1       1       0       0
  // (14)  1.00000000000000E+04
  // (15)       6     706       0       0
  // (16)  1.00000000000000E+02  1.00000000000000E+02  1.00000000000000E+02
  // (17)  1.00000000000000E+02  1.00000000000000E+02  1.00000000000000E+02
  // (18)  ...

  _iMed->_cellFields.resize( nbObjects, (SauvUtilities::DoubleField*) 0 );
  for (int object=0; object!=nbObjects; ++object) // pour chaque field
    {
      initIntReading( 4 );
      int i_sub, nb_sub = getIntNext(); // (1) <nb_sub> 2 6 <title length>
      next(); // skip "2"
      next(); // skip "6"
      int title_length = getIntNext(); // <title length>
      if ( nb_sub < 1 )
        THROW_IK_EXCEPTION("Error of field reading: wrong nb of subcomponents " << nb_sub << lineNb() );

      string description;
      if ( title_length )
        {
          if ( isXRD() )
            {
              initNameReading(1, title_length);
              description = getName();
              next();
            }
          else
            {
              char* line;
              getNextLine( line ); // (2) title
              const int len = 72; // line length
              description = string(line + len - title_length, title_length); // title is right justified
            }
        }
      // look for a line starting with '-' : <reference to support>
      if ( isXRD() )
        {
          initIntReading( nb_sub * 9 );
        }
      else
        {
          do {
            initIntReading( nb_sub * 9 );
          } while ( getInt() >= 0 );
        }
      int total_nb_comp = 0;
      vector<Group*> supports( nb_sub );
      vector<int>     nb_comp( nb_sub );
      for ( i_sub = 0; i_sub < nb_sub; ++i_sub )
        {                                    // (3)
          int supportId     = -getIntNext(); // <reference to support>
          next();                            // ignore <address>
          nb_comp [ i_sub ] =  getIntNext(); // <nb of components in the sub>
          for ( int i = 0; i < 6; ++i ) next();  // ignore 6 ints, in example "0 0 0 -2 0 3"

          if ( supportId < 1 || supportId > (int)_iMed->_groups.size() )
            THROW_IK_EXCEPTION("Error of field reading: wrong mesh reference "<< supportId << lineNb() );
          if ( nb_comp[ i_sub ] < 0 )
            THROW_IK_EXCEPTION("Error of field reading: wrong nb of components " <<nb_comp[ i_sub ] << lineNb() );

          supports[ i_sub ] = &_iMed->_groups[ supportId-1 ];
          total_nb_comp += nb_comp[ i_sub ];
        }
      // (4) dummy strings
      for ( initNameReading( nb_sub, 17 ); more(); next() );
      // (5) dummy strings
      for ( initNameReading( nb_sub ); more(); next() );

      // loop on subcomponents of a field, each of which refers to
      // a certain support and has its own number of components;
      // read component values
      SauvUtilities::DoubleField* fdouble = 0;
      for ( i_sub = 0; i_sub < nb_sub; ++ i_sub )
        {
          vector<string> comp_names( nb_comp[ i_sub ]), comp_type( nb_comp[ i_sub ]);
          // (6) nb_comp addresses of MELVAL structure
          for ( initIntReading( nb_comp[ i_sub ] ); more(); next() );
          // (7) component names
          for ( initNameReading( nb_comp[ i_sub ] ); more(); next() )
            comp_names[ index() ] = getName();
          // (8) component type
          for ( initNameReading( nb_comp[ i_sub ], 17 ); more(); next() ) // 17 is name width
            {
              comp_type[ index() ] = getName();
              // component types must be the same
              if ( index() > 0 && comp_type[ index() ] != comp_type[ index() - 1] )
                THROW_IK_EXCEPTION( "Error of field reading: diff component types <"
                                    << comp_type[ index() ] << "> != <" << comp_type[ index() - 1 ]
                                    << ">" << lineNb() );
            }
          // now type is known, create a field, one for all subs
          bool isReal = (nb_comp[i_sub] > 0) ? (comp_type[0] == "REAL*8") : true;
          if ( !fdouble && total_nb_comp )
            {
              if ( !isReal )
                cout << "Warning: read NOT REAL field, type <" << comp_type[0] << ">" << lineNb() << endl;
              _iMed->_cellFields[ object ] = fdouble = new SauvUtilities::DoubleField( nb_sub, total_nb_comp );
              fdouble->_description = description;
            }
          // store support id and nb components of a sub
          if ( fdouble )
            fdouble->_sub[ i_sub ].setData( nb_comp[ i_sub ], supports[ i_sub ]);
          // loop on components: read values
          for ( int i_comp = 0; i_comp < nb_comp[ i_sub ]; ++i_comp )
            {
              // (9) nb of values
              initIntReading( 4 );
              int nb_val_by_elem = getIntNext();
              int nb_values      = getIntNext();
              next();
              next();
              fdouble->_sub[ i_sub ]._nb_gauss[ i_comp ] = nb_val_by_elem;

              // (10) values
              nb_values *= nb_val_by_elem;
              if ( fdouble )
                {
                  vector<double> & vals = fdouble->addComponent( nb_values );
                  for ( isReal ? initDoubleReading( nb_values ) : initIntReading( nb_values ); more(); next())
                    vals[ index() ] = getDouble();
                  // store component name
                  fdouble->_sub[ i_sub ].compName( i_comp ) = comp_names[ i_comp ];
                }
              else
                {
                  for ( isReal ? initDoubleReading( nb_values ) : initIntReading( nb_values ); more(); next() ) ;
                }
            }
        } // loop on subcomponents of a field

      // set id of a group including all sub supports but only
      // if all subs have the same nb of components
      if ( fdouble && fdouble->hasSameComponentsBySupport() )
        setFieldSupport( supports, fdouble );
      else
        for ( i_sub = 0; i_sub < nb_sub; ++i_sub )
          fdouble->_sub[ i_sub ]._support->_isProfile = true;

    } // end loop on field objects

  // set field names
  setFieldNames( _iMed->_cellFields, objectNames, nameIndices );

} // read_PILE_FIELD()

//================================================================================
/*!
 * \brief Read "PILE NUMERO  10": TABLES
 */
//================================================================================

void SauvReader::read_PILE_TABLES (const int                 nbObjects,
                                   std::vector<std::string>& objectNames,
                                   std::vector<int>&         nameIndices)
{
  // IMP 0020434: mapping GIBI names to MED names

  string table_med_mail = "MED_MAIL";
  string table_med_cham = "MED_CHAM";
  string table_med_comp = "MED_COMP";
  int table_med_mail_id = -1;
  int table_med_cham_id = -1;
  int table_med_comp_id = -1;
  for (size_t iname = 0; iname < objectNames.size(); iname++)
    if      (objectNames[iname] == table_med_mail) table_med_mail_id = nameIndices[iname];
    else if (objectNames[iname] == table_med_cham) table_med_cham_id = nameIndices[iname];
    else if (objectNames[iname] == table_med_comp) table_med_comp_id = nameIndices[iname];

  if ( isASCII() )
    if (table_med_mail_id < 0 && table_med_cham_id < 0 && table_med_comp_id < 0)
      return;

  for (int itable = 1; itable <= nbObjects; itable++)
    {
      // read tables "MED_MAIL", "MED_CHAM" and "MED_COMP", that keeps correspondence
      // between GIBI names (8 symbols if any) and MED names (possibly longer)
      initIntReading(1);
      int nb_table_vals = getIntNext();
      if (nb_table_vals < 0)
        THROW_IK_EXCEPTION("Error of reading PILE NUMERO  10" << lineNb() );

      int name_i_med_pile;
      initIntReading(nb_table_vals);
      for (int i = 0; i < nb_table_vals/4; i++)
        {
          if (itable == table_med_mail_id ||
              itable == table_med_cham_id ||
              itable == table_med_comp_id)
            {
              nameGIBItoMED name_i;
              name_i_med_pile  = getIntNext();
              name_i.med_id    = getIntNext();
              name_i.gibi_pile = getIntNext();
              name_i.gibi_id   = getIntNext();

              if (name_i_med_pile != PILE_STRINGS)
                {
                  // error: med(long) names are always kept in PILE_STRINGS
                }
              if (itable == table_med_mail_id)
                {
                  if (name_i.gibi_pile != PILE_SOUS_MAILLAGE) {
                    // error: must be PILE_SOUS_MAILLAGE
                  }
                  _iMed->_listGIBItoMED_mail.push_back(name_i);
                }
              else if (itable == table_med_cham_id)
                {
                  if (name_i.gibi_pile != PILE_FIELD &&
                      name_i.gibi_pile != PILE_NODES_FIELD)
                    {
                      // error: must be PILE_FIELD or PILE_NODES_FIELD
                    }
                  _iMed->_listGIBItoMED_cham.push_back(name_i);
                }
              else if (itable == table_med_comp_id)
                {
                  if (name_i.gibi_pile != PILE_STRINGS)
                    {
                      // error: gibi(short) names of components are kept in PILE_STRINGS
                    }
                  _iMed->_listGIBItoMED_comp.push_back(name_i);
                }
            }
          else
            {
              // pass table
              for ( int ii = 0; ii < 4; ++ii ) next();
            }
        }
    } // for (int itable = 0; itable < nbObjects; itable++)
}

//================================================================================
/*!
 * \brief Read "PILE NUMERO  27"
 */
//================================================================================

void SauvReader::read_PILE_STRINGS (const int                 nbObjects,
                                    std::vector<std::string>& objectNames,
                                    std::vector<int>&         nameIndices)
{
  // IMP 0020434: mapping GIBI names to MED names
  initIntReading(2);
  int stringLen    = getIntNext();
  int nbSubStrings = getIntNext();
  if (nbSubStrings != nbObjects)
    THROW_IK_EXCEPTION("Error of reading PILE NUMERO  27" << lineNb() );

  string aWholeString;
  if ( isXRD() )
    {
      const int fixedLength = 71;
      while ((int)aWholeString.length() < stringLen)
        {
          int remainLen = stringLen - aWholeString.length();
          int len;
          if ( remainLen > fixedLength )
            {
              len = fixedLength;
            }
          else
            {
              len = remainLen;
            }
          initNameReading(1, len);
          aWholeString += getName();
          next();
        }
    }
  else
    {
      char* line;
      const int fixedLength = 71;
      while ((int)aWholeString.length() < stringLen)
        {
          getNextLine( line );
          int remainLen = stringLen - aWholeString.length();
          if ( remainLen > fixedLength )
            {
              aWholeString += line + 1;
            }
          else
            {
              aWholeString += line + ( 72 - remainLen );
            }
        }
    }
  int prevOffset = 0;
  int currOffset = 0;
  initIntReading(nbSubStrings);
  for (int istr = 1; istr <= nbSubStrings; istr++, next())
    {
      currOffset = getInt();
      // fill mapStrings
      _iMed->_mapStrings[istr] = aWholeString.substr(prevOffset, currOffset - prevOffset);
      prevOffset = currOffset;
    }
}
