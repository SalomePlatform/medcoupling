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
// File      : SauvWriter.cxx
// Created   : Wed Aug 24 12:55:55 2011
// Author    : Edward AGAPOV (eap)

#include "SauvWriter.hxx"

#include "InterpKernelException.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileField.hxx"
#include "MEDFileData.hxx"
#include "CellModel.hxx"

#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdlib>
#include <iomanip>

using namespace MEDCoupling;
using namespace SauvUtilities;
using namespace std;

#define INFOS_MED(txt) cout << txt << endl;

namespace
{
  const char* zeroI8 = "       0"; // FORMAT(I8)

  // ============================================================
  // the class writes endl to the file as soon as <limit> fields
  // have been written after the last endl
  // ============================================================

  class TFieldCounter
  {
    fstream& _file;
    int _count, _limit;
  public:
    TFieldCounter(fstream& f, int limit=0): _file(f), _limit(limit) { init(); }
    void init(int limit=0) // init, is done by stop() as well
    { if (limit) _limit = limit; _count = 0; }
    void operator++(int) // next
    { if ( ++_count == _limit ) { _file << endl; init(); }}
    void stop() // init() and write endl if there was no endl after the last written field
    { if ( _count ) _file << endl; init(); }
    ~TFieldCounter() { stop(); }
  };

  //================================================================================
  /*!
   * \brief Return a name of a field support on all elements
   */
  //================================================================================

  string noProfileName( INTERP_KERNEL::NormalizedCellType type )
  {
    return "INTERP_KERNEL::NormalizedCellType_" + SauvUtilities::toString( type );
  }

  //================================================================================
  /*!
   * \brief Remove white spaces from the head and tail
   */
  //================================================================================

  string cleanName( const string& theName )
  {
    string name = theName;
    if ( !name.empty() )
      {
        // cut off leading white spaces
        string::size_type firstChar = name.find_first_not_of(" \t");
        if (firstChar < name.length())
          {
            name = name.substr(firstChar);
          }
        else
          {
            name = ""; // only whitespaces there - remove them
          }
        // cut off trailing white spaces
        string::size_type lastChar = name.find_last_not_of(" \t");
        if (lastChar < name.length())
          name = name.substr(0, lastChar + 1);
      }
    return name;
  }

  //================================================================================
  /*!
   * \brief Converts MED long names into SAUVE short ones, returnes a healed long name
   */
  //================================================================================

  string addName (map<string,int>& nameMap,
                  map<string,int>& namePrefixesMap,
                  const string&    theName,
                  const int        index)
  {
    // Converts names like:
    // MED:                       GIBI:     
    //   TEMPERATURE_FLUIDE   ->    TEMPE001
    //   TEMPERATURE_SOLIDE   ->    TEMPE002
    //   PRESSION             ->    PRESSION
    //   NU                   ->    NU      
    //   VOLUM001             ->    VOLUM001
    //   VOLUMOFOBJECT        ->    VOLUM003
    //   VOLUM002             ->    VOLUM002
    string healedName = cleanName(theName);
    int ind = index;

    if (!healedName.empty())
      {
        string name = healedName;
        int len = name.length();
        for (int i = 0; i < len; ++i)
          name[i] = toupper(name[i]);

        bool doResave = false; // only for tracing

        // I. Save a short name as it is
        if (len <= 8)
          {
            INFOS_MED("Save <" << theName << "> as <" << name << ">");

            map<string,int>::iterator it = nameMap.find(name);
            if (it != nameMap.end())
              {
                // There is already such name in the map.

                // a. Replace in the map the old pair by the current one
                int old_ind = nameMap[name];
                nameMap[name] = ind;
                // b. Rebuild the old pair (which was in the map,
                //    it seems to be built automatically by step II)
                ind = old_ind;
                // continue with step II
                doResave = true; // only for tracing
              }
            else
              {
                // Save in the map
                nameMap.insert(make_pair(name, ind));

                // Update loc_index for this name (if last free characters represents a number)
                // to avoid conflicts with long names, same in first 5 characters
                if (len == 8)
                  {
                    int new_loc_index = atoi(name.c_str() + 5);
                    if (new_loc_index > 0)
                      {
                        // prefix
                        string str = name.substr(0,5);
                        if (namePrefixesMap.find(str) != namePrefixesMap.end())
                          {
                            int old_loc_index = namePrefixesMap[str];
                            if (new_loc_index < old_loc_index) new_loc_index = old_loc_index;
                          }
                        namePrefixesMap[str] = new_loc_index;
                      }
                  }
                return healedName;
              }
          } // if (len <= 8)

        // II. Cut long name and add a numeric suffix

        // first 5 or less characters of the name
        if (len > 5) name = name.substr(0,5);

        // numeric suffix
        map<string,int>::iterator name2ind = namePrefixesMap.insert( make_pair( name, 0 )).first;
        string numSuffix = SauvUtilities::toString( ++(name2ind->second) );

        if ( numSuffix.size() + name.size() > 8 )
          THROW_IK_EXCEPTION("Can't write not unique name: " << healedName);

        if ( numSuffix.size() < 3 )
          numSuffix.insert( 0, 3 - numSuffix.size(), '0' );

        name += numSuffix;
        nameMap.insert(make_pair(name, ind));

        if (doResave)
          {
            INFOS_MED("Resave previous <" << healedName << "> as <" << name << ">");
          }
        else
          {
            INFOS_MED("Save <" << theName << "> as <" << name << ">");
          }
      }
    return healedName;
  }
}

SauvWriter::SauvWriter():_cpy_grp_if_on_single_family(false)
{
}

SauvWriter* SauvWriter::New()
{
  return new SauvWriter;
}

std::size_t SauvWriter::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> SauvWriter::getDirectChildrenWithNull() const
{
  return std::vector<const BigMemoryObject *>();
}

void SauvWriter::setCpyGrpIfOnASingleFamilyStatus(bool status)
{
  _cpy_grp_if_on_single_family=status;
}

bool SauvWriter::getCpyGrpIfOnASingleFamilyStatus() const
{
  return _cpy_grp_if_on_single_family;
}

//================================================================================
/*!
 * \brief Fills own DS by MEDFileData
 */
//================================================================================

void SauvWriter::setMEDFileDS(const MEDFileData* medData,
                              unsigned           meshIndex)
{
  if ( !medData) THROW_IK_EXCEPTION("NULL MEDFileData");

  MEDFileMeshes * meshes = medData->getMeshes();
  MEDFileFields * fields = medData->getFields();
  if ( !meshes) THROW_IK_EXCEPTION("No meshes in MEDFileData");

  _fileMesh = meshes->getMeshAtPos( meshIndex );
  _fileMesh->incrRef();

  if ( fields )
    for ( int i = 0; i < fields->getNumberOfFields(); ++i )
      {
        MEDFileAnyTypeFieldMultiTS * fB = fields->getFieldAtPos(i);
        MEDFileFieldMultiTS * f = dynamic_cast<MEDFileFieldMultiTS *>(fB);
        if(!f)
          continue;// fields on int32 not managed
        if ( f->getMeshName() == _fileMesh->getName() )
          {
            vector< vector<TypeOfField> > fTypes = f->getTypesOfFieldAvailable();
            if ( fTypes[0].size() == 1 && fTypes[0][0] == ON_NODES )
              _nodeFields.push_back( f );
            else
              _cellFields.push_back( f );
          }
      }
}

//================================================================================
/*!
 * \brief Adds a submesh
 */
//================================================================================

SauvWriter::SubMesh* SauvWriter::addSubMesh(const std::string& name, int dimRelExt)
{
  if ( _subs.capacity() < _subs.size() + 1 )
    THROW_IK_EXCEPTION("SauvWriter: INTERNAL error, wrong evaluation of nb of sub-meshes");
  _subs.resize( _subs.size() + 1 );
  SubMesh& sm = _subs.back();
  sm._name = name;
  sm._dimRelExt = dimRelExt;
  return &sm;
}
//================================================================================
/*!
 * \brief Returns nb of cell types
 */
//================================================================================

int SauvWriter::SubMesh::nbTypes() const
{
  int nb = 0;
  for (int i = 0; i < cellIDsByTypeSize(); ++i )
    nb += int( !_cellIDsByType[i].empty() );
  return nb;
}

//================================================================================
/*!
 * \brief Fill _subs
 */
//================================================================================

void SauvWriter::fillSubMeshes( int& nbSauvObjects, map<string,int>& nameNbMap )
{
  // evaluate nb of _subs in order to avoid re-allocation of _subs
  int nbSubs = 1; // for the very mesh
  nbSubs += _fileMesh->getFamilyInfo().size() + 4;  // + 4 zero families (for each dimRelExt)
  nbSubs += _fileMesh->getGroupInfo().size();
  nbSubs += evaluateNbProfileSubMeshes();
  _subs.clear();
  _subs.reserve( nbSubs );

  fillFamilySubMeshes();
  fillGroupSubMeshes();
  fillProfileSubMeshes();

  // fill names of SubMesh'es and count nb of sauv sub-meshes they will be stored into
  nbSauvObjects = 0;
  map<string,int> namePrefixMap;
  for ( size_t i = 0; i < _subs.size(); ++i )
    {
      SubMesh& sm = _subs[i];

      sm._nbSauvObjects = 0;
      if ( sm._subs.empty() )
        {
          sm._nbSauvObjects = sm.nbTypes();
        }
      else
        {
          sm._nbSauvObjects = 1;
        }

      sm._id = nbSauvObjects+1;
      nbSauvObjects += sm._nbSauvObjects;

      if ( sm._nbSauvObjects )
        sm._name = addName( nameNbMap, namePrefixMap, sm._name, sm._id );

      if ( sm._nbSauvObjects && !sm._name.empty() )
        {
          nameGIBItoMED aMEDName;
          aMEDName.gibi_pile = PILE_SOUS_MAILLAGE;
          aMEDName.gibi_id   = sm._id;
          aMEDName.med_name  = sm._name;
          _longNames[ LN_MAIL ].push_back(aMEDName);
        }
    }
}

//================================================================================
/*!
 * \brief fill sub-meshes of families
 */
//================================================================================

void SauvWriter::fillFamilySubMeshes()
{
  SubMesh* nilSm = (SubMesh*) 0;
  std::vector<int> dims = _fileMesh->getNonEmptyLevelsExt();
  for ( size_t iDim = 0; iDim < dims.size(); ++iDim )
    {
      int dimRelExt = dims[ iDim ];
      MCAuto< MEDCouplingMesh > mesh = _fileMesh->getMeshAtLevel(dimRelExt);
      const DataArrayInt * famIds = _fileMesh->getFamilyFieldAtLevel(dimRelExt);
      if ( !famIds ) continue;

      int curFamID = 0;
      SubMesh* curSubMesh = addSubMesh( "", dimRelExt ); // submesh of zero family
      _famIDs2Sub[0] = curSubMesh;
      int sub0Index = _subs.size()-1;

      const int * famID = famIds->begin(), * famIDEnd = famIds->end();
      for ( int cellID = 0; famID < famIDEnd; ++famID, cellID++ )
        {
          if ( *famID != curFamID )
            {
              curFamID = *famID;
              map< int, SubMesh* >::iterator f2s = _famIDs2Sub.insert( make_pair( curFamID, nilSm )).first;
              if ( !f2s->second )
                f2s->second = addSubMesh( "", dimRelExt ); // no names for families
              curSubMesh = f2s->second;
            }
          INTERP_KERNEL::NormalizedCellType cellType =
            dimRelExt == 1 ? INTERP_KERNEL::NORM_POINT1 : mesh->getTypeOfCell( cellID );
          curSubMesh->_cellIDsByType[ cellType ].push_back( cellID );
        }

      if ( dimRelExt == 1 )
        {
          // clear submesh of nodal zero family
          _famIDs2Sub[0]->_cellIDsByType[ INTERP_KERNEL::NORM_POINT1 ].clear();
        }
      else if ( dimRelExt == 0 )
        {
          // make a submesh including all cells
          if ( sub0Index == (int)(_subs.size()-1) )
            {
              _famIDs2Sub[0]->_name = _fileMesh->getName(); // there is the zero family only
            }
          else
            {
              curSubMesh = addSubMesh( _fileMesh->getName(), dimRelExt );
              if ( _famIDs2Sub[0]->nbTypes() == 0 )
                sub0Index++; // skip an empty zero family
              for ( size_t i = sub0Index; i < _subs.size()-1; ++i )
                curSubMesh->_subs.push_back( & _subs[i] );
            }
        }
    }
}

//================================================================================
/*!
 * \brief fill sub-meshes of groups
 */
//================================================================================

void SauvWriter::fillGroupSubMeshes()
{
  const map<string, vector<string> >& grpFams = _fileMesh->getGroupInfo();
  map<string, vector<string> >::const_iterator g2ff = grpFams.begin();
  for ( ; g2ff != grpFams.end(); ++g2ff )
    {
      const string&        groupName = g2ff->first;
      const vector<string>& famNames = g2ff->second;
      if ( famNames.empty() ) continue;
      std::vector<SubMesh*> famSubMeshes( famNames.size() );
      std::size_t k = 0;
      for ( size_t i = 0; i < famNames.size(); ++i )
        {
          int famID = _fileMesh->getFamilyId( famNames[i].c_str() );
          map< int, SubMesh* >::iterator i2f = _famIDs2Sub.find( famID );
          if ( i2f != _famIDs2Sub.end() )
            {
              famSubMeshes[ k ] = i2f->second;
              ++k;
            }
        }
      if ( k == 0 )
        continue;
      // if a family exists but has no element, no submesh has been found for this family
      // => we have to resize famSubMeshes with the number of submeshes stored
      if (k != famNames.size())
        famSubMeshes.resize(k);
      SubMesh* grpSubMesh = addSubMesh( groupName, famSubMeshes[0]->_dimRelExt );
      if( ! _cpy_grp_if_on_single_family )
        {
          grpSubMesh->_subs.swap( famSubMeshes );
        }
      else
        {
          /* If a group sub mesh consists of only one family, the group is written as
           * a copy of this family.
           * A mesh composed of only one submesh may cause an issue with some Gibi operators.*/
          if (famSubMeshes.size() == 1)
            {
              for(int i = 0; i < famSubMeshes[0]->cellIDsByTypeSize() ; i++)
                {
                  grpSubMesh->_cellIDsByType[i] = famSubMeshes[0]->_cellIDsByType[i];
                }
            }
          else
            grpSubMesh->_subs.swap( famSubMeshes );
        }
    }
}


//================================================================================
/*!
 * \brief fill sub-meshes of profiles
 */
//================================================================================

void SauvWriter::fillProfileSubMeshes()
{
  _profile2Sub.clear();
  SubMesh* nilSm = (SubMesh*) 0;
  for ( int isOnNodes = 0; isOnNodes < 2; ++isOnNodes )
    {
      vector< MCAuto< MEDFileFieldMultiTS > >
        fields = isOnNodes ? _nodeFields : _cellFields;
      for ( size_t i = 0; i < fields.size(); ++i )
        {
          vector< pair<int,int> > iters = fields[i]->getIterations();

          vector<INTERP_KERNEL::NormalizedCellType> types;
          vector< vector<TypeOfField> > typesF;
          vector< vector<string> > pfls, locs;
          fields[i]->getFieldSplitedByType( iters[0].first, iters[0].second,
                                            _fileMesh->getName().c_str(), types, typesF, pfls, locs);
          int dimRelExt;
          for ( size_t iType = 0; iType < types.size(); ++iType )
            {
              if ( types[iType] == INTERP_KERNEL::NORM_ERROR )
                dimRelExt = 1; // on nodes
              else
                dimRelExt = getDimension( types[iType] ) - _fileMesh->getMeshDimension();
              for ( size_t iPfl = 0; iPfl < pfls[iType].size(); ++iPfl )
                {
                  bool isOnAll = pfls[iType][iPfl].empty();
                  if ( isOnAll ) pfls[iType][iPfl] = noProfileName( types[iType] );
                  map< string, SubMesh* >::iterator pfl2sm =
                    _profile2Sub.insert( make_pair( pfls[iType][iPfl], nilSm )).first;
                  if ( !pfl2sm->second )
                    {
                      SubMesh* sm = pfl2sm->second = addSubMesh( "", dimRelExt ); // no names for profiles
                      const DataArrayInt * pfl = isOnAll ? 0 : fields[i]->getProfile( pfls[iType][iPfl].c_str() );
                      makeProfileIDs( sm, types[iType], pfl );
                    }
                }
            }
        }
    }
}

//================================================================================
/*!
 * \brief Return max possible nb of sub-meshes to decsribe field supports
 */
//================================================================================

int SauvWriter::evaluateNbProfileSubMeshes() const
{
  int nb = 0;
  for ( size_t i = 0; i < _nodeFields.size(); ++i )
    nb += 1 + _nodeFields[i]->getPflsReallyUsed().size();

  for ( size_t i = 0; i < _cellFields.size(); ++i )
    {
      nb += _cellFields[i]->getPflsReallyUsed().size();

      vector< pair<int,int> > iters = _cellFields[i]->getIterations();

      vector<INTERP_KERNEL::NormalizedCellType> types;
      vector< vector<TypeOfField> > typesF;
      vector< vector<string> > pfls, locs;
      _cellFields[i]->getFieldSplitedByType( iters[0].first, iters[0].second,
                                             _fileMesh->getName().c_str(), types, typesF, pfls, locs);
      nb += 2 * types.size(); // x 2 - a type can be on nodes and on cells at the same time
    }

  return nb;
}

//================================================================================
/*!
 * \brief Transorm a profile into ids of mesh elements
 */
//================================================================================

void SauvWriter::makeProfileIDs( SubMesh*                          sm,
                                 INTERP_KERNEL::NormalizedCellType type,
                                 const DataArrayInt*               profile )
{
  MCAuto< MEDCouplingMesh >
    mesh = _fileMesh->getMeshAtLevel(sm->_dimRelExt);
  const MEDCouplingUMesh* uMesh = dynamic_cast< const MEDCouplingUMesh* > ((const MEDCouplingMesh*) mesh );

  if ( sm->_dimRelExt == 1 ) type = INTERP_KERNEL::NORM_POINT1;
  vector< int >& ids = sm->_cellIDsByType[ type ];

  if ( sm->_dimRelExt == 1 || !uMesh )
    {
      // profile on nodes or mesh is CARTESIAN
      if ( profile )
        {
          ids.assign( profile->begin(), profile->end() );
        }
      else // on all
        {
          ids.resize( sm->_dimRelExt == 1 ? mesh->getNumberOfNodes() : mesh->getNumberOfCells() );
          for ( size_t i = 0; i < ids.size(); ++i )
            ids[i]=i;
        }
    }
  else
    {
      // profile on cells
      vector<int> code(3);
      code[0] = type;
      if ( profile ) // on profile
        {
          code[1] = profile->getNumberOfTuples();
          code[2] = 0;
        }
      else // on all cells
        {
          code[1] = mesh->getNumberOfCellsWithType( type );
          code[2] = -1;
        }
      vector<const DataArrayInt *> idsPerType( 1, profile );
      MCAuto<DataArrayInt>
        resIDs = uMesh->checkTypeConsistencyAndContig( code, idsPerType );
      if (( const DataArrayInt *) resIDs )
      {
        ids.assign( resIDs->begin(), resIDs->end() );
      }
      else // mesh includes only one type
      {
        int nbE = code[1];
        for ( ids.resize( nbE ); nbE; --nbE )
          ids[ nbE-1 ] = nbE-1;
      }
    }
}

//================================================================================
/*!
 * \brief Write its data into the SAUVE file
 */
//================================================================================

void SauvWriter::write(const std::string& fileName)
{
  std::fstream fileStream;
  fileStream.open( fileName.c_str(), ios::out);
  if
#ifdef WIN32
    ( !fileStream || !fileStream.is_open() )
#else
    ( !fileStream || !fileStream.rdbuf()->is_open() )
#endif
      THROW_IK_EXCEPTION("Can't open the file |"<<fileName<<"|");
  _sauvFile = &fileStream;

  _subs.clear();
  _famIDs2Sub.clear();
  _profile2Sub.clear();
  _longNames[ LN_MAIL ].clear();
  _longNames[ LN_CHAM ].clear();
  _longNames[ LN_COMP ].clear();

  map<string,int> fldNamePrefixMap;

  writeFileHead();
  writeSubMeshes();
  writeNodes();
  writeNodalFields(fldNamePrefixMap);
  writeElemFields(fldNamePrefixMap);
  writeLongNames();
  writeLastRecord();

  _sauvFile->close();
}
//================================================================================
/*!
 * \brief Writes "ENREGISTREMENT DE TYPE" 4 and 7
 */
//================================================================================

void SauvWriter::writeFileHead()
{
  MCAuto< MEDCouplingMesh > mesh = _fileMesh->getMeshAtLevel(0);

  *_sauvFile
    << " ENREGISTREMENT DE TYPE   4" << endl
    << " NIVEAU  16 NIVEAU ERREUR   0 DIMENSION   " << mesh->getSpaceDimension() <<endl
    << " DENSITE 0.00000E+00" << endl
    << " ENREGISTREMENT DE TYPE   7" << endl
    << " NOMBRE INFO CASTEM2000   8" <<endl
    << " IFOUR  -1 NIFOUR   0 IFOMOD  -1 IECHO   1 IIMPI   0 IOSPI   0 ISOTYP   1" << endl
    << " NSDPGE     0" << endl;
}

//================================================================================
/*!
 * \brief Writes names of objects
 */
//================================================================================

void SauvWriter::writeNames( const map<string,int>& nameNbMap )
{
  if ( !nameNbMap.empty() )
  {
    // write names of objects
    // * 8001       FORMAT(8(1X,A8))
    TFieldCounter fcount( *_sauvFile, 8 );
    *_sauvFile << left;
    map<string,int>::const_iterator nameNbIt = nameNbMap.begin();
    for ( ; nameNbIt != nameNbMap.end(); nameNbIt++, fcount++ )
      *_sauvFile << " " << setw(8) << nameNbIt->first;
    fcount.stop();
    *_sauvFile << right;

    // write IDs of named objects in the pile
    // *  8000 FORMAT(10I8)
    nameNbIt = nameNbMap.begin();
    for ( fcount.init(10); nameNbIt != nameNbMap.end(); nameNbIt++, fcount++ )
      *_sauvFile << setw(8) << nameNbIt->second;
  }
}

//================================================================================
/*!
 * \brief Writes "PILE NUMERO   1"
 */
//================================================================================

void SauvWriter::writeSubMeshes()
{
  int nbSauvObjects;
  map<string,int> nameNbMap;
  fillSubMeshes( nbSauvObjects, nameNbMap );

  // * 800   FORMAT (' ENREGISTREMENT DE TYPE', I4)
  *_sauvFile << " ENREGISTREMENT DE TYPE   2" << endl;
  // * 801     FORMAT(' PILE NUMERO',I4,'NBRE OBJETS NOMMES',I8,'NBRE OBJETS',I8)
  *_sauvFile << " PILE NUMERO   1NBRE OBJETS NOMMES" << setw(8) << nameNbMap.size() <<
    "NBRE OBJETS" << setw(8) << nbSauvObjects <<endl;

  writeNames( nameNbMap );

  TFieldCounter fcount( *_sauvFile, 10 ); // 10 intergers per line

  for ( size_t iSub = 0; iSub < _subs.size(); ++iSub )
    {
      SubMesh& sm = _subs[iSub];
      if ( sm._nbSauvObjects < 1 ) continue;

      // The first record of each sub-mesh writes
      // - type of cells; zero means a compound object whose the 2nd record enumerates its components
      // - number of components of a compound object
      // - number of references; each reference means a "pointer" to this sub-mesh
      // - number of nodes per cell
      // - number of cells

      if ( !sm._subs.empty() )
        {
          writeCompoundSubMesh(iSub);
        }
      else
        {
          // write each sub-type as a SAUV sub-mesh
          MCAuto< MEDCouplingMesh >
            mesh = _fileMesh->getMeshAtLevel( sm._dimRelExt );
          MCAuto< MEDCouplingUMesh>
            umesh = mesh->buildUnstructured();

          for ( int iType=0; iType < sm.cellIDsByTypeSize(); ++iType )
            {
              const vector<int>& cellIDs = sm._cellIDsByType[iType];
              if ( cellIDs.empty() ) continue;

              INTERP_KERNEL::NormalizedCellType
                cellType = INTERP_KERNEL::NormalizedCellType( iType );
              const INTERP_KERNEL::CellModel &
                cell = INTERP_KERNEL::CellModel::GetCellModel( cellType );
              int castemType       = SauvUtilities::med2gibiGeom( cellType );
              unsigned nbElemNodes = cell.getNumberOfNodes();
              unsigned nbElems     = cellIDs.size();

              *_sauvFile << setw(8) << castemType
                        << zeroI8
                        << zeroI8
                        << setw(8) << nbElemNodes
                        << setw(8) << nbElems << endl;

              // write color of each element
              // * 8000 FORMAT(10I8)
              for ( size_t i = 0; i < nbElems; ++i, fcount++ ) *_sauvFile << zeroI8;
              fcount.stop();

              // write connectivity
              // gibi IDs are in FORTRAN mode while MEDCoupling IDs are in C mode
              if ( sm._dimRelExt == 1 ) // nodes
                {
                  for ( size_t i = 0; i < nbElems; ++i, fcount++ )
                    *_sauvFile << setw(8) << ( cellIDs[i] + 1 );
                }
              else
                {
                  // indices to transform MED connectivity to GIBI one
                  const int * toMedConn = getGibi2MedQuadraticInterlace( cellType );

                  vector< int > cellConn( nbElemNodes ), transformedConn( nbElemNodes );
                  for ( size_t i = 0; i < nbElems; ++i )
                    {
                      cellConn.clear();
                      umesh->getNodeIdsOfCell( cellIDs[i], cellConn );
                      if ( toMedConn )
                        {
                          for ( unsigned j = 0; j < nbElemNodes; ++j )
                            transformedConn[ toMedConn[ j ]] = cellConn[ j ];
                          cellConn.swap( transformedConn );
                        }
                      for ( unsigned j = 0; j < nbElemNodes; ++j, fcount++ )
                        *_sauvFile << setw(8) << ( cellConn[j] + 1 );
                    }
                }
              fcount.stop();

            } // loop on cell types
        } // not a compound object
    } // loop on sub-meshes
}

//================================================================================
/*!
 * \brief Writes a sum-mesh composed of other sum-meshes
 * This submesh corresponds to a med mesh or group composed of families
 */
//================================================================================

void SauvWriter::writeCompoundSubMesh(int iSub)
{
  SubMesh& sm = _subs[iSub];
  if ( sm._nbSauvObjects < 1 || sm._subs.empty()) return;

  vector< int > subIDs;
  for ( size_t i = 0; i < sm._subs.size(); ++i ) // loop on sub-meshes of families
    for ( int j = 0; j < sm._subs[i]->_nbSauvObjects; ++j )
      subIDs.push_back( sm._subs[i]->_id + j );
      
  *_sauvFile << zeroI8
             << setw(8) << subIDs.size()
             << zeroI8
             << zeroI8
             << zeroI8 << endl;

  TFieldCounter fcount( *_sauvFile, 10 ); // 10 intergers per line
  for ( size_t i = 0; i < subIDs.size(); ++i, fcount++ )
    *_sauvFile << setw(8) << subIDs[i];
}

//================================================================================
/*!
 * \brief Write piles relating to nodes
 */
//================================================================================

void SauvWriter::writeNodes()
{
  MCAuto< MEDCouplingMesh > mesh = _fileMesh->getMeshAtLevel( 1 );
  MCAuto< MEDCouplingUMesh > umesh = mesh->buildUnstructured();

  // write the index connecting nodes with their coodrinates

  const int nbNodes = umesh->getNumberOfNodes();
  *_sauvFile << " ENREGISTREMENT DE TYPE   2" << endl
             << " PILE NUMERO  32NBRE OBJETS NOMMES       0NBRE OBJETS" << setw(8) << nbNodes << endl;
  *_sauvFile << setw(8) << nbNodes << endl;
  //
  TFieldCounter fcount( *_sauvFile, 10 );// * 8000 FORMAT(10I8)
  for ( int i = 0; i < nbNodes; ++i, fcount++ )
    *_sauvFile << setw(8) << i + 1; 
  fcount.stop();

  // write coordinates and density of nodes

  *_sauvFile << " ENREGISTREMENT DE TYPE   2" << endl;
  *_sauvFile << " PILE NUMERO  33NBRE OBJETS NOMMES       0NBRE OBJETS       1" << endl;
  // 
  const int dim = umesh->getSpaceDimension();
  const int nbValues = nbNodes * ( dim + 1 );
  *_sauvFile << setw(8) << nbValues << endl;

  // * 8003   FORMAT(1P,3E22.14)
  const char* density = "  0.00000000000000E+00";
  fcount.init(3);
  _sauvFile->precision(14);
  _sauvFile->setf( ios_base::scientific, ios_base::floatfield );
  _sauvFile->setf( ios_base::uppercase );
  MCAuto< DataArrayDouble> coordArray = umesh->getCoordinatesAndOwner();
  const double precision = 1.e-99; // PAL12077
  for ( int i = 0; i < nbNodes; ++i)
  {
    for ( int j = 0; j < dim; ++j, fcount++ )
      {
        double coo = coordArray->getIJ( i, j );
        bool  zero = ( -precision < coo && coo < precision );
        *_sauvFile << setw(22) << ( zero ? 0.0 : coo );
      }
    *_sauvFile << density;
    fcount++;
  }
}

//================================================================================
/*!
 * \brief Store correspondence between GIBI (short) and MED (long) names
 *
 * IMP 0020434: mapping GIBI names to MED names
 * Store correspondence between GIBI and MED names as one PILE_STRINGS and one
 * PILE_TABLES (in three tables: MED_MAIL, MED_CHAM and MED_COMP)
 */
//================================================================================

void SauvWriter::writeLongNames()
{
  int nbTables =
    3 - _longNames[ LN_MAIL ].empty() - _longNames[ LN_CHAM ].empty() - _longNames[ LN_COMP ].empty();
  if (nbTables == 0) return;

  // ---------------------
  // Write the TABLE pile
  // ---------------------

  *_sauvFile << " ENREGISTREMENT DE TYPE   2" << endl
        << " PILE NUMERO  10NBRE OBJETS NOMMES" << setw(8) << nbTables
        << "NBRE OBJETS" << setw(8) << nbTables << endl;
  // table names
  if (!_longNames[ LN_MAIL ].empty()) *_sauvFile << " MED_MAIL";
  if (!_longNames[ LN_CHAM ].empty()) *_sauvFile << " MED_CHAM";
  if (!_longNames[ LN_COMP ].empty()) *_sauvFile << " MED_COMP";
  *_sauvFile << endl;
  // table indices
  for ( int i = 0; i < nbTables; ++i ) *_sauvFile << setw(8) << i+1;
  *_sauvFile << endl;

  string theWholeString; // concatenated long names
  vector<int> theOffsets;
  int iStr = 1;
  TFieldCounter fcount (*_sauvFile, 10);

  for ( int iTbl = 0; iTbl < LN_NB; ++iTbl )
    {
      vector<nameGIBItoMED>& longNames = _longNames[ iTbl ];
      if ( longNames.empty() ) continue;
      const bool isComp = ( iTbl == LN_COMP);

      // to assure unique MED names
      set<string> medUniqueNames;

      *_sauvFile << setw(8) << longNames.size()*4 << endl; // Nb of table values

      vector<nameGIBItoMED>::iterator itGIBItoMED = longNames.begin();
      for (; itGIBItoMED != longNames.end(); itGIBItoMED++, iStr++)
        {
          // PILE of i-th key (med name)
          *_sauvFile << setw(8) << PILE_STRINGS;
          fcount++;
          // ID of i-th key (med name)
          *_sauvFile << setw(8) << iStr;
          fcount++;
          // PILE of i-th value (gibi name)
          *_sauvFile << setw(8) << itGIBItoMED->gibi_pile;
          fcount++;
          // ID of i-th value (gibi name)
          *_sauvFile << setw(8) << ( isComp ? ++iStr : itGIBItoMED->gibi_id );
          fcount++;

          // add a MED name to the string (while making it be unique for sub-meshes and fields)
          string aMedName = itGIBItoMED->med_name;
          if ( !isComp )
            for (int ind = 1; !medUniqueNames.insert(aMedName).second; ++ind )
              aMedName = itGIBItoMED->med_name + "_" + SauvUtilities::toString( ind );
          theWholeString += aMedName;

          // add an offset
          theOffsets.push_back( theWholeString.size() );
          if ( isComp )
            {
              theWholeString += itGIBItoMED->gibi_name;
              theOffsets.push_back( theWholeString.size() );
            }
        }
      fcount.stop();
    }

  // ----------------------
  // Write the STRING pile
  // ----------------------

  const int nbNames = theOffsets.size();
  *_sauvFile << " ENREGISTREMENT DE TYPE   2" << endl
        << " PILE NUMERO  27NBRE OBJETS NOMMES" << zeroI8 << "NBRE OBJETS" << setw(8) << nbNames << endl
        << setw(8) << theWholeString.length() << setw(8) << nbNames << endl;

  // write the whole string
  const int fixedLength = 71;
  for ( string::size_type aPos = 0; aPos < theWholeString.length(); aPos += fixedLength)
    *_sauvFile << setw(72) << theWholeString.substr(aPos, fixedLength) << endl;

  // write the offsets
  for ( size_t i = 0; i < theOffsets.size(); ++i, fcount++ )
    *_sauvFile << setw(8) << theOffsets[i];
}

//================================================================================
/*!
 * \brief Write beginning of field record
 */
//================================================================================

void SauvWriter::writeFieldNames( const bool                 isNodal,
                                  std::map<std::string,int>& fldNamePrefixMap)
{
  vector< MCAuto< MEDFileFieldMultiTS > >&
    flds = isNodal ? _nodeFields : _cellFields;
  map<string,int> nameNbMap;

  for ( size_t iF = 0; iF < flds.size(); ++iF )
    {
      string name = addName( nameNbMap, fldNamePrefixMap, flds[iF]->getName(), iF+1 );
      nameGIBItoMED aMEDName;
      aMEDName.gibi_pile = isNodal ? PILE_NODES_FIELD : PILE_FIELD;
      aMEDName.gibi_id   = iF+1;
      aMEDName.med_name  = name;
      _longNames[ LN_CHAM ].push_back(aMEDName);
    }

  *_sauvFile << " ENREGISTREMENT DE TYPE   2" << endl
             << ( isNodal ? " PILE NUMERO   2" : " PILE NUMERO  39")
             << "NBRE OBJETS NOMMES" << setw(8) << nameNbMap.size()
             << "NBRE OBJETS"        << setw(8) << flds.size() << endl;
  writeNames( nameNbMap );
}

//================================================================================
/*!
 * \brief Make short names of field components
 *
 * IMP 0020434: mapping GIBI names to MED names
 */
//================================================================================

void SauvWriter::makeCompNames(const string&         fieldName,
                               const vector<string>& compInfo,
                               map<string, string>&  mapMedToGibi)
{
  for ( size_t i = 0; i < compInfo.size(); ++i )
    mapMedToGibi[compInfo[i]] = cleanName( compInfo[i] );

  int compIndex = 1;
  map<string, string>::iterator namesIt = mapMedToGibi.begin();
  for (; namesIt != mapMedToGibi.end(); namesIt++)
    {
      string & compGibiName = (*namesIt).second;
      if (compGibiName.size() > 4) {
        // use new name in form "CXXX", where "XXX" is a number
        do
          {
            compGibiName = SauvUtilities::toString( compIndex++ );
            if ( compGibiName.size() < 3 )
              compGibiName.insert( 0, 3 - compGibiName.size(), '0' );
            compGibiName = "C" + compGibiName;
          }
        while (mapMedToGibi.count(compGibiName) > 0); // real component name could be CXXX
      }

      string compMedName = fieldName + "." + namesIt->first;
      nameGIBItoMED aMEDName;
      aMEDName.med_name  = compMedName;
      aMEDName.gibi_pile = PILE_STRINGS;
      aMEDName.gibi_name = compGibiName;
      _longNames[ LN_COMP ].push_back(aMEDName);
    }
}

//================================================================================
/*!
 * \brief Writes "PILE NUMERO   2": fields on nodes
 */
//================================================================================

void SauvWriter::writeNodalFields(map<string,int>& fldNamePrefixMap)
{
  writeFieldNames( /*isNodal=*/true, fldNamePrefixMap );

  TFieldCounter fcount (*_sauvFile, 10);

  // EXAMPLE ( with no values )

  // (1)       4       7       2       1
  // (2)     -88       0       3     -89       0       1     -90       0       2     -91
  // (2)       0       1
  // (3) FX   FY   FZ   FZ   FX   FY   FLX
  // (4)       0       0       0       0       0       0       0
  // (5)           cree  par  muc pri
  // (6)
  // (7)       2
  for ( size_t iF = 0; iF < _nodeFields.size(); ++iF )
    {
      // (1) write nb subcomponents, nb components(total)
      vector< pair<int,int> >  iters = _nodeFields[iF]->getIterations();
      const vector<string>& compInfo = _nodeFields[iF]->getInfo();
      const int nbSub = iters.size();
      const int nbComp = compInfo.size();
      const int totalNbComp = nbSub * nbComp;
      *_sauvFile << setw(8) << nbSub
                 << setw(8) << totalNbComp
                 << setw(8) << -1         // IFOUR
                 << setw(8) << 0 << endl; // nb attributes

      // (2) for each sub-component (iteration)
      // write support, number of values and number of components
      fcount.init(10);
      vector< int > vals(3);
      for ( std::size_t iIt = 0; iIt < iters.size(); ++iIt )
        {
          pair<int,int> it = iters[iIt];

          vector<INTERP_KERNEL::NormalizedCellType> types;
          vector< vector<TypeOfField> > typesF;
          vector< vector<string> > pfls, locs;
          vector< vector< std::pair<int,int> > > valsVec;
          valsVec=_nodeFields[iF]->getFieldSplitedByType( it.first, it.second, _fileMesh->getName().c_str(),
                                                          types, typesF, pfls, locs);
          // believe that there can be only one type in a nodal field,
          // so do not use a loop on types
          if ( pfls[0][0].empty() ) pfls[0][0] = noProfileName( types[0] );
          map< string, SubMesh* >::iterator pfl2Sub = _profile2Sub.find( pfls[0][0] );
          if ( pfl2Sub == _profile2Sub.end() )
            THROW_IK_EXCEPTION( "SauvWriter::writeNodalFields(): no sub-mesh for profile |"
                                << pfls[0][0] << "|");
          vals[0] = -pfl2Sub->second->_id;
          vals[1] = (valsVec[0][0].second-valsVec[0][0].first);
          vals[2] = compInfo.size();
          for ( size_t i = 0; i < vals.size(); ++i, fcount++ )
            *_sauvFile << setw(8) << vals[i];
        }
      fcount.stop();

      // (3) Write names of components
      map<string, string> mapMedToGibi;
      makeCompNames( _nodeFields[iF]->getName(), compInfo, mapMedToGibi );
      fcount.init(8);
      *_sauvFile << left;
      for ( std::size_t iIt = 0; iIt < iters.size(); ++iIt )
        for ( size_t i = 0; i < compInfo.size(); ++i, fcount++ )
          *_sauvFile << " "  << setw(4) << mapMedToGibi[compInfo[i]];
      *_sauvFile << right;
      fcount.stop();

      // (4) nb harmonics
      fcount.init(10);
      for ( size_t i = 0; i < (std::size_t)totalNbComp; ++i, fcount++ )
        *_sauvFile << " "  << setw(8) << 0;
      fcount.stop();

      string description = _nodeFields[iF]->getName();
      *_sauvFile << endl;                                         // (5) TYPE
      *_sauvFile << setw(72) << description.substr(0,71) << endl; // (6) TITRE
      //*_sauvFile << endl;                                         // (7) 0 attributes

      // write values of each component
      fcount.init( 3 ); // 3 values per a line
      for ( std::size_t iIt = 0; iIt < iters.size(); ++iIt )
        {
          pair<int,int> it = iters[iIt];

          vector<INTERP_KERNEL::NormalizedCellType> types;
          vector< vector<TypeOfField> > typesF;
          vector< vector<string> > pfls, locs;
          vector< vector< std::pair<int,int> > > valsVec;
          valsVec = _nodeFields[iF]->getFieldSplitedByType( it.first, it.second, _fileMesh->getName().c_str(),
                                                            types, typesF, pfls, locs);
          // believe that there can be only one type in a nodal field,
          // so do not perform a loop on types
          const DataArrayDouble* valsArray = _nodeFields[iF]->getUndergroundDataArray(it.first, it.second);
          for ( size_t j = 0; j < compInfo.size(); ++j )
            {
              for ( size_t i = valsVec[0][0].first; i < (std::size_t)valsVec[0][0].second; ++i, fcount++ )
                *_sauvFile << setw(22) << valsArray->getIJ( i, j );
              fcount.stop();
            }
        }
    } // loop on fiels
}

//================================================================================
/*!
 * \brief Writes "PILE NUMERO  39": fields on cells
 */
//================================================================================

void SauvWriter::writeElemFields(map<string,int>& fldNamePrefixMap)
{
  writeFieldNames( /*isNodal=*/false, fldNamePrefixMap );

  TFieldCounter fcount (*_sauvFile, 10);

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
  // (9)       1       1       0       0
  // (10)  3.30000000000000E-01
  // (9)       1       1       0       0
  // (10)  1.00000000000000E+04
  // (9)       6     706       0       0
  // (10)  1.00000000000000E+02  1.00000000000000E+02  1.00000000000000E+02
  // (10)  1.00000000000000E+02  1.00000000000000E+02  1.00000000000000E+02
  // (10)  ...

  for ( size_t iF = 0; iF < _cellFields.size(); ++iF )
    {
      // count nb of sub-components
      int iSub, nbSub = 0;
      vector< pair<int,int> >  iters = _cellFields[iF]->getIterations();
      for ( std::size_t iIt = 0; iIt < iters.size(); ++iIt )
        {
          pair<int,int> it = iters[iIt];

          vector<INTERP_KERNEL::NormalizedCellType> types;
          vector< vector<TypeOfField> > typesF;
          vector< vector<string> > pfls, locs;
          vector< vector< std::pair<int,int> > > valsVec;
          valsVec = _cellFields[iF]->getFieldSplitedByType( it.first, it.second, _fileMesh->getName().c_str(),
                                                            types, typesF, pfls, locs);
          for ( size_t i = 0; i < valsVec.size(); ++i )
            nbSub += valsVec[i].size();
        }
      // (1) write nb sub-components, title length
      *_sauvFile << setw(8) << nbSub
                 << setw(8) << -1 // whatever
                 << setw(8) << 6  // whatever
                 << setw(8) << 72 << endl; // title length
      // (2) title
      string title = _cellFields[iF]->getName();
      *_sauvFile << setw(72) << title.substr(0,71) << endl;
      *_sauvFile << setw(72) << " " << endl;

      // (3) support, nb components
      vector<int> vals(9, 0);
      const vector<string>& compInfo = _cellFields[iF]->getInfo();
      vals[2] = compInfo.size();
      fcount.init(10);
      for ( std::size_t iIt = 0; iIt < iters.size(); ++iIt )
        {
          pair<int,int> it = iters[iIt];

          vector<INTERP_KERNEL::NormalizedCellType> types;
          vector< vector<TypeOfField> > typesF;
          vector< vector<string> > pfls, locs;
          _cellFields[iF]->getFieldSplitedByType( it.first, it.second, _fileMesh->getName().c_str(),
                                                  types, typesF, pfls, locs);
          for ( size_t iType = 0; iType < pfls.size(); ++iType )
            for ( size_t iP = 0; iP < pfls[iType].size(); ++iP )
              {
                if ( pfls[iType][iP].empty() ) pfls[iType][iP] = noProfileName( types[iType] );
                map< string, SubMesh* >::iterator pfl2Sub = _profile2Sub.find( pfls[iType][iP] );
                if ( pfl2Sub == _profile2Sub.end() )
                  THROW_IK_EXCEPTION( "SauvWriter::writeElemFields(): no sub-mesh for profile |"
                                      << pfls[iType][iP] << "|");
                const int supportID = pfl2Sub->second->_id;
                vals[0] = -supportID;

                for ( size_t i = 0; i < vals.size(); ++i, fcount++ )
                  *_sauvFile << setw(8) << vals[ i ];
              }
        }
      fcount.stop();

      // (4) dummy strings
      for ( fcount.init(4), iSub = 0; iSub < nbSub; ++iSub, fcount++ )
        *_sauvFile << "                  ";
      fcount.stop();

      // (5) dummy strings
      for ( fcount.init(8), iSub = 0; iSub < nbSub; ++iSub, fcount++ )
        *_sauvFile << "         ";
      fcount.stop();

      // loop on sub-components of a field, each of which refers to
      // a certain support and has its own number of components
      for ( std::size_t iIt = 0; iIt < iters.size(); ++iIt )
        {
          pair<int,int> it = iters[iIt];
          writeElemTimeStamp( iF, it.first, it.second );
        }
    } // loop on cell fields
}

//================================================================================
/*!
 * \brief Write one elemental time stamp
 */
//================================================================================

void SauvWriter::writeElemTimeStamp(int iF, int iter, int order)
{
  // (6)   317767  317761  317755  317815
  // (7)  YOUN     NU       H        SIGY
  // (8)  REAL*8            REAL*8            REAL*8            REAL*8
  // (9)        1       1       0       0
  // (10)  2.00000000000000E+05
  // (9)       1       1       0       0
  // (10)  3.30000000000000E-01
  // (9)       1       1       0       0
  // (10)  1.00000000000000E+04
  // (9)       6     706       0       0
  // (10)  1.00000000000000E+02  1.00000000000000E+02  1.00000000000000E+02
  // (10)  1.00000000000000E+02  1.00000000000000E+02  1.00000000000000E+02

  TFieldCounter fcount (*_sauvFile, 10);

  vector<INTERP_KERNEL::NormalizedCellType> types;
  vector< vector<TypeOfField> > typesF;
  vector< vector<string> > pfls, locs;
  vector< vector< std::pair<int,int> > > valsVec;
  valsVec = _cellFields[iF]->getFieldSplitedByType( iter, order, _fileMesh->getName().c_str(),
                                                    types, typesF, pfls, locs);
  for ( size_t iType = 0; iType < pfls.size(); ++iType )
    for ( size_t iP = 0; iP < pfls[iType].size(); ++iP )
      {
        const vector<string>& compInfo = _cellFields[iF]->getInfo();

        // (6) component addresses
        int iComp = 0, nbComp = compInfo.size();
        for ( fcount.init(10); iComp < nbComp; ++iComp, fcount++ )
          *_sauvFile << setw(8) << 777; // a good number
        fcount.stop();

        // (7) component names
        map<string, string> mapMedToGibi;
        makeCompNames( _cellFields[iF]->getName(), compInfo, mapMedToGibi );
        *_sauvFile << left;
        for ( fcount.init(8), iComp = 0; iComp < nbComp; ++iComp, fcount++ )
          *_sauvFile << " "  << setw(8) << mapMedToGibi[compInfo[iComp]];
        fcount.stop();

        // (8) component types
        for ( fcount.init(4), iComp = 0; iComp < nbComp; ++iComp, fcount++ )
          *_sauvFile << " "  << setw(17) << "REAL*8";
        fcount.stop();
        *_sauvFile << right;

        // (9) nb values per element, nb of elements
        int nbPntPerCell = 1;
        if ( !locs[iType][iP].empty() )
          {
            int locID = _cellFields[iF]->getLocalizationId( locs[iType][iP].c_str() );
            nbPntPerCell = _cellFields[iF]->getNbOfGaussPtPerCell( locID );
          }
        else if ( typesF[iType][iP] == ON_GAUSS_NE )
          {
            nbPntPerCell = INTERP_KERNEL::CellModel::GetCellModel(types[iType]).getNumberOfNodes();
          }

        // (10) values
        const std::pair<int,int>& bgEnd = valsVec[iType][iP];
        const DataArrayDouble* valArray = _cellFields[iF]->getUndergroundDataArray(iter, order);
        for ( iComp = 0; iComp < nbComp; ++iComp )
          {
            *_sauvFile << setw(8) << nbPntPerCell
                       << setw(8) << (bgEnd.second-bgEnd.first) / nbPntPerCell
                       << setw(8) << 0
                       << setw(8) << 0
                       << endl;
            fcount.init(3);
            for ( size_t i = bgEnd.first; i < (size_t) bgEnd.second; ++i, fcount++ )
              *_sauvFile << setw(22) << valArray->getIJ( i, iComp );
            fcount.stop();
          }
      }
}

//================================================================================
/*!
 * \brief Write the last record of the SAUV file
 */
//================================================================================

void SauvWriter::writeLastRecord()
{
  *_sauvFile << " ENREGISTREMENT DE TYPE   5" << endl;
  *_sauvFile << "LABEL AUTOMATIQUE :   1" << endl;
}
