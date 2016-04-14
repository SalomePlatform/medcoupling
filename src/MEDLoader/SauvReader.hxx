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
// File      : SauvReader.hxx
// Created   : Tue Aug 16 13:04:25 2011
// Author    : Edward AGAPOV (eap)
//
#ifndef __SAUVREADER_HXX__
#define __SAUVREADER_HXX__

#include "MEDLoaderDefines.hxx"
#include "InterpKernelException.hxx"
#include "SauvUtilities.hxx"
#include "MEDCouplingRefCountObject.hxx"

#include <vector>
#include <string>
#include <set>

namespace SauvUtilities
{
  class FileReader;
  struct IntermediateMED;
  struct Group;
  struct DoubleField;
}
namespace MEDCoupling
{
  class MEDFileData;

class SauvReader : public MEDCoupling::RefCountObject
{
 public:
  MEDLOADER_EXPORT static SauvReader* New(const std::string& fileName);
  MEDLOADER_EXPORT MEDCoupling::MEDFileData * loadInMEDFileDS();
  MEDLOADER_EXPORT ~SauvReader();

 private:
  std::size_t getHeapMemorySizeWithoutChildren() const;
  std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
  void readRecord2();
  void readRecord4();
  void readRecord7();

  int readPileNumber(int& nbNamedObjects, int& nbObjects);
  void read_PILE_SOUS_MAILLAGE(const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_NODES_FIELD  (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_TABLES       (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_LREEL        (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_LOGIQUES     (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_FLOATS       (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_INTEGERS     (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_STRINGS      (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_LMOTS        (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_NOEUDS       (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_COORDONNEES  (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_MODL         (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);
  void read_PILE_FIELD        (const int nbObjects, std::vector<std::string>& objectNames, std::vector<int>& nameIndices);

  void setFieldSupport(const std::vector<SauvUtilities::Group*>& supports,
                       SauvUtilities::DoubleField*               field);
  void setFieldNames(const std::vector<SauvUtilities::DoubleField*>& fields,
                     const std::vector<std::string>& objectNames,
                     const std::vector<int>& nameIndices);

  bool isASCII() const                                   { return _fileReader->isASCII(); }
  bool isXRD() const                                     { return !isASCII(); }
  bool getNextLine (char* & line, bool raiseOEF = true ) { return _fileReader->getNextLine( line, raiseOEF ); }
  void initNameReading(int nbValues, int width = 8)      { _fileReader->initNameReading( nbValues, width ); }
  void initIntReading(int nbValues)                      { _fileReader->initIntReading( nbValues ); }
  void initDoubleReading(int nbValues)                   { _fileReader->initDoubleReading( nbValues ); }
  bool more() const                                      { return _fileReader->more(); }
  void next()                                            { _fileReader->next(); }
  int  index() const                                     { return _fileReader->index(); }
  int    getInt() const                                  { return _fileReader->getInt(); }
  int    getIntNext()                                    { int i = getInt(); next(); return i; }
  float  getFloat() const                                { return _fileReader->getFloat(); }
  double getDouble() const                               { return _fileReader->getDouble(); }
  std::string getName() const                            { return _fileReader->getName(); }
  std::string lineNb() const;
  

  std::set<int> _encounteredPiles;

  SauvUtilities::FileReader*      _fileReader;
  SauvUtilities::IntermediateMED* _iMed;
};
}
#endif
