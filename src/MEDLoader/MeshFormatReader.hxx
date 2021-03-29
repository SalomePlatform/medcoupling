// Copyright (C) 2021  CEA/DEN, EDF R&D
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

#ifndef MESHFORMATREADER_HXX
#define MESHFORMATREADER_HXX

#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <utility>
#include <iterator>
#include "MCAuto.hxx"
#include "InterpKernelException.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "MCType.hxx"
#include "MEDMESHConverterUtilities.hxx"
#include "libmesh5.hxx"

#include <fstream>

#ifndef WIN32
#include <features.h>
#endif

namespace MEDCoupling
{
  class DataArrayDouble;
  class DataArrayInt;
  class MEDFileData;
  class MEDFileFields;
  class MEDFileFieldMultiTS;
  class MEDFileUMesh;
  class MEDCouplingUMesh;

  struct MeshFormatElement
  {
    MeshFormatElement(int type, int id = 0) : _type(type), _id(id) {}
    int _type;
    int _id;

    inline friend bool operator==(const MeshFormatElement &a, const MeshFormatElement &b)
    {
      return ((a._type == b._type) && (a._id == b._id));
    }
  };

  struct MeshFormatFamily
  {

    std::map<int, std::vector<MeshFormatElement> *> _meshFormatFams;
    std::map<int, std::vector<MeshFormatElement> *> _meshFormatFams_0;
    std::map<int, std::vector<MeshFormatElement> *> _meshFormatFams_1;
    std::map<int, std::vector<MeshFormatElement> *> _meshFormatFams_2;
    std::map<int, std::vector<MeshFormatElement> *> _meshFormatFamsNodes;

  public:
    void insert(std::pair<int, MeshFormatElement> addToFamily, int dimRelMax)
    {
      insertPairInMap(_meshFormatFams, addToFamily);
      insertPairInMap(getMapAtLevel(dimRelMax), addToFamily);
    }

  private:
    void insertPairInMap(std::map<int, std::vector<MeshFormatElement> *> &aMap, std::pair<int, MeshFormatElement> addToFamily)
    {
      std::map<int, std::vector<MeshFormatElement> *>::iterator it = aMap.find(addToFamily.first);
      if (it != aMap.end())
      {
        aMap[addToFamily.first]->push_back(addToFamily.second);
      }
      else
      {
        std::vector<MeshFormatElement> *tmpVec = new std::vector<MeshFormatElement>;
        tmpVec->push_back(addToFamily.second);
        aMap.insert(std::pair<int, std::vector<MeshFormatElement> *>(addToFamily.first, tmpVec));
      }
    }

  public:
    void remove(std::pair<int, MeshFormatElement> removeFromFamily, int dimRelMax)
    {
      removePairFromMap(_meshFormatFams, removeFromFamily);
      removePairFromMap(getMapAtLevel(dimRelMax), removeFromFamily);
    }

  private:
    void removePairFromMap(std::map<int, std::vector<MeshFormatElement> *> &aMap, const std::pair<int, MeshFormatElement> removeFromFamily)
    {

      if (!aMap.size())
        return;
      std::map<int, std::vector<MeshFormatElement> *>::iterator itTmp = aMap.find(removeFromFamily.first);
      if (itTmp == aMap.end())
        return;
      else
      {
        std::vector<MeshFormatElement> *tmpVec2 = aMap[removeFromFamily.first];
        std::vector<MeshFormatElement>::iterator itt2;
        const MeshFormatElement e = removeFromFamily.second;
        itt2 = std::find(tmpVec2->begin(), tmpVec2->end(), e);
        if (itt2 != tmpVec2->end())
          tmpVec2->erase(itt2);

        if (!tmpVec2->size())
        {
          delete tmpVec2;
          aMap.erase(removeFromFamily.first);
        }
      }
    }

  public:
    std::map<int, std::vector<MeshFormatElement> *> &getMapAtLevel(int dimRelMax)
    {
      switch (dimRelMax)
      {
      case 0:
        return _meshFormatFams_0;
        break;
      case -1:
        return _meshFormatFams_1;
        break;
      case -2:
        return _meshFormatFams_2;
        break;
      case 1:
        return _meshFormatFamsNodes;
        break;
      }
    }

  public:
    ~MeshFormatFamily()
    {

      freeMap(_meshFormatFams);
      freeMap(_meshFormatFams_0);
      freeMap(_meshFormatFams_1);
      freeMap(_meshFormatFams_2);
      freeMap(_meshFormatFamsNodes);
    }

  private:
    void freeMap(std::map<int, std::vector<MeshFormatElement> *> &aMap)
    {
      std::map<int, std::vector<MeshFormatElement> *>::iterator it = aMap.begin();
      for (; it != aMap.end(); ++it)
        delete it->second;
    }
  };
  class MeshFormatReader
  {
  public:
    MEDLOADER_EXPORT MeshFormatReader();
    MEDLOADER_EXPORT MeshFormatReader(const std::string &meshFileName, const std::vector<std::string> &fieldFileName);
    MEDLOADER_EXPORT ~MeshFormatReader();

    MEDLOADER_EXPORT MEDCoupling::MCAuto<MEDCoupling::MEDFileData> loadInMedFileDS();
    MEDLOADER_EXPORT void setMeshName(const std::string &theMeshName);
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT void setFile(const std::string &theFileName);
    MEDLOADER_EXPORT void setFieldFileNames(const std::vector<std::string> &theFieldFileNames);
    MEDLOADER_EXPORT std::vector<std::string> getFieldFileNames() const;
    MEDLOADER_EXPORT std::vector<std::string> getErroMessage() const;

  private:
    MeshFormat::Status addMessage(const std::string &msg, const bool isFatal = false);
    MeshFormat::Status perform();
    MeshFormat::Status performFields();
    MeshFormat::Status setNodes(MEDCoupling::DataArrayDouble *coordArray);
    void setEdges(MEDCoupling::MEDCouplingUMesh *dimMesh1, int nbEdges);
    void setTriangles(MEDCoupling::MEDCouplingUMesh *dimMesh2, int nbTria);
    void setQuadrangles(MEDCoupling::MEDCouplingUMesh *dimMesh2, int nbQuad);
    void setTetrahedras(MEDCoupling::MEDCouplingUMesh *dimMesh3, int nbTet);
    void setPyramids(MEDCoupling::MEDCouplingUMesh *dimMesh3, int nbPyr);
    void setHexahedras(MEDCoupling::MEDCouplingUMesh *dimMesh3, int nbHex);
    void setPrisms(MEDCoupling::MEDCouplingUMesh *dimMesh3, int nbPrism);
    void callParserGetLin(MeshFormat::GmfKwdCod kwd, double *val, int valSize, int *ref);
    void setTypeOfFieldAndDimRel(MeshFormat::GmfKwdCod kwd, MEDCoupling::TypeOfField *typeOfField, int *dimRel);
    void backward_shift(MEDCoupling::mcIdType *, int size);
    void setFields(MeshFormat::GmfKwdCod kwd, int nmbSol, int nbComp);
    INTERP_KERNEL::NormalizedCellType toMedType(MeshFormat::GmfKwdCod kwd);
    void buildFamilies();
    void buildNodesFamilies();
    void buildCellsFamilies();

    std::string _myFile;
    MeshFormat::MeshFormatParser _reader;
    std::string _myCurrentOpenFile;
    int _myCurrentFileId;
    std::string _myMeshName;
    std::vector<std::string> _myFieldFileNames;
    int _dim, _version;
    std::vector<std::string> _myErrorMessages;
    MeshFormat::Status _myStatus;
    MEDCoupling::MCAuto<MEDCoupling::MEDFileData> _myMed;
    MEDCoupling::MCAuto<MEDCoupling::MEDFileUMesh> _uMesh;
    MEDCoupling::MCAuto<MEDCoupling::MEDFileFields> _fields;
    // map id family --element
    MeshFormatFamily _fams;

    int _dim1NbEl;
    int _dim2NbEl;
    int _dim3NbEl;
  };
} // namespace MEDCoupling
#endif //MESHFORMATREADER_HXX
