// Copyright (C) 2021-2025  CEA, EDF
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

#ifndef MESHFORMATWRITER_HXX
#define MESHFORMATWRITER_HXX

#include <vector>
#include <string>
#include <string>
#include "MCAuto.hxx"
#include "InterpKernelException.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "MCType.hxx"
#include "MEDMESHConverterUtilities.hxx"
#include "libmesh5.hxx"

#include "MEDFileMesh.hxx"
#include <fstream>

namespace MEDCoupling
{
class DataArrayDouble;
class DataArrayInt;
class MEDFileData;
class MEDFileFields;
class MEDFileFieldMultiTS;
class MEDFileUMesh;

struct MeshFormatNode
{
    MeshFormatNode(double x, double y, double z = 0, int id = 0, int famId = 0) : _id(id), _famId(famId)
    {
        xyz[0] = x;
        xyz[1] = y;
        xyz[2] = z;
    }
    int _id;
    int _famId;
    double xyz[3];

    //~inline friend bool operator==(const MeshFormatNode& a, const MeshFormatNode&b ) { return ((a._type == b._type) &&
    //(a._id == b._id));}
};

struct MeshFormatCell
{
    MeshFormatCell(INTERP_KERNEL::NormalizedCellType type, int id = 0, int famId = 0)
        : _type(type), _id(id), _famId(famId)
    {
    }
    INTERP_KERNEL::NormalizedCellType _type;
    int _id;
    int _famId;
    void setConn(const std::vector<MEDCoupling::mcIdType> cpy)
    {
        conn.clear();
        std::copy(cpy.begin(), cpy.end(), std::back_inserter(conn));
    }
    std::vector<MEDCoupling::mcIdType> conn;
};

class MeshFormatWriter
{
   public:
    MEDLOADER_EXPORT MeshFormatWriter();
    MEDLOADER_EXPORT MeshFormatWriter(const std::string &meshFileName, const std::vector<std::string> &fieldFileNames);
    MEDLOADER_EXPORT ~MeshFormatWriter();
    MEDLOADER_EXPORT void setMeshFileName(const std::string &meshFileName);
    MEDLOADER_EXPORT std::string getMeshFileName() const;
    MEDLOADER_EXPORT void setFieldFileNames(const std::vector<std::string> &fieldFileNames);
    MEDLOADER_EXPORT std::vector<std::string> getFieldFileNames() const;
    MEDLOADER_EXPORT void setMEDFileDS(MEDCoupling::MEDFileData *mfd);

    MEDLOADER_EXPORT void write();

   private:
    void getNodes(MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh0);
    void getNSEG2(MEDCoupling::mcIdType nbEdgesNSEG2, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh1);
    void getNSEG3(MEDCoupling::mcIdType nbEdgesNSEG2, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh1);
    void getTRI3(MEDCoupling::mcIdType nbTRI3, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh2);
    void getTRI6(MEDCoupling::mcIdType nbTRI6, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh2);
    void getQUAD4(MEDCoupling::mcIdType nbQUAD4, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh2);
    void getQUAD8(MEDCoupling::mcIdType nbQUAD8, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh2);
    void getQUAD9(MEDCoupling::mcIdType nbQUAD9, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh2);
    void getTETRA4(MEDCoupling::mcIdType nbTETRA4, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh3);
    void getTETRA10(MEDCoupling::mcIdType nbTETRA10, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh3);
    void getPYRA5(MEDCoupling::mcIdType nbPYRA5, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh3);
    void getHEXA8(MEDCoupling::mcIdType nbHEXA8, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh3);
    void getHEXA20(MEDCoupling::mcIdType nbHEXA20, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh3);
    void getHEXA27(MEDCoupling::mcIdType nbHEXA27, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh3);
    void getPENTA6(MEDCoupling::mcIdType nbPENTA6, MEDCoupling::MCAuto<MEDCoupling::MEDCouplingUMesh> umesh3);
    int getGmfSolKwd(const int nbComp, const int dim);
    MeshFormat::Status setFieldOnNodes(
        MEDCoupling::MEDFileFieldMultiTS *f, int iteration, int order, size_t compInfoSize
    );
    MeshFormat::Status setFieldOnCells(
        MEDCoupling::MEDFileFieldMultiTS *f, int iteration, int order, std::vector<int> levs
    );
    void extractSymetricTensor(double fullTensor[], double *&symTensor);
    void linkFamilyToNodes();
    void linkFamilyToCells();
    void writeCells();

    bool checkFileName();
    bool checkFieldFileName();
    void forward_shift(std::vector<MEDCoupling::mcIdType> &conn);
    MeshFormat::Status perform();
    MeshFormat::Status performFields();
    MeshFormat::Status addMessage(const std::string &msg, const bool isFatal = false);

    std::string _meshFileName;
    MeshFormat::MeshFormatParser _writer;
    std::vector<std::string> _fieldFileNames;
    MEDCoupling::MEDFileData *_mfd;
    MEDCoupling::MCAuto<MEDCoupling::MEDFileMesh> _mesh;
    std::vector<MEDCoupling::MCAuto<MEDCoupling::MEDFileFieldMultiTS> > _fields;
    std::vector<std::string> _myErrorMessages;
    MeshFormat::Status _myStatus;
    int _myCurrentFileId, _dim, _version;
    std::string _myCurrentOpenFile;
    std::map<int, MeshFormatNode> _idNodeToNode;
    std::map<INTERP_KERNEL::NormalizedCellType, std::map<int, MeshFormatCell> > _typeToIdCellToCell;
};
}  // namespace MEDCoupling
#endif  // MESHFORMATWRITER_HXX
