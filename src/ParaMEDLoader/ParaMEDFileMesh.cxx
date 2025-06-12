// Copyright (C) 2007-2025  CEA, EDF
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
// Author : Anthony Geay (EDF R&D)

#include "ParaMEDFileMesh.hxx"
#include "MCAuto.hxx"
#include "MEDFileMesh.hxx"
#include "MEDFileMeshLL.hxx"
#include "MEDLoader.hxx"
#include "MEDFileField1TS.hxx"
#include "MEDFileUtilities.hxx"
#include "MEDFileEntities.hxx"
#include <iostream>
#include <fstream>

// From MEDLOader.cxx TU
extern med_geometry_type typmai3[INTERP_KERNEL::NORM_MAXTYPE];

using namespace MEDCoupling;

MEDFileMesh *
ParaMEDFileMesh::New(
    int iPart,
    int nbOfParts,
    const std::string &fileName,
    const std::string &mName,
    int dt,
    int it,
    MEDFileMeshReadSelector *mrs
)
{
    MEDFileUtilities::CheckFileForRead(fileName);
    MEDCoupling::MEDCouplingMeshType meshType;
    MEDFileUtilities::AutoFid fid(MEDfileOpen(fileName.c_str(), MED_ACC_RDONLY));
    int dummy0, dummy1;
    std::string dummy2;
    MEDCoupling::MEDCouplingAxisType dummy3;
    MEDFileMeshL2::GetMeshIdFromName(fid, mName, meshType, dummy3, dummy0, dummy1, dummy2);
    switch (meshType)
    {
        case UNSTRUCTURED:
        {
            return ParaMEDFileUMesh::New(iPart, nbOfParts, fileName, mName, dt, it, mrs);
        }
        default:
            throw INTERP_KERNEL::Exception("ParaMEDFileMesh::New : only unstructured mesh supported for the moment !");
    }
}

MEDFileMesh *
ParaMEDFileMesh::ParaNew(
    int iPart,
    int nbOfParts,
    const MPI_Comm &com,
    const MPI_Info &nfo,
    const std::string &fileName,
    const std::string &mName,
    int dt,
    int it,
    MEDFileMeshReadSelector *mrs
)
{
    MEDFileUtilities::CheckFileForRead(fileName);
    MEDCoupling::MEDCouplingMeshType meshType;
    MEDFileUtilities::AutoFid fid(MEDfileOpen(fileName.c_str(), MED_ACC_RDONLY));
    int dummy0, dummy1;
    std::string dummy2;
    MEDCoupling::MEDCouplingAxisType dummy3;
    MEDFileMeshL2::GetMeshIdFromName(fid, mName, meshType, dummy3, dummy0, dummy1, dummy2);
    switch (meshType)
    {
        case UNSTRUCTURED:
        {
            return ParaMEDFileUMesh::ParaNew(iPart, nbOfParts, com, nfo, fileName, mName, dt, it, mrs);
        }
        default:
            throw INTERP_KERNEL::Exception(
                "ParaMEDFileMesh::ParaNew : only unstructured mesh supported for the moment !"
            );
    }
}

MEDFileUMesh *
ParaMEDFileUMesh::New(
    int iPart,
    int nbOfParts,
    const std::string &fileName,
    const std::string &mName,
    int dt,
    int it,
    MEDFileMeshReadSelector *mrs
)
{
    MEDFileUtilities::CheckFileForRead(fileName);
    MEDFileUtilities::AutoFid fid(MEDfileOpen(fileName.c_str(), MED_ACC_RDONLY));
    return ParaMEDFileUMesh::NewPrivate(fid, iPart, nbOfParts, fileName, mName, dt, it, mrs);
}

/*!
 * Opens the given file in parallel so that each processor can load a specific part of the mesh \a mName.
 * Each processor will load the cells contained in the vector \a distrib (only nodes lying on those cells will be
 * loaded), in order to read the entire mesh in parallel (memory consumption is thus distributed among all the
 * processes).
 * \param [in] distrib - map defining for each geometric type, the corresponding vector of cells we want to load with
 * c-type indexing (starting from zero). Each vector in this map has an independant numerotation, which means on one
 * processor, vectors of different types may contain the same numbers, they will not refer to the same cells (the i-th
 * cell of a type A does not correspond to the i-th cell of type B) However they have to differ from one processor to
 * another, as to ensure that: 1) each processor only loads a unique part of the mesh 2) the combined distribution
 * vectors cover the entire mesh
 * \param [in] com - group of MPI processes that will read the file
 * \param [in] nfo- MPI info object (used to manage MPI routines)
 * \param [in] filename - name of the file we want to read
 * \param [in] mName - name of the mesh we want to read
 * \param [in] dt - order at which to read the mesh
 * \param [in] it - iteration at which to read the mesh
 * \param [in] mrs - object used to read additional low-level information
 * \return MEDFileUMesh* - a new instance of MEDFileUMesh. The
 *         caller is to delete this mesh using decrRef() as it is no more needed.
 * \throw exception if the partition of the mesh cells defined by \a distrib does not cover the whole mesh
 */
MEDFileUMesh *
ParaMEDFileUMesh::ParaNew(
    const std::map<INTERP_KERNEL::NormalizedCellType, std::vector<mcIdType>> &distrib,
    const MPI_Comm &com,
    const MPI_Info &nfo,
    const std::string &fileName,
    const std::string &mName,
    int dt,
    int it,
    MEDFileMeshReadSelector *mrs
)
{
    MEDFileUtilities::CheckFileForRead(fileName);
#ifdef HDF5_IS_PARALLEL
    MEDFileUtilities::AutoFid fid(MEDparFileOpen(fileName.c_str(), MED_ACC_RDONLY, com, nfo));
#else
    MEDFileUtilities::AutoFid fid(MEDfileOpen(fileName.c_str(), MED_ACC_RDONLY));
#endif
    return ParaMEDFileUMesh::NewPrivate(fid, com, distrib, fileName, mName, dt, it, mrs);
}

/*!
 * Opens the given file in parallel so that each processor can load its part of the mesh \a mName.
 * The mesh will be equally and linearly distributed among all processes:
 * the list of cells will be divided into \a nbOfParts slices and only slice \a iPart (cells and nodes lying on those
 * cells) will be loaded by the current processor. The entire mesh is thus read in parallel and memory consumption is
 * divided among the group of processes.
 * \param [in] iPart - part of the mesh that will be loaded
 * \param [in] nbOfParts - total number of parts in which to divide the mesh
 * \param [in] com - group of MPI processes that will read the file
 * \param [in] nfo- MPI info object (used to manage MPI routines)
 * \param [in] filename - name of the file we want to read
 * \param [in] mName - name of the mesh we want to read
 * \param [in] dt - Time order at which to read the mesh
 * \param [in] it - Time iteration at which to read the mesh
 * \param [in] mrs - object used to read additional low-level information
 * \return MEDFileUMesh* - a new instance of MEDFileUMesh. The
 *         caller is to delete this mesh using decrRef() as it is no more needed.
 */
MEDFileUMesh *
ParaMEDFileUMesh::ParaNew(
    int iPart,
    int nbOfParts,
    const MPI_Comm &com,
    const MPI_Info &nfo,
    const std::string &fileName,
    const std::string &mName,
    int dt,
    int it,
    MEDFileMeshReadSelector *mrs
)
{
    MEDFileUtilities::CheckFileForRead(fileName);
#ifdef HDF5_IS_PARALLEL
    MEDFileUtilities::AutoFid fid(
        MEDparFileOpen(fileName.c_str(), MED_ACC_RDONLY, com, nfo)
    );  // MPI_COMM_WORLD, MPI_INFO_NULL
#else
    MEDFileUtilities::AutoFid fid(MEDfileOpen(fileName.c_str(), MED_ACC_RDONLY));
#endif
    return ParaMEDFileUMesh::NewPrivate(fid, iPart, nbOfParts, fileName, mName, dt, it, mrs);
}

/*!
 * Loads mesh \a mName in parallel using a custom partition of the mesh cells among the processes.
 * See ParaMEDFileUMesh::ParaNew for detailed description.
 */
MEDFileUMesh *
ParaMEDFileUMesh::NewPrivate(
    med_idt fid,
    const MPI_Comm &com,
    const std::map<INTERP_KERNEL::NormalizedCellType, std::vector<mcIdType>> &distrib,
    const std::string &fileName,
    const std::string &mName,
    int dt,
    int it,
    MEDFileMeshReadSelector *mrs
)
{
    MCAuto<MEDFileUMesh> ret;
    ret = MEDFileUMesh::LoadPartOfFromUserDistrib(fid, mName, distrib, dt, it, mrs);
    return ret.retn();
}

/*!
 * Loads mesh \a mName in parallel using a slice partition of the mesh cells among the processes
 * See ParaMEDFileUMesh::ParaNew for detailed description.
 */
MEDFileUMesh *
ParaMEDFileUMesh::NewPrivate(
    med_idt fid,
    int iPart,
    int nbOfParts,
    const std::string &fileName,
    const std::string &mName,
    int dt,
    int it,
    MEDFileMeshReadSelector *mrs
)
{
    MCAuto<MEDFileUMesh> ret;
    int meshDim, spaceDim;
    mcIdType numberOfNodes;
    std::vector<std::vector<std::pair<INTERP_KERNEL::NormalizedCellType, int>>> typesDistrib(
        GetUMeshGlobalInfo(fileName, mName, meshDim, spaceDim, numberOfNodes)
    );
    std::vector<INTERP_KERNEL::NormalizedCellType> types;
    std::vector<mcIdType> distrib;
    for (std::vector<std::vector<std::pair<INTERP_KERNEL::NormalizedCellType, int>>>::const_iterator it0 =
             typesDistrib.begin();
         it0 != typesDistrib.end();
         it0++)
        for (std::vector<std::pair<INTERP_KERNEL::NormalizedCellType, int>>::const_iterator it1 = (*it0).begin();
             it1 != (*it0).end();
             it1++)
        {
            types.push_back((*it1).first);
            mcIdType tmp[3];
            DataArray::GetSlice(0, (*it1).second, 1, iPart, nbOfParts, tmp[0], tmp[1]);
            tmp[2] = 1;
            distrib.insert(distrib.end(), tmp, tmp + 3);
        }
    ret = MEDFileUMesh::LoadPartOf(fid, mName, types, distrib, dt, it, mrs);
    return ret.retn();
}

/*!
 * Loads field \a fName laying on mesh \a mName from the filename \a fileName in parallel:
 * each processor will load their portion of the field (ie the portion laying on the cells/nodes in the vector \a
 * distrib given in the parameters). WARNING : this will only load the array of values of the field, additionnal
 * information of the local field such as the number of its tuples might be incorrect WARNING : this method does not
 * check the distribution given as input !
 * \param [in] com - group of MPI processes that will read the file
 * \param [in] nfo- MPI info object (used to manage MPI routines)
 * \param [in] fileName - name of the file containing the field
 * \param [in] fName - name of the field we want to load
 * \param [in] mName - name of the mesh on which the field is defined
 * \param [in] distrib - vector of cells/nodes on which we want to load the field \a fName (with c-type indexing, so
 * starting from zero).
 * \param [in] loc - localisation of the field
 * \param [in] geoType - geometrical type of the cells on which the field is laying (only needed if the field is located
 * on cells)
 * \param [in] dt - Time order at which to read the field
 * \param [in] it - Time iteration at which to read the field
 * \return MEDFileField1TS* - a new instance of MEDFileField1TS. The
 *         caller is to delete it using decrRef() as it is no more needed.
 * \throw exception if the field is not of type FLOAT64
 */
MEDFileField1TS *
ParaMEDFileField1TS::ParaNew(
    const MPI_Comm &com,
    const MPI_Info &nfo,
    const std::string &fileName,
    const std::string &fName,
    const std::string &mName,
    const std::vector<mcIdType> &distrib,
    TypeOfField loc,
    INTERP_KERNEL::NormalizedCellType geoType,
    int dt,
    int it
)
{
    MEDFileUtilities::CheckFileForRead(fileName);
#ifdef HDF5_IS_PARALLEL
    MEDFileUtilities::AutoFid fid(MEDparFileOpen(fileName.c_str(), MED_ACC_RDONLY, com, nfo));
#else
    MEDFileUtilities::AutoFid fid(MEDfileOpen(fileName.c_str(), MED_ACC_RDONLY));
#endif
    return ParaMEDFileField1TS::NewPrivate(fid, com, fName, mName, distrib, loc, geoType, dt, it);
}

/*!
 * Loads field \a fName in parallel using a custom partition of the mesh entities (cells or nodes) on which the field is
 * defined among the processes. See ParaMEDFileField1TS::ParaNew for detailed description.
 */
MEDFileField1TS *
ParaMEDFileField1TS::NewPrivate(
    med_idt fid,
    const MPI_Comm &com,
    const std::string &fName,
    const std::string &mName,
    const std::vector<mcIdType> &distrib,
    TypeOfField loc,
    INTERP_KERNEL::NormalizedCellType geoType,
    int dt,
    int it
)
{
    std::vector<std::pair<TypeOfField, INTERP_KERNEL::NormalizedCellType>> tmp = {{loc, geoType}};
    INTERP_KERNEL::AutoCppPtr<MEDFileEntities> entities(MEDFileEntities::BuildFrom(&tmp));
    MCAuto<MEDFileAnyTypeField1TS> partFile(MEDFileAnyTypeField1TS::NewAdv(fid, fName, dt, it, entities, distrib));

    MCAuto<MEDFileField1TS> ret(MEDCoupling::DynamicCast<MEDFileAnyTypeField1TS, MEDFileField1TS>(partFile));
    if (ret.isNotNull())
        return ret.retn();
    else
        throw INTERP_KERNEL::Exception("ParaMEDFileField1TS::ParaNew : only FLOAT64 field supported for the moment !");
}

MEDFileMeshes *
ParaMEDFileMeshes::New(int iPart, int nbOfParts, const std::string &fileName)
{
    std::vector<std::string> ms(GetMeshNames(fileName));
    MCAuto<MEDFileMeshes> ret(MEDFileMeshes::New());
    for (std::vector<std::string>::const_iterator it = ms.begin(); it != ms.end(); it++)
    {
        MCAuto<MEDFileMesh> mesh(ParaMEDFileMesh::New(iPart, nbOfParts, fileName, (*it)));
        ret->pushMesh(mesh);
    }
    return ret.retn();
}

MEDFileMeshes *
ParaMEDFileMeshes::ParaNew(
    int iPart, int nbOfParts, const MPI_Comm &com, const MPI_Info &nfo, const std::string &fileName
)
{
    std::vector<std::string> ms(GetMeshNames(fileName));
    MCAuto<MEDFileMeshes> ret(MEDFileMeshes::New());
    for (std::vector<std::string>::const_iterator it = ms.begin(); it != ms.end(); it++)
    {
        MCAuto<MEDFileMesh> mesh(ParaMEDFileMesh::ParaNew(iPart, nbOfParts, com, nfo, fileName, (*it)));
        ret->pushMesh(mesh);
    }
    return ret.retn();
}
