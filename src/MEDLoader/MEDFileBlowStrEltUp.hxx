// Copyright (C) 2007-2017  CEA/DEN, EDF R&D
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

#ifndef __MEDFILEBLOWSTRELTUP_HXX__
#define __MEDFILEBLOWSTRELTUP_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileUtilities.txx"
#include "MEDFileMesh.hxx"
#include "MEDFileField.hxx"
#include "MEDFileStructureElement.hxx"

#include "MEDCouplingRefCountObject.hxx"

namespace MEDCoupling
{ 
  class MEDFileBlowStrEltUp
  {
  public:
    MEDFileBlowStrEltUp(const MEDFileFields *fsOnlyOnSE, const MEDFileMeshes *ms, const MEDFileStructureElements *ses);
    void generate(MEDFileMeshes *msOut, MEDFileFields *allZeOutFields);
    static void DealWithSE(MEDFileFields *fs, MEDFileMeshes *ms, const MEDFileStructureElements *ses);
  private:
    MCAuto<MEDFileFields> splitFieldsPerLoc(const MEDFileFields *fields, const MEDFileUMesh *mesh, MEDFileMeshes *msOut, MEDFileFields *allZeOutFields);
    MCAuto<MEDFileEltStruct4Mesh> dealWithSEInMesh(const std::string& seName, MEDFileUMesh *mesh, MCAuto<MEDFileUMesh>& mOut, MCAuto<MEDFileFields>& fsOut) const;
    MCAuto<MEDFileEltStruct4Mesh> dealWithMEDBALLInMesh(const MEDFileUMesh *mesh, MCAuto<MEDFileUMesh>& mOut, MCAuto<MEDFileFields>& fsOut) const;
    void dealWithSEInFields(const std::string& seName, const MEDFileFields *fs, const MEDFileEltStruct4Mesh *zeStr, const MEDFileFields *varAtt, MEDFileFields *zeOutputs) const;
    void dealWithMEDBALLSInFields(const MEDFileFields *fs, const MEDFileEltStruct4Mesh *zeStr, const MEDFileFields *varAtt, MEDFileFields *zeOutputs) const;
    static std::string BuildNewMeshName(const std::string& meshName, const std::string& seName);
    static std::string BuildVarAttName(std::size_t iPart, std::size_t totINbParts, std::size_t jPart, std::size_t totJNbParts, const std::string& name);
    static void DealWithConflictNames(MEDFileAnyTypeFieldMultiTS *fmtsToAdd, const MEDFileFields *fs);
  public:
    static const char MED_BALL_STR[];
  private:
    std::vector< MCAuto<MEDFileFields> > _elts;// split by pair(meshName,SEName)
    MCConstAuto<MEDFileMeshes> _ms;
    MCConstAuto<MEDFileStructureElements> _ses;
  };
}

#endif
