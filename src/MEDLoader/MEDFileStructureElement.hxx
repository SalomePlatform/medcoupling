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

#ifndef __MEDFILESTRUCTUREELEMENT_HXX__
#define __MEDFILESTRUCTUREELEMENT_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileUtilities.txx"
#include "MEDFileMesh.hxx"

#include "MEDCouplingRefCountObject.hxx"

namespace MEDCoupling
{
  class MEDFileStructureElement;
  class MEDFileMeshSupports;
  class MEDFileUMesh;
  
  class MEDFileSEHolder
  {
  public:
    std::string getModelName() const;
    std::string getName() const;
  protected:
    MEDFileSEHolder(MEDFileStructureElement *father):_father(father) { }
    void setName(const std::string& name);
    std::size_t getHeapMemorySizeLoc() const;
  private:
    MEDFileStructureElement *_father;
    std::string _name;
  };
  
class MEDFileSEConstAtt : public RefCountObject, public MEDFileWritableStandAlone, public MEDFileSEHolder
  {
  public:
    static MEDFileSEConstAtt *New(med_idt fid, MEDFileStructureElement *father, int idCstAtt, const MEDFileUMesh *mesh);
  public:
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    void writeLL(med_idt fid) const;
    void setProfile(const std::string& name);
    std::string getProfile() const;
  private:
    MEDFileSEConstAtt(med_idt fid, MEDFileStructureElement *father, int idCstAtt, const MEDFileUMesh *mesh);
  private:
    std::string _pfl;
    TypeOfField _tof;
    MCAuto<DataArray> _val;
  };
  
  class MEDFileSEVarAtt : public RefCountObject, public MEDFileWritableStandAlone, public MEDFileSEHolder
  {
  public:
    static MEDFileSEVarAtt *New(med_idt fid, MEDFileStructureElement *father, int idVarAtt);
  public:
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    void writeLL(med_idt fid) const;
    int getNbOfComponents() const { return _nb_compo; }
    MCAuto<DataArray> getGenerator() const { return _gen; }
  private:
    MEDFileSEVarAtt(med_idt fid, MEDFileStructureElement *father, int idVarAtt);
  private:
    int _nb_compo;
    MCAuto<DataArray> _gen;
  };
  
  class MEDFileStructureElement : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileStructureElement *New(med_idt fid, int idSE, const MEDFileMeshSupports *ms);
    MEDLOADER_EXPORT std::string getName() const;
    MEDLOADER_EXPORT int getDynGT() const;
    MEDLOADER_EXPORT TypeOfField getEntity() const;
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT std::vector<std::string> getVarAtts() const;
    MEDLOADER_EXPORT const MEDFileSEVarAtt *getVarAtt(const std::string& varName) const;
  public:
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    void writeLL(med_idt fid) const;
  public:
    static MCAuto<DataArray> BuildFrom(med_attribute_type mat);
    static int EffectiveNbCompo(med_attribute_type mat, int nbCompo);
  private:
    MEDFileStructureElement(med_idt fid, int idSE, const MEDFileMeshSupports *ms);
  private:
    int _id_type;
    std::string _name;
    std::string _sup_mesh_name;
    INTERP_KERNEL::NormalizedCellType _geo_type;
    TypeOfField _tof;
    int _dim;
    std::vector< MCAuto<MEDFileSEConstAtt> > _cst_att;
    std::vector< MCAuto<MEDFileSEVarAtt> > _var_att;
  };
  
  class MEDFileStructureElements : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileStructureElements *New(const std::string& fileName, const MEDFileMeshSupports *ms);
    MEDLOADER_EXPORT static MEDFileStructureElements *New(med_idt fid, const MEDFileMeshSupports *ms);
    MEDLOADER_EXPORT static MEDFileStructureElements *New();
    MEDLOADER_EXPORT int getNumberOf() const;
    MEDLOADER_EXPORT std::vector<int> getDynGTAvail() const;
    MEDLOADER_EXPORT const MEDFileStructureElement *getWithGT(int idGT) const;
    MEDLOADER_EXPORT int getNumberOfNodesPerSE(const std::string& seName) const;
    MEDLOADER_EXPORT const MEDFileStructureElement *getSEWithName(const std::string& seName) const;
    MEDLOADER_EXPORT std::vector<std::string> getVarAttsOf(const std::string& seName) const;
    MEDLOADER_EXPORT const MEDFileSEVarAtt *getVarAttOf(const std::string &seName, const std::string& varName) const;
    MEDLOADER_EXPORT const MEDFileUMesh *getSupMeshWithName(const std::string& name) const;
  public:
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    std::size_t getHeapMemorySizeWithoutChildren() const;
    void writeLL(med_idt fid) const;
  private:
    MEDFileStructureElements(med_idt fid, const MEDFileMeshSupports *ms);
    MEDFileStructureElements();
    ~MEDFileStructureElements();
  private:
    std::vector< MCAuto<MEDFileStructureElement> > _elems;
    MCConstAuto<MEDFileMeshSupports> _sup;
  };
}

#endif
