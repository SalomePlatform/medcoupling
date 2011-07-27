//  Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#ifndef __MEDFILEFIELD_HXX__
#define __MEDFILEFIELD_HXX__

#include "MEDCouplingAutoRefCountObjectPtr.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingMemArray.hxx"

#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <vector>
#include <string>

extern "C"
{
#include "med.h"
}

namespace ParaMEDMEM
{
  class MEDFileFieldLoc : public RefCountObject
  {
  public:
    const std::string& getName() const { return _name; } 
  private:
    int _dim;
    int _nb_gauss_pt;
    int _nb_node_per_cell;
    std::string _name;
    INTERP_KERNEL::NormalizedCellType _geo_type;
    std::vector<double> _ref_coo;
    std::vector<double> _gs_coo;
    std::vector<double> _w;
  };

  class MEDFileFieldPerMeshPerType;
  class MEDFileField1TSWithoutDAS;
  class MEDFileFieldPerMesh;

  class MEDFileFieldPerMeshPerTypePerDisc : public RefCountObject
  {
  public:
    static MEDFileFieldPerMeshPerTypePerDisc *New(MEDFileFieldPerMeshPerType *fath, med_idt fid) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerType *getFather() const;
    int getIteration() const;
    int getOrder() const;
    std::string getName() const;
    std::string getMeshName() const;
    TypeOfField getType() const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    const std::vector<std::string>& getInfos() const;
    std::string getProfile() const;
    std::string getLocalization() const;
  private:
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, med_idt fid) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldPerMeshPerType *_father;
    MEDCouplingAutoRefCountObjectPtr< DataArrayDouble > _arr;
    std::string _profile;
    std::string _localization;
    std::string _mesh_name;
  };

  class MEDFileFieldPerMeshPerType : public RefCountObject
  {
  public:
    static MEDFileFieldPerMeshPerType *New(MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMesh *getFather() const;
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    int getIteration() const;
    int getOrder() const;
    std::string getName() const;
    std::string getMeshName() const;
    TypeOfField getType() const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    const std::vector<std::string>& getInfos() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
  private:
    MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldPerMesh *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> > _field_pm_pt_pd;
    TypeOfField _type;
    INTERP_KERNEL::NormalizedCellType _geo_type;
  };

  class MEDFileFieldPerMesh : public RefCountObject
  {
  public:
    static MEDFileFieldPerMesh *New(MEDFileField1TSWithoutDAS *fath, const char *meshName, double time);
    void pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType);
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    double getTime() const;
    int getIteration() const;
    int getOrder() const;
    std::string getName() const;
    std::string getMeshName() const;
    int getNumberOfComponents() const;
    const std::vector<std::string>& getInfos() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
  private:
    MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, const char *meshName, double time);
  private:
    MEDFileField1TSWithoutDAS *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > > _field_pm_pt;
    std::string _mesh_name;
    double _time;
  private:
    mutable std::vector<TypeOfField> _types;
    mutable std::vector<INTERP_KERNEL::NormalizedCellType> _geo_types;
  };

  class MEDFieldFieldGlobs
  {
  public:
    void loadProfileInFile(med_idt fid, int id, const char *pflName, int lgth) throw(INTERP_KERNEL::Exception);
    void loadProfileInFile(med_idt fid, int id);
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _pfls;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> > _locs;
  };

  class MEDFileField1TSWithoutDAS : public RefCountObject
  {
  public:
    static MEDFileField1TSWithoutDAS *New(const char *fieldName, int iteration, int order, const std::vector<std::string>& infos);
    void pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, double time, const char *meshName) const;
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    int getIteration() const { return _iteration; }
    int getOrder() const { return _order; }
    std::string getName() const { return _name; }
    int getNumberOfComponents() const { return _infos.size(); }
    const std::vector<std::string>& getInfos() const { return _infos; }
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
  protected:
    MEDFileField1TSWithoutDAS(const char *fieldName, int iteration, int order, const std::vector<std::string>& infos);
  protected:
    std::string _name;
    std::vector<std::string> _infos;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > > _field_per_mesh;
    int _iteration;
    int _order;
  private:
    mutable std::vector<TypeOfField> _types;
    mutable std::vector<INTERP_KERNEL::NormalizedCellType> _geo_types;
    mutable std::vector<double> _times;
    mutable std::vector<std::string> _meshes;
  };

  /*!
   * User class.
   */
  class MEDFileField1TS : public MEDFileField1TSWithoutDAS, public MEDFieldFieldGlobs
  {
  public:
    static MEDFileField1TS *New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  };
  
  class MEDFileFieldMultiTSWithoutDAS : public RefCountObject
  {
  public:
    static MEDFileFieldMultiTSWithoutDAS *New(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileFieldMultiTSWithoutDAS(const char *fieldName);
    MEDFileFieldMultiTSWithoutDAS(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    void appendTimeStepEntry(med_idt fid, med_entite_maillage entity, int i, int j) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
  protected:
    std::string _name;
    std::vector<std::string> _infos;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileField1TSWithoutDAS>  > _time_steps;
  };

  /*!
   * User class.
   */
  class MEDFileFieldMultiTS : public MEDFileFieldMultiTSWithoutDAS, public MEDFieldFieldGlobs
  {
  public:
    static MEDFileFieldMultiTS *New(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  };

  /*!
   * Use class.
   */
  class MEDFileFields : public RefCountObject, public MEDFieldFieldGlobs
  {
  public:
    static MEDFileFields *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const;
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
  private:
    MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutDAS> > _fields;
  };
}

#endif
