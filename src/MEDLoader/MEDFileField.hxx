//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
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

#include "MEDFileUtilities.hxx"

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
  class MEDFieldFieldGlobs;
  class MEDCouplingMesh;
  class MEDCouplingFieldDouble;

  class MEDFileFieldLoc : public RefCountObject
  {
  public:
    const std::string& getName() const { return _name; }
    static MEDFileFieldLoc *New(med_idt fid, const char *locName);
    void writeLL(med_idt fid) const;
  private:
    MEDFileFieldLoc(med_idt fid, const char *locName);
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

  class MEDFileFieldPerMeshPerTypePerDisc : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerTypePerDisc *New(MEDFileFieldPerMeshPerType *fath, med_idt fid, int profileIt) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerType *getFather() const;
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getName() const;
    std::string getMeshName() const;
    TypeOfField getType() const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    int getNumberOfTuples() const;
    const std::vector<std::string>& getInfo() const;
    std::string getProfile() const;
    std::string getLocalization() const;
  private:
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, med_idt fid, int profileIt) throw(INTERP_KERNEL::Exception);
    
  private:
    MEDFileFieldPerMeshPerType *_father;
    MEDCouplingAutoRefCountObjectPtr< DataArrayDouble > _arr;
    int _nval;
    int _profile_it;
    std::string _profile;
    std::string _localization;
  };

  class MEDFileFieldPerMeshPerType : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMeshPerType *New(MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMesh *getFather() const;
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getName() const;
    std::string getMeshName() const;
    TypeOfField getType() const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    bool isOnNode(int& type, int& number, const DataArrayInt* &arrs) const throw(INTERP_KERNEL::Exception);
    bool isOnCell(int dimDimReq, int& type, int& number, const DataArrayInt* &arrs) const throw(INTERP_KERNEL::Exception);
    static med_entity_type ConvertIntoMEDFileType(TypeOfField ikType, INTERP_KERNEL::NormalizedCellType ikGeoType, med_geometry_type& medfGeoType);
  private:
    MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldPerMesh *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> > _field_pm_pt_pd;
    TypeOfField _type;
    INTERP_KERNEL::NormalizedCellType _geo_type;
  };

  class MEDFileFieldPerMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMesh *New(MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder);
    void pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType);
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    double getTime() const;
    int getIteration() const;
    int getOrder() const;
    std::string getName() const;
    std::string getMeshName() const { return _mesh_name; }
    int getNumberOfComponents() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    std::vector<int> getDistributionOfTypes(int meshDimRelToMaxExt, int mdim, std::vector<const DataArrayInt *>& arrs) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(int meshDimRelToMaxExt, const MEDFieldFieldGlobs *glob, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder);
  private:
    std::string _mesh_name;
    int _mesh_iteration;
    int _mesh_order;
    int _mesh_csit;
    MEDFileField1TSWithoutDAS *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > > _field_pm_pt;
  private:
    mutable std::vector<TypeOfField> _types;
    mutable std::vector<INTERP_KERNEL::NormalizedCellType> _geo_types;
  };

  class MEDFieldFieldGlobs
  {
  public:
    MEDFieldFieldGlobs(const char *fname);
    void loadProfileInFile(med_idt fid, int id, const char *pflName) throw(INTERP_KERNEL::Exception);
    void loadProfileInFile(med_idt fid, int id);
    void loadGlobals(med_idt fid) throw(INTERP_KERNEL::Exception);
    void writeGlobals(med_idt fid, const MEDFileWritable& opt) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPfls() const;
    std::vector<std::string> getLocs() const;
    void setFileName(const char *fileName);
    const char *getFileName() const { return _file_name.c_str(); }
    virtual std::vector<std::string> getPflsReallyUsed() const = 0;
    virtual std::vector<std::string> getLocsReallyUsed() const = 0;
  protected:
    std::vector< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > _pfls;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldLoc> > _locs;
    std::string _file_name;
  };

  class MEDFileField1TSWithoutDAS : public RefCountObject, public MEDFileWritable
  {
  public:
    int getIteration() const { return _iteration; }
    int getOrder() const { return _order; }
    double getTime() const { return _dt; }
    std::string getName() const { return _name; }
    std::string getDtUnit() const { return _dt_unit; }
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const { return _infos.size(); }
    const std::vector<std::string>& getInfo() const { return _infos; }
    //
    static MEDFileField1TSWithoutDAS *New(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
    void pushBack(TypeOfField type, INTERP_KERNEL::NormalizedCellType geoType, double time) const;
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed2() const;
    std::vector<std::string> getLocsReallyUsed2() const;
  protected:
    MEDCouplingFieldDouble *getFieldAtLevel(int meshDimRelToMaxExt, const MEDFieldFieldGlobs *glob) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(int meshDimRelToMaxExt, const MEDFieldFieldGlobs *glob, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileField1TSWithoutDAS(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
  protected:
    std::string _name;
    std::string _dt_unit;
    std::vector<std::string> _infos;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMesh > > _field_per_mesh;
    int _csit;
    int _iteration;
    int _order;
    double _dt;
  private:
    mutable std::vector<TypeOfField> _types;
    mutable std::vector<INTERP_KERNEL::NormalizedCellType> _geo_types;
    mutable std::vector<double> _times;
  };

  /*!
   * User class.
   */
  class MEDFileField1TS : public MEDFileField1TSWithoutDAS, public MEDFieldFieldGlobs
  {
  public:
    static MEDFileField1TS *New(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void setFileName(const char *fileName);
    MEDCouplingFieldDouble *getFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(int meshDimRelToMaxExt, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
  private:
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    MEDFileField1TS(const char *fileName, const char *fieldName, int iteration, int order) throw(INTERP_KERNEL::Exception);
  };
  
  class MEDFileFieldMultiTSWithoutDAS : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldMultiTSWithoutDAS *New(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception);
    int getNumberOfTS() const;
  protected:
    const std::vector<std::string>& getInfo() const throw(INTERP_KERNEL::Exception);
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    std::string getDtUnit() const throw(INTERP_KERNEL::Exception);
    MEDFileFieldMultiTSWithoutDAS(const char *fieldName);
    MEDFileFieldMultiTSWithoutDAS(med_idt fid, const char *fieldName, int id, const std::vector<std::string>& infos, int nbOfStep) throw(INTERP_KERNEL::Exception);
    void finishLoading(med_idt fid, int nbPdt) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
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
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
  private:
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
  private:
    MEDFileFieldMultiTS(const char *fileName, const char *fieldName) throw(INTERP_KERNEL::Exception);
  };

  /*!
   * Use class.
   */
  class MEDFileFields : public RefCountObject, public MEDFieldFieldGlobs, public MEDFileWritable
  {
  public:
    static MEDFileFields *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    int getNumberOfFields() const;
  private:
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
  private:
    MEDFileFields(const char *fileName) throw(INTERP_KERNEL::Exception);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldMultiTSWithoutDAS> > _fields;
  };
}

#endif
