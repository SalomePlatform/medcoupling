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
  class MEDFileMesh;

  class MEDFileFieldLoc : public RefCountObject
  {
  public:
    const std::string& getName() const { return _name; }
    static MEDFileFieldLoc *New(med_idt fid, const char *locName);
    int getNbOfGaussPtPerCell() const { return _nb_gauss_pt; }
    void writeLL(med_idt fid) const;
    bool isName(const char *name) const { return _name==name; }
    const std::vector<double>& getRefCoords() const { return _ref_coo; }
    const std::vector<double>& getGaussCoords() const { return _gs_coo; }
    const std::vector<double>& getGaussWeights() const { return _w; }
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
    static MEDFileFieldPerMeshPerTypePerDisc *New(MEDFileFieldPerMeshPerType *fath, med_idt fid, TypeOfField type, int profileIt) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMeshPerType *getFather() const;
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getName() const;
    std::string getMeshName() const;
    TypeOfField getType() const;
    void setType(TypeOfField newType);
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    int getNumberOfTuples() const;
    const std::vector<std::string>& getInfo() const;
    std::string getProfile() const;
    std::string getLocalization() const;
    void getFieldAtLevel(TypeOfField type, const MEDFieldFieldGlobs *glob, std::vector<const DataArrayDouble *>& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs,
                         std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
  private:
    MEDFileFieldPerMeshPerTypePerDisc(MEDFileFieldPerMeshPerType *fath, med_idt fid, TypeOfField type, int profileIt) throw(INTERP_KERNEL::Exception);
    
  private:
    TypeOfField _type;
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
    static MEDFileFieldPerMeshPerType *New(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
    const MEDFileFieldPerMesh *getFather() const;
    void finishLoading(med_idt fid, TypeOfField type) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    void getDimension(int& dim) const;
    int getIteration() const;
    int getOrder() const;
    double getTime() const;
    std::string getName() const;
    std::string getMeshName() const;
    INTERP_KERNEL::NormalizedCellType getGeoType() const;
    int getNumberOfComponents() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    void getFieldAtLevel(int meshDim, TypeOfField type, const MEDFieldFieldGlobs *glob, std::vector<const DataArrayDouble *>& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs, std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes) const;
    static med_entity_type ConvertIntoMEDFileType(TypeOfField ikType, INTERP_KERNEL::NormalizedCellType ikGeoType, med_geometry_type& medfGeoType);
  private:
    MEDFileFieldPerMeshPerType(MEDFileFieldPerMesh *fath, INTERP_KERNEL::NormalizedCellType geoType) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileFieldPerMesh *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileFieldPerMeshPerTypePerDisc> > _field_pm_pt_pd;
    INTERP_KERNEL::NormalizedCellType _geo_type;
  };

  class MEDFileFieldPerMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileFieldPerMesh *New(MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder);
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    void getDimension(int& dim) const;
    double getTime() const;
    int getIteration() const;
    int getOrder() const;
    int getMeshIteration() const { return _mesh_iteration; }
    int getMeshOrder() const { return _mesh_order; }
    const std::string& getDtUnit() const;
    std::string getName() const;
    std::string getMeshName() const { return _mesh_name; }
    int getNumberOfComponents() const;
    const std::vector<std::string>& getInfo() const;
    std::vector<std::string> getPflsReallyUsed() const;
    std::vector<std::string> getLocsReallyUsed() const;
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDim, const MEDFieldFieldGlobs *glob, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
  private:
    MEDCouplingFieldDouble *finishField(TypeOfField type, const MEDFieldFieldGlobs *glob,
                                        const std::vector<const DataArrayDouble *>& dads, const std::vector<int>& locs, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *finishField2(TypeOfField type, const MEDFieldFieldGlobs *glob,
                                         const std::vector<const DataArrayDouble *>& dads, const std::vector<const DataArrayInt *>& pfls, const std::vector<int>& locs,
                                         const MEDCouplingMesh *mesh, const DataArrayInt *da) const throw(INTERP_KERNEL::Exception);
    static void SortArraysPerType(const MEDFieldFieldGlobs *glob,
                                  std::vector<INTERP_KERNEL::NormalizedCellType>& geoTypes, std::vector<const DataArrayDouble *>& dads, std::vector<const DataArrayInt *>& pfls, std::vector<int>& locs,
                                  std::vector<int>& code, std::vector<DataArrayInt *>& notNullPfls);
    static int ComputeNbOfElems(const MEDFieldFieldGlobs *glob, const std::vector<const DataArrayDouble *>& dads, const std::vector<int>& locs) throw(INTERP_KERNEL::Exception);
    MEDFileFieldPerMesh(MEDFileField1TSWithoutDAS *fath, int meshCsit, int meshIteration, int meshOrder);
  private:
    std::string _mesh_name;
    int _mesh_iteration;
    int _mesh_order;
    int _mesh_csit;
    MEDFileField1TSWithoutDAS *_father;
    std::vector< MEDCouplingAutoRefCountObjectPtr< MEDFileFieldPerMeshPerType > > _field_pm_pt;
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
    int getNbOfGaussPtPerCell(int locId) const throw(INTERP_KERNEL::Exception);
    int getLocalisationId(const char *loc) const throw(INTERP_KERNEL::Exception);
    const char *getFileName() const { return _file_name.c_str(); }
    const MEDFileFieldLoc& getLocalizationFromId(int locId) const throw(INTERP_KERNEL::Exception);
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
    int getDimension() const;
    int getIteration() const { return _iteration; }
    int getOrder() const { return _order; }
    double getTime() const { return _dt; }
    std::string getName() const { return _name; }
    const std::string& getDtUnit() const { return _dt_unit; }
    std::string getMeshName() const throw(INTERP_KERNEL::Exception);
    int getMeshIteration() const throw(INTERP_KERNEL::Exception);
    int getMeshOrder() const throw(INTERP_KERNEL::Exception);
    int getNumberOfComponents() const { return _infos.size(); }
    bool isDealingTS(int iteration, int order) const;
    std::pair<int,int> getDtIt() const;
    void fillIteration(std::pair<int,int>& p) const;
    const std::vector<std::string>& getInfo() const { return _infos; }
    //
    static MEDFileField1TSWithoutDAS *New(const char *fieldName, int csit, int iteration, int order, const std::vector<std::string>& infos);
    void finishLoading(med_idt fid) throw(INTERP_KERNEL::Exception);
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getPflsReallyUsed2() const;
    std::vector<std::string> getLocsReallyUsed2() const;
    static void CheckMeshDimRel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception);
  public:
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFieldFieldGlobs *glob) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFieldFieldGlobs *glob, const MEDCouplingMesh *mesh) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFieldFieldGlobs *glob, const MEDFileMesh *mesh) const throw(INTERP_KERNEL::Exception);
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
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
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
    std::vector< std::pair<int,int> > getIterations() const;
  public:
    std::vector<std::string> getPflsReallyUsed2() const;
    std::vector<std::string> getLocsReallyUsed2() const;
  protected:
    const MEDFileField1TSWithoutDAS& getTimeStepEntry(int iteration, int order) const throw(INTERP_KERNEL::Exception);
    MEDFileField1TSWithoutDAS& getTimeStepEntry(int iteration, int order) throw(INTERP_KERNEL::Exception);
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
    MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
    MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const throw(INTERP_KERNEL::Exception);
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
