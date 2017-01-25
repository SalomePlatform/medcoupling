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
// Author : Anthony Geay (EDF R&D)

#ifndef __MEDFILEEQUIVALENCE_HXX__
#define __MEDFILEEQUIVALENCE_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDFileUtilities.txx"
#include "MCAuto.hxx"

#include "NormalizedGeometricTypes"

#include <vector>

namespace MEDCoupling
{
  class MEDFileEquivalenceCell;
  class MEDFileEquivalenceNode;
  class MEDFileEquivalences;
  class MEDFileMesh;

  class MEDFileEquivalencePair : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    static MEDFileEquivalencePair *Load(MEDFileEquivalences *father, med_idt fid, const std::string& name, const std::string &desc);
    void writeLL(med_idt fid) const;
    const MEDFileEquivalences *getFather() const { return _father; }
    MEDFileEquivalences *getFather() { return _father; }
    const MEDFileMesh *getMesh() const;
    MEDFileMesh *getMesh();
    MEDFileEquivalencePair *deepCopy(MEDFileEquivalences *father) const;
    bool isEqual(const MEDFileEquivalencePair *other, std::string& what) const;
    void getRepr(std::ostream& oss) const;
    static MEDFileEquivalencePair *New(MEDFileEquivalences *father, const std::string& name);
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
  public:
    MEDLOADER_EXPORT std::string getName() const { return _name; }
    MEDLOADER_EXPORT void setName(const std::string& name) { _name=name; }
    MEDLOADER_EXPORT std::string getDescription() const { return _description; }
    MEDLOADER_EXPORT void setDescription(const std::string& descr) { _description=descr; }
    MEDLOADER_EXPORT MEDFileEquivalenceCell *initCell();
    MEDLOADER_EXPORT MEDFileEquivalenceNode *initNode();
    MEDLOADER_EXPORT MEDFileEquivalenceCell *getCell() { return _cell; }
    MEDLOADER_EXPORT MEDFileEquivalenceNode *getNode() { return _node; }
    MEDLOADER_EXPORT void setArray(int meshDimRelToMaxExt, DataArrayInt *da);
  private:
    MEDFileEquivalencePair(MEDFileEquivalences *father, const std::string& name, const std::string& desc):_father(father),_name(name),_description(desc) { }
    void load(med_idt fid);
  private:
    MEDFileEquivalences *_father;
    std::string _name;
    std::string _description;
    MCAuto<MEDFileEquivalenceCell> _cell;
    MCAuto<MEDFileEquivalenceNode> _node;
  };

  class MEDFileEquivalences : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT const MEDFileMesh *getMesh() const { return _owner; }
    MEDLOADER_EXPORT MEDFileMesh *getMesh() { return _owner; }
    void getDtIt(int &dt, int &it) const;
    std::string getMeshName() const;
    void pushEquivalence(MEDFileEquivalencePair *elt);
    static MEDFileEquivalences *New(MEDFileMesh *owner) { return new MEDFileEquivalences(owner); }
    MEDFileEquivalences *deepCopy(MEDFileMesh *owner) const;
    bool isEqual(const MEDFileEquivalences *other, std::string& what) const;
    void getRepr(std::ostream& oss) const;
  public:
    MEDLOADER_EXPORT MEDFileEquivalencePair *getEquivalence(int i);
    MEDLOADER_EXPORT MEDFileEquivalencePair *getEquivalenceWithName(const std::string& name);
    MEDLOADER_EXPORT int size() const;
    MEDLOADER_EXPORT std::vector<std::string> getEquivalenceNames() const;
    MEDLOADER_EXPORT MEDFileEquivalencePair *appendEmptyEquivalenceWithName(const std::string& name);
    MEDLOADER_EXPORT void killEquivalenceWithName(const std::string& name);
    MEDLOADER_EXPORT void killEquivalenceAt(int i);
    MEDLOADER_EXPORT void clear();
  public:
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    static int PresenceOfEquivalences(med_idt fid, const std::string& meshName);
    static MEDFileEquivalences *Load(med_idt fid, int nbOfEq, MEDFileMesh *owner);
    static void CheckDataArray(const DataArrayInt *data);
  private:
    MEDFileEquivalences(MEDFileMesh *owner):_owner(owner) { }
    void deepCpyFrom(const MEDFileEquivalences& other);
  private:
    MEDFileMesh *_owner;
    std::vector< MCAuto<MEDFileEquivalencePair> > _equ;
  };

  class MEDFileEquivalenceBase : public RefCountObject, public MEDFileWritableStandAlone
  {
  protected:
    MEDFileEquivalenceBase(MEDFileEquivalencePair *father);
    const MEDFileEquivalencePair *getFather() const { return _father; }
    MEDFileEquivalencePair *getFather() { return _father; }
    const MEDFileMesh *getMesh() const { return getFather()->getMesh(); }
    MEDFileMesh *getMesh() { return getFather()->getMesh(); }
  protected:
    ~MEDFileEquivalenceBase() { }
  private:
    MEDFileEquivalencePair *_father;
  };

  class MEDFileEquivalenceData : public MEDFileEquivalenceBase
  {
  public:
    MEDFileEquivalenceData(MEDFileEquivalencePair *owner, DataArrayInt *data);
    MEDLOADER_EXPORT void setArray(DataArrayInt *data);
    MEDLOADER_EXPORT const DataArrayInt *getArray() const { return _data; }
    MEDLOADER_EXPORT DataArrayInt *getArray() { return _data; }
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    bool isEqual(const MEDFileEquivalenceData *other, std::string& what) const;
  protected:
    void writeAdvanced(med_idt fid, med_entity_type medtype, med_geometry_type medgt) const;
  protected:
    ~MEDFileEquivalenceData() { }
  protected:
    MCAuto<DataArrayInt> _data;
  };

  class MEDFileEquivalenceCellType : public MEDFileEquivalenceData
  {
  public:
    MEDFileEquivalenceCellType(MEDFileEquivalencePair *owner, INTERP_KERNEL::NormalizedCellType type, DataArrayInt *data):MEDFileEquivalenceData(owner,data),_type(type) { }
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    INTERP_KERNEL::NormalizedCellType getType() const { return _type; }
    MEDFileEquivalenceCellType *deepCopy(MEDFileEquivalencePair *owner) const;
    bool isEqual(const MEDFileEquivalenceCellType *other, std::string& what) const;
    void getRepr(std::ostream& oss) const;
  public:
    void writeLL(med_idt fid) const;
  protected:
    ~MEDFileEquivalenceCellType() { }
  private:
    INTERP_KERNEL::NormalizedCellType _type;
  };

  class MEDFileEquivalenceCell : public MEDFileEquivalenceBase
  {
  public:
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    static MEDFileEquivalenceCell *Load(med_idt fid, MEDFileEquivalencePair *owner);
    void writeLL(med_idt fid) const;
    MEDFileEquivalenceCell *deepCopy(MEDFileEquivalencePair *owner) const;
    bool isEqual(const MEDFileEquivalenceCell *other, std::string& what) const;
    void getRepr(std::ostream& oss) const;
  public:
    MEDLOADER_EXPORT void clear() { _types.clear(); }
    MEDLOADER_EXPORT std::size_t size() const { return _types.size(); }
    MEDLOADER_EXPORT DataArrayInt *getArray(INTERP_KERNEL::NormalizedCellType type);
    MEDLOADER_EXPORT void setArray(int meshDimRelToMax, DataArrayInt *da);
    MEDLOADER_EXPORT void setArrayForType(INTERP_KERNEL::NormalizedCellType type, DataArrayInt *da);
    MEDLOADER_EXPORT std::vector<INTERP_KERNEL::NormalizedCellType> getTypes() const;
  public:
    MEDFileEquivalenceCell(MEDFileEquivalencePair *owner):MEDFileEquivalenceBase(owner) { }
  private:
    ~MEDFileEquivalenceCell() { }
  private:
    void load(med_idt fid);
    std::string getName() const { return getFather()->getName(); }
  private:
    std::vector< MCAuto<MEDFileEquivalenceCellType> > _types;
  };

  class MEDFileEquivalenceNode : public MEDFileEquivalenceData
  {
  public:
    MEDFileEquivalenceNode(MEDFileEquivalencePair *owner, DataArrayInt *data):MEDFileEquivalenceData(owner,data) { }
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    void writeLL(med_idt fid) const;
    MEDFileEquivalenceNode *deepCopy(MEDFileEquivalencePair *owner) const;
    bool isEqual(const MEDFileEquivalenceNode *other, std::string& what) const;
    void getRepr(std::ostream& oss) const;
  protected:
    ~MEDFileEquivalenceNode() { }
  };
}

#endif

