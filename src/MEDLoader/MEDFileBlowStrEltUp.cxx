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

#include "MEDFileBlowStrEltUp.hxx"
#include "MEDCouplingFieldDouble.hxx"
#include "MEDFileFieldVisitor.hxx"
#include "MEDCouplingPartDefinition.hxx"
#include "MCAuto.txx"

using namespace MEDCoupling;

const char MEDFileBlowStrEltUp::MED_BALL_STR[]="MED_BALL";

MEDFileBlowStrEltUp::MEDFileBlowStrEltUp(const MEDFileFields *fsOnlyOnSE, const MEDFileMeshes *ms, const MEDFileStructureElements *ses)
{
  if(!fsOnlyOnSE || !ms || !ses)
    throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp constructor : NULL input pointer !");
  _ms.takeRef(ms); _ses.takeRef(ses);
  std::vector< std::pair<std::string,std::string> > ps;
  fsOnlyOnSE->getMeshSENames(ps);
  std::size_t sz(ps.size());
  _elts.resize(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const std::pair<std::string,std::string>& p(ps[i]);
      MCAuto<MEDFileFields> f(fsOnlyOnSE->partOfThisLyingOnSpecifiedMeshSEName(p.first,p.second));
      _elts[i]=f;
    }
  for(std::size_t i=0;i<sz;i++)
    {
      const std::pair<std::string,std::string>& p(ps[i]);
      MEDFileMesh *mesh(_ms->getMeshWithName(p.first));
      if(!mesh)
        throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp : NULL mesh !");
      MEDFileUMesh *umesh(dynamic_cast<MEDFileUMesh *>(mesh));
      if(!umesh)
        throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp : Blow up of Stru Elt not managed yet for unstructured meshes !");
    }
}

/*!
 * \param [in] mesh - The mesh containing structure element called \a seName. After the call of this method the Structure elements parts will be removed.
 * \param [out] mOut - the physical mesh of the structure element \a seName in mesh \a mesh
 * \param [out] fsOut - the list of var attribute of structure element \a seName - \b WARNING no time steps here
 */
MCAuto<MEDFileEltStruct4Mesh> MEDFileBlowStrEltUp::dealWithSEInMesh(const std::string& seName, MEDFileUMesh *mesh, MCAuto<MEDFileUMesh>& mOut, MCAuto<MEDFileFields>& fsOut) const
{
  if(!mesh)
    throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithSEInMesh : null pointer !");
  if(seName==MED_BALL_STR)
    {
      MCAuto<MEDFileEltStruct4Mesh> ret(dealWithMEDBALLInMesh(mesh,mOut,fsOut));
      mesh->killStructureElements();
      return ret;
    }
  throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithSEInMesh : only MED_BALL is managed for the moment, but if you are interested please send spec to anthony.geay@edf.fr !");
}

MCAuto<MEDFileEltStruct4Mesh> MEDFileBlowStrEltUp::dealWithMEDBALLInMesh(const MEDFileUMesh *mesh, MCAuto<MEDFileUMesh>& mOut, MCAuto<MEDFileFields>& fsOut) const
{
  mOut=MEDFileUMesh::New(); fsOut=MEDFileFields::New();
  const std::vector< MCAuto<MEDFileEltStruct4Mesh> >& strs(mesh->getAccessOfUndergroundEltStrs());
  MCAuto<MEDFileEltStruct4Mesh> zeStr;
  for(std::vector< MCAuto<MEDFileEltStruct4Mesh> >::const_iterator it=strs.begin();it!=strs.end();it++)
    {
      if((*it)->getGeoTypeName()==MED_BALL_STR)
        {
          zeStr=*it;
          break;
        }
    }
  if(zeStr.isNull())
    {
      std::ostringstream oss; oss << "MEDFileBlowStrEltUp::dealWithMEDBALLInMesh : no geo type with name " <<  MED_BALL_STR << " in " << mesh->getName() << " !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
  const DataArrayDouble *coo(mesh->getCoords());
  if(!coo)
    throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithMEDBALLInMesh : null coords !");
  MCAuto<DataArrayInt> conn(zeStr->getConn());
  if(conn.isNull())
    throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithMEDBALLInMesh : null connectivity !");
  conn->checkAllocated();
  if(conn->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithMEDBALLInMesh : excepted to be single compo !");
  int nbCells(conn->getNumberOfTuples());
  MCAuto<DataArrayDouble> connOut(coo->selectByTupleIdSafe(conn->begin(),conn->end()));
  MCAuto<MEDCouplingUMesh> mcOut(MEDCouplingUMesh::Build0DMeshFromCoords(connOut));
  mcOut->setName(BuildNewMeshName(mesh->getName(),MED_BALL_STR));
  mOut->setMeshAtLevel(0,mcOut);
  const DataArrayInt *ff1(mesh->getFamilyFieldAtLevel(1));
  if(ff1)
    {
      MCAuto<DataArrayInt> ff1o(ff1->selectByTupleIdSafe(conn->begin(),conn->end()));
      mOut->setFamilyFieldArr(1,ff1o);
    }
  const DataArrayInt *nf1(mesh->getNumberFieldAtLevel(1));
  if(nf1)
    {
      MCAuto<DataArrayInt> nf1o(nf1->selectByTupleIdSafe(conn->begin(),conn->end()));
      mOut->setRenumFieldArr(1,nf1o);
    }
  MCAuto<MEDFileUMeshPerTypeCommon> md(zeStr->getMeshDef());
  const DataArrayInt *ff0(md->getFam());
  if(ff0)
    mOut->setFamilyFieldArr(0,const_cast<DataArrayInt *>(ff0));
  const DataArrayInt *nf0(md->getNum());
  if(nf0)
    mOut->setRenumFieldArr(0,const_cast<DataArrayInt *>(nf0));
  mOut->copyFamGrpMapsFrom(*mesh);
  const std::vector< MCAuto<DataArray> >& vars(zeStr->getVars());
  for(std::vector< MCAuto<DataArray> >::const_iterator it=vars.begin();it!=vars.end();it++)
    {
      const DataArray *elt(*it);
      if(!elt)
        continue;
      {
        const DataArrayDouble *eltC(dynamic_cast<const DataArrayDouble *>(elt));
        if(eltC)
          {
            MCAuto<MEDFileFieldMultiTS> fmts(MEDFileFieldMultiTS::New());
            MCAuto<MEDFileField1TS> f1ts(MEDFileField1TS::New());
            MCAuto<MEDCouplingFieldDouble> f(MEDCouplingFieldDouble::New(ON_NODES));
            f->setMesh(mcOut);
            f->setArray(const_cast<DataArrayDouble *>(eltC));
            f->setName(eltC->getName());
            f1ts->setFieldNoProfileSBT(f);
            fmts->pushBackTimeStep(f1ts);
            fsOut->pushField(fmts);
          }
      }
    }
  return zeStr;
}

/*!
 * \param [in] fs - fields lying all on same mesh and on same structure element
 * \param [in] zeStr - ze structure of current structure element
 * \param [in] varAtt - fields containing var att of current structure element. WARNING at this stage the number of iteration are equal to one for each field in \a varAtt
 * \param [out] zeOutputs - ze fields that are the concatenation of fields in \a fs transformed and those in \a varAtt normalized in time space
 */
void MEDFileBlowStrEltUp::dealWithSEInFields(const std::string& seName, const MEDFileFields *fs, const MEDFileEltStruct4Mesh *zeStr, const MEDFileFields *varAtt, MEDFileFields *zeOutputs) const
{
  if(!fs)
    throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithSEInFields : null pointer !");
  if(seName==MED_BALL_STR)
    {
      dealWithMEDBALLSInFields(fs,zeStr,varAtt,zeOutputs);
      return ;
    }
  throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithSEInFields : only MED_BALL is managed for the moment, but if you are interested please send spec to anthony.geay@edf.fr !");
}

void MEDFileBlowStrEltUp::dealWithMEDBALLSInFields(const MEDFileFields *fs, const MEDFileEltStruct4Mesh *zeStr, const MEDFileFields *varAtt, MEDFileFields *zeOutputs) const
{
  int nbf(fs->getNumberOfFields());
  std::vector< MCAuto<MEDFileAnyTypeFieldMultiTS> > elts0;
  std::vector< MEDFileAnyTypeFieldMultiTS * > elts1;
  std::string zeMeshName;
  for(int i=0;i<nbf;i++)
    {
      MCAuto<MEDFileAnyTypeFieldMultiTS> elt(fs->getFieldAtPos(i));
      MCAuto<MEDFileAnyTypeFieldMultiTS> eltOut(elt->buildNewEmpty());
      int nbTS(elt->getNumberOfTS());
      for(int j=0;j<nbTS;j++)
        {
          MCAuto<MEDFileAnyTypeField1TS> eltt(elt->getTimeStepAtPos(j));
          MCAuto<MEDFileAnyTypeField1TS> elttOut(eltt->deepCopy());
          std::string meshName(eltt->getMeshName());
          zeMeshName=BuildNewMeshName(meshName,MED_BALL_STR);
          elttOut->setMeshName(zeMeshName);
          elttOut->convertMedBallIntoClassic();
          eltOut->pushBackTimeStep(elttOut);
        }
      elts0.push_back(eltOut); elts1.push_back(eltOut);
    }
  //
  const MEDFileMesh *zeCurrentMesh(_ms->getMeshWithName(zeMeshName));
  //
  std::size_t ii(0);
  std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > sp(MEDFileAnyTypeFieldMultiTS::SplitIntoCommonTimeSeries(elts1));
  for(std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> >::const_iterator it0=sp.begin();it0!=sp.end();it0++,ii++)
    {
      std::vector< MCAuto<MEDFileFastCellSupportComparator> > fsc;
      std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > sp2(MEDFileAnyTypeFieldMultiTS::SplitPerCommonSupport(*it0,zeCurrentMesh,fsc));
      std::size_t jj(0);
      for(std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> >::const_iterator it1=sp2.begin();it1!=sp2.end();it1++,jj++)
        {
          for(std::vector<MEDFileAnyTypeFieldMultiTS *>::const_iterator it2=(*it1).begin();it2!=(*it1).end();it2++)
            zeOutputs->pushField(*it2);
          // The most exciting part. Users that put profiles on struct elements part of fields. Reduce var att.
          if((*it1).size()<1)
            throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithMEDBALLSInFields : take a deep breath !");
          MCAuto<MEDFileAnyTypeField1TS> zeGuideForPfl;// This var is the reference for pfl management.
          {
            if(!(*it1)[0])
              throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithMEDBALLSInFields : take a deep breath 2 !");
            int pdm((*it1)[0]->getNumberOfTS());
            if(pdm<1)
              throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithMEDBALLSInFields : take a deep breath 3 !");
            zeGuideForPfl=(*it1)[0]->getTimeStepAtPos(0);
          }
          if(zeGuideForPfl.isNull())
            throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithMEDBALLSInFields : take a deep breath 4 !");
          std::vector<std::string> pfls(zeGuideForPfl->getPflsReallyUsed());
          if(pfls.size()>=2)
            throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithMEDBALLSInFields : drink less coffee");
          MCAuto<DataArrayInt> pflMyLove;
          if(pfls.size()==1)
            pflMyLove.takeRef(zeGuideForPfl->getProfile(pfls[0]));
          // Yeah we have pfls
          std::vector<double> t2s;
          std::vector< std::pair<int,int> > t1s((*it1)[0]->getTimeSteps(t2s));
          std::size_t nbTS3(t2s.size());
          int nbf2(varAtt->getNumberOfFields());
          for(int i=0;i<nbf2;i++)
            {
              MCAuto<MEDFileAnyTypeFieldMultiTS> elt(varAtt->getFieldAtPos(i));
              int nbTS2(elt->getNumberOfTS());
              if(nbTS2!=1)
                throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::dealWithMEDBALLSInFields : internal error ! The dealWithMEDBALLInMesh is expected to return a single TS !");
              MCAuto<MEDFileAnyTypeField1TS> elt2(elt->getTimeStepAtPos(0));
              MCAuto<MEDFileAnyTypeFieldMultiTS> elt4(elt->buildNewEmpty());
              for(std::size_t j=0;j<nbTS3;j++)
                {
                  MCAuto<MEDFileAnyTypeField1TS> elt3(elt2->deepCopy());
                  elt3->setTime(t1s[j].first,t1s[j].second,t2s[j]);
                  elt3->setName(BuildVarAttName(ii,sp.size(),jj,sp2.size(),elt3->getName()));
                  if(pflMyLove.isNotNull())
                    elt3->makeReduction(INTERP_KERNEL::NORM_ERROR,ON_NODES,pflMyLove);
                  elt4->pushBackTimeStep(elt3);
                }
              zeOutputs->pushField(elt4);
            }
        }
    }
}

void MEDFileBlowStrEltUp::generate(MEDFileMeshes *msOut, MEDFileFields *allZeOutFields)
{
  for(std::vector< MCAuto<MEDFileFields> >::iterator elt=_elts.begin();elt!=_elts.end();elt++)
    {
      std::vector< std::pair<std::string,std::string> > ps;
      (*elt)->getMeshSENames(ps);
      if(ps.size()!=1)
        throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::generateMeshes : internal error !");
      MEDFileMesh *mesh(_ms->getMeshWithName(ps[0].first));
      if(!mesh)
        throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::generateMeshes : NULL mesh !");
      MEDFileUMesh *umesh(dynamic_cast<MEDFileUMesh *>(mesh));
      if(!umesh)
        throw INTERP_KERNEL::Exception("MEDFileBlowStrEltUp::generateMeshes : Blow up of Stru Elt not managed yet for unstructured meshes !");
      //
      MCAuto<MEDFileFields> classicalSEFields(splitFieldsPerLoc(*elt,umesh,msOut,allZeOutFields));
      if(classicalSEFields.isNotNull())
        {
          MCAuto<MEDFileUMesh> mOut;
          MCAuto<MEDFileFields> fsOut1;
          MCAuto<MEDFileEltStruct4Mesh> zeStr(dealWithSEInMesh(ps[0].second,umesh,mOut,fsOut1));
          msOut->pushMesh(mOut);
          dealWithSEInFields(ps[0].second,classicalSEFields,zeStr,fsOut1,allZeOutFields);
        }
    }
}

std::string MEDFileBlowStrEltUp::BuildNewMeshName(const std::string& meshName, const std::string& seName)
{
  std::ostringstream mNameOut;
  mNameOut << meshName << "_" << seName;
  return mNameOut.str();
}

std::string MEDFileBlowStrEltUp::BuildVarAttName(std::size_t iPart, std::size_t totINbParts, std::size_t jPart, std::size_t totJNbParts, const std::string& name)
{
  if(totINbParts==1 && totJNbParts==1)
    return name;
  std::ostringstream oss;
  oss << name << "@" << iPart << "@" << jPart;
  return oss.str();
}

void MEDFileBlowStrEltUp::DealWithSE(MEDFileFields *fs, MEDFileMeshes *ms, const MEDFileStructureElements *ses)
{
  MCAuto<MEDFileFields> fsSEOnly(fs->partOfThisOnStructureElements());
  fs->killStructureElements();
  MEDFileBlowStrEltUp bu(fsSEOnly,ms,ses);
  bu.generate(ms,fs);
  fs->killStructureElementsInGlobs();
}

//

class FieldWalker2
{
public:
  FieldWalker2(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd);
  std::string getLoc() const { return _loc; }
  std::string getPfl() const { return _pfl; }
  bool isClassic() const { return _is_classic; }
  bool operator!=(const FieldWalker2& other) const;
  bool operator==(const FieldWalker2& other) const;
  const SlicePartDefinition *getPartDef() const { return _pd; }
private:
  std::string _loc;
  std::string _pfl;
  bool _is_classic;
  MCAuto<SlicePartDefinition> _pd;
};

class LocInfo
{
public:
  LocInfo() { }
  LocInfo(const std::vector<FieldWalker2>& fw);
  bool operator==(const LocInfo& other) const { return _locs==other._locs && _pfl==other._pfl; }
  void push(const std::string& loc, const std::string& pfl) { checkUniqueLoc(loc); _locs.push_back(loc); _pfl.push_back(pfl); }
  MCAuto<MEDFileUMesh> generateNonClassicalData(int zePos, const MEDFileUMesh *mesh, const MEDFileFieldGlobsReal *globs) const;
  const PartDefinition *getPartDef() const { return _pd; }
private:
  void checkUniqueLoc(const std::string& loc) const;
  static MCAuto<DataArrayDouble> BuildMeshFromAngleVrille(INTERP_KERNEL::NormalizedCellType gt, const DataArrayDouble *angleDeVrille, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs);
  static MCAuto<DataArrayDouble> BuildMeshFromEpaisseur(INTERP_KERNEL::NormalizedCellType gt, const DataArrayDouble *thikness, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs);
  static MCAuto<DataArrayDouble> BuildMeshPipeSEG3(const DataArrayDouble *angle, const DataArrayDouble *scale, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs);
  static MCAuto<MEDCouplingUMesh> BuildMeshCommon(INTERP_KERNEL::NormalizedCellType gt, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs, MCAuto<DataArrayDouble>& ptsForLoc);
  static MCAuto<DataArrayDouble> BuildMeshFromStructure(INTERP_KERNEL::NormalizedCellType gt, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs);
public:
  static const char ANGLE_DE_VRILLE[];
  static const char ANGLE[];
  static const char SCALE[];
  static const char EPAISSEUR[];
private:
  std::vector<std::string> _locs;
  std::vector<std::string> _pfl;
  MCAuto<PartDefinition> _pd;
};

const char LocInfo::ANGLE_DE_VRILLE[]="ANGLE DE VRILLE";

const char LocInfo::ANGLE[]="ANGLE";

const char LocInfo::SCALE[]="SCALE";

const char LocInfo::EPAISSEUR[]="EPAISSEUR";

LocInfo::LocInfo(const std::vector<FieldWalker2>& fw)
{
  std::size_t sz(fw.size());
  _locs.resize(sz); _pfl.resize(sz);
  if(sz>0)
    _pd=fw[0].getPartDef()->deepCopy();
  for(std::size_t i=0;i<sz;i++)
    {
      _locs[i]=fw[i].getLoc();
      _pfl[i]=fw[i].getPfl();
      if(i>0)
        _pd=(*_pd)+(*(fw[i].getPartDef()));
    }
}

void LocInfo::checkUniqueLoc(const std::string& loc) const
{
  if(std::find(_locs.begin(),_locs.end(),loc)!=_locs.end())
    {
      std::ostringstream oss; oss << "LocInfo::checkUniqueLoc : loc \"" << loc << "\" already exists !";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

MCAuto<MEDCouplingUMesh> LocInfo::BuildMeshCommon(INTERP_KERNEL::NormalizedCellType gt, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs, MCAuto<DataArrayDouble>& ptsForLoc)
{
  MCAuto<DataArrayInt> conn(zeStr->getConn());
  conn=conn->deepCopy(); conn->rearrange(1);
  MCAuto<MEDCouplingUMesh> geoMesh;
  {
    MCAuto<MEDCoupling1SGTUMesh> umesh(MEDCoupling1SGTUMesh::New("",gt));
    umesh->setCoords(mesh->getCoords());
    umesh->setNodalConnectivity(conn);
    geoMesh=umesh->buildUnstructured();
  }
  //
  if(!pfl.empty())
    {
      const DataArrayInt *pflArr(globs->getProfile(pfl));
      geoMesh=geoMesh->buildPartOfMySelf(pflArr->begin(),pflArr->end(),true);
    }
  //
  MCAuto<MEDCouplingFieldDouble> fakeF(MEDCouplingFieldDouble::New(ON_GAUSS_PT));
  fakeF->setMesh(geoMesh);
  fakeF->setGaussLocalizationOnType(gt,loc.getRefCoords(),loc.getGaussCoords(),loc.getGaussWeights());
  ptsForLoc=fakeF->getLocalizationOfDiscr();
  //
  return geoMesh;
}

MCAuto<DataArrayDouble> LocInfo::BuildMeshFromAngleVrille(INTERP_KERNEL::NormalizedCellType gt, const DataArrayDouble *angleDeVrille, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs)
{
  MCAuto<DataArrayDouble> ptsForLoc;
  MCAuto<MEDCouplingUMesh> geoMesh(BuildMeshCommon(gt,pfl,loc,zeStr,mesh,section,globs,ptsForLoc));
  //
  MCConstAuto<DataArrayDouble> angleVrille;
  if(!pfl.empty())
    {
      const DataArrayInt *pflArr(globs->getProfile(pfl));
      angleVrille=angleDeVrille->selectByTupleIdSafe(pflArr->begin(),pflArr->end());
    }
  else
    angleVrille.takeRef(angleDeVrille);
  //
  MCAuto<MEDCouplingFieldDouble> dir(geoMesh->buildDirectionVectorField());
  MCAuto<DataArrayDouble> rot(dir->getArray()->fromCartToSpher());
  int nbCompo(ptsForLoc->getNumberOfComponents());
  MCAuto<DataArrayDouble> secPts(section->getCoords()->changeNbOfComponents(nbCompo,0.));
  int nbSecPts(secPts->getNumberOfTuples()),nbCells(geoMesh->getNumberOfCells()),nbg(loc.getGaussWeights().size());
  {
    const int TAB[3]={2,0,1};
    std::vector<int> v(TAB,TAB+3);
    secPts=secPts->keepSelectedComponents(v);
  }
  const double CENTER[3]={0.,0.,0.},AX0[3]={0.,0.,1.};
  double AX1[3]; AX1[2]=0.;
  std::vector< MCAuto<DataArrayDouble> > arrs(nbCells*nbg);
  for(int j=0;j<nbCells;j++)
    {
      MCAuto<DataArrayDouble> p(secPts->deepCopy());
      double ang0(rot->getIJ(j,2));
      DataArrayDouble::Rotate3DAlg(CENTER,AX0,ang0,nbSecPts,p->begin(),p->getPointer());
      AX1[0]=-sin(ang0); AX1[1]=cos(ang0);// rot Oy around OZ
      double ang1(M_PI/2.-rot->getIJ(j,1));
      DataArrayDouble::Rotate3DAlg(CENTER,AX1,-ang1,nbSecPts,p->begin(),p->getPointer());
      DataArrayDouble::Rotate3DAlg(CENTER,dir->getArray()->begin()+j*3,angleVrille->getIJ(j,0),nbSecPts,p->begin(),p->getPointer());
      for(int l=0;l<nbg;l++)
        {
          MCAuto<DataArrayDouble> p2(p->deepCopy());
          for(int k=0;k<nbCompo;k++)
            p2->applyLin(1.,ptsForLoc->getIJ(j*nbg+l,k),k);
          arrs[j*nbg+l]=p2;
        }
    }
  std::vector<const DataArrayDouble *> arrs2(VecAutoToVecOfCstPt(arrs));
  MCAuto<DataArrayDouble> resu(DataArrayDouble::Aggregate(arrs2));
  return resu;
}

MCAuto<DataArrayDouble> LocInfo::BuildMeshFromEpaisseur(INTERP_KERNEL::NormalizedCellType gt, const DataArrayDouble *thikness, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs)
{
  MCAuto<DataArrayDouble> ptsForLoc;
  MCAuto<MEDCouplingUMesh> geoMesh(BuildMeshCommon(gt,pfl,loc,zeStr,mesh,section,globs,ptsForLoc));
  int nbSecPts(section->getNumberOfNodes()),nbCells(geoMesh->getNumberOfCells()),nbg(loc.getGaussWeights().size());
  MCConstAuto<DataArrayDouble> zeThikness;
  if(!pfl.empty())
    {
      const DataArrayInt *pflArr(globs->getProfile(pfl));
      geoMesh=geoMesh->buildPartOfMySelf(pflArr->begin(),pflArr->end(),true);
      zeThikness=thikness->selectByTupleIdSafe(pflArr->begin(),pflArr->end());
    }
  else
    zeThikness.takeRef(thikness);
  MCAuto<DataArrayDouble> orthoArr;
  {
    MCAuto<MEDCouplingFieldDouble> ortho(geoMesh->buildOrthogonalField());
    orthoArr.takeRef(ortho->getArray());
  }
  int nbCompo(orthoArr->getNumberOfComponents());
  MCAuto<DataArrayDouble> secPts(section->getCoords()->duplicateEachTupleNTimes(nbCompo));
  secPts->rearrange(nbCompo);
  std::vector< MCAuto<DataArrayDouble> > arrs(nbCells*nbg);
  for(int j=0;j<nbCells;j++)
    {
      double thck(zeThikness->getIJ(j,0));
      MCAuto<DataArrayDouble> fact(DataArrayDouble::New()); fact->alloc(1,nbCompo);
      std::copy(orthoArr->begin()+j*nbCompo,orthoArr->begin()+(j+1)*nbCompo,fact->getPointer());
      std::transform(fact->begin(),fact->end(),fact->getPointer(),std::bind2nd(std::multiplies<double>(),thck/2.));
      MCAuto<DataArrayDouble> p(DataArrayDouble::Multiply(secPts,fact));
      for(int l=0;l<nbg;l++)
        {
          MCAuto<DataArrayDouble> p2(p->deepCopy());
          for(int k=0;k<nbCompo;k++)
            p2->applyLin(1.,ptsForLoc->getIJ(j*nbg+l,k),k);
          arrs[j*nbg+l]=p2;
        }
    }
  std::vector<const DataArrayDouble *> arrs2(VecAutoToVecOfCstPt(arrs));
  MCAuto<DataArrayDouble> resu(DataArrayDouble::Aggregate(arrs2));
  return resu;
}

MCAuto<DataArrayDouble> LocInfo::BuildMeshPipeSEG3(const DataArrayDouble *angle, const DataArrayDouble *scale, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs)
{
  static const char MSG1[]="BuildMeshPipeSEG3 : not recognized pattern ! Send mail to anthony.geay@edf.fr with corresponding MED file !";
  MCAuto<DataArrayDouble> ptsForLoc;
  MCAuto<MEDCouplingUMesh> geoMesh(BuildMeshCommon(INTERP_KERNEL::NORM_SEG3,pfl,loc,zeStr,mesh,section,globs,ptsForLoc));
  int nbSecPts(section->getNumberOfNodes()),nbCells(geoMesh->getNumberOfCells()),nbg(loc.getGaussWeights().size());
  MCConstAuto<DataArrayDouble> zeAngle,zeScale;
  if(!pfl.empty())
    {
      const DataArrayInt *pflArr(globs->getProfile(pfl));
      geoMesh=geoMesh->buildPartOfMySelf(pflArr->begin(),pflArr->end(),true);
      zeAngle=angle->selectByTupleIdSafe(pflArr->begin(),pflArr->end());
      zeScale=scale->selectByTupleIdSafe(pflArr->begin(),pflArr->end());
    }
  else
    {
      zeAngle.takeRef(angle);
      zeScale.takeRef(scale);
    }
  if(zeAngle->getNumberOfComponents()!=3 || zeScale->getNumberOfComponents()!=2 || nbg!=3)
    throw INTERP_KERNEL::Exception(MSG1);
  MCAuto<MEDCouplingFieldDouble> dir;
  {
    MCAuto<MEDCouplingUMesh> geoMesh2(geoMesh->deepCopy());
    geoMesh2->convertQuadraticCellsToLinear();
    dir=geoMesh2->buildDirectionVectorField();
  }
  MCAuto<DataArrayDouble> rot(dir->getArray()->fromCartToSpher());
  int nbCompo(ptsForLoc->getNumberOfComponents());
  MCAuto<DataArrayDouble> secPts(section->getCoords()->changeNbOfComponents(nbCompo,0.));
  {
    const int TAB[3]={2,0,1};
    std::vector<int> v(TAB,TAB+3);
    secPts=secPts->keepSelectedComponents(v);
  }
  const double CENTER[3]={0.,0.,0.},AX0[3]={0.,0.,1.};
  double AX1[3]; AX1[2]=0.;
  std::vector< MCAuto<DataArrayDouble> > arrs(nbCells*nbg);
  for(int j=0;j<nbCells;j++)
    {
      MCAuto<DataArrayDouble> p(secPts->deepCopy());
      double ang0(rot->getIJ(j,2));
      DataArrayDouble::Rotate3DAlg(CENTER,AX0,ang0,nbSecPts,p->begin(),p->getPointer());
      AX1[0]=-sin(ang0); AX1[1]=cos(ang0);// rot Oy around OZ
      double ang1(M_PI/2.-rot->getIJ(j,1));
      DataArrayDouble::Rotate3DAlg(CENTER,AX1,-ang1,nbSecPts,p->begin(),p->getPointer());
      for(int l=0;l<3;l++)
        {
          MCAuto<DataArrayDouble> p3(p->deepCopy());
          DataArrayDouble::Rotate3DAlg(CENTER,dir->getArray()->begin()+j*3,zeAngle->getIJ(j,l),nbSecPts,p3->begin(),p3->getPointer());
          MCAuto<DataArrayDouble> p2(p3->deepCopy());
          for(int k=0;k<nbCompo;k++)
            p2->applyLin(1.,ptsForLoc->getIJ(j*nbg+l,k),k);
          arrs[j*nbg+l]=p2;
        }
    }
  std::vector<const DataArrayDouble *> arrs2(VecAutoToVecOfCstPt(arrs));
  MCAuto<DataArrayDouble> resu(DataArrayDouble::Aggregate(arrs2));
  return resu;
}

MCAuto<DataArrayDouble> LocInfo::BuildMeshFromStructure(INTERP_KERNEL::NormalizedCellType gt, const std::string& pfl, const MEDFileFieldLoc& loc, const MEDFileEltStruct4Mesh *zeStr, const MEDFileUMesh *mesh, const MEDFileUMesh *section, const MEDFileFieldGlobsReal *globs)
{
  static const char MSG1[]="BuildMeshFromStructure : not recognized pattern ! Send mail to anthony.geay@edf.fr with corresponding MED file !";
  const std::vector< MCAuto<DataArray> >& vars(zeStr->getVars());
  if(vars.size()==1)
    {
      MCAuto<DataArray> zeArr(vars[0]);
      if(zeArr.isNull())
        throw INTERP_KERNEL::Exception(MSG1);
      MCAuto<DataArrayDouble> zeArr2(DynamicCast<DataArray,DataArrayDouble>(zeArr));
      if(zeArr2.isNull())
        throw INTERP_KERNEL::Exception(MSG1);
      if(zeArr2->getName()==ANGLE_DE_VRILLE || zeArr2->getName()==ANGLE)
        return BuildMeshFromAngleVrille(gt,zeArr2,pfl,loc,zeStr,mesh,section,globs);
      if(zeArr2->getName()==EPAISSEUR || zeArr2->getName()==SCALE)
        return BuildMeshFromEpaisseur(gt,zeArr2,pfl,loc,zeStr,mesh,section,globs);
    }
  if(vars.size()==2)
    {
      MCAuto<DataArray> zeArr0(vars[0]),zeArr1(vars[1]);
      if(zeArr0.isNull() || zeArr1.isNull())
        throw INTERP_KERNEL::Exception(MSG1);
      MCAuto<DataArrayDouble> zeArr00(DynamicCastSafe<DataArray,DataArrayDouble>(zeArr0)),zeArr11(DynamicCastSafe<DataArray,DataArrayDouble>(zeArr1));
      switch(gt)
      {
        case INTERP_KERNEL::NORM_SEG3:
          {
            MCAuto<DataArrayDouble> angle,scale;
            if(zeArr00->getName()==ANGLE)
              angle=zeArr00;
            if(zeArr00->getName()==SCALE)
              scale=zeArr00;
            if(zeArr11->getName()==ANGLE)
              angle=zeArr11;
            if(zeArr11->getName()==SCALE)
              scale=zeArr11;
            if(angle.isNull() || scale.isNull())
              throw INTERP_KERNEL::Exception(MSG1);
            return BuildMeshPipeSEG3(angle,scale,pfl,loc,zeStr,mesh,section,globs);
          }
        default:
          throw INTERP_KERNEL::Exception(MSG1);
      }
    }
  throw INTERP_KERNEL::Exception(MSG1);
}

MCAuto<MEDFileUMesh> LocInfo::generateNonClassicalData(int zePos, const MEDFileUMesh *mesh, const MEDFileFieldGlobsReal *globs) const
{
  static const char MSG1[]="LocInfo::generateNonClassicalData : no spec for GAUSS on StructureElement with more than one cell !";
  std::size_t sz(_locs.size());
  std::vector< MCAuto<DataArrayDouble> > arrs(sz);
  for(std::size_t i=0;i<sz;i++)
    {
      const MEDFileFieldLoc& loc(globs->getLocalization(_locs[i]));
      const MEDFileGTKeeper *gtk(loc.getUndergroundGTKeeper());
      const MEDFileGTKeeperDyn *gtk2(dynamic_cast<const MEDFileGTKeeperDyn *>(gtk));
      if(!gtk2)
        throw INTERP_KERNEL::Exception("LocInfo::generateNonClassicalData : internal error !");
      const MEDFileUMesh *meshLoc(gtk2->getMesh()),*section(gtk2->getSection());
      const MEDFileStructureElement *se(gtk2->getSE());
      INTERP_KERNEL::NormalizedCellType gt;
      {
        std::vector<int> nel(meshLoc->getNonEmptyLevels());
        if(nel.size()!=1)
          throw INTERP_KERNEL::Exception(MSG1);
        if(nel[0]!=0)
          throw INTERP_KERNEL::Exception(MSG1);
        MCAuto<MEDCouplingUMesh> um(meshLoc->getMeshAtLevel(0));
        if(um->getNumberOfCells()!=1)
          throw INTERP_KERNEL::Exception(MSG1);
        gt=um->getTypeOfCell(0);
        std::vector<int> v;
        um->getNodeIdsOfCell(0,v);
        std::size_t sz2(v.size());
        for(std::size_t j=0;j<sz2;j++)
          if(v[j]!=j)
            throw INTERP_KERNEL::Exception(MSG1);
      }
      const std::vector< MCAuto<MEDFileEltStruct4Mesh> >& strs(mesh->getAccessOfUndergroundEltStrs());
      MCAuto<MEDFileEltStruct4Mesh> zeStr;
      for(std::vector< MCAuto<MEDFileEltStruct4Mesh> >::const_iterator it=strs.begin();it!=strs.end();it++)
        {
          if((*it)->getGeoTypeName()==se->getName())
            {
              zeStr=*it;
              break;
            }
        }
      if(zeStr.isNull())
        {
          std::ostringstream oss; oss << "LocInfo::generateNonClassicalData :  : no geo type with name " <<  se->getName() << " in " << mesh->getName() << " !";
          throw INTERP_KERNEL::Exception(oss.str());
        }
      arrs[i]=BuildMeshFromStructure(gt,_pfl[i],loc,zeStr,mesh,section,globs);
    }
  std::vector<const DataArrayDouble *> arrs2(VecAutoToVecOfCstPt(arrs));
  MCAuto<DataArrayDouble> resu(DataArrayDouble::Aggregate(arrs2));
  MCAuto<MEDFileUMesh> ret(MEDFileUMesh::New());
  ret->setCoords(resu);
  std::ostringstream meshName; meshName << mesh->getName() << "_on_" << sz << "_sections" << "_" << zePos;
  ret->setName(meshName.str());
  return ret;
}

FieldWalker2::FieldWalker2(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd)
{
  _loc=pmptpd->getLocalization();
  _pfl=pmptpd->getProfile();
  _is_classic=pmptpd->getType()!=ON_GAUSS_PT;
  _pd=SlicePartDefinition::New(pmptpd->getStart(),pmptpd->getEnd(),1);
}

bool FieldWalker2::operator!=(const FieldWalker2& other) const
{
  return !((*this)==other);
}

bool FieldWalker2::operator==(const FieldWalker2& other) const
{
  bool ret2(false);
  {
    std::string tmp;
    ret2=_pd->isEqual(other._pd,tmp);
  }
  bool ret(_loc==other._loc && _pfl==other._pfl && _is_classic==other._is_classic && ret2);
  return ret;
}

class FieldWalker1
{
public:
  FieldWalker1(const MEDFileAnyTypeField1TSWithoutSDA *ts):_ts(ts),_pm_pt(0),_nb_mesh(0) { }
  void newMeshEntry(const MEDFileFieldPerMesh *fpm);
  void endMeshEntry(const MEDFileFieldPerMesh *fpm) { }
  void newPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt);
  void endPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt);
  void newPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd);
  void checkOK(const FieldWalker1& other) const;
  bool isClassical() const;
  std::vector<FieldWalker2> getNonClassicalData() const { return _fw; }
private:
  const MEDFileAnyTypeField1TSWithoutSDA *_ts;
  const MEDFileFieldPerMeshPerTypeDyn *_pm_pt;
  std::vector<FieldWalker2> _fw;
  int _nb_mesh;
};

void FieldWalker1::newMeshEntry(const MEDFileFieldPerMesh *fpm)
{
  if(_nb_mesh++==1)
    throw INTERP_KERNEL::Exception("FieldWalker1::newMeshEntry : multi mesh not supported !");
}

void FieldWalker1::newPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt)
{
  if(_pm_pt)
    throw INTERP_KERNEL::Exception("FieldWalker1::newPerMeshPerTypeEntry : multi SE loc not managed yet !");
  const MEDFileFieldPerMeshPerTypeDyn *pmpt2(dynamic_cast<const MEDFileFieldPerMeshPerTypeDyn *>(pmpt));
  if(!pmpt2)
    throw INTERP_KERNEL::Exception("newPerMeshPerTypeEntry : internal error !");
  _pm_pt=pmpt2;
}

void FieldWalker1::endPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *)
{
  isClassical();
}

void FieldWalker1::newPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd)
{
  _fw.push_back(FieldWalker2(pmptpd));
}

void FieldWalker1::checkOK(const FieldWalker1& other) const
{
  std::size_t sz(_fw.size());
  if(other._fw.size()!=sz)
    throw INTERP_KERNEL::Exception("checkOK : not OK because size are not the same !");
  for(std::size_t i=0;i<sz;i++)
    if(_fw[i]!=other._fw[i])
      throw INTERP_KERNEL::Exception("checkOK : not OK because an element mismatches !");
}

bool FieldWalker1::isClassical() const
{
  if(_fw.empty())
    throw INTERP_KERNEL::Exception("FieldWalker1::endPerMeshPerTypeEntry : internal error !");
  std::size_t ic(0),inc(0);
  for(std::vector<FieldWalker2>::const_iterator it=_fw.begin();it!=_fw.end();it++)
    {
      if((*it).isClassic())
        ic++;
      else
        inc++;
    }
  if(ic!=0 && inc!=0)
    throw INTERP_KERNEL::Exception("FieldWalker1::endPerMeshPerTypeEntry : mix is not allowed yet !");
  return inc==0;
}

class FieldWalker
{
public:
  FieldWalker(const MEDFileAnyTypeFieldMultiTSWithoutSDA *f):_f(f) { }
  void newTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts);
  void endTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts);
  void newMeshEntry(const MEDFileFieldPerMesh *fpm);
  void endMeshEntry(const MEDFileFieldPerMesh *fpm);
  void newPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt);
  void endPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt);
  void newPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd);
public:
  bool isEmpty() const;
  bool isClassical() const;
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *field() const { return _f; }
  std::vector<FieldWalker2> getNonClassicalData() const { return _fw_prev->getNonClassicalData(); }
private:
  const MEDFileAnyTypeFieldMultiTSWithoutSDA *_f;
  mutable INTERP_KERNEL::AutoCppPtr<FieldWalker1> _fw;
  mutable INTERP_KERNEL::AutoCppPtr<FieldWalker1> _fw_prev;
};

bool FieldWalker::isEmpty() const
{
  return _fw_prev.isNull();
}

bool FieldWalker::isClassical() const
{
  return _fw_prev->isClassical();
}

void FieldWalker::newTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts)
{
  _fw=new FieldWalker1(ts);
}

void FieldWalker::endTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts)
{
  if(_fw_prev.isNull())
    _fw_prev=new FieldWalker1(*_fw);
  else
    _fw_prev->checkOK(*_fw);
  _fw=0;
}

void FieldWalker::newMeshEntry(const MEDFileFieldPerMesh *fpm)
{
  _fw->newMeshEntry(fpm);
}

void FieldWalker::endMeshEntry(const MEDFileFieldPerMesh *fpm)
{
  _fw->endMeshEntry(fpm);
}

void FieldWalker::newPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt)
{
  _fw->newPerMeshPerTypeEntry(pmpt);
}

void FieldWalker::endPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt)
{
  _fw->endPerMeshPerTypeEntry(pmpt);
}

void FieldWalker::newPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd)
{
  _fw->newPerMeshPerTypePerDisc(pmptpd);
}

// this class splits fields into same
class LocSpliter : public MEDFileFieldVisitor
{
public:
  LocSpliter(const MEDFileFieldGlobsReal *globs):_globs(globs),_fw(0) { }
  MCAuto<MEDFileFields> getClassical() const { return _classical; }
  void generateNonClassicalData(const MEDFileUMesh *mesh, std::vector< MCAuto<MEDFileFields> >& outFields, std::vector< MCAuto<MEDFileUMesh> >& outMeshes) const;
private:
  void newFieldEntry(const MEDFileAnyTypeFieldMultiTSWithoutSDA *field);
  void endFieldEntry(const MEDFileAnyTypeFieldMultiTSWithoutSDA *field);
  //
  void newTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts);
  void endTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts);
  //
  void newMeshEntry(const MEDFileFieldPerMesh *fpm);
  void endMeshEntry(const MEDFileFieldPerMesh *fpm);
  //
  void newPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt);
  void endPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt);
  //
  void newPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd);
private:
  const MEDFileFieldGlobsReal *_globs;
  std::vector< LocInfo > _locs;
  std::vector< MCAuto<MEDFileFields> > _fields_on_locs;//size of _locs== size of _fields_on_locs
  MCAuto<MEDFileFields> _classical;
private:
  mutable INTERP_KERNEL::AutoCppPtr<FieldWalker> _fw;
};

void LocSpliter::newFieldEntry(const MEDFileAnyTypeFieldMultiTSWithoutSDA *field)
{
  _fw=new FieldWalker(field);
}

void LocSpliter::endFieldEntry(const MEDFileAnyTypeFieldMultiTSWithoutSDA *field)
{
  if(_fw->isEmpty())
    return ;
  MCAuto<MEDFileAnyTypeFieldMultiTS> f(MEDFileAnyTypeFieldMultiTS::BuildNewInstanceFromContent(const_cast<MEDFileAnyTypeFieldMultiTSWithoutSDA *>(field)));
  if(_fw->isClassical())
    {
      if(_classical.isNull())
        {
          _classical=MEDFileFields::New();
          _classical->shallowCpyGlobs(*_globs);
        }
      _classical->pushField(f);
    }
  else
    {
      std::vector<FieldWalker2> fw2(_fw->getNonClassicalData());
      LocInfo elt(fw2);
      std::vector< LocInfo >::iterator it(std::find(_locs.begin(),_locs.end(),elt));
      if(it==_locs.end())
        {
          _locs.push_back(elt);
          MCAuto<MEDFileFields> zeF(MEDFileFields::New());
          zeF->shallowCpyGlobs(*_globs);
          zeF->pushField(f);
          _fields_on_locs.push_back(zeF);
        }
      else
        {
          MCAuto<MEDFileFields> zeF(_fields_on_locs[std::distance(_locs.begin(),it)]);
          zeF->pushField(f);
        }
    }
}

void LocSpliter::generateNonClassicalData(const MEDFileUMesh *mesh, std::vector< MCAuto<MEDFileFields> >& outFields, std::vector< MCAuto<MEDFileUMesh> >& outMeshes) const
{
  int i(0);
  for(std::vector<LocInfo>::const_iterator it=_locs.begin();it!=_locs.end();it++,i++)
    {
      MCAuto<MEDFileUMesh> m((*it).generateNonClassicalData(i,mesh,_globs));
      outMeshes.push_back(m);
      MCAuto<MEDCouplingUMesh> mcm(MEDCouplingUMesh::Build0DMeshFromCoords(m->getCoords()));
      mcm->setName(m->getName());
      MCAuto<MEDFileFields> fs(_fields_on_locs[i]);
      MCAuto<MEDFileFields> outFs(MEDFileFields::New());
      for(int j=0;j<fs->getNumberOfFields();j++)
        {
          MCAuto<MEDFileAnyTypeFieldMultiTS> fmtsNC(fs->getFieldAtPos(j));
          MCAuto<MEDFileFieldMultiTS> fmts(DynamicCastSafe<MEDFileAnyTypeFieldMultiTS,MEDFileFieldMultiTS>(fmtsNC));
          MCAuto<MEDFileFieldMultiTS> outFmts(MEDFileFieldMultiTS::New());
          for(int k=0;k<fmts->getNumberOfTS();k++)
            {
              MCAuto<MEDFileField1TS> outF1t(MEDFileField1TS::New());
              MCAuto<MEDFileField1TS> f1ts(fmts->getTimeStepAtPos(k));
              int t2,t3;
              double t1(f1ts->getTime(t2,t3));
              MCAuto<MEDCouplingFieldDouble> mcf(MEDCouplingFieldDouble::New(ON_NODES));
              MCAuto<DataArrayDouble> arr,arr2;
              arr.takeRef(f1ts->getUndergroundDataArray());
              arr2=arr->selectPartDef((*it).getPartDef());
              mcf->setArray(arr2);
              mcf->setTime(t1,t2,t3);
              mcf->setName(f1ts->getName());
              mcf->setMesh(mcm);
              outF1t->setFieldNoProfileSBT(mcf);
              outFmts->pushBackTimeStep(outF1t);
            }
          outFs->pushField(outFmts);
        }
      outFields.push_back(outFs);
    }
}

void LocSpliter::newTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts)
{
  _fw->newTimeStepEntry(ts);
}

void LocSpliter::endTimeStepEntry(const MEDFileAnyTypeField1TSWithoutSDA *ts)
{
  _fw->endTimeStepEntry(ts);
}

void LocSpliter::newMeshEntry(const MEDFileFieldPerMesh *fpm)
{
  _fw->newMeshEntry(fpm);
}

void LocSpliter::endMeshEntry(const MEDFileFieldPerMesh *fpm)
{
  _fw->endMeshEntry(fpm);
}

void LocSpliter::newPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt)
{
  _fw->newPerMeshPerTypeEntry(pmpt);
}

void LocSpliter::endPerMeshPerTypeEntry(const MEDFileFieldPerMeshPerTypeCommon *pmpt)
{
  _fw->endPerMeshPerTypeEntry(pmpt);
}

void LocSpliter::newPerMeshPerTypePerDisc(const MEDFileFieldPerMeshPerTypePerDisc *pmptpd)
{
  _fw->newPerMeshPerTypePerDisc(pmptpd);
}

void MEDFileBlowStrEltUp::DealWithConflictNames(MEDFileAnyTypeFieldMultiTS *fmtsToAdd, const MEDFileFields *fs)
{
  std::vector<std::string> fnames(fs->getFieldsNames());
  for(int i=0;i<1000;i++)
    {
      std::ostringstream oss; oss << fmtsToAdd->getName();
      if(i>=1)
        oss << "_" << i-1;
      if(std::find(fnames.begin(),fnames.end(),oss.str())==fnames.end())
        {
          fmtsToAdd->setName(oss.str());
          return ;
        }
    }
  throw INTERP_KERNEL::Exception("DealWithConflictNames : Eh eh interesting !");
}

MCAuto<MEDFileFields> MEDFileBlowStrEltUp::splitFieldsPerLoc(const MEDFileFields *fields, const MEDFileUMesh *mesh, MEDFileMeshes *msOut, MEDFileFields *allZeOutFields)
{
  LocSpliter ls(fields);
  fields->accept(ls);
  std::vector< MCAuto<MEDFileFields> > outFields;
  std::vector< MCAuto<MEDFileUMesh> > outMeshes;
  ls.generateNonClassicalData(mesh,outFields,outMeshes);
  for(std::vector< MCAuto<MEDFileFields> >::iterator it=outFields.begin();it!=outFields.end();it++)
    {
      for(int j=0;j<(*it)->getNumberOfFields();j++)
        {
          MCAuto<MEDFileAnyTypeFieldMultiTS> fmts((*it)->getFieldAtPos(j));
          //DealWithConflictNames(fmts,allZeOutFields);// uncomment to have a writable data structure
          allZeOutFields->pushField(fmts);
        }
    }
  for(std::vector< MCAuto<MEDFileUMesh> >::iterator it=outMeshes.begin();it!=outMeshes.end();it++)
    msOut->pushMesh(*it);
  return ls.getClassical();
}
