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
      MCAuto<MEDFileUMesh> mOut;
      MCAuto<MEDFileFields> fsOut1;
      MCAuto<MEDFileEltStruct4Mesh> zeStr(dealWithSEInMesh(ps[0].second,umesh,mOut,fsOut1));
      msOut->pushMesh(mOut);
      dealWithSEInFields(ps[0].second,*elt,zeStr,fsOut1,allZeOutFields);
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
}
