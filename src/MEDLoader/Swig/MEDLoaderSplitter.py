#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2016  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
#
# See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
#
# Author : Anthony GEAY (CEA/DEN/DM2S/STMF/LGLS)

from MEDLoader import *
import os

class MEDLoaderSplitter:
    @classmethod
    def New(cls,mfd,idsLst):
        """ mfd is a MEDFileData instance containing only one mesh. idsLst is a list of DataArrayInt containing each the ids per processor """
        return MEDLoaderSplitter(mfd,idsLst)
        pass

    def __init__(self,mfd,idsLst):
        """ mfd is a MEDFileData instance containing only one mesh. idsLst is a list of DataArrayInt containing each the ids per processor """
        mfmsh=mfd.getMeshes()
        mfflds=mfd.getFields()
        if len(mfmsh)!=1:
            raise InterpKernelException("Works only with one mesh !")
        mfflds=mfflds.partOfThisLyingOnSpecifiedMeshName(mfmsh[0].getName())
        retf=self.__splitFields(mfmsh[0],mfflds,idsLst)
        retm=self.__splitMesh(mfmsh[0],idsLst)
        self._mfd_splitted=[MEDFileData() for i in range(len(idsLst))]
        for a,b,c in zip(self._mfd_splitted,retf,retm):
            a.setFields(b) ; a.setMeshes(c)
            pass
        pass

    def getSplittedInstances(self):
        return self._mfd_splitted
    
    @classmethod
    def __splitMEDFileField1TSNode(cls,f,f1ts,ids):
        fRet=f[ids]
        f1ts.setFieldNoProfileSBT(fRet)
        pass
    
    @classmethod
    def __splitMEDFileField1TSCell(cls,f,f1ts,ids):
        fRet=f[ids]
        m=fRet.getMesh() ; m.zipCoords()
        o2n=m.getRenumArrForMEDFileFrmt() ; fRet.renumberCells(o2n,False)
        f1ts.setFieldNoProfileSBT(fRet)
        pass
    
    def __splitMEDFileField1TS(self,mm,f1ts,idsLst):
        ret=[f1ts.__class__() for i in range(len(idsLst))]
        dico={ON_CELLS:self.__splitMEDFileField1TSCell,
              ON_NODES:self.__splitMEDFileField1TSNode,
              ON_GAUSS_PT:self.__splitMEDFileField1TSCell,
              ON_GAUSS_NE:self.__splitMEDFileField1TSCell}
        for t in f1ts.getTypesOfFieldAvailable():
            f=f1ts.getFieldOnMeshAtLevel(t,0,mm)
            for i,f0 in enumerate(ret):
                dico[t](f,f0,idsLst[i])
                pass
            pass
        return ret
    
    def __splitFields(self,mm,mfflds,idsLst):
        ret0 = [MEDFileFields() for i in range(len(idsLst))]
        for fmts in mfflds:
            if len(fmts.getPflsReallyUsed())!=0:
                print("Field \"%s\" contains profiles ! Not supported yet ! This field will be ignored !" % (fmts.getName()))
                continue
            pass
            ret1=[fmts.__class__() for i in range(len(idsLst))]
            for f1ts in fmts:
                for fmtsPart,f1tsPart in zip(ret1,self.__splitMEDFileField1TS(mm,f1ts,idsLst)):
                    fmtsPart.pushBackTimeStep(f1tsPart)
                    pass
                pass
            for fieldsPart,fmtsPart in zip(ret0,ret1):
                fieldsPart.pushField(fmtsPart);
                pass
            pass
        return ret0

    def __splitMesh(self,mfm,idsLst):
        ret0 = [MEDFileMeshes() for i in range(len(idsLst))]
        m=mfm[0]
        addlevs=list(mfm.getNonEmptyLevels())[1:]
        dAddlevs={k:mfm[k] for k in addlevs}
        for ret,ids in zip(ret0,idsLst):
            mlPart=mfm.createNewEmpty()
            mPart=m[ids] ; trado2n=mPart.zipCoordsTraducer()
            trad=trado2n.invertArrayO2N2N2O(mPart.getNumberOfNodes())
            mlPart[0]=mPart
            if 0 in mfm.getFamArrNonEmptyLevelsExt():
                mlPart.setFamilyFieldArr(0,mfm.getFamilyFieldAtLevel(0)[ids])
                pass
            if 1 in mfm.getFamArrNonEmptyLevelsExt():
                mlPart.setFamilyFieldArr(1,mfm.getFamilyFieldAtLevel(1)[trad])
                pass
            for k,v in dAddlevs.iteritems():
                part=v.getCellIdsFullyIncludedInNodeIds(trad)
                mSubPart=v[part] ; mSubPart.renumberNodesInConn(trado2n) ; mSubPart.setCoords(mPart.getCoords())
                mlPart[k]=mSubPart
                mlPart.setFamilyFieldArr(k,mfm.getFamilyFieldAtLevel(k)[part])
                pass
            mlPart.copyFamGrpMapsFrom(mfm)
            ret.pushMesh(mlPart)
            pass
        return ret0
    pass
