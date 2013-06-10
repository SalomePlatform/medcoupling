#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2007-2013  CEA/DEN, EDF R&D
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License.
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
# Author Anthony GEAY (CEA/DEN/DM2S/STMF/LGLS)

import numpy as np
from MEDLoader import *
from CaseIO import CaseIO
import sys,re

class CaseReader(CaseIO):
    """ Converting a file in the Case format (Ensight) to the MED format.
    A new file with the same base name and the .med extension is created.
    """

    @classmethod
    def New(cls,fileName):
        """ Static constructor. """
        return CaseReader(fileName)
        pass

    def __init__(self,fileName):
        """ Constructor """
        self._fileName=fileName
        self._dirName=os.path.dirname(self._fileName)
        pass

    def __traduceMesh(self,name,typ,coords,cells):
        """ Convert a CASE mesh into a MEDCouplingUMesh. """
        nbCoords=len(coords)
        coo=np.array(coords,dtype="float64") ; coo=coo.reshape(nbCoords,3)
        coo=DataArrayDouble(coo) ; coo=coo.fromNoInterlace()
        ct=self.dictMCTyp2[typ]
        m=MEDCouplingUMesh(name,MEDCouplingUMesh.GetDimensionOfGeometricType(ct))
        m.setCoords(coo)
        nbNodesPerCell=MEDCouplingMesh.GetNumberOfNodesOfGeometricType(ct)
        cI=DataArrayInt(len(cells)+1) ; cI.iota() ; cI*=nbNodesPerCell+1
        #
        cells2=cells.reshape(len(cells),nbNodesPerCell)
        c2=DataArrayInt(cells2)
        c=DataArrayInt(len(cells),nbNodesPerCell+1) ; c[:,0]=ct ; c[:,1:]=c2-1 ; c.rearrange(1)
        m.setConnectivity(c,cI,True)
        m.checkCoherency2()
        return m

    def __traduceMeshForPolyhed(self,name,coords,arr0,arr1,arr2):
        nbCoords=len(coords)
        coo=np.array(coords,dtype="float64") ; coo=coo.reshape(nbCoords,3)
        coo=DataArrayDouble(coo) ; coo=coo.fromNoInterlace()
        m=MEDCouplingUMesh(name,3)
        m.setCoords(coo)
        #
        arr2=arr2[:]-1
        arr0mc0=DataArrayInt(arr0) ; arr0mc0.computeOffsets2()
        arr0mc1=DataArrayInt(arr0).deepCpy()
        arr0mc2=DataArrayInt(len(arr0),2) ; arr0mc2[:,0]=DataArrayInt(arr0)-1 ; arr0mc2[:,1]=1 ; arr0mc2.rearrange(1) ; arr0mc2.computeOffsets2()
        arr0mc3=DataArrayInt.Range(0,2*len(arr0),2).buildExplicitArrByRanges(arr0mc2)
        arr1mc0=DataArrayInt(arr1) ; arr1mc0.computeOffsets2()
        arr1mc1=arr1mc0[arr0mc0] ; arr1mc1[1:]+=arr0mc0[1:] 
        arr1mc2=DataArrayInt(arr1).deepCpy() ; arr1mc2+=1 ; arr1mc2.computeOffsets2()
        arr2mc0=(arr1mc2[1:])[arr0mc3]
        #
        c=DataArrayInt(arr1.size+arr2.size)
        c[arr1mc1[:-1]]=NORM_POLYHED
        c[arr2mc0]=-1
        a=arr2mc0.buildUnion(arr1mc1[:-1]).buildComplement(len(c))
        c[a]=DataArrayInt(arr2)
        #
        m.setConnectivity(c,arr1mc1,True)
        m.checkCoherency2()
        return m

    def __traduceMeshForPolygon(self,name,coords,arr0,arr1):
        nbCoords=len(coords)
        coo=np.array(coords,dtype="float64") ; coo=coo.reshape(nbCoords,3)
        coo=DataArrayDouble(coo) ; coo=coo.fromNoInterlace()
        m=MEDCouplingUMesh(name,2)
        m.setCoords(coo)
        #
        arr0_0=DataArrayInt(arr0+1) ; arr0_0.computeOffsets2()
        arr0_1=DataArrayInt(len(arr0),2) ; arr0_1[:,1]=DataArrayInt(arr0) ; arr0_1[:,0]=1 ; arr0_1.rearrange(1) ; arr0_1.computeOffsets2()
        arr0_2=DataArrayInt.Range(1,2*len(arr0),2).buildExplicitArrByRanges(arr0_1)
        c=DataArrayInt(len(arr0)+len(arr1)) ; c[:]=0 ; c[arr0_0[:-1]]=NORM_POLYGON
        c[arr0_2]=DataArrayInt(arr1-1)
        #
        m.setConnectivity(c,arr0_0,True)
        m.checkCoherency2()
        return m

    def __convertGeo2MED(self,geoFileName):
        """ Convert all the geometry (all the meshes) contained in teh CASE file into MEDCouplingUMesh'es. """
        fd=open(os.path.join(self._dirName,geoFileName),"r+b") ; fd.seek(0,2) ; end=fd.tell() ; fd.seek(0) ; fd.readline() ; fd.readline()
        name=fd.readline().strip() ; fd.readline() ; fd.readline()
        pos=fd.tell()
        mcmeshes=[]
        elt=fd.read(80) ; elt=elt.strip() ; pos+=80
        while pos!=end:
            if elt!="part":
                raise Exception("Error on reading mesh #1 !")
            fd.seek(fd.tell()+4)
            meshName=fd.read(80).strip()
            if fd.read(len("coordinates"))!="coordinates":
                raise Exception("Error on reading mesh #2 !")
            pos=fd.tell()
            typeOfCoo=np.memmap(fd,dtype='byte',mode='r',offset=int(pos),shape=(1)).tolist()[0]
            pos+=1+17*4
            nbNodes=np.memmap(fd,dtype='int32',mode='r',offset=int(pos),shape=(1,)).tolist()[0]
            pos+=4
            coo=np.memmap(fd,dtype='float32',mode='r',offset=int(pos),shape=(nbNodes,3))
            pos+=nbNodes*3*4 ; fd.seek(pos)#np.array(0,dtype='float%i'%(typeOfCoo)).nbytes
            typ=fd.read(80).strip() ; pos=fd.tell()
            mcmeshes2=[]
            while pos!=end and typ!="part":
                mctyp=self.dictMCTyp2[typ]
                nbCellsOfType=np.memmap(fd,dtype='int32',mode='r',offset=int(pos),shape=(1,)).tolist()[0]
                pos+=4
                if mctyp!=NORM_POLYHED and mctyp!=NORM_POLYGON:
                    nbNodesPerCell=MEDCouplingMesh.GetNumberOfNodesOfGeometricType(mctyp)
                    cells=np.memmap(fd,dtype='int32',mode='r',offset=pos,shape=(nbCellsOfType,nbNodesPerCell))
                    pos+=nbCellsOfType*nbNodesPerCell*4
                    fd.seek(pos)
                    mcmeshes2.append(self.__traduceMesh(meshName,typ,coo,cells))
                elif mctyp==NORM_POLYHED:
                    nbOfFacesPerCell=np.memmap(fd,dtype='int32',mode='r',offset=int(pos),shape=(nbCellsOfType,))
                    pos+=nbCellsOfType*4
                    szOfNbOfNodesPerFacePerCellArr=int(nbOfFacesPerCell.sum())
                    arr1=np.memmap(fd,dtype='int32',mode='r',offset=int(pos),shape=(szOfNbOfNodesPerFacePerCellArr,))#arr1 -> nbOfNodesPerFacePerCellArr
                    pos+=szOfNbOfNodesPerFacePerCellArr*4
                    szOfNodesPerFacePerCellArr=arr1.sum()
                    arr2=np.memmap(fd,dtype='int32',mode='r',offset=int(pos),shape=(szOfNodesPerFacePerCellArr,))#arr2 -> nodesPerFacePerCellArr
                    pos+=szOfNodesPerFacePerCellArr*4 ; fd.seek(pos)
                    mcmeshes2.append(self.__traduceMeshForPolyhed(meshName,coo,nbOfFacesPerCell,arr1,arr2))
                    pass
                else:
                    nbOfNodesPerCell=np.memmap(fd,dtype='int32',mode='r',offset=int(pos),shape=(nbCellsOfType,))
                    pos+=nbCellsOfType*4
                    szOfNbOfNodesPerCellArr=int(nbOfNodesPerCell.sum())
                    arr1=np.memmap(fd,dtype='int32',mode='r',offset=int(pos),shape=(szOfNbOfNodesPerCellArr,))
                    pos+=szOfNbOfNodesPerCellArr*4  ; fd.seek(pos)
                    mcmeshes2.append(self.__traduceMeshForPolygon(meshName,coo,nbOfNodesPerCell,arr1))
                if pos!=end:
                    elt=fd.read(80) ; elt=elt.strip() ; typ=elt[:] ; pos+=80
                    pass
                pass
            coo=mcmeshes2[0].getCoords() ; name=mcmeshes2[0].getName()
            for itmesh in mcmeshes2: itmesh.setCoords(coo)
            m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(mcmeshes2) ; m.setName(name)
            mcmeshes.append(m)
            pass
        #
        ms=MEDFileMeshes()
        ms.resize(len(mcmeshes))
        for i,m in enumerate(mcmeshes):
            mlm=MEDFileUMesh()
            mlm.setMeshAtLevel(0,m)
            ms.setMeshAtPos(i,mlm)
            pass
        return mcmeshes,ms
    
    def __convertField(self,mlfields, mcmeshes, fileName, fieldName, discr, nbCompo, locId, it):
        """ Convert the fields. """
        stars=re.search("[\*]+",fileName).group()
        st="%0"+str(len(stars))+"i"
        trueFileName=fileName.replace(stars,st%(it))
        fd=open(os.path.join(self._dirName,trueFileName),"r+b") ; fd.seek(0,2) ; end=fd.tell() ; fd.seek(0)
        name=fd.readline().strip().split(" ")[0]
        if name!=fieldName:
            raise Exception("ConvertField : mismatch")
        pos=fd.tell()
        st=fd.read(80) ; st=st.strip() ; pos=fd.tell()
        while pos!=end:
            if st!="part":
                raise Exception("ConvertField : mismatch #2")
            fdisc=MEDCouplingFieldDiscretization.New(self.discSpatial2[discr])
            meshId=np.memmap(fd,dtype='int32',mode='r',offset=int(pos),shape=(1)).tolist()[0]-1
            nbOfValues=fdisc.getNumberOfTuples(mcmeshes[meshId])
            vals2=DataArrayDouble(nbOfValues,nbCompo)
            fd.seek(pos+4)
            st=fd.read(80).strip() ; pos=fd.tell()
            offset=0
            while pos!=end and st!="part":
                if st!="coordinates":
                    nbOfValsOfTyp=mcmeshes[meshId].getNumberOfCellsWithType(self.dictMCTyp2[st])
                else:
                    nbOfValsOfTyp=nbOfValues
                    pass
                vals=np.memmap(fd,dtype='float32',mode='r',offset=int(pos),shape=(nbOfValsOfTyp,nbCompo))#np.memmap(fd,dtype='int32',mode='r',offset=159,shape=(1))
                vals2[offset:offset+nbOfValsOfTyp]=DataArrayDouble(np.array(vals,dtype='float64')).fromNoInterlace()
                pos+=nbOfValsOfTyp*nbCompo*4 ; fd.seek(pos)
                st=fd.read(80) ; st=st.strip() ; pos=fd.tell()
                offset+=nbOfValsOfTyp
                pass
            f=MEDCouplingFieldDouble(self.discSpatial2[discr],ONE_TIME) ; f.setName("%s_%s"%(fieldName,mcmeshes[meshId].getName()))
            f.setMesh(mcmeshes[meshId]) ; f.setArray(vals2) ; f.setTime(float(it),it,-1)
            f.checkCoherency()
            mlfields[locId+meshId].appendFieldNoProfileSBT(f)
            pass
        pass
    
    def loadInMEDFileDS(self):
        """ Load a CASE file into a MEDFileData object. """
        f=file(self._fileName)
        lines=f.readlines()
        ind=lines.index("GEOMETRY\n")
        if ind==-1:
            raise Exception("Error with file %s"%(fname))
        geoName=re.match("model:([\W]*)([\w\.]+)",lines[ind+1]).group(2)
        m1,m2=self.__convertGeo2MED(geoName)
        ind=lines.index("VARIABLE\n")
        fieldsInfo=[]
        for i in xrange(ind+1,lines.index("TIME\n")):
            m=re.match("^([\w]+)[\s]+\per[\s]+([\w]+)[\s]*\:[\s]*([\w]+)[\s]+([\S]+)$",lines[i])
            if m:
                spatialDisc=m.groups()[1] ; fieldName=m.groups()[2] ; nbOfCompo=self.dictCompo2[m.groups()[0]] ; fieldFileName=m.groups()[3]
                fieldsInfo.append((fieldName,spatialDisc,nbOfCompo,fieldFileName))
                pass
            pass
        
        expr=re.compile("number[\s]+of[\s]+steps[\s]*\:[\s]*([\d]+)")
        nbOfTimeSteps=int(expr.search(filter(expr.search,lines)[0]).group(1))
        
        expr=re.compile("filename[\s]+start[\s]+number[\s]*\:[\s]*([\d]+)")
        startIt=int(expr.search(filter(expr.search,lines)[0]).group(1))
        
        expr=re.compile("filename[\s]+increment[\s]*\:[\s]*([\d]+)")
        incrIt=int(expr.search(filter(expr.search,lines)[0]).group(1))
        
        curIt=startIt
        mlfields=MEDFileFields()
        mlfields.resize(len(fieldsInfo)*len(m1))
        i=0
        for field in fieldsInfo:
            for m in m1:
                mlfields.setFieldAtPos(i,MEDFileFieldMultiTS())
                i+=1
                pass
            pass
        for ts in xrange(nbOfTimeSteps):
            i=0
            for field in fieldsInfo:
                self.__convertField(mlfields,m1,field[3],field[0],field[1],field[2],i,curIt)
                i+=len(m1)
                pass
            curIt+=incrIt
            pass
        ret=MEDFileData()
        ret.setMeshes(m2)
        del mlfields[filter(lambda x: len(mlfields[x])==0,range(len(mlfields)))]
        ret.setFields(mlfields)
        return ret

    pass
