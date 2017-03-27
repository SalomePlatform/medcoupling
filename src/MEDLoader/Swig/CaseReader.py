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

# http://www-vis.lbl.gov/NERSC/Software/ensight/doc/OnlineHelp/UM-C11.pdf
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
        if cells2.dtype=='int32':
            c2=DataArrayInt(cells2)
        else:
            c2=DataArrayInt(np.array(cells2,dtype="int32"))
            pass
        c=DataArrayInt(len(cells),nbNodesPerCell+1) ; c[:,0]=ct ; c[:,1:]=c2-1 ; c.rearrange(1)
        m.setConnectivity(c,cI,True)
        m.checkConsistency()
        return m

    def __traduceMeshForPolyhed(self,name,coords,arr0,arr1,arr2):
        nbCoords=len(coords)
        coo=np.array(coords,dtype="float64") ; coo=coo.reshape(nbCoords,3)
        coo=DataArrayDouble(coo) ; coo=coo.fromNoInterlace()
        m=MEDCouplingUMesh(name,3)
        m.setCoords(coo)
        #
        arr2=arr2[:]-1
        arr0mc0=DataArrayInt(arr0) ; arr0mc0.computeOffsetsFull()
        arr0mc1=DataArrayInt(arr0).deepCopy()
        arr0mc2=DataArrayInt(len(arr0),2) ; arr0mc2[:,0]=DataArrayInt(arr0)-1 ; arr0mc2[:,1]=1 ; arr0mc2.rearrange(1) ; arr0mc2.computeOffsetsFull()
        arr0mc3=DataArrayInt.Range(0,2*len(arr0),2).buildExplicitArrByRanges(arr0mc2)
        arr1mc0=DataArrayInt(arr1) ; arr1mc0.computeOffsetsFull()
        arr1mc1=arr1mc0[arr0mc0] ; arr1mc1[1:]+=arr0mc0[1:] 
        arr1mc2=DataArrayInt(arr1).deepCopy() ; arr1mc2+=1 ; arr1mc2.computeOffsetsFull()
        arr2mc0=(arr1mc2[1:])[arr0mc3]
        #
        c=DataArrayInt(arr1.size+arr2.size)
        c[arr1mc1[:-1]]=NORM_POLYHED
        c[arr2mc0]=-1
        a=arr2mc0.buildUnion(arr1mc1[:-1]).buildComplement(len(c))
        c[a]=DataArrayInt(arr2)
        #
        m.setConnectivity(c,arr1mc1,True)
        m.checkConsistency()
        return m

    def __traduceMeshForPolygon(self,name,coords,arr0,arr1):
        nbCoords=len(coords)
        coo=np.array(coords,dtype="float64") ; coo=coo.reshape(nbCoords,3)
        coo=DataArrayDouble(coo) ; coo=coo.fromNoInterlace()
        m=MEDCouplingUMesh(name,2)
        m.setCoords(coo)
        #
        arr0_0=DataArrayInt(arr0+1) ; arr0_0.computeOffsetsFull()
        arr0_1=DataArrayInt(len(arr0),2) ; arr0_1[:,1]=DataArrayInt(arr0) ; arr0_1[:,0]=1 ; arr0_1.rearrange(1) ; arr0_1.computeOffsetsFull()
        arr0_2=DataArrayInt.Range(1,2*len(arr0),2).buildExplicitArrByRanges(arr0_1)
        c=DataArrayInt(len(arr0)+len(arr1)) ; c[:]=0 ; c[arr0_0[:-1]]=NORM_POLYGON
        c[arr0_2]=DataArrayInt(arr1-1)
        #
        m.setConnectivity(c,arr0_0,True)
        m.checkConsistency()
        return m

    def __convertGeo2MED(self,geoFileName):
        """ Convert all the geometry (all the meshes) contained in teh CASE file into MEDCouplingUMesh'es. """
        fd=open(os.path.join(self._dirName,geoFileName),"r+b") ; fd.seek(0,2) ; end=fd.tell() ; fd.seek(0) ; title=fd.read(80)
        title=title.strip().lower()
        if "binary" not in title:
            raise Exception("Error only binary geo files are supported for the moment !")
            pass
        zeType=True
        if "fortran" in title:
            mcmeshes=self.__convertGeo2MEDFortran(fd,end) ; zeType=False
        else:
            mcmeshes=self.__convertGeo2MEDC(fd,end)
        #
        ms=MEDFileMeshes()
        ms.resize(len(mcmeshes))
        for i,m in enumerate(mcmeshes):
            mlm=MEDFileUMesh()
            mlm.setMeshAtLevel(0,m)
            ms.setMeshAtPos(i,mlm)
            pass
        return mcmeshes,ms,zeType

    def __convertGeo2MEDFortran(self,fd,end):
        mcmeshes=[]
        fd.read(80) # comment 1
        fd.read(80) # comment 2
        fd.read(80) # node id
        fd.read(80) # element id
        pos=fd.tell()
        elt=fd.read(80) ; elt=elt.strip() ; pos=fd.tell()
        mcmeshes2=[]
        typ="part"
        nbOfTurn=0
        while abs(pos-end)>8 and "part" in typ:
            if "part" not in elt:
                raise Exception("Error on reading mesh fortran #1 !")
            fd.seek(fd.tell()+4)# skip #
            tmp=fd.read(80) ; meshName=tmp.split("P")[-1]
            tmp=fd.read(80)
            if "coordinates" not in tmp:
                raise Exception("Error on reading mesh fortran #2 !")
            pos=fd.tell() # 644
            if nbOfTurn==0:
                pos+=76 # what else ?
            else:
                pos+=40
                pass
            nbNodes=np.memmap(fd,dtype='>i4',mode='r',offset=int(pos),shape=(1,)).tolist()[0]
            pos+=12 # what else ?
            a=np.memmap(fd,dtype='>f4',mode='r',offset=int(pos),shape=(nbNodes))
            b=np.memmap(fd,dtype='>f4',mode='r',offset=int(pos+nbNodes*4+2*4),shape=(nbNodes))
            c=np.memmap(fd,dtype='>f4',mode='r',offset=int(pos+nbNodes*2*4+4*4),shape=(nbNodes))
            coo=np.zeros(dtype=">f4",shape=(nbNodes*3))
            coo[:nbNodes]=a ; coo[nbNodes:2*nbNodes]=b ; coo[2*nbNodes:]=c
            coo=coo.reshape(nbNodes,3)
            pos+=nbNodes*3*4 ; fd.seek(pos)#np.array(0,dtype='float%i'%(typeOfCoo)).nbytes
            typ=fd.read(80).strip() ; pos=fd.tell()
            zeK=""
            for k in self.dictMCTyp2:
                if k in typ:
                    zeK=k
                    break
                    pass
                pass
            pos+=8*4 # yeh man !
            nbCellsOfType=np.memmap(fd,dtype='>i4',mode='r',offset=int(pos),shape=(1,)).tolist()[0]
            pos+=4 # for the number of cells
            pos+=2*4 # because it's great !
            nbNodesPerCell=MEDCouplingMesh.GetNumberOfNodesOfGeometricType(self.dictMCTyp2[zeK])
            nodalConn=np.memmap(fd,dtype='>i4',mode='r',offset=pos,shape=(nbCellsOfType,nbNodesPerCell))
            meshName=meshName.strip()
            mcmeshes2.append(self.__traduceMesh(meshName,zeK,coo,nodalConn))
            pos+=nbNodesPerCell*nbCellsOfType*4
            if abs(pos-end)>8:
                fd.seek(pos) ;elt=fd.read(80) ; typ=elt[:] ; pos+=80 
                pass
            nbOfTurn+=1
            pass
        #coo=mcmeshes2[0].getCoords() ; name=mcmeshes2[0].getName()
        #for itmesh in mcmeshes2: itmesh.setCoords(coo)
        #m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(mcmeshes2) ; m.setName(name)
        #mcmeshes.append(m)
        return mcmeshes2

    def __convertGeo2MEDC(self,fd,end):
        fd.readline()
        name=fd.readline().strip() ; fd.readline() ; fd.readline()
        pos=fd.tell()
        mcmeshes=[]
        elt=fd.read(80) ; elt=elt.strip() ; pos+=80
        while pos!=end:
            if "part" not in elt:
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
                if typ[0]=='\0': pos+=1; continue
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
            if mcmeshes2:
                coo=mcmeshes2[0].getCoords() ; name=mcmeshes2[0].getName()
                for itmesh in mcmeshes2: itmesh.setCoords(coo)
                m=MEDCouplingUMesh.MergeUMeshesOnSameCoords(mcmeshes2) ; m.setName(name)
                mcmeshes.append(m)
            pass
        return mcmeshes
        
    
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
            if meshId >= len( mcmeshes ):
                return
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
            f.checkConsistencyLight()
            mlfields[locId+meshId].appendFieldNoProfileSBT(f)
            pass

    def __convertFieldFortran(self,mlfields, mcmeshes, fileName, fieldName, discr, nbCompo, locId, it):
        """ Convert the fields. """
        if re.search("[\*]+",fileName):
            stars=re.search("[\*]+",fileName).group()
            st="%0"+str(len(stars))+"i"
            trueFileName=fileName.replace(stars,st%(it))
            pass
        else:
            trueFileName=fileName
            pass
        fd=open(os.path.join(self._dirName,trueFileName),"r+b") ; fd.seek(0,2) ; end=fd.tell() ; fd.seek(0)
        name=fd.read(80)
        if fieldName not in name:
            raise Exception("ConvertField : mismatch")
        pos=fd.tell()
        st=fd.read(80) ; st=st.strip() ; pos=fd.tell()
        if "part" not in st:
            raise Exception("ConvertField : mismatch #2")
        st=fd.read(80).strip() ; pos=fd.tell()
        pos+=12 # I love it
        offset=0
        nbTurn=0
        while pos!=end and "part" not in st:
            fdisc=MEDCouplingFieldDiscretization.New(self.discSpatial2[discr])
            nbOfValues=fdisc.getNumberOfTuples(mcmeshes[nbTurn])
            vals2=DataArrayDouble(nbOfValues,nbCompo)
            pos+=24 # I love it again !
            nbOfValsOfTyp=np.memmap(fd,dtype='>i4',mode='r',offset=pos,shape=(1)).tolist()[0]/4
            pos+=4
            vals=np.zeros(dtype=">f4",shape=(nbOfValsOfTyp*nbCompo))
            for iii in range(nbCompo):
                valsTmp=np.memmap(fd,dtype='>f4',mode='r',offset=int(pos),shape=(nbOfValsOfTyp))
                vals[iii*nbOfValsOfTyp:(iii+1)*nbOfValsOfTyp]=valsTmp
                pos+=nbOfValsOfTyp*4
                pos+=2*4 ## hey hey, that is the ultimate class !
                vals2.setInfoOnComponent(iii,chr(ord('X')+iii))
                pass
            if pos>end:
                pos=end
                pass
            vals=vals.reshape(nbOfValsOfTyp,nbCompo)
            vals2[offset:offset+nbOfValsOfTyp]=DataArrayDouble(np.array(vals,dtype='float64')).fromNoInterlace()
            if pos!=end:
                fd.seek(pos)
                st=fd.read(80) ; st=st.strip() ; pos=fd.tell()
                st=fd.read(80) ; st=st.strip() ; pos=fd.tell()
                pass
            f=MEDCouplingFieldDouble(self.discSpatial2[discr],ONE_TIME) ; f.setName("%s_%s"%(fieldName,mcmeshes[nbTurn].getName()))
            f.setMesh(mcmeshes[nbTurn]) ; f.setArray(vals2) ; f.setTime(float(it),it,-1)
            f.checkConsistencyLight()
            mlfields[locId+nbTurn].appendFieldNoProfileSBT(f)
            nbTurn+=1
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
        m1,m2,typeOfFile=self.__convertGeo2MED(geoName)
        fieldsInfo=[] ; nbOfTimeSteps=0
        if "VARIABLE\n" in lines:
            ind=lines.index("VARIABLE\n")
            end=len(lines)-1
            if "TIME\n" in lines:
                end=lines.index("TIME\n")
                pass
            for i in range(ind + 1, end):
                m=re.match("^([\w]+)[\s]+\per[\s]+([\w]+)[\s]*\:[\s]*([\w]+)[\s]+([\S]+)$",lines[i])
                if m:
                    if m.groups()[0]=="constant":
                        continue
                    spatialDisc=m.groups()[1] ; fieldName=m.groups()[2] ; nbOfCompo=self.dictCompo2[m.groups()[0]] ; fieldFileName=m.groups()[3]
                    fieldsInfo.append((fieldName,spatialDisc,nbOfCompo,fieldFileName))
                    pass
                pass
            
            expr=re.compile("number[\s]+of[\s]+steps[\s]*\:[\s]*([\d]+)")
            tmp = [line for line in lines if expr.search(line)]
            if tmp:
                nbOfTimeSteps = int(expr.search(tmp[0]).group(1))
                expr=re.compile("filename[\s]+start[\s]+number[\s]*\:[\s]*([\d]+)")
                startIt = int(expr.search([line for line in lines if expr.search(line)][0]).group(1))
                expr=re.compile("filename[\s]+increment[\s]*\:[\s]*([\d]+)")
                incrIt = int(expr.search([line for line in lines if expr.search(line)][0]).group(1))
            else:
                nbOfTimeSteps=1
                startIt=0
                incrIt=1
                pass
            curIt=startIt
            pass
        mlfields=MEDFileFields()
        mlfields.resize(len(fieldsInfo)*len(m1))
        i=0
        for field in fieldsInfo:
            for m in m1:
                mlfields.setFieldAtPos(i,MEDFileFieldMultiTS())
                i+=1
                pass
            pass
        for ts in range(nbOfTimeSteps):
            i=0
            for field in fieldsInfo:
                if typeOfFile:
                    self.__convertField(mlfields,m1,field[3],field[0],field[1],field[2],i,curIt);
                else:
                    self.__convertFieldFortran(mlfields,m1,field[3],field[0],field[1],field[2],i,curIt)
                    pass
                i+=len(m1)
                pass
            curIt+=incrIt
            pass
        ret=MEDFileData()
        ret.setMeshes(m2)
        del mlfields[[x for x in range(len(mlfields)) if len(mlfields[x]) == 0]]
        ret.setFields(mlfields)
        return ret

    pass
