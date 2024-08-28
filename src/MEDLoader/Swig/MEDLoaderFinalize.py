#  -*- coding: iso-8859-1 -*-
# Copyright (C) 2023-2024  CEA, EDF
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

import logging

def MEDFileUMeshFuseNodesAndCells(self, compType = 2 , eps = 1e-6, logLev = logging.INFO):
  """
  [EDF30179] : Method fusing nodes in this, then fusing cells in this. Fusion is done following eps and compType.

  :param compType : see MEDCouplingPointSet.zipConnectivityTraducer method for explanations
  :param eps: see DataArrayDouble.findCommonTuples for explanations.
  :param logLev: Integer specifying log level
  :return: MEDFileUMesh instance containing the result of nodes and cells fusion
  """
  import MEDLoader as ml
  def getLogger( level = logging.INFO ):
    FORMAT = '%(levelname)s : %(asctime)s : [%(filename)s:%(funcName)s:%(lineno)s] : %(message)s'
    logging.basicConfig( format = FORMAT, level = level )
    return logging.getLogger()
  logger = getLogger( logLev )
  def updateMap( mm : ml.MEDFileUMesh, lev : int, famMap : ml.DataArrayInt, famMapI : ml.DataArrayInt):
    """
    mm instance to be updated
    """
    def famIdManager(lev, famId):
      if lev<= 0:
        return -famId
      else:
        return famId
    nbOfPartSetToBeUpdated = len(famMapI) -1
    for partSetId in range( nbOfPartSetToBeUpdated ):
      newFamId = famIdManager( lev, famMap[ famMapI[partSetId] ] )
      newFamName = f"Family_{newFamId}"
      logger.debug(f"For level {lev} new family : {newFamId}")
      mm.addFamily(newFamName,newFamId)
      for famId in famMap[ famMapI[partSetId]+1 :famMapI[partSetId+1] ]:
        grpsToBeUpdated = mm.getGroupsOnFamily( mm.getFamilyNameGivenId( famIdManager( lev, int(famId) ) ) )
        for grpToBeUpdated in grpsToBeUpdated:
          mm.addFamilyOnGrp( grpToBeUpdated, newFamName )
    pass
  getLogger( level = logLev )
  initNbNodes = len( self.getCoords() )
  logger.info(f"Begin merging nodes with eps = {eps}")
  cc,cci = self.getCoords().findCommonTuples( eps )
  logger.info(f"End of merging nodes with eps = {eps} : Nb of nodes groups to be merged : {len(cci)-1} / {self.getNumberOfNodes()}")
  o2n,newNbNodes = ml.DataArrayInt.ConvertIndexArrayToO2N(initNbNodes,cc,cci)
  n2oNodes = o2n.invertArrayO2N2N2O( newNbNodes )
  newCoords = self.getCoords()[n2oNodes]
  # creation of 
  mmOut = ml.MEDFileUMesh()
  mmOut.copyFamGrpMapsFrom( self )

  for lev in self.getNonEmptyLevels():
    logger.debug(f"Begin level {lev}")
    m1 = self[lev].deepCopy()
    logger.debug(f"Begin renumbering connectivity of level {lev}")
    m1.renumberNodesInConn( o2n )
    logger.debug(f"End renumbering connectivity of level {lev}")
    m1.setCoords( newCoords )
    logger.info(f"Begin of finding of same cells of level {lev}")
    cce,ccei = m1.findCommonCells(compType,0)
    logger.info(f"End of finding of same cells of level {lev} : Nb of cells groups to be merged : {len(ccei)-1} / {m1.getNumberOfCells()}")
    famsCell = self.getFamilyFieldAtLevel(lev)
    if famsCell:
      famsCell = -famsCell
      famsMergedCell,famMap,famMapI = famsCell.forThisAsPartitionBuildReduction(cce,ccei) # <- method updating family field array
      updateMap(mmOut,lev,famMap,famMapI)
      famsMergedCell = -famsMergedCell
    o2nCells,newNbCells = ml.DataArrayInt.ConvertIndexArrayToO2N(m1.getNumberOfCells(),cce,ccei)
    n2oCells = o2nCells.invertArrayO2N2N2O( newNbCells )
    m1 = m1[ n2oCells ]
    m1.setCoords( newCoords )
    m1.setName( self.getName() )
    mmOut[lev] = m1
    if famsCell:
      mmOut.setFamilyFieldArr( lev, famsMergedCell )

  famsNode = self.getFamilyFieldAtLevel(1)
  if famsNode:
    famsMergedNode,famMap,famMapI = famsNode.forThisAsPartitionBuildReduction(cc,cci)
    updateMap(mmOut,1,famMap,famMapI)
    mmOut.setFamilyFieldArr(1, famsMergedNode)
  return mmOut

def MEDFileUMeshReduceToCells(self, level, keepCells, removeOrphanNodes=True):
    """
    Method returning a new MEDFileUMesh, restriction of self to level and keepCell cells at this level.
    This method also 

    :param level: Specifies the top level of the returned MEDFileUMesh expected
    :param keepCells: A DataArrayInt specifying cell ids at level level of self
    :param removeOrphanNodes: Specifies if orphan nodes should be removed at the end
    
    see also MEDFileUMesh.extractPart
    """
    import MEDLoader as ml
    subLevs = [l for l in self.getNonEmptyLevels() if l<=level]
    subMeshes = [self[lev] for lev in subLevs]
    allFamilyFields = [self.getFamilyFieldAtLevel(lev) for lev in subLevs]
    allRefMesh = subMeshes[0]
    refMesh = allRefMesh[keepCells]

    mmOut = ml.MEDFileUMesh()
    # level 0
    mmOut[0] = refMesh
    mmOut.setFamilyFieldArr(0,allFamilyFields[0][keepCells])

    # subLevels
    for curLev,meshLev,famFieldLev in zip(subLevs[1:],subMeshes[1:],allFamilyFields[1:]):
        allMeshLev,d,di, rd,rdi = allRefMesh.explodeMeshTo( curLev-level )
        a,b = allMeshLev.areCellsIncludedIn(meshLev,2)
        if not a:
            raise RuntimeError("Error in mesh {}")
        dlev,dlevi = ml.DataArrayInt.ExtractFromIndexedArrays( keepCells, d,di )
        dlev2 = dlev.buildUniqueNotSorted()
        cellsToKeepLev = ml.DataArrayInt.BuildIntersection([dlev2,b])
        cellsToKeepLev = b.indicesOfSubPart(cellsToKeepLev)
        cellsToKeepLev.sort()
        mmOut[curLev] = meshLev[cellsToKeepLev]
        mmOut.setFamilyFieldArr(curLev,famFieldLev[cellsToKeepLev])

    allFamNodes = mmOut.getFamilyFieldAtLevel(1)
    if allFamNodes:
        mmOut.setFamilyFieldArr(1,allFamNodes[:])

    if removeOrphanNodes:
        mmOut.zipCoords()

    mmOut.copyFamGrpMapsFrom(self)
    return mmOut
