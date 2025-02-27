""" Python version of MEDCouplingUMesh::findNodesToDuplicate() and MEDCouplingUMesh::findCellsToRenumber() methods which are at the core of the
    MEDFileUMesh::buildInnerBoundaryAlongM1Group() algorithm.
    This greatly helps algorithm tuning ...
"""

from medcoupling import *

def findNodesToDuplicate(this, otherDimM1OnSameCoords):
  # Checking star-shaped M1 group:
  meshM2, _,_,_,rdit0 = otherDimM1OnSameCoords.buildDescendingConnectivity() # 2D: a mesh of points, 3D: a mesh of segs
  dsi = rdit0.deltaShiftIndex()
  idsTmp0 = dsi.findIdsNotInRange(-1, 3)  # for 2D: if a point is connected to more than 2 segs. For 3D: if a seg is connected to more than two faces.
  if(idsTmp0.getNumberOfTuples()):
    raise ValueError("")

  # Get extreme nodes from the group (they won't be duplicated except if they also lie on bound of M0 -- see below),
  # ie nodes belonging to the boundary "cells" (might be points) of M1
  xtremIdsM2 = dsi.findIdsEqual(1)
  meshM2Part = meshM2[xtremIdsM2]
  xtrem = meshM2Part.computeFetchedNodeIds()
  # Remove from the list points on the boundary of the M0 mesh (those need duplication!)
  m0desc, dt0, dit0, rdt0, rdit0 = this.buildDescendingConnectivity()
  dsi = rdit0.deltaShiftIndex()
  boundSegs = dsi.findIdsEqual(1)  # boundary segs/faces of the M0 mesh
  m0descSkin = m0desc[boundSegs]
  fNodes = m0descSkin.computeFetchedNodeIds()  # fNodes needs dupl
  # In 3D, some points on the boundary of M0 will NOT be duplicated (where as in 2D, points on the boundary of M0 are always duplicated)
  # Think of a partial (plane) crack in a cube: the points at the tip of the crack and not located inside the volume of the cube are not duplicated
  # although they are technically on the skin of the cube.
  if this.getMeshDimension() == 3 :
      m0descSkinDesc, _, _, _, _ = m0descSkin.buildDescendingConnectivity() # all segments of the skin of the 3D (M0) mesh
      _, corresp = meshM2.areCellsIncludedIn(m0descSkinDesc,2)
      # validIds is the list of segments which are on both the skin of *this*, and in the segments of the M1 group
      # In the cube example above, this is a U shape polyline.
      validIds = corresp.findIdsInRange(0, meshM2.getNumberOfCells())
      if validIds.getNumberOfTuples():
          # Build the set of segments which are: in the desc mesh of the skin of the 3D mesh (M0) **and** in the desc mesh of the M1 group:
          # (the U-shaped polyline described above)
          m1IntersecSkin = m0descSkinDesc[validIds]
          # Its boundary nodes should no be duplicated (this is for example the tip of the crack inside the cube described above)
          notDuplSkin = m1IntersecSkin.findBoundaryNodes()
          fNodes1 = fNodes.buildSubstraction(notDuplSkin) # fNodes1 needs dupl

          # Specific logic to handle singular points :
          #   - a point on this U-shape line used in a cell which has no face in common with M1 is deemed singular.
          #   - indeed, if duplicated, such a point would lead to the duplication of a cell which has no face touching M1 ! The
          #   algorithm would be duplicating too much ...
          # This is a costly algorithm so only go into it if a simple (non sufficient) criteria is met: a node connected to more than 3 segs in meshM2:
          meshM2Desc, _, _, _, rdit0 = meshM2.buildDescendingConnectivity()  # a mesh made of node cells
          dsi = rdit0.deltaShiftIndex()
          singPoints = dsi.findIdsNotInRange(-1,4)   # points connected to (strictly) more than 3 segments
          if singPoints.getNumberOfTuples():
              print ("Hitting singular point logic")
              boundNodes = m1IntersecSkin.computeFetchedNodeIds()
              # If a point on this U-shape line is connected to cells which do not share any face with M1, then it
              # should not be duplicated
              #    1. Extract N D cells touching U-shape line:
              cellsAroundBN = this.getCellIdsLyingOnNodes(boundNodes, False)  # false= take cell in, even if not all nodes are in dupl
              mAroundBN = this[cellsAroundBN]
              mAroundBNDesc, descBN,descIBN,revDescBN,revDescIBN=mAroundBN.buildDescendingConnectivity()
              #    2. Identify cells in sub-mesh mAroundBN which have a face in common with M1
              _, idsOfM1BN = mAroundBNDesc.areCellsIncludedIn(otherDimM1OnSameCoords,2)
              nCells, nCellsDesc = mAroundBN.getNumberOfCells(), mAroundBNDesc.getNumberOfCells()
              idsTouch = DataArrayInt.New(); idsTouch.alloc(0,1)
              for v in idsOfM1BN:
                  if v[0] >= nCellsDesc:    # Keep valid match only
                      continue
                  idx0 = revDescIBN[v[0], 0]
                  c1, c2 = revDescBN[idx0, 0], revDescBN[idx0+1,0]
                  idsTouch.pushBackSilent(c1)
                  idsTouch.pushBackSilent(c2)
              #    3. Build complement
              idsTouchCompl = idsTouch.buildComplement(nCells)
              mAroundBNStrict = mAroundBN[idsTouchCompl]
              nod3 = mAroundBNStrict.computeFetchedNodeIds()
              inters = boundNodes.buildIntersection(nod3)
              print("sing,", inters.getValues())
              fNodes1 = fNodes1.buildSubstraction(inters)  # reminder: fNodes1 represent nodes that need dupl.
          notDup = xtrem.buildSubstraction(fNodes1)
      else:  # if validIds ...
        notDup = xtrem.buildSubstraction(fNodes)
  else:  # if 3D ...
    notDup = xtrem.buildSubstraction(fNodes)

  m1Nodes = otherDimM1OnSameCoords.computeFetchedNodeIds()
  dupl = m1Nodes.buildSubstraction(notDup)
  return dupl


def findCellsToRenumber(this, otherDimM1OnSameCoords, dupl):
  """  Find cells to renumber
  """
  # All N D cells touching our group (even when this is just one point touching)
  cellsAroundGroupLarge = this.getCellIdsLyingOnNodes(dupl, False)  # false= take cell in, even if not all nodes are in dupl
  #
  mAroundGrpLarge=this[cellsAroundGroupLarge]
  mArGrpLargeDesc,descL,descIL,revDescL,revDescIL=mAroundGrpLarge.buildDescendingConnectivity()
  mAroundGrpLarge.writeVTK("/tmp/mAr_large.vtu")
  mArGrpLargeDesc.writeVTK("/tmp/mAr_large_desc.vtu")

  # Extract now all N D cells which have a complete face in touch with the group:
  #    1. Identify cells of M1 group in sub-mesh mAroundGrp
  _, idsOfM1Large = mArGrpLargeDesc.areCellsIncludedIn(otherDimM1OnSameCoords,2)
  nL = mArGrpLargeDesc.getNumberOfCells()
  idsStrict = DataArrayInt.New(); idsStrict.alloc(0,1)
  #    2. Build map giving for each cell ID in mAroundGrp (not in mAroundGrpLarge) the corresponding cell
  #       ID on the other side of the crack:
  toOtherSide, pos = {}, {}
  cnt = 0
  for v in idsOfM1Large:
      if v[0] >= nL:    # Keep valid match only
          continue
      idx0 = revDescIL[v[0], 0]
      # Keep the two cells on either side of the face v of M1:
      c1, c2 = revDescL[idx0, 0], revDescL[idx0+1,0]
      if not c1 in idsStrict:
          pos[c1] = cnt
          idsStrict.pushBackSilent(c1)
          cnt += 1
      if not c2 in idsStrict:
          pos[c2] = cnt
          idsStrict.pushBackSilent(c2)
          cnt += 1
      k1, k2 = pos[c1], pos[c2]
      toOtherSide[k1] = k2
      toOtherSide[k2] = k1

  cellsAroundGroup = cellsAroundGroupLarge[idsStrict]
  mAroundGrp = this[cellsAroundGroup]
  nCells, nCellsLarge = cellsAroundGroup.getNumberOfTuples(), cellsAroundGroupLarge.getNumberOfTuples()
  mArGrpDesc,desc,descI,revDesc,revDescI=mAroundGrp.buildDescendingConnectivity()
  _, idsOfM1 = mArGrpDesc.areCellsIncludedIn(otherDimM1OnSameCoords,2) # TODO : could we avoid recomputing this??
  mAroundGrp.writeVTK("/tmp/mAr.vtu")
  mArGrpDesc.writeVTK("/tmp/mAr_desc.vtu")

  # Neighbor information of the mesh WITH the crack (some neighbors are removed):
  #     In the neighbor information remove the connection between high dimension cells and its low level constituents which are part
  #     of the frontier given in parameter (i.e. the cells of low dimension from the group delimiting the crack):
  DataArrayInt.RemoveIdsFromIndexedArrays(idsOfM1,desc,descI)
  #     Compute the neighbor of each cell in mAroundGrp, taking into account the broken link above. Two
  #     cells on either side of the crack (defined by the mesh of low dimension) are not neighbor anymore.
  neigh, neighI = MEDCouplingUMesh.ComputeNeighborsOfCellsAdv(desc,descI,revDesc,revDescI)

  # For each initial connex part of the M1 mesh (or said differently for each independent crack):
  seed, nIter, cnt = 0, 0, 0
  nIterMax = nCells+1 # Safety net for the loop
  hitCells = DataArrayInt.New(); hitCells.alloc(nCells)
  hitCells.fillWithValue(0)  # 0 : not hit, -1: one side of the crack, +1: other side of the crack
  MAX_CP = 10000 # the choices below assume we won't have more than 10000 different connex parts ...
  PING_FULL_init, PING_PART =  0, MAX_CP
  PONG_FULL_init, PONG_PART = -0,-MAX_CP
  while nIter < nIterMax:
#       print("dbg ", hitCells.getValues())
      t = hitCells.findIdsEqual(0)
      if not t.getNumberOfTuples():
        break
      seed = t[0,0]
      done = False
      cnt += 1
      PING_FULL = PING_FULL_init+cnt
      PONG_FULL = PONG_FULL_init-cnt
      while not done and nIter < nIterMax:  # Start of the ping-pong
          nIter += 1
          # Identify connex zone around the seed
          spreadZone, _ = MEDCouplingUMesh.ComputeSpreadZoneGraduallyFromSeed([seed],  neigh,neighI, -1)
          done = True
          for i, s in enumerate(spreadZone.getValues()):
              hitCells[s] = PING_FULL
              if s in toOtherSide:
                  other = toOtherSide[s]
                  if hitCells[other] != PONG_FULL:
                      done = False
                      hitCells[other] = PONG_PART
                      #  Compute next seed, i.e. a cell on the other side of the crack
                      seed = other
          if done:
              # we might have several disjoing PONG parts in front of a single PING connex part:
              idsPong = hitCells.findIdsEqual(PONG_PART)
              if idsPong.getNumberOfTuples():
                  seed = idsPong[0,0]
                  done = False
              continue  # continue without switching side (or break if done remains false)
          else:
              # Go the other side
              PING_FULL, PONG_FULL = PONG_FULL, PING_FULL
              PING_PART, PONG_PART = PONG_PART, PING_PART

      nonHitCells = hitCells.findIdsEqual(0)
      if nonHitCells.getNumberOfTuples():
        seed = nonHitCells[0,0]
      else:
        break

  if nIter >= nIterMax:
    raise ValueError("Too many iterations - should not happen")

  # Now we have handled all N D cells which have a face touching the M1 group. It remains the cells
  # which are just touching the group by one (or several) node(s):
  # All those cells are in direct contact with a cell which is either PING_FULL or PONG_FULL
  # So first reproject the PING/PONG info onto mAroundGrpLarge:
  hitCellsLarge = DataArrayInt.New(); hitCellsLarge.alloc(nCellsLarge)
  hitCellsLarge.fillWithValue(0)
  hitCellsLarge[idsStrict] = hitCells
  nonHitCells = hitCellsLarge.findIdsEqual(0)
  # Neighbor information in mAroundGrpLarge:
  neighL, neighIL = MEDCouplingUMesh.ComputeNeighborsOfCellsAdv(descL,descIL,revDescL,revDescIL)
  for c in nonHitCells:
      assert(False)
      neighs = neighL[neighIL[c[0]]:neighIL[c[0]+1]]
      for n in neighs:
          neighVal = hitCellsLarge[n[0]]
          if neighVal != 0 and abs(neighVal) < MAX_CP:  # (@test_T0) second part of the test to skip cells being assigned and target only cells assigned in the first part of the algo above
              currVal = hitCellsLarge[c[0]]
              if currVal != 0:   # Several neighbors have a candidate number
                  # Unfortunately in some weird cases (see testBuildInnerBoundary8) a cell in mAroundGrpLarge
                  # might have as neighbor two conflicting spread zone ...
                  if currVal*neighVal < 0:
                      # If we arrive here, the cell was already assigned a number and we found a neighbor with
                      # a different sign ... we must swap the whole spread zone!!
                      print("Ouch - must switch spread zones ...")
                      ids1 = hitCellsLarge.findIdsEqual(neighVal)
                      ids1b = hitCellsLarge.findIdsEqual(-neighVal)
                      ids2 = hitCellsLarge.findIdsEqual(MAX_CP*neighVal)
                      ids2b = hitCellsLarge.findIdsEqual(-MAX_CP*neighVal)
                      hitCellsLarge[ids1] *= -1
                      hitCellsLarge[ids1b] *= -1
                      hitCellsLarge[ids2] *= -1
                      hitCellsLarge[ids2b] *= -1
              else:  # First assignation
                  hitCellsLarge[c[0],0] = MAX_CP*neighVal   # Same sign, but different value to preserve PING_FULL and PONG_FULL

###
###  SO FAR THE LOGIC BELOW WAS NOT NEEDED ....
###

#   # Now handling remaining cells not touched by the above process, called "naked" cells (see cell #20 in mArndLarge in testBuildInnerBoundary8() ...)
#   naked = hitCellsLarge.findIdsEqual(0)
#   mLargC, mLargCI = mArGrpLargeDesc.getNodalConnectivity(),mArGrpLargeDesc.getNodalConnectivityIndex()
#   for c in naked:
#       neighs = neighL[neighIL[c[0]]:neighIL[c[0]+1]]  # ExtractFromIndexedArray?
#       nbDup = {}
#       fac1 = descL[descIL[c[0]]:descIL[c[0]+1]]
#       for n in neighs:
#           if hitCellsLarge[n[0]] == 0:
#               continue   # this neighbour is naked too, nothing we can do for now
#           # Among the values found on neighbour cells, take the one from the neighbour which is connected
#           # with the most "economical" face, i.e. the face made of a minimal number of duplicated points.
#           # TODO: this is a shaky criteria ... find sth more robust ...
#           #   1. find face(s) making the link
#           fac2 = descL[descIL[n[0]]:descIL[n[0]+1]]
#           com = fac1.buildIntersection(fac2)
#           if (com.getNumberOfTuples() == 0):
#               raise ValueError("Internal error : no common face ?")
#           #   2. count number of duplicated node for this face.
#           for f in com: # for all common faces
#               faceNodes = mLargC[mLargCI[f[0]]+1:mLargCI[f[0]+1]] # first +1 to skip type
#               comNod = faceNodes.buildIntersection(dupl)
#               # in case the two cells are in contact by multiple faces, take the most conservative value
#               nbDup[n[0]] = max(nbDup.get(n[0],-1), comNod.getNumberOfTuples())
#       # Minimal value in nbDup?
#       cellIdx = min(nbDup, key=nbDup.get)
#       hitCellsLarge[c[0]] = hitCellsLarge[cellIdx]
#
#   cellsToModifyConn0_torenum = hitCellsLarge.findIdsInRange(1,MAX_CP)    # Positive spread zone number
#   cellsToModifyConn1_torenum = hitCellsLarge.findIdsInRange(-MAX_CP, 0)  # Negative spread zone number


###
###  C++ VERSION OF IT
###
# //
# //  // Now handling remaining cells not touched by the for loop above, called "naked" cells (see cell #20 in mArndGrpLarge in testBuildInnerBoundary8() ...)
# //  DAInt naked = hitCellsLarge->findIdsEqual(0)
# //  const mcIdType *mLargCP=mArGrpLargeDesc->getNodalConnectivity()->begin(), *mLargCIP=mArGrpLargeDesc->getNodalConnectivityIndex()->begin()
# //  for (const auto &c: *naked)
# //    {
# //      std::map<mcIdType, mcIdType> nbDup
# //      // Retrieve list of faces of cell c
# //      mcIdType nbFac1=descILP[c+1]-descILP[c]
# //      std::vector<mcIdType> fac1(nbFac1)
# //      std::copy(descLP+descILP[c], descLP+descILP[c+1], fac1.begin())
# //      std::sort(fac1.begin(), fac1.end())
# //      mcIdType cnt00 = neighILP[c]
# //      for (const mcIdType *n=neighLP+cnt00; cnt00 < neighILP[c+1]; n++, cnt00++)
# //        {
# //          if (hitCellsLargeP[*n] == 0)
# //            continue;   // this neighbour is naked too, nothing we can do for now
# //          // Among the values found on neighbour cells, take the one from the neighbour which is connected
# //          // with the most "economical" face, i.e. the face made of a minimal number of duplicated points.
# //          // TODO: this is a shaky criteria ... find sth more robust ...
# //          //   1. find face(s) making the link
# //          mcIdType nbFac2=descILP[*n+1]-descILP[*n]
# //          std::vector<mcIdType> fac2(nbFac2)
# //          std::copy(descLP+descILP[*n], descLP+descILP[*n+1], fac2.begin())
# //          std::sort(fac2.begin(), fac2.end())
# //          std::vector<mcIdType> comFac
# //          std::set_intersection(fac1.begin(), fac1.end(),
# //                                fac2.begin() ,fac2.end(),
# //                                std::back_inserter(comFac))
# //          if (comFac.size() == 0)
# //            throw INTERP_KERNEL::Exception("MEDCouplingUMesh::findCellsToRenumber: internal error no common face between two cells should not happen")
# //          //   2. count number of duplicated node for this face.
# //          for (const auto &f : comFac) // for all common faces
# //            {
# //              std::vector<mcIdType> comNod
# //              std::set_intersection(nodeIdsToDuplicateBg, nodeIdsToDuplicateEnd,
# //                                    mLargCP+mLargCIP[f]+1, mLargCP+mLargCIP[f+1],    // first +1 to skip type in connectivity
# //                                    std::back_inserter(comNod))
# //              // in case the two cells are in contact by multiple faces, take the most conservative value
# //              mcIdType val=-1
# //              if(nbDup.find(*n) != nbDup.end()) val=nbDup[*n]
# //              nbDup[*n] = std::max(val, (mcIdType)comNod.size())
# //            }
# //        }
# //      // Minimal value in nbDup?
# //      using PairId = std::pair<mcIdType, mcIdType>
# //      auto comp_fonc = [](const PairId& p1, const PairId& p2) { return p1.second < p2.second; }
# //      PairId zemin = *min_element(nbDup.begin(), nbDup.end(), comp_fonc)
# //      hitCellsLargeP[c] = hitCellsLargeP[zemin.first]
# //    }


  cellsToModifyConn0_torenum = hitCellsLarge.findIdsInRange(1,MAX_CP*MAX_CP)    # Positive spread zone number
  cellsToModifyConn1_torenum = hitCellsLarge.findIdsInRange(-MAX_CP*MAX_CP, 0)  # Negative spread zone number
  if cellsToModifyConn0_torenum.getNumberOfTuples() + cellsToModifyConn1_torenum.getNumberOfTuples() != cellsAroundGroupLarge.getNumberOfTuples():
      raise ValueError("Some cells not hit - Internal error should not happen")
  cellsToModifyConn0_torenum.transformWithIndArr(cellsAroundGroupLarge)
  cellsToModifyConn1_torenum.transformWithIndArr(cellsAroundGroupLarge)
  #
  cellIdsNeededToBeRenum=cellsToModifyConn0_torenum
  cellIdsNotModified=cellsToModifyConn1_torenum

  return cellIdsNeededToBeRenum, cellIdsNotModified

