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
// Author : Anthony Geay (CEA/DEN)
#ifndef __INTERPOLATIONCURVE_TXX__
#define __INTERPOLATIONCURVE_TXX__

#include "InterpolationCurve.hxx"
#include "InterpolationOptions.hxx"
#include "CurveIntersectorP0P0.txx"
#include "CurveIntersectorP1P0.txx"
#include "CurveIntersectorP0P1.txx"
#include "CurveIntersectorP1P1.txx"
#include "CurveIntersectorP1P1PL.txx"
#include "BBTree.txx"

#include <time.h>

namespace INTERP_KERNEL
{
  /**
   * \defgroup interpolationCurve InterpolationCurve
   *
   * \class InterpolationCurve
   * \brief Class used to compute the coefficients of the interpolation matrix between 
   * two local meshes in two dimensions.
   */
  template<class RealCurve>
  InterpolationCurve<RealCurve>::InterpolationCurve()
  {
  }

  template<class RealCurve>
  InterpolationCurve<RealCurve>::InterpolationCurve (const InterpolationOptions& io)
    :Interpolation< InterpolationCurve<RealCurve> >(io)
  {
  }

  /** \brief Main function to interpolate 1D meshes.
      \details  The algorithm proceeds in two steps: first a filtering process reduces the number of pairs of elements for which the
      * calculation must be carried out by eliminating pairs that do not intersect based on their bounding boxes. Then, the 
      * volume of intersection is calculated by an object of type IntersectorPlanar for the remaining pairs, and entered into the
      * intersection matrix. 
      * 
      * The matrix is partially sparse : it is a vector of maps of integer - double pairs. 
      * The length of the vector is equal to the number of target elements - for each target element there is a map, regardless
      * of whether the element intersects any source elements or not. But in the maps there are only entries for those source elements
      * which have a non-zero intersection volume with the target element. The vector has indices running from 
      * 0 to (#target elements - 1), meaning that the map for target element i is stored at index i - 1. In the maps, however,
      * the indexing is more natural : the intersection volume of the target element i with source element j is found at matrix[i-1][j].
      * 
   
      * @param myMeshS  Planar source mesh
      * @Param myMeshT  Planar target mesh
      * @return vector containing for each element i of the source mesh, a map giving for each element j
      *         of the target mesh which i intersects, the area of the intersection
      *
      */
  template<class RealCurve>
  template<class MyMeshType, class MatrixType>
  int InterpolationCurve<RealCurve>::interpolateMeshes (const MyMeshType&  myMeshS,
                                                        const MyMeshType&  myMeshT,
                                                        MatrixType&        result,
                                                        const std::string& method)
  {
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol = MyMeshType::My_numPol;

    long global_start = clock();
    int counter=0;   

    long nbMailleS = myMeshS.getNumberOfElements();
    long nbMailleT = myMeshT.getNumberOfElements();
    
    CurveIntersector<MyMeshType,MatrixType>* intersector=0;
    if(method=="P0P0")
      {
        switch (InterpolationOptions::getIntersectionType())
          {
            case Triangulation:
              {
                intersector = new CurveIntersectorP0P0<MyMeshType,MatrixType>(myMeshT, myMeshS,
                                                                              InterpolationOptions::getPrecision(),
                                                                              InterpolationOptions::getBoundingBoxAdjustmentAbs(),
                                                                              InterpolationOptions::getMedianPlane(),
                                                                              InterpolationOptions::getPrintLevel());
                break;
              }
            default:
              throw INTERP_KERNEL::Exception("For P0P0 in 1D or 2D curve only Triangulation supported for the moment !");
          }
      }
    else if(method=="P0P1")
      {
        switch (InterpolationOptions::getIntersectionType())
          {
            case Triangulation:
              {
                intersector = new CurveIntersectorP0P1<MyMeshType,MatrixType>(myMeshT, myMeshS,
                                                                              InterpolationOptions::getPrecision(),
                                                                              InterpolationOptions::getBoundingBoxAdjustmentAbs(),
                                                                              InterpolationOptions::getMedianPlane(),
                                                                              InterpolationOptions::getPrintLevel());
                break;
              }
            default:
              throw INTERP_KERNEL::Exception("For P0P1 in 1D or 2D curve only Triangulation supported for the moment !");
          }
      }
    else if(method=="P1P0")
      {
        switch (InterpolationOptions::getIntersectionType())
          {
            case Triangulation:
              {
                intersector = new CurveIntersectorP1P0<MyMeshType,MatrixType>(myMeshT, myMeshS,
                                                                              InterpolationOptions::getPrecision(),
                                                                              InterpolationOptions::getBoundingBoxAdjustmentAbs(),
                                                                              InterpolationOptions::getMedianPlane(),
                                                                              InterpolationOptions::getPrintLevel());
                break;
              }
            default:
              throw INTERP_KERNEL::Exception("For P1P0 in 1D or 2D curve only Triangulation supported for the moment !");
          }
      }
    else if(method=="P1P1")
      {
        switch (InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector = new CurveIntersectorP1P1<MyMeshType,MatrixType>
              (myMeshT, myMeshS,
               InterpolationOptions::getPrecision(),
               InterpolationOptions::getBoundingBoxAdjustmentAbs(),
               InterpolationOptions::getMedianPlane(),
               InterpolationOptions::getPrintLevel());
            break;
          case PointLocator:
            intersector = new CurveIntersectorP1P1PL<MyMeshType,MatrixType>
              (myMeshT, myMeshS,
               InterpolationOptions::getPrecision(),
               InterpolationOptions::getBoundingBoxAdjustmentAbs(),
               InterpolationOptions::getMedianPlane(),
               InterpolationOptions::getPrintLevel());
            break;
          default:
            throw INTERP_KERNEL::Exception("For P1P1 in 1D or 2D curve only Triangulation and PointLocator supported !");
          }
      }
    else
      throw INTERP_KERNEL::Exception("Invalid method specified ! Must be in : \"P0P0\" \"P0P1\" \"P1P0\" or \"P1P1\"");
    /****************************************************************/
    /* Create a search tree based on the bounding boxes             */
    /* Instanciate the intersector and initialise the result vector */
    /****************************************************************/
 
    long start_filtering=clock();
 
    std::vector<double> bbox;
    intersector->createBoundingBoxes(myMeshS,bbox); // create the bounding boxes
    intersector->adjustBoundingBoxes(bbox, InterpolationOptions::getBoundingBoxAdjustmentAbs());
    BBTree<SPACEDIM,ConnType> my_tree(&bbox[0], 0, 0, nbMailleS);//creating the search structure 

    long end_filtering = clock();

    result.resize(intersector->getNumberOfRowsOfResMatrix());//on initialise.

    /****************************************************/
    /* Loop on the target cells - core of the algorithm */
    /****************************************************/
    long start_intersection = clock();
    const ConnType *connIndxT = myMeshT.getConnectivityIndexPtr();
    for(int iT=0; iT<nbMailleT; iT++)
      {
        int nb_nodesT = connIndxT[iT+1] - connIndxT[iT];
        std::vector<int> intersecting_elems;
        double bb[2*SPACEDIM];
        intersector->getElemBB(bb,myMeshT,OTT<ConnType,numPol>::indFC(iT),nb_nodesT);
        my_tree.getIntersectingElems(bb, intersecting_elems);
        intersector->intersectCells(iT,intersecting_elems,result);
        counter += intersecting_elems.size();
      }
    int ret = intersector->getNumberOfColsOfResMatrix();
    delete intersector;
    
    if (InterpolationOptions::getPrintLevel() >= 1)
      {
        long end_intersection=clock();
        std::cout << "Filtering time= " << end_filtering-start_filtering << std::endl;
        std::cout << "Intersection time= " << end_intersection-start_intersection << std::endl;
        long global_end =clock();    
        std::cout << "Number of computed intersections = " << counter << std::endl;
        std::cout << "Global time= " << global_end - global_start << std::endl;
      }
    return ret;
  }
}

#endif
