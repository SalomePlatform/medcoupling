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

#ifndef __INTERPOLATIONPLANAR_TXX__
#define __INTERPOLATIONPLANAR_TXX__

#include "InterpolationPlanar.hxx"
#include "Interpolation.txx"
#include "InterpolationOptions.hxx"
#include "PlanarIntersector.hxx"
#include "PlanarIntersector.txx"
#include "TriangulationIntersector.hxx"
#include "TriangulationIntersector.txx"
#include "ConvexIntersector.hxx"
#include "ConvexIntersector.txx"
#include "Geometric2DIntersector.hxx"
#include "Geometric2DIntersector.txx"
#include "PointLocator2DIntersector.hxx"
#include "PointLocator2DIntersector.txx"
#include "PlanarIntersectorP0P1PL.hxx"
#include "PlanarIntersectorP0P1PL.txx"
#include "PlanarIntersectorP1P0PL.hxx"
#include "PlanarIntersectorP1P0PL.txx"
#include "PlanarIntersectorP1P1PL.hxx"
#include "PlanarIntersectorP1P1PL.txx"
#include "MappedBarycentric2DIntersectorP1P1.hxx"
#include "MappedBarycentric2DIntersectorP1P1.txx"
#include "VectorUtils.hxx"
#include "BBTree.txx"

#include <limits>
#include <time.h>

namespace INTERP_KERNEL
{
  /**
   * \defgroup interpolationPlanar InterpolationPlanar
   *
   * \class InterpolationPlanar
   * \brief Class used to compute the coefficients of the interpolation matrix between 
   * two local meshes in two dimensions. Meshes can contain mixed triangular and quadrangular elements.
   */
  template<class RealPlanar>
  InterpolationPlanar<RealPlanar>::InterpolationPlanar():_dim_caracteristic(1)
                                                         
  {
  }

  template<class RealPlanar>
  InterpolationPlanar<RealPlanar>::InterpolationPlanar(const InterpolationOptions& io):Interpolation< InterpolationPlanar<RealPlanar> >(io),_dim_caracteristic(1)
                                                         
  {
  }

  /**
   *  \brief  Function used to set the options for the intersection calculation
   * \details The following options can be modified:
   *  -# Intersection_type: the type of algorithm to be used in the computation of the cell-cell intersections.
   *   - Values: Triangle, Convex.
   *   - Default: Triangle.
   *  -# Precision: Level of precision of the computations is precision times the characteristic size of the mesh.
   *   - Values: positive real number.
   *   - Default: 1.0E-12.
   *  -# PrintLevel: Level of verboseness during the computations.
   *   - Values: interger between 0 and 3.
   *   - Default: 0.
   */
  template<class RealPlanar>
  void InterpolationPlanar<RealPlanar>::setOptions(double precision, int printLevel, IntersectionType intersectionType, int orientation)
  {
    InterpolationOptions::setPrecision(precision);
    InterpolationOptions::setPrintLevel(printLevel);
    InterpolationOptions::setIntersectionType(intersectionType);
    InterpolationOptions::setOrientation(orientation);
  }
  
  
  /** \brief Main function to interpolate triangular or quadrangular meshes.
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
      * @return            vector containing for each element i of the source mesh, a map giving for each element j
      *                    of the target mesh which i intersects, the area of the intersection
      *
      */
  template<class RealPlanar>
  template<class MyMeshType, class MatrixType>
  int InterpolationPlanar<RealPlanar>::interpolateMeshes(const MyMeshType& myMeshS, const MyMeshType& myMeshT, MatrixType& result, const std::string& method)
  {
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;

    long global_start =clock();
    int counter=0;   
    /***********************************************************/
    /* Check both meshes are made of triangles and quadrangles */
    /***********************************************************/

    long nbMailleS=myMeshS.getNumberOfElements();
    long nbMailleT=myMeshT.getNumberOfElements();
    
    /**************************************************/
    /* Search the characteristic size of the meshes   */
    /**************************************************/
    
    double BoxS[2*SPACEDIM]; myMeshS.getBoundingBox(BoxS);
    double BoxT[2*SPACEDIM]; myMeshT.getBoundingBox(BoxT);
    double diagonalS,dimCaracteristicS=std::numeric_limits<double>::max();
    if(nbMailleS!=0)
      {
        diagonalS=getDistanceBtw2Pts<SPACEDIM>(BoxS+SPACEDIM,BoxS);
        dimCaracteristicS=diagonalS/nbMailleS;
      }
    double diagonalT,dimCaracteristicT=std::numeric_limits<double>::max();
    if(nbMailleT!=0)
      {
        diagonalT=getDistanceBtw2Pts<SPACEDIM>(BoxT+SPACEDIM,BoxT);
        dimCaracteristicT=diagonalT/nbMailleT;
      }
    
    _dim_caracteristic=std::min(dimCaracteristicS, dimCaracteristicT);
    if (InterpolationOptions::getPrintLevel()>=1)
      {
        std::cout << "  - Characteristic size of the source mesh : " << dimCaracteristicS << std::endl;
        std::cout << "  - Characteristic size of the target mesh: " << dimCaracteristicT << std::endl;
        std::cout << "InterpolationPlanar::computation of the intersections" << std::endl;
      }
    
    PlanarIntersector<MyMeshType,MatrixType>* intersector=0;
    std::string meth = InterpolationOptions::filterInterpolationMethod(method);
    if(meth=="P0P0")
      {
        switch (InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector=new TriangulationIntersector<MyMeshType,MatrixType,PlanarIntersectorP0P0>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                                  InterpolationOptions::getPrecision(),
                                                                                                  InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                  InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                  InterpolationOptions::getMedianPlane(),
                                                                                                  InterpolationOptions::getOrientation(),
                                                                                                  InterpolationOptions::getPrintLevel());
            break;
          case Convex:
            intersector=new ConvexIntersector<MyMeshType,MatrixType,PlanarIntersectorP0P0>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                           InterpolationOptions::getPrecision(),
                                                                                           InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                           InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                           InterpolationOptions::getMedianPlane(),
                                                                                           InterpolationOptions::getDoRotate(),
                                                                                           InterpolationOptions::getOrientation(),
                                                                                           InterpolationOptions::getPrintLevel());
            break;
          case Geometric2D:
            intersector=new Geometric2DIntersector<MyMeshType,MatrixType,PlanarIntersectorP0P0>(myMeshT, myMeshS, _dim_caracteristic,
                                                                                                InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                InterpolationOptions::getMedianPlane(),
                                                                                                InterpolationOptions::getPrecision(),
                                                                                                InterpolationOptions::getOrientation());
            break;
          case PointLocator:
            intersector=new PointLocator2DIntersector<MyMeshType,MatrixType,PlanarIntersectorP0P0>(myMeshT, myMeshS, _dim_caracteristic,
                                                                                                   InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                   InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                   InterpolationOptions::getMedianPlane(),
                                                                                                   InterpolationOptions::getPrecision(),
                                                                                                   InterpolationOptions::getOrientation());
            break;
          default:
            throw INTERP_KERNEL::Exception("For P0P0 planar interpolation possibities are : Triangulation, Convex, Geometric2D, PointLocator !");
          }
      }
    else if(meth=="P0P1")
      {
        switch (InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector=new TriangulationIntersector<MyMeshType,MatrixType,PlanarIntersectorP0P1>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                                  InterpolationOptions::getPrecision(),
                                                                                                  InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                  InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                  InterpolationOptions::getMedianPlane(),
                                                                                                  InterpolationOptions::getOrientation(),
                                                                                                  InterpolationOptions::getPrintLevel());
            break;
          case Convex:
            intersector=new ConvexIntersector<MyMeshType,MatrixType,PlanarIntersectorP0P1>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                           InterpolationOptions::getPrecision(),
                                                                                           InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                           InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                           InterpolationOptions::getMedianPlane(),
                                                                                           InterpolationOptions::getDoRotate(),
                                                                                           InterpolationOptions::getOrientation(),
                                                                                           InterpolationOptions::getPrintLevel());
            break;
          case Geometric2D:
            intersector=new Geometric2DIntersector<MyMeshType,MatrixType,PlanarIntersectorP0P1>(myMeshT, myMeshS, _dim_caracteristic,
                                                                                                InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                InterpolationOptions::getMedianPlane(),
                                                                                                InterpolationOptions::getPrecision(),
                                                                                                InterpolationOptions::getOrientation());
            break;
          case PointLocator:
            intersector=new PlanarIntersectorP0P1PL<MyMeshType,MatrixType>(myMeshT, myMeshS, _dim_caracteristic,
                                                                           InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                           InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                           InterpolationOptions::getMedianPlane(),
                                                                           InterpolationOptions::getPrecision(),
                                                                           InterpolationOptions::getOrientation());
            break;
          case Barycentric:
            intersector=new TriangulationIntersector<MyMeshType,MatrixType,PlanarIntersectorP0P1Bary>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                                      InterpolationOptions::getPrecision(),
                                                                                                      InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                      InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                      InterpolationOptions::getMedianPlane(),
                                                                                                      InterpolationOptions::getOrientation(),
                                                                                                      InterpolationOptions::getPrintLevel());
            break;
          case BarycentricGeo2D:
            intersector=new Geometric2DIntersector<MyMeshType,MatrixType,PlanarIntersectorP0P1Bary>(myMeshT, myMeshS, _dim_caracteristic,
                                                                                                    InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                    InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                    InterpolationOptions::getMedianPlane(),
                                                                                                    InterpolationOptions::getPrecision(),
                                                                                                    InterpolationOptions::getOrientation());
            break;
          default:
            throw INTERP_KERNEL::Exception("For P0P1 planar interpolation possibities are : Triangulation, Convex, Geometric2D, PointLocator, Barycentric, BarycentricGeo2D !");
          }
      }
    else if(meth=="P1P0")
      {
        switch (InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector=new TriangulationIntersector<MyMeshType,MatrixType,PlanarIntersectorP1P0>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                                  InterpolationOptions::getPrecision(),
                                                                                                  InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                  InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                  InterpolationOptions::getMedianPlane(),
                                                                                                  InterpolationOptions::getOrientation(),
                                                                                                  InterpolationOptions::getPrintLevel());
            break;
          case Convex:
            intersector=new ConvexIntersector<MyMeshType,MatrixType,PlanarIntersectorP1P0>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                           InterpolationOptions::getPrecision(),
                                                                                           InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                           InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                           InterpolationOptions::getMedianPlane(),
                                                                                           InterpolationOptions::getDoRotate(),
                                                                                           InterpolationOptions::getOrientation(),
                                                                                           InterpolationOptions::getPrintLevel());
            break;
          case Geometric2D:
            intersector=new Geometric2DIntersector<MyMeshType,MatrixType,PlanarIntersectorP1P0>(myMeshT, myMeshS, _dim_caracteristic,
                                                                                                InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                InterpolationOptions::getMedianPlane(),
                                                                                                InterpolationOptions::getPrecision(),
                                                                                                InterpolationOptions::getOrientation());
            break;
          case PointLocator:
            intersector=new PlanarIntersectorP1P0PL<MyMeshType,MatrixType>(myMeshT, myMeshS, _dim_caracteristic,
                                                                           InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                           InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                           InterpolationOptions::getMedianPlane(),
                                                                           InterpolationOptions::getPrecision(),
                                                                           InterpolationOptions::getOrientation());
            break;
          case Barycentric:
             intersector=new TriangulationIntersector<MyMeshType,MatrixType,PlanarIntersectorP1P0Bary>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                                       InterpolationOptions::getPrecision(),
                                                                                                       InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                       InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                       InterpolationOptions::getMedianPlane(),
                                                                                                       InterpolationOptions::getOrientation(),
                                                                                                       InterpolationOptions::getPrintLevel());
             break;
          case BarycentricGeo2D:
            intersector=new Geometric2DIntersector<MyMeshType,MatrixType,PlanarIntersectorP1P0Bary>(myMeshT, myMeshS, _dim_caracteristic,
                                                                                                    InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                    InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                    InterpolationOptions::getMedianPlane(),
                                                                                                    InterpolationOptions::getPrecision(),
                                                                                                    InterpolationOptions::getOrientation());
            break;
          default:
            throw INTERP_KERNEL::Exception("For P1P0 planar interpolation possibities are : Triangulation, Convex, Geometric2D, PointLocator, BarycentricGeo2D or Barycentric!");
          }
      }
    else if(meth=="P1P1")
      {
        switch (InterpolationOptions::getIntersectionType())
          {
          case Triangulation:
            intersector=new TriangulationIntersector<MyMeshType,MatrixType,PlanarIntersectorP1P1>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                                  InterpolationOptions::getPrecision(),
                                                                                                  InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                  InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                  InterpolationOptions::getMedianPlane(),
                                                                                                  InterpolationOptions::getOrientation(),
                                                                                                  InterpolationOptions::getPrintLevel());
            break;
          case Convex:
            intersector=new ConvexIntersector<MyMeshType,MatrixType,PlanarIntersectorP1P1>(myMeshT,myMeshS,_dim_caracteristic,
                                                                                           InterpolationOptions::getPrecision(),
                                                                                           InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                           InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                           InterpolationOptions::getMedianPlane(),
                                                                                           InterpolationOptions::getDoRotate(),
                                                                                           InterpolationOptions::getOrientation(),
                                                                                           InterpolationOptions::getPrintLevel());
            break;
          case Geometric2D:
            intersector=new Geometric2DIntersector<MyMeshType,MatrixType,PlanarIntersectorP1P1>(myMeshT, myMeshS, _dim_caracteristic,
                                                                                                InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                InterpolationOptions::getMedianPlane(),
                                                                                                InterpolationOptions::getPrecision(),
                                                                                                InterpolationOptions::getOrientation());
            break;
          case PointLocator:
            intersector=new PlanarIntersectorP1P1PL<MyMeshType,MatrixType>(myMeshT, myMeshS, _dim_caracteristic,
                                                                           InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                           InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                           InterpolationOptions::getMedianPlane(),
                                                                           InterpolationOptions::getPrecision(),
                                                                           InterpolationOptions::getOrientation());
            break;
          case MappedBarycentric:
            intersector=new MappedBarycentric2DIntersectorP1P1<MyMeshType,MatrixType>(myMeshT, myMeshS, _dim_caracteristic,
                                                                                     InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                     InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                     InterpolationOptions::getMedianPlane(),
                                                                                     InterpolationOptions::getPrecision(),
                                                                                     InterpolationOptions::getOrientation());
            break;
          default:
            throw INTERP_KERNEL::Exception("For P1P1 planar interpolation possibities are : Triangulation, Convex, Geometric2D, PointLocator, MappedBarycentric !");
          }
      }
    else
      throw INTERP_KERNEL::Exception("Invalid method specified or intersection type ! Must be in : \"P0P0\" \"P0P1\" \"P1P0\" or \"P1P1\"");
    /****************************************************************/
    /* Create a search tree based on the bounding boxes             */
    /* Instanciate the intersector and initialise the result vector */
    /****************************************************************/
 
    long start_filtering=clock();
 
    std::vector<double> bbox;
    intersector->createBoundingBoxes(myMeshS,bbox); // create the bounding boxes
    performAdjustmentOfBB(intersector,bbox);
    const double *bboxPtr=0;
    if(nbMailleS>0)
      bboxPtr=&bbox[0];
    BBTree<SPACEDIM,ConnType> my_tree(bboxPtr, 0, 0,nbMailleS);//creating the search structure 

    long end_filtering=clock();

    result.resize(intersector->getNumberOfRowsOfResMatrix());//on initialise.

    /****************************************************/
    /* Loop on the target cells - core of the algorithm */
    /****************************************************/
    long start_intersection=clock();
    long nbelem_type=myMeshT.getNumberOfElements();
    const ConnType *connIndxT=myMeshT.getConnectivityIndexPtr();
    for(int iT=0; iT<nbelem_type; iT++)
      {
        int nb_nodesT=connIndxT[iT+1]-connIndxT[iT];
        std::vector<int> intersecting_elems;
        double bb[2*SPACEDIM];
        intersector->getElemBB(bb,myMeshT,OTT<ConnType,numPol>::indFC(iT),nb_nodesT);
        my_tree.getIntersectingElems(bb, intersecting_elems);
        intersector->intersectCells(iT,intersecting_elems,result);
        counter+=intersecting_elems.size();
        intersecting_elems.clear();
      }
    int ret=intersector->getNumberOfColsOfResMatrix();
    delete intersector;

    if (InterpolationOptions::getPrintLevel() >=1)
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
