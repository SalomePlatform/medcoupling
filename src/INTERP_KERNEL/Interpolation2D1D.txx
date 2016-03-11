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
#ifndef __INTERPOLATION2D1D_TXX__
#define __INTERPOLATION2D1D_TXX__

#include "Interpolation2D1D.hxx"

namespace INTERP_KERNEL
{

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
  template<class MyMeshType, class MatrixType>
  int Interpolation2D1D::interpolateMeshes(const MyMeshType& myMeshS, const MyMeshType& myMeshT, MatrixType& result, const std::string& method)
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

    /**************************************************/
    /* Search the characteristic size of the meshes   */
    /**************************************************/

    int printLevel = InterpolationOptions::getPrintLevel();
    _dim_caracteristic = CalculateCharacteristicSizeOfMeshes(myMeshS, myMeshT, printLevel);
    if (printLevel>=1)
      {
        std::cout << "Interpolation2D1D::computation of the intersections" << std::endl;
      }

    PlanarIntersector<MyMeshType,MatrixType>* intersector=0;
    std::string meth = InterpolationOptions::filterInterpolationMethod(method);
    if(meth=="P0P0")
      {
        switch (InterpolationOptions::getIntersectionType())
          {
          case Geometric2D:
            intersector=new Geometric2DIntersector<MyMeshType,MatrixType,Planar2D1DIntersectorP0P0>(myMeshT, myMeshS, _dim_caracteristic,
                                                                                                    InterpolationOptions::getMaxDistance3DSurfIntersect(),
                                                                                                    InterpolationOptions::getMinDotBtwPlane3DSurfIntersect(),
                                                                                                    InterpolationOptions::getMedianPlane(),
                                                                                                    InterpolationOptions::getPrecision(),
                                                                                                    InterpolationOptions::getOrientation());
            break;
          default:
            throw INTERP_KERNEL::Exception("Invalid intersection type ! Must be : Geometric2D");
          }
      }
    else
      throw INTERP_KERNEL::Exception("Invalid method specified or intersection type ! Must be : \"P0P0\"");

    /****************************************************************/
    /* Create a search tree based on the bounding boxes             */
    /* Instanciate the intersector and initialise the result vector */
    /****************************************************************/

    long start_filtering=clock();

    std::vector<double> bbox;
    intersector->createBoundingBoxes(myMeshS,bbox); // create the bounding boxes
    const double *bboxPtr=0;
    if(nbMailleS>0)
      bboxPtr=&bbox[0];
    BBTree<SPACEDIM,ConnType> my_tree(bboxPtr, 0, 0,nbMailleS, -getPrecision());//creating the search structure

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

    const DuplicateFacesType& intersectFaces = *intersector->getIntersectFaces();
    DuplicateFacesType::const_iterator iter;
    for (iter = intersectFaces.begin(); iter != intersectFaces.end(); ++iter)
      {
        if (iter->second.size() > 1)
          {
            _duplicate_faces.insert(std::make_pair(iter->first, iter->second));
          }
      }

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
