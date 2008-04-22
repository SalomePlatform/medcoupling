#ifndef __INTERPOLATIONPLANAR_TXX__
#define __INTERPOLATIONPLANAR_TXX__

#include "InterpolationPlanar.hxx"
#include "PlanarIntersector.hxx"
#include "PlanarIntersector.txx"
#include "TriangulationIntersector.hxx"
#include "TriangulationIntersector.txx"
#include "ConvexIntersector.hxx"
#include "ConvexIntersector.txx"
#include "Geometric2DIntersector.hxx"
#include "Geometric2DIntersector.txx"
//#include "GenericIntersector.hxx"
//#include "GenericIntersector.cxx"
#include "VectorUtils.hxx"
#include "BBTree.txx"
#include<time.h>

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
  InterpolationPlanar<RealPlanar>::InterpolationPlanar():_printLevel(0),_precision(DEFAULT_PRECISION),_dimCaracteristic(1),
                                                         _intersectionType(Triangulation)
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
  void InterpolationPlanar<RealPlanar>::setOptions(double precision, int printLevel, IntersectionType intersectionType)
  {
    _precision=precision;
    _printLevel=printLevel;
    _intersectionType=intersectionType;
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
   
      * @param myMesh_S  Planar source mesh
      * @Param myMesh_P  Planar target mesh
      * @return            vector containing for each element i of the source mesh, a map giving for each element j
      *                    of the target mesh which i intersects, the area of the intersection
      *
      */
  template<class RealPlanar>
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MatrixType, class MyMeshType>
  void InterpolationPlanar<RealPlanar>::interpolateMeshes(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& myMesh_S,
                                                          const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& myMesh_P,
                                                          MatrixType& result)
  {
    long global_start =clock();
    int counter=0;   
    /***********************************************************/
    /* Check both meshes are made of triangles and quadrangles */
    /***********************************************************/

    long nbMaille_S =myMesh_S.getNumberOfElements();
    long nbMaille_P =myMesh_P.getNumberOfElements();
    
    /**************************************************/
    /* Search the characteristic size of the meshes   */
    /**************************************************/
    
    double BoxS[2*SPACEDIM]; myMesh_S.getBoundingBox(BoxS);
    double BoxP[2*SPACEDIM]; myMesh_P.getBoundingBox(BoxP);
    double diagonalS=getDistanceBtw2Pts<SPACEDIM>(BoxS+SPACEDIM,BoxS);
    double DimCaracteristic_S=diagonalS/nbMaille_S;
    double diagonalP=getDistanceBtw2Pts<SPACEDIM>(BoxP+SPACEDIM,BoxP);
    double DimCaracteristic_P=diagonalP/nbMaille_P;
    
    _dimCaracteristic=std::min(DimCaracteristic_S, DimCaracteristic_P);
    if (_printLevel>=1)
      {
        std::cout << "  - Characteristic size of the source mesh : " << DimCaracteristic_S << std::endl;
        std::cout << "  - Characteristic size of the target mesh: " << DimCaracteristic_P << std::endl;
        std::cout << "InterpolationPlanar::computation of the intersections" << std::endl;
      }
    
    PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>* intersector;
    
    switch (_intersectionType)
      {
      case Triangulation:
        intersector=new TriangulationIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>(myMesh_P,myMesh_S, _dimCaracteristic,_precision, 
                                                                                              medianPlane(), _printLevel);
        break;
      case Convex:
        intersector=new ConvexIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>(myMesh_P,myMesh_S, _dimCaracteristic, _precision,
                                                                                       medianPlane(), doRotate(), _printLevel);
        break;
      case Geometric2D:
        intersector=new Geometric2DIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>(myMesh_P, myMesh_S, _dimCaracteristic, _precision);
        break;
        // case MEDMEM::Generic:
        //intersector=new GenericIntersector<SPACEDIM>(myMesh_P,myMesh_S, _DimCaracteristic,_Precision,
        //                                         0, 0, _PrintLevel);
        //break;
      }

    /****************************************************************/
    /* Create a search tree based on the bounding boxes             */
    /* Instanciate the intersector and initialise the result vector */
    /****************************************************************/
 
    long start_filtering=clock();
 
    std::vector<double> bbox;
    intersector->createBoundingBoxes(myMesh_S,bbox); // create the bounding boxes
    performAdjustmentOfBB(intersector,bbox);
    BBTree<SPACEDIM> my_tree(&bbox[0], 0, 0,nbMaille_S);//creating the search structure 

    long end_filtering=clock();

    result.resize(nbMaille_P);//on initialise.

    /****************************************************/
    /* Loop on the target cells - core of the algorithm */
    /****************************************************/
    int i_P=0;//global index of cell

    long start_intersection=clock();
    const ConnType *connIndxP=myMesh_P.getConnectivityIndexPtr();
    const ConnType *connIndxS=myMesh_S.getConnectivityIndexPtr();
    long nbelem_type=myMesh_P.getNumberOfElements();
    for(i_P=0; i_P<nbelem_type; i_P++)
      {
        int nb_nodesP=connIndxP[i_P+1]-connIndxP[i_P];
        std::vector<int> intersecting_elems;
        double bb[2*SPACEDIM];
        intersector->getElemBB(bb,myMesh_P,OTT<ConnType,numPol>::indFC(i_P),nb_nodesP);
        my_tree.getIntersectingElems(bb, intersecting_elems);
        int nb_intersecting_elems = intersecting_elems.size();           
        //browsing all the i_S (from mesh S) elems that can 
        //intersect elem i_P (from mesh P)
        for(int ielem=0; ielem<nb_intersecting_elems;ielem++)
          {
            //BBTree structure returns numbers between 0 and n-1
            int i_S=intersecting_elems[ielem]; //MN: Global number of cell ?
            int nb_nodesS=connIndxS[i_S+1]-connIndxS[i_S];
            double surf=intersector->intersectCells(OTT<ConnType,numPol>::indFC(i_P),OTT<ConnType,numPol>::indFC(i_S),nb_nodesP,nb_nodesS);
            result[i_P].insert(std::make_pair(OTT<ConnType,numPol>::indFC(i_S),surf));
            counter++;
            //rempli_vect_d_intersect(i_P+1,i_S+1,MaP,MaS,result_areas);
          }
        intersecting_elems.clear();
      }
    delete intersector;

    /***********************************/
    /*        DEBUG prints             */
    /***********************************/

    if (_printLevel >=1)
      {
        long end_intersection=clock();
//         if (_printLevel >=2)
//           {
//             std::cout << std::endl << "Printing intersection areas:" << std::endl << std::endl;
//             std::cout << "(source cell, target cell): intersection areas" << std::endl;
//             double total=0.0;
//             double total_interm=0.0;
//             int nb_result_areas = result.size();
//             for(int i=0; i< nb_result_areas;i++)
//               { 
//                 std::map<int,double>::iterator surface;
//                 total_interm=0.0;
//                 for( surface=result[i].begin();surface!=result[i].end();surface++)
//                   {
//                     std::cout<< "    ("<<i+1<<" , " << (*surface).first<<")" << " : " << (*surface).second << std::endl;
//                     total_interm +=(*surface).second;
//                   }
//                 std::cout<< " elem " << i+1 << " area= " << total_interm << std::endl;
//                 total+=total_interm;
//               }
//             std::cout << "total area " << total << std::endl;
//           }
        std::cout << "Filtering time= " << end_filtering-start_filtering << std::endl;
        std::cout << "Intersection time= " << end_intersection-start_intersection << std::endl;
        long global_end =clock();    
        std::cout << "Number of computed intersections = " << counter << std::endl;
        std::cout << "Global time= " << global_end - global_start << std::endl;
      }
  }
}

#endif
