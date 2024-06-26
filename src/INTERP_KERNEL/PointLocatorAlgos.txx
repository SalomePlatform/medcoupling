// Copyright (C) 2007-2024  CEA, EDF
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
#ifndef __POINTLOCATORALGOS_TXX__
#define __POINTLOCATORALGOS_TXX__

#include "InterpolationUtils.hxx"
#include "CellModel.hxx"
#include "BBTree.txx"
#include "VectorUtils.hxx"

#include "InterpKernelGeo2DNode.hxx"
#include "InterpKernelGeo2DQuadraticPolygon.hxx"

#include <list>
#include <vector>
#include <set>
#include <limits>
#include <memory>

namespace INTERP_KERNEL
{
  class GenericPointLocatorAlgos
  {
  public:
    virtual ~GenericPointLocatorAlgos() { }
    virtual std::list<mcIdType> locates(const double* x, double eps) = 0;
  };
        
  template<class MyMeshType>
  class PointLocatorAlgos: public GenericPointLocatorAlgos
  {
  private : 
    double* _bb;
    BBTree<MyMeshType::MY_SPACEDIM,typename MyMeshType::MyConnType>* _tree;
    const MyMeshType& _mesh;
  public:
    PointLocatorAlgos(const MyMeshType& mesh):_mesh(mesh)
    {
      typedef typename MyMeshType::MyConnType ConnType;
      const int SPACEDIM=MyMeshType::MY_SPACEDIM;
      const NumberingPolicy numPol=MyMeshType::My_numPol;
      ConnType nelem = _mesh.getNumberOfElements();
      _bb = new double[SPACEDIM*2*nelem];
      const ConnType* conn = _mesh.getConnectivityPtr();
      const ConnType* conn_index = _mesh.getConnectivityIndexPtr();
      const double* coords=_mesh.getCoordinatesPtr();
      for (ConnType i=0; i<nelem; i++)
        {
          for (int idim=0; idim<SPACEDIM; idim++)
            {
              _bb[2*(i*SPACEDIM+idim)]=std::numeric_limits<double>::max();
              _bb[2*(i*SPACEDIM+idim)+1]=-std::numeric_limits<double>::max();
            }
          for (ConnType index= conn_index[i]; index < conn_index[i+1];index++)
            {
              //coordelem points to the coordinates of the current node of the i-th element
              const double* coordelem = coords+OTT<ConnType,numPol>::ind2C(conn[OTT<ConnType,numPol>::ind2C(index)])*SPACEDIM;

              //the bounding box is updated by checking wheher the node is at the min/max in exach dimension
              for (int idim=0; idim<SPACEDIM;idim++)
                {
                  _bb[2*(i*SPACEDIM+idim)]=(coordelem[idim]<_bb[2*(i*SPACEDIM+idim)])?coordelem[idim]:_bb[2*(i*SPACEDIM+idim)];
                  _bb[2*(i*SPACEDIM+idim)+1]=(coordelem[idim]>_bb[2*(i*SPACEDIM+idim)+1])?coordelem[idim]:_bb[2*(i*SPACEDIM+idim)+1];
                }
            }
        }
      _tree=new BBTree<SPACEDIM,typename MyMeshType::MyConnType>(_bb,0,0,nelem);
    }

    ~PointLocatorAlgos()
    {
      delete[] _bb;
      delete _tree;
    }
        
    //returns the list of elements that contains 
    //the point pointed to by x
    std::list<typename MyMeshType::MyConnType> locates(const double* x, double eps)
    {
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;
      std::vector<ConnType> candidates;
      _tree->getElementsAroundPoint(x,candidates);
      std::list<ConnType> retlist;
      for(unsigned int i=0; i< candidates.size(); i++)
        {
          ConnType ielem=candidates[i];
          if (elementContainsPoint(ielem,x,eps))
            retlist.push_back(OTT<ConnType,numPol>::indFC(ielem));
        }
      return retlist;
    }

    static bool isElementContainsPointAlg2DSimple(const double *ptToTest, const double *cellPts, mcIdType nbEdges, double eps)
    {
      /* with dimension 2, it suffices to check all the edges
         and see if the sign of double products from the point
         is always the same.
                         C
                        / \
                       /   \
             Xo       /     \ 
                     A-------B
       
         here XA^XC and XC^XB have different signs*/
      const int SPACEDIM=MyMeshType::MY_SPACEDIM;
      std::unique_ptr<char[]> sign( new char[nbEdges] );
      for (mcIdType iedge=0; iedge<nbEdges; iedge++)
        {
          const double* A=cellPts+SPACEDIM*iedge;
          const double* B=cellPts+SPACEDIM*((iedge+1)%nbEdges);
          double a=mon_determinant(ptToTest, A, B);
          if(a<-eps)
            sign[iedge]=-1;
          else if(a>eps)
            sign[iedge]=1;
          else
            sign[iedge]=0;
        }
      bool ret=decideFromSign(sign.get(), nbEdges);
      return ret;
    }

    // Same as isElementContainsPointAlg2DSimple() with a different input format ...
    static bool isElementContainsPointAlgo2DSimple2(const double *ptToTest, NormalizedCellType type,
                                                    const double *coords, const typename MyMeshType::MyConnType *conn_elem,
                                                    typename MyMeshType::MyConnType conn_elem_sz, double eps)
    {
      const int SPACEDIM=MyMeshType::MY_SPACEDIM;
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;

      const CellModel& cmType=CellModel::GetCellModel(type);
      bool ret=false;

      int nbEdges=cmType.getNumberOfSons();
      double *pts = new double[nbEdges*SPACEDIM];
      for (int iedge=0; iedge<nbEdges; iedge++)
        {
          const double* a=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[iedge]));
          std::copy(a,a+SPACEDIM,pts+iedge*SPACEDIM);
        }
      ret=isElementContainsPointAlg2DSimple(ptToTest,pts,nbEdges,eps);
      delete [] pts;
      return ret;
    }

    static bool isElementContainsPointAlgo2DPolygon(const double *ptToTest, NormalizedCellType type,
                                                    const double *coords, const typename MyMeshType::MyConnType *conn_elem,
                                                    typename MyMeshType::MyConnType conn_elem_sz, double eps)
    {
      // Override precision for this method only:
      INTERP_KERNEL::QuadraticPlanarPrecision prec(eps);

      const int SPACEDIM=MyMeshType::MY_SPACEDIM;
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;

      std::vector<INTERP_KERNEL::Node *> nodes(conn_elem_sz);
      INTERP_KERNEL::QuadraticPolygon *pol(0);
      for(mcIdType j=0;j<conn_elem_sz;j++)
        {
          mcIdType nodeId(OTT<ConnType,numPol>::ind2C(conn_elem[j]));
          nodes[j]=new INTERP_KERNEL::Node(coords[nodeId*SPACEDIM],coords[nodeId*SPACEDIM+1]);
        }
      if(!INTERP_KERNEL::CellModel::GetCellModel(type).isQuadratic())
        pol=INTERP_KERNEL::QuadraticPolygon::BuildLinearPolygon(nodes);
      else
        pol=INTERP_KERNEL::QuadraticPolygon::BuildArcCirclePolygon(nodes);
      INTERP_KERNEL::Node *n(new INTERP_KERNEL::Node(ptToTest[0],ptToTest[1]));
      double a(0.),b(0.),c(0.);
      a=pol->normalizeMe(b,c); n->applySimilarity(b,c,a);
      bool ret=pol->isInOrOut2(n);
      delete pol; n->decrRef();
      return ret;
    }

    static bool isElementContainsPointAlg3D(const double *ptToTest, const typename MyMeshType::MyConnType *conn_elem, typename MyMeshType::MyConnType conn_elem_sz, const double *coords, const CellModel& cmType, double eps)
    {
      const int SPACEDIM=MyMeshType::MY_SPACEDIM;
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;
      
      int nbfaces = cmType.getNumberOfSons2(conn_elem,conn_elem_sz);
      std::unique_ptr<char[]> sign( new char[nbfaces] );
      std::unique_ptr<ConnType[]> connOfSon( new ConnType[conn_elem_sz] );
      std::unique_ptr<double[]> ptsOfTetrahedrizedPolyhedron( new double[3*3*nbfaces] );
      for (int iface=0; iface<nbfaces; iface++)
        {
          NormalizedCellType typeOfSon;
          unsigned nb_face_nodes = cmType.fillSonCellNodalConnectivity2(iface,conn_elem,conn_elem_sz,connOfSon.get(),typeOfSon);
          // in case of polyhedron, points can be colinear on some faces
          // => loop until non colinear points are found
          // or element type is not polyhedron
          // or max number of nodes is reached to avoid infinite loop on this face
          bool non_colinear_points_found = false;
          unsigned i_permutation = 0;
          while (!non_colinear_points_found)
          {
            const double* AA=coords+SPACEDIM*(OTT<ConnType,numPol>::coo2C(connOfSon[(0+i_permutation) % nb_face_nodes]));
            const double* BB=coords+SPACEDIM*(OTT<ConnType,numPol>::coo2C(connOfSon[(1+i_permutation) % nb_face_nodes]));
            const double* CC=coords+SPACEDIM*(OTT<ConnType,numPol>::coo2C(connOfSon[(2+i_permutation) % nb_face_nodes]));
            if (cmType.getEnum() == NORM_POLYHED && isColinear3DPts(AA, BB, CC, eps) && i_permutation < nb_face_nodes)
              // points are colinear => use next 3 other points
              i_permutation+=1;
            else
              {
              std::copy(AA,AA+3,ptsOfTetrahedrizedPolyhedron.get()+9*iface);
              std::copy(BB,BB+3,ptsOfTetrahedrizedPolyhedron.get()+9*iface+3);
              std::copy(CC,CC+3,ptsOfTetrahedrizedPolyhedron.get()+9*iface+6);
              non_colinear_points_found = true;
            }
          }
        }
      double centerPt[3];
      double normalizeFact = NormalizeTetrahedrizedPolyhedron(ptsOfTetrahedrizedPolyhedron.get(),nbfaces,centerPt);
      double ptToTestNorm[3] = {(ptToTest[0]-centerPt[0])/normalizeFact,
      (ptToTest[1]-centerPt[1])/normalizeFact,
      (ptToTest[2]-centerPt[2])/normalizeFact};
      for (int iface=0; iface<nbfaces; iface++)
        {
          double lengthNorm = TripleProduct(ptsOfTetrahedrizedPolyhedron.get() + 9*iface + 3,
                                            ptsOfTetrahedrizedPolyhedron.get() + 9*iface + 6,
                                            ptToTestNorm,
                                            ptsOfTetrahedrizedPolyhedron.get() + 9*iface);
          // assigne sign[iface] : 0 means on the face within eps, 1 or -1 means far from tetra representing face
          if( lengthNorm<-eps )
            sign[iface]=-1;
          else if( lengthNorm>eps )
            sign[iface]=1;
          else
            sign[iface]=0;
        }
      bool ret=decideFromSign(sign.get(), nbfaces);
      return ret;
    }

    /*!
     * Precondition : spacedim==meshdim. To be checked upstream to this call.
     */
    static bool isElementContainsPoint(const double *ptToTest, NormalizedCellType type, const double *coords,
                                       const typename MyMeshType::MyConnType *conn_elem,
                                       typename MyMeshType::MyConnType conn_elem_sz, double eps)
    {
      const int SPACEDIM=MyMeshType::MY_SPACEDIM;
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;

      const CellModel& cmType=CellModel::GetCellModel(type);
      //
      if (SPACEDIM==2)
        {
          if(type != INTERP_KERNEL::NORM_POLYGON && !cmType.isQuadratic())
            return isElementContainsPointAlgo2DSimple2(ptToTest, type, coords, conn_elem, conn_elem_sz, eps);
          else
            return isElementContainsPointAlgo2DPolygon(ptToTest, type, coords, conn_elem, conn_elem_sz, eps);
        }
      if (SPACEDIM==3)
        {
          return isElementContainsPointAlg3D(ptToTest,conn_elem,conn_elem_sz,coords,cmType,eps);
        }

      if(SPACEDIM==1)
        {
          double p1=coords[(OTT<ConnType,numPol>::ind2C(conn_elem[0]))];
          double p2=coords[(OTT<ConnType,numPol>::ind2C(conn_elem[1]))];
          double delta=fabs(p1-p2)+eps;
          double val=*ptToTest-std::min(p1,p2);
          return val>-eps && val<delta;
        }
      throw INTERP_KERNEL::Exception("Invalid spacedim detected ! Managed spaceDim are 2 and 3 !");
    }
        
    bool elementContainsPoint(typename MyMeshType::MyConnType i, const double* x, double eps)
    {
      //as i is extracted from the BBTRee, it is already in C numbering
      //it is not necessary to convert it from F to C
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;
      
      const double* coords= _mesh.getCoordinatesPtr();
      const ConnType* conn=_mesh.getConnectivityPtr();
      const ConnType* conn_index= _mesh.getConnectivityIndexPtr();
      const ConnType* conn_elem=conn+OTT<ConnType,numPol>::ind2C(conn_index[i]);
      int conn_elem_sz=conn_index[i+1]-conn_index[i];
      NormalizedCellType type=_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(i));
      return isElementContainsPoint(x,type,coords,conn_elem,conn_elem_sz,eps);
    }
                
    static bool decideFromSign(const char *sign, mcIdType nbelem)
    {
      char min_sign = 1;
      char max_sign = -1;
      for (int i=0; i<nbelem;i++)
        {
          min_sign=(sign[i]<min_sign)?sign[i]:min_sign;
          max_sign=(sign[i]>max_sign)?sign[i]:max_sign;
        }
      return (min_sign!=-1 || max_sign!=1);     
    }
  };

  template<class MyMeshType>
  class PointLocatorInSimplex : public PointLocatorAlgos<MyMeshType>
  {
    const MyMeshType& _mesh;
  public:
    PointLocatorInSimplex(const MyMeshType& mesh)
      :PointLocatorAlgos<MyMeshType>(mesh),_mesh(mesh)
    {
    }

    //================================================================================
    /*!
     * \brief Returns nodes composing the simplex the point x is in
     */
    //================================================================================

    virtual std::list<typename MyMeshType::MyConnType> locates(const double* x, double eps)
    {
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;

      std::list<ConnType> simplexNodes;
      std::list<ConnType> candidates = PointLocatorAlgos<MyMeshType>::locates(x,eps);
      typename std::list<ConnType>::iterator eIt = candidates.begin();
      for ( ; eIt != candidates.end(); ++eIt )
        {
          const ConnType i = OTT<ConnType,numPol>::ind2C( *eIt );
          const double* coords= _mesh.getCoordinatesPtr();
          const ConnType* conn=_mesh.getConnectivityPtr();
          const ConnType* conn_index= _mesh.getConnectivityIndexPtr();
          const ConnType* conn_elem=conn+OTT<ConnType,numPol>::ind2C(conn_index[i]);
          int conn_elem_sz=conn_index[i+1]-conn_index[i];
          NormalizedCellType type=_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(i));
          CellModel cell = CellModel::GetCellModel(type);

          if ( cell.isQuadratic() )
            throw Exception("P2 not implemented yet");

          if ( cell.isSimplex())
            {
              for ( int n = 0; n < conn_elem_sz; ++n )
                simplexNodes.push_back( conn_elem[ n ]);
            }
          else
            {
              NormalizedCellType simlexType = cell.getDimension()==3 ? NORM_TETRA4 : NORM_TRI3;
              std::vector<mcIdType> sonNodes;
              NormalizedCellType sonType;
              const unsigned nbSons = cell.getNumberOfSons2( conn_elem, conn_elem_sz );
              for ( unsigned s = 0; s < nbSons; ++s )
                {
                  sonNodes.resize( cell.getNumberOfNodesConstituentTheSon2( s, conn_elem, conn_elem_sz ));
                  cell.fillSonCellNodalConnectivity2( s, conn_elem, conn_elem_sz, &sonNodes[0], sonType );
                  std::set<mcIdType> sonNodesSet( sonNodes.begin(), sonNodes.end() );

                  std::set< std::set< ConnType > > checkedSonSimplex;
                  for ( unsigned sn = 0; sn < sonNodes.size(); ++sn )
                    {
                      std::vector< ConnType > simplexConn( cell.getDimension() + 1 );
                      unsigned n;
                      for ( n = 0; n < cell.getDimension()-1; ++n )
                        simplexConn[n] = sonNodes[ (sn+n) % sonNodes.size() ];

                      for ( unsigned n2 = 0; n2 < sonNodes.size()-cell.getDimension()+1; ++n2 )
                        {
                          simplexConn[n] = sonNodes[ (sn+n+n2) % sonNodes.size() ];
                          std::set< ConnType > sonSimplex( simplexConn.begin(), --simplexConn.end());
                          if ( checkedSonSimplex.insert( sonSimplex ).second )
                            {
                              for ( unsigned cn = 0; cn < conn_elem_sz; ++cn )
                                if ( !sonNodesSet.count( conn_elem[cn] ))
                                  {
                                    simplexConn.back() = conn_elem[cn];
                                    if ( this->isElementContainsPoint( x, simlexType, coords,
                                                                       &simplexConn[0], simplexConn.size(), eps ))
                                      {
                                        simplexNodes.insert( simplexNodes.end(),
                                                             simplexConn.begin(), simplexConn.end());
                                        return simplexNodes;
                                      }
                                  }
                            }
                        }
                    }
                }
            }
        }
      return simplexNodes;
    }

  };
}

#endif
