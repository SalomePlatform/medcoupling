//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef __POINTLOCATORALGOS_TXX__
#define __POINTLOCATORALGOS_TXX__

#include "InterpolationUtils.hxx"
#include "CellModel.hxx"
#include "BBTree.txx"

#include <list>
#include <limits>

namespace INTERP_KERNEL
{
  class GenericPointLocatorAlgos
  {
  public:
    virtual ~GenericPointLocatorAlgos() { }
    virtual std::list<int> locates(const double* x) = 0;     
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
      int nelem = _mesh.getNumberOfElements();
      _bb = new double[SPACEDIM*2*nelem];
      const ConnType* conn = _mesh.getConnectivityPtr();
      const ConnType* conn_index = _mesh.getConnectivityIndexPtr();
      const double* coords=_mesh.getCoordinatesPtr();
      for (int i=0; i<nelem; i++)
        {
          for (int idim=0; idim<SPACEDIM; idim++)
            {
              _bb[2*(i*SPACEDIM+idim)]=std::numeric_limits<double>::max();
              _bb[2*(i*SPACEDIM+idim)+1]=-std::numeric_limits<double>::max();
            }
          for (int index= conn_index[i]; index < conn_index[i+1];index++)
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
    std::list<typename MyMeshType::MyConnType> locates(const double* x)
    {
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;
      std::vector<ConnType> candidates;
      _tree->getElementsAroundPoint(x,candidates);
      std::list<ConnType> retlist;
      for(int i=0; i< candidates.size(); i++)
        {
          int ielem=candidates[i];
          if (elementContainsPoint(ielem,x))
            retlist.push_back(OTT<ConnType,numPol>::indFC(ielem));
        }
      return retlist;
    }

    static bool isElementContainsPointAlg2D(const double *ptToTest, const double *cellPts, int nbEdges, double eps)
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
      int* sign = new int[nbEdges];
      for (int iedge=0; iedge<nbEdges; iedge++)
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
      bool ret=decideFromSign(sign, nbEdges);
      delete [] sign;
      return ret;
    }

    static bool isElementContainsPointAlg3D(const double *ptToTest, const int *conn_elem, const double *coords, const CellModel& cmType, double eps)
    {
      const int SPACEDIM=MyMeshType::MY_SPACEDIM;
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;
      
      int nbfaces = cmType.getNumberOfSons();
      int* sign = new int[nbfaces];
      for (int iface=0; iface<nbfaces; iface++)
        {
          const unsigned* connface=cmType.getNodesConstituentTheSon(iface);
          const double* AA=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[connface[0]]));
          const double* BB=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[connface[1]]));
          const double* CC=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[connface[2]]));
                                                        
          double Vol=triple_product(AA,BB,CC,ptToTest);
          if (Vol<-eps)
            sign[iface]=-1;
          else if (Vol>eps)
            sign[iface]=1;
              else
                sign[iface]=0;
        }
      bool ret=decideFromSign(sign, nbfaces);
      delete [] sign;
      return ret;
    }

    static bool isElementContainsPoint(const double *ptToTest, NormalizedCellType type, const double *coords, const int *conn_elem)
    {
      const int SPACEDIM=MyMeshType::MY_SPACEDIM;
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;

      const CellModel& cmType=CellModel::getCellModel(type);
      //
      if (SPACEDIM==2)
        {
          int nbEdges=cmType.getNumberOfSons();
          double *pts = new double[nbEdges*SPACEDIM];
          for (int iedge=0; iedge<nbEdges; iedge++)
            {
              const double* a=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[iedge]));
              std::copy(a,a+SPACEDIM,pts+iedge*SPACEDIM);
            }
          bool ret=isElementContainsPointAlg2D(ptToTest,pts,nbEdges,1e-12);
          delete [] pts;
          return ret;
        }
                        
      if (SPACEDIM==3)
        {
          return isElementContainsPointAlg3D(ptToTest,conn_elem,coords,cmType,1e-12);
        }
    }
        
    bool elementContainsPoint(typename MyMeshType::MyConnType i, const double* x)
    {
      //as i is extracted from the BBTRee, it is already in C numbering
      //it is not necessary to convert it from F to C
      typedef typename MyMeshType::MyConnType ConnType;
      const NumberingPolicy numPol=MyMeshType::My_numPol;
      
      const double* coords= _mesh.getCoordinatesPtr();
      const ConnType* conn=_mesh.getConnectivityPtr();
      const ConnType* conn_index= _mesh.getConnectivityIndexPtr();
      const ConnType* conn_elem=conn+OTT<ConnType,numPol>::ind2C(conn_index[i]);
      NormalizedCellType type=_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(i));
      return isElementContainsPoint(x,type,coords,conn_elem);
    }
                
    static bool decideFromSign(const int* sign, int nbelem)
    {
      int min_sign = 1;
      int max_sign = -1;
      for (int i=0; i<nbelem;i++)
        {
          min_sign=(sign[i]<min_sign)?sign[i]:min_sign;
          max_sign=(sign[i]>max_sign)?sign[i]:max_sign;
        }
      return (min_sign!=-1 || max_sign!=1);     
    }
  };
}

#endif
