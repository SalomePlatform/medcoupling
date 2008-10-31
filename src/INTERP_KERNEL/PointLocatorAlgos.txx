#ifndef _POINT_LOCATOR_ALGOS_TXX_
#define _POINT_LOCATOR_ALGOS_TXX_

#include "InterpolationUtils.hxx"
#include "CellModel.hxx"
#include "BBTree.txx"
#include <list>
#include <limits>

namespace INTERP_KERNEL{

  class GenericPointLocatorAlgos
  {
  public:
    virtual ~GenericPointLocatorAlgos(){};
    virtual std::list<int> locates(const double* x)=0;
         
  };
        
  template<class MyMeshType>
  class PointLocatorAlgos: public GenericPointLocatorAlgos
  {
  private : 
    double* _bb;
    BBTree<MyMeshType::MY_SPACEDIM>* _tree;
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
                        
      _tree=new BBTree<SPACEDIM>(_bb,0,0,nelem);
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
        
    bool elementContainsPoint(typename MyMeshType::MyConnType i, const double*x)
    {
      //as i is extracted from the BBTRee, it is already in C numbering
      //it is not necessary to convert it from F to C
      typedef typename MyMeshType::MyConnType ConnType;
      const int SPACEDIM=MyMeshType::MY_SPACEDIM;
      const NumberingPolicy numPol=MyMeshType::My_numPol;
      
      const double* coords= _mesh.getCoordinatesPtr();
      const ConnType* conn=_mesh.getConnectivityPtr();
      const ConnType* conn_index= _mesh.getConnectivityIndexPtr();
      const ConnType* conn_elem=conn+OTT<ConnType,numPol>::ind2C(conn_index[i]);
      NormalizedCellType type=_mesh.getTypeOfElement(OTT<ConnType,numPol>::indFC(i));
      const CellModel& cmType=CellModel::getCellModel(type);

      int nbnodes = cmType.getNumberOfNodes();//conn_index[i+1]-conn_index[i];
                
      // with dimension 2, it suffices to check all the edges
      // and see if the sign of double products from the point
      //is always the same.
      //                 C
      //                / \
      //               /   \
      //     Xo       /     \ 
      //             A-------B
      //
      //here XA^XC and XC^XB have different signs
      //
      if (SPACEDIM==2)
        {
          //in 2D, nbedges==nbnodes
          int nbedges=nbnodes;
          int* sign = new int[nbedges];
          for (int iedge=0; iedge<nbedges; iedge++)
            {
              const double* A=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[iedge]));
              const double* B;
              if (iedge+1< nbedges)
                B=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[iedge+1]));
              else
                B=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[0]));
                                                        
              double a=mon_determinant(x, A, B);
              if (a<-1e-12)
                sign[iedge]=-1;
              else if (a>1e-12)
                sign[iedge]=1;
              else
                sign[iedge]=0;
            }
          return decide_from_sign(sign, nbedges);
        }
                        
      if (SPACEDIM==3)
        {
          int nbfaces = cmType.getNumberOfSons();
          int* sign = new int[nbfaces];
          for (int iface=0; iface<nbfaces; iface++)
            {
              const unsigned* connface=cmType.getNodesConstituentTheSon(iface);
              const double* AA=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[connface[0]]));
              const double* BB=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[connface[1]]));
              const double* CC=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[connface[2]]));
                                                        
              double Vol=triple_product(AA,BB,CC,x);
              if (Vol<-1e-12)
                sign[iface]=-1;
              else if (Vol>1e-12)
                sign[iface]=1;
              else
                sign[iface]=0;
            }
          return decide_from_sign(sign, nbfaces);
        }
                        
    }
                
    bool decide_from_sign (const int* sign, int nbelem)
    {
      int min_sign =1;
      int max_sign =-1;
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
