#ifndef _POINT_LOCATOR_ALGOS_TXX_
#define _POINT_LOCATOR_ALGOS_TXX_

#include "MEDMEM_Exception.hxx"
#include "InterpolationUtils.hxx"
#include "MEDMEM_CellModel.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "BBTree.txx"
#include <list>
namespace INTERP_KERNEL{

	class GenericPointLocatorAlgos
	{
	public:
		virtual ~GenericPointLocatorAlgos(){};
		virtual std::list<int> locates(const double* x)=0;
	 
	};
	
	template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType> class PointLocatorAlgos: public GenericPointLocatorAlgos
	{
	private : 
		double* _bb;
		BBTree<SPACEDIM>* _tree;
		const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& _mesh;

	public:
		PointLocatorAlgos<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh):_mesh(mesh)
		{
			int nelem = _mesh.getNumberOfElements();
			_bb = new double[SPACEDIM*2*nelem];
			const ConnType* conn = _mesh.getConnectivityPtr();
			const ConnType* conn_index = _mesh.getConnectivityIndexPtr();
			const double* coords=_mesh.getCoordinatesPtr();
			for (int i=0; i<nelem; i++)
				{
					for (int idim=0; idim<SPACEDIM; idim++)
						{
							_bb[2*(i*SPACEDIM+idim)]=HUGE;
							_bb[2*(i*SPACEDIM+idim)+1]=-HUGE;
						}
					for (int index= conn_index[i]; index < conn_index[i+1];index++)
						{
							const double* coordelem = coords+OTT<ConnType,numPol>::ind2C(conn[OTT<ConnType,numPol>::ind2C(index)]);
							for (int idim=0; idim<SPACEDIM;idim++)
								{
									_bb[2*(i*SPACEDIM+idim)]=(coordelem[idim]<_bb[2*(i*SPACEDIM+idim)])?coordelem[idim]:_bb[2*(i*SPACEDIM+idim)];
									_bb[2*(i*SPACEDIM+idim)+1]=(coordelem[idim]>_bb[2*(i*SPACEDIM+idim)+1])?coordelem[idim]:_bb[2*(i*SPACEDIM+idim)+1];
								}
						}
				}
			
			_tree=new BBTree<SPACEDIM>(_bb,0,0,nelem);
		}
		~PointLocatorAlgos<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>()
		{
			delete[] _bb;
			delete _tree;
		}
	
		std::list<int> locates(const double* x)
		{
			vector<int> candidates;
			_tree->getElementsAroundPoint(x,candidates);
			list<int> retlist;
			for (int i=0; i< candidates.size(); i++)
				{
					if (elementContainsPoint(i,x))
						retlist.push_back(OTT<ConnType,numPol>::indFC(i));
				}
			return retlist;
		}
	
		bool elementContainsPoint(int i, const double*x)
		{
			//as i is extracted from the BBTRee, it is already in C numbering
			//it is not necessary to convert it from F to C

			const double* coords= _mesh.getCoordinatesPtr();
			const int* conn=_mesh.getConnectivityPtr();
			const int* conn_index= _mesh.getConnectivityIndexPtr();
			const int* conn_elem=conn+OTT<ConnType,numPol>::ind2C(conn_index[i]);

		
			int nbnodes = conn_index[i+1]-conn_index[i];
		
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
					int sign[nbedges];
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
					MED_EN::medGeometryElement elem_type;
					switch (nbnodes) {
					case 4 :
						elem_type=MED_EN::MED_TETRA4;
						break;
					case 5:
						elem_type=MED_EN::MED_PYRA5;
						break;
					case 6 :
						elem_type=MED_EN::MED_PENTA6;
						break;
					case 8: 
						elem_type=MED_EN::MED_HEXA8;
						break;
					default:
						throw MEDMEM::MEDEXCEPTION("PointLocatorAlgos : bad number of nodes in 3D locator");
					}
					const MEDMEM::CELLMODEL& model=MEDMEM::CELLMODEL_Map::retrieveCellModel(elem_type);
					int nbfaces = model.getNumberOfConstituents(1);
					int sign[nbfaces];
					for (int iface=0; iface<nbfaces; iface++)
						{
							int* connface=model.getNodesConstituent(1,iface+1);
							const double* AA=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[connface[0]-1]));
							const double* BB=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[connface[1]-1]));
							const double* CC=coords+SPACEDIM*(OTT<ConnType,numPol>::ind2C(conn_elem[connface[2]-1]));
							
							double Vol=triple_product(AA,BB,CC,x);
							if (Vol<-1e-12)
								sign[iface]=-1;
							else if (Vol>1e-12)
								sign[iface]=1;
							else
								sign[iface]=0;
						}
					return  decide_from_sign(sign, nbfaces);
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
			if (min_sign==-1 && max_sign==1)
				return false;
			else
				return true;	
		}
	};
}
#endif
