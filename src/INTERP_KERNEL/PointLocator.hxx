#ifndef _POINT_LOCATOR_HXX_
#define _POINT_LOCATOR_HXX_
#include <list>
#include "NormalizedUnstructuredMesh.hxx"
#include "MEDNormalizedUnstructuredMesh.hxx"
#include "PointLocatorAlgos.txx"
namespace INTERP_KERNEL
{
class INTERPKERNEL_EXPORT PointLocator
{
public:
	PointLocator(const MEDMEM::MESH& mesh);
	virtual ~PointLocator();
	std::list<int> locate(const double*x);
	
private:
	//MEDNormalizedUnstructuredMesh<2,2>* _medmesh;
	//GenericPointLocatorAlgos*< 2,2,int,ALL_FORTRAN_MODE,MEDNormalizedUnstructuredMesh<2,2> >* _point_locator;

	GenericMesh* _medmesh;
	GenericPointLocatorAlgos* _point_locator;
};
}
#endif
