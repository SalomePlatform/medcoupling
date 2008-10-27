#include <list>
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Exception.hxx"
#include "PointLocatorAlgos.txx"
#include "PointLocator.hxx"

namespace INTERP_KERNEL {
  PointLocator::PointLocator(const MEDMEM::MESH& mesh)
  {
    int meshdim=mesh.getMeshDimension();
    int spacedim=mesh.getSpaceDimension();
    if (meshdim != spacedim) throw MEDMEM::MEDEXCEPTION("Locator is not implemented for meshdim != spacedim");
    switch (meshdim)
      {
      case 2:
        _medmesh = new MEDNormalizedUnstructuredMesh<2,2> (&mesh);
        _point_locator=new PointLocatorAlgos<MEDNormalizedUnstructuredMesh<2,2> >(*(static_cast<MEDNormalizedUnstructuredMesh<2,2>* >(_medmesh)));
        break;
      case 3:
        _medmesh = new MEDNormalizedUnstructuredMesh<3,3> (&mesh);
        _point_locator=new PointLocatorAlgos<MEDNormalizedUnstructuredMesh<3,3> >(*(static_cast<MEDNormalizedUnstructuredMesh<3,3>* >(_medmesh)));
        break;
      }
  }
  PointLocator::~PointLocator()
  {
    delete _medmesh;
    delete _point_locator;
  }

std::list<int> PointLocator::locate(const double* x)
  {
    return _point_locator->locates(x);
  }
}
