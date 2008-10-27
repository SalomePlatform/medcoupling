#ifndef _POINT_LOCATOR_HXX_
#define _POINT_LOCATOR_HXX_

#include <list>
#include "MEDNormalizedUnstructuredMesh.hxx"
#include "PointLocatorAlgos.txx"

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT PointLocator
  {
  public:
    PointLocator(const MEDMEM::MESH& mesh);
    virtual ~PointLocator();
    std::list<int> locate(const double* x);
  private:
    GenericMesh* _medmesh;
    GenericPointLocatorAlgos* _point_locator;
  };
}

#endif
