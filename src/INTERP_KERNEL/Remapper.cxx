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
#include "Remapper.hxx"
#include "MEDMEM_Exception.hxx"
#include "Interpolation.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3D.txx"
#include "Interpolation3DSurf.txx"
#include "MEDNormalizedUnstructuredMesh.txx"

namespace INTERP_KERNEL
{
	
// 	int InterpolationOptions::_printLevel=0;
// 	IntersectionType  InterpolationOptions::_intersectionType=Triangulation;
// 	double  InterpolationOptions::_precision=1e-12;;
// 	double  InterpolationOptions::_medianPlane =0.5;
// 	bool  InterpolationOptions::_doRotate =true;
// 	double  InterpolationOptions::_boundingBoxAdjustment =0.1;
// 	int  InterpolationOptions::_orientation =0;
// 	SplittingPolicy  InterpolationOptions::_splittingPolicy =GENERAL_48;

  Remapper::Remapper():_matrix(0)
  {
  }

  Remapper::~Remapper()
  {
    delete _matrix;
  }
  /*! This method computes the intersection matrix between 
   * source \a mesh_source and \a mesh_target. It is a preliminary step 
   * that is necessary before calling the \c transfer() method.
   * The method analyses the dimensions of the meshes and checks for compatibility.
   * 
   */
  void Remapper::prepare(const MEDMEM::MESH& mesh_source, const MEDMEM::MESH& mesh_target)
  {
    const int sm_spacedim = mesh_source.getSpaceDimension();
    const int tm_spacedim = mesh_target.getSpaceDimension();
    int sm_meshdim = mesh_source.getMeshDimension();
    int tm_meshdim = mesh_target.getMeshDimension();

    if (tm_spacedim!=sm_spacedim || tm_meshdim!=sm_meshdim)
      throw MEDEXCEPTION("incompatible mesh and/or space dimensions in meshes");
    _matrix= new Matrix<double,ALL_FORTRAN_MODE>;
    if ((sm_spacedim==2)&&(sm_meshdim==2))
      {
        MEDNormalizedUnstructuredMesh<2,2> source_mesh_wrapper(&mesh_source);
        MEDNormalizedUnstructuredMesh<2,2> target_mesh_wrapper(&mesh_target);
        Interpolation2D interpolation;
        interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,*_matrix);
      }
    else if ((sm_spacedim==3)&&(sm_meshdim==3))
      {
        MEDNormalizedUnstructuredMesh<3,3> source_mesh_wrapper(&mesh_source);
        MEDNormalizedUnstructuredMesh<3,3> target_mesh_wrapper(&mesh_target);
        Interpolation3D interpolation;
        interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,*_matrix);
      }
    else if ((sm_spacedim==3)&&(sm_meshdim==2))
      {
        MEDNormalizedUnstructuredMesh<3,2> source_mesh_wrapper(&mesh_source);
        MEDNormalizedUnstructuredMesh<3,2> target_mesh_wrapper(&mesh_target);
        Interpolation3DSurf interpolation;
        interpolation.interpolateMeshes(source_mesh_wrapper,target_mesh_wrapper,*_matrix);
      }
    else
      throw MEDEXCEPTION("no Interpolation exists for the given mesh and space dimensions");
  }

  void Remapper::transfer(const MEDMEM::FIELD<double>& field_source,MEDMEM::FIELD<double>& field_target)
  {
    int source_nbcomp=field_source.getNumberOfComponents();
    int target_nbcomp=field_target.getNumberOfComponents();
    if (source_nbcomp != target_nbcomp)
      throw MEDMEM::MEDEXCEPTION("incoherent number of components for source and target fields");
    if (source_nbcomp>1)
      throw MEDMEM::MEDEXCEPTION("interpolations with more than one component are not yet handled");
    MEDMEM::FIELD<double>* target_volumes = getSupportVolumes(*field_target.getSupport());
    int nbelem_target=field_target.getSupport()->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS);

    double* value_target = const_cast<double*> (field_target.getValue());
    const double* value_source = field_source.getValue();

    _matrix->multiply(value_source, value_target);
    for (int i=0; i< nbelem_target; i++)
      value_target[i]/=target_volumes->getValueIJ(i+1,1);

    delete target_volumes;
  }

  void Remapper::setOptionDouble(const std::string& key, double value)
  {
  }

  void Remapper::setOptionInt(const std::string& key, int value)
  {
  }

  /*!
    \brief returns the volumes of the cells underlying the field \a field

    For 2D geometries, the returned field contains the areas.
    For 3D geometries, the returned field contains the volumes.

    \param field field on which cells the volumes are required
    \return field containing the volumes
  */
  MEDMEM::FIELD<double>* Remapper::getSupportVolumes(const MEDMEM::SUPPORT& support)
  {
    const MEDMEM::MESH* mesh=support.getMesh();
    int dim = mesh->getMeshDimension();
    switch (dim)
      {
      case 2:
        return mesh->getArea(&support);
      case 3:
        return mesh->getVolume(&support);
      default:
        throw MEDMEM::MEDEXCEPTION("interpolation is not available for this dimension");
      }
  }


}
