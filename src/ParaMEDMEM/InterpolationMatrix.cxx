#include "ParaMESH.hxx"
#include "ProcessorGroup.hxx"
#include "MxN_Mapping.hxx"
#include "InterpolationMatrix.hxx"
#include "TranslationRotationMatrix.hxx"
#include "Interpolation.hxx"
#include "Interpolation2D.txx"
#include "Interpolation3DSurf.txx"
#include "Interpolation3D.txx"
#include "MEDNormalizedUnstructuredMesh.hxx"
#include "InterpolationOptions.hxx"

/*! \class InterpolationMatrix

This class enables the storage of an interpolation matrix Wij mapping 
source field Sj to target field Ti via Ti=Vi^(-1).Wij.Sj.
The matrix is built and stored on the processors belonging to the source
group. 

*/

namespace ParaMEDMEM
{

  /*!
    Creates an empty matrix structure linking two distributed supports.
    The method must be called by all processors belonging to source and target groups.
    \param source_support local support
    \param source_group processor group containing the local processors
    \param target_group processor group containing the distant processors
    \param method interpolation method
  */
  InterpolationMatrix::InterpolationMatrix(
                                           const ParaMEDMEM::ParaMESH& source_support, 
                                           const ProcessorGroup& source_group,
                                           const ProcessorGroup& target_group,
                                           const DECOptions& dec_options,
                                           const INTERP_KERNEL::InterpolationOptions& interp_options):
    _source_support(*source_support.getMesh()),
    _mapping(source_group, target_group, dec_options),
    _source_group(source_group),
    _target_group(target_group),
    _source_volume(0),
    DECOptions(dec_options),
    INTERP_KERNEL::InterpolationOptions(interp_options)
  {
    int nbelems = _source_support.getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);
    _row_offsets.resize(nbelems+1);
    for (int i=0; i<nbelems+1; i++)
      _row_offsets[i]=0;

    _coeffs.resize(nbelems);
  }

  InterpolationMatrix::~InterpolationMatrix()
  {
  }

  /*!
    \brief Adds the contribution of a distant subdomain to the interpolation matrix.
    The method adds contribution to the interpolation matrix.
    For each row of the matrix, elements are addded as
    a (column, coeff) pair in the _coeffs array. This column number refers
    to an element on the target side via the _col_offsets array. It is made of a series 
    of (iproc, ielem) pairs. 
    The number of elements per row is stored in the row_offsets array.

    \param distant_support local representation of the distant subdomain
    \param iproc_distant id of the distant subdomain (in the distant group)
    \param distant_elems mapping between the local representation of the subdomain and the actual elem ids on the distant subdomain
  */
  void InterpolationMatrix::addContribution (MEDMEM::MESH& distant_support,
                                             int iproc_distant, int* distant_elems)
  {
    if (distant_support.getMeshDimension() != _source_support.getMeshDimension() ||
        distant_support.getMeshDimension() != _source_support.getMeshDimension() )
      throw MEDMEM::MEDEXCEPTION("local and distant meshes do not have the same space and mesh dimensions");

    //creating the interpolator structure
    vector<map<int,double> > surfaces;
    //computation of the intersection volumes between source and target elements
    int source_size = _source_support.getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);

    if (distant_support.getMeshDimension()==2 && distant_support.getSpaceDimension()==3)
    {
      MEDNormalizedUnstructuredMesh<3,2> target_wrapper (&distant_support);
      MEDNormalizedUnstructuredMesh<3,2> source_wrapper (&_source_support);
      INTERP_KERNEL::Interpolation3DSurf interpolator (*this);
      interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces);
    }
    else if (distant_support.getMeshDimension()==2 && distant_support.getSpaceDimension()==2)
    {
      MEDNormalizedUnstructuredMesh<2,2> target_wrapper(&distant_support);
      MEDNormalizedUnstructuredMesh<2,2> source_wrapper(&_source_support);

      INTERP_KERNEL::Interpolation2D interpolator (*this);
      interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces);
    }
    else if (distant_support.getMeshDimension()==3 && distant_support.getSpaceDimension()==3)
    {
      MEDNormalizedUnstructuredMesh<3,3> target_wrapper (&distant_support);
      MEDNormalizedUnstructuredMesh<3,3> source_wrapper (&_source_support);

      INTERP_KERNEL::Interpolation3D interpolator (*this);
      interpolator.interpolateMeshes(target_wrapper,source_wrapper,surfaces);
    }
    else
    {
      throw MEDMEM::MEDEXCEPTION("no interpolator exists for these mesh and space dimensions ");
    }

    if (surfaces.size() != source_size)
    {
      cout<<"surfaces.size()="<<surfaces.size()<<" source_size="<<source_size<<endl;  
      throw MEDEXCEPTION("uncoherent number of rows in interpolation matrix");
    }

    //computing the vectors containing the source and target element volumes
    MEDMEM::SUPPORT target_support (&distant_support, "all cells", MED_EN::MED_CELL);
    MEDMEM::FIELD<double>* target_triangle_surf = getSupportVolumes(target_support);
    MEDMEM::SUPPORT source_support (const_cast<MEDMEM::MESH*>(&_source_support),
                                    "all cells", MED_EN::MED_CELL);
    MEDMEM::FIELD<double>* source_triangle_surf = getSupportVolumes(source_support);

    //storing the source volumes
    _source_volume.resize(source_size);
    for (int i=0; i<source_size; i++)
      _source_volume[i] = source_triangle_surf->getValueIJ(i+1,1);

    //loop over the elements to build the interpolation
    //matrix structures
    for (int ielem=0; ielem < surfaces.size(); ielem++) 
    {
      _row_offsets[ielem+1] += surfaces[ielem].size();
      //_source_indices.push_back(make_pair(iproc_distant, ielem));
      for (map<int,double>::const_iterator iter = surfaces[ielem].begin();
           iter != surfaces[ielem].end();
           iter++)
      {
        //surface of the target triangle
        double surf = target_triangle_surf->getValueIJ(iter->first,1);

        //locating the (iproc, itriangle) pair in the list of columns
        vector<pair<int,int> >::iterator iter2 =
          find(_col_offsets.begin(), _col_offsets.end(),make_pair(iproc_distant,iter->first));
        int col_id;

        if (iter2 == _col_offsets.end())
        {
          //(iproc, itriangle) is not registered in the list
          //of distant elements

          _col_offsets.push_back(make_pair(iproc_distant,iter->first));
          col_id =_col_offsets.size();
          _mapping.addElementFromSource(iproc_distant,distant_elems[iter->first-1]);
          _target_volume.push_back(surf);
        }
        else 
        {
          col_id = iter2 - _col_offsets.begin() + 1;
        }

        //the non zero coefficient is stored 
        //ielem is the row,
        //col_id is the number of the column
        //iter->second is the value of the coefficient

        _coeffs[ielem].push_back(make_pair(col_id,iter->second));
      }
    }
    delete source_triangle_surf;
    delete target_triangle_surf;
  }
	
  /*! The call to this method updates the arrays on the target side
    so that they know which amount of data from which processor they 
    should expect. 
    That call makes actual interpolations via multiply method 
    available.
  */
  void InterpolationMatrix::prepare()
  {
    int nbelems = _source_support.getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);
    for (int ielem=0; ielem < nbelems; ielem++)
    {
      _row_offsets[ielem+1]+=_row_offsets[ielem];
    }  
    _mapping.prepareSendRecv();
  }

  /*!
    \brief performs t=Ws, where t is the target field, s is the source field

    The call to this method must be called both on the working side 
    and on the idle side. On the working side, the vector  T=VT^(-1).(W.S)
    is computed and sent. On the idle side, no computation is done, but the 
    result from the working side is received and the field is updated.

    \param field source field on processors involved on the source side, target field on processors on the target side
  */
  void InterpolationMatrix::multiply(MEDMEM::FIELD<double>& field) const
  {
    vector<double> target_value(_col_offsets.size()* field.getNumberOfComponents(),0.0);

    //computing the matrix multiply on source side
    if (_source_group.containsMyRank())
    {
      int nbcomp = field.getNumberOfComponents();
      int nbrows = _source_support.getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);

      // performing W.S
      // W is the intersection matrix
      // S is the source vector

      for (int irow=0; irow<nbrows; irow++)
        for (int icomp=0; icomp< nbcomp; icomp++)
        {
          double coeff_row = field.getValueIJ(irow+1,icomp+1);
          for (int icol=_row_offsets[irow]; icol< _row_offsets[irow+1];icol++)
          {
            int colid= _coeffs[irow][icol-_row_offsets[irow]].first;
            double value = _coeffs[irow][icol-_row_offsets[irow]].second;
            target_value[(colid-1)*nbcomp+icomp]+=value*coeff_row;
          }
        }

      // performing VT^(-1).(W.S)
      // where VT^(-1) is the inverse of the diagonal matrix containing 
      // the volumes of target cells

      for (int i=0; i<_col_offsets.size();i++)
        for (int icomp=0; icomp<nbcomp; icomp++)
          target_value[i*nbcomp+icomp] /= _target_volume[i];
    }

    if (_target_group.containsMyRank())
    {
      int nbelems = field.getSupport()->getNumberOfElements(MED_EN::MED_ALL_ELEMENTS);
      int nbcomp = field.getNumberOfComponents();
      double* value = const_cast<double*> (field.getValue());
      for (int i=0; i<nbelems*nbcomp; i++)
        value[i]=0.0;
    }
    //on source side : sending  T=VT^(-1).(W.S)
    //on target side :: receiving T and storing it in field
    _mapping.sendRecv(&target_value[0],field);
  }
  
  /*!
    \brief performs s=WTt, where t is the target field, s is the source field, WT is the transpose matrix from W

    The call to this method must be called both on the working side 
    and on the idle side. On the working side, the target vector T is received and the
    vector  S=VS^(-1).(WT.T) is computed to update the field. 
    On the idle side, no computation is done, but the field is sent.

    \param field source field on processors involved on the source side, target field on processors on the target side
  */
  void InterpolationMatrix::transposeMultiply(MEDMEM::FIELD<double>& field) const
  {
    int nbcomp = field.getNumberOfComponents();
    vector<double> source_value(_col_offsets.size()* nbcomp,0.0);
    _mapping.reverseSendRecv(&source_value[0],field);

    //treatment of the transpose matrix multiply on the source side
    if (_source_group.containsMyRank())
    {
      int nbrows = _source_support.getNumberOfElements(MED_EN::MED_CELL,MED_EN::MED_ALL_ELEMENTS);

      vector<double> target_value(nbrows*nbcomp,0.0);

      //performing WT.T
      //WT is W transpose
      //T is the target vector
      for (int irow = 0; irow < nbrows; irow++)
        for (int icol = _row_offsets[irow]; icol < _row_offsets[irow+1]; icol++)
        {
          int colid = _coeffs[irow][icol-_row_offsets[irow]].first;
          double value = _coeffs[irow][icol-_row_offsets[irow]].second;
          
          for (int icomp=0; icomp<nbcomp; icomp++)
          {
            double coeff_row = source_value[(colid-1)*nbcomp+icomp];
            target_value[irow*nbcomp+icomp] += value*coeff_row;
          }
        }

      //performing VS^(-1).(WT.T)
      //VS^(-1) is the inverse of the diagonal matrix storing
      //volumes of the source cells
      for (int i=0; i<nbrows; i++)
        for (int icomp=0; icomp<nbcomp; icomp++)
          target_value[i*nbcomp+icomp] /= _source_volume[i];

      //storing S= VS^(-1).(WT.T)
      for (int irow=0; irow<nbrows; irow++)
        for (int icomp=0; icomp<nbcomp; icomp++)
          field.setValueIJ(irow+1,icomp+1,target_value[irow*nbcomp+icomp]);
    }
  }

	/*!
		\brief returns the volumes of the cells underlying the field \a field

		For 2D geometries, the returned field contains the areas.
		For 3D geometries, the returned field contains the volumes.

		\param field field on which cells the volumes are required
		\return field containing the volumes
	*/
	MEDMEM::FIELD<double>* InterpolationMatrix::getSupportVolumes(const MEDMEM::SUPPORT& support)
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
