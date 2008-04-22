#ifndef INTERPOLATIONMATRIX_HXX_
#define INTERPOLATIONMATRIX_HXX_

#include "MEDMEM_Field.hxx"
#include "MPI_AccessDEC.hxx"
#include "MxN_Mapping.hxx"

namespace ParaMEDMEM
{
  class InterpolationMatrix
  {
  public:
    
    InterpolationMatrix(const ParaMEDMEM::ParaMESH& source_support, 
			const ProcessorGroup& local_group,
			const ProcessorGroup& distant_group, 
			const string& method);
    
    virtual ~InterpolationMatrix();
    void addContribution(MEDMEM::MESH& distant_support, int iproc_distant, int* distant_elems);
    void multiply(MEDMEM::FIELD<double>&) const;
		void transposeMultiply(MEDMEM::FIELD<double>&)const;
    void prepare();
    int getNbRows() const {return _row_offsets.size();}
    void setAllToAllMethod(const AllToAllMethod& method)
        { _mapping.setAllToAllMethod(method);}
		MPI_AccessDEC* getAccessDEC(){return _mapping.getAccessDEC();}

  private:

		MEDMEM::FIELD<double>* getSupportVolumes(const MEDMEM::SUPPORT& field);

	private:
    vector<int> _row_offsets;
    vector<pair<int,int> > _col_offsets;
    const MEDMEM::MESH& _source_support; 
    MxN_Mapping _mapping;
    string _method;
    const ProcessorGroup& _source_group;
    const ProcessorGroup& _target_group;
    vector<double> _target_volume;
		vector<double> _source_volume;
    vector<vector<pair<int,double> > > _coeffs;
  };
  
}

#endif /*INTERPOLATIONMATRIX_HXX_*/
