#ifndef REMAPPER_HXX_
#define REMAPPER_HXX_

#include "InterpKernelMatrix.hxx"
#include "MEDMEM_Mesh.hxx"
#include "MEDMEM_Support.hxx"
#include "MEDMEM_Field.hxx"

namespace INTERP_KERNEL
{

class INTERPKERNEL_EXPORT Remapper
{
public:
	Remapper();
	virtual ~Remapper();
	void prepare(const MEDMEM::MESH& mesh_source, const MEDMEM::MESH& mesh_target);
	void transfer(const MEDMEM::FIELD<double>& field_source, MEDMEM::FIELD<double>& field_target);
	void setOptionDouble(const std::string& key, double value);
	void setOptionInt(const std::string& key, int value);
private :
	Matrix<double,ALL_FORTRAN_MODE>* _matrix;
	MEDMEM::FIELD<double>* getSupportVolumes(const MEDMEM::SUPPORT& support);

};

}

#endif /*REMAPPER_HXX_*/
