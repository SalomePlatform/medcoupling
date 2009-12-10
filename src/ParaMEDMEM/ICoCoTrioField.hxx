//////////////////////////////////////////////////////////////////////////////
//
// File:        ICoCoTrioField.h
// Directory:   $TRIO_U_ROOT/Kernel/Framework
// Version:     
//
//////////////////////////////////////////////////////////////////////////////

#ifndef _ICOCOTRIOFIELD_HXX_
#define _ICOCOTRIOFIELD_HXX_

#include <ICoCoField.hxx>
namespace ICoCo
{
        /*!
                \brief structure for coupling Trio codes via the ICoCo interface

                This structure contains all the necessary information 
                for constructing a ParaMEDMEM::ParaFIELD (with the addition of the MPI
                communicator). The ICoCo API specifies two kinds of calls for
                the ICoCo::Field : either with the mesh only or with the entire information (mesh and field).
                This structure can therefore be left without _time, _nb_field_components, _field
                information, which are related to the field values.

                _coords and _connectivity tables are always owned by the TrioField.

         */
  class TrioField:public Field
  {
  public:
    
    TrioField();
    ~TrioField();
    void clear();
    void print();
    void set_standalone();
    void dummy_geom();
    TrioField& operator=(const TrioField& NewField);
    void save(std::ostream& os) const;
    void restore(std::istream& in);
    int nb_values() const ;
  public:
    int _type ; // 0 elem 1 nodes
    int _mesh_dim;
    int _space_dim;
    int _nbnodes;
    int _nodes_per_elem;
    int _nb_elems;
    int _itnumber;
    int* _connectivity;
    double* _coords;
    
    double _time1,_time2;
    int _nb_field_components;
    double* _field;
    bool _has_field_ownership;
  }; 
}

#endif
