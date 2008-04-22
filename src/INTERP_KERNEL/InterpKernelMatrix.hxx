#ifndef __INTERPKERNELMATRIX_HXX_
#define __INTERPKERNELMATRIX_HXX__

#include "InterpolationUtils.hxx"

#include <vector>
#include <iostream>
#include <ostream>
#include <istream>
#include <map>

namespace INTERP_KERNEL
{
  template<class T, NumberingPolicy type>
  class Matrix;
  
  template<class U, NumberingPolicy type>
  std::ostream& operator<<(std::ostream& in, const Matrix<U,type>& m);
  template<class U, NumberingPolicy type>
  std::istream& operator>>(std::istream& in, Matrix<U,type>& m);
	
  template<class T, NumberingPolicy type=ALL_C_MODE>
  class Matrix
  {
    class Row : public std::vector< typename std::pair<int,T> >
    {
		public:
			Row():std::vector< typename std::pair<int,T> >(){};
			Row (const Row& row)
			{
				resize(row.size());
				for (int i=0; i<this->size(); i++)
					(*this)[i]=row[i];
			}
			Row& operator= (const Row& row)
			{
				this->resize(row.size());
				for (int i=0; i<this->size(); i++)
					(*this)[i]=row[i];
				return *this;
			}
			
    
      void insert(const std::pair<int,T>& myPair) { push_back(myPair); }
    };
    
  private:
    uint _nb_rows;
    T* _coeffs;
    uint* _cols;
    std::vector<uint> _ncols_offset;
    std::vector< Row > _auxiliary_matrix;
    friend std::ostream& operator<<<>(std::ostream& in, const Matrix<T,type>& m);
    friend std::istream& operator>><>(std::istream& in, Matrix<T,type>& m);
    bool _is_configured;
  public:
    Matrix():_coeffs(0), _cols(0), _nb_rows(0), _is_configured(false)
    { }
    Matrix(int nbrows):_coeffs(0), _cols(0), _is_configured(false)
    { _nb_rows=nbrows; }
    Matrix(std::vector<std::map<int,T> > & matrix) :
      _coeffs(0), _cols(0), _is_configured(false)
    {
      _nb_rows=matrix.size();
      for (int i=0; i<_nb_rows; i++)
        {
          _auxiliary_matrix[i].resize(matrix[i].size());
          typename std::map<int,T>::iterator it;
          for (it=matrix[i].begin(); it != matrix[i].end(); it++)
            _auxiliary_matrix[i].push_back(*it);
        }
      
    }
    /*!Copy constructor
     */
    Matrix(const Matrix & m)
    {
      _is_configured=m._is_configured;
      _nb_rows=m._nb_rows;
      _auxiliary_matrix=m._auxiliary_matrix;
      _ncols_offset=m._ncols_offset;
      if (_is_configured)
        {
          int size=_ncols_offset[_nb_rows];
          _coeffs = new double[size];
          _cols = new uint[size];
          memcpy(_coeffs, m._coeffs, size*sizeof(double));
          memcpy(_cols, m._cols, size*sizeof(int));
        }
    }
    
    ~Matrix()
    {
      delete[] _coeffs;
      delete[] _cols;
    }

		Matrix& operator=(const Matrix& m)
		{
			 _is_configured=m._is_configured;
      _nb_rows=m._nb_rows;
      _auxiliary_matrix=m._auxiliary_matrix;
      _ncols_offset=m._ncols_offset;
      if (_is_configured)
        {
          int size=_ncols_offset[_nb_rows];
          _coeffs = new double[size];
          _cols = new uint[size];
          memcpy(_coeffs, m._coeffs, size*sizeof(double));
          memcpy(_cols, m._cols, size*sizeof(int));
        }
			return this;
		}

    /*! declares a method that specifies the number of rows */
    void resize(uint nbrows)
    {
      _nb_rows=nbrows;
      _auxiliary_matrix.resize(nbrows);
    }
    
    /*! sets (i,j) coefficient to \a value */
    void setIJ(int irow, int icol, T value)
    {
      if (_is_configured)
        throw Exception("filling a configured matrix");
      if (_auxiliary_matrix.empty())
        _auxiliary_matrix.resize(_nb_rows);
      
      for (uint i=0; i< _auxiliary_matrix[OTT<int,type>::ind2C(irow)].size(); i++)
        if (_auxiliary_matrix[OTT<int,type>::ind2C(irow)][i].first == icol)
          {
            _auxiliary_matrix[OTT<int,type>::ind2C(irow)][i].second = value;
            return;
          }
      _auxiliary_matrix[OTT<int,type>::ind2C(irow)].push_back(std::make_pair(icol, value));
    }
    
    /*!
      
      Matrix multiplies vector \a input and stores the result in 
      vector \a output.
      The vector pointed by \a input must be dimensioned
      to the number of columns while the vector pointed by output must be
      dimensioned to the number of rows.
    */
    void multiply(const T* const input, T* const output)
    {
      if (!_is_configured)
        configure();
      
      for (int i=0; i< _nb_rows; i++)
        {
          output[i]=0;
          for (int j=_ncols_offset[i]; j< _ncols_offset[i+1]; j++) {
            int icol = _cols[j];
            output[i]+=input[icol]*_coeffs[j];
          }
        }
    }
    
    /*! This operation freezes the profile of the matrix
      and puts it under a CSR form so that it becomes
      efficient both in terms of memory occupation and
      in terms of multiplication */
    
    void configure()
    {
      _ncols_offset.resize(_nb_rows+1);
      _ncols_offset[0]=0;
      for (int i=0; i<_nb_rows; i++)
        _ncols_offset[i+1]=_ncols_offset[i]+_auxiliary_matrix[i].size();
      int nbcoeffs= _ncols_offset[_nb_rows];
      _cols=new uint[nbcoeffs];
      _coeffs=new T[nbcoeffs];
      uint* cols_ptr=_cols;
      T* coeffs_ptr=_coeffs;
      for (uint i=0; i<_nb_rows; i++)
        {
          for (uint j=0; j<_auxiliary_matrix[i].size(); j++)
            {
              *cols_ptr++ = OTT<int,type>::ind2C(_auxiliary_matrix[i][j].first);
              *coeffs_ptr++ = _auxiliary_matrix[i][j].second;
            }
        }
      _auxiliary_matrix.clear();
      _is_configured=true;
    }

    /*! 
     * 0 <= irow < n
     */
    Row &operator [] (uint irow)
    {
      return _auxiliary_matrix[irow];
    }
    
  };
  
  /*! output to an ascii file
    only nonzero elements are written
    - the first line contains the indexing (0 or 1)
    - the second line contains the number of rows.
    - for each row, a line contains: 
    - the number of nonzero coeffs 
    - and for each coeff : icol, value
    
    for instance, matrix 
    | 1.0 0.0 0.5 |
    | 0.0 1.0 0.0 |
    | 0.2 0.0 1.0 |
    will be displayed in 0-indexing as
    0
    3
    2 0 1.0 2 0.5
    1 1 1.0
    2 0 0.2 2 1.0
  */
  
  template<class T, NumberingPolicy type>
  std::ostream& operator<<(std::ostream& out, const Matrix<T, type>& m)
  {
    if (m._is_configured)
      {
        out << OTT<uint,type>::indFC(0) <<std::endl;
        out << m._nb_rows<<std::endl;
        for (uint i=0; i<m._nb_rows; i++)
          {
            out << m._ncols_offset[i+1]-m._ncols_offset[i];
            for (uint j=m._ncols_offset[i]; j<m._ncols_offset[i+1]; j++)
              out <<"\t"<< OTT<uint,type>::indFC(m._cols[j]) <<"\t"<<m._coeffs[j];
            out<<std::endl;
          }
      }
    else
      {
        out << OTT<uint,type>::indFC(0) <<"\n";
        out << m._nb_rows <<"\n";
        for (uint i=0; i<m._nb_rows; i++)
          {
            out<< m._auxiliary_matrix[i].size();
            for (uint j=0; j<m._auxiliary_matrix[i].size(); j++)
              out << "\t" <<m._auxiliary_matrix[i][j].first <<"\t"
                <<m._auxiliary_matrix[i][j].second;
            out <<"\n";
          }
      }
    return out;
  }
  
  template<class T, NumberingPolicy type>
  std::istream& operator>>(std::istream& in, Matrix<T,type>& m)
  {
    int index_base_test;
    in >> index_base_test;
    if (index_base_test!=OTT<uint,type>::indFC(0))
      {
        std::cerr << "file index is "<<index_base_test<<std::endl;
        throw Exception("incompatible indexing reading matrix");
      }
    in >> m._nb_rows;
    m._auxiliary_matrix.resize(m._nb_rows);
    for (uint i=0; i<m._nb_rows; i++)
      {
        uint ncols;
        in >> ncols;
        m._auxiliary_matrix[i].resize(ncols);
        double value;
        uint col;
        for (uint j=0; j<ncols; j++)
          {
            in>>col;
            in>>value;
            m._auxiliary_matrix[i].push_back(std::make_pair(col, value));
          }
      }
    return in;
  }
}

#endif
