// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

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

    class KeyComparator
    {
    public:
      KeyComparator(int val):_val(val) { }
      bool operator()(const typename std::pair<int,T>& val) { return val.first==_val; }
    protected:
      int _val;
    };

    class Row : public std::vector< std::pair<int,T> >
    {
    public:
      Row():std::vector< std::pair<int,T> >(){};
      Row (const Row& row)
      {
        this->resize(row.size());
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
      typename std::vector< std::pair<int,T> >::const_iterator find(int elem) const
      {
        return std::find_if(std::vector< std::pair<int,T> >::begin(),std::vector< std::pair<int,T> >::end(),KeyComparator(elem));
      }

      void erase(int elem) { std::vector< std::pair<int,T> >::erase(std::find_if(std::vector< std::pair<int,T> >::begin(),std::vector< std::pair<int,T> >::end(),KeyComparator(elem))); }

      void insert(const std::pair<int,T>& myPair) { vector<std::pair<int,T> >::push_back(myPair); }
    };
    
  private:
    unsigned int _nb_rows;
    T* _coeffs;
    unsigned int* _cols;
    std::vector<unsigned int> _ncols_offset;
    std::vector< Row > _auxiliary_matrix;
    friend std::ostream& operator<<<>(std::ostream& in, const Matrix<T,type>& m);
    friend std::istream& operator>><>(std::istream& in, Matrix<T,type>& m);
    bool _is_configured;
  public:
    typedef Row value_type;
  public:
    Matrix():_coeffs(0), _cols(0), _nb_rows(0), _is_configured(false)
    { }
    Matrix(int nbrows):_coeffs(0), _cols(0), _is_configured(false)
    { _nb_rows=nbrows; }
    Matrix(std::vector<std::map<int,T> > & matrix) :
      _coeffs(0), _cols(0), _is_configured(false)
    {
      _nb_rows=matrix.size();
      _auxiliary_matrix.resize(_nb_rows);
      for (int i=0; i<_nb_rows; i++)
        {
          typename std::map<int,T>::iterator it;
          for (it=matrix[i].begin(); it != matrix[i].end(); it++)
            _auxiliary_matrix[i].push_back(*it);//MN: pq push_back plutot que simple affectation?
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
          _cols = new unsigned int[size];
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
          _cols = new unsigned int[size];
          memcpy(_coeffs, m._coeffs, size*sizeof(double));
          memcpy(_cols, m._cols, size*sizeof(int));
        }
      return this;
    }

    /*! declares a method that specifies the number of rows */
    void resize(unsigned int nbrows)
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
      
      for (unsigned int i=0; i< _auxiliary_matrix[OTT<int,type>::ind2C(irow)].size(); i++)
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
          output[i]=0.;
          for (unsigned int j=_ncols_offset[i]; j< _ncols_offset[i+1]; j++)
            {
              int icol = _cols[j];
              output[i]+=input[icol]*_coeffs[j];
            }
        }
    }

    /*!      
      Matrix multiplies vector \a input and stores the result in 
      vector \a output.
      input and output are supposed to represent the same field 
      discretised on two different on meshes.
      nb_comp is the number of components of the fields input and output
      The vector pointed by \a input must be dimensioned
      to the number of columns times nb_comp while the vector pointed by output must be
      dimensioned to the number of rows times nb_comp.
    */
    void multiply(const T* const input, T* const output, int nb_comp)
    {
      if (!_is_configured)
        configure();
      
      for (int i=0; i< _nb_rows; i++)
        {
          for(int comp = 0; comp < nb_comp; comp++)
            output[i*nb_comp+comp]=0.;
          for (unsigned int j=_ncols_offset[i]; j< _ncols_offset[i+1]; j++)
            {
              int icol = _cols[j];
              for(int comp = 0; comp < nb_comp; comp++)
                output[i*nb_comp+comp]+=input[icol*nb_comp+comp]*_coeffs[j];
            }
        }
    }   
    /*!      
      Transpose-multiplies vector \a input and stores the result in 
      vector \a output.
      nb_cols is the number of columns of the matrix, (it is not an attribute of the class) 
      The vector pointed by \a input must be dimensioned
      to the number of lines _nb_rows while the vector pointed by output must be
      dimensioned to the number of columns nb_cols.
    */
    void transposeMultiply(const T* const input, T* const output, int nb_cols)
    {
      if (!_is_configured)
        configure();
      
      for (int icol=0; icol< nb_cols; icol++)
        output[icol]=0.;
      for (int i=0; i< _nb_rows; i++)
        {
          for (unsigned int j=_ncols_offset[i]; j< _ncols_offset[i+1]; j++)
            {
              int icol = _cols[j];
              output[icol]+=input[i]*_coeffs[j];
            }
        }
    }
    /*!      
      Transpose-multiplies vector \a input and stores the result in 
      vector \a output.
      input and output are supposed to represent the same field 
      discretised on two different on meshes.
      nb_comp is the number of components of the fields input and output
      nb_cols is the number of columns of the matrix, (it is not an attribute of the class) 
      The vector pointed by \a input must be dimensioned
      to _nb_rows*nb_comp while the vector pointed by output must be
      dimensioned to nb_cols*nb_comp.
    */
    void transposeMultiply(const T* const input, T* const output, int nb_cols, int nb_comp)
    {
      if (!_is_configured)
        configure();
      
      for (int icol=0; icol< nb_cols; icol++)
        for(int comp = 0; comp < nb_comp; comp++)
          output[icol*nb_comp+comp]=0.;

      for (int i=0; i< _nb_rows; i++)
        {
          for (unsigned int j=_ncols_offset[i]; j< _ncols_offset[i+1]; j++)
            {
              int icol = _cols[j];
              for(int comp = 0; comp < nb_comp; comp++)
                output[icol*nb_comp+comp]+=input[i*nb_comp+comp]*_coeffs[j];
            }
        }
    }
    /*
      Sums the coefficients of each column of the matrix
      nb_cols is the number of columns of the matrix, (it is not an attribute of the class) 
      The vector output must be dimensioned to nb_cols
    */
    void colSum(std::vector< T >& output, int nb_cols)
    {
      if (!_is_configured)
        configure();
      for (int icol=0; icol< nb_cols; icol++)
        output[icol]=0.;
      for (int i=0; i< _nb_rows; i++)
        {
          for (unsigned int j=_ncols_offset[i]; j< _ncols_offset[i+1]; j++)
            {
              int icol = _cols[j];
              output[icol]+=_coeffs[j];
            }
        }
    }

    /*
      Sums the coefficients of each row of the matrix
      The vector output must be dimensioned to _nb_rows
    */
    void rowSum(std::vector< T >& output)
    {
      if (!_is_configured)
        configure();
      for (int i=0; i< _nb_rows; i++)
        {
          output[i]=0;
          for (unsigned int j=_ncols_offset[i]; j< _ncols_offset[i+1]; j++) 
            output[i]+=_coeffs[j];
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
      for (unsigned int i=0; i<_nb_rows; i++)
        _ncols_offset[i+1]=_ncols_offset[i]+_auxiliary_matrix[i].size();
      int nbcoeffs= _ncols_offset[_nb_rows];
      _cols=new unsigned int[nbcoeffs];
      _coeffs=new T[nbcoeffs];
      unsigned int* cols_ptr=_cols;
      T* coeffs_ptr=_coeffs;
      for (unsigned int i=0; i<_nb_rows; i++)
        {
          for (unsigned int j=0; j<_auxiliary_matrix[i].size(); j++)
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
    Row &operator [] (unsigned int irow)
    {
      return _auxiliary_matrix[irow];
    }

    int getNbRows()
    {
      return _nb_rows;
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
        out << OTT<unsigned int,type>::indFC(0) <<std::endl;
        out << m._nb_rows<<std::endl;
        for (unsigned int i=0; i<m._nb_rows; i++)
          {
            out << m._ncols_offset[i+1]-m._ncols_offset[i];
            for (unsigned int j=m._ncols_offset[i]; j<m._ncols_offset[i+1]; j++)
              out <<"\t"<< OTT<unsigned int,type>::indFC(m._cols[j]) <<"\t"<<m._coeffs[j];
            out<<std::endl;
          }
      }
    else
      {
        out << OTT<unsigned int,type>::indFC(0) <<"\n";
        out << m._nb_rows <<"\n";
        for (unsigned int i=0; i<m._nb_rows; i++)
          {
            out<< m._auxiliary_matrix[i].size();
            for (unsigned int j=0; j<m._auxiliary_matrix[i].size(); j++)
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
    if (index_base_test!=OTT<unsigned int,type>::indFC(0))
      {
        std::cerr << "file index is "<<index_base_test<<std::endl;
        throw Exception("incompatible indexing reading matrix");
      }
    in >> m._nb_rows;
    m._auxiliary_matrix.resize(m._nb_rows);
    for (unsigned int i=0; i<m._nb_rows; i++)
      {
        unsigned int ncols;
        in >> ncols;
        m._auxiliary_matrix[i].resize(ncols);
        double value;
        unsigned int col;
        for (unsigned int j=0; j<ncols; j++)
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
