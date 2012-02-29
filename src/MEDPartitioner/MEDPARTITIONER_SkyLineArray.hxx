//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
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

# ifndef __MEDPARTITIONER_MEDSKYLINEARRAY_H__
# define __MEDPARTITIONER_MEDSKYLINEARRAY_H__

#include "MEDPARTITIONER.hxx"

#include <vector>

namespace MEDPARTITIONER {
class MEDPARTITIONER_EXPORT MEDSKYLINEARRAY
{
private :
  std::vector<int> _index;
  std::vector<int> _value;

public :
  // Attention, avec ce constructeur, il n'est possible de remplir le MEDSKYLINEARRAY 
  MEDSKYLINEARRAY();

  // Constructeur par recopie
  MEDSKYLINEARRAY( const MEDSKYLINEARRAY &myArray );

  // Avec ce constructeur la m�moire pour le tableau  de valeur et le
  // tableau d'index est r�serv�e. Il suffit d'effectuer les s�quences
  // d'appels suivantes pour initialiser le MEDSKYLINEARRAY
  // 1) setIndex(index) puis <count> fois setI(i,&listValeurN�I) avec i dans 1..count
  //    rem :   listValeurN�I est dupliqu�e
  // 2) appeler <length> fois setIJ(i,j,valeur) avec i dans 1..count et avec j dans 1..count
  MEDSKYLINEARRAY( const std::vector<int>& index, const std::vector<int>& value );

  ~MEDSKYLINEARRAY();
  //void setMEDSKYLINEARRAY( const int count, const int length, int* index , int* value ) ;

  inline int getNumberOf() const;
  inline int getLength() const;
  inline const int* getIndex() const;
  inline const int* getValue() const;
};

// ---------------------------------------
//              Methodes Inline
// ---------------------------------------
inline int MEDSKYLINEARRAY::getNumberOf() const
{
  return _index.size()-1;
}
inline int MEDSKYLINEARRAY::getLength() const
{
  return _value.size() ;
}
inline const int*  MEDSKYLINEARRAY::getIndex() const
{
  return (const int*)(&_index[0]) ;
} 
inline const int*  MEDSKYLINEARRAY::getValue() const
{
  return (const int*)(&_value[0]) ;
} 
}
# endif
