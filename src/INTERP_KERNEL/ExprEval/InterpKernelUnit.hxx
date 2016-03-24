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
// Author : Anthony Geay (CEA/DEN)

#ifndef __INTERPKERNELUNIT_HXX__
#define __INTERPKERNELUNIT_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpKernelException.hxx"

#include <map>
#include <sstream>

namespace INTERP_KERNEL
{
  class UnitDataBase
  {
  public:
    INTERPKERNEL_EXPORT UnitDataBase();
    INTERPKERNEL_EXPORT const short *getInfoForUnit(const std::string& unit, double& addFact, double& mFact) const;
    INTERPKERNEL_EXPORT static UnitDataBase _uniqueMapForExpr;
    INTERPKERNEL_EXPORT static const int SIZE_OF_UNIT_BASE=5;
  private:
    std::map<std::string,double> _prefix_pow_10;
    std::map<std::string,const short *> _units_semantic;
    std::map<std::string,double> _units_mul;
    std::map<std::string,double> _units_add;
  private:
    static const int NB_OF_PREF_POW10=22;
    static const char *PREF_POW10[NB_OF_PREF_POW10];
    static const double POW10[NB_OF_PREF_POW10];
    static const int NB_OF_UNITS_RECOGN=29;
    static const char *UNITS_RECOGN[NB_OF_UNITS_RECOGN];
    static const short PROJ_IN_BASE[NB_OF_UNITS_RECOGN][SIZE_OF_UNIT_BASE];
    static const double MUL_COEFF[NB_OF_UNITS_RECOGN];
    static const double ADD_COEFF[NB_OF_UNITS_RECOGN];
  };

  class DecompositionInUnitBase
  {
  public:
    INTERPKERNEL_EXPORT DecompositionInUnitBase();
    INTERPKERNEL_EXPORT void setInfo(const short *vals, double addFact, double mFact);
    INTERPKERNEL_EXPORT short operator[](int i) const { return _value[i]; }
    INTERPKERNEL_EXPORT bool operator==(const DecompositionInUnitBase& other) const;
    INTERPKERNEL_EXPORT void getTranslationParams(const DecompositionInUnitBase& other, double& mul, double& add) const;
    INTERPKERNEL_EXPORT bool isEqual(short mass, short lgth, short time, short intensity, short temp,
                                     double add, double mult);
    INTERPKERNEL_EXPORT bool isUnitary() const;
    //! \b WARNING no test is done on the fact that unit is adimensionnal.
    INTERPKERNEL_EXPORT void negate();
    INTERPKERNEL_EXPORT bool isAdimensional() const;
    INTERPKERNEL_EXPORT void tryToConvertInUnit(double val);
    INTERPKERNEL_EXPORT DecompositionInUnitBase &operator*(const DecompositionInUnitBase& other);
    INTERPKERNEL_EXPORT DecompositionInUnitBase &operator/(const DecompositionInUnitBase& other);
    INTERPKERNEL_EXPORT DecompositionInUnitBase &operator^(const DecompositionInUnitBase& other);
  private:
    void dealWithAddFactor(const DecompositionInUnitBase& other);
    static int couldItBeConsideredAsInt(double val);
    static bool areDoubleEquals(double a, double b);
    static double powInt(double val, int exp);
  private:
    short _value[UnitDataBase::SIZE_OF_UNIT_BASE];
    double _add_to_base;
    double _mult_fact_to_base;
  };

  /*!
   * This class deals with units.
   * This class has two main responsabilities :
   *      - interprete units by giving simply their representation in string type.
   *      - performing operations on these units.
   *
   * All the possible units are represented with a unique tuple with 5 elements
   * representing the unique decomposition of a unit in the following base.
   *
   * dimension 0 stands for mass in g (\b NOT kg to simplify parsing).
   * dimension 1 stands for length in m.
   * dimension 2 stands for time in s.
   * dimension 3 stands for elec intensity A.
   * dimension 4 stands for temperature in K.
   */
  class Unit
  {
  public:
    INTERPKERNEL_EXPORT Unit(const char *reprC, bool tryToInterp=true);
    INTERPKERNEL_EXPORT Unit(const char *reprFortran, int sizeOfRepr, bool tryToInterp=true);
    INTERPKERNEL_EXPORT void tryToInterprate() const;
    INTERPKERNEL_EXPORT bool isInterpretationOK() const;
    INTERPKERNEL_EXPORT bool isCompatibleWith(const Unit& other) const;
    INTERPKERNEL_EXPORT double convert(const Unit& target, double sourceVal) const;
    INTERPKERNEL_EXPORT std::string getCoarseRepr() const;
  private:
    std::string _coarse_repr;
    mutable bool _is_interpreted;
    mutable bool _is_interpretation_ok;
    mutable DecompositionInUnitBase _decomp_in_base;
  };
}

#endif
