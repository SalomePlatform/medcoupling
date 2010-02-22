#ifndef __INTERPKERNELUNIT_HXX__
#define __INTERPKERNELUNIT_HXX__

#include "INTERPKERNELEXPREVALDefines.hxx"
#include "InterpKernelException.hxx"

#include <map>
#include <sstream>

namespace INTERP_KERNEL
{
  class INTERPKERNELEXPREVAL_EXPORT UnitDataBase
  {
  public:
    UnitDataBase();
    const short *getInfoForUnit(const std::string& unit,
                                double& addFact, double& mFact) const throw(INTERP_KERNEL::Exception);
    static UnitDataBase _uniqueMapForExpr;
    static const int SIZE_OF_UNIT_BASE=5;
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

  class INTERPKERNELEXPREVAL_EXPORT DecompositionInUnitBase
  {
  public:
    DecompositionInUnitBase();
    void setInfo(const short *vals, double addFact, double mFact);
    short operator[](int i) const { return _value[i]; }
    bool operator==(const DecompositionInUnitBase& other) const;
    void getTranslationParams(const DecompositionInUnitBase& other, double& mul, double& add) const;
    bool isEqual(short mass, short lgth, short time, short intensity, short temp,
                 double add, double mult);
    bool isUnitary() const;
    //! \b WARNING no test is done on the fact that unit is adimensionnal.
    void negate();
    bool isAdimensional() const;
    void tryToConvertInUnit(double val) throw(INTERP_KERNEL::Exception);
    DecompositionInUnitBase &operator*(const DecompositionInUnitBase& other);
    DecompositionInUnitBase &operator/(const DecompositionInUnitBase& other);
    DecompositionInUnitBase &operator^(const DecompositionInUnitBase& other) throw(INTERP_KERNEL::Exception);
  private:
    void dealWithAddFactor(const DecompositionInUnitBase& other);
    static int couldItBeConsideredAsInt(double val) throw(INTERP_KERNEL::Exception);
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
  class INTERPKERNELEXPREVAL_EXPORT Unit
  {
  public:
    Unit(const char *reprC, bool tryToInterp=true);
    Unit(const char *reprFortran, int sizeOfRepr, bool tryToInterp=true);
    void tryToInterprate() const;
    bool isInterpretationOK() const;
    bool isCompatibleWith(const Unit& other) const;
    double convert(const Unit& target, double sourceVal) const;
    std::string getCoarseRepr() const;
  private:
    std::string _coarse_repr;
    mutable bool _is_interpreted;
    mutable bool _is_interpretation_ok;
    mutable DecompositionInUnitBase _decomp_in_base;
  };
}

#endif
