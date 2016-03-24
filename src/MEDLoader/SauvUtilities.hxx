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
// File      : SauvUtilities.hxx
// Created   : Mon Aug 22 18:27:34 2011
// Author    : Edward AGAPOV (eap)

#ifndef __SAUVUTILITIES_HXX__
#define __SAUVUTILITIES_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "NormalizedUnstructuredMesh.hxx"

#include <string>
#include <sstream>

#define THROW_IK_EXCEPTION(text)                        \
  {                                                     \
    std::ostringstream oss; oss << text;                \
    throw INTERP_KERNEL::Exception(oss.str().c_str());  \
  }

namespace SauvUtilities
{
  INTERP_KERNEL::NormalizedCellType MEDLOADER_EXPORT gibi2medGeom( size_t gibiType );
  int med2gibiGeom( INTERP_KERNEL::NormalizedCellType medGeomType );
  const int * getGibi2MedQuadraticInterlace( INTERP_KERNEL::NormalizedCellType type );
  unsigned getDimension( INTERP_KERNEL::NormalizedCellType type );

  enum Readable_Piles
    {
      PILE_SOUS_MAILLAGE=1 ,
      PILE_NODES_FIELD  =2 ,
      PILE_TABLES       =10,
      PILE_LREEL        =18,
      PILE_LOGIQUES     =24,
      PILE_FLOATS       =25,
      PILE_INTEGERS     =26,
      PILE_STRINGS      =27,
      PILE_LMOTS        =29,
      PILE_NOEUDS       =32,
      PILE_COORDONNEES  =33,
      PILE_MODL         =38,
      PILE_FIELD        =39,
      PILE_LAST_READABLE=39
    };

  //================================================================================
  /*!
   * \brief Converts anything to string
   */
  //================================================================================

  template<class T> std::string toString(const T& anything)
  {
    std::ostringstream s; s << anything; return s.str();
  }

  // ==============================================================================
  // IMP 0020434: mapping GIBI names to MED names
  struct nameGIBItoMED
  {
    // GIBI value
    int gibi_pile;    // PILE_SOUS_MAILLAGE or PILE_FIELD/PILE_NODES_FIELD, or PILE_STRINGS(for components)
    int gibi_id;
    std::string gibi_name; // used only for components
    // MED value
    // med_pile = 27; // PILE_STRINGS
    int         med_id;    // used only on reading
    std::string med_name;  // used only on writing
  };

  // ==============================================================================
  /*!
   * \brief Base class for ASCII and XDR file readers
   */
  class FileReader : public MEDCoupling::RefCountObject
  {
  public:
    FileReader(const char* fileName);
    virtual ~FileReader() {}
    virtual bool isASCII() const = 0;

    virtual bool open() = 0;
    virtual bool getNextLine (char* & line, bool raiseOEF = true ) = 0;
    virtual void initNameReading(int nbValues, int width = 8) = 0;
    virtual void initIntReading(int nbValues) = 0;
    virtual void initDoubleReading(int nbValues) = 0;
    virtual bool more() const = 0;
    virtual void next() = 0;
    virtual int  index() const { return _iRead; }
    virtual int    getInt() const = 0;
    virtual float  getFloat() const = 0;
    virtual double getDouble() const = 0;
    virtual std::string getName() const = 0;
  protected:
    std::size_t getHeapMemorySizeWithoutChildren() const { return 0; }
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const { return std::vector<const BigMemoryObject *>(); }
  protected:
    std::string _fileName, _curLocale;
    int _iRead, _nbToRead;
  };
}
#endif
