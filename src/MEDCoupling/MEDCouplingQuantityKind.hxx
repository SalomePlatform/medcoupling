// Copyright (C) 2026  CEA, EDF
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

#pragma once

#include "MEDCouplingRefCountObject.hxx"
#include "MCAuto.hxx"

#include <string>
#include <vector>

namespace MEDCoupling
{
class QuantityKindAbstract : public RefCountObject
{
   protected:
    virtual ~QuantityKindAbstract();

   public:
    MEDCOUPLING_EXPORT virtual std::string repr() const = 0;

    MEDCOUPLING_EXPORT virtual MCAuto<QuantityKindAbstract> clone() const = 0;

    MEDCOUPLING_EXPORT virtual std::string getClassName() const = 0;

    MEDCOUPLING_EXPORT virtual bool isEqual(const QuantityKindAbstract *other, std::string &reason) const = 0;

    MEDCOUPLING_EXPORT virtual std::string serialize() const = 0;

    MEDCOUPLING_EXPORT static MCAuto<QuantityKindAbstract> Deserialize(const std::string &s);

    std::vector<const MEDCoupling::BigMemoryObject *> getDirectChildrenWithNull() const override;

   protected:
    static std::string hexEncode(const std::string &input);
    static std::string hexDecode(const std::string &input);

   private:
    static int hexValue(char c);
};

class QuantityKindUnDef : public QuantityKindAbstract
{
   private:
    QuantityKindUnDef() = default;
    ~QuantityKindUnDef() = default;

   public:
    MEDCOUPLING_EXPORT std::string repr() const override;
    MEDCOUPLING_EXPORT virtual MCAuto<QuantityKindAbstract> clone() const override;
    MEDCOUPLING_EXPORT std::string getClassName() const override { return "QuantityKindUnDef"; }
    MEDCOUPLING_EXPORT bool isEqual(const QuantityKindAbstract *other, std::string &reason) const override;
    MEDCOUPLING_EXPORT static MCAuto<QuantityKindUnDef> New();
    std::size_t getHeapMemorySizeWithoutChildren() const override;
    MEDCOUPLING_EXPORT std::string serialize() const override;
};

class QuantityKindEnum : public QuantityKindAbstract
{
   private:
    explicit QuantityKindEnum(std::string value);

   public:
    MEDCOUPLING_EXPORT std::string repr() const override;

    MEDCOUPLING_EXPORT virtual MCAuto<QuantityKindAbstract> clone() const override;

    MEDCOUPLING_EXPORT std::string getClassName() const override { return "QuantityKindEnum"; }

    MEDCOUPLING_EXPORT bool isEqual(const QuantityKindAbstract *other, std::string &reason) const override;

    MEDCOUPLING_EXPORT static MCAuto<QuantityKindEnum> New(std::string value);

    MEDCOUPLING_EXPORT const std::string &value() const;

    MEDCOUPLING_EXPORT static const std::vector<std::string> &AllowedValues();

    MEDCOUPLING_EXPORT static int MCIdOfValue(const std::string &value);

    MEDCOUPLING_EXPORT static bool IsValidValue(const std::string &value);

    MEDCOUPLING_EXPORT static std::string AllowedValuesAt(int id);

    MEDCOUPLING_EXPORT static std::string AllowedValuesAsString();

    MEDCOUPLING_EXPORT std::string serialize() const override;

    std::size_t getHeapMemorySizeWithoutChildren() const override;

   private:
    std::string _value;
};

class QuantityKindUser : public QuantityKindAbstract
{
   private:
    explicit QuantityKindUser(std::string value);
    ~QuantityKindUser() = default;

   public:
    MEDCOUPLING_EXPORT std::string repr() const override;

    MEDCOUPLING_EXPORT virtual MCAuto<QuantityKindAbstract> clone() const override;

    MEDCOUPLING_EXPORT std::string getClassName() const override { return "QuantityKindUser"; }

    MEDCOUPLING_EXPORT bool isEqual(const QuantityKindAbstract *other, std::string &reason) const override;

    MEDCOUPLING_EXPORT static MCAuto<QuantityKindUser> New(std::string value);

    MEDCOUPLING_EXPORT const std::string &value() const;

    MEDCOUPLING_EXPORT std::string serialize() const override;

    std::size_t getHeapMemorySizeWithoutChildren() const override;

   private:
    std::string _value;
};
}  // namespace MEDCoupling
