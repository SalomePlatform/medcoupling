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

#include "MEDCouplingQuantityKind.hxx"
#include "InterpKernelException.hxx"

using namespace MEDCoupling;

QuantityKindAbstract::~QuantityKindAbstract() = default;

std::vector<const MEDCoupling::BigMemoryObject *>
QuantityKindAbstract::getDirectChildrenWithNull() const
{
    return {};
}

std::string
QuantityKindAbstract::hexEncode(const std::string &input)
{
    constexpr char digits[] = "0123456789ABCDEF";

    std::string out;
    out.reserve(input.size() * 2);

    for (unsigned char c : input)
    {
        out.push_back(digits[c >> 4]);
        out.push_back(digits[c & 0x0F]);
    }

    return out;
}

std::string
QuantityKindAbstract::hexDecode(const std::string &input)
{
    if (input.size() % 2 != 0)
    {
        THROW_IK_EXCEPTION("Invalid hex string length");
    }

    std::string out;
    out.reserve(input.size() / 2);

    for (std::size_t i = 0; i < input.size(); i += 2)
    {
        int hi = hexValue(input[i]);
        int lo = hexValue(input[i + 1]);

        if (hi < 0 || lo < 0)
        {
            THROW_IK_EXCEPTION("Invalid hex character");
        }

        out.push_back(static_cast<char>((hi << 4) | lo));
    }

    return out;
}

int
QuantityKindAbstract::hexValue(char c)
{
    if (c >= '0' && c <= '9')
    {
        return c - '0';
    }

    if (c >= 'A' && c <= 'F')
    {
        return c - 'A' + 10;
    }

    if (c >= 'a' && c <= 'f')
    {
        return c - 'a' + 10;
    }

    return -1;
}

MCAuto<QuantityKindAbstract>
QuantityKindAbstract::Deserialize(const std::string &s)
{
    const std::string prefix = "QK1|";

    if (s.compare(0, prefix.size(), prefix) != 0)
    {
        THROW_IK_EXCEPTION("Invalid QuantityKind format");
    }

    const std::size_t typeStart = prefix.size();
    const std::size_t typeEnd = s.find('|', typeStart);

    if (typeEnd == std::string::npos)
    {
        THROW_IK_EXCEPTION("Missing type separator");
    }

    const std::string type = s.substr(typeStart, typeEnd - typeStart);

    const std::string payload = s.substr(typeEnd + 1);

    if (type == "UNDEF")
    {
        if (!payload.empty())
        {
            THROW_IK_EXCEPTION("UNDEF must not contain payload");
        }

        return StaticCast<QuantityKindUnDef, QuantityKindAbstract>(QuantityKindUnDef::New());
    }

    if (type == "ENUM")
    {
        return StaticCast<QuantityKindEnum, QuantityKindAbstract>(QuantityKindEnum::New(hexDecode(payload)));
    }

    if (type == "USER")
    {
        return StaticCast<QuantityKindUser, QuantityKindAbstract>(QuantityKindUser::New(hexDecode(payload)));
    }

    THROW_IK_EXCEPTION("Unknown QuantityKind type");
}

namespace
{
template <class IN>
const IN *
IsEqualCommon(const IN *zeThis, const QuantityKindAbstract *other, std::string &reason)
{
    if (!other)
    {
        reason = "other QK is null";
        return nullptr;
    }
    const IN *otherC(dynamic_cast<const IN *>(other));
    if (!otherC)
    {
        std::ostringstream oss;
        oss << "other QK is " << other->getClassName() << " this is " << zeThis->getClassName();
        reason = oss.str();
    }
    return otherC;
}
}  // namespace

std::string
QuantityKindUnDef::repr() const
{
    return std::string("Undef QK");
}

MCAuto<QuantityKindAbstract>
QuantityKindUnDef::clone() const
{
    return MCAuto<QuantityKindAbstract>(QuantityKindUnDef::New().retn());
}

bool
QuantityKindUnDef::isEqual(const QuantityKindAbstract *other, std::string &reason) const
{
    const QuantityKindUnDef *otherC(IsEqualCommon<QuantityKindUnDef>(this, other, reason));
    if (otherC)
        return true;
    else
        return false;
}

MCAuto<QuantityKindUnDef>
QuantityKindUnDef::New()
{
    return MCAuto<QuantityKindUnDef>(new QuantityKindUnDef);
}

std::size_t
QuantityKindUnDef::getHeapMemorySizeWithoutChildren() const
{
    return sizeof(QuantityKindUnDef);
}

std::string
QuantityKindUnDef::serialize() const
{
    return "QK1|UNDEF|";
}

std::string
QuantityKindEnum::repr() const
{
    std::ostringstream oss;
    oss << "Enum QK with value \"" << _value << "\"";
    return oss.str();
}

MCAuto<QuantityKindAbstract>
QuantityKindEnum::clone() const
{
    return MCAuto<QuantityKindAbstract>(QuantityKindEnum::New(this->_value).retn());
}

QuantityKindEnum::QuantityKindEnum(std::string value) : _value(std::move(value))
{
    if (!IsValidValue(_value))
    {
        THROW_IK_EXCEPTION(
            "Invalid QuantityKindEnum value \"" << _value << "\". Allowed are : " << AllowedValuesAsString()
        );
    }
}

bool
QuantityKindEnum::isEqual(const QuantityKindAbstract *other, std::string &reason) const
{
    const QuantityKindEnum *otherC(IsEqualCommon<QuantityKindEnum>(this, other, reason));
    if (otherC)
    {
        if (_value != otherC->_value)
        {
            reason = "QK different";
            return false;
        }
        return true;
    }
    else
        return false;
}

MCAuto<QuantityKindEnum>
QuantityKindEnum::New(std::string value)
{
    return MCAuto(new QuantityKindEnum(value));
}

const std::string &
QuantityKindEnum::value() const
{
    return _value;
}

const std::vector<std::string> &
QuantityKindEnum::AllowedValues()
{
    static const std::vector<std::string> values = {
        "Displacement",
        "PositionVector",
        "Speed",
        "Velocity",
        "Acceleration",
        "Strain",
        "Stress",
        "Pressure",
        "Force",
        "Torque",
        "MomentOfForce",
        "Temperature",
        "HeatFlowRate",
        "HeatFlux",
        "HeatFluxDensity",
        "ThermalConductivity",
        "HeatTransferCoefficient",
        "ThermalConductance",
        "Time",
        "Length",
        "Area",
        "Volume",
        "Mass",
        "Density",
        "Energy",
        "Power",
        "ElectricPotential",
        "ElectricCurrent",
        "Frequency",
        "AngularFrequency",
        "AngularVelocity",
        "SoundPressure",
        "SoundPressureLevel",
        "SoundIntensity",
        "Impedance",
        "AcousticImpedance",
        "DiffusionCoefficient",
        "RelativeHumidity",
        "ParticleFluence",
        "ParticleFluenceRate",
        "Irradiance",
        "StrainEnergyReleaseRate",
        "StressIntensityFactor",
        "DynamicViscosity",
        "KinematicViscosity",
        "SpecificHeatCapacityAtConstantPressure",
        "SpecificEnthalpy",
        "SpecificEnergy",
        "ThermalDiffusivity"
    };
    return values;
}

int
QuantityKindEnum::MCIdOfValue(const std::string &value)
{
    const auto &values(AllowedValues());
    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (values[i] == value)
        {
            return (int)i;
        }
    }
    THROW_IK_EXCEPTION("MCIdOfValue : value \"" << value << "\" not found. Allowed are : " << AllowedValuesAsString());
}

std::string
QuantityKindEnum::AllowedValuesAt(int id)
{
    int len((int)AllowedValues().size());
    if (id < 0 || id >= len)
    {
        THROW_IK_EXCEPTION("AllowedValuesAt : " << id << " requested. Must be in [0," << len << ") !");
    }
    return AllowedValues()[id];
}

std::string
QuantityKindEnum::AllowedValuesAsString()
{
    std::ostringstream oss;

    const auto &values(AllowedValues());

    for (std::size_t i = 0; i < values.size(); ++i)
    {
        if (i != 0)
        {
            oss << ", ";
        }

        oss << values[i];
    }

    return oss.str();
}

bool
QuantityKindEnum::IsValidValue(const std::string &value)
{
    const auto &values(AllowedValues());
    return std::find(values.begin(), values.end(), value) != values.end();
}

std::string
QuantityKindEnum::serialize() const
{
    return "QK1|ENUM|" + hexEncode(_value);
}

std::size_t
QuantityKindEnum::getHeapMemorySizeWithoutChildren() const
{
    return sizeof(QuantityKindEnum) + _value.capacity();
}

std::string
QuantityKindUser::repr() const
{
    std::ostringstream oss;
    oss << "User QK with value \"" << _value << "\"";
    return oss.str();
}

MCAuto<QuantityKindUser>
QuantityKindUser::New(std::string value)
{
    return MCAuto<QuantityKindUser>(new QuantityKindUser(value));
}

MCAuto<QuantityKindAbstract>
QuantityKindUser::clone() const
{
    return MCAuto<QuantityKindAbstract>(QuantityKindUser::New(this->_value).retn());
}

bool
QuantityKindUser::isEqual(const QuantityKindAbstract *other, std::string &reason) const
{
    const QuantityKindUser *otherC(IsEqualCommon<QuantityKindUser>(this, other, reason));
    if (otherC)
    {
        if (_value != otherC->_value)
        {
            reason = "QK different";
            return false;
        }
        return true;
    }
    else
    {
        return false;
    }
}

QuantityKindUser::QuantityKindUser(std::string value) : _value(std::move(value)) {}

const std::string &
QuantityKindUser::value() const
{
    return _value;
}

std::string
QuantityKindUser::serialize() const
{
    return "QK1|USER|" + hexEncode(_value);
}

std::size_t
QuantityKindUser::getHeapMemorySizeWithoutChildren() const
{
    return sizeof(QuantityKindUser) + _value.capacity();
}
