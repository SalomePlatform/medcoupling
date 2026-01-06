// Copyright (C) 2007-2026  CEA, EDF
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

#pragma once

#include "MEDCouplingRefCountObject.hxx"
#include "InterpKernelException.hxx"

#include <vector>
#include <algorithm>

namespace MEDCoupling
{
template <class T>
class MCAuto
{
   public:
    static MCAuto TakeRef(T *ptr)
    {
        MCAuto ret;
        ret.takeRef(ptr);
        return ret;
    }
    MCAuto(const MCAuto &other) : _ptr(nullptr) { referPtr(other._ptr); }
    MCAuto(T *ptr = nullptr) : _ptr(ptr) {}
    ~MCAuto() { destroyPtr(); }
    void checkNotNull() const
    {
        if (!_ptr)
            throw INTERP_KERNEL::Exception("Pointer is nullptr !");
    }
    bool isNull() const { return _ptr == nullptr; }
    bool isNotNull() const { return !isNull(); }
    void nullify()
    {
        destroyPtr();
        _ptr = nullptr;
    }
    bool operator==(const MCAuto &other) const { return _ptr == other._ptr; }
    bool operator==(const T *other) const { return _ptr == other; }
    MCAuto &operator=(const MCAuto &other)
    {
        if (_ptr != other._ptr)
        {
            destroyPtr();
            referPtr(other._ptr);
        }
        return *this;
    }
    MCAuto &operator=(T *ptr)
    {
        if (_ptr != ptr)
        {
            destroyPtr();
            _ptr = ptr;
        }
        return *this;
    }
    void takeRef(T *ptr)
    {
        if (_ptr != ptr)
        {
            destroyPtr();
            referPtr(ptr);
        }
    }
    T *operator->() { return _ptr; }
    const T *operator->() const { return _ptr; }
    T &operator*() { return *_ptr; }
    const T &operator*() const { return *_ptr; }
    operator T *() { return _ptr; }
    operator const T *() const { return _ptr; }
    T *retn()
    {
        if (_ptr)
            _ptr->incrRef();
        return _ptr;
    }
    T *retnConstCast() const
    {
        if (_ptr)
            _ptr->incrRef();
        return _ptr;
    }
    T *iAmATrollConstCast() const { return _ptr; }

   private:
    void referPtr(T *ptr)
    {
        _ptr = ptr;
        if (_ptr)
            _ptr->incrRef();
    }
    void destroyPtr()
    {
        if (_ptr)
            _ptr->decrRef();
    }

   private:
    T *_ptr;
};

template <class T>
std::vector<const T *>
FromVecAutoToVecOfConst(const std::vector<MCAuto<T>> &inputVec)
{
    std::size_t size(inputVec.size());
    std::vector<const T *> ret(size);
    typename std::vector<const T *>::iterator itArrays(ret.begin());
    std::for_each(inputVec.begin(), inputVec.end(), [&itArrays](MCAuto<T> elt) { *itArrays++ = elt; });
    return ret;
}

template <class T>
std::vector<MCAuto<T>>
FromVecToVecAuto(const std::vector<T *> &inputVec)
{
    std::size_t size(inputVec.size());
    std::vector<MCAuto<T>> ret(size);
    typename std::vector<MCAuto<T>>::iterator itArrays(ret.begin());
    std::for_each(
        inputVec.begin(),
        inputVec.end(),
        [&itArrays](T *elt)
        {
            (*itArrays).takeRef(elt);
            itArrays++;
        }
    );
    return ret;
}

template <class T>
std::vector<MCAuto<T>>
FromVecToVecAutoStealRef(const std::vector<T *> &inputVec)
{
    std::size_t size(inputVec.size());
    std::vector<MCAuto<T>> ret(size);
    typename std::vector<MCAuto<T>>::iterator itArrays(ret.begin());
    std::for_each(
        inputVec.begin(),
        inputVec.end(),
        [&itArrays](T *elt)
        {
            *itArrays = elt;
            itArrays++;
        }
    );
    return ret;
}

template <class T>
std::vector<MCAuto<T>>
FromVecConstToVecAuto(const std::vector<const T *> &inputVec)
{
    std::size_t size(inputVec.size());
    std::vector<MCAuto<T>> ret(size);
    typename std::vector<MCAuto<T>>::iterator itArrays(ret.begin());
    std::for_each(
        inputVec.begin(),
        inputVec.end(),
        [&itArrays](const T *elt)
        {
            (*itArrays).takeRef(const_cast<T *>(elt));
            itArrays++;
        }
    );
    return ret;
}

template <class T, class U>
typename MEDCoupling::MCAuto<U>
DynamicCast(typename MEDCoupling::MCAuto<T> &autoSubPtr) noexcept(true)
{
    T *subPtr(autoSubPtr);
    U *ptr(dynamic_cast<U *>(subPtr));
    typename MEDCoupling::MCAuto<U> ret(ptr);
    if (ptr)
        ptr->incrRef();
    return ret;
}

template <class T, class U>
typename MEDCoupling::MCAuto<U>
DynamicCastSafe(typename MEDCoupling::MCAuto<T> &autoSubPtr)
{
    T *subPtr(autoSubPtr);
    U *ptr(dynamic_cast<U *>(subPtr));
    if (subPtr && !ptr)
        throw INTERP_KERNEL::Exception("DynamicCastSafe : U is not a subtype of T !");
    typename MEDCoupling::MCAuto<U> ret(ptr);
    if (ptr)
        ptr->incrRef();
    return ret;
}

template <class T>
typename std::vector<const T *>
ToConstVect(const typename std::vector<MCAuto<T>> &vec)
{
    std::size_t sz(vec.size());
    std::vector<const T *> ret(sz);
    for (std::size_t i = 0; i < sz; i++) ret[i] = (const T *)vec[i];
    return ret;
}

template <class T>
class MCConstAuto
{
   public:
    MCConstAuto(const MCConstAuto &other) : _ptr(0) { referPtr(other._ptr); }
    MCConstAuto(const typename MEDCoupling::MCAuto<T> &other) : _ptr(0) { referPtr((const T *)other); }
    MCConstAuto(const T *ptr = 0) : _ptr(ptr) {}
    ~MCConstAuto() { destroyPtr(); }
    bool isNull() const { return _ptr == 0; }
    bool isNotNull() const { return !isNull(); }
    void nullify()
    {
        destroyPtr();
        _ptr = 0;
    }
    bool operator==(const MCConstAuto &other) const { return _ptr == other._ptr; }
    bool operator==(const T *other) const { return _ptr == other; }
    MCConstAuto &operator=(const MCConstAuto &other)
    {
        if (_ptr != other._ptr)
        {
            destroyPtr();
            referPtr(other._ptr);
        }
        return *this;
    }
    MCConstAuto &operator=(const typename MEDCoupling::MCAuto<T> &other)
    {
        if (_ptr != (const T *)other)
        {
            destroyPtr();
            referPtr((const T *)other);
        }
        return *this;
    }
    MCConstAuto &operator=(const T *ptr)
    {
        if (_ptr != ptr)
        {
            destroyPtr();
            _ptr = ptr;
        }
        return *this;
    }
    void takeRef(const T *ptr)
    {
        if (_ptr != ptr)
        {
            destroyPtr();
            _ptr = ptr;
            if (_ptr)
                _ptr->incrRef();
        }
    }
    const T *operator->() { return _ptr; }
    const T *operator->() const { return _ptr; }
    const T &operator*() { return *_ptr; }
    const T &operator*() const { return *_ptr; }
    operator const T *() const { return _ptr; }
    T *shameOnMeConstCast() const { return const_cast<T *>(_ptr); }

   private:
    void referPtr(const T *ptr)
    {
        _ptr = ptr;
        if (_ptr)
            _ptr->incrRef();
    }
    void destroyPtr()
    {
        if (_ptr)
            _ptr->decrRef();
    }

   private:
    const T *_ptr;
};

template <class T>
std::vector<const T *>
FromVecAutoToVecOfConst(const std::vector<MCConstAuto<T>> &inputVec)
{
    std::size_t size(inputVec.size());
    std::vector<const T *> ret(size);
    typename std::vector<const T *>::iterator itArrays(ret.begin());
    std::for_each(inputVec.begin(), inputVec.end(), [&itArrays](MCConstAuto<T> elt) { *itArrays++ = elt; });
    return ret;
}

}  // namespace MEDCoupling
