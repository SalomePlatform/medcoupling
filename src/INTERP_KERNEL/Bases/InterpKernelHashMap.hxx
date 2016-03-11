// Copyright (C) 2001-2016  CEA/DEN, EDF R&D
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

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/*
 * Copyright (c) 1996
 * Silicon Graphics Computer Systems, Inc.
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Silicon Graphics makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 *
 * Copyright (c) 1994
 * Hewlett-Packard Company
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies and
 * that both that copyright notice and this permission notice appear
 * in supporting documentation.  Hewlett-Packard Company makes no
 * representations about the suitability of this software for any
 * purpose.  It is provided "as is" without express or implied warranty.
 *
 */
#ifndef __INTERPKERNELHASHMAP__
#define __INTERPKERNELHASHMAP__

#include "InterpKernelStlExt.hxx"
#include "InterpKernelHashTable.hxx"

namespace INTERP_KERNEL
{
  template<class _Key, class _Tp, class _HashFn = hash<_Key>,
           class _EqualKey = std::equal_to<_Key>, class _Alloc = std::allocator<_Tp> >
  class HashMap
  {
  private:
    typedef hashtable<std::pair<const _Key, _Tp>,_Key, _HashFn,
                      STLEXT::Select1st<std::pair<const _Key, _Tp> >,
                      _EqualKey, _Alloc> _Ht;

    _Ht _M_ht;
    
  public:
    typedef typename _Ht::key_type key_type;
    typedef _Tp data_type;
    typedef _Tp mapped_type;
    typedef typename _Ht::value_type value_type;
    typedef typename _Ht::hasher hasher;
    typedef typename _Ht::key_equal key_equal;
    
    typedef typename _Ht::size_type size_type;
    typedef typename _Ht::difference_type difference_type;
    typedef typename _Ht::pointer pointer;
    typedef typename _Ht::const_pointer const_pointer;
    typedef typename _Ht::reference reference;
    typedef typename _Ht::const_reference const_reference;
    
    typedef typename _Ht::iterator iterator;
    typedef typename _Ht::const_iterator const_iterator;
    
    typedef typename _Ht::allocator_type allocator_type;
      
    hasher hash_funct() const { return _M_ht.hash_funct(); }

    key_equal key_eq() const { return _M_ht.key_eq(); }
    
    allocator_type get_allocator() const { return _M_ht.get_allocator(); }

    HashMap() : _M_ht(100, hasher(), key_equal(), allocator_type()) {}
  
    explicit HashMap(size_type __n) : _M_ht(__n, hasher(), key_equal(), allocator_type()) {}

    HashMap(size_type __n, const hasher& __hf) : _M_ht(__n, __hf, key_equal(), allocator_type()) {}

    HashMap(size_type __n, const hasher& __hf, const key_equal& __eql,
            const allocator_type& __a = allocator_type()) : _M_ht(__n, __hf, __eql, __a) {}
    
    template<class _InputIterator>
    HashMap(_InputIterator __f, _InputIterator __l) : _M_ht(100, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
    
    template<class _InputIterator>
    HashMap(_InputIterator __f, _InputIterator __l, size_type __n) : _M_ht(__n, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }

    template<class _InputIterator>
    HashMap(_InputIterator __f, _InputIterator __l, size_type __n, const hasher& __hf)
      : _M_ht(__n, __hf, key_equal(), allocator_type())
    { _M_ht.insert_unique(__f, __l); }
    
    template<class _InputIterator>
    HashMap(_InputIterator __f, _InputIterator __l, size_type __n,
            const hasher& __hf, const key_equal& __eql,
            const allocator_type& __a = allocator_type()) : _M_ht(__n, __hf, __eql, __a)
    { _M_ht.insert_unique(__f, __l); }
    
    size_type size() const { return _M_ht.size(); }
    
    size_type max_size() const { return _M_ht.max_size(); }
    
    bool empty() const { return _M_ht.empty(); }
    
    void swap(HashMap& __hs) { _M_ht.swap(__hs._M_ht); }
    
    template<class _K1, class _T1, class _HF, class _EqK, class _Al>
    friend bool operator== (const HashMap<_K1, _T1, _HF, _EqK, _Al>&,
                            const HashMap<_K1, _T1, _HF, _EqK, _Al>&);
    
    iterator begin() { return _M_ht.begin(); }
    
    iterator end() { return _M_ht.end(); }
    
    const_iterator begin() const { return _M_ht.begin(); }
    
    const_iterator end() const { return _M_ht.end(); }
    
    std::pair<iterator, bool> insert(const value_type& __obj) { return _M_ht.insert_unique(__obj); }
    
    template<class _InputIterator>
    void insert(_InputIterator __f, _InputIterator __l) { _M_ht.insert_unique(__f, __l); }
    
    std::pair<iterator, bool>
    insert_noresize(const value_type& __obj) { return _M_ht.insert_unique_noresize(__obj); }
    
    iterator find(const key_type& __key) { return _M_ht.find(__key); }
    
    const_iterator find(const key_type& __key) const { return _M_ht.find(__key); }
    
    _Tp& operator[](const key_type& __key) { return _M_ht.find_or_insert(value_type(__key, _Tp())).second; }
    
    size_type count(const key_type& __key) const { return _M_ht.count(__key); }
    
    std::pair<iterator, iterator> equal_range(const key_type& __key) { return _M_ht.equal_range(__key); }
    
    std::pair<const_iterator, const_iterator> equal_range(const key_type& __key) const { return _M_ht.equal_range(__key); }
    
    size_type erase(const key_type& __key) { return _M_ht.erase(__key); }
    
    void erase(iterator __it) { _M_ht.erase(__it); }
    
    void erase(iterator __f, iterator __l) { _M_ht.erase(__f, __l); }

    void clear() { _M_ht.clear(); }

    void resize(size_type __hint) { _M_ht.resize(__hint); }
    
    size_type bucket_count() const { return _M_ht.bucket_count(); }

    size_type max_bucket_count() const { return _M_ht.max_bucket_count(); }
    
    size_type elems_in_bucket(size_type __n) const { return _M_ht.elems_in_bucket(__n); }
  };
  
  template<class _Key, class _Tp, class _HashFn, class _EqlKey, class _Alloc>
  inline bool operator==(const HashMap<_Key, _Tp, _HashFn, _EqlKey, _Alloc>& __hm1,
                         const HashMap<_Key, _Tp, _HashFn, _EqlKey, _Alloc>& __hm2)
  { return __hm1._M_ht == __hm2._M_ht; }
  
  template<class _Key, class _Tp, class _HashFn, class _EqlKey, class _Alloc>
  inline bool operator!=(const HashMap<_Key, _Tp, _HashFn, _EqlKey, _Alloc>& __hm1,
                         const HashMap<_Key, _Tp, _HashFn, _EqlKey, _Alloc>& __hm2)
  { return !(__hm1 == __hm2); }
  
  template<class _Key, class _Tp, class _HashFn, class _EqlKey, class _Alloc>
  inline void swap(HashMap<_Key, _Tp, _HashFn, _EqlKey, _Alloc>& __hm1,
                   HashMap<_Key, _Tp, _HashFn, _EqlKey, _Alloc>& __hm2)
  { __hm1.swap(__hm2); }

  template<class _Key, class _Tp,
           class _HashFn = hash<_Key>,
           class _EqualKey = std::equal_to<_Key>,
           class _Alloc = std::allocator<_Tp> >
  class HashMultiMap
  {
  private:
    typedef hashtable<std::pair<const _Key, _Tp>, _Key, _HashFn,
                      STLEXT::Select1st<std::pair<const _Key, _Tp> >, _EqualKey, _Alloc>
    _Ht;
    _Ht _M_ht;
  public:
    typedef typename _Ht::key_type key_type;
    typedef _Tp data_type;
    typedef _Tp mapped_type;
    typedef typename _Ht::value_type value_type;
    typedef typename _Ht::hasher hasher;
    typedef typename _Ht::key_equal key_equal;
    
    typedef typename _Ht::size_type size_type;
    typedef typename _Ht::difference_type difference_type;
    typedef typename _Ht::pointer pointer;
    typedef typename _Ht::const_pointer const_pointer;
    typedef typename _Ht::reference reference;
    typedef typename _Ht::const_reference const_reference;
    
    typedef typename _Ht::iterator iterator;
    typedef typename _Ht::const_iterator const_iterator;
    
    typedef typename _Ht::allocator_type allocator_type;
  
    hasher hash_funct() const { return _M_ht.hash_funct(); }
    
    key_equal key_eq() const { return _M_ht.key_eq(); }
    
    allocator_type get_allocator() const { return _M_ht.get_allocator(); }
    
    HashMultiMap() : _M_ht(100, hasher(), key_equal(), allocator_type()) { }
    
    explicit HashMultiMap(size_type __n) : _M_ht(__n, hasher(), key_equal(), allocator_type()) {}
    
    HashMultiMap(size_type __n, const hasher& __hf) : _M_ht(__n, __hf, key_equal(), allocator_type()) {}
    
    HashMultiMap(size_type __n, const hasher& __hf, const key_equal& __eql,
                 const allocator_type& __a = allocator_type()) : _M_ht(__n, __hf, __eql, __a) {}
    
    template<class _InputIterator>
    HashMultiMap(_InputIterator __f, _InputIterator __l) : _M_ht(100, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
    
    template<class _InputIterator>
    HashMultiMap(_InputIterator __f, _InputIterator __l, size_type __n) : _M_ht(__n, hasher(), key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
    
    template<class _InputIterator>
    HashMultiMap(_InputIterator __f, _InputIterator __l, size_type __n, const hasher& __hf)
      : _M_ht(__n, __hf, key_equal(), allocator_type())
    { _M_ht.insert_equal(__f, __l); }
    
    template<class _InputIterator>
    HashMultiMap(_InputIterator __f, _InputIterator __l, size_type __n,
                 const hasher& __hf, const key_equal& __eql,
                 const allocator_type& __a = allocator_type())
      : _M_ht(__n, __hf, __eql, __a)
    { _M_ht.insert_equal(__f, __l); }
    
    size_type size() const { return _M_ht.size(); }
    
    size_type max_size() const { return _M_ht.max_size(); }
    
    bool empty() const { return _M_ht.empty(); }
    
    void swap(HashMultiMap& __hs) { _M_ht.swap(__hs._M_ht); }
    
    template<class _K1, class _T1, class _HF, class _EqK, class _Al>
    friend bool operator==(const HashMultiMap<_K1, _T1, _HF, _EqK, _Al>&,
                           const HashMultiMap<_K1, _T1, _HF, _EqK, _Al>&);
    
    iterator begin() { return _M_ht.begin(); }
    
    iterator end() { return _M_ht.end(); }
    
    const_iterator begin() const { return _M_ht.begin(); }
    
    const_iterator end() const { return _M_ht.end(); }
    
    iterator insert(const value_type& __obj) { return _M_ht.insert_equal(__obj); }
    
    template<class _InputIterator>
    void insert(_InputIterator __f, _InputIterator __l) { _M_ht.insert_equal(__f,__l); }
    
    iterator insert_noresize(const value_type& __obj) { return _M_ht.insert_equal_noresize(__obj); }
    
    iterator find(const key_type& __key) { return _M_ht.find(__key); }
    
    const_iterator find(const key_type& __key) const { return _M_ht.find(__key); }
    
    size_type count(const key_type& __key) const { return _M_ht.count(__key); }
    
    std::pair<iterator, iterator> equal_range(const key_type& __key) { return _M_ht.equal_range(__key); }
    
    std::pair<const_iterator, const_iterator> equal_range(const key_type& __key) const { return _M_ht.equal_range(__key); }
    
    size_type erase(const key_type& __key) { return _M_ht.erase(__key); }
    
    void erase(iterator __it) { _M_ht.erase(__it); }
    
    void erase(iterator __f, iterator __l) { _M_ht.erase(__f, __l); }
    
    void clear() { _M_ht.clear(); }
    
    void resize(size_type __hint) { _M_ht.resize(__hint); }
    
    size_type bucket_count() const { return _M_ht.bucket_count(); }
    
    size_type max_bucket_count() const { return _M_ht.max_bucket_count(); }
    
    size_type elems_in_bucket(size_type __n) const { return _M_ht.elems_in_bucket(__n); }
  };
  
  template<class _Key, class _Tp, class _HF, class _EqKey, class _Alloc>
  inline bool operator==(const HashMultiMap<_Key, _Tp, _HF, _EqKey, _Alloc>& __hm1,
                         const HashMultiMap<_Key, _Tp, _HF, _EqKey, _Alloc>& __hm2)
  { return __hm1._M_ht == __hm2._M_ht; }
  
  template<class _Key, class _Tp, class _HF, class _EqKey, class _Alloc>
  inline bool operator!=(const HashMultiMap<_Key, _Tp, _HF, _EqKey, _Alloc>& __hm1,
                         const HashMultiMap<_Key, _Tp, _HF, _EqKey, _Alloc>& __hm2)
  { return !(__hm1 == __hm2); }
  
  template<class _Key, class _Tp, class _HashFn, class _EqlKey, class _Alloc>
  inline void swap(HashMultiMap<_Key, _Tp, _HashFn, _EqlKey, _Alloc>& __hm1,
                   HashMultiMap<_Key, _Tp, _HashFn, _EqlKey, _Alloc>& __hm2)
  { __hm1.swap(__hm2); }
  
}

namespace std
{
  // Specialization of insert_iterator so that it will work for HashMap
  // and HashMultiMap.
  template<class _Key, class _Tp, class _HashFn,  class _EqKey, class _Alloc>
  class insert_iterator<INTERP_KERNEL::HashMap<_Key, _Tp, _HashFn, 
                                               _EqKey, _Alloc> >
  {
  protected:
    typedef INTERP_KERNEL::HashMap<_Key, _Tp, _HashFn, _EqKey, _Alloc>
    _Container;
    _Container* container;
  public:
    typedef _Container          container_type;
    typedef output_iterator_tag iterator_category;
    typedef void                value_type;
    typedef void                difference_type;
    typedef void                pointer;
    typedef void                reference;
      
    insert_iterator(_Container& __x) : container(&__x) {}
    
    insert_iterator(_Container& __x, typename _Container::iterator) : container(&__x) {}
    
    insert_iterator<_Container>& operator=(const typename _Container::value_type& __value__)
    {
      container->insert(__value__);
      return *this;
    }
    
    insert_iterator<_Container>& operator*() { return *this; }
    
    insert_iterator<_Container>& operator++() { return *this; }

    insert_iterator<_Container>& operator++(int) { return *this; }
  };

  template<class _Key, class _Tp, class _HashFn,  class _EqKey, class _Alloc>
  class insert_iterator<INTERP_KERNEL::HashMultiMap<_Key, _Tp, _HashFn,
                                                    _EqKey, _Alloc> >
  {
  protected:
    typedef INTERP_KERNEL::HashMultiMap<_Key, _Tp, _HashFn, _EqKey, _Alloc>
    _Container;
    _Container* container;
    typename _Container::iterator iter;
    
  public:
    typedef _Container          container_type;
    typedef output_iterator_tag iterator_category;
    typedef void                value_type;
    typedef void                difference_type;
    typedef void                pointer;
    typedef void                reference;
    
    insert_iterator(_Container& __x) : container(&__x) {}

    insert_iterator(_Container& __x, typename _Container::iterator) : container(&__x) {}

    insert_iterator<_Container>& operator=(const typename _Container::value_type& __value__)
    {
      container->insert(__value__);
      return *this;
    }

    insert_iterator<_Container>& operator*() { return *this; }

    insert_iterator<_Container>& operator++() { return *this; }

    insert_iterator<_Container>& operator++(int) { return *this; }
  };
}

#endif
