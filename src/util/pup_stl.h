/*
Pup routines for STL classes.

After including this header, you can parameter-marshall
a variable consisting of STL containers such as vectors,
lists, maps, strings, or pairs.

This includes variables of type "std::list<int>", or even
"std::map<double, std::vector<std::string> >".

NOT included are the rarer types like valarray or slice.

Orion Sky Lawlor, olawlor@acm.org, 7/22/2002
*/
#ifndef _UIUC_CHARM_PUP_STL_H
#define _UIUC_CHARM_PUP_STL_H

#include <conv-config.h>

/*It's kind of annoying that we have to drag all these headers in
  just so the std:: parameter declarations will compile.
 */
#include <algorithm>
#include <array>
#include <set>
#include <vector>
#include <deque>
#include <forward_list>
#include <list>
#include <map>
#include <memory>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <complex>
#include <utility> /*for std::pair*/
#include "pup.h"

#include <cstddef>

namespace PUP {
  /*************** Simple classes ***************/
  // Non-const version is required for puping std::pair
  template <class A,class B>
  inline void operator|(er &p,typename std::pair<A,B> &v);
  template <class A,class B>
  inline void operator|(er &p,typename std::pair<const A,B> &v);
  template <class T>
  inline void operator|(er &p,std::complex<T> &v);
  template <class T>
  inline void operator|(er &p, std::unique_ptr<T, std::default_delete<T>> &ptr);
  template <class charType>
  inline void operator|(er &p,typename std::basic_string<charType> &v);
  inline void operator|(er &p,std::string &v);
  template <class container>
  inline size_t PUP_stl_container_size(er &p,container &c);
  template <class container, class dtype>
  inline void PUP_stl_container_items(er &p, container &c, size_t nElem);
  template <> inline void PUP_stl_container_items<std::vector<bool>,bool>(er &p, std::vector<bool> &c, size_t nElem);
  template <class container,class dtype>
  inline void PUP_stl_container(er &p,container &c);
  template <class container,class dtype>
  inline void PUP_stl_map(er &p,container &c);
  template <class T>
  inline void operator|(er &p,typename std::vector<T> &v);
  template <class T>
  inline void operator|(er &p,typename std::deque<T> &d);
  template <class T>
  inline void operator|(er &p,typename std::list<T> &v);
  template <class T>
  inline void operator|(er &p,typename std::forward_list<T> &fl);
  template <class V,class T,class Cmp>
  inline void operator|(er &p,typename std::map<V,T,Cmp> &m);
  template <class V,class T,class Cmp>
  inline void operator|(er &p,typename std::multimap<V,T,Cmp> &m);
  template <class T>
  inline void operator|(er &p,typename std::set<T> &m);
  template <class T,class Cmp>
  inline void operator|(er &p,typename std::multiset<T,Cmp> &m);
  template <> inline void operator|(er &p,std::vector<bool> &v);

  template <class A,class B>
  inline void operator|(er &p,typename std::pair<A,B> &v)
  {
    p.syncComment(sync_index);
    p|v.first;
    p.syncComment(sync_item);
    p|v.second;
  }
  // Const version is required for puping std::map
  template <class A,class B> 
  inline void operator|(er &p,typename std::pair<const A,B> &v)
  {
    p.syncComment(sync_index);
    p|*(A *)&v.first; /* cast away constness on A */
    p.syncComment(sync_item);
    p|v.second;
  }
  template <class T>
  inline void operator|(er &p,std::complex<T> &v)
  {
    T re=v.real(), im=v.imag();
    p|re; p|im;
    v=std::complex<T>(re,im);
  }
  template <class charType> 
  inline void operator|(er &p,typename std::basic_string<charType> &v)
  {
    size_t nChar=v.length();
    p|nChar;
    if (p.isUnpacking()) { //Unpack to temporary buffer
      charType *buf=new charType[nChar];
      p(buf,nChar);
      v=std::basic_string<charType>(buf,nChar);
      delete[] buf;
    }
    else /*packing*/ { //Do packing in-place from data
      //Have to cast away constness here
      p((charType *)v.data(),nChar);
    }
  }
  inline void operator|(er &p,std::string &v)
  {
    p.syncComment(sync_begin_object,"std::string");
    size_t nChar=v.length();
    p|nChar;
    if (p.isUnpacking()) { //Unpack to temporary buffer
      char *buf=new char[nChar];
      p(buf,nChar);
      v=std::basic_string<char>(buf,nChar);
      delete[] buf;
    }
    else /*packing*/ { //Do packing in-place from data
      //Have to cast away constness here
      p((char *)v.data(),nChar);
    }
    p.syncComment(sync_end_object);
  }

  /**************** Containers *****************/

  template <class container>
  inline void reserve_if_applicable(container &c, size_t nElem)
  {
    c.clear();
    c.reserve(nElem);
  }
  template <class dtype>
  inline void reserve_if_applicable(std::deque<dtype> &c, size_t nElem)
  {
    c.clear();
  }
  template <class dtype>
  inline void reserve_if_applicable(std::list<dtype> &c, size_t nElem)
  {
    c.clear();
  }
  template <class dtype>
  inline void reserve_if_applicable(std::set<dtype> &c, size_t nElem)
  {
    c.clear();
  }
  template <class dtype>
  inline void reserve_if_applicable(std::multiset<dtype> &c, size_t nElem)
  {
    c.clear();
  }
  template <class K, class V>
  inline void reserve_if_applicable(std::map<K, V> &c, size_t nElem)
  {
    c.clear();
  }
  template <class K, class V>
  inline void reserve_if_applicable(std::multimap<K, V> &c, size_t nElem)
  {
    c.clear();
  }

  template <class container, class... Args>
  inline void emplace(container &c, Args&&... args)
  {
    c.emplace(std::forward<Args>(args)...);
  }
  template <class dtype, class... Args>
  inline void emplace(std::vector<dtype> &c, Args&&... args)
  {
    c.emplace_back(std::forward<Args>(args)...);
  }
  template <class dtype, class... Args>
  inline void emplace(std::deque<dtype> &c, Args&&... args)
  {
    c.emplace_back(std::forward<Args>(args)...);
  }
  template <class dtype, class... Args>
  inline void emplace(std::list<dtype> &c, Args&&... args)
  {
    c.emplace_back(std::forward<Args>(args)...);
  }

  //Impl. util: pup the length of a container
  template <class container>
  inline size_t PUP_stl_container_size(er &p,container &c) {
    size_t nElem=c.size();
    p|nElem;
    return nElem; 
  }

  //Impl. util: pup each current item of a container (no allocation)
  template <class container, class dtype>
  inline void PUP_stl_container_items(er &p, container &c, size_t nElem)
  {
    if (p.isUnpacking())
    {
      reserve_if_applicable(c, nElem);
      for (size_t i = 0; i < nElem; ++i)
      {
        p.syncComment(sync_item);
        detail::TemporaryObjectHolder<dtype> n;
        p|n;
        emplace(c, std::move(n.t));
      }
    }
    else
    {
      for (typename container::iterator it=c.begin(); it!=c.end(); ++it)
      {
        p.syncComment(sync_item);
        // Cast away the constness (needed for std::set)
        p|*(dtype *)&(*it);
      }
    }
  }

  // Specialized to work with vector<bool>
  template<>
  inline void PUP_stl_container_items<std::vector<bool>, bool>(er &p, std::vector<bool> &c, size_t nElem)
  {
    if (p.isUnpacking())
      c.resize(nElem);

    std::deque<bool> q(c.begin(), c.end());
    
    for (std::deque<bool>::iterator it = q.begin(); it != q.end(); it++)
    {
      p.syncComment(sync_item);
      p|*it;
    }
  }

  template <class container,class dtype>
  inline void PUP_stl_container(er &p,container &c) {
    p.syncComment(sync_begin_array);
    size_t nElem=PUP_stl_container_size(p,c);
    PUP_stl_container_items<container, dtype>(p, c, nElem);
    p.syncComment(sync_end_array);
  }

  // forward_list does not have: .size(), .emplace(), .emplace_back()
  template <class dtype>
  inline void PUP_stl_forward_list(er &p,std::forward_list<dtype> &c) {
    p.syncComment(sync_begin_array);
    size_t nElem;
    if (p.isUnpacking())
    {
      p | nElem;
      auto iter = c.before_begin();
      for (size_t i = 0; i < nElem; ++i)
      {
        p.syncComment(sync_item);
        detail::TemporaryObjectHolder<dtype> n;
        p|n;
        iter = c.emplace_after(iter, std::move(n.t));
      }
    }
    else
    {
      nElem = 0;
      for (auto& n: c)
      {
        ++nElem;
      }
      p | nElem;
      for (auto& n : c)
      {
        p.syncComment(sync_item);
        p | n;
      }
    }
    p.syncComment(sync_end_array);
  }

  template <class container, class K, class V>
  inline void PUP_stl_map(er &p,container &c) {
    p.syncComment(sync_begin_list);
    size_t nElem=PUP_stl_container_size(p,c);
    if (p.isUnpacking())
      { //Unpacking: Extract each element and insert:
        reserve_if_applicable(c, nElem);
        for (size_t i=0;i<nElem;i++)
        {
          detail::TemporaryObjectHolder<K> k;
          detail::TemporaryObjectHolder<V> v;

          // keep in sync with std::pair
          p.syncComment(sync_index);
          p | k;
          p.syncComment(sync_item);
          p | v;

          c.emplace(std::piecewise_construct, std::forward_as_tuple(std::move(k.t)), std::forward_as_tuple(std::move(v.t)));
        }
      }
    else
    {
      for (auto& kv : c)
      {
        p | kv;
      }
    }
    p.syncComment(sync_end_list);
  }

  template <class T>
  inline void operator|(er &p, typename std::vector<T> &v) {
    if (PUP::as_bytes<T>::value) {
      size_t nElem = PUP_stl_container_size(p, v);
      if (p.isUnpacking()) {
        v.resize(nElem);
        v.shrink_to_fit();
      }
      PUParray(p, v.data(), nElem);
    } else {
      PUP_stl_container<std::vector<T>, T>(p, v);
    }
  }

  template <class T>
  inline void operator|(er &p,typename std::deque<T> &d)
  { PUP_stl_container<std::deque<T>,T>(p,d); }
  template <class T> 
  inline void operator|(er &p,typename std::list<T> &v)
  { PUP_stl_container<std::list<T>,T>(p,v); }
  template <class T>
  inline void operator|(er &p,typename std::forward_list<T> &fl)
  { PUP_stl_forward_list<T>(p,fl); }

  template <class V,class T,class Cmp> 
  inline void operator|(er &p,typename std::map<V,T,Cmp> &m)
  { PUP_stl_map<std::map<V,T,Cmp>,V,T >(p,m); }
  template <class V,class T,class Cmp> 
  inline void operator|(er &p,typename std::multimap<V,T,Cmp> &m)
  { PUP_stl_map<std::multimap<V,T,Cmp>,V,T >(p,m); }
  /// \warning This does not work with custom hash functions that have state
  template <class V,class T,class Cmp>
  inline void operator|(er &p,typename std::unordered_map<V,T,Cmp> &m)
  { PUP_stl_map<std::unordered_map<V,T,Cmp>,V,T >(p,m); }
  template <class V,class T,class Cmp>
  inline void operator|(er &p,typename std::unordered_multimap<V,T,Cmp> &m)
  { PUP_stl_map<std::unordered_multimap<V,T,Cmp>,V,T >(p,m); }

  template <class T>
  inline void operator|(er &p,typename std::set<T> &m)
  { PUP_stl_container<std::set<T>,T >(p,m); }
  template <class T,class Cmp>
  inline void operator|(er &p,typename std::multiset<T,Cmp> &m)
  { PUP_stl_container<std::multiset<T,Cmp>,T >(p,m); }
  template <class T>
  inline void operator|(er &p,typename std::unordered_set<T> &m)
  { PUP_stl_container<std::unordered_set<T>,T >(p,m); }
  template <class T,class Cmp>
  inline void operator|(er &p,typename std::unordered_multiset<T,Cmp> &m)
  { PUP_stl_container<std::unordered_multiset<T,Cmp>,T >(p,m); }

  // Specialized to work with vector<bool>, which doesn't
  // have data() or shrink_to_fit() members
  template <>
  inline void operator|(er &p,std::vector<bool> &v) {
    PUP_stl_container<std::vector<bool>, bool>(p, v);
  }

// Distributed under the MIT License.
// The following allows for pupping STL structures with pointers to abstract
// base classes. Requires is used in place of enable_if_t to enforce
// requirements on template parameters for the following PUP methods.
template <bool B>
struct requires_impl {
  using template_error_type_failed_to_meet_requirements_on_template_parameters
      = std::nullptr_t;
};

template <>
struct requires_impl<false> {};

template <bool B>
using Requires = typename requires_impl<
    B>::template_error_type_failed_to_meet_requirements_on_template_parameters;

  template <typename T, std::size_t N,
            Requires<!PUP::as_bytes<T>::value> = nullptr>
  inline void pup(PUP::er& p, std::array<T, N>& a) {
    std::for_each(a.begin(), a.end(), [&p](T& t) { p | t; });
  }

  template <typename T, std::size_t N,
            Requires<PUP::as_bytes<T>::value> = nullptr>
  inline void pup(PUP::er& p, std::array<T, N>& a) {
    PUParray(p, a.data(), N);
  }

  template <typename T, std::size_t N>
  inline void operator|(er& p, std::array<T, N>& a) {
    pup(p, a);
  }

  template <typename T, Requires<std::is_enum<T>::value> = nullptr>
  inline void operator|(PUP::er& p, T& s) {
    pup_bytes(&p, static_cast<void*>(&s), sizeof(T));
  }

  template <size_t N = 0, typename... Args,
            Requires<0 == sizeof...(Args)> = nullptr>
  void pup_tuple_impl(PUP::er& /* p */, std::tuple<Args...>& /* t */) {
  }

  template <size_t N = 0, typename... Args,
            Requires<(0 < sizeof...(Args) && 0 == N)> = nullptr>
  void pup_tuple_impl(PUP::er& p, std::tuple<Args...>& t) {
    p | std::get<N>(t);
  }

  template <size_t N, typename... Args,
            Requires<(sizeof...(Args) > 0 && N > 0)> = nullptr>
  void pup_tuple_impl(PUP::er& p, std::tuple<Args...>& t) {
    p | std::get<N>(t);
    pup_tuple_impl<N - 1>(p, t);
  }

  template <typename... Args>
  inline void pup(PUP::er& p, std::tuple<Args...>& t) {
    pup_tuple_impl<sizeof...(Args) - 1>(p, t);
  }

  template <typename... Args>
  inline void operator|(PUP::er& p, std::tuple<Args...>& t) {
    pup(p, t);
  }

  template <typename T,
            Requires<!std::is_base_of<PUP::able, T>::value> = nullptr>
  inline void pup(PUP::er& p, std::unique_ptr<T>& t) {
    bool is_nullptr = nullptr == t;
    p | is_nullptr;
    if (!is_nullptr) {
      T* t1;
      if (p.isUnpacking()) {
        t1 = new T;
      } else {
        t1 = t.get();
      }
      p | *t1;
      if (p.isUnpacking()) {
        t.reset(t1);
      }
    }
  }

  template <typename T, Requires<std::is_base_of<PUP::able, T>::value> = nullptr>
  inline void pup(PUP::er& p, std::unique_ptr<T>& t) {
    T* t1 = nullptr;
    if (p.isUnpacking()) {
      p | t1;
      t = std::unique_ptr<T>(t1);
    } else {
      t1 = t.get();
      p | t1;
    }
  }

  template <typename T>
  inline void operator|(PUP::er& p, std::unique_ptr<T>& t) {
    pup(p, t);
  }

} // end of namespace PUP

#endif
