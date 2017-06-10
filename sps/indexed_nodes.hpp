/**
 * @file   indexed_nodes.hpp
 * @author Jens Munk Hansen <jens.munk.hansen@gmail.com>
 * @date   Fri Jan 10 21:58:00 2014
 *
 * @brief Indexed nodes
 *
 *
 */

/*
 *  This file is part of SOFUS.
 *
 *  SOFUS is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  SOFUS is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with SOFUS.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

// TODO: Make specialization of garbage_collected type, indexed and indexed nodes

//#include <sps/env.h>
#include <sps/cstdio>
#include <sps/cstdlib>
#include <cstdint>

#include <vector>
#include <algorithm>
#include <stdexcept>

#include <typeinfo>
#include <cstring>

#if (defined(__GNUC__) && __GNUC__ >= 3)                                \
  || (defined(_WIN32))                                                  \
  || (defined(linux) && defined(__INTEL_COMPILER) && defined(__ICC))
# define SPS_TYPE_ID_NAME
#endif

namespace sps {
  class base_node {
  public:
    base_node() : children(NULL) {}
  protected:
    virtual ~base_node() {}

    /// Two-dimensional array of nodes (multiple instances of different types)
    base_node** children;
  };

  /// Forward declaration of indexed_auto_ptr_collector
  template<typename T> class indexed_node_collector;

  // Suggestion, include template int specifying how many can co-exist
  template <typename T>
  class indexed_node : public base_node {
  public:
    indexed_node() : type("None"), t(NULL) {}

    /**
     * Ctor
     *
     * @param ptr
     */
#ifdef SPS_TYPE_ID_NAME
    indexed_node(T*& ptr) : type(typeid(T).name()), t(ptr)
#else
    indexed_node(T*& ptr) : type(&typeid(T)), t(ptr)
#endif
    {
      signature = this;
      indexed_node_collector<T>::register_handle(this);
      ptr = 0;
    }

    /**
     * Dtor
     *
     */
    ~indexed_node()
    {
      delete t;        // destroy object
      signature= NULL; // destroy signature
    }

    /**
     * Convert index to an node<T>.
     *
     * @param index
     *
     * @return
     */
    static indexed_node* from_index( const uint64_t index ) throw(std::runtime_error);

    /**
     * Convert node<T> to an index.
     *
     *
     * @return
     */
    uint64_t to_index();

    /**
     * Get the actual object contained by handle
     *
     *
     * @return
     */
    T& get_object() const
    {
      return *t;
    }

  private:

    /// used as a unique object signature
    indexed_node* signature;
    /// type checkig information
    base_id_t type;
    /// object pointer
    T *t;

    /// Collector
    friend class indexed_node_collector<T>;
  };

  namespace nodes {

    /**
       \defgroup IndexedNode helper functions
       @{
    */

    /**
     * Create integer handle for object
     *
     * @param t The object
     *
     * @return
     */
    template <typename T> uint64_t create_handle(T* t)
    {
      indexed_node<T>* handle = new indexed_node<T>(t);
      return handle->to_index();
    }

    /**
     * Obtain object represented by index.
     *
     * @param index
     *
     * @return
     */
    template <typename T> T& get_object(const uint64_t index) throw(std::runtime_error)
    {
      indexed_node<T>* handle= indexed_node<T>::from_index(index);
      return handle->get_object();
    }

    /**
     * Obtain index_node wrapper by index
     *
     * @param index
     *
     * @return
     */
    template <typename T> indexed_node<T>& get_node(const uint64_t index) throw(std::runtime_error)
    {
      indexed_node<T>* handle= indexed_node<T>::from_index(index);
      return *handle;
    }

    /**
     * Destroy object. If deleting object, rather than leaving it to
     * garbage collection, you must delete it via the handle; do not
     * delete T* directly.
     *
     * @param index
     */
    template <typename T> void destroy_object(const uint64_t index) throw(std::runtime_error)
    {
      indexed_node<T>* handle= indexed_node<T>::from_index(index);
      delete handle;
    }

    /**
     * Clone object using index
     *
     * @param index
     *
     * @return
     */
    template <typename T> uint64_t clone_object(const uint64_t index) throw(std::runtime_error)
    {
      indexed_node<T>* hCurrent = indexed_node<T>::from_index(index);
      T* clone = new T(hCurrent->get_object());
      indexed_node<T>* hClone = new indexed_node<T>(clone);
      return hClone->to_index();
    }

//@}

  };

  template <typename T>
  class indexed_node_collector {
  public:
    static std::vector<indexed_node<T>*> objvector;

    /**
     * Dtor
     *
     */
    ~indexed_node_collector()
    {
      size_t nObjectsCleared = 0;

      typename std::vector<indexed_node<T>*>::iterator i;
      typename std::vector<indexed_node<T>*>::iterator end = objvector.end();
      for (i = objvector.begin(); i!=end; ++i) {
        // check for valid signature
        if ((*i)->signature == *i) {
          delete *i;
          nObjectsCleared++;
        }
      }
      if (nObjectsCleared) {
        fprintf(stderr,"Garbage collector: ");
        fprintf(stderr,"Cleared %zu %s item(s)\n",
                nObjectsCleared, this->objname);
      }
    }

    static void register_handle (indexed_node<T>* obj)
    {
      static indexed_node_collector singleton(obj);

      typename std::vector<indexed_node<T>*>::iterator i;
      typename std::vector<indexed_node<T>*>::iterator end = objvector.end();

      i = std::find(objvector.begin(),objvector.end(),nullptr);
      if (i!=end) {
        (*i) = obj;
      } else {
        singleton.objvector.push_back(obj);
      }
    }
  private:
    /// Name
    char objname[256];

    /**
     * Private ctor (prevent construction)
     *
     * @param obj
     */
    indexed_node_collector(indexed_node<T>* obj)
    {
      char buffer[128];
#ifdef SPS_TYPE_ID_NAME
      sprintf(buffer, "%zu", strlen(obj->type));
#else
      sprintf(buffer, "%zu", strlen(obj->type.name()));
#endif

#if defined(__GNUC__)
      // GNU use the length followed by the name
      strcpy(objname,obj->type+strlen(buffer));
#elif defined(_MSC_VER)
      // Microsoft have 6 extra characters in names
      strcpy(objname,obj->type+6);
#endif
    }
  };

/// Explicit instantiation
  template <typename T>
  std::vector<indexed_node<T>*> indexed_node_collector<T>::objvector;

  template <typename T>
  uint64_t indexed_node<T>::to_index()
  {
    // Find pointer in vector
    auto it = find(indexed_node_collector<T>::objvector.begin(),indexed_node_collector<T>::objvector.end(),this);
    if (it!=indexed_node_collector<T>::objvector.end()) {
      return uint64_t(it - indexed_node_collector<T>::objvector.begin());
    } else {
      fprintf(stderr, "Invalid index.\n");
      throw std::runtime_error("indexed_node<T>::to_index");
    }
  }
};
