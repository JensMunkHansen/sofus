#include <sps/config.h>
#include <sps/stdlib.h>
#include <sps/cstdio>

#include <list>
#include <memory>

#define _INDEXED_TYPES_CLONEABLE 1

#include <sps/indexed_types.hpp>
#include <sps/indexed_types_custom.hpp>
#include <sps/indexed_nodes.hpp>
#include <sps/aligned_allocator.hpp>
#include <iostream>
#include <vector>

class A {
public:
  A()
  {
    data = std::shared_ptr<int>(new int(2));
  }
#ifdef _INDEXED_TYPES_CLONEABLE
  A* clone() const
  {
    A* out = new A(*this);
    out->data = std::shared_ptr<int>(new int(*this->data.get()));
    return out;
  }
#endif
private:
public:
  std::shared_ptr<int> data;
};



using namespace sps;

int main(const int argc, const char* argv[])
{
  std::vector<int,aligned_allocator<int> > vec = std::vector<int, aligned_allocator<int> >(4);

  A* pA = new A();
  uint64_t i = types::create_handle<A>(pA);
  fprintf(stdout, "index of new object is: %zu\n",i);

  pA = &(types::get_object<A>(i));

  *pA->data.get() = 4;

  printf("%d\n",*pA->data.get());

  *pA->data = 5;

  printf("%d\n",*pA->data.get());

#ifdef _INDEXED_TYPES_CLONEABLE
  uint64_t j = types::clone_object<A>(i);
  fprintf(stdout, "index of new object is: %zu\n",j);
  types::destroy_object<A>(j);
#endif

  // Nodes not working
  A* pA1 = new A();
  uint64_t k = types::create_handle(pA1);
  fprintf(stdout, "index of new object is: %zu\n",k);
  // assert equal to 1

  types::destroy_object<A>(i);
  A* pA2 = new A();
  uint64_t l = types::create_handle(pA2);
  fprintf(stdout, "index of new object is: %zu\n",l);
  // assert equal to 0

  //throw sps::bad_indexed_auto_ptr();
  return EXIT_SUCCESS;
}


