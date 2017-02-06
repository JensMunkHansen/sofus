#include <sps/memory>

#include <gtest/gtest.h>

class A {
public:
  A()
  {
    m_data = sps::deleted_aligned_array_create<float>(20);
  }
private:
  sps::deleted_aligned_array<float> m_data;
};

TEST(memory_test, aligned_unique_array)
{
  auto a = sps::deleted_aligned_array_create<float>(0);

  auto b = sps::deleted_aligned_array_create<float>(35);

  auto c = sps::deleted_aligned_multi_array_create<double,2>(1,1);

  auto d = sps::deleted_aligned_multi_array_create<double,2>(0,0);

  auto f = sps::deleted_aligned_multi_array_create<double,2>(1,0);

  auto g = sps::deleted_aligned_multi_array_create<double,2>(0,1);

  ASSERT_TRUE(true);
}

TEST(memory_test, unique_member)
{
  A a;
}

int main(int argc, char* argv[])
{
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}

// GNU extension
//
