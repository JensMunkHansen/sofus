#include <sps/progress.hpp>
#include <iostream>
namespace sps {

  void ProgressBarStdOut::show(float percent)
  {
    std::cout << "Progress: " << percent << " %" << std::endl;
  }

  // Why is this needed??
  void ProgressBarInterface::show(float percent)
  {
    std::cout << "Progress: " << percent << " %" << std::endl;
  }
}
