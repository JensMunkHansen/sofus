#include <sps/cenv.h>
#include <sps/progress.hpp>
#include <iostream>

namespace sps {

  void ProgressBar::show(float percent)
  {
    // Is this needed???
    SPS_UNREFERENCED_PARAMETERS(percent);
  }

#if 0
  void ProgressBar::cancel()
  {
    std::unique_lock<std::mutex> lock(mtx);
    // Order is insignificant so no need for std::atomic<bool>
    bStop = true;
    cnd.notify_all();
  }
#endif

  // std::chrono::duration::zero()
  template <typename P>
  void ProgressBar::wait(const P& period)
  {
    std::unique_lock<std::mutex> lock(mtx);
    if (bStop || cnd.wait_for(lock, period) == std::cv_status::no_timeout) {
      bStop = false;
      //throw 1;//cancelled_error();
    }
  }

  void ProgressBarStdOut::show(float percent)
  {
    std::cout << "Progress: " << percent << " %" << std::endl;
  }
}
