#pragma once

#include <sps/sps_export.h>
#include <sps/progress_if.hpp>

#include <mutex>
#include <condition_variable>

// TODO: Export only interface
namespace sps {
  class ProgressBar : public ProgressBarInterface {
  public:
    virtual void show(float percent) = 0;
    virtual ~ProgressBar() = default;

    // Base functions
    //void cancel();

    template <typename P>
    void wait(const P& period);
  private:
    bool bStop;
    std::mutex mtx;
    std::condition_variable cnd;
  };

  class SPS_EXPORT ProgressBarStdOut : public ProgressBarInterface {
      // Add stuff like above
  public:
    virtual void show(float percent);
    virtual ~ProgressBarStdOut() = default;
  };
}
