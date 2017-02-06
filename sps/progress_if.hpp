#pragma once

#include <sps/sps_export.h>

namespace sps {
  class SPS_EXPORT ProgressBarInterface {
  public:
    virtual void show(float percent) = 0;
    virtual ~ProgressBarInterface() = default;
    // Base functions
    //virtual void cancel();

    template <typename P>
    void wait(const P& period);
  };
}
