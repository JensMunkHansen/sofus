#pragma once

namespace sps {
  class ProgressBarInterface {
  public:
    virtual void show(float percent);
    virtual ~ProgressBarInterface() = default;
  };

  class ProgressBarStdOut : public ProgressBarInterface {
  public:
    virtual void show(float percent);
    virtual ~ProgressBarStdOut() = default;
  };
}
