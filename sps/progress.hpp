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
