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
#include <queue>
#include <thread>
#include <mutex>
#include <condition_variable>

namespace sps {
  template <typename T>
  class queue {
  public:

    T pop()
    {
      std::unique_lock<std::mutex> mlock(mutex_);
      while (queue_.empty()) {
        cond_.wait(mlock);
      }
      auto item = queue_.front();
      queue_.pop();
      return item;
    }

    void pop(T& item)
    {
      std::unique_lock<std::mutex> mlock(mutex_);
      while (queue_.empty()) {
        cond_.wait(mlock);
      }
      item = queue_.front();
      queue_.pop();
    }

    void push(const T& item)
    {
      std::unique_lock<std::mutex> mlock(mutex_);
      queue_.push(item);
      mlock.unlock();
      cond_.notify_one();
    }

    void push(T&& item)
    {
      std::unique_lock<std::mutex> mlock(mutex_);
      queue_.push(std::move(item));
      mlock.unlock();
      cond_.notify_one();
    }

  private:
    std::queue<T> queue_;
    std::mutex mutex_;
    std::condition_variable cond_;
  };
}
