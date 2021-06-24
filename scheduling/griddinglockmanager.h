#ifndef GRIDDING_LOCK_MANAGER_H
#define GRIDDING_LOCK_MANAGER_H

#include <mutex>

/**
 * @brief Abstract interface class providing access to the
 * locking mechanism of the GriddingTaskManager or childs thereof.
 */
class GriddingLockManager {
 protected:
  class WriterLock {
   public:
    virtual void lock() = 0;
    virtual void unlock() = 0;
  };

 public:
  virtual ~GriddingLockManager(){};

  using WriterGroupLockGuard = std::lock_guard<WriterLock>;

  virtual WriterGroupLockGuard LockWriterGroup(size_t writerGroupIndex) = 0;
};
#endif