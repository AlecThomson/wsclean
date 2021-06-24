#ifndef GRIDDING_LOCK_MANAGER_H
#define GRIDDING_LOCK_MANAGER_H

/**
 * @brief Abstract interface class providing access to the
 * locking mechanism of the GriddingTaskManager or childs thereof.
 */
class GriddingLockManager {
 protected:
  class WriterLockBase {
   public:
    virtual void lock() = 0;
    virtual void unlock() = 0;
  };

 public:
  virtual ~GriddingLockManager(){};

  // Becomes obsolete --> std::lock_guard<WriterLock>
  class WriterGroupLockGuard {
   public:
    WriterGroupLockGuard(WriterLockBase& lock, size_t& counter)
        : _lock(lock), _writerGroupCounter(counter) {
      _lock.lock();
    }
    ~WriterGroupLockGuard() {
      ++_writerGroupCounter;
      _lock.unlock();
    }
    size_t GetCounter() const { return _writerGroupCounter; }

   private:
    WriterLockBase& _lock;
    size_t& _writerGroupCounter;
  };

  virtual WriterGroupLockGuard LockWriterGroup(
      size_t writerGroupIndex) const = 0;
};
#endif