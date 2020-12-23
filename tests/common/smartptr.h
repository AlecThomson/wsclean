#include <boost/test/unit_test.hpp>

namespace test {

/**
 * Helper class for using using unique pointers in tests.
 * @tparam T Object type for the unique pointer.
 */
template <class T>
class UniquePtr {
 public:
  /**
   * Constructor. It uses perfect forwarding for calling T's constructor.
   */
  template <typename... Args>
  UniquePtr(Args&&... args)
      : _unique(new T(std::forward<Args>(args)...)), _raw(_unique.get()) {}

  /**
   * Extract the unique_ptr from the helper class.
   * The test code may only call this function once.
   */
  std::unique_ptr<T> take() {
    BOOST_TEST_REQUIRE(static_cast<bool>(_unique));
    return std::move(_unique);
  }

  T& operator*() { return *_raw; }
  T* operator->() { return _raw; }
  T* get() { return _raw; }

 private:
  std::unique_ptr<T> _unique;
  T* _raw;
};

}  // namespace test
