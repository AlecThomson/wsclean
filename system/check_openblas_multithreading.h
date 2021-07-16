#include <cstddef>
#include <dlfcn.h>

// Detect linkage to a multi-threaded version of OpenBLAS
// OpenBLAS multi-threading interfers with multi-threading in wsclean and OpenMP
class OpenBLASMultithreadingCheck {
 public:
  OpenBLASMultithreadingCheck() {
    // Ask the dynamic linker to lookup the openblas_set_num_threads function
    void (*openblas_set_num_threads)(int) = reinterpret_cast<void (*)(int)>(
        dlsym(RTLD_DEFAULT, "openblas_set_num_threads"));

    // If openblas_set_num_threads is present, the executable is linked to a
    // multithreaded version of OpenBLAS
    if (openblas_set_num_threads) {
      _multithreaded_openblas_linkage = openblas_set_num_threads;

      // Read the OPENBLAS_NUM_THREADS environment variable
      _openblas_num_threads = 0;
      char *openblas_num_threads_env_var = getenv("OPENBLAS_NUM_THREADS");
      if (openblas_num_threads_env_var != nullptr) {
        _openblas_num_threads = atoi(openblas_num_threads_env_var);
      }

      if (_openblas_num_threads != 1) {
        setenv("OPENBLAS_NUM_THREADS", "1", 1);

        std::clog << "WSClean has been linked to a multi-threaded version of "
                     "OpenBLAS. \n"
                     "OpenBLAS multi-threading interfers with other "
                     "multi-threaded parts "
                     "of WCSlean.\n"
                     "This has a severe impact on performance.\n"
                     "Temporary setting OPENBLAS_NUM_THREADS=1\n\n";
      }
    }
  }

  ~OpenBLASMultithreadingCheck() {
    if (_openblas_num_threads == 0) {
      unsetenv("OPENBLAS_NUM_THREADS");
    } else if (_openblas_num_threads != 1) {
      int n = 10;
      char openblas_num_threads_str[n];
      snprintf(openblas_num_threads_str, n, "%d", _openblas_num_threads);
      setenv("OPENBLAS_NUM_THREADS", openblas_num_threads_str, 1);
      std::clog << "Restoring OPENBLAS_NUM_THREADS=" << _openblas_num_threads
                << "\n";
    }
  }

 private:
  bool _multithreaded_openblas_linkage = false;
  int _openblas_num_threads = 1;
};
