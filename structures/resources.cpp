#include "resources.h"

#include <atomic>
#include <cassert>
#include <cmath>
#include <unistd.h>  // for sysconf

#include <aocommon/logger.h>

using aocommon::Logger;

const Resources Resources::GetPart(size_t part_size) const {
  assert(part_size != 0);
  const size_t n_cpus =
      std::max<size_t>(1, (n_cpus_ + part_size - 1) / part_size);
  const int64_t memory = memory_ / part_size;
  return Resources(n_cpus, memory);
}

int64_t GetAvailableMemory(double memory_fraction, double abs_memory_limit) {
  assert(memory_fraction > 0.0 && memory_fraction <= 1.0);
  // During the first run of this function, some information is reported
  // This needs to be thread safe because gridders can call this function in
  // parallel.
  static std::atomic<bool> isFirst(true);
  const bool print_output = isFirst.exchange(false);
  const long page_count = sysconf(_SC_PHYS_PAGES);
  const long page_size = sysconf(_SC_PAGE_SIZE);
  int64_t memory = (int64_t)page_count * (int64_t)page_size;
  const double memory_size_in_gb = (double)memory / (1024.0 * 1024.0 * 1024.0);
  if (memory_fraction == 1.0 && abs_memory_limit == 0.0 && print_output) {
    Logger::Info << "Detected " << std::round(memory_size_in_gb * 10.0) / 10.0
                 << " GB of system memory, usage not limited.\n";
  } else {
    double limit_in_gb = memory_size_in_gb * memory_fraction;
    if (abs_memory_limit != 0.0)
      limit_in_gb = std::min(limit_in_gb, abs_memory_limit);
    if (print_output) {
      Logger::Info << "Detected " << std::round(memory_size_in_gb * 10.0) / 10.0
                   << " GB of system memory, usage limited to "
                   << std::round(limit_in_gb * 10.0) / 10.0 << " GB (frac="
                   << std::round(memory_fraction * 1000.0) / 10.0 << "%, ";
      if (abs_memory_limit == 0.0)
        Logger::Info << "no limit)\n";
      else
        Logger::Info << "limit=" << std::round(abs_memory_limit * 10.0) / 10.0
                     << "GB)\n";
    }

    memory = int64_t((double)page_count * (double)page_size * memory_fraction);
    if (abs_memory_limit != 0.0 &&
        double(memory) > double(1024.0 * 1024.0 * 1024.0) * abs_memory_limit)
      memory =
          int64_t(double(abs_memory_limit) * double(1024.0 * 1024.0 * 1024.0));
  }
  return memory;
}
