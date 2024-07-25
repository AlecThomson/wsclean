#ifndef DISTRIBUTED_WORKER_H_
#define DISTRIBUTED_WORKER_H_

#include "../main/settings.h"

#include "../scheduling/mpiworkerscheduler.h"

namespace wsclean {

class Worker {
 public:
  Worker(const Settings& settings) : scheduler_{settings} {}

  void Run();

 private:
  MpiWorkerScheduler scheduler_;
};

}  // namespace wsclean

#endif
