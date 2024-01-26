#ifndef DISTRIBUTED_WORKER_H_
#define DISTRIBUTED_WORKER_H_

#include "../main/settings.h"

#include "../scheduling/mpiworkerscheduler.h"

class Worker {
 public:
  Worker(const Settings& settings) : scheduler_{settings} {}

  void Run();

 private:
  void grid(size_t body_size);

  MpiWorkerScheduler scheduler_;
};

#endif
