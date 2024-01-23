#ifndef DISTRIBUTED_WORKER_H_
#define DISTRIBUTED_WORKER_H_

#include "../main/settings.h"

class Worker {
 public:
  Worker(const Settings& settings);

  void Run();

 private:
  void grid(size_t body_size);

  const Settings settings_;
  int rank_;
};

#endif
