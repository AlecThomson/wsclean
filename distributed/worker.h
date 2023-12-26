#ifndef DISTRIBUTED_WORKER_H_
#define DISTRIBUTED_WORKER_H_

#include "../main/settings.h"
#include "../scheduling/griddingtaskmanager.h"

class Worker {
 public:
  Worker(const Settings& settings)
      : scheduler_(GriddingTaskManager::Make(settings)) {}

  void Run();

 private:
  void grid(size_t bodySize);

  std::unique_ptr<GriddingTaskManager> scheduler_;
};

#endif
