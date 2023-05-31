#ifndef SLAVE_H
#define SLAVE_H

#include "../main/settings.h"

class Worker {
 public:
  Worker(const Settings& settings) : _settings(settings) {}

  void Run();

 private:
  void grid(size_t bodySize);

  const Settings _settings;
};

#endif
