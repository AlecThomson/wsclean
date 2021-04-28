#ifndef MSREADERS_TIMESTEPBUFFERREADER_H
#define MSREADERS_TIMESTEPBUFFERREADER_H

#include "msreader.h"

class TimestepBufferReader final : public MSReader {
 public:
  TimestepBufferReader(MSProvider* msProvider) : MSReader(msProvider){};
  virtual ~TimestepBufferReader(){};
};

#endif