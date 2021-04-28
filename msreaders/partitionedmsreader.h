#ifndef MSREADERS_PARTITIONEDMSREADER_H
#define MSREADERS_PARTITIONEDMSREADER_H

#include "msreader.h"

class PartitionedMSReader final : public MSReader {
 public:
  PartitionedMSReader(MSProvider* msProvider) : MSReader(msProvider){};
  virtual ~PartitionedMSReader(){};
};

#endif