#include "timestepbuffer.h"
#include "../msreaders/timestepbufferreader.h"

std::unique_ptr<MSReader> TimestepBuffer::GetReader() {
    std::cout << "Should return a reader?!"<<std::endl;
    std::unique_ptr<MSReader> reader(new TimestepBufferReader(this));
    std::cout << "Created a reader, only return it!"<<std::endl;
    return reader;
}