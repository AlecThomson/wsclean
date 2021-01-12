#include "griddingresult.h"

#include "../io/serialostream.h"
#include "../io/serialistream.h"

void GriddingResult::Serialize(SerialOStream& stream) const {
  imageRealResult.Serialize(stream);
  imageImaginaryResult.Serialize(stream);
  stream.Double(beamSize)
      .Double(imageWeight)
      .Double(normalizationFactor)
      .UInt64(actualWGridSize)
      .UInt64(griddedVisibilityCount)
      .Double(effectiveGriddedVisibilityCount)
      .Double(visibilityWeightSum)
      .UInt64(actualInversionWidth)
      .UInt64(actualInversionHeight)
      .Ptr(cache);
}

void GriddingResult::Unserialize(SerialIStream& stream) {
  imageRealResult.Unserialize(stream);
  imageImaginaryResult.Unserialize(stream);
  stream.Double(beamSize)
      .Double(imageWeight)
      .Double(normalizationFactor)
      .UInt64(actualWGridSize)
      .UInt64(griddedVisibilityCount)
      .Double(effectiveGriddedVisibilityCount)
      .Double(visibilityWeightSum)
      .UInt64(actualInversionWidth)
      .UInt64(actualInversionHeight)
      .Ptr(cache);
}
