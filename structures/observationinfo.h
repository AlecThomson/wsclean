#ifndef OBSERVATION_INFO_H
#define OBSERVATION_INFO_H

#include <aocommon/io/serialstreamfwd.h>

#include <string>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

struct ObservationInfo {
  double phaseCentreRA = 0.0;
  double phaseCentreDec = 0.0;
  bool hasShiftedPhaseCentre = false;
  double shiftL = 0.0;
  double shiftM = 0.0;
  double facetShiftL = 0.0;
  double facetShiftM = 0.0;
  std::string telescopeName;
  std::string observer;
  std::string fieldName;

  void Serialize(aocommon::SerialOStream& stream) const;
  void Unserialize(aocommon::SerialIStream& stream);
};

/// Generates observation info from the measurement set tables.
struct ObservationInfo ReadObservationInfo(casacore::MeasurementSet& ms,
                                           size_t fieldId);

#endif
