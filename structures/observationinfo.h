#ifndef OBSERVATION_INFO_H
#define OBSERVATION_INFO_H

#include <string>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

struct ObservationInfo {
  double phaseCentreRA = 0.0, phaseCentreDec = 0.0;
  double startTime = 0.0;
  bool hasDenormalPhaseCentre = false;
  double phaseCentreDL = 0.0, phaseCentreDM = 0.0;
  std::string telescopeName;
  std::string observer;
  std::string fieldName;

  void Serialize(class SerialOStream& stream) const;
  void Unserialize(class SerialIStream& stream);

  static ObservationInfo ReadObservationInfo(casacore::MeasurementSet& ms,
                                             size_t fieldId);
};

#endif
