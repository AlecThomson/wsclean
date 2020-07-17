#ifndef LOFAR_BEAM_TERM_H
#define LOFAR_BEAM_TERM_H

#include <thread>

#include <casacore/ms/MeasurementSets/MeasurementSet.h>

#include <complex>

#include "atermstub.h"
#include "atermbeam.h"

#ifdef HAVE_EVERYBEAM

#include <EveryBeam/load.h>
#include <EveryBeam/coords/coord_utils.h>

class LofarBeamTerm : public ATermBeam
{
public:
	LofarBeamTerm(casacore::MeasurementSet& ms, const CoordinateSystem& coordinateSystem, const std::string& dataColumnName);
	
	void SetUseDifferentialBeam(bool useDiffBeam)
	{
		_useDifferentialBeam = useDiffBeam;
	}
	
	void SetUseChannelFrequency(bool useChannelFrequency)
	{
		_useChannelFrequency = useChannelFrequency;
	}
	
private:
	bool calculateBeam(std::complex<float>* buffer, double time, double frequency, size_t fieldId) final override;

	// std::vector<LOFAR::StationResponse::Station::Ptr> _stations;
	std::unique_ptr<everybeam::telescope::Telescope> telescope_;
	size_t _width, _height;
	double _subbandFrequency, _phaseCentreRA, _phaseCentreDec, _dl, _dm, _phaseCentreDL, _phaseCentreDM;
	// casacore::MDirection _delayDir, _preappliedBeamDir, _tileBeamDir;
	// casacore::MPosition _arrayPos;
	casacore::MDirection _preappliedBeamDir;
	bool _useDifferentialBeam, _useChannelFrequency;

	everybeam::coords::CoordinateSystem _coordinate_system;

	// LOFAR::StationResponse::vector3r_t _l_vector_itrf;
	// LOFAR::StationResponse::vector3r_t _m_vector_itrf;
	// LOFAR::StationResponse::vector3r_t _n_vector_itrf;
	// std::vector<aocommon::MC2x2F> _inverseCentralGain;
	// LOFAR::StationResponse::vector3r_t _station0, _tile0;
	
	// aocommon::Lane<size_t> *_lane;
	// size_t _nThreads;
	// std::vector<std::thread> _threads;
};

#else
using LofarBeamTerm = ATermStub;
#endif // HAVE_EVERYBEAM

#endif 
