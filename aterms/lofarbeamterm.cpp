#include "lofarbeamterm.h"

#include <aocommon/banddata.h>
#include <aocommon/matrix2x2.h>

#include "../system.h"

#include "../lofar/lofarbeamkeywords.h"

#include <aocommon/imagecoordinates.h>

#include "../wsclean/logger.h"

#include <EveryBeam/options.h>
#include <EveryBeam/element_response.h>

#include <EveryBeam/coords/ITRFDirection.h>
#include <EveryBeam/coords/ITRFConverter.h>
#include <EveryBeam/gridded_response/griddedresponse.h>

#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/tables/Tables/TableRecord.h>

#include <aocommon/lane.h>

#include <algorithm>

using namespace everybeam;
using namespace aocommon;

LofarBeamTerm::LofarBeamTerm(casacore::MeasurementSet& ms, const CoordinateSystem& coordinateSystem, const std::string& dataColumnName) :
	_width(coordinateSystem.width),
	_height(coordinateSystem.height),
	_dl(coordinateSystem.dl), _dm(coordinateSystem.dm),
	_phaseCentreDL(coordinateSystem.phaseCentreDL),
	_phaseCentreDM(coordinateSystem.phaseCentreDM),
	_useDifferentialBeam(false),
	_useChannelFrequency(true)
{	
	// 

	// Default options 
	Options options;

	// Response model
	ElementResponseModel response_model = ElementResponseModel::kHamaker;
	
	// Load telescope
	telescope_ = Load(ms, response_model, options);

	casacore::MSAntenna aTable(ms.antenna());

	casacore::MSField fieldTable(ms.field());
	if(fieldTable.nrow() != 1)
		throw std::runtime_error("Set has multiple fields");

	casacore::MPosition::ScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	casacore::MPosition arrayPos = antPosColumn(0);
	
	// Compute Ra/Dec pointing direction
	casacore::MEpoch::ScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	casacore::MDirection::ScalarColumn phaseDirColumn(fieldTable, fieldTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
	casacore::MDirection phaseDir = phaseDirColumn(0);
	casacore::MEpoch curtime = timeColumn(0);
	casacore::MeasFrame frame(arrayPos, curtime);
	casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	casacore::MDirection j2000 = casacore::MDirection::Convert(phaseDir, j2000Ref)();
	casacore::Vector<casacore::Double> j2000Val = j2000.getValue().get();
	_phaseCentreRA = j2000Val[0];
	_phaseCentreDec = j2000Val[1];

	_coordinate_system = {.width=_width, .height=_height, .ra=_phaseCentreRA, .dec=_phaseCentreDec, .dl=_dl, .dm=_dm, .phase_centre_dl=_phaseCentreDL, .phase_centre_dm=_phaseCentreDM};

	// casacore::MSAntenna aTable(ms.antenna());

	// size_t nStations = aTable.nrow();
	// size_t nCPUs = System::ProcessorCount();
	// _nThreads = std::min(nCPUs, nStations);
	// _threads.resize(_nThreads);

	// casacore::MPosition::ScalarColumn antPosColumn(aTable, aTable.columnName(casacore::MSAntennaEnums::POSITION));
	// _arrayPos = antPosColumn(0);
	// _stations.resize(aTable.nrow());
	// readStations(ms, _stations.begin());
	
	// BandData band(ms.spectralWindow());
	// _subbandFrequency = band.CentreFrequency();
	
	// casacore::MSField fieldTable(ms.field());
	// if(fieldTable.nrow() != 1)
	// 	throw std::runtime_error("Set has multiple fields");
	
	// casacore::MEpoch::ScalarColumn timeColumn(ms, ms.columnName(casacore::MSMainEnums::TIME));
	// casacore::MDirection::ScalarColumn phaseDirColumn(fieldTable, fieldTable.columnName(casacore::MSFieldEnums::PHASE_DIR));
	// casacore::MDirection phaseDir = phaseDirColumn(0);
	// casacore::MEpoch curtime = timeColumn(0);
	// casacore::MeasFrame frame(_arrayPos, curtime);
	// casacore::MDirection::Ref j2000Ref(casacore::MDirection::J2000, frame);
	// casacore::MDirection j2000 = casacore::MDirection::Convert(phaseDir, j2000Ref)();
	// casacore::Vector<casacore::Double> j2000Val = j2000.getValue().get();
	// _phaseCentreRA = j2000Val[0];
	// _phaseCentreDec = j2000Val[1];

	// casacore::ScalarMeasColumn<casacore::MDirection> delayDirColumn(fieldTable, casacore::MSField::columnName(casacore::MSFieldEnums::DELAY_DIR));
	// _delayDir = delayDirColumn(0);
	
	// if(fieldTable.tableDesc().isColumn("LOFAR_TILE_BEAM_DIR")) {
	// 	casacore::ArrayMeasColumn<casacore::MDirection> tileBeamDirColumn(fieldTable, "LOFAR_TILE_BEAM_DIR");
	// 	_tileBeamDir = *(tileBeamDirColumn(0).data());
	// } else {
	// 	throw std::runtime_error("LOFAR_TILE_BEAM_DIR column not found");
	// }
	
	// _useDifferentialBeam = LOFARBeamKeywords::GetPreappliedBeamDirection(ms, dataColumnName, _useDifferentialBeam, _preappliedBeamDir);
	// Logger::Debug << "Tile direction: " << RaDecCoord::RaDecToString(_tileBeamDir.getAngle().getValue()[0], _tileBeamDir.getAngle().getValue()[1]) << '\n';
	// Logger::Debug << "Delay direction: " << RaDecCoord::RaDecToString(_delayDir.getAngle().getValue()[0], _delayDir.getAngle().getValue()[1]) << '\n';
}

bool LofarBeamTerm::calculateBeam(std::complex<float>* buffer, double time, double frequency, size_t)
{
	// Get the gridded response
	std::unique_ptr<gridded_response::GriddedResponse> grid_response = telescope_->GetGriddedResponse(_coordinate_system);

	// // Compute the beam (AllStations)
	bool grid_response_pass = grid_response->CalculateAllStations(buffer, time,
                                                                  frequency);
	return grid_response_pass;

	// 
	// aocommon::Lane<size_t> lane(_nThreads);
	// _lane = &lane;

	// LOFAR::StationResponse::ITRFConverter itrfConverter(time);
	// setITRFVector(itrfConverter.toDirection(_delayDir), _station0);
	// setITRFVector(itrfConverter.toDirection(_tileBeamDir), _tile0);

	// const casacore::Unit radUnit("rad");

	// casacore::MDirection lDir(casacore::MVDirection(
	// 	casacore::Quantity(_phaseCentreRA + M_PI/2, radUnit),
	// 	casacore::Quantity(0, radUnit)),
	// 	casacore::MDirection::J2000);
	// setITRFVector(itrfConverter.toDirection(lDir), _l_vector_itrf);

	// casacore::MDirection mDir(casacore::MVDirection(
	// 	casacore::Quantity(_phaseCentreRA, radUnit),
	// 	casacore::Quantity(_phaseCentreDec + M_PI/2, radUnit)),
	// 	casacore::MDirection::J2000);
	// setITRFVector(itrfConverter.toDirection(mDir), _m_vector_itrf);

	// casacore::MDirection nDir(casacore::MVDirection(
	// 	casacore::Quantity(_phaseCentreRA, radUnit),
	// 	casacore::Quantity(_phaseCentreDec, radUnit)),
	// 	casacore::MDirection::J2000);
	// setITRFVector(itrfConverter.toDirection(nDir), _n_vector_itrf);

	// vector3r_t diffBeamCentre;
	// setITRFVector(itrfConverter.toDirection(_preappliedBeamDir), diffBeamCentre);
	// _inverseCentralGain.resize(_stations.size());
	// for(size_t a=0; a!=_stations.size(); ++a)
	// {
	// 	double sbFrequency = _useChannelFrequency ? frequency : _subbandFrequency;
	// 	matrix22c_t gainMatrix = _stations[a]->response(time, frequency, diffBeamCentre, sbFrequency, _station0, _tile0);
	// 	if(_useDifferentialBeam)
	// 	{
	// 		_inverseCentralGain[a][0] = gainMatrix[0][0];
	// 		_inverseCentralGain[a][1] = gainMatrix[0][1];
	// 		_inverseCentralGain[a][2] = gainMatrix[1][0];
	// 		_inverseCentralGain[a][3] = gainMatrix[1][1];
	// 		if(!_inverseCentralGain[a].Invert())
	// 		{
	// 			_inverseCentralGain[a] = MC2x2F::Zero();
	// 		}
	// 	}
	// }

	// for(size_t i=0; i!=_nThreads; ++i)
	// {
	// 	_threads[i] = std::thread(&LofarBeamTerm::calcThread, this, buffer, time, frequency);
	// }

	// for(size_t y=0; y!=_height; ++y)
	// {
	// 	for(size_t antennaIndex=0; antennaIndex!=_stations.size(); ++antennaIndex)
	// 	{
	// 		size_t job_id = y*_stations.size() + antennaIndex;
	// 		lane.write(job_id);
	// 	}
	// }

	// lane.write_end();
	// for(size_t i=0; i!=_nThreads; ++i)
	// 	_threads[i].join();
	
	// saveATermsIfNecessary(buffer, _stations.size(), _width, _height);
}