#include "wsclean.h"

#include "areaset.h"
#include "fitswriter.h"
#include "imageweights.h"
#include "inversionalgorithm.h"
#include "modelrenderer.h"
#include "msselection.h"
#include "wsinversion.h"

#include "model/model.h"

#include "cleanalgorithms/cleanalgorithm.h"
#include "cleanalgorithms/joinedclean.h"
#include "cleanalgorithms/simpleclean.h"
#include "cleanalgorithms/multiscaleclean.h"

#include "msproviders/contiguousms.h"

#include "model/areaparser.h"

#include "uvector.h"
#include "gaussianfitter.h"
#include "angle.h"
#include "dftpredictionalgorithm.h"

#include "lofar/lmspredicter.h"
#include "progressbar.h"

#include <iostream>
#include <memory>

std::string commandLine;

WSClean::WSClean() :
	_imgWidth(2048), _imgHeight(2048), _channelsOut(1),
	_pixelScaleX(0.01 * M_PI / 180.0), _pixelScaleY(0.01 * M_PI / 180.0),
	_threshold(0.0), _gain(0.1), _mGain(1.0), _cleanBorderRatio(0.05),
	_manualBeamMajorSize(0.0), _manualBeamMinorSize(0.0),
	_manualBeamPA(0.0), _fittedBeam(false), _circularBeam(false),
	_memFraction(1.0), _absMemLimit(0.0),
	_minUVInLambda(0.0), _maxUVInLambda(0.0), _wLimit(0.0),
	_multiscaleThresholdBias(0.7), _multiscaleScaleBias(0.6),
	_nWLayers(0), _nIter(0), _antialiasingKernelSize(7), _overSamplingFactor(63),
	_threadCount(sysconf(_SC_NPROCESSORS_ONLN)),
	_globalSelection(),
	_columnName(), _cleanAreasFilename(),
	_polarizations(),
	_weightMode(WeightMode::UniformWeighted),
	_prefixName("wsclean"),
	_allowNegative(true), _smallPSF(false), _smallInversion(true), _stopOnNegative(false), _makePSF(false), _isGriddingImageSaved(false),
	_dftPrediction(false), _dftWithBeam(false),
	_temporaryDirectory(),
	_forceReorder(false), _forceNoReorder(false), _joinedPolarizationCleaning(false), _joinedFrequencyCleaning(false),
	_mfsWeighting(false), _multiscale(false),
	_gridMode(LayeredImager::KaiserBessel),
	_filenames(),
	_commandLine(),
	_inversionWatch(false), _predictingWatch(false), _cleaningWatch(false),
	_isFirstInversion(true), _doReorder(false),
	_majorIterationNr(0)
{
	_polarizations.insert(Polarization::StokesI);
}

WSClean::~WSClean()
{
	freeCleanAlgorithms();
}

void WSClean::initFitsWriter(FitsWriter& writer)
{
	double
		ra = _inversionAlgorithm->PhaseCentreRA(),
		dec = _inversionAlgorithm->PhaseCentreDec(),
		pixelScaleX = _inversionAlgorithm->PixelSizeX(),
		pixelScaleY = _inversionAlgorithm->PixelSizeY(),
		freqHigh = _inversionAlgorithm->HighestFrequencyChannel(),
		freqLow = _inversionAlgorithm->LowestFrequencyChannel(),
		freqCentre = (freqHigh + freqLow) * 0.5,
		bandwidth = _inversionAlgorithm->BandEnd() - _inversionAlgorithm->BandStart(),
		beamSize = _inversionAlgorithm->BeamSize(),
		dateObs = _inversionAlgorithm->StartTime();
		
	writer.SetImageDimensions(_inversionAlgorithm->ImageWidth(), _inversionAlgorithm->ImageHeight(), ra, dec, pixelScaleX, pixelScaleY);
	writer.SetFrequency(freqCentre, bandwidth);
	writer.SetDate(dateObs);
	writer.SetPolarization(_inversionAlgorithm->Polarization());
	writer.SetOrigin("WSClean", "W-stacking imager written by Andre Offringa");
	writer.AddHistory(commandLine);
	if(_manualBeamMajorSize != 0.0) {
		writer.SetBeamInfo(_manualBeamMajorSize, _manualBeamMinorSize, _manualBeamPA);
	}
	else {
		writer.SetBeamInfo(beamSize, beamSize, 0.0);
	}
	if(_inversionAlgorithm->HasDenormalPhaseCentre())
		writer.SetPhaseCentreShift(_inversionAlgorithm->PhaseCentreDL(), _inversionAlgorithm->PhaseCentreDM());
	
	writer.SetExtraKeyword("WSCIMGWG", _inversionAlgorithm->ImageWeight());
	writer.SetExtraKeyword("WSCNWLAY", _inversionAlgorithm->WGridSize());
	writer.SetExtraKeyword("WSCDATAC", _inversionAlgorithm->DataColumnName());
	writer.SetExtraKeyword("WSCWEIGH", _inversionAlgorithm->Weighting().ToString());
	writer.SetExtraKeyword("WSCGKRNL", _inversionAlgorithm->AntialiasingKernelSize());
	if(_inversionAlgorithm->Selection().HasChannelRange())
	{
		writer.SetExtraKeyword("WSCCHANS", _inversionAlgorithm->Selection().ChannelRangeStart());
		writer.SetExtraKeyword("WSCCHANE", _inversionAlgorithm->Selection().ChannelRangeEnd());
	}
	if(_inversionAlgorithm->Selection().HasInterval())
	{
		writer.SetExtraKeyword("WSCTIMES", _inversionAlgorithm->Selection().IntervalStart());
		writer.SetExtraKeyword("WSCTIMEE", _inversionAlgorithm->Selection().IntervalEnd());
	}
	writer.SetExtraKeyword("WSCFIELD", _inversionAlgorithm->Selection().FieldId());
}

void WSClean::copyWSCleanKeywords(FitsReader& reader, FitsWriter& writer)
{
	const size_t
		N_STRKEYWORDS=2, N_DBLKEYWORDS=17;
	const char* strKeywords[N_STRKEYWORDS] =
		{ "WSCDATAC", "WSCWEIGH" };
	const char* dblKeywords[N_DBLKEYWORDS] =
		{ "WSCIMGWG", "WSCNWLAY", "WSCGKRNL", "WSCCHANS", "WSCCHANE", "WSCTIMES", "WSCTIMEE", "WSCFIELD",
			"WSCNITER", "WSCTHRES", "WSCGAIN", "WSCMGAIN", "WSCNEGCM", "WSCNEGST", "WSCSMPSF",
			"WSCMINOR", "WSCMAJOR"
		};
	for(size_t i=0; i!=N_STRKEYWORDS; ++i)
		writer.CopyStringKeywordIfExists(reader, strKeywords[i]);
	for(size_t i=0; i!=N_DBLKEYWORDS; ++i)
		writer.CopyDoubleKeywordIfExists(reader, dblKeywords[i]);
}

void WSClean::setCleanParameters(FitsWriter& writer, const CleanAlgorithm& clean)
{
	writer.SetExtraKeyword("WSCNITER", clean.MaxNIter());
	writer.SetExtraKeyword("WSCTHRES", clean.Threshold());
	writer.SetExtraKeyword("WSCGAIN", clean.SubtractionGain());
	writer.SetExtraKeyword("WSCMGAIN", clean.StopGain());
	writer.SetExtraKeyword("WSCNEGCM", clean.AllowNegativeComponents());
	writer.SetExtraKeyword("WSCNEGST", clean.StopOnNegativeComponents());
	writer.SetExtraKeyword("WSCSMPSF", clean.ResizePSF());
}

void WSClean::updateCleanParameters(FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr)
{
	writer.SetExtraKeyword("WSCMINOR", minorIterationNr);
	writer.SetExtraKeyword("WSCMAJOR", majorIterationNr);
}

void WSClean::imagePSF(size_t currentChannelIndex, size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Constructing PSF ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoImagePSF(true);
	_inversionAlgorithm->SetVerbose(_isFirstInversion);
	_inversionAlgorithm->Invert();
		
	CleanAlgorithm::RemoveNaNsInPSF(_inversionAlgorithm->ImageRealResult(), _imgWidth, _imgHeight);
	initFitsWriter(_fitsWriter);
	_psfImages.SetFitsWriter(_fitsWriter);
	_psfImages.Store(_inversionAlgorithm->ImageRealResult(), *_polarizations.begin(), joinedChannelIndex, false);
	_inversionWatch.Pause();
	
	_isFirstInversion = false;
	if(_fittedBeam)
	{
		double bMaj, bMin, bPA;
		GaussianFitter beamFitter;
		std::cout << "Fitting beam... " << std::flush;
		beamFitter.Fit2DGaussianCentred(
			_inversionAlgorithm->ImageRealResult(),
			_imgWidth, _imgHeight,
			_inversionAlgorithm->BeamSize()*2.0/(_pixelScaleX+_pixelScaleY),
			bMaj, bMin, bPA);
		if(bMaj < 1.0) bMaj = 1.0;
		if(bMin < 1.0) bMin = 1.0;
		bMaj = bMaj*0.5*(_pixelScaleX+_pixelScaleY);
		bMin = bMin*0.5*(_pixelScaleX+_pixelScaleY);
		std::cout << "major=" << Angle::ToNiceString(bMaj) << ", minor=" <<
		Angle::ToNiceString(bMin) << ", PA=" << Angle::ToNiceString(bPA) << ", theoretical=" <<
		Angle::ToNiceString(_inversionAlgorithm->BeamSize())<< ".\n";
		
		_manualBeamMajorSize = bMaj;
		if(_circularBeam)
		{
			_manualBeamMinorSize = _manualBeamMajorSize;
			_manualBeamPA = 0.0;
		}
		else {
			_manualBeamMinorSize = bMin;
			_manualBeamPA = bPA;
		}
		initFitsWriter(_fitsWriter);
	}
	else {
		std::cout << "Beam size is " << _inversionAlgorithm->BeamSize()*(180.0*60.0/M_PI) << " arcmin.\n";
	}
		
	std::cout << "Writing psf image... " << std::flush;
	const std::string name(getPSFPrefix(currentChannelIndex) + "-psf.fits");
	_fitsWriter.Write(name, _inversionAlgorithm->ImageRealResult());
	std::cout << "DONE\n";
}

void WSClean::imageGridding()
{
	std::cout << "Writing gridding correction image... " << std::flush;
	double* gridding = _imageAllocator.Allocate(_imgWidth * _imgHeight);
	_inversionAlgorithm->GetGriddingCorrectionImage(&gridding[0]);
	FitsWriter fitsWriter;
	initFitsWriter(fitsWriter);
	fitsWriter.SetImageDimensions(_inversionAlgorithm->ActualInversionWidth(), _inversionAlgorithm->ActualInversionHeight());
	fitsWriter.Write(_prefixName + "-gridding.fits", &gridding[0]);
	_imageAllocator.Free(gridding);
	std::cout << "DONE\n";
}

void WSClean::imageMainFirst(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Constructing image ==\n";
	_inversionWatch.Start();
	if(_nWLayers != 0)
		_inversionAlgorithm->SetWGridSize(_nWLayers);
	else
		_inversionAlgorithm->SetNoWGridSize();
	_inversionAlgorithm->SetDoImagePSF(false);
	_inversionAlgorithm->SetVerbose(_isFirstInversion);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
	_inversionAlgorithm->SetVerbose(false);
	
	storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
		storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, true, _inversionAlgorithm->ImageImaginaryResult());
}

void WSClean::imageMainNonFirst(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Constructing image ==\n";
	_inversionWatch.Start();
	_inversionAlgorithm->SetDoSubtractModel(true);
	_inversionAlgorithm->Invert();
	_inversionWatch.Pause();
	
	storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, false, _inversionAlgorithm->ImageRealResult());
	if(Polarization::IsComplex(polarization))
		storeAndCombineXYandYX(_residualImages, polarization, joinedChannelIndex, true, _inversionAlgorithm->ImageImaginaryResult());
}

void WSClean::storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image)
{
	if(polarization == Polarization::YX && _polarizations.count(Polarization::XY)!=0)
	{
		std::cout << "Adding XY and YX together...\n";
		double
			*xyImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
		dest.Load(xyImage, Polarization::XY, joinedChannelIndex, isImaginary);
		size_t count = _imgWidth*_imgHeight;
		if(isImaginary)
		{
			for(size_t i=0; i!=count; ++i)
				xyImage[i] = (xyImage[i]-image[i])*0.5;
		}
		else {
			for(size_t i=0; i!=count; ++i)
				xyImage[i] = (xyImage[i]+image[i])*0.5;
		}
		dest.Store(xyImage, Polarization::XY, joinedChannelIndex, isImaginary);
		_imageAllocator.Free(xyImage);
	}
	else {
		dest.Store(image, polarization, joinedChannelIndex, isImaginary);
	}
}

void WSClean::predict(PolarizationEnum polarization, size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Converting model image to visibilities ==\n";
	const size_t size = _imgWidth*_imgHeight;
	double
		*modelImageReal = _imageAllocator.Allocate(size),
		*modelImageImaginary = 0;
		
	if(polarization == Polarization::YX)
	{
		_modelImages.Load(modelImageReal, Polarization::XY, joinedChannelIndex, false);
		modelImageImaginary = _imageAllocator.Allocate(size);
		_modelImages.Load(modelImageImaginary, Polarization::XY, joinedChannelIndex, true);
		for(size_t i=0; i!=size; ++i)
			modelImageImaginary[i] = -modelImageImaginary[i];
	}
	else {
		_modelImages.Load(modelImageReal, polarization, joinedChannelIndex, false);
		if(Polarization::IsComplex(polarization))
		{
			modelImageImaginary = _imageAllocator.Allocate(size);
			_modelImages.Load(modelImageImaginary, polarization, joinedChannelIndex, true);
		}
	}
	
	_predictingWatch.Start();
	_inversionAlgorithm->SetAddToModel(false);
	if(Polarization::IsComplex(polarization))
		_inversionAlgorithm->Predict(modelImageReal, modelImageImaginary);
	else
		_inversionAlgorithm->Predict(modelImageReal);
	_predictingWatch.Pause();
	_imageAllocator.Free(modelImageReal);
	_imageAllocator.Free(modelImageImaginary);
}

void WSClean::dftPredict(size_t joinedChannelIndex)
{
	std::cout << std::flush << " == Predicting visibilities ==\n";
	const size_t size = _imgWidth*_imgHeight;
	double
		*modelImageReal = _imageAllocator.Allocate(size),
		*modelImageImaginary = 0;
		
	std::unique_ptr<DFTPredictionImage> image(new DFTPredictionImage(_imgWidth, _imgHeight, _imageAllocator));
		
	for(PolarizationEnum curPol : _polarizations)
	{
		if(curPol == Polarization::YX)
		{
			_modelImages.Load(modelImageReal, Polarization::XY, joinedChannelIndex, false);
			modelImageImaginary = _imageAllocator.Allocate(size);
			_modelImages.Load(modelImageImaginary, Polarization::XY, joinedChannelIndex, true);
			for(size_t i=0; i!=size; ++i)
				modelImageImaginary[i] = -modelImageImaginary[i];
			image->Add(curPol, modelImageReal, modelImageImaginary);
			_imageAllocator.Free(modelImageReal);
		}
		else {
			_modelImages.Load(modelImageReal, curPol, joinedChannelIndex, false);
			if(Polarization::IsComplex(curPol))
			{
				modelImageImaginary = _imageAllocator.Allocate(size);
				_modelImages.Load(modelImageImaginary, curPol, joinedChannelIndex, true);
				image->Add(curPol, modelImageReal, modelImageImaginary);
				_imageAllocator.Free(modelImageReal);
			}
			else {
				image->Add(curPol, modelImageReal);
			}
		}
	}
	_imageAllocator.Free(modelImageReal);
	
	casa::MeasurementSet firstMS(_filenames.front());
	BandData firstBand(firstMS.spectralWindow());
	DFTPredictionInput input;
	image->FindComponents(input, _inversionAlgorithm->PhaseCentreRA(), _inversionAlgorithm->PhaseCentreDec(), _pixelScaleX, _pixelScaleY, _inversionAlgorithm->PhaseCentreDL(), _inversionAlgorithm->PhaseCentreDM(), firstBand.ChannelCount());
	// Free the input model images
	image.reset();
	std::cout << "Number of components to be predicted: " << input.ComponentCount() << '\n';
	
	_predictingWatch.Start();
	
	for(size_t filenameIndex=0; filenameIndex!=_filenames.size(); ++filenameIndex)
	{
		const std::string& msName = _filenames[filenameIndex];
		
		size_t polIndex = 0;
		std::vector<MSProvider*> msProviders(_polarizations.size());
		for(PolarizationEnum pol : _polarizations)
		{
			msProviders[polIndex] = initializeMSProvider(filenameIndex, joinedChannelIndex, pol);
			++polIndex;
		}
		casa::MeasurementSet ms(msName);
		
		size_t nRow = ms.nrow();
		LMSPredicter predicter(ms, _threadCount);
		predicter.SetApplyBeam(_dftWithBeam);
		predicter.Input() = input;
		
		if(_dftWithBeam)
		{
			std::cout << "Converting model to absolute values...\n";
			predicter.Input().ConvertApparentToAbsolute(ms);
		}
		
		std::cout << "Creating row mapping...\n";
		std::vector<size_t> msToRowId;
		msProviders[0]->MakeMSRowToRowIdMapping(msToRowId, _globalSelection);
		
		ProgressBar progress("Predicting visibilities for " + msName);
		BandData band(ms.spectralWindow());
		predicter.Start();
		LMSPredicter::RowData row;
		ao::uvector<std::complex<float>> buffer[4];
		for(size_t p=0; p!=_polarizations.size(); ++p)
			buffer[p].assign(band.ChannelCount(), std::complex<float>(0.0));
		
		while(predicter.GetNextRow(row))
		{
			// Write to MS provider(s)
			polIndex = 0;
			for(PolarizationEnum pol : _polarizations)
			{
				for(size_t ch=0; ch!=band.ChannelCount(); ++ch)
				{
					switch(pol)
					{
						case Polarization::XX:
							buffer[polIndex][ch] = row.modelData[ch][0];
							break;
						case Polarization::XY:
							buffer[polIndex][ch] = row.modelData[ch][1];
							break;
						case Polarization::YX:
							buffer[polIndex][ch] = row.modelData[ch][2];
							break;
						case Polarization::YY:
							buffer[polIndex][ch] = row.modelData[ch][3];
							break;
						case Polarization::StokesI:
							buffer[polIndex][ch] =
								(row.modelData[ch][0] +
								row.modelData[ch][3])*0.5;
							break;
						default:
							throw std::runtime_error("Can't predict for this polarization at the moment");
					}
				}
				++polIndex;
			}
			
			boost::mutex::scoped_lock lock(predicter.IOMutex());
			for(size_t polIndex=0; polIndex!=_polarizations.size(); ++polIndex)
			{
				msProviders[polIndex]->WriteModel(msToRowId[row.rowIndex], buffer[polIndex].data());
			}
			lock.unlock();
				
			predicter.FinishRow(row);
			progress.SetProgress(row.rowIndex+1, nRow);
		}
		for(MSProvider* provider : msProviders)
			delete provider;
	}
	
	_predictingWatch.Pause();
}

void WSClean::initializeImageWeights(const MSSelection& partSelection)
{
	if(_mfsWeighting)
	{
		std::cout << "Reusing MFS weights for " << _weightMode.ToString() << " weighting.\n";
	}
	else {
		std::cout << "Precalculating weights for " << _weightMode.ToString() << " weighting... " << std::flush;
		_imageWeights.reset(new ImageWeights(_weightMode, _imgWidth, _imgHeight, _pixelScaleX, _pixelScaleY, _weightMode.SuperWeight()));
		for(size_t i=0; i!=_inversionAlgorithm->MeasurementSetCount(); ++i)
		{
			_imageWeights->Grid(_inversionAlgorithm->MeasurementSet(i), partSelection);
			if(_inversionAlgorithm->MeasurementSetCount() > 1)
				std::cout << i << ' ' << std::flush;
		}
		_imageWeights->FinishGridding();
		initializeWeightTapers();
		std::cout << "DONE\n";
	}
	_inversionAlgorithm->SetPrecalculatedWeightInfo(_imageWeights.get());
}

void WSClean::initializeWeightTapers()
{
	if(_minUVInLambda!=0.0)
		_imageWeights->SetMinUVRange(_minUVInLambda);
	if(_maxUVInLambda!=0.0)
		_imageWeights->SetMaxUVRange(_maxUVInLambda);
}

void WSClean::initializeMFSImageWeights()
{
	std::cout << "Precalculating MFS weights for " << _weightMode.ToString() << " weighting...\n";
	_imageWeights.reset(new ImageWeights(_weightMode, _imgWidth, _imgHeight, _pixelScaleX, _pixelScaleY, _weightMode.SuperWeight()));
	for(size_t i=0; i!=_filenames.size(); ++i)
	{
		if(_doReorder)
		{
			for(size_t ch=0; ch!=_channelsOut; ++ch)
			{
				MSSelection partSelection(_globalSelection);
				selectChannels(partSelection, ch, _channelsOut);
				PartitionedMS msProvider(_partitionedMSHandles[i], ch, *_polarizations.begin());
				_imageWeights->Grid(msProvider, partSelection);
			}
		}
		else {
			ContiguousMS msProvider(_filenames[i], _columnName, _globalSelection, *_polarizations.begin(), _mGain != 1.0);
			_imageWeights->Grid(msProvider,  _globalSelection);
			std::cout << '.' << std::flush;
		}
	}
	_imageWeights->FinishGridding();
	initializeWeightTapers();
}

void WSClean::freeCleanAlgorithms()
{
	for(std::vector<CleanAlgorithm*>::iterator caPtr=_cleanAlgorithms.begin(); caPtr!=_cleanAlgorithms.end(); ++caPtr)
		delete *caPtr;
	_cleanAlgorithms.clear();
}

void WSClean::initializeCleanAlgorithm()
{
	freeCleanAlgorithms();
	size_t count = 0;
	double beamSize = _inversionAlgorithm->BeamSize();
	if(_joinedPolarizationCleaning)
	{
		size_t polCount;
		if(Polarization::HasFullPolarization(_polarizations))
			polCount = 4;
		else if(Polarization::HasDualLinearPolarization(_polarizations))
			polCount = 2;
		else
			throw std::runtime_error("Joined polarization cleaning was requested, but can't find a compatible set of 2 or 4 pols to clean");
			
		_cleanAlgorithms.resize(1);
		if(_joinedFrequencyCleaning)
		{
			if(_multiscale)
			{
				if(polCount == 4)
				{
					_cleanAlgorithms[0] =
					new MultiScaleClean
					<clean_algorithms::MultiImageSet
					<clean_algorithms::PolarizedImageSet<4>>>(beamSize, _pixelScaleX, _pixelScaleY);
				}
				else {
					_cleanAlgorithms[0] =
					new MultiScaleClean
					<clean_algorithms::MultiImageSet<clean_algorithms::PolarizedImageSet<2>>>(beamSize, _pixelScaleX, _pixelScaleY);
				}
			}
			else {
				if(polCount == 4)
					_cleanAlgorithms[0] = new JoinedClean<clean_algorithms::MultiImageSet<clean_algorithms::PolarizedImageSet<4>>>();
				else
					_cleanAlgorithms[0] = new JoinedClean<clean_algorithms::MultiImageSet<clean_algorithms::PolarizedImageSet<2>>>();
			}
		}
		else {
			if(_multiscale)
			{
				if(polCount == 4)
					_cleanAlgorithms[0] = new MultiScaleClean<clean_algorithms::PolarizedImageSet<4>>(beamSize, _pixelScaleX, _pixelScaleY);
				else
					_cleanAlgorithms[0] = new MultiScaleClean<clean_algorithms::PolarizedImageSet<2>>(beamSize, _pixelScaleX, _pixelScaleY);
			}
			else
			{
				if(polCount == 4)
					_cleanAlgorithms[0] = new JoinedClean<clean_algorithms::PolarizedImageSet<4>>();
				else
					_cleanAlgorithms[0] = new JoinedClean<clean_algorithms::PolarizedImageSet<2>>();
			}
		}
		count = 1;
	}
	else {
		_cleanAlgorithms.resize(_polarizations.size());
		for(size_t p=0; p!=_polarizations.size(); ++p)
		{
			if(_joinedFrequencyCleaning)
			{
				if(_multiscale)
					_cleanAlgorithms[p] = new MultiScaleClean<clean_algorithms::MultiImageSet<clean_algorithms::SingleImageSet>>(beamSize, _pixelScaleX, _pixelScaleY);
				else
					_cleanAlgorithms[p] = new JoinedClean<clean_algorithms::MultiImageSet<clean_algorithms::SingleImageSet>>();
			}
			else {
				if(_multiscale)
					_cleanAlgorithms[p] = new MultiScaleClean<clean_algorithms::SingleImageSet>(beamSize, _pixelScaleX, _pixelScaleY);
				else
					_cleanAlgorithms[p] = new SimpleClean();
			}
		}
		count = _polarizations.size();
	}
	for(size_t p=0; p!=count; ++p)
	{
		_cleanAlgorithms[p]->SetMaxNIter(_nIter);
		_cleanAlgorithms[p]->SetThreshold(_threshold);
		_cleanAlgorithms[p]->SetSubtractionGain(_gain);
		_cleanAlgorithms[p]->SetStopGain(_mGain);
		_cleanAlgorithms[p]->SetCleanBorderRatio(_cleanBorderRatio);
		_cleanAlgorithms[p]->SetAllowNegativeComponents(_allowNegative);
		_cleanAlgorithms[p]->SetStopOnNegativeComponents(_stopOnNegative);
		_cleanAlgorithms[p]->SetResizePSF(_smallPSF);
		_cleanAlgorithms[p]->SetThreadCount(_threadCount);
		_cleanAlgorithms[p]->SetMultiscaleScaleBias(_multiscaleScaleBias);
		_cleanAlgorithms[p]->SetMultiscaleThresholdBias(_multiscaleThresholdBias);
		if(!_cleanAreasFilename.empty())
		{
			_cleanAreas.reset(new AreaSet());
			AreaParser parser;
			std::ifstream caFile(_cleanAreasFilename.c_str());
			parser.Parse(*_cleanAreas, caFile);
			_cleanAreas->SetImageProperties(_pixelScaleX, _pixelScaleY, _inversionAlgorithm->PhaseCentreRA(), _inversionAlgorithm->PhaseCentreDec(), _imgWidth, _imgHeight);
			_cleanAlgorithms[p]->SetCleanAreas(*_cleanAreas);
		}
	}
}

void WSClean::prepareInversionAlgorithm(PolarizationEnum polarization)
{
	static_cast<WSInversion&>(*_inversionAlgorithm).SetGridMode(_gridMode);
	_inversionAlgorithm->SetImageWidth(_imgWidth);
	_inversionAlgorithm->SetImageHeight(_imgHeight);
	_inversionAlgorithm->SetPixelSizeX(_pixelScaleX);
	_inversionAlgorithm->SetPixelSizeY(_pixelScaleY);
	if(_nWLayers != 0)
		_inversionAlgorithm->SetWGridSize(_nWLayers);
	else
		_inversionAlgorithm->SetNoWGridSize();
	_inversionAlgorithm->SetAntialiasingKernelSize(_antialiasingKernelSize);
	_inversionAlgorithm->SetOverSamplingFactor(_overSamplingFactor);
	_inversionAlgorithm->SetPolarization(polarization);
	_inversionAlgorithm->SetIsComplex(polarization == Polarization::XY || polarization == Polarization::YX);
	_inversionAlgorithm->SetDataColumnName(_columnName);
	_inversionAlgorithm->SetWeighting(_weightMode);
	_inversionAlgorithm->SetSelection(_currentPartSelection);
	_inversionAlgorithm->SetWLimit(_wLimit/100.0);
	_inversionAlgorithm->SetSmallInversion(_smallInversion);
}

void WSClean::checkPolarizations()
{
	if(_joinedPolarizationCleaning)
	{
		if(Polarization::HasFullPolarization(_polarizations))
		{
			if(_polarizations.size() != 4)
				throw std::runtime_error("Joined polarization cleaning requested. I found four suitable polarizations, but additional polarizations were specified");
		}
		else if(Polarization::HasDualPolarization(_polarizations))
		{
			if(_polarizations.size() != 2)
				throw std::runtime_error("Joined polarization cleaning requested. I found two suitable polarizations, but additional polarizations were specified that did not match");
		}
		else 
			throw std::runtime_error("Joined polarization cleaning requested, but neither 2 or 4 polarizations are imaged that are suitable for this");
	}
	else {
		if((_polarizations.count(Polarization::XY)!=0 || _polarizations.count(Polarization::YX)!=0) && _nIter!=0)
			throw std::runtime_error("You are imaging XY and/or YX polarizations and have enabled cleaning (niter!=0). This is not possible -- you have to specify '-joinpolarizations' or disable cleaning.");
	}
}

void WSClean::performReordering(bool isPredictMode)
{
	for(std::vector<std::string>::const_iterator i=_filenames.begin(); i != _filenames.end(); ++i)
	{
		_partitionedMSHandles.push_back(PartitionedMS::Partition(*i, _channelsOut, _globalSelection, _columnName, true, _mGain != 1.0 || isPredictMode, _polarizations, _temporaryDirectory));
	}
}

void WSClean::RunClean()
{
	// If no column specified, determine column to use
	if(_columnName.empty())
	{
		casa::MeasurementSet ms(_filenames.front());
		bool hasCorrected = ms.tableDesc().isColumn("CORRECTED_DATA");
		if(hasCorrected) {
			std::cout << "First measurement set has corrected data: tasks will be applied on the corrected data column.\n";
			_columnName = "CORRECTED_DATA";
		} else {
			std::cout << "No corrected data in first measurement set: tasks will be applied on the data column.\n";
			_columnName= "DATA";
		}
	}

	checkPolarizations();
	
	_doReorder = preferReordering();
	
	if(_doReorder) performReordering(false);
	
	_firstMSBand = MultiBandData(casa::MeasurementSet(_filenames[0]).spectralWindow(), casa::MeasurementSet(_filenames[0]).dataDescription());
	_weightPerChannel.assign(_channelsOut, 0);
	
	if(_mfsWeighting)
		initializeMFSImageWeights();
	
	if(_joinedFrequencyCleaning)
	{
		// run all frequencies
		runIndependentChannel(0);
	}
	else {
		for(size_t outChannelIndex=0; outChannelIndex!=_channelsOut; ++outChannelIndex)
		{
			runIndependentChannel(outChannelIndex);
		}
	}
	
	if(_channelsOut > 1)
	{
		for(std::set<PolarizationEnum>::const_iterator pol=_polarizations.begin(); pol!=_polarizations.end(); ++pol)
		{
			if(!(*pol == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
			{
				makeMFSImage("image.fits", *pol, false);
				if(_nIter > 0)
				{
					makeMFSImage("residual.fits", *pol, false);
					makeMFSImage("model.fits", *pol, false);
				}
				if(Polarization::IsComplex(*pol))
				{
					makeMFSImage("image.fits", *pol, true);
					if(_nIter > 0)
					{
					  makeMFSImage("residual.fits", *pol, true);
					  makeMFSImage("model.fits", *pol, true);
					}
				}
			}
		}
	}
}

void WSClean::RunPredict()
{
	if(_joinedFrequencyCleaning)
		throw std::runtime_error("Joined frequency cleaning specified for prediction: prediction doesn't clean, parameter invalid");
	if(_joinedPolarizationCleaning)
		throw std::runtime_error("Joined polarization cleaning specified for prediction: prediction doesn't clean, parameter invalid");
	
	_columnName = "DATA";
	
	checkPolarizations();
	
	_doReorder = preferReordering();
	
	if(_doReorder) performReordering(true);
	
	_firstMSBand = MultiBandData(casa::MeasurementSet(_filenames[0]).spectralWindow(), casa::MeasurementSet(_filenames[0]).dataDescription());
	
	for(size_t outChannelIndex=0; outChannelIndex!=_channelsOut; ++outChannelIndex)
	{
		predictChannel(outChannelIndex);
	}
}

void WSClean::selectChannels(MSSelection& selection, size_t outChannelIndex, size_t channelsOut)
{
	size_t startCh, endCh;
	if(selection.HasChannelRange())
	{
		startCh = selection.ChannelRangeStart();
		endCh = selection.ChannelRangeEnd();
	}
	else {
		startCh = 0;
		endCh = _firstMSBand.MaxChannels();
	}
	size_t newStart = startCh + (endCh - startCh) * outChannelIndex / channelsOut;
	size_t newEnd = startCh + (endCh - startCh) * (outChannelIndex+1) / channelsOut;
	selection.SetChannelRange(newStart, newEnd);
}

void WSClean::runIndependentChannel(size_t outChannelIndex)
{
	_inversionAlgorithm.reset(new WSInversion(&_imageAllocator, _threadCount, _memFraction, _absMemLimit));
	
	size_t  joinedChannelsOut;
	if(_joinedFrequencyCleaning)
		joinedChannelsOut = _channelsOut;
	else
		joinedChannelsOut = 1;
	
	_currentPartSelection = _globalSelection;
	if(!_joinedFrequencyCleaning)
		selectChannels(_currentPartSelection, outChannelIndex, _channelsOut);
	
	_modelImages.Initialize(_fitsWriter, _polarizations.size(), joinedChannelsOut, _prefixName + "-model", _imageAllocator);
	_residualImages.Initialize(_fitsWriter, _polarizations.size(), joinedChannelsOut, _prefixName + "-residual", _imageAllocator);
	_psfImages.Initialize(_fitsWriter, 1, joinedChannelsOut, _prefixName + "-psf", _imageAllocator);
	
	const std::string rootPrefix = _prefixName;
		
	for(size_t ch=0; ch!=joinedChannelsOut; ++ch)
	{
		const MSSelection selectionBefore = _currentPartSelection;
		size_t currentChannelIndex;
		if(_joinedFrequencyCleaning)
		{
			selectChannels(_currentPartSelection, ch, _channelsOut);
			currentChannelIndex = ch;
		}
		else {
			currentChannelIndex = outChannelIndex;
		}
		for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
		{
			runFirstInversion(currentChannelIndex, *curPol, ch);
		}
		_currentPartSelection = selectionBefore;
	}
	
	initializeCleanAlgorithm();

	initFitsWriter(_fitsWriter);
	setCleanParameters(_fitsWriter, *_cleanAlgorithms[0]);
	updateCleanParameters(_fitsWriter, 0, 0);
		
	if(_nIter > 0)
	{
		// Start major cleaning loop
		_majorIterationNr = 1;
		bool reachedMajorThreshold = false;
		do {
			performClean(outChannelIndex, reachedMajorThreshold, _majorIterationNr);
			
			if(_mGain != 1.0)
			{
				for(size_t ch=0; ch!=joinedChannelsOut; ++ch)
				{
					const MSSelection selectionBefore = _currentPartSelection;
					size_t currentChannelIndex;
					if(_joinedFrequencyCleaning)
					{
						selectChannels(_currentPartSelection, ch, _channelsOut);
						currentChannelIndex = ch;
					}
					else {
						currentChannelIndex = outChannelIndex;
					}
					if(_dftPrediction)
					{
						dftPredict(ch);
						for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
						{
							prepareInversionAlgorithm(*curPol);
							initializeCurMSProviders(currentChannelIndex, *curPol);
							if(*curPol == *_polarizations.begin())
								initializeImageWeights(_currentPartSelection);
		
							imageMainNonFirst(*curPol, ch);
							clearCurMSProviders();
						}
					}
					else {
						for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
						{
							prepareInversionAlgorithm(*curPol);
							
							initializeCurMSProviders(currentChannelIndex, *curPol);
							if(*curPol == *_polarizations.begin())
								initializeImageWeights(_currentPartSelection);
		
							predict(*curPol, ch);
							
							imageMainNonFirst(*curPol, ch);
							clearCurMSProviders();
						} // end of polarization loop
					}
					_currentPartSelection = selectionBefore;
				} // end of joined channels loop
				
				++_majorIterationNr;
			}
			
		} while(reachedMajorThreshold);
		
		std::cout << _majorIterationNr << " major iterations were performed.\n";
	}
	
	_inversionAlgorithm->FreeImagingData();
	
	// Restore model to residuals and save all images
	for(size_t ch=0; ch!=joinedChannelsOut; ++ch)
	{
		size_t currentChannelIndex = _joinedFrequencyCleaning ? ch : outChannelIndex;
		MSSelection curSelection(_globalSelection);
		selectChannels(curSelection, currentChannelIndex, _channelsOut);
		
		double freqLow, freqHigh;
		if(_firstMSBand.BandCount() == 1)
		{
			const BandData curBand = BandData(_firstMSBand.FirstBand(), curSelection.ChannelRangeStart(), curSelection.ChannelRangeEnd());
			freqLow = curBand.LowestFrequency();
			freqHigh = curBand.HighestFrequency();
		}
		else {
			freqLow = _firstMSBand.LowestFrequency();
			freqHigh = _firstMSBand.HighestFrequency();
		}
						
		for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
		{
			if(!(*curPol == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
			{
				std::vector<bool> imaginaryStates(1, false);
				if(*curPol == Polarization::XY || *curPol == Polarization::YX)
					imaginaryStates.push_back(true);
				
				for(std::vector<bool>::const_iterator isImaginary=imaginaryStates.begin(); isImaginary!=imaginaryStates.end(); ++isImaginary)
				{
					double* restoredImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
					_residualImages.Load(restoredImage, *curPol, currentChannelIndex, *isImaginary);
					if(_nIter != 0)
						writeFits("residual.fits", restoredImage, *curPol, currentChannelIndex, *isImaginary);
					double* modelImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
					_modelImages.Load(modelImage, *curPol, currentChannelIndex, *isImaginary);
					ModelRenderer renderer(_fitsWriter.RA(), _fitsWriter.Dec(), _pixelScaleX, _pixelScaleY, _fitsWriter.PhaseCentreDL(), _fitsWriter.PhaseCentreDM());
					if(_multiscale)
					{
						std::cout << "Rendering sources to restored image... " << std::flush;
						renderer.Restore(restoredImage, modelImage, _imgWidth, _imgHeight, _fitsWriter.BeamSizeMajorAxis(), _fitsWriter.BeamSizeMinorAxis(), _fitsWriter.BeamPositionAngle(), Polarization::StokesI);
						std::cout << "DONE\n";
					}
					else {
						Model model;
						// A model cannot hold instrumental pols (xx/xy/yx/yy), hence always use Stokes I here
						CleanAlgorithm::GetModelFromImage(model, modelImage, _imgWidth, _imgHeight, _fitsWriter.RA(), _fitsWriter.Dec(), _pixelScaleX, _pixelScaleY, _fitsWriter.PhaseCentreDL(), _fitsWriter.PhaseCentreDM(), 0.0, _fitsWriter.Frequency(), Polarization::StokesI);
						
						double beamMaj = _fitsWriter.BeamSizeMajorAxis();
						double beamMin = _fitsWriter.BeamSizeMinorAxis();
						double beamPA = _fitsWriter.BeamPositionAngle();
						if(beamMaj == beamMin) {
							std::cout << "Rendering " << model.SourceCount() << " circular sources to restored image... " << std::flush;
							renderer.Restore(restoredImage, _imgWidth, _imgHeight, model, beamMaj, freqLow, freqHigh, Polarization::StokesI);
						}
						else {
							std::cout << "Rendering " << model.SourceCount() << " elliptical sources to restored image... " << std::flush;
							renderer.Restore(restoredImage, _imgWidth, _imgHeight, model, beamMaj, beamMin, beamPA, freqLow, freqHigh, Polarization::StokesI);
						}
						std::cout << "DONE\n";
					}
					std::cout << "Writing restored image... " << std::flush;
					writeFits("image.fits", restoredImage, *curPol, currentChannelIndex, *isImaginary);
					std::cout << "DONE\n";
					_imageAllocator.Free(restoredImage);
					_imageAllocator.Free(modelImage);
				}
			}
		}
	}
	
	_imageAllocator.ReportStatistics();
	std::cout << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", cleaning: " << _cleaningWatch.ToString() << '\n';
	
	_prefixName = rootPrefix;
	
	// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
	_inversionAlgorithm.reset();
}

void WSClean::predictChannel(size_t outChannelIndex)
{
	_inversionAlgorithm.reset(new WSInversion(&_imageAllocator, _threadCount, _memFraction, _absMemLimit));
	
	_currentPartSelection = _globalSelection;
	selectChannels(_currentPartSelection, outChannelIndex, _channelsOut);
	
	_modelImages.Initialize(_fitsWriter, _polarizations.size(), 1, _prefixName + "-model", _imageAllocator);
	
	const std::string rootPrefix = _prefixName;
		
	for(std::set<PolarizationEnum>::const_iterator curPol=_polarizations.begin(); curPol!=_polarizations.end(); ++curPol)
	{
		// load image(s) from disk
		if(*curPol != Polarization::YX || _polarizations.count(Polarization::XY)==0)
		{
			FitsReader reader(getPrefix(*curPol, outChannelIndex, false) + "-model.fits");
			_fitsWriter = FitsWriter(reader);
			_modelImages.SetFitsWriter(_fitsWriter);
			std::cout << "Reading " << reader.Filename() << "...\n";
			double* buffer = _imageAllocator.Allocate(_imgWidth*_imgHeight);
			if(reader.ImageWidth()!=_imgWidth || reader.ImageHeight()!=_imgHeight)
				throw std::runtime_error("Inconsistent image size: input image did not match with specified dimensions.");
			reader.Read(buffer);
			for(size_t i=0; i!=_imgWidth*_imgHeight; ++i)
				if(!std::isfinite(buffer[i]))
					throw std::runtime_error("The input image contains non-finite values -- can't predict from an image with non-finite values");
			_modelImages.Store(buffer, *curPol, 0, false);
			if(*curPol == Polarization::XY)
			{
				reader = FitsReader(getPrefix(Polarization::XY, outChannelIndex, true) + "-model.fits");
				std::cout << "Reading " << reader.Filename() << "...\n";
				if(reader.ImageWidth()!=_imgWidth || reader.ImageHeight()!=_imgHeight)
					throw std::runtime_error("Inconsistent image size: input image did not match with specified dimensions.");
				reader.Read(buffer);
				_modelImages.Store(buffer, Polarization::XY, 0, true);
			}
			_imageAllocator.Free(buffer);
		}
		
		prepareInversionAlgorithm(*curPol);
		
		initializeCurMSProviders(outChannelIndex, *curPol);
		if(*curPol == *_polarizations.begin())
			initializeImageWeights(_currentPartSelection);

		predict(*curPol, 0);
		
		clearCurMSProviders();
	} // end of polarization loop
	
	_imageAllocator.ReportStatistics();
	std::cout << "Inversion: " << _inversionWatch.ToString() << ", prediction: " << _predictingWatch.ToString() << ", cleaning: " << _cleaningWatch.ToString() << '\n';
	
	_prefixName = rootPrefix;
	
	// Needs to be destructed before image allocator, or image allocator will report error caused by leaked memory
	_inversionAlgorithm.reset();
}

MSProvider* WSClean::initializeMSProvider(size_t filenameIndex, size_t currentChannelIndex, PolarizationEnum polarization)
{
	if(_doReorder)
		return new PartitionedMS(_partitionedMSHandles[filenameIndex], currentChannelIndex, polarization);
	else
		return new ContiguousMS(_filenames[filenameIndex], _columnName, _currentPartSelection, polarization, _mGain != 1.0);
}

void WSClean::initializeCurMSProviders(size_t currentChannelIndex, PolarizationEnum polarization)
{
	_inversionAlgorithm->ClearMeasurementSetList();
	for(size_t i=0; i != _filenames.size(); ++i)
	{
		MSProvider* msProvider = initializeMSProvider(i, currentChannelIndex, polarization);
		_inversionAlgorithm->AddMeasurementSet(msProvider);
		_currentPolMSes.push_back(msProvider);
	}
}

void WSClean::clearCurMSProviders()
{
	for(std::vector<MSProvider*>::iterator i=_currentPolMSes.begin(); i != _currentPolMSes.end(); ++i)
		delete *i;
	_currentPolMSes.clear();
}

void WSClean::runFirstInversion(size_t currentChannelIndex, PolarizationEnum polarization, size_t joinedChannelIndex)
{
	initializeCurMSProviders(currentChannelIndex, polarization);
	
	if(polarization == *_polarizations.begin())
		initializeImageWeights(_currentPartSelection);
	
	prepareInversionAlgorithm(polarization);
	
	const bool firstBeforePSF = _isFirstInversion;

	if((_nIter > 0 || _makePSF) && polarization == *_polarizations.begin())
		imagePSF(currentChannelIndex, joinedChannelIndex);
	
	initFitsWriter(_fitsWriter);
	_modelImages.SetFitsWriter(_fitsWriter);
	_residualImages.SetFitsWriter(_fitsWriter);
	
	imageMainFirst(polarization, joinedChannelIndex);
	
	_weightPerChannel[currentChannelIndex] = _inversionAlgorithm->ImageWeight();
	
	if(_isGriddingImageSaved && firstBeforePSF && _inversionAlgorithm->HasGriddingCorrectionImage())
		imageGridding();
	
	_isFirstInversion = false;
	
	// Set model to zero: already done if this is YX of XY/YX imaging combi
	if(!(polarization == Polarization::YX && _polarizations.count(Polarization::XY)!=0))
	{
		double* modelImage = _imageAllocator.Allocate(_imgWidth * _imgHeight);
		memset(modelImage, 0, _imgWidth * _imgHeight * sizeof(double));
		_modelImages.Store(modelImage, polarization, joinedChannelIndex, false);
		if(Polarization::IsComplex(polarization))
			_modelImages.Store(modelImage, polarization, joinedChannelIndex, true);
		_imageAllocator.Free(modelImage);
	}
	
	if(polarization == Polarization::XY && _polarizations.count(Polarization::YX)!=0)
	{ // Skip saving XY of XY/YX combi
	}
	else {
		PolarizationEnum savedPol = polarization;
		if(savedPol == Polarization::YX && _polarizations.count(Polarization::XY)!=0)
			savedPol = Polarization::XY;
		double* dirtyImage = _imageAllocator.Allocate(_imgWidth * _imgHeight);
		_residualImages.Load(dirtyImage, savedPol, joinedChannelIndex, false);
		std::cout << "Writing dirty image...\n";
		writeFits("dirty.fits", dirtyImage, savedPol, currentChannelIndex, false);
		if(Polarization::IsComplex(polarization))
		{
			_residualImages.Load(dirtyImage, savedPol, joinedChannelIndex, true);
			writeFits("dirty.fits", dirtyImage, savedPol, currentChannelIndex, true);
		}
		_imageAllocator.Free(dirtyImage);
	}
	
	clearCurMSProviders();
}

void WSClean::performClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr)
{
	std::cout << std::flush << " == Cleaning (" << majorIterationNr << ") ==\n";
	
	if(_joinedFrequencyCleaning)
	{
		if(_joinedPolarizationCleaning)
		{
			if(Polarization::HasFullPolarization(_polarizations))
				performJoinedPolFreqClean<4>(reachedMajorThreshold, majorIterationNr);
			else if(Polarization::HasDualLinearPolarization(_polarizations))
				performJoinedPolFreqClean<2>(reachedMajorThreshold, majorIterationNr);
			else throw std::runtime_error("Incompatible polarization combination for joined polarization cleaning");
		}
		else {
			if(_polarizations.size() == 1)
				performJoinedPolFreqClean<1>(reachedMajorThreshold, majorIterationNr);
			else
				throw std::runtime_error("Can only joinedly clean frequencies when cleaning the imaged polarizations simultaneously as well. Either provide the '-joinpolarizations' option to enable this, or image each polarization separately with separate WSClean commands.");
		}
	}
	else if(_joinedPolarizationCleaning) {
		if(Polarization::HasFullPolarization(_polarizations))
			performJoinedPolClean<4>(currentChannelIndex, reachedMajorThreshold, majorIterationNr);
		else if(Polarization::HasDualLinearPolarization(_polarizations))
			performJoinedPolClean<2>(currentChannelIndex, reachedMajorThreshold, majorIterationNr);
		else throw std::runtime_error("Incompatible polarization combination for joined polarization cleaning");
	}
	else {
		std::set<PolarizationEnum>::const_iterator pol = _polarizations.begin();
		for(size_t pIndex=0; pIndex!=_polarizations.size(); ++pIndex)
		{
			performSimpleClean(*_cleanAlgorithms[pIndex], currentChannelIndex, reachedMajorThreshold, majorIterationNr, *pol);
			++pol;
		}
	}
}

void WSClean::performSimpleClean(CleanAlgorithm& cleanAlgorithm, size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr, PolarizationEnum polarization)
{
	clean_algorithms::SingleImageSet
		residualImage(_imgWidth*_imgHeight, _imageAllocator),
		modelImage(_imgWidth*_imgHeight, _imageAllocator);
	double
		*psfImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
		
	_residualImages.Load(residualImage.Data(), polarization, 0, false);
	_modelImages.Load(modelImage.Data(), polarization, 0, false);
	_psfImages.Load(psfImage, polarization, 0, false);
	
	_cleaningWatch.Start();

	std::vector<double*> psfs(1, psfImage);
	TypedCleanAlgorithm<clean_algorithms::SingleImageSet>& tAlgorithm =
		static_cast<TypedCleanAlgorithm<clean_algorithms::SingleImageSet>&>(cleanAlgorithm);
	tAlgorithm.ExecuteMajorIteration(residualImage, modelImage, psfs, _imgWidth, _imgHeight, reachedMajorThreshold);
	_cleaningWatch.Pause();
	
	_modelImages.Store(modelImage.Data(), polarization, 0, false);
	_residualImages.Store(residualImage.Data(), polarization, 0, false);
	
	updateCleanParameters(_fitsWriter, cleanAlgorithm.IterationNumber(), majorIterationNr);
	
	if(majorIterationNr == 1)
	{
		if(_mGain != 1.0)
		{
			std::cout << "Writing first iteration image... " << std::flush;
			writeFits("first-residual.fits", residualImage.Data(), polarization, currentChannelIndex, false);
			std::cout << "DONE\n";
		}
	}
	if(!reachedMajorThreshold)
	{
		std::cout << "Writing model image... " << std::flush;
		writeFits("model.fits", modelImage.Data(), polarization, currentChannelIndex, false);
		std::cout << "DONE\n";
	}
	
	_imageAllocator.Free(psfImage);
}

template<size_t PolCount>
void WSClean::performJoinedPolClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr)
{
	typename JoinedClean<clean_algorithms::PolarizedImageSet<PolCount>>::ImageSet
		modelSet(_imgWidth*_imgHeight, _imageAllocator),
		residualSet(_imgWidth*_imgHeight, _imageAllocator);
	
	std::vector<double*> psfImages;
	double* psfImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
	_psfImages.Load(psfImage, *_polarizations.begin(), 0, false);
	psfImages.push_back(psfImage);
	
	bool hasStokesPols = (PolCount==4) && Polarization::HasFullStokesPolarization(_polarizations);
	if(hasStokesPols)
	{
		modelSet.LoadStokes(_modelImages, 0);
		residualSet.LoadStokes(_residualImages, 0);
	}
	else {
		modelSet.LoadLinear(_modelImages, 0);
		residualSet.LoadLinear(_residualImages, 0);
	}

	_cleaningWatch.Start();
	static_cast<TypedCleanAlgorithm<clean_algorithms::PolarizedImageSet<PolCount>>&>(*_cleanAlgorithms[0]).ExecuteMajorIteration(residualSet, modelSet, psfImages, _imgWidth, _imgHeight, reachedMajorThreshold);
	_cleaningWatch.Pause();
	
	_imageAllocator.Free(psfImage);
	if(hasStokesPols)
	{
		modelSet.StoreStokes(_modelImages, 0);
		residualSet.StoreStokes(_residualImages, 0);
	}
	else {
		modelSet.StoreLinear(_modelImages, 0);
		residualSet.StoreLinear(_residualImages, 0);
	}
	
	updateCleanParameters(_fitsWriter, _cleanAlgorithms[0]->IterationNumber(), majorIterationNr);
	
	const PolarizationEnum
		linPolsFour[4] = { Polarization::XX, Polarization::XY, Polarization::XY, Polarization::YY },
		linPolsTwo[2] = { Polarization::XX, Polarization::YY },
		stokesPols[4] = { Polarization::StokesI, Polarization::StokesQ, Polarization::StokesU, Polarization::StokesV },
		*pols = hasStokesPols ? stokesPols : (PolCount==4 ? linPolsFour : linPolsTwo);
	for(size_t i=0; i!=PolCount; ++i)
	{
		PolarizationEnum polarization = pols[i];
		bool isImaginary = (i==2) && !hasStokesPols;
		
		if(majorIterationNr == 1)
		{
			if(_mGain != 1.0)
			{
				std::cout << "Writing first iteration image... " << std::flush;
				writeFits("first-residual.fits", residualSet.GetImage(i), polarization, currentChannelIndex, isImaginary);
				std::cout << "DONE\n";
			}
		}
		if(!reachedMajorThreshold)
		{
			std::cout << "Writing model image... " << std::flush;
			writeFits("model.fits", modelSet.GetImage(i), polarization, currentChannelIndex, isImaginary);
			std::cout << "DONE\n";
		}
	}
}

template<size_t PolCount>
void WSClean::performJoinedPolFreqClean(bool& reachedMajorThreshold, size_t majorIterationNr)
{
	typename JoinedClean<clean_algorithms::MultiImageSet<clean_algorithms::PolarizedImageSet<PolCount>>>::ImageSet
		modelSet(_imgWidth*_imgHeight, _channelsOut, _imageAllocator),
		residualSet(_imgWidth*_imgHeight, _channelsOut, _imageAllocator);
	
	bool hasStokesPols = (PolCount==4) && Polarization::HasFullStokesPolarization(_polarizations);
	
	std::vector<double*> psfImages;
	for(size_t ch=0; ch!=_channelsOut; ++ch)
	{
		double* psfImage = _imageAllocator.Allocate(_imgWidth*_imgHeight);
		_psfImages.Load(psfImage, *_polarizations.begin(), ch, false);
		psfImages.push_back(psfImage);
		
		if(PolCount==1)
		{
			modelSet.Load(_modelImages, *_polarizations.begin(), ch);
			residualSet.Load(_residualImages, *_polarizations.begin(), ch);
		}
		else if(hasStokesPols)
		{
			modelSet.LoadStokes(_modelImages, ch);
			residualSet.LoadStokes(_residualImages, ch);
		}
		else {
			modelSet.LoadLinear(_modelImages, ch);
			residualSet.LoadLinear(_residualImages, ch);
		}
	}

	_cleaningWatch.Start();
	static_cast<TypedCleanAlgorithm<clean_algorithms::MultiImageSet<clean_algorithms::PolarizedImageSet<PolCount>>>&>(*_cleanAlgorithms[0]).ExecuteMajorIteration(residualSet, modelSet, psfImages, _imgWidth, _imgHeight, reachedMajorThreshold);
	_cleaningWatch.Pause();
	
	for(size_t ch=0; ch!=_channelsOut; ++ch)
	{
		_imageAllocator.Free(psfImages[ch]);
		if(PolCount == 1)
		{
			modelSet.Store(_modelImages, *_polarizations.begin(), ch);
			residualSet.Store(_residualImages, *_polarizations.begin(), ch);
		}
		else if(hasStokesPols)
		{
			modelSet.StoreStokes(_modelImages, ch);
			residualSet.StoreStokes(_residualImages, ch);
		}
		else {
			modelSet.StoreLinear(_modelImages, ch);
			residualSet.StoreLinear(_residualImages, ch);
		}
	}
	
	updateCleanParameters(_fitsWriter, _cleanAlgorithms[0]->IterationNumber(), majorIterationNr);
	
	const PolarizationEnum
		linPolsFour[4] = { Polarization::XX, Polarization::XY, Polarization::XY, Polarization::YY },
		linPolsTwo[2] = { Polarization::XX, Polarization::YY },
		stokesPols[4] = { Polarization::StokesI, Polarization::StokesQ, Polarization::StokesU, Polarization::StokesV },
		*pols = hasStokesPols ? stokesPols : (PolCount==4 ? linPolsFour : linPolsTwo);
	for(size_t ch=0; ch!=_channelsOut; ++ch)
	{
		for(size_t i=0; i!=PolCount; ++i)
		{
			PolarizationEnum polarization = pols[i];
			bool isImaginary = (i==2) && !hasStokesPols;
			
			if(majorIterationNr == 1)
			{
				if(_mGain != 1.0)
				{
					std::cout << "Writing first iteration image... " << std::flush;
					writeFits("first-residual.fits", residualSet.GetImage(i,ch), polarization, ch, isImaginary);
					std::cout << "DONE\n";
				}
			}
			if(!reachedMajorThreshold)
			{
				std::cout << "Writing model image... " << std::flush;
				writeFits("model.fits", modelSet.GetImage(i,ch), polarization, ch, isImaginary);
				std::cout << "DONE\n";
			}
		}
	}
}

void WSClean::makeMFSImage(const string& suffix, PolarizationEnum pol, bool isImaginary)
{
	double lowestFreq = 0.0, highestFreq = 0.0;
	const size_t size = _imgWidth * _imgHeight;
	ao::uvector<double> mfsImage(size, 0.0), addedImage(size), weightImage(size, 0.0);
	double weightSum = 0.0;
	FitsWriter writer;
	for(size_t ch=0; ch!=_channelsOut; ++ch)
	{
		const std::string name(getPrefix(pol, ch, isImaginary) + '-' + suffix);
		FitsReader reader(name);
		if(ch == 0)
		{
			writer = FitsWriter(reader);
			copyWSCleanKeywords(reader, writer);
			lowestFreq = reader.Frequency() - reader.Bandwidth()*0.5;
			highestFreq = reader.Frequency() + reader.Bandwidth()*0.5;
		}
		else {
			lowestFreq = std::min(lowestFreq, reader.Frequency() - reader.Bandwidth()*0.5);
			highestFreq = std::max(highestFreq, reader.Frequency() + reader.Bandwidth()*0.5);
		}
		double weight;
		if(!reader.ReadDoubleKeyIfExists("WSCIMGWG", weight))
			weight = 0.0;
		weightSum += weight;
		reader.Read(addedImage.data());
		for(size_t i=0; i!=size; ++i)
		{
			if(std::isfinite(addedImage[i])) {
				mfsImage[i] += addedImage[i] * weight;
				weightImage[i] += weight;
			}
		}
	}
	for(size_t i=0; i!=size; ++i)
		mfsImage[i] /= weightImage[i];
	writer.SetFrequency((lowestFreq+highestFreq)*0.5, highestFreq-lowestFreq);
	writer.SetNoBeamInfo();
	writer.SetExtraKeyword("WSCIMGWG", weightSum);
	writer.RemoveExtraKeyword("WSCCHANS");
	writer.RemoveExtraKeyword("WSCCHANE");
	writer.Write(getMFSPrefix(pol, isImaginary) + '-' + suffix, mfsImage.data());
}

void WSClean::writeFits(const string& suffix, const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary)
{
	double centreFrequency, bandwidth;
	if(_firstMSBand.BandCount() == 1)
	{
		MSSelection selection(_globalSelection);
		selectChannels(selection, channelIndex, _channelsOut);
		BandData band(_firstMSBand.FirstBand(), selection.ChannelRangeStart(), selection.ChannelRangeEnd());
		centreFrequency = band.CentreFrequency();
		bandwidth = band.Bandwidth();
	}
	else {
		centreFrequency = _firstMSBand.CentreFrequency();
		bandwidth = _firstMSBand.Bandwidth();
	}
	const std::string name(getPrefix(pol, channelIndex, isImaginary) + '-' + suffix);
	initFitsWriter(_fitsWriter);
	_fitsWriter.SetPolarization(pol);
	_fitsWriter.SetFrequency(centreFrequency, bandwidth);
	_fitsWriter.SetExtraKeyword("WSCIMGWG", _weightPerChannel[channelIndex]);
	size_t polIndex;
	if(_joinedPolarizationCleaning)
		polIndex = 0;
	else
		Polarization::TypeToIndex(pol, _polarizations, polIndex);
	if(_cleanAlgorithms.size()>polIndex && _cleanAlgorithms[polIndex] != 0)
	{
		setCleanParameters(_fitsWriter, *_cleanAlgorithms[polIndex]);
		updateCleanParameters(_fitsWriter, _cleanAlgorithms[polIndex]->IterationNumber(), _majorIterationNr);
	}
	_fitsWriter.Write(name, image);
}
