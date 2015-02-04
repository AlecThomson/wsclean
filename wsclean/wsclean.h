#ifndef WSCLEAN_H
#define WSCLEAN_H

#include "msproviders/msprovider.h"
#include "msproviders/partitionedms.h"

#include "layeredimager.h"
#include "msselection.h"
#include "polarizationenum.h"
#include "weightmode.h"
#include "imagebufferallocator.h"
#include "stopwatch.h"
#include "cachedimageset.h"

#include <set>

class WSClean
{
public:
	WSClean();
	~WSClean();
	
	void SetImageSize(size_t width, size_t height) { _imgWidth = width; _imgHeight = height; }
	void SetPixelScale(double pixelScale) { _pixelScaleX = pixelScale; _pixelScaleY = pixelScale; }
	void SetNWlayers(size_t nWLayers) { _nWLayers = nWLayers; }
	void SetCleanGain(double gain) { _gain = gain; }
	void SetCleanMGain(double mGain) { _mGain = mGain; }
	void SetNIter(size_t nIter) { _nIter = nIter; }
	void SetCleanBorderRatio(double borderRatio) { _cleanBorderRatio = borderRatio; }
	void SetFitsMask(const std::string& fitsMask) { _fitsMask = fitsMask; }
	void SetCASAMask(const std::string& casaMask) { _casaMask = casaMask; }
	void SetThreshold(double threshold) { _threshold = threshold; }
	void SetColumnName(const std::string& columnName) { _columnName = columnName; }
	void SetPolarizations(const std::set<PolarizationEnum>& polarizations) { _polarizations = polarizations; }
	void SetAllowNegative(bool allowNegative) { _allowNegative = allowNegative; }
	void SetStopOnNegative(bool stopOnNegative) { _stopOnNegative = stopOnNegative; }
	void SetMakePSF(bool makePSF) { _makePSF = makePSF; }
	void SetPrefixName(const std::string& prefixName) { _prefixName = prefixName; }
	void SetGridMode(LayeredImager::GridModeEnum gridMode) { _gridMode = gridMode; }
	void SetSmallPSF(bool smallPSF) { _smallPSF = smallPSF; }
	void SetSmallInversion(bool smallInversion) { _smallInversion = smallInversion; }
	void SetIntervalSelection(size_t startTimestep, size_t endTimestep) {
		_globalSelection.SetInterval(startTimestep, endTimestep);
	}
	void SetChannelSelection(size_t startChannel, size_t endChannel) {
		_globalSelection.SetChannelRange(startChannel, endChannel);
	}
	void SetFieldSelection(size_t fieldId) {
		_globalSelection.SetFieldId(fieldId);
	}
	void SetChannelsOut(size_t channelsOut) { _channelsOut = channelsOut; }
	void SetJoinPolarizations(bool joinPolarizations) { _joinedPolarizationCleaning = joinPolarizations; }
	void SetJoinChannels(bool joinChannels) { _joinedFrequencyCleaning = joinChannels; }
	void SetIntervalCount(size_t intervalCount) { _intervalCount = intervalCount; }
	void SetMultiscale(bool multiscale) { _multiscale = multiscale; }
	void SetMultiscaleThresholdBias(double thresholdBias)
	{ _multiscaleThresholdBias = thresholdBias; }
	void SetMultiscaleScaleBias(double scaleBias)
	{ _multiscaleScaleBias = scaleBias; }
	void SetMFSWeighting(bool mfsWeighting) { _mfsWeighting = mfsWeighting; }
	void SetWeightMode(enum WeightMode::WeightingEnum weighting) {
		_weightMode.SetMode(WeightMode(weighting));
	}
	void SetBriggsWeighting(double robustness) {
		_weightMode.SetMode(WeightMode::Briggs(robustness));
	}
	void SetSuperWeight(double superWeight) { _weightMode.SetSuperWeight(superWeight); }
	void SetBeamSize(double major, double minor, double positionAngle) {
		_manualBeamMajorSize = major;
		_manualBeamMinorSize = minor;
		_manualBeamPA = positionAngle;
	}
	void SetFittedBeam(bool fittedBeam) { _fittedBeam = fittedBeam; }
	void SetCircularBeam(bool circularBeam) { _circularBeam = circularBeam; }
	void SetAntialiasingKernelSize(size_t kernelSize) { _antialiasingKernelSize = kernelSize; }
	void SetOversamplingFactor(size_t oversampling) { _overSamplingFactor = oversampling; }
	void SetThreadCount(size_t threadCount) { _threadCount = threadCount; }
	void SetTemporaryDirectory(const std::string& tempDir) { _temporaryDirectory = tempDir; }
	void SetForceReorder(bool forceReorder) { _forceReorder = forceReorder; }
	void SetForceNoReorder(bool forceNoReorder) { _forceNoReorder = forceNoReorder; }
	void SetMemFraction(double memFraction) { _memFraction = memFraction; }
	void SetMemAbsLimit(double absMemLimit) { _absMemLimit = absMemLimit; }
	void SetMinUVWInM(double minUVW) { _globalSelection.SetMinUVWInM(minUVW); }
	void SetMaxUVWInM(double maxUVW) { _globalSelection.SetMaxUVWInM(maxUVW); }
	void SetMinUVInLambda(double lambda) { _minUVInLambda = lambda; }
	void SetMaxUVInLambda(double lambda) { _maxUVInLambda = lambda; }
	void SetWLimit(double wLimit) { _wLimit = wLimit; }
	void SetCommandLine(const std::string& cmdLine) { _commandLine = cmdLine; }
	void SetSaveGriddingImage(bool isGriddingImageSaved) { _isGriddingImageSaved = isGriddingImageSaved; }
	void SetDFTPrediction(bool dftPrediction) { _dftPrediction = dftPrediction; }
	void SetDFTWithBeam(bool applyBeam) { _dftWithBeam = applyBeam; }
	void SetUseMoreSane(bool useMoreSane) { _useMoreSane = useMoreSane; }
	void SetMoreSaneLocation(const std::string& location) { _moreSaneLocation = location; }
	
	void AddInputMS(const std::string& msPath) { _filenames.push_back(msPath); }
	
	void RunClean();
	
	void RunPredict();
	
	bool JoinedFrequencyCleaning() const { return _joinedFrequencyCleaning; }
private:
	void runIndependentChannel(size_t outChannelIndex);
	void predictChannel(size_t outChannelIndex);
	
	void runFirstInversion(size_t outChannelIndex, PolarizationEnum polarization, size_t joinedChannelIndex);
	void performClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr);
	void performSimpleClean(class CleanAlgorithm& cleanAlgorithm, size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr, PolarizationEnum polarization);
	template<size_t PolCount>
	void performJoinedPolClean(size_t currentChannelIndex, bool& reachedMajorThreshold, size_t majorIterationNr);
	template<size_t PolCount>
	void performJoinedPolFreqClean(bool& reachedMajorThreshold, size_t majorIterationNr);
	void prepareInversionAlgorithm(PolarizationEnum polarization);
	
	void checkPolarizations();
	void performReordering(bool isPredictMode);
	
	void initFitsWriter(class FitsWriter& writer);
	void copyWSCleanKeywords(FitsReader& reader, FitsWriter& writer);
	//void copyDoubleKeywordIfExists(FitsReader& reader, FitsWriter& writer, const char* keywordName);
	void setCleanParameters(class FitsWriter& writer, const class CleanAlgorithm& clean);
	void updateCleanParameters(class FitsWriter& writer, size_t minorIterationNr, size_t majorIterationNr);
	void initializeWeightTapers();
	void initializeImageWeights(const MSSelection& partSelection);
	void initializeMFSImageWeights();
	void initializeCleanAlgorithm();
	void freeCleanAlgorithms();
	MSProvider* initializeMSProvider(size_t filenameIndex, size_t currentChannelIndex, PolarizationEnum polarization);
	void initializeCurMSProviders(size_t currentChannelIndex, PolarizationEnum polarization);
	void clearCurMSProviders();
	void storeAndCombineXYandYX(CachedImageSet& dest, PolarizationEnum polarization, size_t joinedChannelIndex, bool isImaginary, const double* image);
	void selectChannels(MSSelection& selection, size_t outChannelIndex, size_t channelsOut);
	MSSelection selectInterval(MSSelection& fullSelection);
	
	void imagePSF(size_t currentChannelIndex, size_t joinedChannelIndex);
	void imageGridding();
	void imageMainFirst(PolarizationEnum polarization, size_t joinedChannelIndex);
	void imageMainNonFirst(PolarizationEnum polarization, size_t joinedChannelIndex);
	void predict(PolarizationEnum polarization, size_t joinedChannelIndex);
	void dftPredict(size_t joinedChannelIndex);
	
	void makeMFSImage(const string& suffix, PolarizationEnum pol, bool isImaginary);
	void writeFits(const string& suffix, const double* image, PolarizationEnum pol, size_t channelIndex, bool isImaginary);
	
	std::string fourDigitStr(size_t val) const
	{
		std::ostringstream str;
		if(val < 1000) str << '0';
		if(val < 100) str << '0';
		if(val < 10) str << '0';
		str << val;
		return str.str();
	}
	
	std::string getPSFPrefix(size_t channelIndex) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _prefixName;
		if(_intervalCount != 1)
			partPrefixNameStr << "-t" << fourDigitStr(_currentIntervalIndex);
		if(_channelsOut != 1)
			partPrefixNameStr << '-' << fourDigitStr(channelIndex);
		return partPrefixNameStr.str();
	}
	
	std::string getPrefix(PolarizationEnum polarization, size_t channelIndex, bool isImaginary) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _prefixName;
		if(_intervalCount != 1)
			partPrefixNameStr << "-t" << fourDigitStr(_currentIntervalIndex);
		if(_channelsOut != 1)
			partPrefixNameStr << '-' << fourDigitStr(channelIndex);
		if(_polarizations.size() != 1)
		{
			partPrefixNameStr << '-' << Polarization::TypeToShortString(polarization);
			if(isImaginary)
				partPrefixNameStr << 'i';
		}
		return partPrefixNameStr.str();
	}
	
	std::string getMFSPrefix(PolarizationEnum polarization, bool isImaginary) const
	{
		std::ostringstream partPrefixNameStr;
		partPrefixNameStr << _prefixName;
		if(_intervalCount != 1)
			partPrefixNameStr << "-t" << fourDigitStr(_currentIntervalIndex);
		if(_channelsOut != 1)
			partPrefixNameStr << "-MFS";
		if(_polarizations.size() != 1)
		{
			partPrefixNameStr << '-' << Polarization::TypeToShortString(polarization);
			if(isImaginary)
				partPrefixNameStr << 'i';
		}
		return partPrefixNameStr.str();
	}
	
	bool preferReordering() const
	{
		return ((_channelsOut != 1) || (_polarizations.size()>=4) || _forceReorder) && !_forceNoReorder;
	}
	
	size_t _imgWidth, _imgHeight, _channelsOut, _intervalCount;
	double _pixelScaleX, _pixelScaleY, _threshold, _gain, _mGain, _cleanBorderRatio;
	std::string _fitsMask, _casaMask;
	double _manualBeamMajorSize, _manualBeamMinorSize, _manualBeamPA;
	bool _fittedBeam, _circularBeam;
	double _memFraction, _absMemLimit, _minUVInLambda, _maxUVInLambda, _wLimit, _multiscaleThresholdBias, _multiscaleScaleBias;
	size_t _nWLayers, _nIter, _antialiasingKernelSize, _overSamplingFactor, _threadCount;
	MSSelection _globalSelection, _currentPartSelection;
	std::string _columnName;
	std::set<PolarizationEnum> _polarizations;
	WeightMode _weightMode;
	std::string _prefixName;
	bool _allowNegative, _smallPSF, _smallInversion, _stopOnNegative, _useMoreSane, _makePSF, _isGriddingImageSaved, _dftPrediction, _dftWithBeam;
	std::string _temporaryDirectory, _moreSaneLocation;
	bool _forceReorder, _forceNoReorder, _joinedPolarizationCleaning, _joinedFrequencyCleaning, _mfsWeighting, _multiscale;
	enum LayeredImager::GridModeEnum _gridMode;
	std::vector<std::string> _filenames;
	std::string _commandLine;
	std::vector<size_t> _weightPerChannel;
	
	std::unique_ptr<class InversionAlgorithm> _inversionAlgorithm;
	std::unique_ptr<class ImageWeights> _imageWeights;
	std::vector<class CleanAlgorithm*> _cleanAlgorithms;
	ao::uvector<bool> _cleanMask;
	ImageBufferAllocator<double> _imageAllocator;
	Stopwatch _inversionWatch, _predictingWatch, _deconvolutionWatch;
	bool _isFirstInversion, _doReorder;
	size_t _currentIntervalIndex, _majorIterationNr;
	CachedImageSet _psfImages, _modelImages, _residualImages;
	std::vector<PartitionedMS::Handle> _partitionedMSHandles;
	FitsWriter _fitsWriter;
	std::vector<MSProvider*> _currentPolMSes;
	MultiBandData _firstMSBand;
};

#endif
