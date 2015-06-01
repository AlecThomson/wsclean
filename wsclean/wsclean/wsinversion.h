#ifndef WS_INVERSION
#define WS_INVERSION

#include "inversionalgorithm.h"
#include "layeredimager.h"

#include "../lane.h"
#include "../multibanddata.h"

#include <complex>
#include <memory>

#include <casa/Arrays/Array.h>
#include <tables/Tables/ArrayColumn.h>

namespace casa {
	class MeasurementSet;
}
template<typename NumType>
class ImageBufferAllocator;

class WSInversion : public InversionAlgorithm
{
	public:
		WSInversion(class ImageBufferAllocator<double>* imageAllocator, size_t threadCount, double memFraction, double absMemLimit);
	
		virtual void Invert();
		
		virtual void Predict(double* image) { Predict(image, 0); }
		virtual void Predict(double* real, double* imaginary);
		
		virtual double *ImageRealResult() const { return _imager->RealImage(); }
		virtual double *ImageImaginaryResult() const {
			if(!IsComplex())
				throw std::runtime_error("No imaginary result available for non-complex inversion");
			return _imager->ImaginaryImage();
		}
		virtual double PhaseCentreRA() const { return _phaseCentreRA; }
		virtual double PhaseCentreDec() const { return _phaseCentreDec; }
		virtual double HighestFrequencyChannel() const { return _freqHigh; }
		virtual double LowestFrequencyChannel() const { return _freqLow; }
		virtual double BandStart() const { return _bandStart; }
		virtual double BandEnd() const { return _bandEnd; }
		virtual double BeamSize() const { return _beamSize; }
		virtual double StartTime() const { return _startTime; }
		virtual bool HasDenormalPhaseCentre() const { return _denormalPhaseCentre; }
		virtual double PhaseCentreDL() const { return _phaseCentreDL; }
		virtual double PhaseCentreDM() const { return _phaseCentreDM; }
		virtual double ImageWeight() const { return _totalWeight; }
		
		enum LayeredImager::GridModeEnum GridMode() const { return _gridMode; }
		void SetGridMode(LayeredImager::GridModeEnum gridMode) { _gridMode = gridMode; }
		
		virtual bool HasGriddingCorrectionImage() const { return _gridMode == LayeredImager::KaiserBessel; }
		virtual void GetGriddingCorrectionImage(double *image) const { _imager->GetGriddingCorrectionImage(image); }
		
		size_t ActualInversionWidth() const { return _actualInversionWidth; }
		size_t ActualInversionHeight() const { return _actualInversionHeight; }
		
		virtual void FreeImagingData()
		{
			_imager.reset();
		}
	private:
		struct InversionWorkItem
		{
			double u, v, w;
			size_t dataDescId;
			std::complex<float> *data;
		};
		struct InversionWorkSample
		{
			double uInLambda, vInLambda, wInLambda;
			std::complex<float> sample;
		};
		struct PredictionWorkItem
		{
			double u, v, w;
			std::complex<float> *data;
			size_t rowId, dataDescId;
		};
		
		struct MSData
		{
			public:
				MSData();
				~MSData();
				class MSProvider *msProvider;
				MultiBandData bandData;
				size_t startChannel, endChannel;
				size_t matchingRows, totalRowsProcessed;
				double minW, maxW;
				size_t rowStart, rowEnd;
			
				MultiBandData SelectedBand() const { return MultiBandData(bandData, startChannel, endChannel); }
			private:
				MSData(const MSData &source);
				
				void operator=(const MSData &source);
		};
		
		void initializeMeasurementSet(size_t msIndex, MSData &msData);
		void gridMeasurementSet(MSData &msData);
		void countSamplesPerLayer(MSData &msData);

		void predictMeasurementSet(MSData &msData);

		void workThread(ao::lane<InversionWorkItem>* workLane)
		{
			InversionWorkItem workItem;
			while(workLane->read(workItem))
			{
				_imager->AddData(workItem.data, workItem.dataDescId, workItem.u, workItem.v, workItem.w);
				delete[] workItem.data;
			}
		}
		
		void workThreadParallel(const MultiBandData* selectedBand);
		void workThreadPerSample(ao::lane<InversionWorkSample>* workLane);
		
		void predictCalcThread(ao::lane<PredictionWorkItem>* inputLane, ao::lane<PredictionWorkItem>* outputLane);
		void predictWriteThread(ao::lane<PredictionWorkItem>* samplingWorkLane, const MSData* msData);
		static void rotateVisibilities(const BandData &bandData, double shiftFactor, std::complex<float>* dataIter);

		std::unique_ptr<LayeredImager> _imager;
		std::unique_ptr<ao::lane<InversionWorkItem>> _inversionWorkLane;
		double _phaseCentreRA, _phaseCentreDec, _phaseCentreDL, _phaseCentreDM;
		bool _denormalPhaseCentre, _hasFrequencies;
		double _freqHigh, _freqLow;
		double _bandStart, _bandEnd;
		double _beamSize;
		double _totalWeight;
		double _startTime;
		LayeredImager::GridModeEnum _gridMode;
		size_t _cpuCount, _laneBufferSize;
		int64_t _memSize;
		ImageBufferAllocator<double>* _imageBufferAllocator;
		size_t _actualInversionWidth, _actualInversionHeight;
		double _actualPixelSizeX, _actualPixelSizeY;
};

#endif
