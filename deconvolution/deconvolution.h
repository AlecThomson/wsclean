#ifndef DECONVOLUTION_H
#define DECONVOLUTION_H

#include "paralleldeconvolution.h"

#include <aocommon/polarization.h>
#include <aocommon/uvector.h>

#include <cstring>

class DeconvolutionTable;
struct DeconvolutionTableEntry;

class Deconvolution {
 public:
  explicit Deconvolution(const class Settings& settings);
  ~Deconvolution();

  void Perform(const DeconvolutionTable& groupTable,
               bool& reachedMajorThreshold, size_t majorIterationNr);

  void InitializeDeconvolutionAlgorithm(
      const DeconvolutionTable& groupTable,
      aocommon::PolarizationEnum psfPolarization, double beamSize,
      size_t threadCount);

  void InitializeImages(class CachedImageSet& residuals, CachedImageSet& models,
                        CachedImageSet& psfs) {
    _residualImages = &residuals;
    _modelImages = &models;
    _psfImages = &psfs;
  }

  void FreeDeconvolutionAlgorithms() {
    _parallelDeconvolution.FreeDeconvolutionAlgorithms();
  }

  bool IsInitialized() const { return _parallelDeconvolution.IsInitialized(); }

  /// Return IterationNumber of the underlying \c DeconvolutionAlgorithm
  size_t IterationNumber() const;

  void SaveSourceList(const DeconvolutionTable& table,
                      long double phaseCentreRA, long double phaseCentreDec) {
    _parallelDeconvolution.SaveSourceList(*_modelImages, table, phaseCentreRA,
                                          phaseCentreDec);
  }

  void SavePBSourceList(const DeconvolutionTable& table,
                        long double phaseCentreRA, long double phaseCentreDec) {
    _parallelDeconvolution.SavePBSourceList(*_modelImages, table, phaseCentreRA,
                                            phaseCentreDec);
  }

 private:
  void correctChannelForPB(class ComponentList& list,
                           const DeconvolutionTableEntry& entry) const;

  void readMask(const DeconvolutionTable& groupTable);

  const class Settings& _settings;

  ParallelDeconvolution _parallelDeconvolution;

  aocommon::UVector<bool> _cleanMask;

  bool _autoMaskIsFinished;
  aocommon::UVector<double> _channelFrequencies;
  aocommon::UVector<float> _channelWeights;
  std::set<aocommon::PolarizationEnum> _polarizations;
  aocommon::PolarizationEnum _psfPolarization;
  size_t _imgWidth, _imgHeight;
  CachedImageSet *_psfImages, *_modelImages, *_residualImages;
  aocommon::UVector<bool> _autoMask;
  double _beamSize, _pixelScaleX, _pixelScaleY;
};

#endif
