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

  ComponentList GetComponentList() const {
    return _parallelDeconvolution.GetComponentList(*_table, *_modelImages);
  }

  /**
   * @brief Exposes a const reference to either the first algorithm, or - in
   * case of a multiscale clean - the algorithm with the maximum number of scale
   * counts.
   */
  const DeconvolutionAlgorithm& MaxScaleCountAlgorithm() const {
    return _parallelDeconvolution.MaxScaleCountAlgorithm();
  }

  void Perform(bool& reachedMajorThreshold, size_t majorIterationNr);

  void InitializeDeconvolutionAlgorithm(
      std::unique_ptr<DeconvolutionTable> table,
      aocommon::PolarizationEnum psfPolarization, double beamSize,
      size_t threadCount);

  void FreeDeconvolutionAlgorithms();

  bool IsInitialized() const { return _parallelDeconvolution.IsInitialized(); }

  /// Return IterationNumber of the underlying \c DeconvolutionAlgorithm
  size_t IterationNumber() const;

  void SaveSourceList(const DeconvolutionTable& table,
                      long double phaseCentreRA, long double phaseCentreDec) {
    _parallelDeconvolution.SaveSourceList(table, phaseCentreRA, phaseCentreDec);
  }

  void SavePBSourceList(const DeconvolutionTable& table,
                        long double phaseCentreRA, long double phaseCentreDec) {
    _parallelDeconvolution.SavePBSourceList(table, phaseCentreRA,
                                            phaseCentreDec);
  }

 private:
  void readMask(const DeconvolutionTable& groupTable);

  const class Settings& _settings;

  std::unique_ptr<DeconvolutionTable> _table;

  ParallelDeconvolution _parallelDeconvolution;

  aocommon::UVector<bool> _cleanMask;

  bool _autoMaskIsFinished;
  aocommon::UVector<double> _channelFrequencies;
  aocommon::UVector<float> _channelWeights;
  aocommon::PolarizationEnum _psfPolarization;
  size_t _imgWidth, _imgHeight;
  aocommon::UVector<bool> _autoMask;
  double _beamSize, _pixelScaleX, _pixelScaleY;
};

#endif
