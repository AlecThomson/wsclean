#ifndef PARALLEL_DECONVOLUTION_H
#define PARALLEL_DECONVOLUTION_H

#include "../system/fftwmanager.h"

#include "../structures/primarybeamimageset.h"

#include "controllablelog.h"

#include <aocommon/image.h>
#include <aocommon/uvector.h>

#include <memory>
#include <mutex>
#include <vector>

class ComponentList;

class ParallelDeconvolution {
 public:
  ParallelDeconvolution(const class Settings& settings);

  ~ParallelDeconvolution();

  class DeconvolutionAlgorithm& FirstAlgorithm() {
    return *_algorithms.front();
  }
  const class DeconvolutionAlgorithm& FirstAlgorithm() const {
    return *_algorithms.front();
  }

  ComponentList GetComponentList(const DeconvolutionTable& table,
                                 const CachedImageSet& modelImages) const;

  /**
   * @brief Same as @c FirstAlgorithm , except that for a multi-scale clean
   * the algorithm with the maximum number of scale counts is returned.
   */
  const DeconvolutionAlgorithm& MaxScaleCountAlgorithm() const;

  void SetAllocator(class ImageBufferAllocator* allocator) {
    _allocator = allocator;
  }

  void SetAlgorithm(std::unique_ptr<class DeconvolutionAlgorithm> algorithm);

  void SetRMSFactorImage(aocommon::Image&& image);

  void SetThreshold(double threshold);

  bool IsInitialized() const { return !_algorithms.empty(); }

  void SetAutoMaskMode(bool trackPerScaleMasks, bool usePerScaleMasks);

  void SetCleanMask(const bool* mask);

  void SetSpectrallyForcedImages(std::vector<aocommon::Image>&& images);

  void ExecuteMajorIteration(class ImageSet& dataImage,
                             class ImageSet& modelImage,
                             const aocommon::UVector<const float*>& psfImages,
                             bool& reachedMajorThreshold);

  void FreeDeconvolutionAlgorithms() {
    _algorithms.clear();
    _mask = nullptr;
  }

  void SaveSourceList(const DeconvolutionTable& table,
                      long double phaseCentreRA, long double phaseCentreDec);

  void SavePBSourceList(const DeconvolutionTable& table,
                        long double phaseCentreRA,
                        long double phaseCentreDec) const;

  class FFTWManager& GetFFTWManager() {
    return _fftwManager;
  }

 private:
  void executeParallelRun(class ImageSet& dataImage, class ImageSet& modelImage,
                          const aocommon::UVector<const float*>& psfImages,
                          bool& reachedMajorThreshold);

  struct SubImage {
    size_t index, x, y, width, height;
    // Mask to be used during deconvoution (combines user mask with the
    // boundary mask)
    aocommon::UVector<bool> mask;
    // Selects the pixels inside this subimage
    aocommon::UVector<bool> boundaryMask;
    double peak;
    bool reachedMajorThreshold;
  };

  void runSubImage(SubImage& subImg, ImageSet& dataImage,
                   const ImageSet& modelImage, ImageSet& resultModel,
                   const aocommon::UVector<const float*>& psfImages,
                   double majorIterThreshold, bool findPeakOnly,
                   std::mutex* mutex);

  void correctChannelForPB(class ComponentList& list,
                           const struct DeconvolutionTableEntry& entry) const;

  PrimaryBeamImageSet loadAveragePrimaryBeam(
      size_t imageIndex, const class DeconvolutionTable& table) const;

  void writeSourceList(const ComponentList& componentList,
                       const std::string& filename, long double phaseCentreRA,
                       long double phaseCentreDec) const;

  FFTWManager _fftwManager;
  std::vector<std::unique_ptr<class DeconvolutionAlgorithm>> _algorithms;
  FacetLogSet _logs;
  size_t _horImages, _verImages;
  const Settings& _settings;
  ImageBufferAllocator* _allocator;
  const bool* _mask;
  std::vector<aocommon::Image> _spectrallyForcedImages;
  bool _trackPerScaleMasks, _usePerScaleMasks;
  std::vector<aocommon::UVector<bool>> _scaleMasks;
  std::unique_ptr<class ComponentList> _componentList;
  aocommon::Image _rmsImage;
};

#endif
