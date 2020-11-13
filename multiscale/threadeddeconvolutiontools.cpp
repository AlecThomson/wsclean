#include "threadeddeconvolutiontools.h"
#include "multiscaletransforms.h"

#include "../deconvolution/peakfinder.h"
#include "../deconvolution/simpleclean.h"

ThreadedDeconvolutionTools::ThreadedDeconvolutionTools(size_t threadCount)
    : _taskLanes(threadCount),
      _resultLanes(threadCount),
      _threadCount(threadCount) {
  for (size_t i = 0; i != _threadCount; ++i) {
    _taskLanes[i] = new aocommon::Lane<ThreadTask*>(1);
    _resultLanes[i] = new aocommon::Lane<ThreadResult*>(1);
    _threadGroup.emplace_back(&ThreadedDeconvolutionTools::threadFunc, this,
                              _taskLanes[i], _resultLanes[i]);
  }
}

ThreadedDeconvolutionTools::~ThreadedDeconvolutionTools() {
  for (size_t i = 0; i != _threadCount; ++i) {
    _taskLanes[i]->write_end();
  }

  for (std::thread& t : _threadGroup) t.join();

  for (size_t i = 0; i != _threadCount; ++i) {
    delete _taskLanes[i];
    delete _resultLanes[i];
  }
}

void ThreadedDeconvolutionTools::SubtractImage(float* image, const float* psf,
                                               size_t width, size_t height,
                                               size_t x, size_t y,
                                               float factor) {
  for (size_t thr = 0; thr != _threadCount; ++thr) {
    SubtractionTask* task = new SubtractionTask();
    task->image = image;
    task->psf = psf;
    task->width = width;
    task->height = height;
    task->x = x;
    task->y = y;
    task->factor = factor;
    task->startY = height * thr / _threadCount;
    task->endY = height * (thr + 1) / _threadCount;
    if (thr == _threadCount - 1) {
      ThreadResult* result = (*task)();
      delete result;
    } else {
      _taskLanes[thr]->write(task);
    }
  }
  for (size_t thr = 0; thr != _threadCount - 1; ++thr) {
    ThreadResult* result = 0;
    _resultLanes[thr]->read(result);
    delete result;
  }
}

ThreadedDeconvolutionTools::ThreadResult*
ThreadedDeconvolutionTools::SubtractionTask::operator()() {
  SimpleClean::PartialSubtractImage(image, psf, width, height, x, y, factor,
                                    startY, endY);
  return 0;
}

void ThreadedDeconvolutionTools::MultiScaleTransform(
    MultiScaleTransforms* msTransforms, std::vector<ImageF>& images,
    ImageF& scratch, float scale) {
  size_t imageIndex = 0;
  size_t nextThread = 0;
  msTransforms->PrepareTransform(scratch.data(), scale);
  while (imageIndex < images.size()) {
    FinishMultiScaleTransformTask* task = new FinishMultiScaleTransformTask();
    task->msTransforms = msTransforms;
    task->image = &images[imageIndex];
    task->kernel = &scratch;
    _taskLanes[nextThread]->write(task);

    ++nextThread;
    if (nextThread == _threadCount) {
      for (size_t thr = 0; thr != nextThread; ++thr) {
        ThreadResult* result = 0;
        _resultLanes[thr]->read(result);
        delete result;
      }
      nextThread = 0;
    }
    ++imageIndex;
  }
  for (size_t thr = 0; thr != nextThread; ++thr) {
    ThreadResult* result = nullptr;
    _resultLanes[thr]->read(result);
    delete result;
  }
}

ThreadedDeconvolutionTools::ThreadResult*
ThreadedDeconvolutionTools::FinishMultiScaleTransformTask::operator()() {
  msTransforms->FinishTransform(image->data(), kernel->data());
  return nullptr;
}

void ThreadedDeconvolutionTools::MultiScaleTransform(
    MultiScaleTransforms* msTransforms, std::vector<ImageF>& images,
    aocommon::UVector<float> scales) {
  size_t imageIndex = 0;
  size_t nextThread = 0;

  size_t scratchCount = std::min(images.size(), _threadCount);
  std::unique_ptr<ImageF::Ptr[]> scratchImages(new ImageF::Ptr[scratchCount]);
  for (size_t i = 0; i != scratchCount; ++i)
    scratchImages[i] =
        ImageF::Make(msTransforms->Width(), msTransforms->Height());

  while (imageIndex < images.size()) {
    MultiScaleTransformTask* task = new MultiScaleTransformTask();
    task->msTransforms = msTransforms;
    task->image = &images[imageIndex];
    task->scratch = scratchImages[nextThread].get();
    task->scale = scales[imageIndex];
    _taskLanes[nextThread]->write(task);

    ++nextThread;
    if (nextThread == _threadCount) {
      for (size_t thr = 0; thr != nextThread; ++thr) {
        ThreadResult* result = 0;
        _resultLanes[thr]->read(result);
        delete result;
      }
      nextThread = 0;
    }
    ++imageIndex;
  }
  for (size_t thr = 0; thr != nextThread; ++thr) {
    ThreadResult* result = nullptr;
    _resultLanes[thr]->read(result);
    delete result;
  }
}

ThreadedDeconvolutionTools::ThreadResult*
ThreadedDeconvolutionTools::MultiScaleTransformTask::operator()() {
  msTransforms->Transform(*image, *scratch, scale);
  return nullptr;
}

void ThreadedDeconvolutionTools::FindMultiScalePeak(
    MultiScaleTransforms* msTransforms, const ImageF& image,
    const aocommon::UVector<float>& scales,
    std::vector<ThreadedDeconvolutionTools::PeakData>& results,
    bool allowNegativeComponents, const bool* mask,
    const std::vector<aocommon::UVector<bool>>& scaleMasks, float borderRatio,
    const ImageF& rmsFactorImage, bool calculateRMS) {
  size_t imageIndex = 0;
  size_t nextThread = 0;
  size_t resultIndex = 0;

  results.resize(scales.size());

  size_t size = std::min(scales.size(), _threadCount);
  std::unique_ptr<ImageF::Ptr[]> imageData(new ImageF::Ptr[size]);
  std::unique_ptr<ImageF::Ptr[]> scratchData(new ImageF::Ptr[size]);
  for (size_t i = 0; i != size; ++i) {
    imageData[i] = ImageF::Make(msTransforms->Width(), msTransforms->Height());
    scratchData[i] =
        ImageF::Make(msTransforms->Width(), msTransforms->Height());
  }

  while (imageIndex < scales.size()) {
    FindMultiScalePeakTask* task = new FindMultiScalePeakTask();
    task->msTransforms = msTransforms;
    (*imageData[nextThread]) = image;
    task->image = imageData[nextThread].get();
    task->scratch = scratchData[nextThread].get();
    task->scale = scales[imageIndex];
    task->allowNegativeComponents = allowNegativeComponents;
    if (scaleMasks.empty())
      task->mask = mask;
    else
      task->mask = scaleMasks[imageIndex].data();
    task->borderRatio = borderRatio;
    task->calculateRMS = calculateRMS;
    task->rmsFactorImage = &rmsFactorImage;
    _taskLanes[nextThread]->write(task);

    ++nextThread;
    if (nextThread == _threadCount) {
      for (size_t thr = 0; thr != nextThread; ++thr) {
        ThreadResult* result = nullptr;
        _resultLanes[thr]->read(result);
        results[resultIndex].normalizedValue =
            static_cast<FindMultiScalePeakResult*>(result)->normalizedValue;
        results[resultIndex].unnormalizedValue =
            static_cast<FindMultiScalePeakResult*>(result)->unnormalizedValue;
        results[resultIndex].x =
            static_cast<FindMultiScalePeakResult*>(result)->x;
        results[resultIndex].y =
            static_cast<FindMultiScalePeakResult*>(result)->y;
        results[resultIndex].rms =
            static_cast<FindMultiScalePeakResult*>(result)->rms;
        delete result;
        ++resultIndex;
      }
      nextThread = 0;
    }
    ++imageIndex;
  }
  for (size_t thr = 0; thr != nextThread; ++thr) {
    ThreadResult* result = nullptr;
    _resultLanes[thr]->read(result);
    results[resultIndex].unnormalizedValue =
        static_cast<FindMultiScalePeakResult*>(result)->unnormalizedValue;
    results[resultIndex].normalizedValue =
        static_cast<FindMultiScalePeakResult*>(result)->normalizedValue;
    results[resultIndex].x = static_cast<FindMultiScalePeakResult*>(result)->x;
    results[resultIndex].y = static_cast<FindMultiScalePeakResult*>(result)->y;
    results[resultIndex].rms =
        static_cast<FindMultiScalePeakResult*>(result)->rms;
    delete result;
    ++resultIndex;
  }
}

ThreadedDeconvolutionTools::ThreadResult*
ThreadedDeconvolutionTools::FindMultiScalePeakTask::operator()() {
  msTransforms->Transform(*image, *scratch, scale);
  const size_t width = msTransforms->Width(), height = msTransforms->Height(),
               scaleBorder = size_t(ceil(scale * 0.5)),
               horBorderSize =
                   std::max<size_t>(round(width * borderRatio), scaleBorder),
               vertBorderSize =
                   std::max<size_t>(round(height * borderRatio), scaleBorder);
  FindMultiScalePeakResult* result = new FindMultiScalePeakResult();
  if (calculateRMS)
    result->rms = RMS(*image, width * height);
  else
    result->rms = -1.0;
  if (rmsFactorImage->empty()) {
    if (mask == 0)
      result->unnormalizedValue = PeakFinder::Find(
          image->data(), width, height, result->x, result->y,
          allowNegativeComponents, 0, height, horBorderSize, vertBorderSize);
    else
      result->unnormalizedValue =
          PeakFinder::FindWithMask(image->data(), width, height, result->x,
                                   result->y, allowNegativeComponents, 0,
                                   height, mask, horBorderSize, vertBorderSize);

    result->normalizedValue = result->unnormalizedValue;
  } else {
    for (size_t i = 0; i != rmsFactorImage->size(); ++i)
      (*scratch)[i] = (*image)[i] * (*rmsFactorImage)[i];

    if (mask == 0)
      result->unnormalizedValue = PeakFinder::Find(
          scratch->data(), width, height, result->x, result->y,
          allowNegativeComponents, 0, height, horBorderSize, vertBorderSize);
    else
      result->unnormalizedValue =
          PeakFinder::FindWithMask(scratch->data(), width, height, result->x,
                                   result->y, allowNegativeComponents, 0,
                                   height, mask, horBorderSize, vertBorderSize);

    if (result->unnormalizedValue) {
      result->normalizedValue =
          (*result->unnormalizedValue) /
          (*rmsFactorImage)[result->x + result->y * width];
    } else {
      result->normalizedValue = boost::optional<float>();
    }
  }
  return result;
}
