#include "imageset.h"
#include "spectralfitter.h"

#include "../math/nlplfitter.h"

#include "../io/cachedimageset.h"
#include "../io/logger.h"

#include "../structures/primarybeam.h"
#include "../structures/primarybeamimageset.h"

#include <aocommon/staticfor.h>

ImageSet::ImageSet(const ImagingTable* table, const class Settings& settings)
    : _images(),
      _width(0),
      _height(0),
      _channelsInDeconvolution((settings.deconvolutionChannelCount == 0)
                                   ? table->SquaredGroups().size()
                                   : settings.deconvolutionChannelCount),
      _squareJoinedChannels(settings.squaredJoins),
      _imagingTable(*table),
      _imageIndexToPSFIndex(),
      _linkedPolarizations(settings.linkedPolarizations),
      _settings(settings) {
  size_t nPol = table->SquaredGroups().front().size();
  size_t nImages = nPol * _channelsInDeconvolution;
  _images.resize(nImages);
  _imageIndexToPSFIndex.resize(nImages);

  initializePolFactor();
  initializeIndices();
  CalculateDeconvolutionFrequencies(*table, _frequencies, _weights,
                                    _channelsInDeconvolution);
}

ImageSet::ImageSet(const ImagingTable* table, const class Settings& settings,
                   size_t width, size_t height)
  : ImageSet(table, settings) {
  _width = width;
  _height = height;
  allocateImages();
}

void ImageSet::initializeIndices() {
  size_t lastDeconvolutionChannel = 0;
  size_t deconvolutionChannelStartIndex = 0, lastOutChannel = 0;
  size_t imgIndex = 0;
  for (const ImagingTableEntry& entry : _imagingTable) {
    size_t outChannel = entry.outputChannelIndex;
    size_t chIndex = (outChannel * _channelsInDeconvolution) /
                     _imagingTable.SquaredGroups().size();
    if (outChannel != lastOutChannel && chIndex == lastDeconvolutionChannel) {
      // New output channel maps to an earlier deconvolution channel:
      // start at index of previous deconvolution channel
      imgIndex = deconvolutionChannelStartIndex;
    }
    if (chIndex != lastDeconvolutionChannel) {
      deconvolutionChannelStartIndex = imgIndex;
    }
    _entryIndexToImageIndex.emplace(entry.index, imgIndex);
    lastOutChannel = outChannel;
    lastDeconvolutionChannel = chIndex;
    ++imgIndex;
  }
  for (size_t channel = 0; channel != _channelsInDeconvolution; ++channel) {
    size_t sqIndex = channelToSqIndex(channel);
    const ImagingTable::Group& sqGroup = _imagingTable.SquaredGroups()[sqIndex];
    for (const ImagingTable::EntryPtr& entry : sqGroup) {
      size_t imageIndex = _entryIndexToImageIndex.find(entry->index)->second;
      _imageIndexToPSFIndex[imageIndex] = channel;
    }
  }
}

void ImageSet::LoadAndAverage(const CachedImageSet& imageSet) {
  for (size_t i = 0; i != _images.size(); ++i) _images[i] = 0.0;

  ImageF scratch(_width, _height);

  aocommon::UVector<size_t> averagedWeights(_images.size(), 0.0);
  size_t imgIndex = 0;
  for (size_t sqIndex = 0; sqIndex != _imagingTable.SquaredGroups().size();
       ++sqIndex) {
    // The next loop iterates over the polarizations. The logic in the next loop
    // makes sure that images of the same polarizations and that belong to the
    // same deconvolution channel are averaged together.
    const size_t imgIndexForChannel = imgIndex;
    const ImagingTable::Group& sqGroup = _imagingTable.SquaredGroups()[sqIndex];
    for (const ImagingTable::EntryPtr& e : sqGroup) {
      for (size_t i = 0; i != e->imageCount; ++i) {
        imageSet.Load(scratch.data(), e->polarization, e->outputChannelIndex,
                      i == 1);
        _images[imgIndex].AddWithFactor(scratch, e->imageWeight);
        averagedWeights[imgIndex] += e->imageWeight;
        ++imgIndex;
      }
    }
    const size_t thisChannelIndex = (sqIndex * _channelsInDeconvolution) /
                                    _imagingTable.SquaredGroups().size();
    const size_t nextChannelIndex = ((sqIndex + 1) * _channelsInDeconvolution) /
                                    _imagingTable.SquaredGroups().size();
    // If the next loaded image belongs to the same deconvolution channel as the
    // previously loaded, they need to be averaged together.
    if (thisChannelIndex == nextChannelIndex) imgIndex = imgIndexForChannel;
  }

  for (size_t i = 0; i != _images.size(); ++i)
    _images[i] *= 1.0 / double(averagedWeights[i]);
}

void ImageSet::LoadAndAveragePSFs(
    const CachedImageSet& psfSet,
    std::vector<aocommon::UVector<float>>& psfImages,
    aocommon::PolarizationEnum psfPolarization) {
  for (size_t chIndex = 0; chIndex != _channelsInDeconvolution; ++chIndex)
    psfImages[chIndex].assign(_width * _height, 0.0);

  ImageF scratch(_width, _height);

  aocommon::UVector<size_t> averagedWeights(_channelsInDeconvolution, 0.0);
  for (size_t sqIndex = 0; sqIndex != _imagingTable.SquaredGroups().size();
       ++sqIndex) {
    size_t chIndex = (sqIndex * _channelsInDeconvolution) /
                     _imagingTable.SquaredGroups().size();
    const ImagingTable::Group& sqGroup = _imagingTable.SquaredGroups()[sqIndex];
    const ImagingTableEntry& entry = *sqGroup.front();
    const double inputChannelWeight = entry.imageWeight;
    psfSet.Load(scratch.data(), psfPolarization, entry.outputChannelIndex, 0);
     for (size_t i = 0; i != _width * _height; ++i) {
      psfImages[chIndex][i] += scratch[i] * inputChannelWeight;
    }
    averagedWeights[chIndex] += inputChannelWeight;
  }

  for (size_t chIndex = 0; chIndex != ChannelsInDeconvolution(); ++chIndex) {
    for (size_t i = 0; i != _width * _height; ++i) {
      psfImages[chIndex][i] *= 1.0 / double(averagedWeights[chIndex]);
    }
  }
}

void ImageSet::InterpolateAndStore(CachedImageSet& imageSet,
                                   const SpectralFitter& fitter) {
  if (_channelsInDeconvolution == _imagingTable.SquaredGroups().size()) {
    directStore(imageSet);
  } else {
    Logger::Info << "Interpolating from " << _channelsInDeconvolution << " to "
                 << _imagingTable.SquaredGroups().size() << " channels...\n";

    // TODO should use spectralimagefitter to do the interpolation of images;
    // here we should just unpack the data structure

    // The following loop will make an 'image' with at each pixel
    // the terms of the fit. By doing this first, it is not necessary
    // to have all channel images in memory at the same time.
    // TODO: this assumes that polarizations are not joined!
    size_t nTerms = fitter.NTerms();
    aocommon::UVector<float> termsImage(_width * _height * nTerms);
    aocommon::StaticFor<size_t> loop(_settings.threadCount);
    loop.Run(0, _width * _height, [&](size_t pxStart, size_t pxEnd) {
      aocommon::UVector<float> spectralPixel(_channelsInDeconvolution);
      aocommon::UVector<float> termsPixel(nTerms);
      for (size_t px = pxStart; px != pxEnd; ++px) {
        bool isZero = true;
        for (size_t s = 0; s != _images.size(); ++s) {
          float value = _images[s][px];
          spectralPixel[s] = value;
          isZero = isZero && (value == 0.0);
        }
        float* termsPtr = &termsImage[px * nTerms];
        // Skip fitting if it is zero; most of model images will be zero, so
        // this can save a lot of time.
        if (isZero) {
          for (float* p = termsPtr; p != termsPtr + nTerms; ++p) *p = 0.0;
        } else {
          fitter.Fit(termsPixel, spectralPixel.data());
          for (size_t i = 0; i != nTerms; ++i) termsPtr[i] = termsPixel[i];
        }
      }
    });

    // Now that we know the fit for each pixel, evaluate the function for each
    // pixel of each output channel.
    Image scratch(_width, _height);
    size_t imgIndex = 0;
    for (const ImagingTableEntry& e : _imagingTable) {
      double freq = e.CentralFrequency();
      loop.Run(0, _width * _height, [&](size_t pxStart, size_t pxEnd) {
        aocommon::UVector<float> termsPixel(nTerms);
        for (size_t px = pxStart; px != pxEnd; ++px) {
          const float* termsPtr = &termsImage[px * nTerms];
          for (size_t i = 0; i != nTerms; ++i) termsPixel[i] = termsPtr[i];
          scratch[px] = fitter.Evaluate(termsPixel, freq);
        }
      });

      imageSet.Store(scratch.data(), e.polarization, e.outputChannelIndex,
                     false);
      ++imgIndex;
    }
  }
}

void ImageSet::AssignAndStore(CachedImageSet& imageSet) {
  if (_channelsInDeconvolution == _imagingTable.SquaredGroups().size()) {
    directStore(imageSet);
  } else {
    Logger::Info << "Assigning from " << _channelsInDeconvolution << " to "
                 << _imagingTable.SquaredGroups().size() << " channels...\n";
    size_t imgIndex = 0;
    size_t sqIndex = 0;
    for (const ImagingTable::Group& sqGroup : _imagingTable.SquaredGroups()) {
      size_t imgIndexForChannel = imgIndex;
      for (const ImagingTable::EntryPtr& e : sqGroup) {
        for (size_t i = 0; i != e->imageCount; ++i) {
          imageSet.Store(_images[imgIndex].data(), e->polarization,
                         e->outputChannelIndex, i == 1);
          ++imgIndex;
        }
      }
      size_t thisChannelIndex = (sqIndex * _channelsInDeconvolution) /
                                _imagingTable.SquaredGroups().size();
      size_t nextChannelIndex = ((sqIndex + 1) * _channelsInDeconvolution) /
                                _imagingTable.SquaredGroups().size();
      if (thisChannelIndex == nextChannelIndex) imgIndex = imgIndexForChannel;
      ++sqIndex;
    }
  }
}

void ImageSet::directStore(CachedImageSet& imageSet) {
  size_t imgIndex = 0;
  for (const ImagingTableEntry& e : _imagingTable) {
    for (size_t i = 0; i != e.imageCount; ++i) {
      imageSet.Store(_images[imgIndex].data(), e.polarization,
                     e.outputChannelIndex, i == 1);
      ++imgIndex;
    }
  }
}

void ImageSet::getSquareIntegratedWithNormalChannels(ImageF& dest,
                                                     ImageF& scratch) const {
  // In case only one frequency channel is used, we do not have to use
  // 'scratch', which saves copying and normalizing the data.
  if (_channelsInDeconvolution == 1) {
    const ImagingTable::Group& sqGroup = _imagingTable.SquaredGroups().front();
    if (sqGroup.size() == 1) {
      const ImagingTable::EntryPtr& entry = sqGroup.front();
      dest = entryToImage(entry);
    } else {
      const bool useAllPolarizations = _linkedPolarizations.empty();
      bool isFirst = true;
      for (const ImagingTable::EntryPtr& entry : sqGroup) {
        if (useAllPolarizations ||
            _linkedPolarizations.count(entry->polarization) != 0) {
          if (isFirst) {
            dest = entryToImage(entry);
            dest.Square();
            isFirst = false;
          } else {
            dest.AddSquared(entryToImage(entry));
          }
        }
      }
      squareRootMultiply(dest, std::sqrt(_polarizationNormalizationFactor));
    }
  } else {
    double weightSum = 0.0;
    for (size_t chIndex = 0; chIndex != _channelsInDeconvolution; ++chIndex) {
      size_t sqIndex = channelToSqIndex(chIndex);
      const ImagingTable::Group& sqGroup =
          _imagingTable.SquaredGroups()[sqIndex];
      const double groupWeight = _weights[chIndex];
      weightSum += groupWeight;
      if (sqGroup.size() == 1) {
        const ImagingTable::EntryPtr& entry = sqGroup.front();
        scratch = entryToImage(entry);
      } else {
        const bool useAllPolarizations = _linkedPolarizations.empty();
        bool isFirst = true;
        for (const ImagingTable::EntryPtr& entry : sqGroup) {
          if (useAllPolarizations ||
              _linkedPolarizations.count(entry->polarization) != 0) {
            if (isFirst) {
              scratch = entryToImage(entry);
              scratch.Square();
              isFirst = false;
            } else {
              scratch.AddSquared(entryToImage(entry));
            }
          }
        }
        scratch.Sqrt();
      }

      if (chIndex == 0)
        assignMultiply(dest, scratch, groupWeight);
      else
        dest.AddWithFactor(scratch, groupWeight);
    }
    if (_channelsInDeconvolution > 0)
      dest *= std::sqrt(_polarizationNormalizationFactor) / weightSum;
    else
      dest = 0.0;
  }
}

void ImageSet::getSquareIntegratedWithSquaredChannels(ImageF& dest) const {
  bool isFirst = true;
  const bool useAllPolarizations = _linkedPolarizations.empty();
  for (size_t channel = 0; channel != _channelsInDeconvolution; ++channel) {
    // TODO this should be weighted
    size_t sqIndex = channelToSqIndex(channel);
    const ImagingTable::Group& sqGroup = _imagingTable.SquaredGroups()[sqIndex];
    for (const ImagingTable::EntryPtr& entry : sqGroup) {
      if (useAllPolarizations ||
          _linkedPolarizations.count(entry->polarization) != 0) {
        if (isFirst) {
          dest = entryToImage(entry);
          dest.Square();
          isFirst = false;
        } else {
          dest.AddSquared(entryToImage(entry));
        }
      }
    }
  }
  double factor = _channelsInDeconvolution > 0
                      ? sqrt(_polarizationNormalizationFactor) /
                            double(_channelsInDeconvolution)
                      : 0.0;
  squareRootMultiply(dest, factor);
}

void ImageSet::getLinearIntegratedWithNormalChannels(ImageF& dest) const {
  const bool useAllPolarizations = _linkedPolarizations.empty();
  if (_channelsInDeconvolution == 1 &&
      _imagingTable.SquaredGroups().front().size() == 1) {
    const ImagingTable::Group& sqGroup = _imagingTable.SquaredGroups().front();
    const ImagingTable::EntryPtr& entry = sqGroup.front();
    dest = entryToImage(entry);
  } else {
    bool isFirst = true;
    double weightSum = 0.0;
    for (size_t channel = 0; channel != _channelsInDeconvolution; ++channel) {
      size_t sqIndex = channelToSqIndex(channel);
      const ImagingTable::Group& sqGroup =
          _imagingTable.SquaredGroups()[sqIndex];
      const double groupWeight = _weights[channel];
      weightSum += groupWeight;
      for (const ImagingTable::EntryPtr& entry : sqGroup) {
        if (useAllPolarizations ||
            _linkedPolarizations.count(entry->polarization) != 0) {
          if (isFirst) {
            assignMultiply(dest, entryToImage(entry), groupWeight);
            isFirst = false;
          } else {
            dest.AddWithFactor(entryToImage(entry), groupWeight);
          }
        }
      }
    }
    if (weightSum > 0.0)
      dest *= _polarizationNormalizationFactor / weightSum;
    else
      dest = 0.0;
  }
}

void ImageSet::CalculateDeconvolutionFrequencies(
    const ImagingTable& groupTable, aocommon::UVector<double>& frequencies,
    aocommon::UVector<float>& weights, size_t nDeconvolutionChannels) {
  size_t nInputChannels = groupTable.SquaredGroups().size();
  if (nDeconvolutionChannels == 0) nDeconvolutionChannels = nInputChannels;
  frequencies.assign(nDeconvolutionChannels, 0.0);
  weights.assign(nDeconvolutionChannels, 0.0);
  aocommon::UVector<double> weightSums(nDeconvolutionChannels, 0);
  for (size_t i = 0; i != nInputChannels; ++i) {
    const ImagingTableEntry& entry = *groupTable.SquaredGroups()[i].front();
    double freq = entry.CentralFrequency(), weight = entry.imageWeight;
    size_t deconvolutionChannel = i * nDeconvolutionChannels / nInputChannels;
    frequencies[deconvolutionChannel] += freq * weight;
    weights[deconvolutionChannel] += weight;
  }
  for (size_t i = 0; i != nDeconvolutionChannels; ++i)
    frequencies[i] /= weights[i];
}

void ImageSet::GetIntegratedPSF(ImageF& dest,
                                const aocommon::UVector<const float*>& psfs) {
  // TODO should use weighting!
  std::copy_n(psfs[0], _width * _height, dest.data());
  for (size_t img = 1; img != PSFCount(); ++img) {
    for (size_t i = 0; i != _width * _height; ++i) dest[i] += psfs[img][i];
  }
  if (PSFCount() != 1) {
    for (size_t i = 0; i != _width * _height; ++i)
      dest[i] *= 1.0 / float(PSFCount());
  }
}
