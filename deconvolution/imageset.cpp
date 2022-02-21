#include "imageset.h"

#include "spectralfitter.h"

#include "../math/nlplfitter.h"

#include "../structures/primarybeam.h"
#include "../structures/primarybeamimageset.h"

#include <aocommon/logger.h>
#include <aocommon/staticfor.h>

#include <cassert>

using aocommon::Image;
using aocommon::Logger;

ImageSet::ImageSet(const DeconvolutionTable& table,
                   const class Settings& settings)
    : _images(),
      _width(0),
      _height(0),
      _squareJoinedChannels(settings.squaredJoins),
      _deconvolutionTable(table),
      _imageIndexToPSFIndex(),
      _linkedPolarizations(settings.linkedPolarizations),
      _settings(settings) {
  size_t nPol = table.OriginalGroups().front().size();
  size_t nImages = nPol * ChannelsInDeconvolution();
  _images.resize(nImages);
  _imageIndexToPSFIndex.resize(nImages);

  initializePolFactor();
  initializeIndices();
  aocommon::UVector<double> frequencies;
  CalculateDeconvolutionFrequencies(table, frequencies, _weights,
                                    ChannelsInDeconvolution());
}

ImageSet::ImageSet(const DeconvolutionTable& table,
                   const class Settings& settings, size_t width, size_t height)
    : ImageSet(table, settings) {
  _width = width;
  _height = height;
  allocateImages();
}

void ImageSet::initializeIndices() {
  _entryIndexToImageIndex.reserve(_deconvolutionTable.Size());
  size_t imgIndex = 0;
  for (const std::vector<int>& group :
       _deconvolutionTable.DeconvolutionGroups()) {
    const size_t deconvolutionChannelStartIndex = imgIndex;
    for (const int originalIndex : group) {
      imgIndex = deconvolutionChannelStartIndex;

      for (const DeconvolutionTableEntry* entry :
           _deconvolutionTable.OriginalGroups()[originalIndex]) {
        assert(entry->index == _entryIndexToImageIndex.size());
        _entryIndexToImageIndex.push_back(imgIndex);

        ++imgIndex;
      }
    }
  }

  for (size_t chIndex = 0; chIndex != ChannelsInDeconvolution(); ++chIndex) {
    const DeconvolutionTable::Group& originalGroup =
        _deconvolutionTable.FirstOriginalGroup(chIndex);
    for (const DeconvolutionTableEntry* entry : originalGroup) {
      const size_t imageIndex = _entryIndexToImageIndex[entry->index];
      _imageIndexToPSFIndex[imageIndex] = chIndex;
    }
  }
}

void ImageSet::SetImages(ImageSet&& source) {
  for (size_t imageIndex = 0; imageIndex != _images.size(); ++imageIndex) {
    _images[imageIndex] = std::move(source._images[imageIndex]);
  }
  _width = source._width;
  _height = source._height;
  source._width = 0;
  source._height = 0;
}

void ImageSet::LoadAndAverage(bool use_residual_image) {
  for (Image& image : _images) {
    image = 0.0;
  }

  Image scratch(_width, _height);

  aocommon::UVector<double> averagedWeights(_images.size(), 0.0);

  size_t imgIndex = 0;
  for (const std::vector<int>& group :
       _deconvolutionTable.DeconvolutionGroups()) {
    const size_t deconvolutionChannelStartIndex = imgIndex;
    for (const int originalIndex : group) {
      // The next loop iterates over the polarizations. The logic in the next
      // loop makes sure that images of the same polarizations and that belong
      // to the same deconvolution channel are averaged together.
      imgIndex = deconvolutionChannelStartIndex;
      for (const DeconvolutionTableEntry* entry_ptr :
           _deconvolutionTable.OriginalGroups()[originalIndex]) {
        if (use_residual_image) {
          entry_ptr->residual_accessor->Load(scratch);
        } else {
          entry_ptr->model_accessor->Load(scratch);
        }
        _images[imgIndex].AddWithFactor(scratch, entry_ptr->image_weight);
        averagedWeights[imgIndex] += entry_ptr->image_weight;
        ++imgIndex;
      }
    }
  }

  for (size_t i = 0; i != _images.size(); ++i) {
    _images[i] *= 1.0 / averagedWeights[i];
  }
}

void ImageSet::LoadAndAveragePSFs(
    std::vector<aocommon::UVector<float>>& psfImages,
    aocommon::PolarizationEnum psfPolarization) {
  for (size_t chIndex = 0; chIndex != ChannelsInDeconvolution(); ++chIndex)
    psfImages[chIndex].assign(_width * _height, 0.0);

  Image scratch(_width, _height);

  aocommon::UVector<double> averagedWeights(ChannelsInDeconvolution(), 0.0);
  for (size_t groupIndex = 0;
       groupIndex != _deconvolutionTable.OriginalGroups().size();
       ++groupIndex) {
    size_t chIndex = (groupIndex * ChannelsInDeconvolution()) /
                     _deconvolutionTable.OriginalGroups().size();
    const DeconvolutionTable::Group& channelGroup =
        _deconvolutionTable.OriginalGroups()[groupIndex];
    const DeconvolutionTableEntry& entry = *channelGroup.front();
    const double inputChannelWeight = entry.image_weight;
    entry.psf_accessor->Load(scratch);
    for (size_t i = 0; i != _width * _height; ++i) {
      psfImages[chIndex][i] += scratch[i] * inputChannelWeight;
    }
    averagedWeights[chIndex] += inputChannelWeight;
  }

  for (size_t chIndex = 0; chIndex != ChannelsInDeconvolution(); ++chIndex) {
    const double factor =
        averagedWeights[chIndex] == 0.0 ? 0.0 : 1.0 / averagedWeights[chIndex];
    for (size_t i = 0; i != _width * _height; ++i) {
      psfImages[chIndex][i] *= factor;
    }
  }
}

void ImageSet::InterpolateAndStoreModel(const SpectralFitter& fitter) {
  if (ChannelsInDeconvolution() ==
      _deconvolutionTable.OriginalGroups().size()) {
    size_t imgIndex = 0;
    for (const DeconvolutionTableEntry& e : _deconvolutionTable) {
      e.model_accessor->Store(_images[imgIndex]);
      ++imgIndex;
    }
  } else {
    Logger::Info << "Interpolating from " << ChannelsInDeconvolution() << " to "
                 << _deconvolutionTable.OriginalGroups().size()
                 << " channels...\n";

    // TODO should use spectralimagefitter to do the interpolation of images;
    // here we should just unpack the data structure

    // The following loop will make an 'image' with at each pixel
    // the terms of the fit. By doing this first, it is not necessary
    // to have all channel images in memory at the same time.
    // TODO: this assumes that polarizations are not joined!
    size_t nTerms = fitter.NTerms();
    aocommon::UVector<float> termsImage(_width * _height * nTerms);
    aocommon::StaticFor<size_t> loop(_settings.threadCount);
    loop.Run(0, _height, [&](size_t yStart, size_t yEnd) {
      aocommon::UVector<float> spectralPixel(ChannelsInDeconvolution());
      aocommon::UVector<float> termsPixel(nTerms);
      for (size_t y = yStart; y != yEnd; ++y) {
        size_t px = y * _width;
        for (size_t x = 0; x != _width; ++x) {
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
            fitter.Fit(termsPixel, spectralPixel.data(), x, y);
            for (size_t i = 0; i != nTerms; ++i) termsPtr[i] = termsPixel[i];
          }
          ++px;
        }
      }
    });

    // Now that we know the fit for each pixel, evaluate the function for each
    // pixel of each output channel.
    Image scratch(_width, _height);
    for (const DeconvolutionTableEntry& e : _deconvolutionTable) {
      double freq = e.CentralFrequency();
      loop.Run(0, _width * _height, [&](size_t pxStart, size_t pxEnd) {
        aocommon::UVector<float> termsPixel(nTerms);
        for (size_t px = pxStart; px != pxEnd; ++px) {
          const float* termsPtr = &termsImage[px * nTerms];
          for (size_t i = 0; i != nTerms; ++i) termsPixel[i] = termsPtr[i];
          scratch[px] = fitter.Evaluate(termsPixel, freq);
        }
      });

      e.model_accessor->Store(scratch);
    }
  }
}

void ImageSet::AssignAndStoreResidual() {
  if (ChannelsInDeconvolution() ==
      _deconvolutionTable.OriginalGroups().size()) {
    size_t imgIndex = 0;
    for (const DeconvolutionTableEntry& e : _deconvolutionTable) {
      e.residual_accessor->Store(_images[imgIndex]);
      ++imgIndex;
    }
  } else {
    Logger::Info << "Assigning from " << ChannelsInDeconvolution() << " to "
                 << _deconvolutionTable.OriginalGroups().size()
                 << " channels...\n";
    size_t imgIndex = 0;
    size_t groupIndex = 0;
    for (const DeconvolutionTable::Group& channelGroup :
         _deconvolutionTable.OriginalGroups()) {
      size_t imgIndexForChannel = imgIndex;
      for (const DeconvolutionTableEntry* e : channelGroup) {
        e->residual_accessor->Store(_images[imgIndex]);
        ++imgIndex;
      }
      size_t thisChannelIndex = (groupIndex * ChannelsInDeconvolution()) /
                                _deconvolutionTable.OriginalGroups().size();
      size_t nextChannelIndex = ((groupIndex + 1) * ChannelsInDeconvolution()) /
                                _deconvolutionTable.OriginalGroups().size();
      if (thisChannelIndex == nextChannelIndex) imgIndex = imgIndexForChannel;
      ++groupIndex;
    }
  }
}

void ImageSet::getSquareIntegratedWithNormalChannels(Image& dest,
                                                     Image& scratch) const {
  // In case only one frequency channel is used, we do not have to use
  // 'scratch', which saves copying and normalizing the data.
  if (ChannelsInDeconvolution() == 1) {
    const DeconvolutionTable::Group& originalGroup =
        _deconvolutionTable.OriginalGroups().front();
    if (originalGroup.size() == 1) {
      const DeconvolutionTableEntry& entry = *originalGroup.front();
      dest = entryToImage(entry);
    } else {
      const bool useAllPolarizations = _linkedPolarizations.empty();
      bool isFirst = true;
      for (const DeconvolutionTableEntry* entry_ptr : originalGroup) {
        if (useAllPolarizations ||
            _linkedPolarizations.count(entry_ptr->polarization) != 0) {
          if (isFirst) {
            dest = entryToImage(*entry_ptr);
            dest.Square();
            isFirst = false;
          } else {
            dest.AddSquared(entryToImage(*entry_ptr));
          }
        }
      }
      squareRootMultiply(dest, std::sqrt(_polarizationNormalizationFactor));
    }
  } else {
    double weightSum = 0.0;
    bool isFirstChannel = true;

    for (size_t chIndex = 0; chIndex != ChannelsInDeconvolution(); ++chIndex) {
      const DeconvolutionTable::Group& originalGroup =
          _deconvolutionTable.FirstOriginalGroup(chIndex);
      const double groupWeight = _weights[chIndex];
      // if the groupWeight is zero, the image might contain NaNs, so we
      // shouldn't add it to the total in that case.
      if (groupWeight != 0.0) {
        weightSum += groupWeight;
        if (originalGroup.size() == 1) {
          const DeconvolutionTableEntry& entry = *originalGroup.front();
          scratch = entryToImage(entry);
        } else {
          const bool useAllPolarizations = _linkedPolarizations.empty();
          bool isFirstPolarization = true;
          for (const DeconvolutionTableEntry* entry_ptr : originalGroup) {
            if (useAllPolarizations ||
                _linkedPolarizations.count(entry_ptr->polarization) != 0) {
              if (isFirstPolarization) {
                scratch = entryToImage(*entry_ptr);
                scratch.Square();
                isFirstPolarization = false;
              } else {
                scratch.AddSquared(entryToImage(*entry_ptr));
              }
            }
          }
          scratch.Sqrt();
        }
      }

      if (isFirstChannel) {
        assignMultiply(dest, scratch, groupWeight);
        isFirstChannel = false;
      } else {
        dest.AddWithFactor(scratch, groupWeight);
      }
    }
    dest *= std::sqrt(_polarizationNormalizationFactor) / weightSum;
  }
}

void ImageSet::getSquareIntegratedWithSquaredChannels(Image& dest) const {
  bool isFirst = true;
  const bool useAllPolarizations = _linkedPolarizations.empty();
  double weightSum = 0.0;

  for (size_t chIndex = 0; chIndex != ChannelsInDeconvolution(); ++chIndex) {
    const double groupWeight = _weights[chIndex];
    weightSum += groupWeight;
    if (groupWeight != 0.0) {
      const DeconvolutionTable::Group& originalGroup =
          _deconvolutionTable.FirstOriginalGroup(chIndex);
      for (const DeconvolutionTableEntry* entry_ptr : originalGroup) {
        if (useAllPolarizations ||
            _linkedPolarizations.count(entry_ptr->polarization) != 0) {
          if (isFirst) {
            dest = entryToImage(*entry_ptr);
            dest.SquareWithFactor(groupWeight);
            isFirst = false;
          } else {
            dest.AddSquared(entryToImage(*entry_ptr), groupWeight);
          }
        }
      }
    }
  }
  double factor = weightSum > 0.0
                      ? std::sqrt(_polarizationNormalizationFactor) / weightSum
                      : 0.0;
  squareRootMultiply(dest, factor);
}

void ImageSet::getLinearIntegratedWithNormalChannels(Image& dest) const {
  const bool useAllPolarizations = _linkedPolarizations.empty();
  if (_deconvolutionTable.DeconvolutionGroups().size() == 1 &&
      _deconvolutionTable.OriginalGroups().front().size() == 1) {
    const DeconvolutionTable::Group& originalGroup =
        _deconvolutionTable.OriginalGroups().front();
    const DeconvolutionTableEntry& entry = *originalGroup.front();
    dest = entryToImage(entry);
  } else {
    bool isFirst = true;
    double weightSum = 0.0;
    for (size_t chIndex = 0; chIndex != ChannelsInDeconvolution(); ++chIndex) {
      const DeconvolutionTable::Group& originalGroup =
          _deconvolutionTable.FirstOriginalGroup(chIndex);
      const double groupWeight = _weights[chIndex];
      // if the groupWeight is zero, the image might contain NaNs, so we
      // shouldn't add it to the total in that case.
      if (groupWeight != 0.0) {
        weightSum += groupWeight;
        for (const DeconvolutionTableEntry* entry_ptr : originalGroup) {
          if (useAllPolarizations ||
              _linkedPolarizations.count(entry_ptr->polarization) != 0) {
            if (isFirst) {
              assignMultiply(dest, entryToImage(*entry_ptr), groupWeight);
              isFirst = false;
            } else {
              dest.AddWithFactor(entryToImage(*entry_ptr), groupWeight);
            }
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
    const DeconvolutionTable& groupTable,
    aocommon::UVector<double>& frequencies, aocommon::UVector<float>& weights,
    size_t nDeconvolutionChannels) {
  const size_t nInputChannels = groupTable.OriginalGroups().size();
  if (nDeconvolutionChannels == 0) nDeconvolutionChannels = nInputChannels;
  frequencies.assign(nDeconvolutionChannels, 0.0);
  weights.assign(nDeconvolutionChannels, 0.0);
  std::vector<double> unweightedFrequencies(nDeconvolutionChannels, 0.0);
  std::vector<size_t> counts(nDeconvolutionChannels, 0);
  for (size_t i = 0; i != nInputChannels; ++i) {
    const DeconvolutionTableEntry& entry =
        *groupTable.OriginalGroups()[i].front();
    const double freq = entry.CentralFrequency();
    const double weight = entry.image_weight;
    const size_t deconvolutionChannel =
        i * nDeconvolutionChannels / nInputChannels;

    frequencies[deconvolutionChannel] += freq * weight;
    weights[deconvolutionChannel] += weight;

    unweightedFrequencies[deconvolutionChannel] += freq;
    ++counts[deconvolutionChannel];
  }
  for (size_t i = 0; i != nDeconvolutionChannels; ++i) {
    // Even when there is no data for a given frequency and the weight
    // is zero, it is still desirable to have a proper value for the frequency
    // (e.g. for extrapolating flux).
    if (weights[i] > 0.0)
      frequencies[i] /= weights[i];
    else
      frequencies[i] = unweightedFrequencies[i] / counts[i];
  }
}

void ImageSet::GetIntegratedPSF(Image& dest,
                                const aocommon::UVector<const float*>& psfs) {
  if (ChannelsInDeconvolution() == 1)
    std::copy_n(psfs[0], _width * _height, dest.Data());
  else {
    bool isFirst = true;
    double weightSum = 0.0;
    for (size_t channel = 0; channel != ChannelsInDeconvolution(); ++channel) {
      const double groupWeight = _weights[channel];
      // if the groupWeight is zero, the image might contain NaNs, so we
      // shouldn't add it to the total in that case.
      if (groupWeight != 0.0) {
        weightSum += groupWeight;
        if (isFirst) {
          for (size_t i = 0; i != _width * _height; ++i)
            dest[i] = psfs[channel][i] * groupWeight;
          isFirst = false;
        } else {
          for (size_t i = 0; i != _width * _height; ++i)
            dest[i] += psfs[channel][i] * groupWeight;
        }
      }
    }
    const double factor = weightSum == 0.0 ? 0.0 : 1.0 / weightSum;
    for (size_t i = 0; i != _width * _height; ++i) dest[i] *= factor;
  }
}
