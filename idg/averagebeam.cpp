#include "averagebeam.h"

#include "../io/cachedimageset.h"

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

void AverageBeam::Serialize(aocommon::SerialOStream& stream) const {
  if (_scalarBeam) {
    stream.Bool(true);
    stream.Vector(*_scalarBeam);
  } else {
    stream.Bool(false);
  }

  if (_matrixInverseBeam) {
    stream.Bool(true);
    stream.Vector(*_matrixInverseBeam);
  } else {
    stream.Bool(false);
  }
}

void AverageBeam::Unserialize(aocommon::SerialIStream& stream) {
  bool hasScalar = stream.Bool();
  if (hasScalar) {
    _scalarBeam.reset(new std::vector<float>());
    stream.Vector(*_scalarBeam);
  } else
    _scalarBeam.reset();

  bool hasMatrixInverse = stream.Bool();
  if (hasMatrixInverse) {
    _matrixInverseBeam.reset(new std::vector<std::complex<float>>());
    stream.Vector(*_matrixInverseBeam);
  } else
    _matrixInverseBeam.reset();
}

std::unique_ptr<AverageBeam> AverageBeam::Load(
    const CachedImageSet& scalar_cache, const CachedImageSet& matrix_cache,
    size_t frequency_index, size_t n_pixels) {
  std::unique_ptr<AverageBeam> result;
  if (!scalar_cache.Empty()) {
    result.reset(new AverageBeam());
    result->_scalarBeam->resize(n_pixels);
    scalar_cache.Load(result->_scalarBeam->data(),
                      aocommon::PolarizationEnum::StokesI, frequency_index,
                      false);
    constexpr size_t kNPolarizations = 4;
    constexpr aocommon::PolarizationEnum kPolarizations[kNPolarizations] = {
        aocommon::PolarizationEnum::XX,
        aocommon::PolarizationEnum::XY,
        aocommon::PolarizationEnum::YX,
        aocommon::PolarizationEnum::YY,
    };

    result->_matrixInverseBeam->resize(kNPolarizations * n_pixels);
    aocommon::UVector<float> real_image(n_pixels);
    aocommon::UVector<float> imaginary_image(n_pixels);
    for (size_t p = 0; p != kNPolarizations; ++p) {
      matrix_cache.Load(real_image.data(), kPolarizations[p], frequency_index,
                        false);
      matrix_cache.Load(imaginary_image.data(), kPolarizations[p],
                        frequency_index, false);
      for (size_t i = 0; i != real_image.size(); ++i) {
        (*result->_matrixInverseBeam)[i * kNPolarizations] =
            std::complex<float>(real_image[i], imaginary_image[i]);
      }
    }
  }
  return result;
}

void AverageBeam::Store(CachedImageSet& scalar_cache,
                        CachedImageSet& matrix_cache,
                        size_t frequency_index) const {
  scalar_cache.Store(_scalarBeam->data(), aocommon::PolarizationEnum::StokesI,
                     frequency_index, false);

  constexpr size_t kNPolarizations = 4;
  constexpr aocommon::PolarizationEnum kPolarizations[kNPolarizations] = {
      aocommon::PolarizationEnum::XX,
      aocommon::PolarizationEnum::XY,
      aocommon::PolarizationEnum::YX,
      aocommon::PolarizationEnum::YY,
  };
  aocommon::UVector<float> real_image(_matrixInverseBeam->size() /
                                      kNPolarizations);
  aocommon::UVector<float> imaginary_image(_matrixInverseBeam->size() /
                                           kNPolarizations);
  for (size_t p = 0; p != kNPolarizations; ++p) {
    for (size_t i = 0; i != real_image.size(); ++i) {
      real_image[i] = (*_matrixInverseBeam)[i * kNPolarizations].real();
      imaginary_image[i] = (*_matrixInverseBeam)[i * kNPolarizations].imag();
    }
    matrix_cache.Store(real_image.data(), kPolarizations[p], frequency_index,
                       false);
    matrix_cache.Store(imaginary_image.data(), kPolarizations[p],
                       frequency_index, false);
  }
}
