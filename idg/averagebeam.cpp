#include "averagebeam.h"

#include "../io/cachedimageset.h"

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

void AverageBeam::Serialize(aocommon::SerialOStream& stream) const {
  if (scalar_beam_) {
    stream.Bool(true);
    stream.Vector(*scalar_beam_);
  } else {
    stream.Bool(false);
  }

  if (matrix_inverse_beam_) {
    stream.Bool(true);
    stream.Vector(*matrix_inverse_beam_);
  } else {
    stream.Bool(false);
  }
}

void AverageBeam::Unserialize(aocommon::SerialIStream& stream) {
  bool hasScalar = stream.Bool();
  if (hasScalar) {
    scalar_beam_.reset(new std::vector<float>());
    stream.Vector(*scalar_beam_);
  } else
    scalar_beam_.reset();

  bool hasMatrixInverse = stream.Bool();
  if (hasMatrixInverse) {
    matrix_inverse_beam_.reset(new std::vector<std::complex<float>>());
    stream.Vector(*matrix_inverse_beam_);
  } else
    matrix_inverse_beam_.reset();
}

std::unique_ptr<AverageBeam> AverageBeam::Load(
    const CachedImageSet& scalar_cache, const CachedImageSet& matrix_cache,
    size_t frequency_index) {
  std::unique_ptr<AverageBeam> result;
  if (!scalar_cache.Empty()) {
    Logger::Debug << "Loading average beam from cache.\n";
    result.reset(new AverageBeam());
    const size_t n_scalar_pixels =
        scalar_cache.Writer().Width() * scalar_cache.Writer().Height();
    result->scalar_beam_.reset(new std::vector<float>(n_scalar_pixels));
    scalar_cache.Load(result->scalar_beam_->data(),
                      aocommon::PolarizationEnum::StokesI, frequency_index,
                      false);
    constexpr size_t kNPolarizations = 4;
    constexpr aocommon::PolarizationEnum kPolarizations[kNPolarizations] = {
        aocommon::PolarizationEnum::XX,
        aocommon::PolarizationEnum::XY,
        aocommon::PolarizationEnum::YX,
        aocommon::PolarizationEnum::YY,
    };

    assert(!matrix_cache.Empty());
    const size_t n_matrix_pixels =
        matrix_cache.Writer().Width() * matrix_cache.Writer().Height();
    result->matrix_inverse_beam_.reset(new std::vector<std::complex<float>>(
        kNPolarizations * n_matrix_pixels));
    aocommon::UVector<float> real_image(n_matrix_pixels);
    aocommon::UVector<float> imaginary_image(n_matrix_pixels);
    for (size_t p = 0; p != kNPolarizations; ++p) {
      matrix_cache.Load(real_image.data(), kPolarizations[p], frequency_index,
                        false);
      matrix_cache.Load(imaginary_image.data(), kPolarizations[p],
                        frequency_index, true);
      for (size_t i = 0; i != real_image.size(); ++i) {
        (*result->matrix_inverse_beam_)[i * kNPolarizations + p] =
            std::complex<float>(real_image[i], imaginary_image[i]);
      }
    }
  }
  return result;
}

void AverageBeam::Store(CachedImageSet& scalar_cache,
                        CachedImageSet& matrix_cache,
                        size_t frequency_index) const {
  if (!Empty()) {
    scalar_cache.Store(scalar_beam_->data(),
                       aocommon::PolarizationEnum::StokesI, frequency_index,
                       false);

    constexpr size_t kNPolarizations = 4;
    constexpr aocommon::PolarizationEnum kPolarizations[kNPolarizations] = {
        aocommon::PolarizationEnum::XX,
        aocommon::PolarizationEnum::XY,
        aocommon::PolarizationEnum::YX,
        aocommon::PolarizationEnum::YY,
    };
    matrix_cache.Writer().SetImageDimensions(matrix_width_, matrix_height_);
    aocommon::UVector<float> real_image(matrix_width_ * matrix_height_);
    aocommon::UVector<float> imaginary_image(matrix_width_ * matrix_height_);
    for (size_t p = 0; p != kNPolarizations; ++p) {
      for (size_t i = 0; i != real_image.size(); ++i) {
        real_image[i] = (*matrix_inverse_beam_)[i * kNPolarizations + p].real();
        imaginary_image[i] =
            (*matrix_inverse_beam_)[i * kNPolarizations + p].imag();
      }
      matrix_cache.Store(real_image.data(), kPolarizations[p], frequency_index,
                         false);
      matrix_cache.Store(imaginary_image.data(), kPolarizations[p],
                         frequency_index, true);
    }
  }
}
