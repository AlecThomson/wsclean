#include "averagebeam.h"

#include "../io/cachedimageset.h"

#include <aocommon/io/serialistream.h>
#include <aocommon/io/serialostream.h>

void AverageBeam::Serialize(aocommon::SerialOStream& stream) const {
  if (scalar_beam_) {
    stream.Bool(true);
    stream.Vector(*scalar_beam_).UInt64(scalar_width_).UInt64(scalar_height_);
  } else {
    stream.Bool(false);
  }

  if (matrix_inverse_beam_) {
    stream.Bool(true);
    stream.Vector(*matrix_inverse_beam_)
        .UInt64(matrix_width_)
        .UInt64(matrix_height_);
  } else {
    stream.Bool(false);
  }
}

void AverageBeam::Unserialize(aocommon::SerialIStream& stream) {
  const bool hasScalar = stream.Bool();
  if (hasScalar) {
    scalar_beam_.reset(new std::vector<float>());
    stream.Vector(*scalar_beam_).UInt64(scalar_width_).UInt64(scalar_height_);
  } else
    scalar_beam_.reset();

  const bool hasMatrixInverse = stream.Bool();
  if (hasMatrixInverse) {
    matrix_inverse_beam_.reset(new std::vector<std::complex<float>>());
    stream.Vector(*matrix_inverse_beam_)
        .UInt64(matrix_width_)
        .UInt64(matrix_height_);
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

    // Scalar beam
    result->scalar_width_ = scalar_cache.Writer().Width();
    result->scalar_height_ = scalar_cache.Writer().Height();
    const size_t n_scalar_pixels =
        result->scalar_width_ * result->scalar_height_;
    result->scalar_beam_ = std::make_shared<std::vector<float>>(n_scalar_pixels);
    scalar_cache.Load(result->scalar_beam_->data(),
                      aocommon::PolarizationEnum::StokesI, frequency_index,
                      false);

    // Inverse matrix beam
    // It is stored as one real and one imaginary fits image. The 16 elements
    // are consecutively stored along the X axis, so its width is 4x its height.
    // This makes these images not easy to interpret, but it avoids having to
    // reshuffle the data, and they are only for temporary storage, not for
    // interpretation.
    constexpr size_t kNPolarizations = 4;
    assert(!matrix_cache.Empty());
    result->matrix_width_ =
        matrix_cache.Writer().Width() / (kNPolarizations * kNPolarizations);
    result->matrix_height_ = matrix_cache.Writer().Height();
    const size_t n_matrix_elements =
        matrix_cache.Writer().Width() * matrix_cache.Writer().Height();
    result->matrix_inverse_beam_ =
      std::make_shared<std::vector<std::complex<float>>>(n_matrix_elements);
    aocommon::UVector<float> real_image(n_matrix_elements);
    aocommon::UVector<float> imaginary_image(n_matrix_elements);
    matrix_cache.Load(real_image.data(), aocommon::PolarizationEnum::StokesI,
                      frequency_index, false);
    matrix_cache.Load(imaginary_image.data(),
                      aocommon::PolarizationEnum::StokesI, frequency_index,
                      true);
    for (size_t i = 0; i != n_matrix_elements; ++i) {
      (*result->matrix_inverse_beam_)[i] =
          std::complex<float>(real_image[i], imaginary_image[i]);
    }
  }
  return result;
}

void AverageBeam::Store(CachedImageSet& scalar_cache,
                        CachedImageSet& matrix_cache,
                        size_t frequency_index) const {
  if (!Empty()) {
    // Scalar beam
    scalar_cache.Writer().SetImageDimensions(scalar_width_, scalar_height_);
    scalar_cache.Store(scalar_beam_->data(),
                       aocommon::PolarizationEnum::StokesI, frequency_index,
                       false);

    // Matrix beam
    constexpr size_t kNPolarizations = 4;
    const size_t n_pixels =
        matrix_width_ * matrix_height_ * kNPolarizations * kNPolarizations;
    matrix_cache.Writer().SetImageDimensions(
        matrix_width_ * kNPolarizations * kNPolarizations, matrix_height_);
    aocommon::UVector<float> real_image(n_pixels);
    aocommon::UVector<float> imaginary_image(n_pixels);
    for (size_t i = 0; i != n_pixels; ++i) {
      real_image[i] = (*matrix_inverse_beam_)[i].real();
      imaginary_image[i] = (*matrix_inverse_beam_)[i].imag();
    }
    matrix_cache.Store(real_image.data(), aocommon::PolarizationEnum::StokesI,
                       frequency_index, false);
    matrix_cache.Store(imaginary_image.data(),
                       aocommon::PolarizationEnum::StokesI, frequency_index,
                       true);
  }
}
