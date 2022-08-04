#ifndef CACHED_IMAGE_ACCESSOR_H_
#define CACHED_IMAGE_ACCESSOR_H_

#include <aocommon/dataaccessor.h>

#include "cachedimageset.h"

/**
 * @brief DataAccessor implementation that internally uses CachedImageSet.
 */
class CachedImageAccessor : public aocommon::DataAccessor<float> {
 public:
  /**
   * @brief Construct a new CachedImageAccessor object
   *
   * @param image_set The CachedImageSet for loading and storing the image.
   * @param polarization Image polarization.
   * @param frequency_index Frequency index of the image.
   * @param is_imaginary False: Image has real values. True: Image has imaginary
   * values.
   */
  CachedImageAccessor(CachedImageSet& image_set,
                      aocommon::PolarizationEnum polarization,
                      size_t frequency_index, bool is_imaginary)
      : image_set_(image_set),
        polarization_(polarization),
        frequency_index_(frequency_index),
        is_imaginary_(is_imaginary) {}

  std::size_t Size() const override {
    return image_set_.Writer().Width() * image_set_.Writer().Height();
  }

  void Load(float* data) const override {
    image_set_.Load(data, polarization_, frequency_index_, is_imaginary_);
  }

  void Store(const float* data) override {
    image_set_.Store(data, polarization_, frequency_index_, is_imaginary_);
  }

  /**
   * Get functions, mainly for testing purposes.
   * @{
   */
  const CachedImageSet& GetImageSet() const { return image_set_; }
  aocommon::PolarizationEnum GetPolarization() const { return polarization_; }
  size_t GetFrequencyIndex() const { return frequency_index_; }
  bool GetIsImaginary() const { return is_imaginary_; }
  /** @} */

 private:
  CachedImageSet& image_set_;
  aocommon::PolarizationEnum polarization_;
  size_t frequency_index_;
  bool is_imaginary_;
};

#endif