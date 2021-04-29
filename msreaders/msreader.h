#ifndef MSREADERS_MSREADER_H
#define MSREADERS_MSREADER_H

#include "../msproviders/msprovider.h"

/**
 * An MSReader provides the capability to read (meta)data
 * independent from the location of the reader position in the MSProvider
 */
class MSReader {
 public:
  MSReader(MSProvider* msProvider) : _msProvider(msProvider){};

  virtual ~MSReader(){};

  /**
   * This provides a unique, consecutive number that corresponds to
   * the current reading position. Note that this number does not have
   * to map directly to measurement set row indices, because unselected
   * data does not affect the RowId. @ref MakeIdToMSRowMapping()
   * can be used to convert this Id to a measurement row number.
   */
  virtual size_t RowId() const = 0;

  /**
   * Returns true as long as there is more data available for reading.
   */
  virtual bool CurrentRowAvailable() = 0;

  /**
   * Move the reading position to the next row.
   */
  virtual void NextInputRow() = 0;

  /**
   * @{
   * Read meta data from the current reading position.
   */
  virtual void ReadMeta(double& u, double& v, double& w,
                        size_t& dataDescId) = 0;

  virtual void ReadMeta(MSProvider::MetaData& metaData) = 0;
  /** @} */

  /**
   * Read visibility data from current reading position.
   */
  virtual void ReadData(std::complex<float>* buffer) = 0;

  /**
   * Read the model visibilities from the current reading position.
   */
  virtual void ReadModel(std::complex<float>* buffer) = 0;

  virtual void ReadWeights(float* buffer) = 0;

  virtual void ReadWeights(std::complex<float>* buffer) = 0;

  /**
   * Write imaging weights to the current READING position.
   * Note that despite this is a write operation, the reading position is
   * used nevertheless. This is because it is written while reading the meta
   * data inside WSClean, hence it would be inconvenient if the writing position
   * would be used.
   */
  virtual void WriteImagingWeights(const float* buffer) = 0;

 protected:
  MSProvider* _msProvider;

  static void copyData(std::complex<float>* dest, size_t startChannel,
                       size_t endChannel,
                       const std::vector<aocommon::PolarizationEnum>& polsIn,
                       const casacore::Array<std::complex<float>>& data,
                       aocommon::PolarizationEnum polOut);

  template <typename NumType>
  static bool isCFinite(const std::complex<NumType>& c) {
    return std::isfinite(c.real()) && std::isfinite(c.imag());
  }
};

#endif
