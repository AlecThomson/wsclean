#ifndef WGRIDDING_GRIDDER_SIMPLE_H
#define WGRIDDING_GRIDDER_SIMPLE_H

#include <complex>
#include <cstddef>
#include <vector>

#include <aocommon/banddata.h>
#include "ducc0/wgridder/wgridder.h"

#include "../gridding/gainmode.h"

namespace wsclean {

class MsGridder;

class WGriddingGridderBase {
 public:
  virtual ~WGriddingGridderBase() = default;
  virtual size_t ConstantMemoryUsage() const = 0;
  virtual size_t PerVisibilityMemoryUsage() const = 0;
  virtual void InitializeInversion() = 0;
  virtual void AddInversionData(size_t n_rows, size_t n_chan, const double *uvw,
                                const double *freq,
                                const std::complex<float> *vis) = 0;
  virtual void AddInversionDataWithCorrectionCallback(
      GainMode mode, size_t n_polarizations, size_t n_rows, const double *uvws,
      const double *frequencies, size_t n_channels,
      const aocommon::BandData &selected_band,
      const std::pair<size_t, size_t> *antennas,
      const std::complex<float> *visibilities, const size_t *time_offsets,
      MsGridder *gridder, size_t n_antenna) = 0;
  virtual void FinalizeImage(double multiplicationFactor) = 0;
  virtual std::vector<float> RealImage() = 0;
  virtual void InitializePrediction(const float *image_data) = 0;
  virtual void PredictVisibilities(size_t n_rows, size_t n_chan,
                                   const double *uvw, const double *freq,
                                   std::complex<float> *vis) const = 0;
};

/* Memory usage of this gridder is:
   between calls:
     width*height*4  between calls (dirty image buffer)
   during gridding/degridding calls:
     width*height*(
       4      +    (dirty image buffer)
       8*2*2  +    (padded complex uv grid)
       8*2*2  )    (second uv grid used during FFT, not really necessary)
     + nvis_unflagged*8 (index arrays, rough guess)
*/
template <typename NumT>
class WGriddingGridder_Simple final : public WGriddingGridderBase {
 private:
  static constexpr double sigma_min = 1.1;
  static constexpr double sigma_max = 2.0;
  size_t width_, height_, width_t_, height_t_, nthreads_;
  double pixelSizeX_, pixelSizeY_;
  double l_shift_, m_shift_;
  double epsilon_;
  std::vector<NumT> img;
  size_t verbosity_;
  bool tuning_;

 public:
  /** Construct a new gridder with given settings.
   * @param width The width of the untrimmed image in pixels
   * @param height The height of the untrimmed image in pixels.
   * @param width_t The width of the trimmed image in pixels
   * @param height_t The height of the trimmed image in pixels.
   * @param pixelSizeX The angular width of a pixel in radians.
   * @param pixelSizeY The angular height of a pixel in radians.
   * @param nthreads The number of threads to use
   * @param epsilon The requested accuracy of the gridding process.
   *   Affects the support of the employed kernel. Useful values
   *   range between 1e-2 and 1e-6 (for single-precision visibilities).
   * @param verbosity The amount of diagnostic output printed
   *   0: no output
   *   1: print short overview for every inversion/prediction
   *   2: print information for every processed w-plane
   */
  WGriddingGridder_Simple(size_t width, size_t height, size_t width_t,
                          size_t height_t, double pixelSizeX, double pixelSizeY,
                          double l_shift, double m_shift, size_t nthreads,
                          double epsilon = 1e-4, size_t verbosity = 0,
                          bool tuning_ = false);

  WGriddingGridder_Simple(const WGriddingGridder_Simple &) = delete;
  WGriddingGridder_Simple &operator=(const WGriddingGridder_Simple &) = delete;

  /**
   * @return The constant base memory usage of the object in bytes
   */
  size_t ConstantMemoryUsage() const override;
  /**
   * @return Additional memory required per gridded visibility in bytes.
   */
  size_t PerVisibilityMemoryUsage() const override;

  /**
   * Initialize a new inversion gridding pass. This just
   * intializes the accumulated dirty image with zero.
   */
  void InitializeInversion() override;

  /** Add more data to the current inversion operation.
   * The visibilities will be gridded, and the dirty image
   * will be updated accordingly.
   * visibilities with value 0 will be skipped entirely.
   * @param n_chan The number of frequency channels
   * @param uvw pointer to n_rows*3 doubles containing UVW in m.
   *        U(row) := uvw[3*row  ]
   *        V(row) := uvw[3*row+1]
   *        W(row) := uvw[3*row+2]
   * @param freq pointer to n_chan doubles containing channel frequencies
   * @param vis pointer to nrow*n_chan complex<float> containing weighted and
   * corrected visibilities: visibility(row, chan) := vis[row*n_chan + chan]
   */
  void AddInversionData(size_t n_rows, size_t n_chan, const double *uvw,
                        const double *freq,
                        const std::complex<float> *vis) override;
  /** Equivalent to @ref AddInversionData() but without facet solutions
   * pre-applied and with additional paramaters to allow the creation of a
   * callback that can apply solutions "on the fly" as required
   *
   * It is expected that corrections have already been summed via @ref
   * ApplyCorrections<ModifierBehaviour::kSum>()
   *
   * @param n_polarizations The number of polarizations per visibility in @ref
   * visibilities
   * @param antennas Pointer to n_rows `std::pair<size_t, size_t>` containing
   * the antenna pair for each row of visibilities.
   * @param uvws Pointer to n_rows*3 doubles containing UVW in meters.
   *        U(row) := uvws[3*row  ]
   *        V(row) := uvws[3*row+1]
   *        W(row) := uvws[3*row+2]
   * @param frequencies Pointer to n_channels doubles containing channel
   * frequencies
   * @param visibilities Pointer to n_rows*n_channels*n_polarizations
   * `complex<float>` containing weighted but uncorrected and not yet collapsed
   * visibilities: visibility(row, chan) := vis[row*n_chan + chan]
   * @param time_offsets Pointer to n_rows `size_t` containing the time offset
   * as calculated by @ref CacheParmResponse() for the corresponding visibility
   * row when applying @ref ApplyCorrections<ModifierBehaviour::kSum>() on it
   * For further explanation see @ref VisibilityCallbackBuffer::time_offsets_
   * @param gridder Pointer to a gridder that can be called back into in order
   * to apply solutions
   */
  void AddInversionDataWithCorrectionCallback(
      GainMode mode, size_t n_polarizations, size_t n_rows, const double *uvws,
      const double *frequencies, size_t n_channels,
      const aocommon::BandData &selected_band,
      const std::pair<size_t, size_t> *antennas,
      const std::complex<float> *visibilities, const size_t *time_offsets,
      MsGridder *gridder, size_t n_antenna) override;

  /**
   * Finalize inversion once all passes are performed.
   * @param multiplicationFactor Apply this factor to all pixels. This can be
   * used to normalize the image for the weighting scheme.
   */
  void FinalizeImage(double multiplicationFactor) override;

  /**
   * Get the untrimmed image result of inversion. This is an array of size width
   * x height, and can be indexed with [x + y*width]. It is allowed to change
   * this image, e.g. set the horizon to zero before saving to fits. This call
   * is only valid once @ref FinalizeImage() has been called.
   */
  std::vector<float> RealImage() override;

  /**
   * Initialize gridder for prediction and specify image to predict for.
   * @param image The (untrimmed) model image that is to be predicted for. This
   * is an array of width * height size, index by (x + width*y).
   */
  void InitializePrediction(const float *image_data) override;

  /** Predicts visibilities from the current dirty image.
   * FIXME: how do we indicate flagged visibilities that do not
   *        need to be computed? Some special value on input?
   * @param n_rows The number of MS rows being passed
   * @param n_chan The number of frequency channels
   * @param uvw pointer to n_rows*3 doubles containing UVW in m.
   *        U(row) := uvw[3*row  ]
   *        V(row) := uvw[3*row+1]
   *        W(row) := uvw[3*row+2]
   * @param freq pointer to n_chan doubles containing channel frequencies
   * @param vis pointer to nrow*n_chan complex<float> containing weighted
   *        visibilities
   *        visibility(row, chan) := vis[row*n_chan + chan]
   */
  void PredictVisibilities(size_t n_rows, size_t n_chan, const double *uvw,
                           const double *freq,
                           std::complex<float> *vis) const override;

 private:
  /** Internal helper to handle the processing of
   * AddInversionData/AddInversionDataWithCorrectionCallback
   * @tparam TMs Container type for visibilities, must be derived from or
   * compatible with cmav interface
   * @param freq Buffer via which all frequencies for this inversion can be
   * accessed, indexing must be in ascending order
   * @param ms Virtual or actual buffer via which all visibilities for this
   * inversion can be accessed, indexing must be in ascending order
   */
  template <typename Tms>
  void AddInversionMs(size_t n_rows, const double *uvw,
                      const ducc0::cmav<double, 1> &freq, Tms &ms);

  // Helper function to convert mode to a template paramater
  template <typename... Params>
  void AddInversionMs(GainMode mode, Params... params);
  // Helper function to convert polarisations to a template paramater
  template <GainMode Mode, typename... Params>
  void AddInversionMs(size_t n_polarizations, Params... params);
  // Helper function to finally construct the templated object now that we have
  // mode and polarisations as template paramaters
  template <GainMode Mode, size_t NPolarizations, typename... Params>
  void AddInversionMs(size_t n_rows, const double *uvw,
                      const ducc0::cmav<double, 1> &freq, Params... params);
};

#endif

/*
Usage scenario:

WGriddingGridder_Simple gridder(width, height, pixelSizeX, pixelSizeY, nthreads,
1e-5);
// determine number of visibilities that can be gridded in one go, using
// gridder.memUsage() and information about available memory.
// Making the chunks as large as posssible will improve perfrmance.
// Making chunks compact in w will also help a lot.
gridder.InitializeInversion();
for (auto &chunk: chunks)
  gridder.AddInversionData(...)
auto res = gridder.RealImage();
*/

extern template class WGriddingGridder_Simple<float>;
extern template class WGriddingGridder_Simple<double>;

}  // namespace wsclean
