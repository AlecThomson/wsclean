#ifndef GRID_MODE_ENUM_H
#define GRID_MODE_ENUM_H

namespace wsclean {

/** Gridding modes that are supported for interpolating samples on the uv-grid.
 */
enum class GriddingKernelMode {

  /** Simple method that places/samples a visibility on the nearest uv-cell. */
  NearestNeighbour,

  /** Interpolate with a Kaiser-Bessel kernel. This attenuates aliasing. The
   * Kaiser-Bessel window is very similar to the prolate spheroidal kernel,
   * which is considered the optimal gridding window.
   * When this mode is selected, the kernel size and oversampling factor can be
   * specified. This is the recommended and default mode.
   */
  KaiserBessel,

  /** Like KB, but it is not multiplied with a low-pass filtering sinc function.
   */
  KaiserBesselWithoutSinc,

  /** Interpolate with a rectangular window. This will give the sharpest
   * transition at the edge of the image, so will maximally attenuate objects
   * just outside the image. However, objects further from the edge will not be
   * as much attenuated compared to the KB window, which has much deeper
   * sidelobes further out.
   */
  Rectangular,

  /** Window the low-pass filter with a Gaussian trimmed at 3 sigma.
   */
  Gaussian,
  GaussianWithoutSinc,

  /** Blackman-nutall window.
   */
  BlackmanNuttall,
  BlackmanNuttallWithoutSinc,

  /** Blackman-Harris window.
   */
  BlackmanHarris
};

}  // namespace wsclean

#endif
