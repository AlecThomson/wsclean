#ifndef GAUSSIAN_FITTER_H
#define GAUSSIAN_FITTER_H

#include <cstring>
#include <iostream>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>

class GaussianFitter {
 public:
  GaussianFitter() : _posConstrained(0.0) {}

  void Fit2DGaussianCentred(const float* image, size_t width, size_t height,
                            double beamEst, double& beamMaj, double& beamMin,
                            double& beamPA, double boxScaleFactor = 10.0,
                            bool verbose = false);

  void Fit2DCircularGaussianCentred(const float* image, size_t width,
                                    size_t height, double& beamSize,
                                    double boxScaleFactor = 10.0);

  void Fit2DGaussianFull(const float* image, size_t width, size_t height,
                         double& val, double& posX, double& posY,
                         double& beamMaj, double& beamMin, double& beamPA,
                         double* floorLevel = nullptr);

  const float* Image() const { return _image; }
  size_t Width() const { return _width; }
  size_t Height() const { return _height; }
  size_t ScaleFactor() const { return _scaleFactor; }

 private:
  const float* _image;
  size_t _width, _height, _scaleFactor;

  void fit2DGaussianCentredInBox(const float* image, size_t width,
                                 size_t height, double beamEst, double& beamMaj,
                                 double& beamMin, double& beamPA,
                                 size_t boxWidth, size_t boxHeight,
                                 bool verbose);

  void fit2DCircularGaussianCentredInBox(const float* image, size_t width,
                                         size_t height, double& beamSize,
                                         size_t boxWidth, size_t boxHeight);

  /**
   * This function performs a single fit of a Gaussian. The position of the
   * Gaussian is constrained to be in the centre of the image. The Gaussian is
   * fitted such that the squared residuals (data - model) are minimal.
   *
   * This function is typically used to find the beam-shape of the point-spread
   * function. The beam estimate is used as initial value for the minor and
   * major shape.
   */
  void fit2DGaussianCentred(const float* image, size_t width, size_t height,
                            double beamEst, double& beamMaj, double& beamMin,
                            double& beamPA, bool verbose);

  void fit2DCircularGaussianCentred(const float* image, size_t width,
                                    size_t height, double& beamSize);

  void fit2DGaussianWithAmplitudeInBox(const float* image, size_t width,
                                       size_t height, double& val, double& posX,
                                       double& posY, double& beamMaj,
                                       double& beamMin, double& beamPA,
                                       double* floorLevel, size_t xStart,
                                       size_t xEnd, size_t yStart, size_t yEnd);

  /**
   * Fits the position, size and amplitude of a Gaussian. If floorLevel is not
   * a nullptr, the floor (background level, or zero level) is fitted too.
   */
  void fit2DGaussianWithAmplitude(const float* image, size_t width,
                                  size_t height, double& val, double& posX,
                                  double& posY, double& beamMaj,
                                  double& beamMin, double& beamPA,
                                  double* floorLevel);

  /**
   * Like fit2DGaussianCentred(), but includes Gaussian centre X and Y position
   * and amplitude in the fitted parameters.
   *
   * This function can typically be used for source fitting.
   */
  void fit2DGaussianWithAmplitude(double& val, double& posX, double& posY,
                                  double& beamMaj, double& beamMin,
                                  double& beamPA);

  /**
   * Like fit2DGaussianWithAmplitude(), but includes floorLevel as fitted
   * parameter. Floor is the background/zero level on which the Gaussian
   * resides.
   */
  void fit2DGaussianWithAmplitudeWithFloor(double& val, double& posX,
                                           double& posY, double& beamMaj,
                                           double& beamMin, double& beamPA,
                                           double& floorLevel);

  static int fitting_func_with_amplitude(const gsl_vector* xvec, void* data,
                                         gsl_vector* f);

  static int fitting_deriv_with_amplitude(const gsl_vector* xvec, void* data,
                                          gsl_matrix* J);

  static int fitting_both_with_amplitude(const gsl_vector* x, void* data,
                                         gsl_vector* f, gsl_matrix* J) {
    fitting_func_with_amplitude(x, data, f);
    fitting_deriv_with_amplitude(x, data, J);
    return GSL_SUCCESS;
  }

  static int fitting_func_with_amplitude_and_floor(const gsl_vector* xvec,
                                                   void* data, gsl_vector* f);

  static int fitting_deriv_with_amplitude_and_floor(const gsl_vector* xvec,
                                                    void* data, gsl_matrix* J);

  static int fitting_both_with_amplitude_and_floor(const gsl_vector* x,
                                                   void* data, gsl_vector* f,
                                                   gsl_matrix* J) {
    fitting_func_with_amplitude_and_floor(x, data, f);
    fitting_deriv_with_amplitude_and_floor(x, data, J);
    return GSL_SUCCESS;
  }

  double _xInit, _yInit, _posConstrained;
};

#endif
