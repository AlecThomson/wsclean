#ifndef GAUSSIAN_FITTER_H
#define GAUSSIAN_FITTER_H

#include <cstring>
#include <iostream>

class GaussianFitter {
 public:
  GaussianFitter()
      : image_(nullptr),
        width_(0),
        height_(0),
        scale_factor_(0),
        x_init_(0.0),
        y_init_(0.0),
        pos_constrained_(0.0) {}

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

  const float* Image() const { return image_; }
  size_t Width() const { return width_; }
  size_t Height() const { return height_; }
  size_t ScaleFactor() const { return scale_factor_; }
  double XInit() const { return x_init_; };
  double YInit() const { return y_init_; };
  double PosConstrained() const { return pos_constrained_; };

 private:
  void Fit2DGaussianCentredInBox(const float* image, size_t width,
                                 size_t height, double beamEst, double& beamMaj,
                                 double& beamMin, double& beamPA,
                                 size_t boxWidth, size_t boxHeight,
                                 bool verbose);

  void Fit2DCircularGaussianCentredInBox(const float* image, size_t width,
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
  void Fit2DGaussianCentred(const float* image, size_t width, size_t height,
                            double beamEst, double& beamMaj, double& beamMin,
                            double& beamPA, bool verbose);

  void Fit2DCircularGaussianCentred(const float* image, size_t width,
                                    size_t height, double& beamSize);

  void Fit2DGaussianWithAmplitudeInBox(const float* image, size_t width,
                                       size_t height, double& val, double& posX,
                                       double& posY, double& beamMaj,
                                       double& beamMin, double& beamPA,
                                       double* floorLevel, size_t xStart,
                                       size_t xEnd, size_t yStart, size_t yEnd);

  /**
   * Fits the position, size and amplitude of a Gaussian. If floorLevel is not
   * a nullptr, the floor (background level, or zero level) is fitted too.
   */
  void Fit2DGaussianWithAmplitude(const float* image, size_t width,
                                  size_t height, double& val, double& posX,
                                  double& posY, double& beamMaj,
                                  double& beamMin, double& beamPA,
                                  double* floorLevel);

  /**
   * Like Fit2DGaussianCentred(), but includes Gaussian centre X and Y position
   * and amplitude in the fitted parameters.
   *
   * This function can typically be used for source fitting.
   */
  void Fit2DGaussianWithAmplitude(double& val, double& posX, double& posY,
                                  double& beamMaj, double& beamMin,
                                  double& beamPA);

  /**
   * Like Fit2DGaussianWithAmplitude(), but includes floorLevel as fitted
   * parameter. Floor is the background/zero level on which the Gaussian
   * resides.
   */
  void Fit2DGaussianWithAmplitudeWithFloor(double& val, double& posX,
                                           double& posY, double& beamMaj,
                                           double& beamMin, double& beamPA,
                                           double& floorLevel);

  const float* image_;
  size_t width_;
  size_t height_;
  size_t scale_factor_;
  double x_init_;
  double y_init_;
  double pos_constrained_;
};

#endif
