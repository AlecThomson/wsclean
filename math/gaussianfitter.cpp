#include "gaussianfitter.h"

#include <cmath>

#include <aocommon/matrix2x2.h>
#include <aocommon/uvector.h>

using aocommon::Matrix2x2;

namespace {
void ToAnglesAndFWHM(double sx, double sy, double beta, double& ellipseMaj,
                     double& ellipseMin, double& ellipsePA) {
  const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
  const double betaFact = 1.0 - beta * beta;
  double cov[4];
  cov[0] = sx * sx / betaFact;
  cov[1] = beta * sx * sy / betaFact;
  cov[2] = cov[1];
  cov[3] = sy * sy / betaFact;

  double e1, e2, vec1[2], vec2[2];
  Matrix2x2::EigenValuesAndVectors(cov, e1, e2, vec1, vec2);
  if (std::isfinite(e1)) {
    ellipseMaj = sqrt(std::fabs(e1)) * sigmaToBeam;
    ellipseMin = sqrt(std::fabs(e2)) * sigmaToBeam;
    if (ellipseMaj < ellipseMin) {
      std::swap(ellipseMaj, ellipseMin);
      vec1[0] = vec2[0];
      vec1[1] = vec2[1];
    }
    ellipsePA = -atan2(vec1[0], vec1[1]);
  } else {
    ellipseMaj = sqrt(std::fabs(sx)) * sigmaToBeam;
    ellipseMin = sqrt(std::fabs(sx)) * sigmaToBeam;
    ellipsePA = 0.0;
  }
}

/**
 * Calculates the difference between a gaussian with the specified parameters
 * at position x,y and the given value.
 */
double ErrCentered(double val, double x, double y, double sx, double sy,
                   double beta) {
  return std::exp(-x * x / (2.0 * sx * sx) + beta * x * y / (sx * sy) -
                  y * y / (2.0 * sy * sy)) -
         val;
}

/**
 * Calculates the difference between a gaussian with the specified parameters
 * at position x,y and the given value.
 */
double ErrCircularCentered(double val, double x, double y, double s) {
  return std::exp((-x * x - y * y) / (2.0 * s * s)) - val;
}

double ErrFull(double val, double v, double x, double y, double sx, double sy,
               double beta) {
  return std::exp(-x * x / (2.0 * sx * sx) + beta * x * y / (sx * sy) -
                  y * y / (2.0 * sy * sy)) *
             v -
         val;
}

/**
 * Fitting function for fit2DGaussianCentred(). Calculates the sum of the
 * squared errors(/residuals).
 */
int FittingFuncCentered(const gsl_vector* xvec, void* data, gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double sx = gsl_vector_get(xvec, 0);
  const double sy = gsl_vector_get(xvec, 1);
  const double beta = gsl_vector_get(xvec, 2);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int xMid = width / 2;
  const int yMid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t dataIndex = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - yMid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - xMid) * scale;
      double e = ErrCentered(fitter.Image()[dataIndex], x, y, sx, sy, beta);
      gsl_vector_set(f, dataIndex, e);
      ++dataIndex;
    }
  }
  return GSL_SUCCESS;
}

int FittingFuncCircularCentered(const gsl_vector* xvec, void* data,
                                gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double s = gsl_vector_get(xvec, 0);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int xMid = width / 2;
  const int yMid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t dataIndex = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - yMid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - xMid) * scale;
      double e = ErrCircularCentered(fitter.Image()[dataIndex], x, y, s);
      gsl_vector_set(f, dataIndex, e);
      ++dataIndex;
    }
  }
  return GSL_SUCCESS;
}

/**
 * Derivative function belong with fit2DGaussianCentred().
 */
int FittingDerivCentered(const gsl_vector* xvec, void* data, gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double sx = gsl_vector_get(xvec, 0);
  const double sy = gsl_vector_get(xvec, 1);
  const double beta = gsl_vector_get(xvec, 2);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int xMid = width / 2;
  const int yMid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t dataIndex = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - yMid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - xMid) * scale;
      double expTerm =
          std::exp(-x * x / (2.0 * sx * sx) + beta * x * y / (sx * sy) -
                   y * y / (2.0 * sy * sy));
      double dsx =
          (beta * x * y / (sx * sx * sy) + x * x / (sx * sx * sx)) * expTerm;
      double dsy =
          (beta * x * y / (sy * sy * sx) + y * y / (sy * sy * sy)) * expTerm;
      double dbeta = x * y / (sx * sy) * expTerm;
      gsl_matrix_set(J, dataIndex, 0, dsx);
      gsl_matrix_set(J, dataIndex, 1, dsy);
      gsl_matrix_set(J, dataIndex, 2, dbeta);
      ++dataIndex;
    }
  }
  return GSL_SUCCESS;
}

int FittingDerivCircularCentered(const gsl_vector* xvec, void* data,
                                 gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double s = gsl_vector_get(xvec, 0);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int xMid = width / 2, yMid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t dataIndex = 0;
  for (size_t yi = 0; yi != height; ++yi) {
    double y = (yi - yMid) * scale;
    for (size_t xi = 0; xi != width; ++xi) {
      double x = (xi - xMid) * scale;
      double expTerm = std::exp((-x * x - y * y) / (2.0 * s * s));
      // derivative of exp((-x*x - y*y)/(2.0*s*s)) to s
      // = (-x*x - y*y)/2.0*-2/(s*s*s)
      // = (-x*x - y*y)/(-s*s*s)
      // = (x*x + y*y)/(s*s*s)
      double ds = ((x * x + y * y) / (s * s * s)) * expTerm;
      gsl_matrix_set(J, dataIndex, 0, ds);
      ++dataIndex;
    }
  }
  return GSL_SUCCESS;
}

/**
 * Squared error and derivative function together.
 */
int FittingBothCentered(const gsl_vector* x, void* data, gsl_vector* f,
                        gsl_matrix* J) {
  FittingFuncCentered(x, data, f);
  FittingDerivCentered(x, data, J);
  return GSL_SUCCESS;
}

int FittingBothCircularCentered(const gsl_vector* x, void* data, gsl_vector* f,
                                gsl_matrix* J) {
  FittingFuncCircularCentered(x, data, f);
  FittingDerivCircularCentered(x, data, J);
  return GSL_SUCCESS;
}

int fitting_func_with_amplitude(const gsl_vector* xvec, void* data,
                                gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int xMid = width / 2;
  const int yMid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t dataIndex = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double yS = yc + (yi - yMid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      double xS = xc + (xi - xMid) * scale;
      double e = ErrFull(fitter.Image()[dataIndex], v, xS, yS, sx, sy, beta);
      gsl_vector_set(f, dataIndex, e);
      ++dataIndex;
    }
  }
  return GSL_SUCCESS;
}

int fitting_deriv_with_amplitude(const gsl_vector* xvec, void* data,
                                 gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double scale = 1.0 / fitter.ScaleFactor();
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  if (fitter.PosConstrained() != 0.0 &&
      (std::fabs(xc - fitter.XInit()) > fitter.PosConstrained() * scale ||
       std::fabs(yc - fitter.YInit()) > fitter.PosConstrained() * scale)) {
    std::cout << "GSL_EDOM\n";
    return GSL_EDOM;
  }
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int xMid = width / 2;
  const int yMid = height / 2;

  size_t dataIndex = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double y = yc + (yi - yMid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      // TODO I need to go over the signs -- ds, dy, dsx, dsy in particular
      double x = xc + (xi - xMid) * scale;
      double expTerm =
          std::exp(-x * x / (2.0 * sx * sx) + beta * x * y / (sx * sy) -
                   y * y / (2.0 * sy * sy));
      double dv = expTerm;
      expTerm *= v;
      double dx = (-beta * y / (sx * sy) - x / (sx * sx)) * expTerm;
      double dy = (-beta * x / (sy * sx) - y / (sy * sy)) * expTerm;
      double dsx =
          (beta * x * y / (sx * sx * sy) + x * x / (sx * sx * sx)) * expTerm;
      double dsy =
          (beta * x * y / (sy * sy * sx) + y * y / (sy * sy * sy)) * expTerm;
      double dbeta = x * y / (sx * sy) * expTerm;
      gsl_matrix_set(J, dataIndex, 0, dv);
      gsl_matrix_set(J, dataIndex, 1, dx);
      gsl_matrix_set(J, dataIndex, 2, dy);
      gsl_matrix_set(J, dataIndex, 3, dsx);
      gsl_matrix_set(J, dataIndex, 4, dsy);
      gsl_matrix_set(J, dataIndex, 5, dbeta);
      ++dataIndex;
    }
  }
  return GSL_SUCCESS;
}

int fitting_both_with_amplitude(const gsl_vector* x, void* data, gsl_vector* f,
                                gsl_matrix* J) {
  fitting_func_with_amplitude(x, data, f);
  fitting_deriv_with_amplitude(x, data, J);
  return GSL_SUCCESS;
}
int fitting_func_with_amplitude_and_floor(const gsl_vector* xvec, void* data,
                                          gsl_vector* f) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double scale = 1.0 / fitter.ScaleFactor();
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  const double fl = gsl_vector_get(xvec, 6);
  if (fitter.PosConstrained() != 0.0 &&
      (std::fabs(xc - fitter.XInit()) > fitter.PosConstrained() * scale ||
       std::fabs(yc - fitter.YInit()) > fitter.PosConstrained() * scale))
    return GSL_EDOM;
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int xMid = width / 2;
  const int yMid = height / 2;

  size_t dataIndex = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double yS = yc + (yi - yMid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      double xS = xc + (xi - xMid) * scale;
      double e =
          ErrFull(fitter.Image()[dataIndex], v, xS, yS, sx, sy, beta) + fl;
      gsl_vector_set(f, dataIndex, e);
      ++dataIndex;
    }
  }
  return GSL_SUCCESS;
}

int fitting_deriv_with_amplitude_and_floor(const gsl_vector* xvec, void* data,
                                           gsl_matrix* J) {
  const GaussianFitter& fitter = *static_cast<const GaussianFitter*>(data);
  const double v = gsl_vector_get(xvec, 0);
  const double xc = gsl_vector_get(xvec, 1);
  const double yc = gsl_vector_get(xvec, 2);
  const double sx = gsl_vector_get(xvec, 3);
  const double sy = gsl_vector_get(xvec, 4);
  const double beta = gsl_vector_get(xvec, 5);
  const size_t width = fitter.Width();
  const size_t height = fitter.Height();
  const int xMid = width / 2;
  const int yMid = height / 2;
  const double scale = 1.0 / fitter.ScaleFactor();

  size_t dataIndex = 0;
  for (int yi = 0; yi != int(height); ++yi) {
    double y = yc + (yi - yMid) * scale;
    for (int xi = 0; xi != int(width); ++xi) {
      double x = xc + (xi - xMid) * scale;
      double expTerm = exp(-x * x / (2.0 * sx * sx) + beta * x * y / (sx * sy) -
                           y * y / (2.0 * sy * sy));
      double dv = expTerm;
      expTerm *= v;
      double dx = (-beta * y / (sx * sy) - x / (sx * sx)) * expTerm;
      double dy = (-beta * x / (sy * sx) - y / (sy * sy)) * expTerm;
      double dsx =
          (beta * x * y / (sx * sx * sy) + x * x / (sx * sx * sx)) * expTerm;
      double dsy =
          (beta * x * y / (sy * sy * sx) + y * y / (sy * sy * sy)) * expTerm;
      double dbeta = x * y / (sx * sy) * expTerm;
      double dfl = 1.0;
      gsl_matrix_set(J, dataIndex, 0, dv);
      gsl_matrix_set(J, dataIndex, 1, dx);
      gsl_matrix_set(J, dataIndex, 2, dy);
      gsl_matrix_set(J, dataIndex, 3, dsx);
      gsl_matrix_set(J, dataIndex, 4, dsy);
      gsl_matrix_set(J, dataIndex, 5, dbeta);
      gsl_matrix_set(J, dataIndex, 6, dfl);
      ++dataIndex;
    }
  }
  return GSL_SUCCESS;
}

int fitting_both_with_amplitude_and_floor(const gsl_vector* x, void* data,
                                          gsl_vector* f, gsl_matrix* J) {
  fitting_func_with_amplitude_and_floor(x, data, f);
  fitting_deriv_with_amplitude_and_floor(x, data, J);
  return GSL_SUCCESS;
}

}  // namespace

void GaussianFitter::Fit2DGaussianCentred(const float* image, size_t width,
                                          size_t height, double beamEst,
                                          double& beamMaj, double& beamMin,
                                          double& beamPA, double boxScaleFactor,
                                          bool verbose) {
  size_t prefSize = std::max<size_t>(std::ceil(boxScaleFactor),
                                     std::ceil(beamEst * boxScaleFactor));
  if (prefSize % 2 != 0) ++prefSize;
  if (prefSize < width || prefSize < height) {
    size_t nIter = 0;
    bool boxWasLargeEnough;
    do {
      size_t boxWidth = std::min(prefSize, width);
      size_t boxHeight = std::min(prefSize, height);
      if (verbose) std::cout << "Fit initial value:" << beamEst << "\n";
      fit2DGaussianCentredInBox(image, width, height, beamEst, beamMaj, beamMin,
                                beamPA, boxWidth, boxHeight, verbose);
      if (verbose)
        std::cout << "Fit result:" << beamMaj << " x " << beamMin << " px, "
                  << beamPA << " (box was " << boxWidth << " x " << boxHeight
                  << ")\n";

      boxWasLargeEnough =
          (beamMaj * boxScaleFactor * 0.8 < boxWidth || boxWidth >= width) &&
          (beamMaj * boxScaleFactor * 0.8 < boxHeight || boxHeight >= height);
      if (!boxWasLargeEnough) {
        prefSize = std::max<size_t>(std::ceil(boxScaleFactor),
                                    std::ceil(beamMaj * boxScaleFactor));
        if (prefSize % 2 != 0) ++prefSize;
        beamEst = std::max(beamMaj, beamEst);
      }
      ++nIter;
    } while (!boxWasLargeEnough && nIter < 5);
  } else {
    if (verbose) std::cout << "Image is as large as the fitting box.\n";
    fit2DGaussianCentred(image, width, height, beamEst, beamMaj, beamMin,
                         beamPA, verbose);
  }
}

void GaussianFitter::Fit2DCircularGaussianCentred(const float* image,
                                                  size_t width, size_t height,
                                                  double& beamSize,
                                                  double boxScaleFactor) {
  double initialValue = beamSize;
  size_t prefSize = std::max<size_t>(std::ceil(boxScaleFactor),
                                     std::ceil(beamSize * boxScaleFactor));
  if (prefSize % 2 != 0) ++prefSize;
  if (prefSize < width || prefSize < height) {
    size_t boxWidth = std::min(prefSize, width);
    size_t boxHeight = std::min(prefSize, height);
    size_t nIter = 0;
    bool boxWasLargeEnough;
    do {
      fit2DCircularGaussianCentredInBox(image, width, height, beamSize,
                                        boxWidth, boxHeight);

      boxWasLargeEnough =
          (beamSize * boxScaleFactor * 0.8 < boxWidth || width >= boxWidth) &&
          (beamSize * boxScaleFactor * 0.8 < boxHeight || height >= boxHeight);
      if (!boxWasLargeEnough) {
        prefSize = std::max<size_t>(std::ceil(boxScaleFactor),
                                    std::ceil(beamSize * boxScaleFactor));
        if (prefSize % 2 != 0) ++prefSize;
        beamSize = std::max(initialValue, beamSize);
      }
      ++nIter;
    } while (!boxWasLargeEnough && nIter < 5);
  } else {
    fit2DCircularGaussianCentred(image, width, height, beamSize);
  }
}

void GaussianFitter::Fit2DGaussianFull(const float* image, size_t width,
                                       size_t height, double& val, double& posX,
                                       double& posY, double& beamMaj,
                                       double& beamMin, double& beamPA,
                                       double* floorLevel) {
  size_t prefSize = std::max<size_t>(10, std::ceil(beamMaj * 10.0));
  if (prefSize % 2 != 0) ++prefSize;
  if (prefSize < width || prefSize < height) {
    size_t xStart = std::max<int>(0, int(round(posX)) - int(prefSize) / 2);
    size_t xEnd = std::min(width, size_t(round(posX)) + prefSize / 2);
    size_t yStart = std::max<int>(0, int(round(posY)) - int(prefSize) / 2);
    size_t yEnd = std::min(height, size_t(round(posY)) + prefSize / 2);
    size_t nIter = 0;
    bool boxWasLargeEnough;
    do {
      fit2DGaussianWithAmplitudeInBox(image, width, height, val, posX, posY,
                                      beamMaj, beamMin, beamPA, floorLevel,
                                      xStart, xEnd, yStart, yEnd);

      size_t boxWidth = xEnd - xStart;
      size_t boxHeight = yEnd - yStart;
      boxWasLargeEnough = (beamMaj * 4.0 < boxWidth || width >= boxWidth) &&
                          (beamMaj * 4.0 < boxHeight || height >= boxHeight);
      if (!boxWasLargeEnough) {
        prefSize = std::max<size_t>(10, std::ceil(beamMaj * 10.0));
        if (prefSize % 2 != 0) ++prefSize;
      }
      ++nIter;
    } while (!boxWasLargeEnough && nIter < 5);
  } else {
    fit2DGaussianWithAmplitude(image, width, height, val, posX, posY, beamMaj,
                               beamMin, beamPA, floorLevel);
  }
}

void GaussianFitter::fit2DGaussianCentredInBox(const float* image, size_t width,
                                               size_t height, double beamEst,
                                               double& beamMaj, double& beamMin,
                                               double& beamPA, size_t boxWidth,
                                               size_t boxHeight, bool verbose) {
  size_t startX = (width - boxWidth) / 2;
  size_t startY = (height - boxHeight) / 2;
  aocommon::UVector<float> smallImage(boxWidth * boxHeight);
  for (size_t y = startY; y != (height + boxHeight) / 2; ++y) {
    std::copy_n(&image[y * width + startX], boxWidth,
                &smallImage[(y - startY) * boxWidth]);
  }

  fit2DGaussianCentred(&smallImage[0], boxWidth, boxHeight, beamEst, beamMaj,
                       beamMin, beamPA, verbose);
}

void GaussianFitter::fit2DCircularGaussianCentredInBox(
    const float* image, size_t width, size_t height, double& beamSize,
    size_t boxWidth, size_t boxHeight) {
  size_t startX = (width - boxWidth) / 2;
  size_t startY = (height - boxHeight) / 2;
  aocommon::UVector<float> smallImage(boxWidth * boxHeight);
  for (size_t y = startY; y != (height + boxHeight) / 2; ++y) {
    std::copy_n(&image[y * width + startX], boxWidth,
                &smallImage[(y - startY) * boxWidth]);
  }

  fit2DCircularGaussianCentred(&smallImage[0], boxWidth, boxHeight, beamSize);
}

void GaussianFitter::fit2DGaussianCentred(const float* image, size_t width,
                                          size_t height, double beamEst,
                                          double& beamMaj, double& beamMin,
                                          double& beamPA, bool verbose) {
  _width = width;
  _height = height;
  _image = image;
  _scaleFactor = (width + height) / 2;

  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, _width * _height, 3);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingFuncCentered;
  fdf.df = &FittingDerivCentered;
  fdf.fdf = &FittingBothCentered;
  fdf.n = _width * _height;
  fdf.p = 3;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
  double initialValsArray[3] = {beamEst / (_scaleFactor * double(sigmaToBeam)),
                                beamEst / (_scaleFactor * double(sigmaToBeam)),
                                0.0};
  gsl_vector_view initialVals = gsl_vector_view_array(initialValsArray, 3);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initialVals.vector);

  int status;
  size_t iter = 0;
  do {
    if (verbose) std::cout << "Iteration " << iter << ": ";
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  double sx = gsl_vector_get(solver->x, 0), sy = gsl_vector_get(solver->x, 1),
         beta = gsl_vector_get(solver->x, 2);

  gsl_multifit_fdfsolver_free(solver);

  ToAnglesAndFWHM(sx, sy, beta, beamMaj, beamMin, beamPA);
  beamMaj *= _scaleFactor;
  beamMin *= _scaleFactor;
}

void GaussianFitter::fit2DCircularGaussianCentred(const float* image,
                                                  size_t width, size_t height,
                                                  double& beamSize) {
  _width = width;
  _height = height;
  _image = image;
  _scaleFactor = (width + height) / 2;

  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, _width * _height, 1);

  gsl_multifit_function_fdf fdf;
  fdf.f = &FittingFuncCircularCentered;
  fdf.df = &FittingDerivCircularCentered;
  fdf.fdf = &FittingBothCircularCentered;
  fdf.n = _width * _height;
  fdf.p = 1;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  const long double kSigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
  double initialValsArray[1] = {beamSize /
                                (_scaleFactor * double(kSigmaToBeam))};
  gsl_vector_view initialVals = gsl_vector_view_array(initialValsArray, 1);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initialVals.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  const double s = gsl_vector_get(solver->x, 0);
  gsl_multifit_fdfsolver_free(solver);

  beamSize = s * kSigmaToBeam * _scaleFactor;
}

void GaussianFitter::fit2DGaussianWithAmplitudeInBox(
    const float* image, size_t width, size_t /*height*/, double& val,
    double& posX, double& posY, double& beamMaj, double& beamMin,
    double& beamPA, double* floorLevel, size_t xStart, size_t xEnd,
    size_t yStart, size_t yEnd) {
  size_t boxWidth = xEnd - xStart;
  size_t boxHeight = yEnd - yStart;
  aocommon::UVector<float> smallImage(boxWidth * boxHeight);
  for (size_t y = yStart; y != yEnd; ++y) {
    std::copy_n(&image[y * width + xStart], boxWidth,
                &smallImage[(y - yStart) * boxWidth]);
  }

  posX -= xStart;
  posY -= yStart;
  fit2DGaussianWithAmplitude(&smallImage[0], boxWidth, boxHeight, val, posX,
                             posY, beamMaj, beamMin, beamPA, floorLevel);
  posX += xStart;
  posY += yStart;
}

/**
 * Fits the position, size and amplitude of a Gaussian. If floorLevel is not
 * a nullptr, the floor (background level, or zero level) is fitted too.
 */
void GaussianFitter::fit2DGaussianWithAmplitude(const float* image,
                                                size_t width, size_t height,
                                                double& val, double& posX,
                                                double& posY, double& beamMaj,
                                                double& beamMin, double& beamPA,
                                                double* floorLevel) {
  _width = width;
  _height = height;
  _image = image;
  _scaleFactor = (width + height) / 2;

  if (floorLevel == nullptr)
    fit2DGaussianWithAmplitude(val, posX, posY, beamMaj, beamMin, beamPA);
  else
    fit2DGaussianWithAmplitudeWithFloor(val, posX, posY, beamMaj, beamMin,
                                        beamPA, *floorLevel);
}

void GaussianFitter::fit2DGaussianWithAmplitude(double& val, double& posX,
                                                double& posY, double& beamMaj,
                                                double& beamMin,
                                                double& beamPA) {
  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, _width * _height, 6);

  gsl_multifit_function_fdf fdf;
  fdf.f = &fitting_func_with_amplitude;
  fdf.df = &fitting_deriv_with_amplitude;
  fdf.fdf = &fitting_both_with_amplitude;
  fdf.n = _width * _height;
  fdf.p = 6;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
  _xInit = -(posX - _width / 2) / _scaleFactor;
  _yInit = -(posY - _height / 2) / _scaleFactor;
  double initialValsArray[6] = {val,
                                _xInit,
                                _yInit,
                                beamMaj / (_scaleFactor * double(sigmaToBeam)),
                                beamMaj / (_scaleFactor * double(sigmaToBeam)),
                                0.0};
  gsl_vector_view initialVals = gsl_vector_view_array(initialValsArray, 6);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initialVals.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  val = gsl_vector_get(solver->x, 0);
  posX = -1.0 * gsl_vector_get(solver->x, 1) * _scaleFactor + _width / 2;
  posY = -1.0 * gsl_vector_get(solver->x, 2) * _scaleFactor + _height / 2;
  double sx = gsl_vector_get(solver->x, 3), sy = gsl_vector_get(solver->x, 4),
         beta = gsl_vector_get(solver->x, 5);

  gsl_multifit_fdfsolver_free(solver);

  ToAnglesAndFWHM(sx, sy, beta, beamMaj, beamMin, beamPA);
  beamMaj *= _scaleFactor;
  beamMin *= _scaleFactor;
}

void GaussianFitter::fit2DGaussianWithAmplitudeWithFloor(
    double& val, double& posX, double& posY, double& beamMaj, double& beamMin,
    double& beamPA, double& floorLevel) {
  const gsl_multifit_fdfsolver_type* T = gsl_multifit_fdfsolver_lmsder;
  gsl_multifit_fdfsolver* solver =
      gsl_multifit_fdfsolver_alloc(T, _width * _height, 7);

  gsl_multifit_function_fdf fdf;
  fdf.f = &fitting_func_with_amplitude_and_floor;
  fdf.df = &fitting_deriv_with_amplitude_and_floor;
  fdf.fdf = &fitting_both_with_amplitude_and_floor;
  fdf.n = _width * _height;
  fdf.p = 7;
  fdf.params = this;

  // Using the FWHM formula for a Gaussian:
  const long double sigmaToBeam = 2.0L * sqrtl(2.0L * logl(2.0L));
  _xInit = -(posX - _width / 2) / _scaleFactor;
  _yInit = -(posY - _height / 2) / _scaleFactor;
  double initialValsArray[7] = {val,
                                _xInit,
                                _yInit,
                                beamMaj / (_scaleFactor * double(sigmaToBeam)),
                                beamMaj / (_scaleFactor * double(sigmaToBeam)),
                                0.0,
                                0.0};
  gsl_vector_view initialVals = gsl_vector_view_array(initialValsArray, 7);
  gsl_multifit_fdfsolver_set(solver, &fdf, &initialVals.vector);

  int status;
  size_t iter = 0;
  do {
    iter++;
    status = gsl_multifit_fdfsolver_iterate(solver);

    if (status) break;

    status = gsl_multifit_test_delta(solver->dx, solver->x, 1e-7, 1e-7);

  } while (status == GSL_CONTINUE && iter < 500);

  val = gsl_vector_get(solver->x, 0);
  posX = -1.0 * gsl_vector_get(solver->x, 1) * _scaleFactor + _width / 2;
  posY = -1.0 * gsl_vector_get(solver->x, 2) * _scaleFactor + _height / 2;
  double sx = gsl_vector_get(solver->x, 3), sy = gsl_vector_get(solver->x, 4),
         beta = gsl_vector_get(solver->x, 5);
  floorLevel = gsl_vector_get(solver->x, 6);

  gsl_multifit_fdfsolver_free(solver);

  ToAnglesAndFWHM(sx, sy, beta, beamMaj, beamMin, beamPA);
  beamMaj *= _scaleFactor;
  beamMin *= _scaleFactor;
}
