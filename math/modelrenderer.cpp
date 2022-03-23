
#include "modelrenderer.h"

#include <aocommon/imagecoordinates.h>
#include <aocommon/uvector.h>

#include <schaapcommon/fft/convolution.h>
#include <schaapcommon/fft/restoreimage.h>

#include <cmath>

#include <boost/algorithm/clamp.hpp>

using aocommon::ImageCoordinates;
using boost::algorithm::clamp;

namespace {
long double Gaussian(long double x, long double sigma) {
  const long double xi = x / sigma;
  // TODO: ask Andre: why is the second part not evaluated? Current formulation
  // doesn't entirely match definition of a Gaussian.
  return std::exp(static_cast<long double>(-0.5) * xi * xi);
  // / (sigma * sqrt(static_cast<long double>(2.0) * M_PIl));
}

void RenderGaussianComponent(float* imageData, size_t imageWidth,
                             size_t imageHeight, long double phaseCentreRA,
                             long double phaseCentreDec,
                             long double pixelScaleL, long double pixelScaleM,
                             long double phaseCentreDL,
                             long double phaseCentreDM, long double posRA,
                             long double posDec, long double gausMaj,
                             long double gausMin, long double gausPA,
                             long double flux) {
  // Using the FWHM formula for a Gaussian:
  const long double fwhm_constant = (2.0L * std::sqrt(2.0L * std::log(2.0L)));
  const long double sigmaMaj = gausMaj / fwhm_constant;
  const long double sigmaMin = gausMin / fwhm_constant;
  // TODO this won't work for non-equally spaced dimensions
  const long double minPixelScale = std::min(pixelScaleL, pixelScaleM);

  // Make rotation matrix
  long double transf[4];
  // Position angle is angle from North:
  const long double angle = gausPA + 0.5 * M_PI;
  transf[2] = std::sin(angle);
  transf[0] = std::cos(angle);
  transf[1] = -transf[2];
  transf[3] = transf[0];
  const double sigmaMax = std::max(std::fabs(sigmaMaj * transf[0]),
                                   std::fabs(sigmaMaj * transf[1]));
  // Multiply with scaling matrix to make variance 1.
  transf[0] = transf[0] / sigmaMaj;
  transf[1] = transf[1] / sigmaMaj;
  transf[2] = transf[2] / sigmaMin;
  transf[3] = transf[3] / sigmaMin;
  const int boundingBoxSize = std::ceil(sigmaMax * 20.0 / minPixelScale);
  long double sourceL, sourceM;
  ImageCoordinates::RaDecToLM(posRA, posDec, phaseCentreRA, phaseCentreDec,
                              sourceL, sourceM);

  // Calculate the bounding box
  int sourceX, sourceY;
  ImageCoordinates::LMToXY<long double>(
      sourceL - phaseCentreDL, sourceM - phaseCentreDM, pixelScaleL,
      pixelScaleM, imageWidth, imageHeight, sourceX, sourceY);
  const int xLeft = clamp(sourceX - boundingBoxSize, 0, int(imageWidth));
  const int xRight = clamp(sourceX + boundingBoxSize, xLeft, int(imageWidth));
  const int yTop = clamp(sourceY - boundingBoxSize, 0, int(imageHeight));
  const int yBottom = clamp(sourceY + boundingBoxSize, yTop, int(imageHeight));

  aocommon::UVector<double> values;
  double fluxSum = 0.0;
  for (int y = yTop; y != yBottom; ++y) {
    for (int x = xLeft; x != xRight; ++x) {
      long double l, m;
      ImageCoordinates::XYToLM<long double>(x, y, pixelScaleL, pixelScaleM,
                                            imageWidth, imageHeight, l, m);
      l += phaseCentreDL;
      m += phaseCentreDM;
      const long double lTransf =
          (l - sourceL) * transf[0] + (m - sourceM) * transf[1];
      const long double mTransf =
          (l - sourceL) * transf[2] + (m - sourceM) * transf[3];
      const long double dist = std::sqrt(lTransf * lTransf + mTransf * mTransf);
      long double v = Gaussian(dist, 1.0L) * flux;
      fluxSum += double(v);
      values.emplace_back(v);
    }
  }
  const double* iter = values.data();
  const double factor = flux / fluxSum;
  for (int y = yTop; y != yBottom; ++y) {
    float* imageDataPtr = imageData + y * imageWidth + xLeft;
    for (int x = xLeft; x != xRight; ++x) {
      (*imageDataPtr) += *iter * factor;
      ++imageDataPtr;
      ++iter;
    }
  }
}
}  // namespace

/** Restore a circular beam*/
void ModelRenderer::Restore(double* imageData, size_t imageWidth,
                            size_t imageHeight, const Model& model,
                            long double beamSize, long double startFrequency,
                            long double endFrequency,
                            aocommon::PolarizationEnum polarization) const {
  // Using the FWHM formula for a Gaussian:
  const long double sigma = beamSize / (2.0L * sqrtl(2.0L * logl(2.0L)));

  const int boundingBoxSize =
      std::ceil(sigma * 20.0 / std::min(_pixelScaleL, _pixelScaleM));
  for (const ModelSource& source : model) {
    for (const ModelComponent& component : source) {
      const long double posRA = component.PosRA();
      const long double posDec = component.PosDec();
      long double sourceL;
      long double sourceM;
      ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA,
                                  _phaseCentreDec, sourceL, sourceM);
      const SpectralEnergyDistribution& sed = component.SED();
      const long double intFlux =
          sed.IntegratedFlux(startFrequency, endFrequency, polarization);

      int sourceX, sourceY;
      ImageCoordinates::LMToXY<long double>(
          sourceL - _phaseCentreDL, sourceM - _phaseCentreDM, _pixelScaleL,
          _pixelScaleM, imageWidth, imageHeight, sourceX, sourceY);

      const int xLeft = clamp(sourceX - boundingBoxSize, 0, int(imageWidth));
      const int xRight =
          clamp(sourceX + boundingBoxSize, xLeft, int(imageWidth));
      const int yTop = clamp(sourceY - boundingBoxSize, 0, int(imageHeight));
      const int yBottom =
          clamp(sourceY + boundingBoxSize, yTop, int(imageHeight));

      for (int y = yTop; y != yBottom; ++y) {
        double* imageDataPtr = imageData + y * imageWidth + xLeft;
        for (int x = xLeft; x != xRight; ++x) {
          long double l, m;
          ImageCoordinates::XYToLM<long double>(
              x, y, _pixelScaleL, _pixelScaleM, imageWidth, imageHeight, l, m);
          l += _phaseCentreDL;
          m += _phaseCentreDM;
          const long double dist = std::sqrt((l - sourceL) * (l - sourceL) +
                                             (m - sourceM) * (m - sourceM));
          const long double g = Gaussian(dist, sigma);
          (*imageDataPtr) += double(g * intFlux);
          ++imageDataPtr;
        }
      }
    }
  }
}

void ModelRenderer::renderPointComponent(float* imageData, size_t imageWidth,
                                         size_t imageHeight, long double posRA,
                                         long double posDec,
                                         long double flux) const {
  long double sourceL, sourceM;
  ImageCoordinates::RaDecToLM(posRA, posDec, _phaseCentreRA, _phaseCentreDec,
                              sourceL, sourceM);
  sourceL -= _phaseCentreDL;
  sourceM -= _phaseCentreDM;

  int sourceX, sourceY;
  ImageCoordinates::LMToXY<long double>(sourceL, sourceM, _pixelScaleL,
                                        _pixelScaleM, imageWidth, imageHeight,
                                        sourceX, sourceY);

  if (sourceX >= 0 && sourceX < (int)imageWidth && sourceY >= 0 &&
      sourceY < (int)imageHeight) {
    float* imageDataPtr = imageData + sourceY * imageWidth + sourceX;
    (*imageDataPtr) += double(flux);
  }
}

/** Restore a model with an elliptical beam */
void ModelRenderer::Restore(float* imageData, size_t imageWidth,
                            size_t imageHeight, const Model& model,
                            long double beamMaj, long double beamMin,
                            long double beamPA, long double startFrequency,
                            long double endFrequency,
                            aocommon::PolarizationEnum polarization,
                            size_t threadCount) const {
  aocommon::UVector<float> renderedWithoutBeam(imageWidth * imageHeight, 0.0);
  RenderModel(renderedWithoutBeam.data(), imageWidth, imageHeight, model,
              startFrequency, endFrequency, polarization);
  // Restore(imageData, renderedWithoutBeam.data(), imageWidth, imageHeight,
  //         beamMaj, beamMin, beamPA, _pixelScaleL, _pixelScaleM, threadCount);
  schaapcommon::fft::RestoreImage(
      imageData, renderedWithoutBeam.data(), imageWidth, imageHeight, beamMaj,
      beamMin, beamPA, _pixelScaleL, _pixelScaleM, threadCount);
}

/**
 * Render without beam convolution, such that each point-source is one pixel.
 */
void ModelRenderer::RenderModel(float* imageData, size_t imageWidth,
                                size_t imageHeight, const Model& model,
                                long double startFrequency,
                                long double endFrequency,
                                aocommon::PolarizationEnum polarization) const {
  for (const ModelSource& source : model) {
    for (const ModelComponent& component : source) {
      const long double posRA = component.PosRA(), posDec = component.PosDec();
      const long double intFlux = component.SED().IntegratedFlux(
          startFrequency, endFrequency, polarization);

      if (component.Type() == ModelComponent::GaussianSource) {
        const long double gausMaj = component.MajorAxis(),
                          gausMin = component.MinorAxis();
        const long double gausPA = component.PositionAngle();
        RenderGaussianComponent(
            imageData, imageWidth, imageHeight, _phaseCentreRA, _phaseCentreDec,
            _pixelScaleL, _pixelScaleM, _phaseCentreDL, _phaseCentreDM, posRA,
            posDec, gausMaj, gausMin, gausPA, intFlux);
      } else
        renderPointComponent(imageData, imageWidth, imageHeight, posRA, posDec,
                             intFlux);
    }
  }
}
