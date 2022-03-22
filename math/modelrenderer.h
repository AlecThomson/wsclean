#ifndef MODELRENDERER_H
#define MODELRENDERER_H

#include <cstring>

#include <aocommon/polarization.h>

#include "../model/model.h"

class ModelRenderer {
 public:
  ModelRenderer(long double phaseCentreRA, long double phaseCentreDec,
                long double pixelScaleL, long double pixelScaleM,
                long double phaseCentreDL, long double phaseCentreDM)
      : _phaseCentreRA(phaseCentreRA),
        _phaseCentreDec(phaseCentreDec),
        _pixelScaleL(pixelScaleL),
        _pixelScaleM(pixelScaleM),
        _phaseCentreDL(phaseCentreDL),
        _phaseCentreDM(phaseCentreDM) {}

  ModelRenderer(const ModelRenderer&) = delete;
  ModelRenderer& operator=(const ModelRenderer&) = delete;

  /**
   * Restore with circular beam
   */
  void Restore(double* imageData, size_t imageWidth, size_t imageHeight,
               const Model& model, long double beamSize,
               long double startFrequency, long double endFrequency,
               aocommon::PolarizationEnum polarization);

  /**
   * Restore a model with an elliptical beam
   */
  void Restore(float* imageData, size_t imageWidth, size_t imageHeight,
               const Model& model, long double beamMaj, long double beamMin,
               long double beamPA, long double startFrequency,
               long double endFrequency,
               aocommon::PolarizationEnum polarization, size_t threadCount);

  /**
   * Restore elliptical beam using a FFT deconvolution
   */
  void Restore(float* imageData, const float* modelData, size_t imageWidth,
               size_t imageHeight, long double beamMaj, long double beamMin,
               long double beamPA, size_t threadCount = 1) {
    Restore(imageData, modelData, imageWidth, imageHeight, beamMaj, beamMin,
            beamPA, _pixelScaleL, _pixelScaleM, threadCount);
  }

  /**
   * Restore elliptical beam using a FFT deconvolution (static version).
   */
  static void Restore(float* imageData, const float* modelData,
                      size_t imageWidth, size_t imageHeight,
                      long double beamMaj, long double beamMin,
                      long double beamPA, long double pixelScaleL,
                      long double pixelScaleM, size_t threadCount);

  /**
   * Render without beam convolution, such that each point-source is one pixel.
   */
  void RenderModel(float* imageData, size_t imageWidth, size_t imageHeight,
                   const Model& model, long double startFrequency,
                   long double endFrequency,
                   aocommon::PolarizationEnum polarization);

  static void RenderGaussianComponent(
      float* imageData, size_t imageWidth, size_t imageHeight,
      long double phaseCentreRA, long double phaseCentreDec,
      long double pixelScaleL, long double pixelScaleM,
      long double phaseCentreDL, long double phaseCentreDM, long double posRA,
      long double posDec, long double gausMaj, long double gausMin,
      long double gausPA, long double flux);

 private:
  void renderPointComponent(float* imageData, size_t imageWidth,
                            size_t imageHeight, long double posRA,
                            long double posDec, long double flux);

  const long double _phaseCentreRA;
  const long double _phaseCentreDec;
  const long double _pixelScaleL;
  const long double _pixelScaleM;
  const long double _phaseCentreDL;
  const long double _phaseCentreDM;
};

#endif
