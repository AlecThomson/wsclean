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
               aocommon::PolarizationEnum polarization) const;

  /**
   * Restore a model with an elliptical beam
   */
  void Restore(float* imageData, size_t imageWidth, size_t imageHeight,
               const Model& model, long double beamMaj, long double beamMin,
               long double beamPA, long double startFrequency,
               long double endFrequency,
               aocommon::PolarizationEnum polarization,
               size_t threadCount) const;

  // /**
  //  * Restore elliptical beam using a FFT deconvolution
  //  */
  // void Restore(float* imageData, const float* modelData, size_t imageWidth,
  //              size_t imageHeight, long double beamMaj, long double beamMin,
  //              long double beamPA, size_t threadCount = 1) {
  //   Restore(imageData, modelData, imageWidth, imageHeight, beamMaj, beamMin,
  //           beamPA, _pixelScaleL, _pixelScaleM, threadCount);
  // }

 private:
  /**
   * Render without beam convolution, such that each point-source is one pixel.
   */
  void RenderModel(float* imageData, size_t imageWidth, size_t imageHeight,
                   const Model& model, long double startFrequency,
                   long double endFrequency,
                   aocommon::PolarizationEnum polarization) const;

  void renderPointComponent(float* imageData, size_t imageWidth,
                            size_t imageHeight, long double posRA,
                            long double posDec, long double flux) const;

  const long double _phaseCentreRA;
  const long double _phaseCentreDec;
  const long double _pixelScaleL;
  const long double _pixelScaleM;
  const long double _phaseCentreDL;
  const long double _phaseCentreDM;
};

// class Renderer {
//  public:
//   Renderer(const Model& model, long double phaseCentreRA,
//            long double phaseCentreDec, long double pixelScaleL,
//            long double pixelScaleM, long double phaseCentreDL,
//            long double phaseCentreDM)
//       : model_(model),
//         _phaseCentreRA(phaseCentreRA),
//         _phaseCentreDec(phaseCentreDec),
//         _pixelScaleL(pixelScaleL),
//         _pixelScaleM(pixelScaleM),
//         _phaseCentreDL(phaseCentreDL),
//         _phaseCentreDM(phaseCentreDM) {}

//   Renderer(const Renderer&) = delete;
//   Renderer& operator=(const Renderer&) = delete;

//   /**
//    * Restore model component using a circular beam, goverened by
//    * \c beamSize .
//    */
//   void RestoreWithCircularBeam(double* imageData, size_t imageWidth,
//                                size_t imageHeight, const Model& model,
//                                long double beamSize, long double
//                                startFrequency, long double endFrequency,
//                                aocommon::PolarizationEnum polarization)
//                                const;

//   /**
//    * Restore a model with an elliptical beam
//    */
//   void RestoreWithEllipticalBeam(float* imageData, size_t imageWidth,
//                                  size_t imageHeight, const Model& model,
//                                  long double beamMaj, long double beamMin,
//                                  long double beamPA, long double
//                                  startFrequency, long double endFrequency,
//                                  aocommon::PolarizationEnum polarization,
//                                  size_t threadCount) const;

//  private:
//   /**
//    * Render without beam convolution, such that each point-source is one
//    pixel.
//    */
//   void RenderModel(float* imageData, size_t imageWidth, size_t imageHeight,
//                    const Model& model, long double startFrequency,
//                    long double endFrequency,
//                    aocommon::PolarizationEnum polarization) const;

//   void renderPointComponent(float* imageData, size_t imageWidth,
//                             size_t imageHeight, long double posRA,
//                             long double posDec, long double flux) const;

//   const Model& model_;
//   const long double _phaseCentreRA;
//   const long double _phaseCentreDec;
//   const long double _pixelScaleL;
//   const long double _pixelScaleM;
//   const long double _phaseCentreDL;
//   const long double _phaseCentreDM;
// };

#endif
