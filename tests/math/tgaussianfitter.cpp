#include <boost/test/unit_test.hpp>

#include "../../math/gaussianfitter.h"
#include "../../math/renderer.h"

#include "../../model/model.h"
#include "../../model/powerlawsed.h"

#include <aocommon/image.h>

#include <schaapcommon/fft/restoreimage.h>

BOOST_AUTO_TEST_SUITE(gaussian_fitter)

BOOST_AUTO_TEST_CASE(fit) {
  for (size_t beamPAindex = 0; beamPAindex != 10; ++beamPAindex) {
    const size_t width = 512, height = 512;
    aocommon::Image model(width, height, 0.0);
    aocommon::Image restored(width, height, 0.0);
    model[((height / 2) * width) + (width / 2)] = 1.0;
    const long double pixelScale = 1.0L /*amin*/ * (M_PI / 180.0 / 60.0);
    const long double beamMaj = 20.0L * pixelScale;
    const long double beamMin = 5.0L * pixelScale;
    const long double beamPA = beamPAindex * M_PI / 10.0;
    const size_t threadCount = 1;
    schaapcommon::fft::RestoreImage(restored.Data(), model.Data(), width,
                                    height, beamMaj, beamMin, beamPA,
                                    pixelScale, pixelScale, threadCount);

    GaussianFitter fitter;
    double fitMaj, fitMin, fitPA;
    fitter.Fit2DGaussianCentred(restored.Data(), width, height, 5.0, fitMaj,
                                fitMin, fitPA, 10.0, false);
    fitPA = fmod((fitPA + 2.0 * M_PI), M_PI);

    BOOST_CHECK_CLOSE_FRACTION(fitMaj, 20.0, 1e-3);
    BOOST_CHECK_CLOSE_FRACTION(fitMin, 5.0, 1e-3);
    BOOST_CHECK_CLOSE_FRACTION(fitPA, beamPA, 1e-3);
  }
}

BOOST_AUTO_TEST_CASE(fit_with_bad_initial_value) {
  const size_t threadCount = 1;
  const size_t width = 64;
  const size_t height = 64;
  aocommon::Image restored(width, height, 0.0);
  PowerLawSED sed(150.0e6, 1.0);
  ModelComponent component;
  component.SetPosDec(0.0);
  component.SetPosRA(0.0);
  component.SetSED(sed);
  ModelSource source;
  source.AddComponent(component);
  Model model;
  model.AddSource(source);
  // long double pixelScale = 1 /*amin*/ * (M_PI / 180.0 / 60.0);
  // const long double phaseCentreDL = 0.0;
  // const long double phaseCentreDM = 0.0;

  renderer::ImageCoordinateSettings imageSettings;
  imageSettings.phase_centre_ra = 0.0;
  imageSettings.phase_centre_dec = 0.0;
  imageSettings.pixel_scale_l = 1 /*amin*/ * (M_PI / 180.0 / 60.0);
  imageSettings.pixel_scale_m = imageSettings.pixel_scale_l;
  imageSettings.phase_centre_dl = 0.0;
  imageSettings.phase_centre_dm = 0.0;

  const long double beamMaj = 4.0L * imageSettings.pixel_scale_l;
  const long double beamMin = 4.0L * imageSettings.pixel_scale_l;
  const long double beamPA = 0.0;
  const long double estimatedBeamPx = 1.0;  // this is on purpose way off

  // Renderer renderer(model, 0.0, 0.0, pixelScale, pixelScale, phaseCentreDL,
  //                   phaseCentreDM);
  // renderer.RestoreWithEllipticalBeam(
  //     restored.Data(), width, height, beamMaj, beamMin, beamPA, 100e6, 200e6,
  //     aocommon::Polarization::StokesI, threadCount);
  renderer::RestoreWithEllipticalBeam(
      restored, imageSettings, model, beamMaj, beamMin, beamPA, 100e6, 200e6,
      aocommon::Polarization::StokesI, threadCount);

  GaussianFitter fitter;
  double fitMajor, fitMinor, fitPA;
  fitter.Fit2DGaussianCentred(restored.Data(), width, height, estimatedBeamPx,
                              fitMajor, fitMinor, fitPA, 10.0, false);

  BOOST_CHECK_CLOSE_FRACTION(fitMajor, 4.0, 1e-4);
  BOOST_CHECK_CLOSE_FRACTION(fitMinor, 4.0, 1e-4);
}

BOOST_AUTO_TEST_CASE(fit_circular) {
  const size_t width = 64;
  const size_t height = 64;
  const size_t threadCount = 1;
  aocommon::Image restored(width, height, 0.0);
  PowerLawSED sed(150.0e6, 1.0);
  ModelComponent component;
  component.SetPosDec(0.0);
  component.SetPosRA(0.0);
  component.SetSED(sed);
  ModelSource source;
  source.AddComponent(component);
  Model model;
  model.AddSource(source);

  renderer::ImageCoordinateSettings imageSettings;
  imageSettings.phase_centre_ra = 0.0;
  imageSettings.phase_centre_dec = 0.0;
  imageSettings.pixel_scale_l = 1 /*amin*/ * (M_PI / 180.0 / 60.0);
  imageSettings.pixel_scale_m = imageSettings.pixel_scale_l;
  imageSettings.phase_centre_dl = 0.0;
  imageSettings.phase_centre_dm = 0.0;

  const long double beamMaj = 4.0L * imageSettings.pixel_scale_l;
  const long double beamMin = 4.0L * imageSettings.pixel_scale_l;
  const long double beamPA = 0.0;
  const long double estimatedBeamPx = 1.0;  // this is on purpose way off

  // const long double pixelScale = 1 /*amin*/ * (M_PI / 180.0 / 60.0);
  // const long double beamMaj = 4 * pixelScale;
  // const long double beamMin = 4 * pixelScale;
  // const long double beamPA = 0.0;
  // const long double estimatedBeamPx = 1.0;  // this is on purpose way off
  // const long double phaseCentreDL = 0.0;
  // const long double phaseCentreDM = 0.0;
  // Renderer renderer(model, 0.0, 0.0, pixelScale, pixelScale, phaseCentreDL,
  //                   phaseCentreDM);
  // renderer.RestoreWithEllipticalBeam(
  //     restored.Data(), width, height, beamMaj, beamMin, beamPA, 100e6, 200e6,
  //     aocommon::Polarization::StokesI, threadCount);

  renderer::RestoreWithEllipticalBeam(
      restored, imageSettings, model, beamMaj, beamMin, beamPA, 100e6, 200e6,
      aocommon::Polarization::StokesI, threadCount);

  GaussianFitter fitter;
  double fitMajor = estimatedBeamPx;
  fitter.Fit2DCircularGaussianCentred(restored.Data(), width, height, fitMajor);

  BOOST_CHECK_CLOSE_FRACTION(fitMajor, 4.0, 1e-4);
}

BOOST_AUTO_TEST_CASE(fit_small_beam) {
  const size_t width = 64;
  const size_t height = 64;
  const size_t threadCount = 1;
  aocommon::Image restored(width, height, 0.0);
  PowerLawSED sed(150.0e6, 1.0);
  ModelComponent component;
  component.SetPosDec(0.0);
  component.SetPosRA(0.0);
  component.SetSED(sed);
  ModelSource source;
  source.AddComponent(component);
  Model model;
  model.AddSource(source);
  // const long double pixelScale = 1 /*amin*/ * (M_PI / 180.0 / 60.0);
  // const long double beamMaj = 4 * pixelScale, beamMin = 0.5 * pixelScale;
  // const long double beamPA = 0.0;
  // const long double estimatedBeamPx = 1.0;  // this is on purpose way off
  // const long double phaseCentreDL = 0.0;
  // const long double phaseCentreDM = 0.0;

  renderer::ImageCoordinateSettings imageSettings;
  imageSettings.phase_centre_ra = 0.0;
  imageSettings.phase_centre_dec = 0.0;
  imageSettings.pixel_scale_l = 1 /*amin*/ * (M_PI / 180.0 / 60.0);
  imageSettings.pixel_scale_m = imageSettings.pixel_scale_l;
  imageSettings.phase_centre_dl = 0.0;
  imageSettings.phase_centre_dm = 0.0;

  const long double beamMaj = 4.0L * imageSettings.pixel_scale_l;
  const long double beamMin = 0.5L * imageSettings.pixel_scale_l;
  const long double beamPA = 0.0;
  const long double estimatedBeamPx = 1.0;  // this is on purpose way off

  renderer::RestoreWithEllipticalBeam(
      restored, imageSettings, model, beamMaj, beamMin, beamPA, 100e6, 200e6,
      aocommon::Polarization::StokesI, threadCount);

  // Renderer renderer(model, 0.0, 0.0, pixelScale, pixelScale, phaseCentreDL,
  //                   phaseCentreDM);
  // renderer.RestoreWithEllipticalBeam(
  //     restored.Data(), width, height, beamMaj, beamMin, beamPA, 100e6, 200e6,
  //     aocommon::Polarization::StokesI, threadCount);

  GaussianFitter fitter;
  double fitMajor = estimatedBeamPx, fitMinor = estimatedBeamPx, fitPA = 0.0;
  fitter.Fit2DGaussianCentred(restored.Data(), width, height, estimatedBeamPx,
                              fitMajor, fitMinor, fitPA, 10.0, false);

  BOOST_CHECK_CLOSE_FRACTION(fitMinor, 0.5, 1e-4);
}

BOOST_AUTO_TEST_SUITE_END()
