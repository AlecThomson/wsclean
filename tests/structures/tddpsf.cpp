#include "../../structures/ddpsf.h"

#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(dd_psf)

BOOST_AUTO_TEST_CASE(create_single_psf) {
  const double phaseCentreRA = 0;
  const double phaseCentreDec = 0;
  const double pixelScaleX = 1;
  const double pixelScaleY = 1;

  const double psfsGridWidth = 1;
  const double psfsGridHeight = 1;
  size_t trimmedImageHeight = 100;
  size_t trimmedImageWidth = 100;

  const std::vector<std::shared_ptr<schaapcommon::facets::Facet>> ddpsfs =
      CreateRectangularPsfs(phaseCentreRA, phaseCentreDec, pixelScaleX,
                            pixelScaleY, trimmedImageHeight, trimmedImageWidth,
                            psfsGridWidth, psfsGridHeight);

  BOOST_CHECK_EQUAL(ddpsfs.size(), 1);
}

BOOST_AUTO_TEST_CASE(create_multiple_psfs) {
  const double phaseCentreRA = 0;
  const double phaseCentreDec = 0;
  const double pixelScaleX = 1;
  const double pixelScaleY = 1;

  const double psfsGridWidth = 5;
  const double psfsGridHeight = 5;
  size_t trimmedImageHeight = 100;
  size_t trimmedImageWidth = 100;

  const std::vector<std::shared_ptr<schaapcommon::facets::Facet>> ddpsfs =
      CreateRectangularPsfs(phaseCentreRA, phaseCentreDec, pixelScaleX,
                            pixelScaleY, trimmedImageHeight, trimmedImageWidth,
                            psfsGridWidth, psfsGridHeight);

  BOOST_REQUIRE_EQUAL(ddpsfs.size(), 25);
}

BOOST_AUTO_TEST_SUITE_END()