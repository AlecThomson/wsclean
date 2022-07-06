#include "../../structures/ddpsf.h"

#include <boost/test/unit_test.hpp>

namespace {
const double kPhaseCentreRa = 0;
const double kPhaseCentreDec = 0;
const double kPixelScaleX = 1;
const double kPixelScaleY = 1;
const size_t kTrimmedImageHeight = 100;
const size_t kTrimmedImageWidth = 100;
}  // namespace

BOOST_AUTO_TEST_SUITE(dd_psf)

BOOST_AUTO_TEST_CASE(create_single_psf) {
  const double kPsfGridWidth = 1;
  const double kPsfGridHeight = 1;

  const std::vector<schaapcommon::facets::Facet> dd_psfs =
      CreateRectangularPsfs(kPhaseCentreRa, kPhaseCentreDec, kPixelScaleX,
                            kPixelScaleY, kTrimmedImageHeight,
                            kTrimmedImageWidth, kPsfGridWidth, kPsfGridHeight);

  BOOST_CHECK_EQUAL(dd_psfs.size(), 1);
}

BOOST_AUTO_TEST_CASE(create_multiple_psfs) {
  const double kPsfGridWidth = 5;
  const double kPsfGridHeight = 5;

  const std::vector<schaapcommon::facets::Facet> dd_psfs =
      CreateRectangularPsfs(kPhaseCentreRa, kPhaseCentreDec, kPixelScaleX,
                            kPixelScaleY, kTrimmedImageHeight,
                            kTrimmedImageWidth, kPsfGridWidth, kPsfGridHeight);

  BOOST_REQUIRE_EQUAL(dd_psfs.size(), 25);
}

BOOST_AUTO_TEST_SUITE_END()