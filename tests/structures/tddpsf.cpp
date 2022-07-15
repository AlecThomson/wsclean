#include "../../structures/ddpsf.h"

#include <boost/test/unit_test.hpp>

namespace {

struct SettingsFixture {
  SettingsFixture() {
    settings.pixelScaleX = 0.01;
    settings.pixelScaleY = 0.01;
    settings.trimmedImageWidth = 100;
    settings.trimmedImageHeight = 100;
  }

  Settings settings;
  ObservationInfo observation_info;  // Use default values only.
};

}  // namespace

BOOST_AUTO_TEST_SUITE(dd_psf)

BOOST_FIXTURE_TEST_CASE(create_single_psf, SettingsFixture) {
  settings.psfsGridWidth = 1;
  settings.psfsGridHeight = 1;

  const std::vector<schaapcommon::facets::Facet> dd_psfs =
      CreateRectangularPsfs(settings, observation_info);

  BOOST_CHECK_EQUAL(dd_psfs.size(), 1);
}

BOOST_FIXTURE_TEST_CASE(create_multiple_psfs, SettingsFixture) {
  settings.psfsGridWidth = 5;
  settings.psfsGridHeight = 5;

  const std::vector<schaapcommon::facets::Facet> dd_psfs =
      CreateRectangularPsfs(settings, observation_info);

  BOOST_REQUIRE_EQUAL(dd_psfs.size(), 25);
}

BOOST_AUTO_TEST_SUITE_END()