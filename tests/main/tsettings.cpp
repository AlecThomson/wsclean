#include <boost/test/unit_test.hpp>

#include "../../main/settings.h"

BOOST_AUTO_TEST_SUITE(settings)

BOOST_AUTO_TEST_CASE(get_feather_size_default) {
  Settings settings;
  settings.trimmedImageWidth = 10000;
  settings.trimmedImageHeight = 15625;
  BOOST_CHECK_EQUAL(settings.GetFeatherSize(),
                    125);  // 1% of sqrt(10000*15625)
}

BOOST_AUTO_TEST_CASE(get_feather_size_override) {
  Settings settings;
  settings.featherSize = 100;
  settings.trimmedImageWidth = 5000;
  settings.trimmedImageHeight = 7000;
  BOOST_CHECK_EQUAL(settings.GetFeatherSize(), 100);
}

BOOST_AUTO_TEST_SUITE_END()
