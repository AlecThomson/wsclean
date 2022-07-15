#include "../../structures/facetutil.h"

#include <boost/test/unit_test.hpp>

using schaapcommon::facets::Facet;

BOOST_AUTO_TEST_SUITE(dd_psf)

BOOST_AUTO_TEST_CASE(create_single_psf) {
  Facet::InitializationData facet_data(0.01, 100);

  const std::vector<schaapcommon::facets::Facet> dd_psfs =
      CreateFacetGrid(facet_data, 1, 1);

  BOOST_CHECK_EQUAL(dd_psfs.size(), 1);
}

BOOST_AUTO_TEST_CASE(create_multiple_psfs) {
  Facet::InitializationData facet_data(0.01, 100);

  const std::vector<schaapcommon::facets::Facet> dd_psfs =
      CreateFacetGrid(facet_data, 5, 5);

  BOOST_REQUIRE_EQUAL(dd_psfs.size(), 25);
}

BOOST_AUTO_TEST_SUITE_END()