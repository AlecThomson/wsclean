#include "../../main/commandline.h"
#include "../../main/wsclean.h"

#include <boost/test/unit_test.hpp>

namespace {
const char* kMWA_MS = "test_data/MWA_MOCK.ms/";
const char* kFacets = "test_data/facets.reg";
}  // namespace

BOOST_AUTO_TEST_SUITE(facet_stitching)

BOOST_AUTO_TEST_CASE(wstacking) {
  BOOST_REQUIRE(boost::filesystem::is_directory(kMWA_MS));
  BOOST_REQUIRE(boost::filesystem::is_regular_file(kFacets));

  WSClean wsclean;
  CommandLine commandLine;
  std::vector<const char*> args = {
      "wsclean", "-quiet", "-size",        "256",   "256",
      "-scale",  "4amin",  "-pol",         "XX,YY", "-facet-regions",
      kFacets,   "-name",  "facet-stitch", kMWA_MS};
  commandLine.Parse(wsclean, args.size(), args.data(), false);
  commandLine.Run(wsclean);

  std::vector<std::string> filesExpected;

  filesExpected.emplace_back("facet-stitch-XX-dirty.fits");
  filesExpected.emplace_back("facet-stitch-YY-dirty.fits");
  filesExpected.emplace_back("facet-stitch-XX-image.fits");
  filesExpected.emplace_back("facet-stitch-YY-image.fits");

  for (const std::string& str : filesExpected)
    BOOST_CHECK(boost::filesystem::is_regular_file(str));

  for (const std::string& str : filesExpected)
    BOOST_CHECK(boost::filesystem::remove(str));

  BOOST_CHECK(boost::filesystem::remove(kFacets));
}

BOOST_AUTO_TEST_SUITE_END()
