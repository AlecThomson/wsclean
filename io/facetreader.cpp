
#include "facetreader.h"

#include <schaapcommon/facets/ds9facetfile.h>

std::vector<schaapcommon::facets::Facet> FacetReader::ReadFacets(
    const Settings& settings, double phaseCentreRA, double phaseCentreDec) {
  std::vector<schaapcommon::facets::Facet> facets;
  if (!settings.facetRegionFilename.empty()) {
    facets = schaapcommon::facets::DS9FacetFile(settings.facetRegionFilename)
                 .Read(phaseCentreRA, phaseCentreDec, settings.pixelScaleX,
                       settings.pixelScaleY, settings.trimmedImageWidth,
                       settings.trimmedImageHeight);

    if (facets.empty())
      throw std::runtime_error("No facets found in " +
                               settings.facetRegionFilename);
  }

  return facets;
}
