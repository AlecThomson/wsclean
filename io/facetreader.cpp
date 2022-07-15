
#include "facetreader.h"

#include <schaapcommon/facets/ds9facetfile.h>

#include "../structures/facetutil.h"

using schaapcommon::facets::DS9FacetFile;
using schaapcommon::facets::Facet;

std::vector<std::shared_ptr<Facet>> FacetReader::ReadFacets(
    std::string filename, double width, double height, double pixelScaleX,
    double pixelScaleY, double phaseCentreRA, double phaseCentreDec,
    double shiftL, double shiftM, double imagePadding, bool make_square) {
  Facet::InitializationData data(pixelScaleX, pixelScaleY, width, height);
  data.phase_centre.ra = phaseCentreRA;
  data.phase_centre.dec = phaseCentreDec;
  data.shift_l = shiftL;
  data.shift_m = shiftM;
  data.padding = imagePadding;
  data.align = 2;
  data.make_square = make_square;

  std::vector<std::shared_ptr<Facet>> facets;
  if (!filename.empty()) {
    facets = DS9FacetFile(filename).ReadShared(data);

    if (facets.empty()) {
      throw std::runtime_error("No facets found in " + filename);
    }
  }

  return facets;
}

std::size_t FacetReader::CountFacets(const std::string& filename) {
  return filename.empty() ? 0 : DS9FacetFile(filename).Count();
}
