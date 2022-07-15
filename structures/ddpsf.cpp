#include "ddpsf.h"

#include <aocommon/imagecoordinates.h>

using schaapcommon::facets::Facet;

Facet::InitializationData CreateFacetInitilizationData(
    const Settings& settings, const ObservationInfo& observation_info) {
  Facet::InitializationData data(settings.pixelScaleX, settings.pixelScaleY,
                                 settings.trimmedImageWidth,
                                 settings.trimmedImageHeight);
  data.phase_centre.ra = observation_info.phaseCentreRA;
  data.phase_centre.dec = observation_info.phaseCentreDec;
  data.shift_l = observation_info.shiftL;
  data.shift_m = observation_info.shiftM;
  data.padding = settings.imagePadding;
  data.align = 2;
  data.make_square = settings.gridderType == GridderType::IDG;
  return data;
}

std::vector<Facet> CreateRectangularPsfs(
    const Settings& settings, const ObservationInfo& observation_info) {
  const Facet::InitializationData facet_data =
      CreateFacetInitilizationData(settings, observation_info);
  const size_t psf_grid_width = settings.psfsGridWidth;
  const size_t psf_grid_height = settings.psfsGridHeight;

  std::vector<Facet> ddpsfs;
  ddpsfs.reserve(psf_grid_height * psf_grid_width);

  double single_psf_height = settings.trimmedImageHeight / psf_grid_height;
  double single_psf_width = settings.trimmedImageWidth / psf_grid_width;

  for (int grid_y = 0; grid_y < static_cast<int>(psf_grid_height); ++grid_y) {
    for (int grid_x = 0; grid_x < static_cast<int>(psf_grid_width); ++grid_x) {
      const int facet_start_x = grid_x * single_psf_width;
      const int facet_start_y = grid_y * single_psf_height;
      const int facet_end_x = (grid_x + 1) * single_psf_width;
      const int facet_end_y = (grid_y + 1) * single_psf_height;
      const schaapcommon::facets::BoundingBox box(
          {{facet_start_x, facet_start_y}, {facet_end_x, facet_end_y}});

      ddpsfs.emplace_back(facet_data, box);
      // add a name label for this box
      ddpsfs.back().SetDirectionLabel(std::to_string(grid_x) + ", " +
                                      std::to_string(grid_y));
    }
  }

  return ddpsfs;
}