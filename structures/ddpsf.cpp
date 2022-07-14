#include "ddpsf.h"

#include <aocommon/imagecoordinates.h>

using schaapcommon::facets::Facet;

std::vector<Facet> CreateRectangularPsfs(
    const double phase_centre_ra, const double phase_centre_dec,
    const double pixel_scale_x, const double pixel_scale_y,
    const size_t trimmed_image_height, const size_t trimmed_image_width,
    const size_t psf_grid_width, const size_t psf_grid_height) {
  Facet::InitializationData facet_data(
      pixel_scale_x, pixel_scale_y, trimmed_image_width, trimmed_image_height);
  facet_data.phase_centre.ra = phase_centre_ra;
  facet_data.phase_centre.dec = phase_centre_dec;

  std::vector<Facet> ddpsfs;
  ddpsfs.reserve(psf_grid_height * psf_grid_width);

  double single_psf_height = trimmed_image_height / psf_grid_height;
  double single_psf_width = trimmed_image_width / psf_grid_width;

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