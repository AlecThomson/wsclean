#include "ddpsf.h"

#include <aocommon/imagecoordinates.h>

std::vector<schaapcommon::facets::Facet> CreateRectangularPsfs(
    const double phase_centre_ra, const double phase_centre_dec,
    const double pixelScaleX, const double pixel_scale_y,
    const size_t trimmed_image_height, const size_t trimmed_image_width,
    const size_t psf_grid_width, const size_t psf_grid_height) {
  std::vector<schaapcommon::facets::Facet> ddpsfs;
  ddpsfs.reserve(psf_grid_height * psf_grid_width);

  double single_psf_height = trimmed_image_height / psf_grid_height;
  double single_psf_width = trimmed_image_width / psf_grid_width;

  for (size_t grid_y = 0; grid_y < psf_grid_height; ++grid_y) {
    for (size_t grid_x = 0; grid_x < psf_grid_width; ++grid_x) {
      schaapcommon::facets::Facet facet;

      const auto add_vertex = [&](size_t x, size_t y) {
        double l = 0.0;
        double m = 0.0;
        double ra = 0.0;
        double dec = 0.0;
        aocommon::ImageCoordinates::XYToLM(x, y, pixelScaleX, pixel_scale_y,
                                           trimmed_image_width,
                                           trimmed_image_height, l, m);
        aocommon::ImageCoordinates::LMToRaDec(l, m, phase_centre_ra,
                                              phase_centre_dec, ra, dec);
        facet.AddVertex(ra, dec);
      };

      const size_t facet_start_x = grid_x * single_psf_width;
      const size_t facet_start_y = grid_y * single_psf_height;
      const size_t facet_end_x = (grid_x + 1) * single_psf_width;
      const size_t facet_end_y = (grid_y + 1) * single_psf_height;
      add_vertex(facet_start_x, facet_start_y);
      add_vertex(facet_start_x, facet_end_y);
      add_vertex(facet_end_x, facet_end_y);
      add_vertex(facet_end_x, facet_start_y);

      // add a name label for this box
      facet.SetDirectionLabel(std::to_string(grid_x) + ", " +
                              std::to_string(grid_y));

      ddpsfs.push_back(std::move(facet));
    }
  }

  return ddpsfs;
}