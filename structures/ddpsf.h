#ifndef STRUCTURES_DDPSF_H_
#define STRUCTURES_DDPSF_H_

#include <vector>

#include <schaapcommon/facets/facet.h>

std::vector<schaapcommon::facets::Facet> CreateRectangularPsfs(
    double phase_centre_ra, double phase_centre_dec, double pixel_scale_x,
    double pixel_scale_y, size_t trimmed_image_height,
    size_t trimmed_image_width, size_t psf_grid_width, size_t psf_grid_height);

#endif  // STRUCTURES_DDPSF_H_