#ifndef STRUCTURES_DDPSF_H_
#define STRUCTURES_DDPSF_H_

#include <vector>

#include <schaapcommon/facets/facet.h>

#include "../main/settings.h"
#include "../structures/observationinfo.h"

schaapcommon::facets::Facet::InitializationData CreateFacetInitilizationData(
    const Settings& settings, const ObservationInfo& observation_info);

std::vector<schaapcommon::facets::Facet> CreateRectangularPsfs(
    const Settings& settings, const ObservationInfo& observation_info);

#endif  // STRUCTURES_DDPSF_H_