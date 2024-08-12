#include "gridder_simple.h"

#include <complex>

#include "ducc0/wgridder/wgridder.h"

#include "../gridding/gainmode.h"

using namespace ducc0;

namespace wsclean {

#include "callback_implementation.h"

template void WGriddingGridder_Simple<float>::AddInversionMs<GainMode::kFull>(
    size_t, size_t, const double*,
    std::reference_wrapper<const cmav<double, 1>>, size_t,
    std::reference_wrapper<const aocommon::BandData>,
    const std::pair<size_t, size_t>*, const std::complex<float>*, const size_t*,
    MsGridder*, size_t);

}  // namespace wsclean
