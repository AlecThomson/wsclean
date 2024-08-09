#include "msgridderdata.h"

#include "msprovidercollection.h"

namespace wsclean {

MsGridderData::MsGridderData(const Settings& settings)
    : settings_(settings),
      visibility_weighting_mode_(settings.visibilityWeightingMode) {
  settings_ = settings;
#ifdef HAVE_EVERYBEAM
  visibility_modifier_.SetBeamInfo(settings.beamMode,
                                   settings.beamNormalisationMode);
#endif
}

void MsGridderData::StartMeasurementSet(
    size_t ms_count, const MsProviderCollection::MsData& ms_data,
    bool is_predict) {
  InitializePointResponse(ms_data);
  if (visibility_modifier_.HasH5Parm()) {
    visibility_modifier_.InitializeCacheParmResponse(ms_data.antenna_names,
                                                     ms_data.SelectedBand(),
                                                     ms_data.original_ms_index);
  }
  original_ms_index_ = ms_data.original_ms_index;
  writer_lock_index_ =
      facet_group_index_ * ms_count + ms_data.original_ms_index;
  n_vis_polarizations_ = ms_data.ms_provider->NPolarizations();
  gain_mode_ = GetGainMode(Polarization(), n_vis_polarizations_);
  const size_t n_channels = ms_data.SelectedBand().ChannelCount();
  scratch_image_weights_.resize(n_channels);
  if (is_predict) {
    scratch_model_data_.resize(n_channels *
                               ms_data.ms_provider->NPolarizations());
    predict_reader_ = ms_data.ms_provider->MakeReader();
  }
}

void MsGridderData::ReadPredictMetaData(MSProvider::MetaData& metadata) {
  predict_reader_->ReadMeta(metadata);
  predict_reader_->NextInputRow();
}

void MsGridderData::ResetVisibilityModifierCache(size_t ms_count) {
  if (visibility_modifier_.HasH5Parm()) {
    visibility_modifier_.ResetCache(ms_count);
  }
  visibility_modifier_.ResetSums();
}

void MsGridderData::InitializePointResponse(
    const MsProviderCollection::MsData& ms_data) {
#ifdef HAVE_EVERYBEAM
  if (settings_.applyFacetBeam || settings_.gridWithBeam) {
    const std::string element_response_string =
        !settings_.beamModel.empty() ? settings_.beamModel : "DEFAULT";
    visibility_modifier_.InitializePointResponse(
        ms_data.ms_provider->MS(), settings_.facetBeamUpdateTime,
        element_response_string, ms_data.band_data.ChannelCount(),
        settings_.dataColumnName, settings_.mwaPath);
  } else {
    visibility_modifier_.SetNoPointResponse();
  }
#else
  if (settings_.applyFacetBeam || settings_.gridWithBeam) {
    throw std::runtime_error(
        "-apply-facet-beam or -grid-with-beam was set, but wsclean was not "
        "compiled "
        "with EveryBeam. Please compile wsclean with EveryBeam to "
        "use the Facet Beam functionality");
  }
#endif
}

}  // namespace wsclean
