#ifndef IDG_MS_GRIDDER_H
#define IDG_MS_GRIDDER_H

#ifdef HAVE_IDG

#include "../gridding/msgridderbase.h"
#include "../structures/resources.h"

#include <idg-api.h>

#include <aocommon/image.h>
#include <aocommon/lane.h>
#include <aocommon/uvector.h>
#include <aocommon/fits/fitswriter.h>

#ifdef HAVE_EVERYBEAM
#include <EveryBeam/aterms/atermbase.h>
#endif  // HAVE_EVERYBEAM

#include "../main/stopwatch.h"

struct ImagingTableEntry;
class ImageFilename;
class Settings;

class IdgMsGridder final : public MSGridderBase {
 public:
  IdgMsGridder(const Settings& settings, const Resources& resources,
               const MSManager& measurement_sets, const size_t gridder_index);

  virtual ~IdgMsGridder() final override;

  virtual void Invert() final override;

  virtual void Predict(std::vector<aocommon::Image>&& images) final override;

  virtual std::vector<aocommon::Image> ResultImages() final override;

  static void SavePBCorrectedImages(aocommon::FitsWriter& writer,
                                    const ImageFilename& filename,
                                    const std::string& filenameKind,
                                    const Settings& settings);

  static void SaveBeamImage(const ImagingTableEntry& entry,
                            ImageFilename& filename, const Settings& settings,
                            double ra, double dec, double pdl, double pdm,
                            const AverageBeam& average_beam);

  void SetAverageBeam(std::unique_ptr<AverageBeam> average_beam);

  std::unique_ptr<AverageBeam> ReleaseAverageBeam();

 private:
  std::unique_ptr<AverageBeam> _averageBeam;

  virtual size_t getSuggestedWGridSize() const override final {
    return 1;  // TODO
  }

  void gridMeasurementSet(const MSManager::Data& msData,
                          const MSManager::FacetData& ms_facet_data);
  void gridThreadFunction();

  void predictMeasurementSet(const MSManager::Data& msData,
                             const MSManager::FacetData& ms_facet_data);
  void readConfiguration();

  void setIdgType();

#ifdef HAVE_EVERYBEAM
  std::unique_ptr<class everybeam::aterms::ATermBase> getATermMaker(
      const MSManager::Data& msData);
  bool prepareForMeasurementSet(
      const MSManager::Data& ms_data, const MSManager::FacetData& ms_facet_data,
      std::unique_ptr<everybeam::aterms::ATermBase>& aTermMaker,
      aocommon::UVector<std::complex<float>>& aTermBuffer,
      idg::api::BufferSetType);
#else
  bool prepareForMeasurementSet(
      const MSManager::Data& ms_data, const MSManager::FacetData& ms_facet_data,
      aocommon::UVector<std::complex<float>>& aTermBuffer,
      idg::api::BufferSetType);
#endif  // HAVE_EVERYBEAM

  struct IDGInversionRow : public MSGridderBase::InversionRow {
    size_t antenna1, antenna2, timeIndex;
  };
  struct IDGPredictionRow {
    double uvw[3];
    size_t antenna1, antenna2, timeIndex, rowId;
  };
  void predictRow(IDGPredictionRow& row,
                  const std::vector<std::string>& antennaNames);
  void computePredictionBuffer(const std::vector<std::string>& antennaNames);

  std::unique_ptr<idg::api::BufferSet> _bufferset;
  size_t _subgridSize;
  aocommon::UVector<double> _image;
  aocommon::UVector<float> _taper_subgrid;
  aocommon::UVector<float> _taper_grid;
  MSProvider* _outputProvider;
  aocommon::BandData _selectedBand;
  idg::api::Type _proxyType;
  int _buffersize;
  idg::api::options_type _options;
  Stopwatch _griddingWatch;
  Stopwatch _degriddingWatch;
  const Resources _resources;
};

void init_optimal_taper_1D(int subgridsize, int gridsize, float kernelsize,
                           float padding, float* taper_subgrid,
                           float* taper_grid);
void init_optimal_gridding_taper_1D(int subgridsize, int gridsize,
                                    float kernelsize, float* taper_subgrid,
                                    float* taper_grid);

#else

#include "../gridding/unavailablegridder.h"

#define IdgMsGridder UnavailableGridder

#endif  // HAVE IDG

#endif
