#include "directmsgridder.h"

#include "msgriddermanager.h"

#include "../main/progressbar.h"
#include "../msproviders/msprovider.h"
#include "../msproviders/msreaders/msreader.h"

#include <thread>
#include <vector>

#include <aocommon/threadpool.h>

namespace wsclean {

template <typename num_t>
DirectMSGridder<num_t>::DirectMSGridder(
    const Settings& settings, const Resources& resources,
    MsProviderCollection& ms_provider_collection)
    : MSGridderBase(settings, ms_provider_collection), _resources(resources) {}

template <typename num_t>
void DirectMSGridder<num_t>::Invert() {
  _sqrtLMTable = GetSqrtLMLookupTable();
  const size_t width = TrimWidth(), height = TrimHeight();

  resetVisibilityCounters();

  ProgressBar progress("Performing direct Fourier transform");

  _inversionLane.resize(_resources.NCpus() * 1024);

  for (size_t t = 0; t != _resources.NCpus(); ++t) {
    _layers.emplace_back(aocommon::ImageBase<num_t>(ImageWidth(), ImageHeight(),
                                                    static_cast<num_t>(0.0)));
  }

  aocommon::ThreadPool& thread_pool = aocommon::ThreadPool::GetInstance();
  thread_pool.SetNThreads(_resources.NCpus() + 1);
  thread_pool.StartParallelExecution([&](size_t thread_index) {
    const size_t layer = thread_index - 1;
    InversionSample sample;
    while (_inversionLane.read(sample)) {
      gridSample(sample, layer);
    }
  });

  for (size_t i = 0; i != GetMsCount(); ++i) {
    invertMeasurementSet(GetMsData(i), progress, i);
  }

  _inversionLane.write_end();
  thread_pool.FinishParallelExecution();

  aocommon::ImageBase<num_t> scratch = std::move(_layers.back());
  _layers.pop_back();

  for (const aocommon::ImageBase<num_t>& layer : _layers) {
    scratch += layer;
  }

  _layers.clear();
  _sqrtLMTable.Reset();

  // Wrap the image correctly and normalize it
  _image = aocommon::Image(TrimWidth(), TrimHeight());
  const double weight_factor = 1.0 / ImageWeight();
  for (size_t y = 0; y != height; ++y) {
    size_t ySrc = (height - y) + height / 2;
    if (ySrc >= height) ySrc -= height;

    for (size_t x = 0; x != width; ++x) {
      size_t xSrc = x + width / 2;
      if (xSrc >= width) xSrc -= width;

      _image[x + y * width] = scratch[xSrc + ySrc * width] * weight_factor;
    }
  }
}

template <typename num_t>
inline void DirectMSGridder<num_t>::gridSample(const InversionSample& sample,
                                               size_t layerIndex) {
  // Contribution of one visibility:
  //  I(l, m) = V(u, v, w) exp (2 pi i (ul + vm + w (sqrt(1 - l^2 - m^2) - 1)))
  //   Since every visibility has a conjugate visibility for (-u, -v, -w), we
  //   can simultaneously add:
  // Ic(l, m) = V^*(u, v, w) exp (-2 pi i (ul + vm + w (sqrt(1 - l^2 - m^2) -
  // 1)))
  //   Adding those together gives one real value:
  //     I+Ic = real(V) 2 cos (2 pi (ul + vm + w (sqrt(1 - l^2 - m^2) - 1))) -
  //            imag(V) 2 sin (2 pi (ul + vm + w (sqrt(1 - l^2 - m^2) - 1)))
  aocommon::ImageBase<num_t>& layer = _layers[layerIndex];
  const std::complex<num_t> val = sample.sample;
  const size_t width = TrimWidth();
  const size_t height = TrimHeight();
  constexpr num_t minTwoPi = num_t(-2.0 * M_PI);
  const num_t u = sample.uInLambda;
  const num_t v = sample.vInLambda;
  const num_t w = sample.wInLambda;

  for (size_t y = 0; y != height; ++y) {
    const size_t yIndex = y * height;

    size_t ySrc = (height - y) + height / 2;
    if (ySrc >= height) ySrc -= height;
    const num_t m =
        num_t(((num_t)ySrc - (height / 2)) * PixelSizeY() + MShift());

    for (size_t x = 0; x != width; ++x) {
      size_t xSrc = x + width / 2;
      if (xSrc >= width) xSrc -= width;
      const num_t l =
          num_t(((width / 2) - (num_t)xSrc) * PixelSizeX() + LShift());

      const size_t index = yIndex + x;
      const num_t angle = minTwoPi * (u * l + v * m + w * _sqrtLMTable[index]);
      layer[index] +=
          val.real() * std::cos(angle) - val.imag() * std::sin(angle);
    }
  }
}

template <typename num_t>
void DirectMSGridder<num_t>::invertMeasurementSet(
    const MsProviderCollection::MsData& msData, ProgressBar& progress,
    size_t msIndex) {
  StartMeasurementSet(msData, false);
  const size_t n_vis_polarizations = msData.ms_provider->NPolarizations();
  const aocommon::BandData selectedBand(msData.SelectedBand());
  aocommon::UVector<std::complex<float>> modelBuffer(
      selectedBand.ChannelCount() * n_vis_polarizations);
  aocommon::UVector<float> weightBuffer(selectedBand.ChannelCount() *
                                        n_vis_polarizations);
  aocommon::UVector<bool> isSelected(selectedBand.ChannelCount(), true);

  InversionRow newItem;
  aocommon::UVector<std::complex<float>> newItemData(
      selectedBand.ChannelCount() * n_vis_polarizations);
  newItem.data = newItemData.data();

  std::vector<size_t> idToMSRow;
  msData.ms_provider->MakeIdToMSRowMapping(idToMSRow);
  size_t rowIndex = 0;
  std::unique_ptr<MSReader> msReader = msData.ms_provider->MakeReader();
  while (msReader->CurrentRowAvailable()) {
    progress.SetProgress(msIndex * idToMSRow.size() + rowIndex,
                         GetMsCount() * idToMSRow.size());

    MSProvider::MetaData metaData;
    msReader->ReadMeta(metaData);
    newItem.uvw[0] = metaData.uInM;
    newItem.uvw[1] = metaData.vInM;
    newItem.uvw[2] = metaData.wInM;

    GetCollapsedVisibilities(*msReader, msData.antenna_names.size(), newItem,
                             selectedBand, weightBuffer.data(),
                             modelBuffer.data(), isSelected.data(), metaData);
    InversionSample sample;
    for (size_t ch = 0; ch != selectedBand.ChannelCount(); ++ch) {
      const double wl = selectedBand.ChannelWavelength(ch);
      sample.uInLambda = newItem.uvw[0] / wl;
      sample.vInLambda = newItem.uvw[1] / wl;
      sample.wInLambda = newItem.uvw[2] / wl;
      sample.sample = newItem.data[ch];
      _inversionLane.write(sample);
    }

    msReader->NextInputRow();
    ++rowIndex;
  }
}

template <typename num_t>
void DirectMSGridder<num_t>::Predict(std::vector<aocommon::Image>&& /*image*/) {
  throw std::runtime_error(
      "Prediction not yet implemented for direct FT gridding");
}

template <typename num_t>
aocommon::ImageBase<num_t> DirectMSGridder<num_t>::GetSqrtLMLookupTable()
    const {
  const size_t width = TrimWidth();
  const size_t height = TrimHeight();
  aocommon::ImageBase<num_t> sqrtLMTable(ImageWidth(), ImageHeight());
  num_t* iter = sqrtLMTable.Data();
  for (size_t y = 0; y != height; ++y) {
    size_t ySrc = (height - y) + height / 2;
    if (ySrc >= height) ySrc -= height;
    num_t m = num_t(((num_t)ySrc - (height / 2)) * PixelSizeY() + MShift());

    for (size_t x = 0; x != width; ++x) {
      size_t xSrc = x + width / 2;
      if (xSrc >= width) xSrc -= width;
      num_t l = num_t(((width / 2) - (num_t)xSrc) * PixelSizeX() + LShift());

      if (l * l + m * m < 1.0)
        *iter = std::sqrt(1.0 - l * l - m * m) - 1.0;
      else
        *iter = 0.0;
      ++iter;
    }
  }
  return sqrtLMTable;
}

template class DirectMSGridder<float>;
template class DirectMSGridder<double>;
template class DirectMSGridder<long double>;

}  // namespace wsclean
