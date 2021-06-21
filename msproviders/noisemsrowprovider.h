#ifndef NOISE_MS_ROW_PROVIDER_H
#define NOISE_MS_ROW_PROVIDER_H

#include "directmsrowprovider.h"

#include <fstream>
#include <map>
#include <random>

class NoiseMSRowProvider : public DirectMSRowProvider {
 public:
  NoiseMSRowProvider(const string& msPath, const MSSelection& selection,
                     const std::map<size_t, size_t>& selectedDataDescIds,
                     const std::string& dataColumnName, bool requireModel)
      : DirectMSRowProvider(msPath, selection, selectedDataDescIds,
                            dataColumnName, requireModel),
        _rng(std::random_device{}()),
        _distribution(0.0, 1.0) {}

  void SetNoiseLevel(double noiseStdDevJy) {
    _distribution = std::normal_distribution<float>(0.0, noiseStdDevJy);
  }

  void SetNoiseBaselineFile(const std::string& filename) {
    std::ifstream file(filename);
    while (file) {
      size_t ant1, ant2;
      float stddev;
      file >> ant1 >> ant2 >> stddev;
      if (file) {
        if (ant1 > ant2) std::swap(ant1, ant2);
        _noiseMap.emplace(std::make_pair(ant1, ant2), stddev);
      }
    }
  }

  virtual void ReadData(DataArray& data, FlagArray& flags, WeightArray& weights,
                        double& u, double& v, double& w, uint32_t& dataDescId,
                        uint32_t& antenna1, uint32_t& antenna2,
                        uint32_t& fieldId, double& time) final override {
    DirectMSRowProvider::ReadData(data, flags, weights, u, v, w, dataDescId,
                                  antenna1, antenna2, fieldId, time);
    float stddev;
    if (_noiseMap.empty())
      stddev = 1.0;
    else {
      size_t a1 = antenna1, a2 = antenna2;
      if (a1 > a2) std::swap(a1, a2);
      auto iter = _noiseMap.find(std::make_pair(a1, a2));
      if (iter == _noiseMap.end())
        throw std::runtime_error(
            "The following baseline was not present in the baseline noise "
            "map: " +
            std::to_string(a1) + " x " + std::to_string(a2));
      stddev = iter->second;
    }
    for (DataArray::contiter iter = data.cbegin(); iter != data.cend();
         ++iter) {
      if (std::isfinite(iter->real()) && std::isfinite(iter->imag())) {
        iter->real(_distribution(_rng) * stddev);
        iter->imag(_distribution(_rng) * stddev);
      } else {
        iter->real(std::numeric_limits<float>::quiet_NaN());
        iter->imag(std::numeric_limits<float>::quiet_NaN());
      }
    }
  }

 private:
  std::mt19937 _rng;
  std::normal_distribution<float> _distribution;
  // Maps an antenna1, antenna2 pair to the stddev level
  std::map<std::pair<size_t, size_t>, float> _noiseMap;
};

#endif
