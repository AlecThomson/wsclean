#ifndef POLYNOMIAL_CHANNEL_FITTER_H
#define POLYNOMIAL_CHANNEL_FITTER_H

#include <vector>

class PolynomialChannelFitter {
 public:
  void Clear() {
    _channels.clear();
    _dataPoints.clear();
  }

  void AddChannel(double startFrequency, double endFrequency) {
    _channels.emplace_back(startFrequency, endFrequency);
  }

  void AddDataPoint(size_t channel, double y) {
    _dataPoints.emplace_back(channel, y);
  }

  void Fit(std::vector<double>& terms, size_t nTerms);

  static double Evaluate(double x, const std::vector<double>& terms);

 private:
  /**
   * Start and end frequencies of the channels
   */
  std::vector<std::pair<double, double>> _channels;
  std::vector<std::pair<size_t, double>> _dataPoints;
};

#endif
