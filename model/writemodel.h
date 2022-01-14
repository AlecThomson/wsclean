#ifndef WSCLEAN_MODEL_WRITEMODEL_H_
#define WSCLEAN_MODEL_WRITEMODEL_H_

#include <fstream>
#include <boost/filesystem/operations.hpp>
#include <aocommon/radeccoord.h>
#include <aocommon/uvector.h>

namespace wsclean {
namespace model {
inline void WriteHeaderForSpectralTerms(std::ostream& stream,
                                        double reference_frequency) {
  stream.precision(15);
  stream << "Format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, "
            "ReferenceFrequency='"
         << reference_frequency << "', MajorAxis, MinorAxis, Orientation\n";
}

inline void AddSiTerms(std::ostream& stream,
                       const aocommon::UVector<float>& si_terms) {
  stream << '[';
  if (!si_terms.empty()) {
    stream << si_terms[0];
    for (size_t i = 1; i != si_terms.size(); ++i) {
      stream << ',' << si_terms[i];
    }
  }
  stream << ']';
}

inline void WritePolynomialPointComponent(
    std::ostream& stream, const std::string& name, long double ra,
    long double dec, double i, bool use_log_si,
    const aocommon::UVector<float>& polarization_terms,
    double reference_frequency_hz) {
  stream << name << ",POINT," << aocommon::RaDecCoord::RAToString(ra, ':')
         << ',' << aocommon::RaDecCoord::DecToString(dec, '.') << ',' << i
         << ',';
  AddSiTerms(stream, polarization_terms);
  stream << "," << (use_log_si ? "true" : "false") << ","
         << reference_frequency_hz << ",,,\n";
}

inline void WritePolynomialGaussianComponent(
    std::ostream& stream, const std::string& name, long double ra,
    long double dec, double i, bool use_log_si,
    const aocommon::UVector<float>& polarization_terms,
    double reference_frequency_hz, double maj, double min,
    double position_angle) {
  stream << name << ",GAUSSIAN," << aocommon::RaDecCoord::RAToString(ra, ':')
         << ',' << aocommon::RaDecCoord::DecToString(dec, '.') << ',' << i
         << ',';
  AddSiTerms(stream, polarization_terms);
  stream << "," << (use_log_si ? "true" : "false") << ","
         << reference_frequency_hz << "," << maj << ',' << min << ','
         << position_angle << "\n";
}
}  // namespace model
}  // namespace wsclean
#endif
