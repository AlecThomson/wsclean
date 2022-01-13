#ifndef SYSTEM_DP3_H
#define SYSTEM_DP3_H

#include <fstream>
#include <boost/filesystem/operations.hpp>
#include <aocommon/radeccoord.h>
#include <aocommon/uvector.h>

class DP3 {
 public:
  static void WriteStandardHeader(std::ostream& stream, double refFrequency) {
    // TODO: not used inside WSClean
    stream.precision(15);
    stream
        << "# (Name, Type, Ra, Dec, I, Q, U, V, ReferenceFrequency='"
        << refFrequency
        << "', SpectralIndex, MajorAxis, MinorAxis, Orientation) = format\n\n";
  }

  static void WriteHeaderForSpectralTerms(std::ostream& stream,
                                          double refFrequency) {
    stream.precision(15);
    stream << "Format = Name, Type, Ra, Dec, I, SpectralIndex, LogarithmicSI, "
              "ReferenceFrequency='"
           << refFrequency << "', MajorAxis, MinorAxis, Orientation\n";
  }

  static void WriteOldHeaderForSpectralTerms(
      std::ostream& stream, double refFrequency,
      const std::string& spectralFunction) {
    stream.precision(15);
    stream << "# (Name, Type, Ra, Dec, SpectralTerms, MajorAxis, MinorAxis, "
              "Orientation) = format\n"
           << "# ReferenceFrequency = " << refFrequency << '\n'
           << "# SpectralFunction = " << spectralFunction << '\n';
  }

  static void addSITerms(std::ostream& stream,
                         const aocommon::UVector<float>& siTerms) {
    stream << '[';
    if (!siTerms.empty()) {
      stream << siTerms[0];
      for (size_t i = 1; i != siTerms.size(); ++i) {
        stream << ',' << siTerms[i];
      }
    }
    stream << ']';
  }

  static void WritePointComponent(std::ostream& stream, const std::string& name,
                                  long double ra, long double dec, double i,
                                  double q, double u, double v, double freq,
                                  const aocommon::UVector<float>& siTerms) {
    stream << name << ",POINT," << aocommon::RaDecCoord::RAToString(ra, ':')
           << ',' << aocommon::RaDecCoord::DecToString(dec, '.') << ',' << i
           << ',' << q << ',' << u << ',' << v << ',' << freq << ",";
    addSITerms(stream, siTerms);
    stream << ",,,\n";
  }

  static void WriteGaussianComponent(std::ostream& stream,
                                     const std::string& name, long double ra,
                                     long double dec, double i, double q,
                                     double u, double v, double freq,
                                     const aocommon::UVector<float>& siTerms,
                                     double maj, double min, double posangle) {
    stream << name << ",GAUSSIAN," << aocommon::RaDecCoord::RAToString(ra, ':')
           << ',' << aocommon::RaDecCoord::DecToString(dec, '.') << ',' << i
           << ',' << q << ',' << u << ',' << v << ',' << freq << ",";
    addSITerms(stream, siTerms);
    stream << "," << maj << ',' << min << ',' << posangle << "\n";
  }

  static void WritePolynomialPointComponent(
      std::ostream& stream, const std::string& name, long double ra,
      long double dec, double i, bool useLogSI,
      const aocommon::UVector<float>& polTerms, double referenceFrequencyHz) {
    stream << name << ",POINT," << aocommon::RaDecCoord::RAToString(ra, ':')
           << ',' << aocommon::RaDecCoord::DecToString(dec, '.') << ',' << i
           << ',';
    addSITerms(stream, polTerms);
    stream << "," << (useLogSI ? "true" : "false") << ","
           << referenceFrequencyHz << ",,,\n";
  }

  static void WritePolynomialGaussianComponent(
      std::ostream& stream, const std::string& name, long double ra,
      long double dec, double i, bool useLogSI,
      const aocommon::UVector<float>& polTerms, double referenceFrequencyHz,
      double maj, double min, double posangle) {
    stream << name << ",GAUSSIAN," << aocommon::RaDecCoord::RAToString(ra, ':')
           << ',' << aocommon::RaDecCoord::DecToString(dec, '.') << ',' << i
           << ',';
    addSITerms(stream, polTerms);
    stream << "," << (useLogSI ? "true" : "false") << ","
           << referenceFrequencyHz << "," << maj << ',' << min << ','
           << posangle << "\n";
  }

  static void WriteOldPolynomialPointComponent(
      std::ostream& stream, const std::string& name, long double ra,
      long double dec, const aocommon::UVector<float>& polTerms) {
    stream << name << ",POINT," << aocommon::RaDecCoord::RAToString(ra, ':')
           << ',' << aocommon::RaDecCoord::DecToString(dec, '.') << ',';
    addSITerms(stream, polTerms);
    stream << ",,,\n";
  }

  static void WriteOldPolynomialGaussianComponent(
      std::ostream& stream, const std::string& name, long double ra,
      long double dec, const aocommon::UVector<float>& polTerms, double maj,
      double min, double posangle) {
    stream << name << ",GAUSSIAN," << aocommon::RaDecCoord::RAToString(ra, ':')
           << ',' << aocommon::RaDecCoord::DecToString(dec, '.') << ',';
    addSITerms(stream, polTerms);
    stream << "," << maj << ',' << min << ',' << posangle << "\n";
  }
};

#endif
