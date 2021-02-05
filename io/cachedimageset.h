#ifndef CACHED_IMAGE_SET_H
#define CACHED_IMAGE_SET_H

#include <aocommon/fits/fitsreader.h>
#include <aocommon/fits/fitswriter.h>

#include "logger.h"

#include <string.h>
#include <set>

using aocommon::FitsReader;
using aocommon::FitsWriter;

class CachedImageSet {
 public:
  CachedImageSet() : _polCount(0), _freqCount(0), _image() {}

  ~CachedImageSet() {
    for (const std::string& filename : _storedNames)
      std::remove(filename.c_str());
  }

  CachedImageSet(const CachedImageSet& source) = delete;
  CachedImageSet& operator=(const CachedImageSet& source) = delete;

  void Initialize(const FitsWriter& writer, size_t polCount, size_t freqCount,
                  const std::string& prefix, size_t facetCount = 0) {
    _writer = writer;
    _polCount = polCount;
    _freqCount = freqCount;
    _facetCount = facetCount;
    _prefix = prefix;
    _image.reset();
  }

  void SetFitsWriter(const FitsWriter& writer) { _writer = writer; }

  template <typename NumT>
  void Load(NumT* image, aocommon::PolarizationEnum polarization,
            size_t freqIndex, bool isImaginary) const {
    if (_writer.Width() == 0 || _writer.Height() == 0)
      throw std::runtime_error("Writer is not set.");
    Logger::Debug << "Loading " << name(polarization, freqIndex, isImaginary)
                  << '\n';
    // TODO: maybe _facetCount <= 1 is fine?
    if (_polCount == 1 && _freqCount == 1 && _facetCount == 0)
      if (_image.empty())
        throw std::runtime_error("Loading image before store");
      else
        std::copy(_image.data(),
                  _image.data() + _writer.Width() * _writer.Height(), image);
    else {
      FitsReader reader(name(polarization, freqIndex, isImaginary));
      reader.Read(image);
    }
  }

  template <typename NumT>
  void LoadFacet(NumT* image, aocommon::PolarizationEnum polarization,
                 size_t freqIndex, bool isImaginary, size_t facetIndex) const {
    if (_facetCount == 0) {
      Load<NumT>(image, polarization, freqIndex, isImaginary);
    } else {
      // No need to check width and height here, because we do not cache
      std::string n =
          name_facet(polarization, freqIndex, facetIndex, isImaginary);
      Logger::Debug << "Loading " << n << '\n';
      FitsReader reader(n);
      reader.Read(image);
    }
  }

  template <typename NumT>
  void Store(const NumT* image, aocommon::PolarizationEnum polarization,
             size_t freqIndex, bool isImaginary) {
    if (_writer.Width() == 0 || _writer.Height() == 0)
      throw std::runtime_error("Writer is not set.");
    Logger::Debug << "Storing " << name(polarization, freqIndex, isImaginary)
                  << '\n';
    // TODO: maybe leave out _facetCount, or use _facetCount <= 1?
    if (_polCount == 1 && _freqCount == 1 && _facetCount == 0) {
      if (_image.empty()) {
        _image = Image(_writer.Width(), _writer.Height());
      }
      std::copy(image, image + _writer.Width() * _writer.Height(),
                _image.data());
    } else {
      std::string n = name(polarization, freqIndex, isImaginary);
      _writer.Write(n, image);
      _storedNames.insert(n);
    }
  }

  template <typename NumT>
  void StoreFacet(const NumT* image, aocommon::PolarizationEnum polarization,
                  size_t freqIndex, bool isImaginary, size_t facetIndex,
                  size_t facet_width, size_t facet_height) {
    if (_facetCount == 0) {
      // If _facetCount 0, use the main "Store" as is
      Store<NumT>(image, polarization, freqIndex, isImaginary);
    } else {
      std::string n =
          name_facet(polarization, freqIndex, facetIndex, isImaginary);
      Logger::Debug << "Storing " << n << '\n';

      // Initialize FacetWriter, using facet_width and facet_height
      FitsWriter facetWriter;
      facetWriter.SetImageDimensions(facet_width, facet_height);
      facetWriter.Write(n, image);
      _storedNames.insert(n);
    }
  }

 private:
  std::string name(aocommon::PolarizationEnum polarization, size_t freqIndex,
                   bool isImaginary) const {
    std::ostringstream str;
    str << name_trunk(polarization, freqIndex, isImaginary);
    str << "-tmp.fits";
    return str.str();
  }

  std::string name_facet(aocommon::PolarizationEnum polarization,
                         size_t freqIndex, size_t facetIndex,
                         bool isImaginary) const {
    std::ostringstream str;
    str << name_trunk(polarization, freqIndex, isImaginary);
    if (_facetCount > 0) {
      str << "-facet-";
      if (facetIndex < 10) str << '0';
      if (freqIndex < 100) str << '0';
      if (freqIndex < 1000) str << '0';
      str << facetIndex;
    }
    str << "-tmp.fits";
    return str.str();
  }

  std::string name_trunk(aocommon::PolarizationEnum polarization,
                         size_t freqIndex, bool isImaginary) const {
    if (_freqCount == 1) {
      if (isImaginary)
        return _prefix + '-' +
               aocommon::Polarization::TypeToShortString(polarization) + "i";
      else
        return _prefix + '-' +
               aocommon::Polarization::TypeToShortString(polarization);
    } else {
      std::ostringstream str;
      str << _prefix + '-' +
                 aocommon::Polarization::TypeToShortString(polarization);
      if (isImaginary) str << 'i';
      str << '-';
      if (freqIndex < 10) str << '0';
      if (freqIndex < 100) str << '0';
      if (freqIndex < 1000) str << '0';
      return str.str();
    }
  }

  FitsWriter _writer;
  size_t _polCount, _freqCount, _facetCount;
  std::string _prefix;

  Image _image;
  std::set<std::string> _storedNames;
};

#endif
