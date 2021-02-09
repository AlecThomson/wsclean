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
            size_t freqIndex, bool isImaginary, size_t facetIndex = 0) const {
    if (_writer.Width() == 0 || _writer.Height() == 0)
      throw std::runtime_error("Writer is not set.");
    Logger::Debug << "Loading "
                  << name(polarization, freqIndex, facetIndex, isImaginary)
                  << '\n';
    // TODO: maybe _facetCount <= 1 is fine?
    if (_polCount == 1 && _freqCount == 1 && _facetCount == 0)
      if (_image.empty())
        throw std::runtime_error("Loading image before store");
      else
        std::copy(_image.data(),
                  _image.data() + _writer.Width() * _writer.Height(), image);
    else {
      FitsReader reader(name(polarization, freqIndex, facetIndex, isImaginary));
      reader.Read(image);
    }
  }

  template <typename NumT>
  void Store(const NumT* image, aocommon::PolarizationEnum polarization,
             size_t freqIndex, bool isImaginary, size_t facetIndex = 0) {
    if (_writer.Width() == 0 || _writer.Height() == 0)
      throw std::runtime_error("Writer is not set.");
    Logger::Debug << "Storing "
                  << name(polarization, freqIndex, facetIndex, isImaginary)
                  << '\n';
    // TODO: maybe _facetCount <= 1 is fine?
    if (_polCount == 1 && _freqCount == 1 && _facetCount == 0) {
      if (_image.empty()) {
        _image = Image(_writer.Width(), _writer.Height());
      }
      std::copy(image, image + _writer.Width() * _writer.Height(),
                _image.data());
    } else {
      std::string n = name(polarization, freqIndex, facetIndex, isImaginary);
      _writer.Write(n, image);
      _storedNames.insert(n);
    }
  }

 private:
  std::string name(aocommon::PolarizationEnum polarization, size_t freqIndex,
                   size_t facetIndex, bool isImaginary) const {
    if (_freqCount == 1) {
      if (isImaginary)
        return _prefix + '-' +
               aocommon::Polarization::TypeToShortString(polarization) +
               "i-tmp.fits";
      else
        return _prefix + '-' +
               aocommon::Polarization::TypeToShortString(polarization) +
               "-tmp.fits";
    } else {
      std::ostringstream str;
      str << _prefix + '-' +
                 aocommon::Polarization::TypeToShortString(polarization);
      if (isImaginary) str << 'i';
      str << '-';
      if (freqIndex < 10) str << '0';
      if (freqIndex < 100) str << '0';
      if (freqIndex < 1000) str << '0';
      str << freqIndex;
      if (_facetCount > 0) {
        str << "-FACET_ID-";
        if (facetIndex < 10) str << '0';
        if (freqIndex < 100) str << '0';
        if (freqIndex < 1000) str << '0';
        str << facetIndex;
      }
      str << "-tmp.fits";
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
