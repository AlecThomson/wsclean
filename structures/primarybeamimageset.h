#ifndef PRIMARY_BEAM_IMAGE_SET
#define PRIMARY_BEAM_IMAGE_SET

#include <vector>

#include <boost/filesystem/operations.hpp>

#include <aocommon/hmatrix4x4.h>
#include <aocommon/matrix2x2.h>
#include <aocommon/polarization.h>

#include <aocommon/fits/fitsreader.h>

#include "../structures/image.h"

class PrimaryBeamImageSet {
 public:
  PrimaryBeamImageSet() : _beamImages(), _width(0), _height(0) {}

  PrimaryBeamImageSet(size_t width, size_t height, size_t nImages)
      : _beamImages(nImages), _width(width), _height(height) {
    for (size_t i = 0; i != nImages; ++i) _beamImages[i] = Image(width, height);
  }

  static aocommon::FitsReader GetAReader(const std::string& beamPrefix,
                                         bool stokesIOnly) {
    if (stokesIOnly) {
      return aocommon::FitsReader(beamPrefix + "-I.fits");
    } else {
      std::string aFilename = beamPrefix + "-XX.fits";
      if (!boost::filesystem::exists(beamPrefix + "-XX.fits"))
        aFilename = beamPrefix + "-0.fits";
      return aocommon::FitsReader(aFilename);
    }
  }

  static PrimaryBeamImageSet Load(const std::string& beamPrefix, size_t width,
                                  size_t height, bool stokesIOnly) {
    if (stokesIOnly) {
      PrimaryBeamImageSet beamImages(width, height, 8);
      // IDG produces only a Stokes I beam, and has already corrected for the
      // rest. Currently we just load that beam into real component of XX and
      // YY, and set the other 6 images to zero. This is a bit wasteful so might
      // require a better strategy for big images.
      aocommon::FitsReader reader(beamPrefix + "-I.fits");
      reader.Read(beamImages[0].data());
      for (size_t i = 0; i != width * height; ++i)
        beamImages[0][i] = std::sqrt(beamImages[0][i]);
      std::copy_n(beamImages[0].data(), width * height, beamImages[6].data());
      for (size_t i = 1; i != 8; ++i) {
        if (i != 6) std::fill_n(beamImages[i].data(), width * height, 0.0);
      }
      return beamImages;
    } else {
      try {
        PrimaryBeamImageSet beamImages(width, height, 8);
        aocommon::PolarizationEnum linPols[4] = {
            aocommon::Polarization::XX, aocommon::Polarization::XY,
            aocommon::Polarization::YX, aocommon::Polarization::YY};
        for (size_t i = 0; i != 8; ++i) {
          aocommon::PolarizationEnum p = linPols[i / 2];
          std::string polStr;
          if (i % 2 == 0)  // real?
            polStr = aocommon::Polarization::TypeToShortString(p);
          else
            polStr = aocommon::Polarization::TypeToShortString(p) + "i";
          aocommon::FitsReader reader(beamPrefix + "-" + polStr + ".fits");
          reader.Read(beamImages[i].data());
        }
        return beamImages;
      } catch (std::exception&) {
        PrimaryBeamImageSet beamImages(width, height, 16);
        for (size_t i = 0; i != 16; ++i) {
          aocommon::FitsReader reader(beamPrefix + "-" + std::to_string(i) +
                                      ".fits");
          reader.Read(beamImages[i].data());
        }
        return beamImages;
      }
    }
  }

  void SetToZero() {
    for (Image& img : _beamImages)
      std::fill_n(img.data(), _width * _height, 0.0);
  }

  double GetUnpolarizedCorrectionFactor(size_t x, size_t y) {
    if (_beamImages.size() == 16) {
      const size_t j = y * _width + x;
      aocommon::HMC4x4 beam = aocommon::HMC4x4::FromData(
          {_beamImages[0][j], _beamImages[1][j], _beamImages[2][j],
           _beamImages[3][j], _beamImages[4][j], _beamImages[5][j],
           _beamImages[6][j], _beamImages[7][j], _beamImages[8][j],
           _beamImages[9][j], _beamImages[10][j], _beamImages[11][j],
           _beamImages[12][j], _beamImages[13][j], _beamImages[14][j],
           _beamImages[15][j]});
      if (!beam.Invert()) beam = aocommon::HMC4x4::Zero();
      aocommon::Vector4 v{0.5, 0.0, 0.0, 0.5};
      v = beam * v;
      return v[0].real() + v[3].real();
    } else {
      throw std::runtime_error(
          "PrimaryBeamImageSet::GetUnpolarizedCorrectionFactor(): Expected 16 "
          "images, but found " +
          std::to_string(_beamImages.size()));
    }
  }

  const Image& operator[](size_t index) const { return _beamImages[index]; }
  Image& operator[](size_t index) { return _beamImages[index]; }

  PrimaryBeamImageSet& operator+=(const PrimaryBeamImageSet& rhs) {
    if (_beamImages.size() != rhs._beamImages.size())
      throw std::runtime_error("Primary beam image sets don't match");
    for (size_t i = 0; i != _beamImages.size(); ++i) {
      for (size_t j = 0; j != _width * _height; ++j)
        _beamImages[i][j] += rhs._beamImages[i][j];
    }
    return *this;
  }

  PrimaryBeamImageSet& operator*=(double factor) {
    for (size_t i = 0; i != _beamImages.size(); ++i) {
      for (size_t j = 0; j != _width * _height; ++j)
        _beamImages[i][j] *= factor;
    }
    return *this;
  }

  void ApplyStokesI(float* stokesI, double beamLimit) const {
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beamLimit = beamLimit * beamLimit;
    const size_t size = _width * _height;
    for (size_t j = 0; j != size; ++j) {
      aocommon::HMC4x4 beam = aocommon::HMC4x4::FromData(
          {_beamImages[0][j], _beamImages[1][j], _beamImages[2][j],
           _beamImages[3][j], _beamImages[4][j], _beamImages[5][j],
           _beamImages[6][j], _beamImages[7][j], _beamImages[8][j],
           _beamImages[9][j], _beamImages[10][j], _beamImages[11][j],
           _beamImages[12][j], _beamImages[13][j], _beamImages[14][j],
           _beamImages[15][j]});
      if (beam.Norm() > beamLimit) {
        if (!beam.Invert()) beam = aocommon::HMC4x4::Zero();
        const float factor = stokesI[j] * 0.5;
        aocommon::Vector4 v{factor, 0.0, 0.0, factor};
        v = beam * v;
        stokesI[j] = v[0].real() + v[3].real();
      } else {
        stokesI[j] = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  void ApplyDiagonal(float* images[2], double beamLimit) {
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beamLimit = beamLimit * beamLimit;
    size_t size = _width * _height;
    for (size_t j = 0; j != size; ++j) {
      aocommon::HMC4x4 beam = aocommon::HMC4x4::FromData(
          {_beamImages[0][j], _beamImages[1][j], _beamImages[2][j],
           _beamImages[3][j], _beamImages[4][j], _beamImages[5][j],
           _beamImages[6][j], _beamImages[7][j], _beamImages[8][j],
           _beamImages[9][j], _beamImages[10][j], _beamImages[11][j],
           _beamImages[12][j], _beamImages[13][j], _beamImages[14][j],
           _beamImages[15][j]});
      if (beam.Norm() > beamLimit) {
        if (!beam.Invert()) beam = aocommon::HMC4x4::Zero();
        aocommon::Vector4 v{images[0][j], 0., 0., images[1][j]};
        // Implementation could be made more efficient - at the expense of
        // explicitness - since only a few entries in the beam matrix are
        // relevant
        v = beam * v;
        images[0][j] = v[0].real();
        images[1][j] = v[3].real();
      } else {
        for (size_t p = 0; p != 2; ++p)
          images[p][j] = std::numeric_limits<float>::quiet_NaN();
      }
    }
  };

  void ApplyFullStokes(float* images[4], double beamLimit) const {
    // The beam will be compared to a squared quantity (matrix norm), so square
    // it:
    beamLimit = beamLimit * beamLimit;
    size_t size = _width * _height;
    for (size_t j = 0; j != size; ++j) {
      aocommon::HMC4x4 beam = aocommon::HMC4x4::FromData(
          {_beamImages[0][j], _beamImages[1][j], _beamImages[2][j],
           _beamImages[3][j], _beamImages[4][j], _beamImages[5][j],
           _beamImages[6][j], _beamImages[7][j], _beamImages[8][j],
           _beamImages[9][j], _beamImages[10][j], _beamImages[11][j],
           _beamImages[12][j], _beamImages[13][j], _beamImages[14][j],
           _beamImages[15][j]});
      if (beam.Norm() > beamLimit) {
        if (!beam.Invert()) beam = aocommon::HMC4x4::Zero();
        double stokesVal[4] = {images[0][j], images[1][j], images[2][j],
                               images[3][j]};
        aocommon::Vector4 v;
        aocommon::Polarization::StokesToLinear(stokesVal, v.data());
        v = beam * v;
        aocommon::Polarization::LinearToStokes(v.data(), stokesVal);
        for (size_t p = 0; p != 4; ++p) images[p][j] = stokesVal[p];
      } else {
        for (size_t p = 0; p != 4; ++p)
          images[p][j] = std::numeric_limits<float>::quiet_NaN();
      }
    }
  }

  size_t NImages() const { return _beamImages.size(); }
  size_t Width() const { return _width; }
  size_t Height() const { return _height; }

 private:
  std::vector<Image> _beamImages;
  size_t _width, _height;
};

#endif
