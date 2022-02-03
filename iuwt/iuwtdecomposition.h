#ifndef IUWT_DECOMPOSITION_H
#define IUWT_DECOMPOSITION_H

#include <aocommon/fits/fitswriter.h>

#include "iuwtmask.h"

#include <aocommon/image.h>
#include <aocommon/staticfor.h>
#include <aocommon/uvector.h>

#include <iostream>
#include <string>
#include <sstream>

class IUWTDecompositionScale {
 public:
  aocommon::Image& Coefficients() { return _coefficients; }
  const aocommon::Image& Coefficients() const { return _coefficients; }
  float& operator[](size_t index) { return _coefficients[index]; }
  const float& operator[](size_t index) const { return _coefficients[index]; }

 private:
  aocommon::Image _coefficients;
};

class IUWTDecomposition {
 public:
  IUWTDecomposition(int scaleCount, size_t width, size_t height)
      : _scales(scaleCount + 1),
        _scaleCount(scaleCount),
        _width(width),
        _height(height) {}

  IUWTDecomposition* CreateTrimmed(int newScaleCount, size_t x1, size_t y1,
                                   size_t x2, size_t y2) const {
    std::unique_ptr<IUWTDecomposition> p(
        new IUWTDecomposition(newScaleCount, x2 - x1, y2 - y1));
    for (int i = 0; i != newScaleCount; ++i) {
      copySmallerPart(_scales[i].Coefficients(), p->_scales[i].Coefficients(),
                      x1, y1, x2, y2);
    }
    p->_scales.back().Coefficients() = aocommon::Image(_width, _height, 0.0);
    return p.release();
  }

  void Convolve(aocommon::Image& image, int toScale) {
    aocommon::Image scratch(image.Width(), image.Height());
    for (int scale = 0; scale != toScale; ++scale) {
      convolve(image.Data(), image.Data(), scratch.Data(), _width, _height,
               scale + 1);
    }
  }

  void DecomposeSimple(aocommon::Image& input) {
    aocommon::Image scratch(input.Width(), input.Height());
    for (int scale = 0; scale != int(_scaleCount); ++scale) {
      aocommon::Image& coefficients = _scales[scale].Coefficients();
      coefficients = aocommon::Image(_width, _height);

      aocommon::Image tmp(_width, _height);
      convolve(tmp.Data(), input.Data(), scratch.Data(), _width, _height,
               scale);
      difference(coefficients.Data(), input.Data(), tmp.Data(), _width,
                 _height);
      memcpy(input.Data(), tmp.Data(), sizeof(float) * _width * _height);
    }
    _scales.back().Coefficients() = input;
  }

  void RecomposeSimple(aocommon::Image& output) {
    output = _scales[0].Coefficients();
    for (size_t scale = 1; scale != _scaleCount - 1; ++scale) {
      for (size_t i = 0; i != _width * _height; ++i)
        output[i] += _scales[scale][i];
    }
  }

  void Decompose(aocommon::StaticFor<size_t>& loop, const float* input,
                 float* scratch, bool includeLargest) {
    DecomposeMT(loop, input, scratch, includeLargest);
  }

  void DecomposeMT(aocommon::StaticFor<size_t>& loop, const float* input,
                   float* scratch, bool includeLargest);

  void DecomposeST(const float* input, float* scratch) {
    aocommon::UVector<float> i0(input, input + _width * _height);
    aocommon::Image i1(_width, _height);
    aocommon::Image i2(_width, _height);
    for (int scale = 0; scale != int(_scaleCount); ++scale) {
      aocommon::Image& coefficients = _scales[scale].Coefficients();
      coefficients = aocommon::Image(_width, _height);
      convolve(i1.Data(), i0.data(), scratch, _width, _height, scale + 1);
      convolve(i2.Data(), i1.Data(), scratch, _width, _height, scale + 1);

      // coefficients = i0 - i2
      difference(coefficients.Data(), i0.data(), i2.Data(), _width, _height);

      // i0 = i1;
      if (scale + 1 != int(_scaleCount))
        memcpy(i0.data(), i1.Data(), sizeof(float) * _width * _height);
    }
    _scales.back().Coefficients() = i1;
  }

  void Recompose(aocommon::Image& output, bool includeLargest) {
    aocommon::Image scratch1(_width, _height);
    aocommon::Image scratch2(_width, _height);
    bool isZero;
    if (includeLargest) {
      output = _scales.back().Coefficients();
      isZero = false;
    } else {
      output = aocommon::Image(_width, _height, 0.0);
      isZero = true;
    }
    for (int scale = int(_scaleCount) - 1; scale != -1; --scale) {
      const aocommon::Image& coefficients = _scales[scale].Coefficients();
      if (isZero) {
        output = coefficients;
        isZero = false;
      } else {
        // output = output (x) IUWT
        convolve(scratch2.Data(), output.Data(), scratch1.Data(), _width,
                 _height, scale + 1);
        output = scratch2;

        // output += coefficients
        for (size_t i = 0; i != output.Size(); ++i)
          output[i] += coefficients[i];
      }
    }
  }

  size_t NScales() const { return _scaleCount; }

  size_t Width() const { return _width; }

  size_t Height() const { return _height; }

  IUWTDecompositionScale& operator[](int scale) { return _scales[scale]; }

  const IUWTDecompositionScale& operator[](int scale) const {
    return _scales[scale];
  }

  void ApplyMask(const IUWTMask& mask) {
    for (size_t scale = 0; scale != _scaleCount; ++scale) {
      for (size_t i = 0; i != _scales[scale].Coefficients().Size(); ++i) {
        if (!mask[scale][i]) _scales[scale][i] = 0.0;
      }
    }
    _scales[_scaleCount].Coefficients() = aocommon::Image(_width, _height, 0.0);
  }

  void Save(const std::string& prefix) {
    std::cout << "Saving scales...\n";
    aocommon::FitsWriter writer;
    writer.SetImageDimensions(_width, _height);
    for (size_t scale = 0; scale != _scales.size(); ++scale) {
      std::ostringstream str;
      str << prefix << "-iuwt-" << scale << ".fits";
      writer.Write(str.str(), _scales[scale].Coefficients().Data());
    }
  }

  static int EndScale(size_t maxImageDimension) {
    return std::max(int(log2(maxImageDimension)) - 3, 2);
  }

  static size_t MinImageDimension(int endScale) {
    return (1 << (endScale + 3));
  }

  std::string Summary() const {
    std::ostringstream str;
    str << "IUWTDecomposition, NScales()=" << NScales()
        << ", MinImageDimension()=" << MinImageDimension(NScales() + 1)
        << ", width=" << _width << ", height=" << _height;
    return str.str();
  }

 private:
  static void convolveComponentHorizontal(const float* input, float* output,
                                          size_t width, size_t height,
                                          float val, int dist) {
    size_t minX = std::max<int>(0, -dist),
           maxX = std::min<int>(width, width - dist);
    for (size_t y = 0; y != height; ++y) {
      float* outputPtr = &output[y * width];
      const float* inputPtr = &input[y * width];
      for (size_t x = minX; x < maxX; ++x) {
        outputPtr[x] += inputPtr[x + dist] * val;
      }
    }
  }

  static void convolveComponentVertical(const float* input, float* output,
                                        size_t width, size_t height, float val,
                                        int dist) {
    size_t minY = std::max<int>(0, -dist),
           maxY = std::min<int>(height, height - dist);
    for (size_t y = minY; y < maxY; ++y) {
      float* outputPtr = &output[y * width];
      const float* inputPtr = &input[(y + dist) * width];
      for (size_t x = 0; x != width; ++x) {
        outputPtr[x] += inputPtr[x] * val;
      }
    }
  }

  static void convolveComponentVerticalPartial(const float* input,
                                               float* output, size_t width,
                                               size_t height, size_t startX,
                                               size_t endX, float val,
                                               int dist) {
    size_t minY = std::max<int>(0, -dist),
           maxY = std::min<int>(height, height - dist);
    for (size_t y = minY; y < maxY; ++y) {
      float* outputPtr = &output[y * width];
      const float* inputPtr = &input[(y + dist) * width];
      for (size_t x = startX; x != endX; ++x) {
        outputPtr[x] += inputPtr[x] * val;
      }
    }
  }

  static void convolve(float* output, const float* image, float* scratch,
                       size_t width, size_t height, int scale) {
    for (size_t i = 0; i != width * height; ++i) scratch[i] = 0.0;
    const size_t H_SIZE = 5;
    const float h[H_SIZE] = {1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0,
                             1.0 / 16.0};
    int scaleDist = (1 << scale);
    for (int hIndex = 0; hIndex != H_SIZE; ++hIndex) {
      int hShift = hIndex - H_SIZE / 2;
      convolveComponentHorizontal(image, scratch, width, height, h[hIndex],
                                  (scaleDist - 1) * hShift);
    }
    for (size_t i = 0; i != width * height; ++i) output[i] = 0.0;
    for (int hIndex = 0; hIndex != H_SIZE; ++hIndex) {
      int hShift = hIndex - H_SIZE / 2;
      convolveComponentVertical(scratch, output, width, height, h[hIndex],
                                (scaleDist - 1) * hShift);
    }
  }

  static void convolveMT(aocommon::StaticFor<size_t>& loop, float* output,
                         const float* image, float* scratch, size_t width,
                         size_t height, int scale);

  static void convolveHorizontalPartial(float* output, const float* image,
                                        size_t width, size_t startY,
                                        size_t endY, int scale) {
    size_t startIndex = startY * width;
    convolveHorizontalFast(&output[startIndex], &image[startIndex], width,
                           endY - startY, scale);
  }

  static void convolveHorizontal(float* output, const float* image,
                                 size_t width, size_t height, int scale) {
    for (size_t i = 0; i != width * height; ++i) output[i] = 0.0;
    const size_t H_SIZE = 5;
    const float h[H_SIZE] = {1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0,
                             1.0 / 16.0};
    int scaleDist = (1 << scale);
    for (int hIndex = 0; hIndex != H_SIZE; ++hIndex) {
      int hShift = hIndex - H_SIZE / 2;
      convolveComponentHorizontal(image, output, width, height, h[hIndex],
                                  (scaleDist - 1) * hShift);
    }
  }

  static void convolveHorizontalFast(float* output, const float* image,
                                     size_t width, size_t height, int scale);

  static void convolveVerticalPartial(float* output, const float* image,
                                      size_t width, size_t height,
                                      size_t startX, size_t endX, int scale) {
    for (size_t i = 0; i != width * height; ++i) output[i] = 0.0;
    const size_t H_SIZE = 5;
    const float h[H_SIZE] = {1.0 / 16.0, 4.0 / 16.0, 6.0 / 16.0, 4.0 / 16.0,
                             1.0 / 16.0};
    int scaleDist = (1 << scale);
    for (int hIndex = 0; hIndex != H_SIZE; ++hIndex) {
      int hShift = hIndex - H_SIZE / 2;
      convolveComponentVerticalPartial(image, output, width, height, startX,
                                       endX, h[hIndex],
                                       (scaleDist - 1) * hShift);
    }
  }

  static void convolveVerticalPartialFast(float* output, const float* image,
                                          size_t width, size_t height,
                                          size_t startX, size_t endX,
                                          int scale);

  static void convolveVerticalPartialFastFailed(float* output,
                                                const float* image,
                                                size_t width, size_t height,
                                                size_t startX, size_t endX,
                                                int scale);

  static void differencePartial(float* dest, const float* lhs, const float* rhs,
                                size_t width, size_t startY, size_t endY) {
    size_t startIndex = startY * width;
    difference(&dest[startIndex], &lhs[startIndex], &rhs[startIndex], width,
               endY - startY);
  }

  static void differenceMT(aocommon::StaticFor<size_t>& loop, float* dest,
                           const float* lhs, const float* rhs, size_t width,
                           size_t height);

  static void difference(float* dest, const float* lhs, const float* rhs,
                         size_t width, size_t height) {
    for (size_t i = 0; i != width * height; ++i) {
      dest[i] = lhs[i] - rhs[i];
    }
  }

  void copySmallerPart(const aocommon::Image& input, aocommon::Image& output,
                       size_t x1, size_t y1, size_t x2, size_t y2) const {
    size_t newWidth = x2 - x1;
    output = aocommon::Image(newWidth, y2 - y1);
    for (size_t y = y1; y != y2; ++y) {
      const float* oldPtr = &input[y * _width];
      float* newPtr = &output[(y - y1) * newWidth];
      for (size_t x = x1; x != x2; ++x) {
        newPtr[x - x1] = oldPtr[x];
      }
    }
  }

 private:
  std::vector<IUWTDecompositionScale> _scales;
  size_t _scaleCount, _width, _height;
};

#endif
