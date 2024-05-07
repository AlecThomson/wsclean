#include <boost/test/unit_test.hpp>

#include "../../math/imageoperations.h"

BOOST_AUTO_TEST_SUITE(image_operations)

using aocommon::HMC4x4;
using aocommon::Image;
using aocommon::MC2x2;

void CheckConstantImage(const Image& image, float expected,
                        const std::string& description) {
  for (float value : image) {
    if (std::fabs(value - expected) > 1e-6) {
      std::ostringstream str;
      str << description << ", expecting " << expected << ", got " << value;
      BOOST_CHECK_EQUAL(str.str(), "");
      break;
    }
    BOOST_CHECK_LT(std::fabs(value - expected), 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(correct_images_for_mueller_matrix) {
  Image i(3, 2, 77.0);
  Image q(3, 2, 0.0);
  Image u(3, 2, 0.0);
  Image v(3, 2, 0.0);
  std::array<aocommon::Image*, 4> images{&i, &q, &u, &v};
  math::CorrectImagesForMuellerMatrix(HMC4x4::Unit(), images);
  CheckConstantImage(i, 77.0, "Stokes I");
  CheckConstantImage(q, 0.0, "Stokes Q");
  CheckConstantImage(u, 0.0, "Stokes U");
  CheckConstantImage(v, 0.0, "Stokes V");

  const HMC4x4 i_to_q =
      HMC4x4::KroneckerProduct(MC2x2{2, 0, 0, 0}, MC2x2{1, 0, 0, 0});
  math::CorrectImagesForMuellerMatrix(i_to_q, images);
  CheckConstantImage(i, 77.0, "Stokes I");
  CheckConstantImage(q, 77.0, "Stokes Q");
  CheckConstantImage(u, 0.0, "Stokes U");
  CheckConstantImage(v, 0.0, "Stokes V");

  i = 10.0;
  q = 0.0;
  u = 0.0;
  v = 0.0;
  constexpr std::complex j(0.0, 1.0);
  const HMC4x4 i_to_v{
      0.0, -j,  j,   0.0,  // XX=0
      j,   0.0, 0.0, 0.0,  // XY=iXX
      -j,  0.0, 0.0, 0.0,  // YX=-iXX
      0.0, 0.0, 0.0, 0.0   // YY=0
  };
  math::CorrectImagesForMuellerMatrix(i_to_v, images);
  CheckConstantImage(i, 0.0, "Stokes I");
  CheckConstantImage(q, 0.0, "Stokes Q");
  CheckConstantImage(u, 0.0, "Stokes U");
  // V = -i(XY - YX)/2 = -i((10i) - (-10i)) = (10 + 10) / 2 = 10
  CheckConstantImage(v, 10.0, "Stokes V");
}

BOOST_AUTO_TEST_SUITE_END()
