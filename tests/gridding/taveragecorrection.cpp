#include <boost/test/unit_test.hpp>

#include <aocommon/hmatrix4x4.h>

#include "../../gridding/averagecorrection.h"

#include <cmath>
#include <complex>
#include <iostream>

using aocommon::HMC4x4;

namespace wsclean {

void CheckMatrix(const HMC4x4& m, const HMC4x4& reference) {
  for (size_t i = 0; i != 16; ++i) {
    std::ostringstream a;
    std::ostringstream b;
    a << "index " << i << ", diff "
      << std::round(std::abs(m[i] - reference[i]) * 100.0) / 100.0 << ", got "
      << m[i] << ", expected " << reference[i];
    b << "index " << i << ", diff 0, got " << m[i] << ", expected "
      << reference[i];
    BOOST_CHECK_EQUAL(a.str(), b.str());
  }
}

BOOST_AUTO_TEST_SUITE(average_correction)

BOOST_AUTO_TEST_CASE(square_root_unit) {
  const HMC4x4 result = PrincipalSquareRoot(HMC4x4::Unit());
  CheckMatrix(result, HMC4x4::Unit());
}

BOOST_AUTO_TEST_CASE(square_root_scalar) {
  const HMC4x4 result = PrincipalSquareRoot(HMC4x4::Unit() * 4.0);
  CheckMatrix(result, HMC4x4::Unit() * 2.0);
}

BOOST_AUTO_TEST_CASE(square_root_diagonal) {
  const HMC4x4 input{1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 4, 0, 0, 0, 0, 16};
  const HMC4x4 result = PrincipalSquareRoot(input.Square());
  CheckMatrix(result, input);
}

BOOST_AUTO_TEST_CASE(square_root_zero_eigenvalue) {
  const std::complex<double> j(0, 1);
  const HMC4x4 input{5, 2.0 + j, 0, 0, 2.0 - j, 5, 0, 0,
                     0, 0,       0, 0, 0,       0, 0, 0};
  const HMC4x4 result = PrincipalSquareRoot(input);
  CheckMatrix(result.Square(), input);
}

BOOST_AUTO_TEST_CASE(square_root_complex) {
  const std::complex<double> j(0, 1);
  const HMC4x4 complex_matrix{2,  j,  j, j, -j, 2,  j, j,
                              -j, -j, 2, 1, -j, -j, 1, 2};
  const HMC4x4 result = PrincipalSquareRoot(complex_matrix);
  CheckMatrix(result.Square(), complex_matrix);
}

BOOST_AUTO_TEST_CASE(square_root_non_positive_semi_definite) {
  const std::complex<double> j(0, 1);
  const HMC4x4 complex_matrix{2,  j,  j, j, -j, 2,  j,  j,
                              -j, -j, 2, j, -j, -j, -j, 2};
  const HMC4x4 result = PrincipalSquareRoot(complex_matrix);
  BOOST_CHECK(!std::isfinite(result[0].real()));
  BOOST_CHECK(!std::isfinite(result[5].real()));
  BOOST_CHECK(!std::isfinite(result[10].real()));
  BOOST_CHECK(!std::isfinite(result[15].real()));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace wsclean
