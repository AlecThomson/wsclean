#include "../../msproviders/reordering.h"
#include "aocommon/polarization.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <tuple>
#include <vector>
#include <complex>
#include <cstddef>

using aocommon::Polarization;
using aocommon::PolarizationEnum;
using wsclean::reordering::ExtractData;
using wsclean::reordering::ExtractWeights;

using ComplexVector = std::vector<std::complex<float>>;

const ComplexVector test_data{
    10.0f, 11.0f, 12.0f, 13.0f, 20.0f, 21.0f, 22.0f, 23.0f,
    30.0f, 31.0f, 32.0f, 33.0f, 40.0f, 41.0f, 42.0f, 43.0f,
};
const std::vector<float> test_weights{0.1f, 0.1f, 0.1f, 0.1f, 0.2f, 0.2f,
                                      0.2f, 0.2f, 0.3f, 0.3f, 0.3f, 0.3f,
                                      0.4f, 0.4f, 0.4f, 0.4f};
// std::vector<bool> does not have a .data() method
const bool test_flags[] = {false, false, false, false, false, false,
                           false, false, true,  true,  true,  true,
                           false, false, false, false};

BOOST_AUTO_TEST_SUITE(reorder_extraction)

BOOST_DATA_TEST_CASE(extract_data_linear_pol_to_stokes,
                     boost::unit_test::data::make({Polarization::StokesI,
                                                   Polarization::StokesQ,
                                                   Polarization::StokesU,
                                                   Polarization::StokesV}),
                     pol_out) {
  constexpr size_t kNChannel = 4;
  std::vector<std::complex<float>> actual(kNChannel);
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::StokesI, ComplexVector{11.5f, 21.5f, 31.5f, 41.5f}},
      {Polarization::StokesQ, ComplexVector{-1.5f, -1.5f, -1.5f, -1.5f}},
      {Polarization::StokesU, ComplexVector{11.5f, 21.5f, 31.5f, 41.5f}},
      {Polarization::StokesV,
       ComplexVector{{0.0f, 0.5f}, {0.0f, 0.5f}, {0.0f, 0.5f}, {0.0f, 0.5f}}},
  };
  const ComplexVector expected_result_for_pol = expected.find(pol_out)->second;

  ExtractData(actual.data(), 0, kNChannel,
              std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                         Polarization::YX, Polarization::YY},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_DATA_TEST_CASE(extract_data_circular_pol_to_stokes,
                     boost::unit_test::data::make({Polarization::StokesI,
                                                   Polarization::StokesQ,
                                                   Polarization::StokesU,
                                                   Polarization::StokesV}),
                     pol_out) {
  constexpr size_t kNChannel = 4;
  std::vector<std::complex<float>> actual(kNChannel);
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::StokesI, ComplexVector{11.5f, 21.5f, 31.5f, 41.5f}},
      {Polarization::StokesQ, ComplexVector{11.5f, 21.5f, 31.5f, 41.5f}},
      {Polarization::StokesU,
       ComplexVector{{0.0f, 0.5}, {0.0f, 0.5f}, {0.0f, 0.5f}, {0.0f, 0.5f}}},
      {Polarization::StokesV, ComplexVector{-1.5f, -1.5f, -1.5f, -1.5f}},
  };
  const ComplexVector expected_result_for_pol = expected.find(pol_out)->second;

  ExtractData(actual.data(), 0, kNChannel,
              std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                         Polarization::LR, Polarization::LL},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_DATA_TEST_CASE(
    extract_data_linear_pol_to_linear,
    boost::unit_test::data::make({Polarization::XX, Polarization::XY,
                                  Polarization::YX, Polarization::YY}),
    pol_out) {
  constexpr size_t kNChannel = 4;
  std::vector<std::complex<float>> actual(kNChannel);
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::XX, ComplexVector{10.0f, 20.0f, 30.0f, 40}},
      {Polarization::XY, ComplexVector{11.0f, 21.0f, 31.0f, 41}},
      {Polarization::YX, ComplexVector{12.0f, 22.0f, 32.0f, 42}},
      {Polarization::YY, ComplexVector{13.0f, 23.0f, 33.0f, 43.0f}},
  };
  const ComplexVector expected_result_for_pol = expected.find(pol_out)->second;

  ExtractData(actual.data(), 0, kNChannel,
              std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                         Polarization::YX, Polarization::YY},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_DATA_TEST_CASE(
    extract_data_circular_pol_to_circular,
    boost::unit_test::data::make({Polarization::RR, Polarization::RL,
                                  Polarization::LR, Polarization::LL}),
    pol_out) {
  constexpr size_t kNChannel = 4;
  std::vector<std::complex<float>> actual(kNChannel);
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::RR, ComplexVector{10.0f, 20.0f, 30.0f, 40}},
      {Polarization::RL, ComplexVector{11.0f, 21.0f, 31.0f, 41}},
      {Polarization::LR, ComplexVector{12.0f, 22.0f, 32.0f, 42}},
      {Polarization::LL,
       ComplexVector{
           13.0f,
           23.0f,
           33.0f,
           43.0f,
       }},
  };
  const ComplexVector expected_result_for_pol = expected.find(pol_out)->second;

  ExtractData(actual.data(), 0, kNChannel,
              std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                         Polarization::LR, Polarization::LL},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_AUTO_TEST_CASE(extract_data_linear_pol_to_instrumental) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 4;
  std::vector<std::complex<float>> actual(kNChannel * kNPolPerFile);
  const ComplexVector expected(test_data);
  constexpr PolarizationEnum pol_out = Polarization::Instrumental;

  ExtractData(actual.data(), 0, kNChannel,
              std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                         Polarization::YX, Polarization::YY},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extract_data_circular_pol_to_instrumental) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 4;
  std::vector<std::complex<float>> actual(kNChannel * kNPolPerFile);
  const ComplexVector expected(test_data);
  constexpr PolarizationEnum pol_out = Polarization::Instrumental;

  ExtractData(actual.data(), 0, kNChannel,
              std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                         Polarization::LR, Polarization::LL},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extract_data_linear_pol_to_diag_instrumental) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 2;
  std::vector<std::complex<float>> actual(kNChannel * kNPolPerFile);
  const ComplexVector expected{
      10.0f, 13.0f, 20.0f, 23.0f, 30.0f, 33.0f, 40.0f, 43.0f,
  };
  constexpr PolarizationEnum pol_out = Polarization::DiagonalInstrumental;

  ExtractData(actual.data(), 0, kNChannel,
              std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                         Polarization::YX, Polarization::YY},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extract_data_circular_pol_to_diag_instrumental) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 2;
  std::vector<std::complex<float>> actual(kNChannel * kNPolPerFile);
  const ComplexVector expected{
      10.0f, 13.0f, 20.0f, 23.0f, 30.0f, 33.0f, 40.0f, 43.0f,
  };
  constexpr PolarizationEnum pol_out = Polarization::DiagonalInstrumental;

  ExtractData(actual.data(), 0, kNChannel,
              std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                         Polarization::LR, Polarization::LL},
              test_data.data(), pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_DATA_TEST_CASE(extract_weights_linear_pol_to_stokes,
                     boost::unit_test::data::make({Polarization::StokesI,
                                                   Polarization::StokesQ,
                                                   Polarization::StokesU,
                                                   Polarization::StokesV}),
                     pol_out) {
  constexpr size_t kNChannel = 4;
  std::vector<float> actual(kNChannel);
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, std::vector<float>> expected{
      {Polarization::StokesI, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesQ, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesU, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesV, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
  };

  const std::vector<float> expected_result_for_pol =
      expected.find(pol_out)->second;

  ExtractWeights(actual.data(), 0, kNChannel,
                 std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                            Polarization::YX, Polarization::YY},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_DATA_TEST_CASE(extract_weights_circular_pol_to_stokes,
                     boost::unit_test::data::make({Polarization::StokesI,
                                                   Polarization::StokesQ,
                                                   Polarization::StokesU,
                                                   Polarization::StokesV}),
                     pol_out) {
  constexpr size_t kNChannel = 4;
  std::vector<float> actual(kNChannel);
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, std::vector<float>> expected{
      {Polarization::StokesI, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesQ, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesU, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
      {Polarization::StokesV, std::vector<float>{0.4f, 0.8f, 0.0f, 1.6f}},
  };
  const std::vector<float> expected_result_for_pol =
      expected.find(pol_out)->second;

  ExtractWeights(actual.data(), 0, kNChannel,
                 std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                            Polarization::LR, Polarization::LL},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_DATA_TEST_CASE(
    extract_weights_linear_pol_to_linear,
    boost::unit_test::data::make({Polarization::XX, Polarization::XY,
                                  Polarization::YX, Polarization::YY}),
    pol_out) {
  constexpr size_t kNChannel = 4;
  std::vector<float> actual(kNChannel);
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, std::vector<float>> expected{
      {Polarization::XX, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::XY, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::YX, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::YY, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
  };
  const std::vector<float> expected_result_for_pol =
      expected.find(pol_out)->second;

  ExtractWeights(actual.data(), 0, kNChannel,
                 std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                            Polarization::YX, Polarization::YY},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_DATA_TEST_CASE(
    extract_weights_circular_pol_to_circular,
    boost::unit_test::data::make({Polarization::RR, Polarization::RL,
                                  Polarization::LR, Polarization::LL}),
    pol_out) {
  constexpr size_t kNChannel = 4;
  std::vector<float> actual(kNChannel);
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, std::vector<float>> expected{
      {Polarization::RR, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::RL, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::LR, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
      {Polarization::LL, std::vector<float>{0.1f, 0.2f, 0.0f, 0.4f}},
  };
  const std::vector<float> expected_result_for_pol =
      expected.find(pol_out)->second;

  ExtractWeights(actual.data(), 0, kNChannel,
                 std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                            Polarization::LR, Polarization::LL},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_AUTO_TEST_CASE(extract_weights_linear_pol_to_instrumental) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 4;
  std::vector<float> actual(kNChannel * kNPolPerFile);
  std::vector<float> expected(test_weights);
  size_t idx = 0;
  for (float& weight : expected) {
    weight *= 4;
    // Channel 3 is flagged
    weight *= idx >= 8 && idx < 12 ? 0.0f : 1.0f;
    idx++;
  }

  constexpr PolarizationEnum pol_out = Polarization::Instrumental;

  ExtractWeights(actual.data(), 0, kNChannel,
                 std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                            Polarization::YX, Polarization::YY},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extract_weights_circular_pol_to_instrumental) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 4;
  std::vector<float> actual(kNChannel * kNPolPerFile);
  std::vector<float> expected(test_weights);
  size_t idx = 0;
  for (float& weight : expected) {
    weight *= 4;
    // Channel 3 is flagged
    weight *= idx >= 8 && idx < 12 ? 0.0f : 1.0f;
    idx++;
  }

  constexpr PolarizationEnum pol_out = Polarization::Instrumental;

  ExtractWeights(actual.data(), 0, kNChannel,
                 std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                            Polarization::LR, Polarization::LL},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extract_weight_linear_pol_to_diag_instrumental) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 2;
  std::vector<float> actual(kNChannel * kNPolPerFile);
  const std::vector<float> expected{0.4f, 0.4f, 0.8f, 0.8f,
                                    0.0f, 0.0f, 1.6f, 1.6f};
  constexpr PolarizationEnum pol_out = Polarization::DiagonalInstrumental;

  ExtractWeights(actual.data(), 0, kNChannel,
                 std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                            Polarization::YX, Polarization::YY},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(extract_weight_circular_pol_to_diag_instrumental) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 2;
  std::vector<float> actual(kNChannel * kNPolPerFile);
  const std::vector<float> expected{0.4f, 0.4f, 0.8f, 0.8f,
                                    0.0f, 0.0f, 1.6f, 1.6f};
  constexpr PolarizationEnum pol_out = Polarization::DiagonalInstrumental;

  ExtractWeights(actual.data(), 0, kNChannel,
                 std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                            Polarization::LR, Polarization::LL},
                 test_data.data(), test_weights.data(), test_flags, pol_out);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_SUITE_END()
