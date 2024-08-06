#include "../../msproviders/reordering.h"
#include "aocommon/polarization.h"

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <stdexcept>
#include <tuple>
#include <vector>
#include <complex>
#include <cstddef>

using aocommon::Polarization;
using aocommon::PolarizationEnum;
using wsclean::reordering::StoreData;
using wsclean::reordering::StoreWeights;

using ComplexVector = std::vector<std::complex<float>>;

const ComplexVector test_data{10.0f, 11.0f, 12.0f, 13.0f, 20.0f, 21.0f,
                              22.0f, 23.0f, 30.0f, 31.0f, 32.0f, 33.0f,
                              40.0f, 41.0f, 42.0f, 43.0f};

const std::vector<float> test_weights{0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f,
                                      0.7f, 0.8f, 0.9f, 1.0f, 1.1f, 1.2f,
                                      1.3f, 1.4f, 1.4f, 1.5f};

BOOST_AUTO_TEST_SUITE(reorder_storeion)

BOOST_DATA_TEST_CASE(store_data_stokes_to_linear_pol,
                     boost::unit_test::data::make({
                         Polarization::StokesI,
                         Polarization::StokesQ,
                         Polarization::StokesU,
                         Polarization::StokesV,
                     }),
                     pol_source) {
  constexpr size_t kNChannel = 4;
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, std::set<PolarizationEnum>> pol_dest_mapping{
      {Polarization::StokesI, {Polarization::XX, Polarization::YY}},
      {Polarization::StokesQ, {Polarization::XX, Polarization::YY}},
      {Polarization::StokesU, {Polarization::XY, Polarization::YX}},
      {Polarization::StokesV, {Polarization::XY, Polarization::YX}}};

  const std::set<PolarizationEnum> pol_dest =
      pol_dest_mapping.find(pol_source)->second;

  std::vector<std::complex<float>> actual(kNChannel * pol_dest.size());
  const std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::StokesI,
       ComplexVector{10.0f, 10.0f, 11.0f, 11.0f, 12.0f, 12.0f, 13.0f, 13.0f}},
      {Polarization::StokesQ, ComplexVector{10.0f, -10.0f, 11.0f, -11.0f, 12.0f,
                                            -12.0f, 13.0f, -13.0f}},
      {Polarization::StokesU,
       ComplexVector{10.0f, 10.0f, 11.0f, 11.0f, 12.0f, 12.0f, 13.0f, 13.0f}},

      {
          Polarization::StokesV,
          ComplexVector{{0.0f, 10.0f},
                        {0.0f, -10.0f},
                        {0.0f, 11.0f},
                        {0.0f, -11.0f},
                        {0.0f, 12.0f},
                        {0.0f, -12.0f},
                        {0.0f, 13.0f},
                        {0.0f, -13.0f}},
      },
  };
  const ComplexVector expected_result_for_pol =
      expected.find(pol_source)->second;

  StoreData<false>(actual.data(), 0, kNChannel, pol_dest, test_data.data(),
                   pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_DATA_TEST_CASE(store_data_stokes_to_circular_pol,
                     boost::unit_test::data::make({
                         Polarization::StokesI,
                         Polarization::StokesQ,
                         Polarization::StokesU,
                         Polarization::StokesV,
                     }),
                     pol_source) {
  constexpr size_t kNChannel = 4;
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, std::set<PolarizationEnum>> pol_dest_mapping{
      {Polarization::StokesI, {Polarization::LL, Polarization::RR}},
      {Polarization::StokesQ, {Polarization::LL, Polarization::RR}},
      {Polarization::StokesU, {Polarization::LR, Polarization::RL}},
      {Polarization::StokesV, {Polarization::LR, Polarization::RL}}};

  const std::set<PolarizationEnum> pol_dest =
      pol_dest_mapping.find(pol_source)->second;

  std::vector<std::complex<float>> actual(kNChannel * pol_dest.size());
  const std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::StokesI,
       ComplexVector{10.0f, 10.0f, 11.0f, 11.0f, 12.0f, 12.0f, 13.0f, 13.0f}},
      {Polarization::StokesQ,
       ComplexVector{10.0f, 0.0f, 11.0f, 0.0f, 12.0f, 0.0f, 13.0f, 0.0f}},
      {
          Polarization::StokesU,
          ComplexVector{{0.0f, 10.0f},
                        {0.0f, -10.0f},
                        {0.0f, 11.0f},
                        {0.0f, -11.0f},
                        {0.0f, 12.0f},
                        {0.0f, -12.0f},
                        {0.0f, 13.0f},
                        {0.0f, -13.0f}},
      },
      {
          Polarization::StokesV,
          ComplexVector{-10.0f, 0.0f, -11.0f, 0.0f, -12.0f, 0.0f, -13.0f, 0.0f},
      },
  };
  const ComplexVector expected_result_for_pol =
      expected.find(pol_source)->second;

  StoreData<false>(actual.data(), 0, kNChannel, pol_dest, test_data.data(),
                   pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_DATA_TEST_CASE(store_data_linear_to_linear,
                     boost::unit_test::data::make({
                         Polarization::XX,
                         Polarization::XY,
                         Polarization::YX,
                         Polarization::YY,
                     }),
                     pol_source) {
  constexpr size_t kNChannel = 2;
  constexpr size_t kNPolDest = 4;
  std::vector<std::complex<float>> actual(kNChannel * kNPolDest);
  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::XX,
       ComplexVector{10.0f, 0.0f, 0.0f, 0.0f, 11.0f, 0.0f, 0.0f, 0.0f}},
      {Polarization::XY,
       ComplexVector{0.0f, 10.0f, 0.0f, 0.0f, 0.0f, 11.0f, 0.0f, 0.0f}},
      {Polarization::YX,
       ComplexVector{0.0f, 0.0f, 10.0f, 0.0f, 0.0f, 0.0f, 11.0f, 0.0f}},
      {Polarization::YY,
       ComplexVector{0.0f, 0.0f, 0.0f, 10.0f, 0.0f, 0.0f, 0.0f, 11.0f}},
  };
  const ComplexVector expected_result_for_pol =
      expected.find(pol_source)->second;

  StoreData<false>(
      actual.data(), 0, kNChannel,
      std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                 Polarization::YX, Polarization::YY},
      test_data.data(), pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_DATA_TEST_CASE(store_data_circular_to_circular,
                     boost::unit_test::data::make({
                         Polarization::RR,
                         Polarization::RL,
                         Polarization::LL,
                         Polarization::LR,
                     }),
                     pol_source) {
  constexpr size_t kNChannel = 2;
  constexpr size_t kNPolDest = 4;
  std::vector<std::complex<float>> actual(kNChannel * kNPolDest);

  // Boost auto test case doesn't support containers inside data::make
  const std::map<PolarizationEnum, ComplexVector> expected{
      {Polarization::RR,
       ComplexVector{10.0f, 0.0f, 0.0f, 0.0f, 11.0f, 0.0f, 0.0f, 0.0f}},
      {Polarization::RL,
       ComplexVector{0.0f, 10.0f, 0.0f, 0.0f, 0.0f, 11.0f, 0.0f, 0.0f}},
      {Polarization::LR,
       ComplexVector{0.0f, 0.0f, 10.0f, 0.0f, 0.0f, 0.0f, 11.0f, 0.0f}},
      {Polarization::LL,
       ComplexVector{0.0f, 0.0f, 0.0f, 10.0f, 0.0f, 0.0f, 0.0f, 11.0f}},
  };

  const ComplexVector expected_result_for_pol =
      expected.find(pol_source)->second;

  StoreData<false>(actual.data(), 0, kNChannel,
                   std::set<PolarizationEnum>{
                       Polarization::RR,
                       Polarization::RL,
                       Polarization::LL,
                       Polarization::LR,
                   },
                   test_data.data(), pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_result_for_pol.begin(),
                                expected_result_for_pol.end());
}

BOOST_AUTO_TEST_CASE(store_data_instrumental_to_linear_pol) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 4;
  std::vector<std::complex<float>> actual(kNChannel * kNPolPerFile);
  const ComplexVector expected(test_data);
  const PolarizationEnum pol_source = Polarization::Instrumental;

  StoreData<false>(
      actual.data(), 0, kNChannel,
      std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                 Polarization::YX, Polarization::YY},
      test_data.data(), pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(store_data_instrumental_to_circular_pol) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 4;
  std::vector<std::complex<float>> actual(kNChannel * kNPolPerFile);
  const ComplexVector expected(test_data);
  const PolarizationEnum pol_source = Polarization::Instrumental;

  StoreData<false>(
      actual.data(), 0, kNChannel,
      std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                 Polarization::LR, Polarization::LL},
      test_data.data(), pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(store_data_diag_instrumental_to_linear_pol) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 4;
  std::vector<std::complex<float>> actual(kNChannel * kNPolPerFile);
  const ComplexVector expected{10.0f, 0, 0, 11.0f, 12.0f, 0, 0, 13.0f,
                               20.0f, 0, 0, 21.0f, 22.0f, 0, 0, 23};
  const PolarizationEnum pol_source = Polarization::DiagonalInstrumental;

  StoreData<false>(
      actual.data(), 0, kNChannel,
      std::set<PolarizationEnum>{Polarization::XX, Polarization::XY,
                                 Polarization::YX, Polarization::YY},
      test_data.data(), pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(store_data_diag_instrumental_to_circular_pol) {
  constexpr size_t kNChannel = 4;
  constexpr size_t kNPolPerFile = 4;
  std::vector<std::complex<float>> actual(kNChannel * kNPolPerFile);
  const ComplexVector expected{10.0f, 0, 0, 11.0f, 12.0f, 0, 0, 13.0f,
                               20.0f, 0, 0, 21.0f, 22.0f, 0, 0, 23};
  const PolarizationEnum pol_source = Polarization::DiagonalInstrumental;

  StoreData<false>(
      actual.data(), 0, kNChannel,
      std::set<PolarizationEnum>{Polarization::RR, Polarization::RL,
                                 Polarization::LR, Polarization::LL},
      test_data.data(), pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(store_weight_stokes_i_to_linear) {
  PolarizationEnum pol_source = Polarization::StokesI;
  constexpr size_t kNChannel = 4;
  const std::set<PolarizationEnum> pol_dest{Polarization::XX, Polarization::YY};
  const ComplexVector expected{0.1f, 0.1f, 0.2f, 0.2f, 0.3f, 0.3f, 0.4f, 0.4f};

  std::vector<float> actual(kNChannel * pol_dest.size());
  StoreWeights(actual.data(), 0, kNChannel, pol_dest, test_weights.data(),
               pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(), expected.begin(),
                                expected.end());
}

BOOST_AUTO_TEST_CASE(store_weight_stokes_linear_to_linear) {
  constexpr size_t kNChannel = 4;
  const std::set<PolarizationEnum> pol_dest{Polarization::XX, Polarization::YY};
  const ComplexVector expected_after_xx{0.1f, 0.0f, 0.2f, 0.0f,
                                        0.3f, 0.0f, 0.4f, 0.0f};

  std::vector<float> actual(kNChannel * pol_dest.size());

  // Extract XX
  StoreWeights(actual.data(), 0, kNChannel, pol_dest, test_weights.data(),
               Polarization::XX);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_after_xx.begin(),
                                expected_after_xx.end());

  ComplexVector expected_after_xx_and_yy{0.1f, 0.1f, 0.2f, 0.2f,
                                         0.3f, 0.3f, 0.4f, 0.4f};

  // Extract YY
  StoreWeights(actual.data(), 0, kNChannel, pol_dest, test_weights.data(),
               Polarization::YY);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                expected_after_xx_and_yy.begin(),
                                expected_after_xx_and_yy.end());
}

BOOST_AUTO_TEST_CASE(store_weight_instrumental_to_linear) {
  const PolarizationEnum pol_source = Polarization::Instrumental;
  constexpr size_t kNChannel = 4;
  const std::set<PolarizationEnum> pol_dest{Polarization::XX, Polarization::XY,
                                            Polarization::YX, Polarization::YY};

  std::vector<float> actual(kNChannel * pol_dest.size());
  StoreWeights(actual.data(), 0, kNChannel, pol_dest, test_weights.data(),
               pol_source);
  BOOST_CHECK_EQUAL_COLLECTIONS(actual.begin(), actual.end(),
                                test_weights.begin(), test_weights.end());
}

BOOST_DATA_TEST_CASE(
    store_weight_throws_error_on_invalid,
    boost::unit_test::data::make({Polarization::StokesQ, Polarization::StokesU,
                                  Polarization::StokesV,
                                  Polarization::DiagonalInstrumental}),
    pol_source) {
  constexpr size_t kNChannel = 4;
  const std::set<PolarizationEnum> pol_dest{Polarization::XX, Polarization::YY};
  const ComplexVector expected{0.1f, 0.1f, 0.2f, 0.2f, 0.3f, 0.3f, 0.4f, 0.4f};

  std::vector<float> actual(kNChannel * pol_dest.size());
  ;

  BOOST_CHECK_THROW(StoreWeights(actual.data(), 0, kNChannel, pol_dest,
                                 test_weights.data(), pol_source),
                    std::runtime_error);
}

BOOST_AUTO_TEST_SUITE_END()
