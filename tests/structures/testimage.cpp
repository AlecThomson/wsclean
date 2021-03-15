#include <boost/test/unit_test.hpp>

#include "../../structures/image.h"

#include <aocommon/uvector.h>

#include <random>

BOOST_AUTO_TEST_SUITE(image_operations)

BOOST_AUTO_TEST_CASE(median_empty) {
  BOOST_CHECK_EQUAL(DImage::Median(nullptr, 0), 0.0);
  BOOST_CHECK_EQUAL(ImageF::Median(nullptr, 0), 0.0);
}

BOOST_AUTO_TEST_CASE(median_single) {
  aocommon::UVector<double> arr(1, 1.0);
  BOOST_CHECK_EQUAL(DImage::Median(arr.data(), arr.size()), 1.0);

  arr[0] = std::numeric_limits<double>::quiet_NaN();
  DImage::Median(arr.data(),
                 arr.size());  // undefined -- just make sure it doesn't crash
}

BOOST_AUTO_TEST_CASE(median_two_elements) {
  aocommon::UVector<double> arr1(2, 1.0);
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr1.data(), arr1.size()), 1.0,
                             1e-5);

  aocommon::UVector<double> arr2(2, 0.0);
  arr2[1] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr2.data(), arr2.size()), 1.0,
                             1e-5);

  aocommon::UVector<double> arr3(2, 1.0);
  arr3[1] = -1.0;
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr3.data(), arr3.size()), 0.0,
                             1e-5);

  aocommon::UVector<double> arr4(2, 13.0);
  arr3[1] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr4.data(), arr4.size()), 13.0,
                             1e-5);
}

BOOST_AUTO_TEST_CASE(median_three_elements) {
  aocommon::UVector<double> arr(3, 1.0);
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr.data(), arr.size()), 1.0, 1e-5);

  arr[0] = 0.0;
  arr[1] = 1.0;
  arr[2] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr.data(), arr.size()), 1.0, 1e-5);

  arr[0] = 3.0;
  arr[1] = -3.0;
  arr[2] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr.data(), arr.size()), 2.0, 1e-5);

  arr[1] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr.data(), arr.size()), 2.5, 1e-5);

  arr[0] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr.data(), arr.size()), 2.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(mad_empty) {
  BOOST_CHECK_EQUAL(DImage::MAD(nullptr, 0), 0.0);
}

BOOST_AUTO_TEST_CASE(mad_single) {
  aocommon::UVector<double> arr(1, 1.0);
  BOOST_CHECK_EQUAL(DImage::MAD(arr.data(), arr.size()), 0.0);
}

BOOST_AUTO_TEST_CASE(mad_two_elements) {
  aocommon::UVector<double> arr(2, 1.0);
  BOOST_CHECK_EQUAL(DImage::MAD(arr.data(), arr.size()), 0.0);

  arr[0] = 0.0;
  arr[1] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr.data(), arr.size()), 1.0, 1e-5);

  arr[0] = 1.0;
  arr[1] = -1.0;
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr.data(), arr.size()), 0.0, 1e-5);

  arr[0] = 13.0;
  arr[1] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr.data(), arr.size()), 13.0,
                             1e-5);
}

BOOST_AUTO_TEST_CASE(mad_three_elements) {
  aocommon::UVector<double> arr(3, 1.0);
  BOOST_CHECK_CLOSE_FRACTION(DImage::MAD(arr.data(), arr.size()), 0.0, 1e-5);

  arr[0] = 0.0;
  arr[1] = 1.0;
  arr[2] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(DImage::MAD(arr.data(), arr.size()), 1.0, 1e-5);

  arr[0] = 3.0;
  arr[1] = -3.0;
  arr[2] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(DImage::MAD(arr.data(), arr.size()), 1.0, 1e-5);

  arr[1] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(DImage::MAD(arr.data(), arr.size()), 0.5, 1e-5);

  arr[0] = std::numeric_limits<double>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(DImage::MAD(arr.data(), arr.size()), 0.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(stddev_from_mad) {
  std::mt19937 rnd;
  std::normal_distribution<double> dist(1.0, 5.0);
  aocommon::UVector<double> data(10000);
  for (size_t i = 0; i != data.size(); ++i) data[i] = dist(rnd);
  BOOST_CHECK_CLOSE_FRACTION(DImage::StdDevFromMAD(data.data(), data.size()),
                             5.0, 0.05);
}

BOOST_AUTO_TEST_SUITE_END()
