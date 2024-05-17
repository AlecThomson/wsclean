#include "averagecorrection.h"

#include <complex>

extern "C" void zheev_(const char* jobz, const char* uplo, const int* n,
                       std::complex<double>* a, const int* lda, double* w,
                       std::complex<double>* work, const int* lwork,
                       double* rwork, int* info);

aocommon::HMC4x4 PrincipalSquareRoot(const aocommon::HMC4x4& matrix) {
  double ev[4];
  constexpr int n = 4;
  constexpr int work_size = 20;
  std::complex<double> work[work_size];
  double rwork[std::max(1, 3 * n - 2)];  // size as required by zheev
  int info = 0;
  std::complex<double> a[16];
  for (size_t col = 0; col != 4; ++col) {
    for (size_t row = 0; row != 4; ++row) {
      // LAPACK uses col-first, HMC4x4 uses row first, so translate:
      a[col * 4 + row] = matrix[row * 4 + col];
    }
  }
  const char job_mode = 'V';        // Get eigenvalues and eigenvectors
  const char upper_or_lower = 'L';  // Lower triangle of A is stored
  // CHEEV computes all eigenvalues and, optionally, eigenvectors of a
  // complex Hermitian matrix.
  zheev_(&job_mode, &upper_or_lower,
         &n,  // Order of A
         a,
         &n,  // leading dimension of the array A
         ev, work, &work_size, rwork, &info);
  if (info == 0) {
    bool is_positive_semi_definite = true;
    for (size_t i = 0; i != 4; ++i) {
      if (ev[i] < 0.0) {
        is_positive_semi_definite = false;
        break;
      }
      ev[i] = std::sqrt(ev[i]);
    }

    if (is_positive_semi_definite) {
      const auto c = [](std::complex<double> z) -> std::complex<double> {
        return std::conj(z);
      };
      // Recompose the full matrix: Calculate M = A Λ^1/2 A^H
      // where A is the 4 × 4 matrix whose ith column is the eigenvector qi of
      // the matrix. Λ is the diagonal matrix whose diagonal elements are the
      // corresponding eigenvalues, Λii = λi.
      //
      // Note that LAPACK uses row-first ordering, HMC4x4 uses col first.
      // Therefore, the LAPACK results are translated. Because the last term
      // (A^H) was already translated, it is translated twice and thus does not
      // need to be translated.
      return aocommon::HMC4x4{
          // Row 1
          std::norm(a[0]) * ev[0] + std::norm(a[4]) * ev[1] +
              std::norm(a[8]) * ev[2] + std::norm(a[12]) * ev[3],
          0, 0, 0,
          // Row 2
          a[1] * ev[0] * c(a[0]) + a[5] * ev[1] * c(a[4]) +
              a[9] * ev[2] * c(a[8]) + a[13] * ev[3] * c(a[12]),
          std::norm(a[1]) * ev[0] + std::norm(a[5]) * ev[1] +
              std::norm(a[9]) * ev[2] + std::norm(a[13]) * ev[3],
          0, 0,
          // Row 3
          a[2] * ev[0] * c(a[0]) + a[6] * ev[1] * c(a[4]) +
              a[10] * ev[2] * c(a[8]) + a[14] * ev[3] * c(a[12]),
          a[2] * ev[0] * c(a[1]) + a[6] * ev[1] * c(a[5]) +
              a[10] * ev[2] * c(a[9]) + a[14] * ev[3] * c(a[13]),
          std::norm(a[2]) * ev[0] + std::norm(a[6]) * ev[1] +
              std::norm(a[10]) * ev[2] + std::norm(a[14]) * ev[3],
          0,
          // Row 4
          a[3] * ev[0] * c(a[0]) + a[7] * ev[1] * c(a[4]) +
              a[11] * ev[2] * c(a[8]) + a[15] * ev[3] * c(a[12]),
          a[3] * ev[0] * c(a[1]) + a[7] * ev[1] * c(a[5]) +
              a[11] * ev[2] * c(a[9]) + a[15] * ev[3] * c(a[13]),
          a[3] * ev[0] * c(a[2]) + a[7] * ev[1] * c(a[6]) +
              a[11] * ev[2] * c(a[10]) + a[15] * ev[3] * c(a[14]),
          std::norm(a[3]) * ev[0] + std::norm(a[7]) * ev[1] +
              std::norm(a[11]) * ev[2] + std::norm(a[15]) * ev[3]};
    }
  }
  // Return a matrix with NaNs on the diagonal.
  return aocommon::HMC4x4::Unit() * std::numeric_limits<double>::quiet_NaN();
}

std::string ToString(const AverageCorrection& average_correction) {
  std::ostringstream str;
  str << "AverageCorrection, ";
  if (average_correction.IsScalar()) {
    str << "scalar: " << average_correction.GetScalarValue();
  } else {
    str << "matrix:\n" << average_correction.GetMatrixValue().String();
  }
  return str.str();
}
