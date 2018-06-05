// Distributed under the MIT License.
// See LICENSE.txt for details.
// Code written by Geoffrey Lovelace (glovelace@fullerton.edu)

/* Test lapack's dggev call to find the generalized eigenvectors and
 * eigenvalues given two square matrices A and B.
 *
 * To compile on macOS, you can use
 *  i) blas/lapack in the Accelerate framework
 *   - In terminal, compile with
 * clang++ -framework Accelerate -o Test_dggev_Accelerate Test_dggev.cpp
 *
 *  ii) openblas
 *   - Install homebrew (https://brew.sh)
 *   - In terminal, do "brew install openblas"
 *   - In terminal, compile with
 * clang++ -L/usr/local/opt/openblas/lib -lopenblas -o Test_dggev Test_dggev.cpp
 */

#include <iostream>
#include <xmmintrin.h>

// This header file includes some hard-coded matrices for the test,
// defined in a header to keep the main code clean.
// Since I exported these matrices in csv format, it was easiest to
// turn them into C++ code for row-major arrays. I will later
// convert to column-major arrays as lapack expects.
#include "Matrices.hpp"

using namespace std;

// The test matrices are 168x168
const int matrix_dimension = 168;

// LAPACK routine to do the generalized eigenvalue problem
extern "C" {
extern void dggev_(char*, char*, int*, double*, int*, double*, int*, double*,
                   double*, double*, double*, int*, double*, int*, double*,
                   int*, int*);
}

int main(void) {
  // Enable floating-point exceptions on overflow, invalid operation, and
  // divide by zero
  _mm_setcsr(_MM_MASK_MASK &
             ~(_MM_MASK_OVERFLOW | _MM_MASK_INVALID | _MM_MASK_DIV_ZERO));

  // Construct column-major order matrices A and B from the
  // row-major matrices in Matrices.hpp
  double matrix_A[matrix_dimension * matrix_dimension];
  double matrix_B[matrix_dimension * matrix_dimension];
  for(int row = 0; row < matrix_dimension; ++row) {
    for(int col = 0; col < matrix_dimension; ++col) {
      matrix_A[matrix_dimension * col + row]
        = matrix_A_rowmajor[matrix_dimension * row + col];
      matrix_B[matrix_dimension * col + row]
        = matrix_B_rowmajor[matrix_dimension * row + col];
    }
  }
  
  // Set up the lapack call
  char compute_left_eigenvectors = 'N'; // 'N' = don't compute
  char compute_right_eigenvectors = 'V'; // 'V' = do compute
  int info = 0; // info = 0 after dggev call: success.
  
  // Lapack expects an int for the matrix size
  int matrix_and_vector_size = static_cast<int>(matrix_dimension);
  
  // Arrays to store the results.
  // Lapack splits the eigenvalues into unnormalized real and imaginary
  // parts, which it calls alphar and alphai, and a normalization,
  // which it calls beta. The real and imaginary parts of the eigenvalues are
  // found by dividing the unnormalized results by the normalization.
  double eigenvalues_real_part[matrix_dimension];
  double eigenvalues_im_part[matrix_dimension];
  double eigenvalue_normalization[matrix_dimension];
  double eigenvectors[matrix_dimension * matrix_dimension];
  for(int i = 0; i < matrix_dimension; ++i) {
    eigenvalues_real_part[i] = 0.0;
    eigenvalues_im_part[i] = 0.0;
    eigenvalue_normalization[i] = 0.0;
    for(int j = 0; j < matrix_dimension; ++j) {
      eigenvectors[j*matrix_dimension + i] = 0.0;
    }
  }
  
  // Lapack uses this array for scratch work
  int work_size = matrix_and_vector_size * 8;
  double lapack_work[work_size];
  for(int i = 0; i < work_size; ++i) {
    lapack_work[i] = 0.0;
  }
  
  // Call lapack's dggev to get the generalized eigenvectors and eigenvalues
  dggev_(&compute_left_eigenvectors, &compute_right_eigenvectors,
         &matrix_and_vector_size, matrix_A, &matrix_and_vector_size,
         matrix_B, &matrix_and_vector_size,
         eigenvalues_real_part, eigenvalues_im_part,
         eigenvalue_normalization, eigenvectors,
         &matrix_and_vector_size, eigenvectors, &matrix_and_vector_size,
         lapack_work, &work_size, &info);
  
  // Print the real part of the first eigenvalue
  cout << "Real_eigenvalue(" << 0 << ") = "
       << eigenvalues_real_part[0]/eigenvalue_normalization[0] << '\n';
  
  return 0;
}
