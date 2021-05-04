#include <gtest/gtest.h>
#include <ode/ode.h>

namespace {

template <typename T>
bool cmpNum(T n1, T n2) {
  auto diff = std::fabs(n1 - n2);
  if (diff < 1e-5)
    return true;
  else
    printf("Diff: %e\n", diff);
  return false;
}

#define HEADER
#define ShowMat0f(MAT, X, Y) printMatrix0f(#MAT, MAT, (X), (Y))
#define ShowMat(MAT, X, Y) printMatrixf(#MAT, MAT, (X), (Y))
#define ShowV3(MAT) ShowMat(MAT, 3, 1)
#define assert ASSERT_TRUE
#define ShowPMat(MAT, X, Y) printPaddedMatrixf(#MAT, MAT, (X), (Y))

// test miscellaneous math functions
TEST(RandomNumberGenerator, 0) { ASSERT_TRUE(dTestRand()); }

TEST(testMatrixMultiply, 0) {
  // A is 2x3, B is 3x4, B2 is B except stored columnwise, C is 2x4
  int i;
  dReal host_A[8];
  dReal host_B[12];
  dReal host_A2[12];
  dReal host_B2[16];
  dReal host_C[8];
  dReal ode_C[8];

  dReal *C = cuda_makeOnDevice(8);

  HEADER;
  dSetZero(host_A, 8);
  for (i = 0; i < 3; i++) host_A[i] = i + 2;
  for (i = 0; i < 3; i++) host_A[i + 4] = i + 3 + 2;
  dReal *A = cuda_copyToDevice(host_A, 8);
  for (i = 0; i < 12; i++) host_B[i] = i + 8;
  dReal *B = cuda_copyToDevice(host_B, 12);
  dSetZero(host_A2, 12);
  for (i = 0; i < 6; i++) host_A2[i + 2 * (i / 2)] = host_A[i + i / 3];
  dReal *A2 = cuda_copyToDevice(host_A2, 12);
  dSetZero(host_B2, 16);
  for (i = 0; i < 12; i++) host_B2[i + i / 3] = host_B[i];
  dReal *B2 = cuda_copyToDevice(host_B2, 16);

  ShowMat0f(host_A, 2, 3);
  ShowMat0f(host_B, 3, 4);
// #define printMatrix printMatrix0f
#define printMatrix printMatrix0f
  dMultiply0(ode_C, host_A, host_B, 2, 3, 4, false);
  printMatrix("C_ODE", ode_C, 2, 4);

  cuda_dMultiply0(C, A, B, 2, 3, 4);
  cuda_copyFromDevice(C, host_C, 8);
  printMatrix("C_CUD", host_C, 2, 4);
  for (i = 0; i < 8; ++i) {
    ASSERT_TRUE(cmpNum(ode_C[i], host_C[i]));
  }
  //   printf ("\tpassed (1)\n");

  printMatrix("A2", host_A2, 3, 2);
  printMatrix("B", host_B, 3, 4);

  dMultiply1(ode_C, host_A2, host_B, 2, 3, 4, false);
  printMatrix("C_ODE", ode_C, 2, 4);
  if (ode_C[0] != 160 || ode_C[1] != 172 || ode_C[2] != 184 ||
      ode_C[3] != 196 || ode_C[4] != 196 || ode_C[5] != 211 ||
      ode_C[6] != 226 || ode_C[7] != 241)
    printf("\tFAILED (2)\n");
  else
    printf("\tpassed (2)\n");

  cuda_dMultiply1(C, A2, B, 2, 3, 4);
  cuda_copyFromDevice(C, host_C, 8);
  printMatrix("C_CUD", host_C, 2, 4);
  for (i = 0; i < 8; ++i) {
    assert(cmpNum(ode_C[i], host_C[i]));
  }
  printf("\tpassed (2)\n");

  dMultiply2(ode_C, host_A, host_B2, 2, 3, 4, false);
  printMatrix("C_ODE", ode_C, 2, 4);

  cuda_dMultiply2(C, A, B2, 2, 3, 4);
  cuda_copyFromDevice(C, host_C, 8);
  printMatrix("C_CUD", host_C, 2, 4);
  for (i = 0; i < 8; ++i) {
    assert(cmpNum(ode_C[i], host_C[i]));
  }
  printf("\tpassed (3)\n");
#undef printMatrix
  cuda_freeFromDevice(A);
  cuda_freeFromDevice(B);
  cuda_freeFromDevice(A2);
  cuda_freeFromDevice(B2);
  cuda_freeFromDevice(C);
}

TEST(testMatrixMultiply, 3) {
  const int dim = 3;
  using Mat3 = dReal[dPAD(dim) * dPAD(dim)];
  Mat3 A, B, C, host_C;
  const dReal *C_on_CUDA = host_C, *C_ODE = C;
  for (int i = 0; i < dim * dim; ++i) {
    A[i + i / dim] = i;
  }
  makeIdMatrix(B, dim + 1, 1);
  dSetZero(C, (dim + 1) * (dim + 1));

  ShowPMat(A, 3, 3);
  ShowPMat(B, 3, 3);
  ShowPMat(C, 3, 3);

  dReal *dev_A = cuda_copyPaddedToDevice(A, dim);
  dReal *dev_B = cuda_copyPaddedToDevice(B, dim);
  dReal *dev_C = cuda_copyPaddedToDevice(C, dim);

  // with ODE
  // C = A * B
  dMultiply0(C, A, B, dim, dim, dim);
  ShowPMat(C, 3, 3);

  // with CUDA
  // dev_C = dev_A * dev_B
  cuda_dMultiply0(dev_C, dev_A, dev_B, dim, dim, dim);
  cuda_copyFromDevice(dev_C, host_C, dim * dim);
  ShowMat(host_C, 3, 3);
  cuda_freeFromDevice(dev_A);
  cuda_freeFromDevice(dev_B);
  cuda_freeFromDevice(dev_C);

  for (int i = 0; i < dim * dim; i++) {
    assert(C[(i / dim) * dPAD(dim) + i % dim] == host_C[i]);
  }
}

TEST(testMatrixMultiply, 4) {
  /* We're computing C = (A^T)*B
   * row and col are the height and width of A, so C is a col*col matrix.
   * Thus row corresponds to q, col to p and r */
  int p = 5;
  int q = 3;
  int r = 7;
  int row = q;
  int col = 2;
  dReal A[q * p], B[q * r], C[p * r], host_C[p * r];
  for (int i = 0; i < q * p; ++i) {
    A[i] = i + 1;
  }
  for (int i = 0; i < q * r; ++i) {
    B[i] = i + 1;
  }
  ShowMat0f(A, q, p);
  ShowMat0f(B, q, r);
  dSetZero(C, p * r);
  ShowMat0f(C, p, r);
  dReal *dev_A = cuda_copyToDevice(A, q * p);
  dReal *dev_B = cuda_copyToDevice(B, q * r);
  dReal *dev_C = cuda_copyToDevice(C, p * r);
  dMultiply2(C, A, B, p, q, r, false);
  printMatrix0f("C_ODE", C, (p), (r));
  cuda_dMultiply2(dev_C, dev_A, dev_B, p, q, r);
  cuda_copyFromDevice(dev_C, host_C, p * r);
  printMatrix0f("C_CUD", host_C, p, r);
  cuda_freeFromDevice(dev_A);
  cuda_freeFromDevice(dev_B);
  cuda_freeFromDevice(dev_C);
  for (int i = 0; i < p * r; i++) {
    assert(C[i] == host_C[i]);
  }
}

TEST(testMatrixMultiply, 5) {
  int row = 5;
  int col = 2;
  dReal A[9 * 17], B[17 * 25], C[9 * 25], host_C[9 * 25];
  for (int i = 0; i < 9 * 17; ++i) {
    A[i] = i + 1;
  }
  printMatrix0f("A", A, 9, 17);
  for (int i = 0; i < 17 * 25; ++i) {
    B[i] = i + 1;
  }
  printMatrix0f("B", B, 17, 25);
  dSetZero(C, 9 * 25);
  printMatrix0f("C", C, 9, 25);
  dReal *dev_A = cuda_copyToDevice(A, 9 * 17);
  dReal *dev_B = cuda_copyToDevice(B, 17 * 25);
  dReal *dev_C = cuda_realMalloc(9 * 25);
  dMultiply0(C, A, B, 9, 17, 25, false);
  ShowMat0f(C, 9, 25);
  cuda_dMultiply0(dev_C, dev_A, dev_B, 9, 17, 25);
  cuda_copyFromDevice(dev_C, host_C, 9 * 25);
  ShowMat0f(host_C, 9, 25);
  cuda_freeFromDevice(dev_A);
  cuda_freeFromDevice(dev_B);
  cuda_freeFromDevice(dev_C);
  for (int i = 0; i < 9 * 25; i++) {
    assert(cmpNum(C[i], host_C[i]));
  }
}

#define fat_matrix_with_dim(dim)             \
  TEST(fat_matrix, dim) {                    \
    dReal *A = cuda_makeOnDevice(dim);       \
    dReal *B = cuda_makeOnDevice(dim);       \
    dReal *C = cuda_makeOnDevice(dim);       \
    cuda_dSetValue(A, 1, dim *dim);          \
    cuda_dSetValue(B, 2, dim *dim);          \
    cuda_dSetZero(C, dim *dim);              \
    cuda_dMultiply0(C, A, B, dim, dim, dim); \
    cuda_dSetValue(C, 1, dim *dim);          \
    cuda_dMultiply1(C, A, B, dim, dim, dim); \
    cuda_dSetValue(C, 2, dim *dim);          \
    cuda_dMultiply2(C, A, B, dim, dim, dim); \
    cuda_dSetValue(C, 3, dim *dim);          \
    cuda_freeFromDevice(A);                  \
    cuda_freeFromDevice(B);                  \
    cuda_freeFromDevice(C);                  \
  }

fat_matrix_with_dim(100);
fat_matrix_with_dim(200);
fat_matrix_with_dim(400);

}  // namespace
