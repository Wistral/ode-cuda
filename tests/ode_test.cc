/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/
#include <gtest/gtest.h>
#include <limits>
#include <type_traits>
// #include <setjmp.h>
#include <ode/ode.h>

#ifdef _MSC_VER
#pragma warning(disable : 4244 4305)  // for VC++, no precision loss complaints
#endif

// build fix defines
#define HEADER
#define dSetDebugHandler(X)
#define setjmp(X)
#define TRAP_MESSAGE(do, ifnomsg, ifmsg) do

//****************************************************************************
// matrix accessors
#define _A(i, j) A[(i)*4 + (j)]
#define _I(i, j) I[(i)*4 + (j)]
#define _R(i, j) R[(i)*4 + (j)]

//****************************************************************************
// tolerances
#ifdef dDOUBLE
const double tol = 1e-10;
#endif
#ifdef dSINGLE
const double tol = 1e-4;
#endif

// USE `IsAlmostEqual()` REPLACE `cmp()`
// Test whether two float or double numbers are equal.
// ulp: units in the last place.
template <typename T, typename T2>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
IsAlmostEqual(T x, T2 y, int ulp = 2, double _EPSILON = 1e-5) {
  // the machine epsilon has to be scaled to the magnitude of the values used
  // and multiplied by the desired precision in ULPs (units in the last place)
  auto cp = std::numeric_limits<T>::epsilon() * std::fabs(x + y) * ulp;
  return std::fabs(x - y) < cp
         // unless the result is subnormal
         || std::fabs(x - y) < std::numeric_limits<T>::min() ||
         std::fabs(x - y) < _EPSILON;
}

#define cmp IsAlmostEqual

//****************************************************************************
// matrix utility stuff

// compare a 3x3 matrix with the identity
int cmpIdentityMat3(dMatrix3 A) {
  return (cmp(_A(0, 0), 1.0) && cmp(_A(0, 1), 0.0) && cmp(_A(0, 2), 0.0) &&
          cmp(_A(1, 0), 0.0) && cmp(_A(1, 1), 1.0) && cmp(_A(1, 2), 0.0) &&
          cmp(_A(2, 0), 0.0) && cmp(_A(2, 1), 0.0) && cmp(_A(2, 2), 1.0));
}

// transpose a 3x3 matrix in-line
void transpose3x3(dMatrix3 A) {
  dReal tmp;
  tmp = A[4];
  A[4] = A[1];
  A[1] = tmp;
  tmp = A[8];
  A[8] = A[2];
  A[2] = tmp;
  tmp = A[9];
  A[9] = A[6];
  A[6] = tmp;
}

//****************************************************************************
// test miscellaneous math functions
TEST(testRandomNumberGenerator, ) {
  HEADER;
  ASSERT_TRUE(dTestRand());
}

TEST(testInfinity, ) {
  HEADER;
  ASSERT_TRUE(1e10 < dInfinity && -1e10 > -dInfinity && -dInfinity < dInfinity);
}

TEST(testPad, ) {
  HEADER;
  char s[100];
  s[0] = 0;
  for (int i = 0; i <= 16; i++) sprintf(s + strlen(s), "%d ", dPAD(i));

  ASSERT_TRUE(!strcmp(s, "0 1 4 4 4 8 8 8 8 12 12 12 12 16 16 16 16 "));
}

TEST(testCrossProduct, ) {
  HEADER;

  dVector3 a1, a2, b, c;
  dMatrix3 B;
  dMakeRandomVector(b, 3, 1.0);
  dMakeRandomVector(c, 3, 1.0);

  dCROSS(a1, =, b, c);

  dSetZero(B, 12);
  dCROSSMAT(B, b, 4, +, -);
  dMultiply0(a2, B, c, 3, 3, 1);

  dReal diff = dMaxDifference(a1, a2, 3, 1);
  ASSERT_TRUE(diff < tol);
}

TEST(testSetZero, ) {
  HEADER;
  dReal a[100];
  dMakeRandomVector(a, 100, 1.0);
  dSetZero(a, 100);
  for (int i = 0; i < 100; i++) ASSERT_TRUE(a[i] == 0.0);
}

TEST(testNormalize, 3) {
  HEADER;
  int i, j, bad = 0;
  dVector3 n1, n2;
  for (i = 0; i < 1000; i++) {
    dMakeRandomVector(n1, 3, 1.0);
    for (j = 0; j < 3; j++) n2[j] = n1[j];
    dNormalize3(n2);
    if (dFabs(dDOT(n2, n2) - 1.0) > tol) bad |= 1;
    if (dFabs(n2[0] / n1[0] - n2[1] / n1[1]) > tol) bad |= 2;
    if (dFabs(n2[0] / n1[0] - n2[2] / n1[2]) > tol) bad |= 4;
    if (dFabs(n2[1] / n1[1] - n2[2] / n1[2]) > tol) bad |= 8;
    if (dFabs(dDOT(n2, n1) - dSqrt(dDOT(n1, n1))) > tol) bad |= 16;

    ASSERT_FALSE(bad);
  }
}

TEST(testPlaneSpace, ) {
  HEADER;
  dVector3 n, p, q;
  int bad = 0;
  for (int i = 0; i < 1000; i++) {
    dMakeRandomVector(n, 3, 1.0);
    dNormalize3(n);
    dPlaneSpace(n, p, q);
    if (fabs(dDOT(n, p)) > tol) bad = 1;
    if (fabs(dDOT(n, q)) > tol) bad = 1;
    if (fabs(dDOT(p, q)) > tol) bad = 1;
    if (fabs(dDOT(p, p) - 1) > tol) bad = 1;
    if (fabs(dDOT(q, q) - 1) > tol) bad = 1;
  }
  ASSERT_FALSE(bad);
}

//****************************************************************************
// test matrix functions

#define MSIZE 21
#define MSIZE4 24  // MSIZE rounded up to 4

void fat_matrix(int dim) {
  dReal *A = (dReal *)malloc(sizeof(dReal) * dim * dim);
  dReal *B = (dReal *)malloc(sizeof(dReal) * dim * dim);
  dReal *C = (dReal *)malloc(sizeof(dReal) * dim * dim);
  ASSERT_TRUE(A && B && C);
  dSetValue(A, 1, dim * dim);
  dSetValue(B, 2, dim * dim);
  dSetZero(C, dim * dim);
  dMultiply0(C, A, B, dim, dim, dim);
  dSetValue(C, 1, dim * dim);
  dMultiply1(C, A, B, dim, dim, dim);
  dSetValue(C, 2, dim * dim);
  dMultiply2(C, A, B, dim, dim, dim);
  dSetValue(C, 3, dim * dim);
  free(A);
  free(B);
  free(C);
}

TEST(testMatrixMultiply, ) {
  // A is 2x3, B is 3x4, B2 is B except stored columnwise, C is 2x4
  dReal A[8], B[12], A2[12], B2[16], C[8];
  int i;

  HEADER;
  dSetZero(A, 8);
  for (i = 0; i < 3; i++) A[i] = i + 2;
  for (i = 0; i < 3; i++) A[i + 4] = i + 3 + 2;
  for (i = 0; i < 12; i++) B[i] = i + 8;
  dSetZero(A2, 12);
  for (i = 0; i < 6; i++) A2[i + 2 * (i / 2)] = A[i + i / 3];
  dSetZero(B2, 16);
  for (i = 0; i < 12; i++) B2[i + i / 3] = B[i];

  dMultiply0(C, A, B, 2, 3, 4);
  ASSERT_FALSE(C[0] != 116 || C[1] != 125 || C[2] != 134 || C[3] != 143 ||
               C[4] != 224 || C[5] != 242 || C[6] != 260 || C[7] != 278);

  dMultiply1(C, A2, B, 2, 3, 4);
  ASSERT_FALSE(C[0] != 160 || C[1] != 172 || C[2] != 184 || C[3] != 196 ||
               C[4] != 196 || C[5] != 211 || C[6] != 226 || C[7] != 241);

  dMultiply2(C, A, B2, 2, 3, 4);
  ASSERT_FALSE(C[0] != 83 || C[1] != 110 || C[2] != 137 || C[3] != 164 ||
               C[4] != 164 || C[5] != 218 || C[6] != 272 || C[7] != 326);
}

TEST(testSmallMatrixMultiply, ) {
  dMatrix3 A, B, C, A2;
  dVector3 a, a2, x;

  HEADER;
  dMakeRandomMatrix(A, 3, 3, 1.0);
  dMakeRandomMatrix(B, 3, 3, 1.0);
  dMakeRandomMatrix(C, 3, 3, 1.0);
  dMakeRandomMatrix(x, 3, 1, 1.0);

  // dMULTIPLY0_331()
  dMULTIPLY0_331(a, B, x);
  dMultiply0(a2, B, x, 3, 3, 1);
  ASSERT_FALSE(dMaxDifference(a, a2, 3, 1) > tol);

  // dMULTIPLY1_331()
  dMULTIPLY1_331(a, B, x);
  dMultiply1(a2, B, x, 3, 3, 1);
  ASSERT_FALSE(dMaxDifference(a, a2, 3, 1) > tol);

  // dMULTIPLY0_133
  dMULTIPLY0_133(a, x, B);
  dMultiply0(a2, x, B, 1, 3, 3);
  ASSERT_FALSE(dMaxDifference(a, a2, 1, 3) > tol);

  // dMULTIPLY0_333()
  dMULTIPLY0_333(A, B, C);
  dMultiply0(A2, B, C, 3, 3, 3);
  ASSERT_FALSE(dMaxDifference(A, A2, 3, 3) > tol);

  // dMULTIPLY1_333()
  dMULTIPLY1_333(A, B, C);
  dMultiply1(A2, B, C, 3, 3, 3);
  ASSERT_FALSE(dMaxDifference(A, A2, 3, 3) > tol);

  // dMULTIPLY2_333()
  dMULTIPLY2_333(A, B, C);
  dMultiply2(A2, B, C, 3, 3, 3);
  ASSERT_FALSE(dMaxDifference(A, A2, 3, 3) > tol);
}

TEST(testCholeskyFactorization, ) {
  dReal A[MSIZE4 * MSIZE], B[MSIZE4 * MSIZE], C[MSIZE4 * MSIZE], diff;
  HEADER;
  dMakeRandomMatrix(A, MSIZE, MSIZE, 1.0);
  dMultiply2(B, A, A, MSIZE, MSIZE, MSIZE);
  memcpy(A, B, MSIZE4 * MSIZE * sizeof(dReal));
  ASSERT_TRUE(dFactorCholesky(B, MSIZE));
  dClearUpperTriangle(B, MSIZE);
  dMultiply2(C, B, B, MSIZE, MSIZE, MSIZE);
  diff = dMaxDifference(A, C, MSIZE, MSIZE);
  ASSERT_FALSE(diff > tol);
}

TEST(testCholeskySolve, ) {
  dReal A[MSIZE4 * MSIZE], L[MSIZE4 * MSIZE], b[MSIZE], x[MSIZE], btest[MSIZE],
      diff;
  HEADER;

  // get A,L = PD matrix
  dMakeRandomMatrix(A, MSIZE, MSIZE, 1.0);
  dMultiply2(L, A, A, MSIZE, MSIZE, MSIZE);
  memcpy(A, L, MSIZE4 * MSIZE * sizeof(dReal));

  // get b,x = right hand side
  dMakeRandomMatrix(b, MSIZE, 1, 1.0);
  memcpy(x, b, MSIZE * sizeof(dReal));

  // factor L
  ASSERT_TRUE(dFactorCholesky(L, MSIZE));
  dClearUpperTriangle(L, MSIZE);

  // solve A*x = b
  dSolveCholesky(L, x, MSIZE);

  // compute A*x and compare it with b
  dMultiply2(btest, A, x, MSIZE, MSIZE, 1);
  diff = dMaxDifference(b, btest, MSIZE, 1);
  printf("\tmaximum difference = %.6e \n", diff);
  ASSERT_FALSE(diff > tol * 5);
}

TEST(testInvertPDMatrix, ) {
  int i, j, ok;
  dReal A[MSIZE4 * MSIZE], Ainv[MSIZE4 * MSIZE], I[MSIZE4 * MSIZE];
  HEADER;

  dAASSERT(cmp(std::numeric_limits<dReal>::min(), 0));
  // dAASSERT(cmp(0., 0.));
  // dAASSERT(cmp(1e-7, 0));
  // dAASSERT(cmp(1e-8, 0));

#define showMat(MAT) printMatrixf(#MAT, MAT, MSIZE, MSIZE4)

  dMakeRandomMatrix(A, MSIZE, MSIZE, 1.0);
  //   showMat(A);
  dMultiply2(Ainv, A, A, MSIZE, MSIZE, MSIZE);
  //   showMat(Ainv);
  memcpy(A, Ainv, MSIZE4 * MSIZE * sizeof(dReal));
  dSetZero(Ainv, MSIZE4 * MSIZE);

  ASSERT_TRUE(dInvertPDMatrix(A, Ainv, MSIZE));
  dMultiply0(I, A, Ainv, MSIZE, MSIZE, MSIZE);

  //   showMat(I);
  // compare with identity
  auto err = 1e-4;
  for (i = 0; i < MSIZE; i++) {
    for (j = 0; j < MSIZE; j++) {
      if (i != j) {
        auto diff = std::abs(I[i * MSIZE4 + j] - 0.);
        if (diff > err) printf("diff OF : %f\n", diff);
        ASSERT_TRUE(IsAlmostEqual(I[i * MSIZE4 + j], 0, 2, err));
      }
    }
  }
  for (i = 0; i < MSIZE; i++) {
    auto diff = std::abs(I[i * MSIZE4 + i] - 1.);
    if (diff > err) printf("diff OF : %f\n", diff);
    ASSERT_TRUE(diff < err);
  }
}

TEST(testIsPositiveDefinite, ) {
  dReal A[MSIZE4 * MSIZE], B[MSIZE4 * MSIZE];
  HEADER;
  dMakeRandomMatrix(A, MSIZE, MSIZE, 1.0);
  dMultiply2(B, A, A, MSIZE, MSIZE, MSIZE);
  ASSERT_FALSE(dIsPositiveDefinite(A, MSIZE));
  ASSERT_TRUE(dIsPositiveDefinite(B, MSIZE));
}

TEST(testFastLDLTFactorization, ) {
  int i, j;
  dReal A[MSIZE4 * MSIZE], L[MSIZE4 * MSIZE], DL[MSIZE4 * MSIZE],
      ATEST[MSIZE4 * MSIZE], d[MSIZE], diff;
  HEADER;
  dMakeRandomMatrix(A, MSIZE, MSIZE, 1.0);
  dMultiply2(L, A, A, MSIZE, MSIZE, MSIZE);
  memcpy(A, L, MSIZE4 * MSIZE * sizeof(dReal));

  dFactorLDLT(L, d, MSIZE, MSIZE4);
  dClearUpperTriangle(L, MSIZE);
  for (i = 0; i < MSIZE; i++) L[i * MSIZE4 + i] = 1.0;

  dSetZero(DL, MSIZE4 * MSIZE);
  for (i = 0; i < MSIZE; i++) {
    for (j = 0; j < MSIZE; j++) DL[i * MSIZE4 + j] = L[i * MSIZE4 + j] / d[j];
  }

  dMultiply2(ATEST, L, DL, MSIZE, MSIZE, MSIZE);
  diff = dMaxDifference(A, ATEST, MSIZE, MSIZE);
  ASSERT_TRUE(diff < tol);
  // printf ("\tmaximum difference = %.6e - %s\n",diff,
  //   diff > tol ? "FAILED" : "passed");
}

TEST(testSolveLDLT, ) {
  dReal A[MSIZE4 * MSIZE], L[MSIZE4 * MSIZE], d[MSIZE], x[MSIZE], b[MSIZE],
      btest[MSIZE], diff;
  HEADER;
  dMakeRandomMatrix(A, MSIZE, MSIZE, 1.0);
  dMultiply2(L, A, A, MSIZE, MSIZE, MSIZE);
  memcpy(A, L, MSIZE4 * MSIZE * sizeof(dReal));

  dMakeRandomMatrix(b, MSIZE, 1, 1.0);
  memcpy(x, b, MSIZE * sizeof(dReal));

  dFactorLDLT(L, d, MSIZE, MSIZE4);
  dSolveLDLT(L, d, x, MSIZE, MSIZE4);

  dMultiply2(btest, A, x, MSIZE, MSIZE, 1);
  diff = dMaxDifference(b, btest, MSIZE, 1);
  if (diff > tol) printf("\tmaximum difference = %.6e\n", diff);
  ASSERT_TRUE(diff < tol * 30);
}

TEST(testLDLTAddTL, ) {
  int i, j;
  dReal A[MSIZE4 * MSIZE], L[MSIZE4 * MSIZE], d[MSIZE], a[MSIZE],
      DL[MSIZE4 * MSIZE], ATEST[MSIZE4 * MSIZE], diff;
  HEADER;

  dMakeRandomMatrix(A, MSIZE, MSIZE, 1.0);
  dMultiply2(L, A, A, MSIZE, MSIZE, MSIZE);
  memcpy(A, L, MSIZE4 * MSIZE * sizeof(dReal));
  dFactorLDLT(L, d, MSIZE, MSIZE4);

  // delete first row and column of factorization
  for (i = 0; i < MSIZE; i++) a[i] = -A[i * MSIZE4];
  a[0] += 1;
  dLDLTAddTL(L, d, a, MSIZE, MSIZE4);
  for (i = 1; i < MSIZE; i++) L[i * MSIZE4] = 0;
  d[0] = 1;

  // get modified L*D*L'
  dClearUpperTriangle(L, MSIZE);
  for (i = 0; i < MSIZE; i++) L[i * MSIZE4 + i] = 1.0;
  dSetZero(DL, MSIZE4 * MSIZE);
  for (i = 0; i < MSIZE; i++) {
    for (j = 0; j < MSIZE; j++) DL[i * MSIZE4 + j] = L[i * MSIZE4 + j] / d[j];
  }
  dMultiply2(ATEST, L, DL, MSIZE, MSIZE, MSIZE);

  // compare it to A with its first row/column removed
  for (i = 1; i < MSIZE; i++) A[i * MSIZE4] = A[i] = 0;
  A[0] = 1;
  diff = dMaxDifference(A, ATEST, MSIZE, MSIZE);
  if (diff > tol) printf("\tmaximum difference = %.6e\n", diff);
  ASSERT_TRUE(diff < tol);
  // printf ("\tmaximum difference = %.6e - %s\n",diff,
  //   diff > tol ? "FAILED" : "passed");
}

TEST(testLDLTRemove, ) {
  int i, j, r, p[MSIZE];
  dReal A[MSIZE4 * MSIZE], L[MSIZE4 * MSIZE], d[MSIZE], L2[MSIZE4 * MSIZE],
      d2[MSIZE], DL2[MSIZE4 * MSIZE], Atest1[MSIZE4 * MSIZE],
      Atest2[MSIZE4 * MSIZE], diff, maxdiff;
  HEADER;

  // make array of A row pointers
  dReal *Arows[MSIZE];
  for (i = 0; i < MSIZE; i++) Arows[i] = A + i * MSIZE4;

  // fill permutation vector
  for (i = 0; i < MSIZE; i++) p[i] = i;

  dMakeRandomMatrix(A, MSIZE, MSIZE, 1.0);
  dMultiply2(L, A, A, MSIZE, MSIZE, MSIZE);
  memcpy(A, L, MSIZE4 * MSIZE * sizeof(dReal));
  dFactorLDLT(L, d, MSIZE, MSIZE4);

  maxdiff = 1e10;
  for (r = 0; r < MSIZE; r++) {
    // get Atest1 = A with row/column r removed
    memcpy(Atest1, A, MSIZE4 * MSIZE * sizeof(dReal));
    dRemoveRowCol(Atest1, MSIZE, MSIZE4, r);

    // test that the row/column removal worked
    int bad = 0;
    for (i = 0; i < MSIZE; i++) {
      for (j = 0; j < MSIZE; j++) {
        if (i != r && j != r) {
          int ii = i;
          int jj = j;
          if (ii >= r) ii--;
          if (jj >= r) jj--;
          if (A[i * MSIZE4 + j] != Atest1[ii * MSIZE4 + jj]) bad = 1;
        }
      }
    }
    if (bad) printf("\trow/col removal FAILED for row %d\n", r);

    // zero out last row/column of Atest1
    for (i = 0; i < MSIZE; i++) {
      Atest1[(MSIZE - 1) * MSIZE4 + i] = 0;
      Atest1[i * MSIZE4 + MSIZE - 1] = 0;
    }

    // get L2*D2*L2' = adjusted factorization to remove that row
    memcpy(L2, L, MSIZE4 * MSIZE * sizeof(dReal));
    memcpy(d2, d, MSIZE * sizeof(dReal));
    dLDLTRemove(/*A*/ Arows, p, L2, d2, MSIZE, MSIZE, r, MSIZE4);

    // get Atest2 = L2*D2*L2'
    dClearUpperTriangle(L2, MSIZE);
    for (i = 0; i < (MSIZE - 1); i++) L2[i * MSIZE4 + i] = 1.0;
    for (i = 0; i < MSIZE; i++) L2[(MSIZE - 1) * MSIZE4 + i] = 0;
    d2[MSIZE - 1] = 1;
    dSetZero(DL2, MSIZE4 * MSIZE);
    for (i = 0; i < (MSIZE - 1); i++) {
      for (j = 0; j < MSIZE - 1; j++)
        DL2[i * MSIZE4 + j] = L2[i * MSIZE4 + j] / d2[j];
    }

    dMultiply2(Atest2, L2, DL2, MSIZE, MSIZE, MSIZE);

    diff = dMaxDifference(Atest1, Atest2, MSIZE, MSIZE);
    if (diff < maxdiff) maxdiff = diff;

    /*
    dPrintMatrix (Atest1,MSIZE,MSIZE);
    printf ("\n");
    dPrintMatrix (Atest2,MSIZE,MSIZE);
    printf ("\n");
    */
  }
  diff = maxdiff;
  if (diff > tol) printf("\tmaximum difference = %.6e\n", diff);
  ASSERT_TRUE(diff < tol);
}

//****************************************************************************
// test mass stuff

#define NUMP 10  // number of particles

void printMassParams(dMass *m) {
  printf("mass = %.4f\n", m->mass);
  printf("com  = (%.4f,%.4f,%.4f)\n", m->c[0], m->c[1], m->c[2]);
  printf(
      "I    = [ %10.4f %10.4f %10.4f ]\n"
      "       [ %10.4f %10.4f %10.4f ]\n"
      "       [ %10.4f %10.4f %10.4f ]\n",
      m->_I(0, 0), m->_I(0, 1), m->_I(0, 2), m->_I(1, 0), m->_I(1, 1),
      m->_I(1, 2), m->_I(2, 0), m->_I(2, 1), m->_I(2, 2));
}

void compareMassParams(dMass *m1, dMass *m2, const char *msg) {
  int i, j, ok = 1;
  ASSERT_TRUE(cmp(m1->mass, m2->mass) && cmp(m1->c[0], m2->c[0]) &&
              cmp(m1->c[1], m2->c[1]) && cmp(m1->c[2], m2->c[2]));
  //       {
  //   printf("m1, m2 ->c mass not equal!\n");
  //   ok = 0;
  // }

  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      if (!cmp(m1->_I(i, j), m2->_I(i, j))) {
        printf("m1, m2 ->I[%d, %d] {%f, %f} mass not equal!\n", i, j,
               m1->_I(i, j), m2->_I(i, j));
        printMassParams(m1);
        printMassParams(m2);
        ok = 0;
      }
  ASSERT_TRUE(ok);
}

// compute the mass parameters of a particle set
void computeMassParams(dMass *m, dReal q[NUMP][3], dReal pm[NUMP]) {
  int i, j;
  dMassSetZero(m);
  for (i = 0; i < NUMP; i++) {
    m->mass += pm[i];
    for (j = 0; j < 3; j++) m->c[j] += pm[i] * q[i][j];
    m->_I(0, 0) += pm[i] * (q[i][1] * q[i][1] + q[i][2] * q[i][2]);
    m->_I(1, 1) += pm[i] * (q[i][0] * q[i][0] + q[i][2] * q[i][2]);
    m->_I(2, 2) += pm[i] * (q[i][0] * q[i][0] + q[i][1] * q[i][1]);
    m->_I(0, 1) -= pm[i] * (q[i][0] * q[i][1]);
    m->_I(0, 2) -= pm[i] * (q[i][0] * q[i][2]);
    m->_I(1, 2) -= pm[i] * (q[i][1] * q[i][2]);
  }
  for (j = 0; j < 3; j++) m->c[j] /= m->mass;
  m->_I(1, 0) = m->_I(0, 1);
  m->_I(2, 0) = m->_I(0, 2);
  m->_I(2, 1) = m->_I(1, 2);
}

TEST(testMassFunctions, ) {
  auto ShowMass = [](dMass *m) {};
  dMass m;
  int i, j;
  dReal q[NUMP][3];  // particle positions
  dReal pm[NUMP];    // particle masses
  dMass m1{}, m2{};
  dMatrix3 R;

  HEADER;

  printf("\t");
  dMassSetZero(&m);
  TRAP_MESSAGE(dMassSetParameters(&m, 10, 0, 0, 0, 1, 2, 3, 4, 5, 6),
               printf(" FAILED (1)\n"), printf(" passed (1)\n"));

  printf("\t");
  dMassSetZero(&m);
  TRAP_MESSAGE(
      dMassSetParameters(&m, 10, 0.1, 0.2, 0.15, 3, 5, 14, 3.1, 3.2, 4),
      printf("passed (2)\n"), printf(" FAILED (2)\n"));
  EXPECT_TRUE(m.mass == 10 && m.c[0] == REAL(0.1) && m.c[1] == REAL(0.2));
  EXPECT_TRUE(m.c[2] == REAL(0.15) && m._I(0, 0) == 3 && m._I(1, 1) == 5);
  EXPECT_TRUE(m._I(2, 2) == 14 && m._I(0, 1) == REAL(3.1));
  EXPECT_TRUE(m._I(0, 2) == REAL(3.2) && m._I(1, 2) == 4);
  EXPECT_TRUE(m._I(1, 0) == REAL(3.1) && m._I(2, 0) == REAL(3.2));
  ASSERT_DOUBLE_EQ(m._I(2, 1), 4);

  dMassSetZero(&m);
  dMassSetSphere(&m, 1.4, 0.86);
  EXPECT_TRUE(cmp(m.mass, 3.73002719949386) && m.c[0] == 0 && m.c[1] == 0 &&
              m.c[2] == 0 && cmp(m._I(0, 0), 1.10349124669826) &&
              cmp(m._I(1, 1), 1.10349124669826) &&
              cmp(m._I(2, 2), 1.10349124669826) && m._I(0, 1) == 0 &&
              m._I(0, 2) == 0 && m._I(1, 2) == 0 && m._I(1, 0) == 0 &&
              m._I(2, 0) == 0 && m._I(2, 1) == 0);

  dMassSetZero(&m);
  dMassSetCapsule(&m, 1.3, 1, 0.76, 1.53);
  EXPECT_TRUE(cmp(m.mass, 5.99961928996029) && m.c[0] == 0 && m.c[1] == 0 &&
              m.c[2] == 0 && cmp(m._I(0, 0), 1.59461986077384) &&
              cmp(m._I(1, 1), 4.21878433864904) &&
              cmp(m._I(2, 2), 4.21878433864904) && m._I(0, 1) == 0 &&
              m._I(0, 2) == 0 && m._I(1, 2) == 0 && m._I(1, 0) == 0 &&
              m._I(2, 0) == 0 && m._I(2, 1) == 0);

  dMassSetZero(&m);
  dMassSetBox(&m, 0.27, 3, 4, 5);
  EXPECT_TRUE(cmp(m.mass, 16.2) && m.c[0] == 0 && m.c[1] == 0 && m.c[2] == 0 &&
              cmp(m._I(0, 0), 55.35) && cmp(m._I(1, 1), 45.9) &&
              cmp(m._I(2, 2), 33.75) && m._I(0, 1) == 0 && m._I(0, 2) == 0 &&
              m._I(1, 2) == 0 && m._I(1, 0) == 0 && m._I(2, 0) == 0 &&
              m._I(2, 1) == 0);

  // test dMassAdjust?

  // make random particles and compute the mass, COM and inertia, then
  // translate and repeat.
  for (i = 0; i < NUMP; i++) {
    pm[i] = dRandReal() + 0.5;
    for (j = 0; j < 3; j++) {
      q[i][j] = 2.0 * (dRandReal() - 0.5);
    }
  }
  computeMassParams(&m1, q, pm);
  memcpy(&m2, &m1, sizeof(dMass));
  dMassTranslate(&m2, 1, 2, -3);
  for (i = 0; i < NUMP; i++) {
    q[i][0] += 1;
    q[i][1] += 2;
    q[i][2] -= 3;
  }
  computeMassParams(&m1, q, pm);
  compareMassParams(&m1, &m2, "7");

  // rotate the masses
  _R(0, 0) = -0.87919618797635;
  _R(0, 1) = 0.15278881840384;
  _R(0, 2) = -0.45129772879842;
  _R(1, 0) = -0.47307856232664;
  _R(1, 1) = -0.39258064912909;
  _R(1, 2) = 0.78871864932708;
  _R(2, 0) = -0.05666336483842;
  _R(2, 1) = 0.90693771059546;
  _R(2, 2) = 0.41743652473765;
  dMassRotate(&m2, R);
  for (i = 0; i < NUMP; i++) {
    dReal a[3];
    dMultiply0(a, &_R(0, 0), &q[i][0], 3, 3, 1);
    q[i][0] = a[0];
    q[i][1] = a[1];
    q[i][2] = a[2];
  }
  computeMassParams(&m1, q, pm);
  compareMassParams(&m1, &m2, "8");
}

// //****************************************************************************
// // test rotation stuff
void makeRandomRotation(dMatrix3 R) {
  dReal *u1 = R, *u2 = R + 4, *u3 = R + 8;
  dMakeRandomVector(u1, 3, 1.0);
  dNormalize3(u1);
  dMakeRandomVector(u2, 3, 1.0);
  dReal d = dDOT(u1, u2);
  u2[0] -= d * u1[0];
  u2[1] -= d * u1[1];
  u2[2] -= d * u1[2];
  dNormalize3(u2);
  dCROSS(u3, =, u1, u2);
}

TEST(testRtoQandQtoR, ) {
  HEADER;
  dMatrix3 R, I, R2;
  dQuaternion q;
  int i;

  // test makeRandomRotation()
  makeRandomRotation(R);
  dMultiply2(I, R, R, 3, 3, 3);
  ASSERT_TRUE(cmpIdentityMat3(I));

  for (i = 0; i < 100; i++) {
    dMakeRandomVector(q, 4, 1.0);
    dNormalize4(q);
    dQtoR(q, R);
    dMultiply2(I, R, R, 3, 3, 3);
    ASSERT_TRUE(cmpIdentityMat3(I));
  }

  // test R -> Q -> R works
  dReal maxdiff = 0;
  for (i = 0; i < 100; i++) {
    makeRandomRotation(R);
    dRtoQ(R, q);
    dQtoR(q, R2);
    dReal diff = dMaxDifference(R, R2, 3, 3);
    if (diff > maxdiff) maxdiff = diff;
  }
  ASSERT_TRUE(maxdiff < tol);
}

TEST(testQuaternionMultiply, ) {
  HEADER;
  dMatrix3 RA, RB, RC, Rtest;
  dQuaternion qa, qb, qc;
  dReal diff, maxdiff = 0;

  for (int i = 0; i < 100; i++) {
    makeRandomRotation(RB);
    makeRandomRotation(RC);
    dRtoQ(RB, qb);
    dRtoQ(RC, qc);

    dMultiply0(RA, RB, RC, 3, 3, 3);
    dQMultiply0(qa, qb, qc);
    dQtoR(qa, Rtest);
    diff = dMaxDifference(Rtest, RA, 3, 3);
    if (diff > maxdiff) maxdiff = diff;

    dMultiply1(RA, RB, RC, 3, 3, 3);
    dQMultiply1(qa, qb, qc);
    dQtoR(qa, Rtest);
    diff = dMaxDifference(Rtest, RA, 3, 3);
    if (diff > maxdiff) maxdiff = diff;

    dMultiply2(RA, RB, RC, 3, 3, 3);
    dQMultiply2(qa, qb, qc);
    dQtoR(qa, Rtest);
    diff = dMaxDifference(Rtest, RA, 3, 3);
    if (diff > maxdiff) maxdiff = diff;

    dMultiply0(RA, RC, RB, 3, 3, 3);
    transpose3x3(RA);
    dQMultiply3(qa, qb, qc);
    dQtoR(qa, Rtest);
    diff = dMaxDifference(Rtest, RA, 3, 3);
    if (diff > maxdiff) maxdiff = diff;
  }
  diff = maxdiff;
  if (diff > tol) printf("\tmaximum difference = %.6e\n", diff);
  ASSERT_TRUE(diff < tol);
}

TEST(testRotationFunctions, ) {
  dMatrix3 R1;
  HEADER;

  printf("\tdRSetIdentity - ");
  dMakeRandomMatrix(R1, 3, 3, 1.0);
  dRSetIdentity(R1);
  ASSERT_TRUE(cmpIdentityMat3(R1));
  //  printf ("passed\n"); else printf ("FAILED\n");

  // printf ("\tdRFromAxisAndAngle - ");

  // printf ("\n");
  // printf ("\tdRFromEulerAngles - ");

  // printf ("\n");
  // printf ("\tdRFrom2Axes - ");

  // printf ("\n");
}

//****************************************************************************

#include "array.h"
#define dDebug(CON, MSG, ...)                   \
  if (CON) {                                    \
    EXPECT_TRUE(false) << "diagnostic message"; \
  }

// matrix header on the stack

class dMatrixComparison {
  struct dMatInfo;
  dArray<dMatInfo *> mat;
  int afterfirst, index;

 public:
  dMatrixComparison();
  ~dMatrixComparison();

  dReal nextMatrix(dReal *A, int n, int m, int lower_tri, const char *name,
                   ...);
  // add a new n*m matrix A to the sequence. the name of the matrix is given
  // by the printf-style arguments (name,...). if this is the first sequence
  // then this object will simply record the matrices and return 0.
  // if this the second or subsequent sequence then this object will compare
  // the matrices with the first sequence, and report any differences.
  // the matrix error will be returned. if `lower_tri' is 1 then only the
  // lower triangle of the matrix (including the diagonal) will be compared
  // (the matrix must be square).

  void end();
  // end a sequence.

  void reset();
  // restarts the object, so the next sequence will be the first sequence.

  void dump();
  // print out info about all the matrices in the sequence
};

struct dMatrixComparison::dMatInfo {
  int n, m;        // size of matrix
  char name[128];  // name of the matrix
  dReal *data;     // matrix data
  int size;        // size of `data'
};

dMatrixComparison::dMatrixComparison() {
  afterfirst = 0;
  index = 0;
}

dMatrixComparison::~dMatrixComparison() { reset(); }

dReal dMatrixComparison::nextMatrix(dReal *A, int n, int m, int lower_tri,
                                    const char *name, ...) {
  if (A == 0 || n < 1 || m < 1 || name == 0)
    dDebug(0, "bad args to nextMatrix");
  int num = n * dPAD(m);

  if (afterfirst == 0) {
    dMatInfo *mi = (dMatInfo *)dAlloc(sizeof(dMatInfo));
    mi->n = n;
    mi->m = m;
    mi->size = num * sizeof(dReal);
    mi->data = (dReal *)dAlloc(mi->size);
    memcpy(mi->data, A, mi->size);

    va_list ap;
    va_start(ap, name);
    vsprintf(mi->name, name, ap);
    if (strlen(mi->name) >= sizeof(mi->name)) dDebug(0, "name too long");

    mat.push(mi);
    return 0;
  } else {
    if (lower_tri && n != m)
      dDebug(0, "dMatrixComparison, lower triangular matrix must be square");
    if (index >= mat.size()) dDebug(0, "dMatrixComparison, too many matrices");
    dMatInfo *mp = mat[index];
    index++;

    dMatInfo mi;
    va_list ap;
    va_start(ap, name);
    vsprintf(mi.name, name, ap);
    if (strlen(mi.name) >= sizeof(mi.name)) dDebug(0, "name too long");

    if (strcmp(mp->name, mi.name) != 0)
      dDebug(0, "dMatrixComparison, name mismatch (\"%s\" and \"%s\")",
             mp->name, mi.name);
    if (mp->n != n || mp->m != m)
      dDebug(0, "dMatrixComparison, size mismatch (%dx%d and %dx%d)", mp->n,
             mp->m, n, m);

    dReal maxdiff;
    if (lower_tri) {
      maxdiff = dMaxDifferenceLowerTriangle(A, mp->data, n);
    } else {
      maxdiff = dMaxDifference(A, mp->data, n, m);
    }
    if (maxdiff > tol)
      dDebug(0,
             "dMatrixComparison, matrix error (size=%dx%d, name=\"%s\", "
             "error=%.4e)",
             n, m, mi.name, maxdiff);
    return maxdiff;
  }
}

void dMatrixComparison::end() {
  if (mat.size() <= 0) dDebug(0, "no matrices in sequence");
  afterfirst = 1;
  index = 0;
}

void dMatrixComparison::reset() {
  for (int i = 0; i < mat.size(); i++) {
    dFree(mat[i]->data, mat[i]->size);
    dFree(mat[i], sizeof(dMatInfo));
  }
  mat.setSize(0);
  afterfirst = 0;
  index = 0;
}

void dMatrixComparison::dump() {
  for (int i = 0; i < mat.size(); i++)
    printf("%d: %s (%dx%d)\n", i, mat[i]->name, mat[i]->n, mat[i]->m);
}

//****************************************************************************
// unit test

TEST(dTestMatrixComparison, ) {
  volatile int i;
  dMessageFunction *orig_debug = dGetDebugHandler();

  dMatrixComparison mc;
  dReal A[50 * 50];

  // make first sequence
  unsigned long seed = dRandGetSeed();
  for (i = 1; i < 49; i++) {
    dMakeRandomMatrix(A, i, i + 1, 1.0);
    mc.nextMatrix(A, i, i + 1, 0, "A%d", i);
  }
  mc.end();

  // mc.dump();

  dRandSetSeed(seed);

  for (i = 1; i < 49; i++) {
    dMakeRandomMatrix(A, i, i + 1, 1.0);
    mc.nextMatrix(A, i, i + 1, 0, "A%d", i);
  }
  mc.end();

  // test broken sequences (with matrix error)
  dRandSetSeed(seed);
  volatile int passcount = 0;
  for (i = 1; i < 49; i++) {
    dMakeRandomMatrix(A, i, i + 1, 1.0);
    A[(i - 1) * dPAD(i + 1) + i] += REAL(0.01);
    mc.nextMatrix(A, i, i + 1, 0, "A%d", i);
    passcount++;
  }
  mc.end();
  EXPECT_TRUE(passcount == 48);

  // test broken sequences (with name error)
  dRandSetSeed(seed);
  passcount = 0;
  for (i = 1; i < 49; i++) {
    dMakeRandomMatrix(A, i, i + 1, 1.0);
    mc.nextMatrix(A, i, i + 1, 0, "B%d", i);
    passcount++;
  }
  mc.end();
  EXPECT_TRUE(passcount == 48);
  dRandSetSeed(seed);

  for (i = 1; i < 49; i++) {
    dMakeRandomMatrix(A, i, i + 1, 1.0);
    mc.nextMatrix(A, i, i + 1, 0, "A%d", i);
  }
  mc.end();
}

#define DO(x)
#define NUM 100

extern void checkWorld(dxWorld *w);

TEST(dTestDataStructures, ) {
  int i;
  DO(printf("testDynamicsStuff()\n"));

  // dRandSetSeed(0);

  dBodyID body[NUM];
  int nb = 0;
  dJointID joint[NUM];
  int nj = 0;

  for (i = 0; i < NUM; i++) body[i] = 0;
  for (i = 0; i < NUM; i++) joint[i] = 0;

  DO(printf("creating world\n"));
  dWorldID w = dWorldCreate();
  checkWorld(w);

  int test_c = 1e3;
  while (test_c-- > 0) {
    if (nb < NUM && dRandReal() > 0.5) {
      DO(printf("creating body\n"));
      body[nb] = dBodyCreate(w);
      DO(printf("\t--> %p\n", body[nb]));
      nb++;
      checkWorld(w);
      DO(printf("%d BODIES, %d JOINTS\n", nb, nj));
    }
    if (nj < NUM && nb > 2 && dRandReal() > 0.5) {
      dBodyID b1 = body[dRand() % nb];
      dBodyID b2 = body[dRand() % nb];
      if (b1 != b2) {
        DO(printf("creating joint, attaching to %p,%p\n", b1, b2));
        joint[nj] = dJointCreateBall(w, nullptr);
        DO(printf("\t-->%p\n", joint[nj]));
        checkWorld(w);
        dJointAttach(joint[nj], b1, b2);
        nj++;
        checkWorld(w);
        DO(printf("%d BODIES, %d JOINTS\n", nb, nj));
      }
    }
    if (nj > 0 && nb > 2 && dRandReal() > 0.5) {
      dBodyID b1 = body[dRand() % nb];
      dBodyID b2 = body[dRand() % nb];
      if (b1 != b2) {
        int k = int(dRand() % nj);
        DO(printf("reattaching joint %p\n", joint[k]));
        dJointAttach(joint[k], b1, b2);
        checkWorld(w);
        DO(printf("%d BODIES, %d JOINTS\n", nb, nj));
      }
    }
    if (nb > 0 && dRandReal() > 0.5) {
      int k = int(dRand() % nb);
      DO(printf("destroying body %p\n", body[k]));
      dBodyDestroy(body[k]);
      checkWorld(w);
      for (; k < (NUM - 1); k++) body[k] = body[k + 1];
      nb--;
      DO(printf("%d BODIES, %d JOINTS\n", nb, nj));
    }
    if (nj > 0 && dRandReal() > 0.5) {
      int k = int(dRand() % nj);
      DO(printf("destroying joint %p\n", joint[k]));
      dJointDestroy(joint[k]);
      checkWorld(w);
      for (; k < (NUM - 1); k++) joint[k] = joint[k + 1];
      nj--;
      DO(printf("%d BODIES, %d JOINTS\n", nb, nj));
    }
  }
}
