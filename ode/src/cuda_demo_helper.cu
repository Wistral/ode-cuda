#include <stdio.h>
#include <ode/cuda_demo_helper.h>

ODE_API void printMatrixPad4(const char *name, const char *fmt, dReal const *a,
                             const int h, const int w) {
  printf("%s:\n", name);
  for (int row = 0; row < h; row++) {
    for (int col = 0; col < w; col++) printf(fmt, a[row * dPAD(w) + col]);
    printf("\n");
  }
  printf("\n");
}

ODE_API void printMatrixBase(const char *name, const char *fmt, dReal const *a,
                             const int h, const int w, bool pad) {
  printf("%s:\n", name);
  for (int row = 0; row < h; row++) {
    for (int col = 0; col < w; col++)
      printf(fmt, a[(pad ? (row * dPAD(w) + col) : (row * w + col))]);
    printf("\n");
  }
  printf("\n");
}

ODE_API void printMatrix(char *name, dReal *a, int h, int w) {
  printMatrixBase(name, "%d, ", a, h, w, 0);
}

ODE_API void printMatrixf(const char *name, dReal const *a, const int h,
                          const int w) {
  printMatrixBase(name, "%f, ", a, h, w, 0);
}

ODE_API void printPaddedMatrixf(const char *name, dReal const *a, const int h,
                                const int w) {
  printMatrixBase(name, "%f, ", a, h, w, 1);
}

ODE_API void printMatrix0f(const char *name, dReal const *a, const int h,
                           const int w) {
  printMatrixBase(name, "%.0f, ", a, h, w, 0);
}

ODE_API void makeIdMatrix(dReal *a, int s, int n)
{
	for (int row=0; row<s; row++) {
		for (int col=0; col<s; col++) {
			if (row==col) { a[row*s+col] = n; }
			else { a[row*s+col] = 0; }
		}
	}
}

