/*
 * File: rdivide.c
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 24-Nov-2015 01:13:59
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "rdivide.h"

/* Function Definitions */

/*
 * Arguments    : const double x_data[]
 *                const int x_size[2]
 *                const double y_data[]
 *                double z_data[]
 *                int z_size[2]
 * Return Type  : void
 */
void b_rdivide(const double x_data[], const int x_size[2], const double y_data[],
               double z_data[], int z_size[2])
{
  int loop_ub;
  int i1;
  z_size[0] = x_size[0];
  z_size[1] = x_size[1];
  loop_ub = x_size[0] * x_size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    z_data[i1] = x_data[i1] / y_data[i1];
  }
}

/*
 * Arguments    : double x
 *                double y
 * Return Type  : double
 */
double rdivide(double x, double y)
{
  return x / y;
}

/*
 * File trailer for rdivide.c
 *
 * [EOF]
 */
