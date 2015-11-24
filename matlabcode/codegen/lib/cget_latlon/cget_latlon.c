/*
 * File: cget_latlon.c
 *
 * MATLAB Coder version            : 2.8
 * C/C++ source code generated on  : 24-Nov-2015 01:13:59
 */

/* Include Files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "cadjlon.h"
#include "rdivide.h"

/* Function Declarations */
static double rt_atan2d_snf(double u0, double u1);

/* Function Definitions */

/*
 * Arguments    : double u0
 *                double u1
 * Return Type  : double
 */
static double rt_atan2d_snf(double u0, double u1)
{
  double y;
  int b_u0;
  int b_u1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else if (rtIsInf(u0) && rtIsInf(u1)) {
    if (u0 > 0.0) {
      b_u0 = 1;
    } else {
      b_u0 = -1;
    }

    if (u1 > 0.0) {
      b_u1 = 1;
    } else {
      b_u1 = -1;
    }

    y = atan2(b_u0, b_u1);
  } else if (u1 == 0.0) {
    if (u0 > 0.0) {
      y = RT_PI / 2.0;
    } else if (u0 < 0.0) {
      y = -(double)(RT_PI / 2.0);
    } else {
      y = 0.0;
    }
  } else {
    y = atan2(u0, u1);
  }

  return y;
}

/*
 * Arguments    : double blat
 *                double blon
 *                double ranges
 *                double azm
 *                double major_axis
 *                double esquared
 *                double *plat
 *                double *plon
 * Return Type  : void
 */
void cget_latlon(double blat, double blon, double ranges, double azm, double
                 major_axis, double esquared, double *plat, double *plon)
{
  double c1;
  double c2;
  double D;
  double P;
  double ss;
  double S;
  double phi1;
  int ellipse;
  double onef;
  double f;
  double f4;
  double al12;
  int ii_size_idx_0;
  int ii_size_idx_1;
  int is1_size_idx_0;
  int is1_size_idx_1;
  signed char iv0[1];
  int loop_ub;
  int i0;
  signed char iv1[1];
  double th1;
  double costh1;
  double sinth1;
  double sina12;
  int im1_size_idx_0;
  int im1_size_idx_1;
  signed char iv2[1];
  signed char iv3[1];
  double cosa12;
  double M;
  double N;
  double dv0[1];
  double tmp_data[1];
  double dv1[1];
  int tmp_size[2];
  double b_tmp_data[1];
  int b_tmp_size[2];
  double c_tmp_data[1];
  double dv2[1];
  double dv3[1];
  double dv4[1];
  double dv5[1];
  int c_tmp_size[2];
  double dv6[1];
  double s1;
  double d0;
  double d;
  double dv7[1];
  double u;
  double V;
  double sind;
  double ds;
  double cosds;
  double sinds;
  double dv8[1];
  double al21;
  double phi2;
  double de;
  double b_ellipse;
  double lam2;

  /* CGET_LATLON        Get latitudes and longitudes along a line */
  /*  */
  /*        function [plat,plon]=cget_latlon(blat,blon,ranges,azm,major_axis); */
  /*  */
  /*  This is a translation to matlab of the USGS PROJ-4.4 geographic */
  /*  projection library to compute positions along a series of ranges */
  /*  along a given azimuth from a reference point. */
  /*  Uses ellipsoidal earth model set to Clark 1966 standard. */
  /*  */
  /*  Inputs: */
  /*        blat,blon       lat,lon coordinates of start point */
  /*                        in degrees (scalar) */
  /*  */
  /*        ranges          range in nmi (scalar or vector) */
  /*        azm             single-valued (same size as range) (degrees) */
  /*        major_axes - skip for an elliptical earth, set=0 for spherical Earth */
  /*  */
  /*  Outputs: */
  /*        plat,plon       lat,lons along path */
  /*  */
  c1 = 0.0;
  c2 = 0.0;
  D = 0.0;
  P = 0.0;
  ss = 0.0;
  if (major_axis == 0.0) {
    major_axis = 6.3782064E+6;
    esquared = 0.0;
  }

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /*  now do the stuff that used to be in .c */
  S = ranges * 1852.0;
  phi1 = blat * 0.017453292519943295;
  ellipse = (esquared != 0.0);
  if (ellipse == 1) {
    onef = sqrt(1.0 - esquared);
    f = 1.0 - onef;
    f4 = (1.0 - onef) / 4.0;
  } else {
    onef = 1.0;
    f = 0.0;
    f4 = 0.0;
  }

  al12 = azm * 0.017453292519943295;
  cadjlon(&al12);
  if (fabs(al12) > 1.5707963267948966) {
    ii_size_idx_0 = 1;
    ii_size_idx_1 = 1;
  } else {
    ii_size_idx_0 = 0;
    ii_size_idx_1 = 0;
  }

  is1_size_idx_0 = ii_size_idx_0;
  is1_size_idx_1 = ii_size_idx_1;
  iv0[0] = 0;
  loop_ub = ii_size_idx_0 * ii_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    iv0[0] = 1;
    i0 = 1;
  }

  if (fabs(al12) <= 1.5707963267948966) {
    ii_size_idx_0 = 1;
    ii_size_idx_1 = 1;
  } else {
    ii_size_idx_0 = 0;
    ii_size_idx_1 = 0;
  }

  iv1[0] = iv0[0];
  loop_ub = ii_size_idx_0 * ii_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    iv1[0] = 0;
    i0 = 1;
  }

  if (ellipse == 1) {
    th1 = atan(onef * tan(phi1));
  } else {
    th1 = phi1;
  }

  costh1 = cos(th1);
  sinth1 = sin(th1);
  sina12 = sin(al12);
  if (fabs(sina12) < 1.0E-10) {
    ii_size_idx_0 = 1;
    ii_size_idx_1 = 1;
  } else {
    ii_size_idx_0 = 0;
    ii_size_idx_1 = 0;
  }

  im1_size_idx_0 = ii_size_idx_0;
  im1_size_idx_1 = ii_size_idx_1;
  iv2[0] = 0;
  loop_ub = ii_size_idx_0 * ii_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    iv2[0] = 1;
    i0 = 1;
  }

  if (fabs(sina12) >= 1.0E-10) {
    ii_size_idx_0 = 1;
    ii_size_idx_1 = 1;
  } else {
    ii_size_idx_0 = 0;
    ii_size_idx_1 = 0;
  }

  iv3[0] = iv2[0];
  loop_ub = ii_size_idx_0 * ii_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    iv3[0] = 0;
    i0 = 1;
  }

  if (iv3[0] == 1) {
    sina12 = 0.0;
    if (fabs(al12) < 1.5707963267948966) {
      cosa12 = 1.0;
    } else {
      cosa12 = -1.0;
    }

    M = 0.0;
  } else {
    cosa12 = cos(al12);
    M = costh1 * sina12;
  }

  N = costh1 * cosa12;
  if (ellipse == 1) {
    /*    if merid(j) == 1 */
    dv0[0] = 0.0;
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv0[0] = f4;
      i0 = 1;
    }

    tmp_data[0] = 0.0;
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 1.0 - dv0[0];
      i0 = 1;
    }

    dv1[0] = tmp_data[0];
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv1[0] = tmp_data[0] * tmp_data[0];
      i0 = 1;
    }

    tmp_size[0] = im1_size_idx_0;
    tmp_size[1] = im1_size_idx_1;
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = dv0[0];
      i0 = 1;
    }

    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      b_tmp_data[0] = dv1[0];
      i0 = 1;
    }

    b_rdivide(tmp_data, tmp_size, b_tmp_data, c_tmp_data, b_tmp_size);
    dv2[0] = 0.0;
    loop_ub = b_tmp_size[0] * b_tmp_size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      dv2[0] = c_tmp_data[i0];
    }

    /*    else  */
    dv3[0] = 0.0;
    loop_ub = ii_size_idx_0 * ii_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv3[0] = f * M;
      i0 = 1;
    }

    c1 = dv3[0];
    dv4[0] = dv0[0];
    loop_ub = ii_size_idx_0 * ii_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv4[0] = f4 * (1.0 - M * M);
      i0 = 1;
    }

    c2 = dv4[0];
    dv5[0] = dv1[0];
    loop_ub = ii_size_idx_0 * ii_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv5[0] = (1.0 - dv4[0]) * ((1.0 - dv4[0]) - dv3[0] * M);
      i0 = 1;
    }

    D = dv5[0];
    c_tmp_size[0] = ii_size_idx_0;
    c_tmp_size[1] = ii_size_idx_1;
    loop_ub = ii_size_idx_0 * ii_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = (1.0 + 0.5 * dv3[0] * M) * dv4[0];
      i0 = 1;
    }

    loop_ub = ii_size_idx_0 * ii_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      b_tmp_data[0] = dv5[0];
      i0 = 1;
    }

    b_rdivide(tmp_data, c_tmp_size, b_tmp_data, c_tmp_data, b_tmp_size);
    dv6[0] = dv2[0];
    loop_ub = b_tmp_size[0] * b_tmp_size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      dv6[0] = c_tmp_data[i0];
    }

    P = dv6[0];

    /*    end */
  }

  if (iv3[0] == 1) {
    s1 = 1.5707963267948966 - th1;
  } else {
    if (fabs(M) >= 1.0) {
      d0 = 0.0;
    } else {
      d0 = acos(M);
    }

    s1 = sinth1 / sin(d0);
    if (fabs(s1) >= 1.0) {
      s1 = 0.0;
    } else {
      s1 = acos(s1);
    }
  }

  if (ellipse == 1) {
    d = rdivide(S, D * major_axis);
    dv7[0] = d;
    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv7[0] = -d;
      i0 = 1;
    }

    u = 2.0 * (s1 - dv7[0]);
    V = cos(u + dv7[0]);
    sind = sin(dv7[0]);
    ds = (dv7[0] + c2 * c2 * sind * cos(dv7[0]) * (2.0 * V * V - 1.0)) - 2.0 * P
      * V * (1.0 - 2.0 * P * cos(u)) * sind;
    ss = (s1 + s1) - ds;
  } else {
    ds = S / major_axis;
    dv7[0] = ds;
    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv7[0] = -ds;
      i0 = 1;
    }

    ds = dv7[0];
  }

  cosds = cos(ds);
  sinds = sin(ds);
  dv8[0] = sinds;
  loop_ub = is1_size_idx_0 * is1_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    dv8[0] = -sinds;
    i0 = 1;
  }

  al21 = N * cosds - sinth1 * dv8[0];
  if (iv3[0] == 1) {
    phi2 = atan(tan((1.5707963267948966 + s1) - ds) / onef);
    if (al21 > 0.0) {
      if (iv1[0] == 1) {
        de = 3.1415926535897931;
      } else {
        phi2 = -phi2;
        de = 0.0;
      }
    } else if (iv1[0] == 1) {
      phi2 = -phi2;
      de = 0.0;
    } else {
      de = 3.1415926535897931;
    }
  } else {
    al21 = atan(M / al21);
    if (al21 > 0.0) {
      al21 += 3.1415926535897931;
    }

    if (al12 < 0.0) {
      al21 -= 3.1415926535897931;
    }

    cadjlon(&al21);
    if (ellipse == 1) {
      b_ellipse = onef * M;
    } else {
      b_ellipse = M;
    }

    phi2 = atan(-(sinth1 * cosds + N * dv8[0]) * sin(al21) / b_ellipse);
    de = rt_atan2d_snf(dv8[0] * sina12, costh1 * cosds - sinth1 * dv8[0] *
                       cosa12);
    if (ellipse == 1) {
      if (iv1[0] == 1) {
        de += c1 * ((1.0 - c2) * ds + c2 * dv8[0] * cos(ss));
      } else {
        de -= c1 * ((1.0 - c2) * ds - c2 * dv8[0] * cos(ss));
      }
    }
  }

  lam2 = blon * 0.017453292519943295 + de;
  cadjlon(&lam2);
  *plon = lam2 * 57.295779513082323;
  *plat = phi2 * 57.295779513082323;
}

/*
 * File trailer for cget_latlon.c
 *
 * [EOF]
 */
