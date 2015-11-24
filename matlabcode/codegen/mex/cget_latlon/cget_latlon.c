/*
 * cget_latlon.c
 *
 * Code generation for function 'cget_latlon'
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "cget_latlon.h"
#include "eml_error.h"
#include "cadjlon.h"
#include "rdivide.h"
#include "acos.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 80, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtRSInfo b_emlrtRSI = { 87, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtRSInfo c_emlrtRSI = { 117, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtRSInfo d_emlrtRSI = { 121, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtRSInfo e_emlrtRSI = { 131, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtRSInfo f_emlrtRSI = { 137, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtRSInfo g_emlrtRSI = { 180, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtRSInfo h_emlrtRSI = { 198, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtRSInfo i_emlrtRSI = { 14, "sqrt",
  "C:\\Program Files\\MATLAB\\MATLAB Production Server\\R2015a\\toolbox\\eml\\lib\\matlab\\elfun\\sqrt.m"
};

static emlrtMCInfo emlrtMCI = { 52, 1, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo emlrtECI = { -1, 156, 1, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo b_emlrtECI = { -1, 153, 4, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo c_emlrtECI = { -1, 121, 7, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo d_emlrtECI = { -1, 117, 32, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo e_emlrtECI = { -1, 144, 4, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo f_emlrtECI = { 2, 121, 16, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo g_emlrtECI = { 2, 121, 20, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo h_emlrtECI = { -1, 120, 7, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo i_emlrtECI = { 2, 120, 16, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo j_emlrtECI = { 2, 120, 31, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo k_emlrtECI = { 2, 120, 42, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo l_emlrtECI = { -1, 119, 7, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo m_emlrtECI = { -1, 119, 26, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo n_emlrtECI = { 2, 119, 44, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo o_emlrtECI = { -1, 117, 7, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo p_emlrtECI = { 2, 117, 16, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtECInfo q_emlrtECI = { -1, 116, 33, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

static emlrtRSInfo k_emlrtRSI = { 52, "cget_latlon",
  "C:\\Users\\Brenton\\Downloads\\Google Drive\\Miscellaneous\\IGPP\\igpp\\matlabcode\\cget_latlon.m"
};

/* Function Declarations */
static void keyboard(const emlrtStack *sp, emlrtMCInfo *location);

/* Function Definitions */
static void keyboard(const emlrtStack *sp, emlrtMCInfo *location)
{
  emlrtCallMATLABR2012b(sp, 0, NULL, 0, NULL, "keyboard", true, location);
}

void cget_latlon(const emlrtStack *sp, real_T blat, real_T blon, real_T ranges,
                 real_T azm, real_T major_axis, real_T esquared, real_T *plat,
                 real_T *plon)
{
  real_T c1;
  real_T c2;
  real_T D;
  real_T P;
  real_T ss;
  real_T S;
  real_T phi1;
  int32_T ellipse;
  real_T onef;
  real_T f;
  real_T f4;
  real_T al12;
  int8_T ii_size[2];
  int32_T is1_size_idx_0;
  int32_T is1_size_idx_1;
  int8_T iv0[1];
  int32_T loop_ub;
  int32_T i0;
  int8_T iv1[1];
  real_T th1;
  real_T costh1;
  real_T sinth1;
  real_T sina12;
  int32_T im1_size_idx_0;
  int32_T im1_size_idx_1;
  int8_T iv2[1];
  int32_T im0_size_idx_0;
  int32_T im0_size_idx_1;
  int8_T iv3[1];
  real_T cosa12;
  real_T M;
  real_T N;
  real_T dv0[1];
  int8_T tmp_data[1];
  real_T b_tmp_data[1];
  int8_T tmp_size[2];
  int32_T ii[2];
  int32_T iv4[2];
  int8_T c_tmp_data[1];
  real_T d_tmp_data[1];
  real_T dv1[1];
  int32_T b_tmp_size[2];
  real_T e_tmp_data[1];
  int32_T c_tmp_size[2];
  int32_T d_tmp_size[2];
  real_T dv2[1];
  real_T dv3[1];
  real_T dv4[1];
  int8_T e_tmp_size[2];
  int8_T f_tmp_size[2];
  real_T dv5[1];
  int32_T g_tmp_size[2];
  int32_T h_tmp_size[2];
  real_T dv6[1];
  real_T s1;
  real_T d;
  real_T dv7[1];
  real_T u;
  real_T V;
  real_T sind;
  real_T ds;
  real_T cosds;
  real_T sinds;
  real_T dv8[1];
  real_T al21;
  real_T phi2;
  real_T de;
  real_T b_ellipse;
  real_T lam2;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;

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

  st.site = &k_emlrtRSI;
  keyboard(&st, &emlrtMCI);

  /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
  /*  now do the stuff that used to be in .c */
  S = ranges * 1852.0;
  phi1 = blat * 0.017453292519943295;
  ellipse = (esquared != 0.0);
  if (ellipse == 1) {
    st.site = &emlrtRSI;
    if (1.0 - esquared < 0.0) {
      b_st.site = &i_emlrtRSI;
      eml_error(&b_st);
    }

    onef = muDoubleScalarSqrt(1.0 - esquared);
    f = 1.0 - onef;
    f4 = (1.0 - onef) / 4.0;
  } else {
    onef = 1.0;
    f = 0.0;
    f4 = 0.0;
  }

  al12 = azm * 0.017453292519943295;
  st.site = &b_emlrtRSI;
  cadjlon(&st, &al12);
  if (muDoubleScalarAbs(al12) > 1.5707963267948966) {
    ii_size[0] = 1;
    ii_size[1] = 1;
  } else {
    ii_size[0] = 0;
    ii_size[1] = 0;
  }

  is1_size_idx_0 = ii_size[0];
  is1_size_idx_1 = ii_size[1];
  iv0[0] = 0;
  loop_ub = ii_size[0] * ii_size[1];
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    iv0[0] = 1;
    i0 = 1;
  }

  if (muDoubleScalarAbs(al12) <= 1.5707963267948966) {
    ii_size[0] = 1;
    ii_size[1] = 1;
  } else {
    ii_size[0] = 0;
    ii_size[1] = 0;
  }

  iv1[0] = iv0[0];
  loop_ub = ii_size[0] * ii_size[1];
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    iv1[0] = 0;
    i0 = 1;
  }

  if (ellipse == 1) {
    th1 = muDoubleScalarAtan(onef * muDoubleScalarTan(phi1));
  } else {
    th1 = phi1;
  }

  costh1 = muDoubleScalarCos(th1);
  sinth1 = muDoubleScalarSin(th1);
  sina12 = muDoubleScalarSin(al12);
  if (muDoubleScalarAbs(sina12) < 1.0E-10) {
    ii_size[0] = 1;
    ii_size[1] = 1;
  } else {
    ii_size[0] = 0;
    ii_size[1] = 0;
  }

  im1_size_idx_0 = ii_size[0];
  im1_size_idx_1 = ii_size[1];
  iv2[0] = 0;
  loop_ub = ii_size[0] * ii_size[1];
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    iv2[0] = 1;
    i0 = 1;
  }

  if (muDoubleScalarAbs(sina12) >= 1.0E-10) {
    ii_size[0] = 1;
    ii_size[1] = 1;
  } else {
    ii_size[0] = 0;
    ii_size[1] = 0;
  }

  im0_size_idx_0 = ii_size[0];
  im0_size_idx_1 = ii_size[1];
  iv3[0] = iv2[0];
  loop_ub = ii_size[0] * ii_size[1];
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    iv3[0] = 0;
    i0 = 1;
  }

  if (iv3[0] == 1) {
    sina12 = 0.0;
    if (muDoubleScalarAbs(al12) < 1.5707963267948966) {
      cosa12 = 1.0;
    } else {
      cosa12 = -1.0;
    }

    M = 0.0;
  } else {
    cosa12 = muDoubleScalarCos(al12);
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

    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 1;
      i0 = 1;
    }

    i0 = im1_size_idx_0 * im1_size_idx_1;
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    if (i0 != loop_ub) {
      emlrtSizeEqCheck1DR2012b(i0, loop_ub, &q_emlrtECI, sp);
    }

    b_tmp_data[0] = 0.0;
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      b_tmp_data[0] = 1.0 - dv0[0];
      i0 = 1;
    }

    ii_size[0] = (int8_T)im1_size_idx_0;
    ii_size[1] = (int8_T)im1_size_idx_1;
    tmp_size[0] = (int8_T)im1_size_idx_0;
    tmp_size[1] = (int8_T)im1_size_idx_1;
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 0;
      i0 = 1;
    }

    for (i0 = 0; i0 < 2; i0++) {
      ii[i0] = ii_size[i0];
      iv4[i0] = tmp_size[i0];
    }

    if ((ii[0] != iv4[0]) || (ii[1] != iv4[1])) {
      emlrtSizeEqCheckNDR2012b(ii, iv4, &p_emlrtECI, sp);
    }

    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      c_tmp_data[0] = 1;
      i0 = 1;
    }

    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      d_tmp_data[0] = b_tmp_data[0] * b_tmp_data[0];
      i0 = 1;
    }

    i0 = im1_size_idx_0 * im1_size_idx_1;
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    if (i0 != loop_ub) {
      emlrtSizeEqCheck1DR2012b(i0, loop_ub, &o_emlrtECI, sp);
    }

    dv1[0] = b_tmp_data[0];
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv1[0] = d_tmp_data[0];
      i0 = 1;
    }

    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 0;
      i0 = 1;
    }

    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      c_tmp_data[0] = 1;
      i0 = 1;
    }

    b_tmp_size[0] = im1_size_idx_0;
    b_tmp_size[1] = im1_size_idx_1;
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      b_tmp_data[0] = dv0[0];
      i0 = 1;
    }

    c_tmp_size[0] = im1_size_idx_0;
    c_tmp_size[1] = im1_size_idx_1;
    loop_ub = im1_size_idx_0 * im1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      e_tmp_data[0] = dv1[0];
      i0 = 1;
    }

    st.site = &c_emlrtRSI;
    b_rdivide(&st, b_tmp_data, b_tmp_size, e_tmp_data, c_tmp_size, d_tmp_data,
              d_tmp_size);
    i0 = im1_size_idx_0 * im1_size_idx_1;
    loop_ub = d_tmp_size[0] * d_tmp_size[1];
    if (i0 != loop_ub) {
      emlrtSizeEqCheck1DR2012b(i0, loop_ub, &d_emlrtECI, sp);
    }

    dv2[0] = 0.0;
    loop_ub = d_tmp_size[0] * d_tmp_size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      dv2[0] = d_tmp_data[i0];
    }

    /*    else  */
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 1;
      i0 = 1;
    }

    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      d_tmp_data[0] = f * M;
      i0 = 1;
    }

    i0 = im0_size_idx_0 * im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    if (i0 != loop_ub) {
      emlrtSizeEqCheck1DR2012b(i0, loop_ub, &l_emlrtECI, sp);
    }

    dv3[0] = 0.0;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv3[tmp_data[0] - 1] = d_tmp_data[0];
      i0 = 1;
    }

    c1 = dv3[0];
    ii_size[0] = (int8_T)im0_size_idx_0;
    ii_size[1] = (int8_T)im0_size_idx_1;
    tmp_size[0] = (int8_T)im0_size_idx_0;
    tmp_size[1] = (int8_T)im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 0;
      i0 = 1;
    }

    for (i0 = 0; i0 < 2; i0++) {
      ii[i0] = ii_size[i0];
      iv4[i0] = tmp_size[i0];
    }

    if ((ii[0] != iv4[0]) || (ii[1] != iv4[1])) {
      emlrtSizeEqCheckNDR2012b(ii, iv4, &n_emlrtECI, sp);
    }

    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      c_tmp_data[0] = 1;
      i0 = 1;
    }

    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      d_tmp_data[0] = f4 * (1.0 - M * M);
      i0 = 1;
    }

    i0 = im0_size_idx_0 * im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    if (i0 != loop_ub) {
      emlrtSizeEqCheck1DR2012b(i0, loop_ub, &m_emlrtECI, sp);
    }

    dv4[0] = dv0[0];
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv4[0] = d_tmp_data[0];
      i0 = 1;
    }

    c2 = dv4[0];
    ii_size[0] = (int8_T)im0_size_idx_0;
    ii_size[1] = (int8_T)im0_size_idx_1;
    tmp_size[0] = (int8_T)im0_size_idx_0;
    tmp_size[1] = (int8_T)im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 0;
      i0 = 1;
    }

    e_tmp_size[0] = (int8_T)im0_size_idx_0;
    e_tmp_size[1] = (int8_T)im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      c_tmp_data[0] = 0;
      i0 = 1;
    }

    f_tmp_size[0] = (int8_T)im0_size_idx_0;
    f_tmp_size[1] = (int8_T)im0_size_idx_1;
    for (i0 = 0; i0 < 2; i0++) {
      iv4[i0] = e_tmp_size[i0];
      ii[i0] = f_tmp_size[i0];
    }

    if ((iv4[0] != ii[0]) || (iv4[1] != ii[1])) {
      emlrtSizeEqCheckNDR2012b(iv4, ii, &k_emlrtECI, sp);
    }

    d_tmp_size[0] = im0_size_idx_0;
    d_tmp_size[1] = im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      d_tmp_data[0] = dv3[0] * M;
      i0 = 1;
    }

    for (i0 = 0; i0 < 2; i0++) {
      iv4[i0] = tmp_size[i0];
      ii[i0] = d_tmp_size[i0];
    }

    if ((iv4[0] != ii[0]) || (iv4[1] != ii[1])) {
      emlrtSizeEqCheckNDR2012b(iv4, ii, &j_emlrtECI, sp);
    }

    for (i0 = 0; i0 < 2; i0++) {
      ii[i0] = ii_size[i0];
      iv4[i0] = tmp_size[i0];
    }

    if ((ii[0] != iv4[0]) || (ii[1] != iv4[1])) {
      emlrtSizeEqCheckNDR2012b(ii, iv4, &i_emlrtECI, sp);
    }

    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      c_tmp_data[0] = 1;
      i0 = 1;
    }

    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      d_tmp_data[0] = (1.0 - dv4[0]) * ((1.0 - dv4[0]) - d_tmp_data[0]);
      i0 = 1;
    }

    i0 = im0_size_idx_0 * im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    if (i0 != loop_ub) {
      emlrtSizeEqCheck1DR2012b(i0, loop_ub, &h_emlrtECI, sp);
    }

    dv5[0] = dv1[0];
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv5[c_tmp_data[0] - 1] = d_tmp_data[0];
      i0 = 1;
    }

    D = dv5[0];
    tmp_size[0] = (int8_T)im0_size_idx_0;
    tmp_size[1] = (int8_T)im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 0;
      i0 = 1;
    }

    d_tmp_size[0] = im0_size_idx_0;
    d_tmp_size[1] = im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      d_tmp_data[0] = 0.5 * dv3[0];
      i0 = 1;
    }

    for (i0 = 0; i0 < 2; i0++) {
      iv4[i0] = d_tmp_size[i0];
      ii[i0] = tmp_size[i0];
    }

    if ((iv4[0] != ii[0]) || (iv4[1] != ii[1])) {
      emlrtSizeEqCheckNDR2012b(iv4, ii, &g_emlrtECI, sp);
    }

    ii_size[0] = (int8_T)im0_size_idx_0;
    ii_size[1] = (int8_T)im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      d_tmp_data[0] *= M;
      i0 = 1;
    }

    for (i0 = 0; i0 < 2; i0++) {
      iv4[i0] = d_tmp_size[i0];
      ii[i0] = ii_size[i0];
    }

    if ((iv4[0] != ii[0]) || (iv4[1] != ii[1])) {
      emlrtSizeEqCheckNDR2012b(iv4, ii, &f_emlrtECI, sp);
    }

    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 0;
      i0 = 1;
    }

    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      c_tmp_data[0] = 1;
      i0 = 1;
    }

    g_tmp_size[0] = im0_size_idx_0;
    g_tmp_size[1] = im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      b_tmp_data[0] = (1.0 + d_tmp_data[0]) * dv4[0];
      i0 = 1;
    }

    h_tmp_size[0] = im0_size_idx_0;
    h_tmp_size[1] = im0_size_idx_1;
    loop_ub = im0_size_idx_0 * im0_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      e_tmp_data[0] = dv5[0];
      i0 = 1;
    }

    st.site = &d_emlrtRSI;
    b_rdivide(&st, b_tmp_data, g_tmp_size, e_tmp_data, h_tmp_size, d_tmp_data,
              d_tmp_size);
    i0 = im0_size_idx_0 * im0_size_idx_1;
    loop_ub = d_tmp_size[0] * d_tmp_size[1];
    if (i0 != loop_ub) {
      emlrtSizeEqCheck1DR2012b(i0, loop_ub, &c_emlrtECI, sp);
    }

    dv6[0] = dv2[0];
    loop_ub = d_tmp_size[0] * d_tmp_size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      dv6[c_tmp_data[i0] - 1] = d_tmp_data[i0];
    }

    P = dv6[0];

    /*    end */
  }

  if (iv3[0] == 1) {
    s1 = 1.5707963267948966 - th1;
  } else {
    if (muDoubleScalarAbs(M) >= 1.0) {
      s1 = 0.0;
    } else {
      s1 = M;
      st.site = &e_emlrtRSI;
      b_acos(&st, &s1);
    }

    s1 = sinth1 / muDoubleScalarSin(s1);
    if (muDoubleScalarAbs(s1) >= 1.0) {
      s1 = 0.0;
    } else {
      st.site = &f_emlrtRSI;
      b_acos(&st, &s1);
    }
  }

  if (ellipse == 1) {
    d = rdivide(S, D * major_axis);
    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 1;
      i0 = 1;
    }

    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      d_tmp_data[0] = -d;
      i0 = 1;
    }

    i0 = is1_size_idx_0 * is1_size_idx_1;
    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    if (i0 != loop_ub) {
      emlrtSizeEqCheck1DR2012b(i0, loop_ub, &e_emlrtECI, sp);
    }

    dv7[0] = d;
    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv7[tmp_data[0] - 1] = d_tmp_data[0];
      i0 = 1;
    }

    u = 2.0 * (s1 - dv7[0]);
    V = muDoubleScalarCos(u + dv7[0]);
    sind = muDoubleScalarSin(dv7[0]);
    ds = (dv7[0] + c2 * c2 * sind * muDoubleScalarCos(dv7[0]) * (2.0 * V * V -
           1.0)) - 2.0 * P * V * (1.0 - 2.0 * P * muDoubleScalarCos(u)) * sind;
    ss = (s1 + s1) - ds;
  } else {
    ds = S / major_axis;
    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      tmp_data[0] = 1;
      i0 = 1;
    }

    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      d_tmp_data[0] = -ds;
      i0 = 1;
    }

    i0 = is1_size_idx_0 * is1_size_idx_1;
    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    if (i0 != loop_ub) {
      emlrtSizeEqCheck1DR2012b(i0, loop_ub, &b_emlrtECI, sp);
    }

    dv7[0] = ds;
    loop_ub = is1_size_idx_0 * is1_size_idx_1;
    i0 = 0;
    while (i0 <= loop_ub - 1) {
      dv7[tmp_data[0] - 1] = d_tmp_data[0];
      i0 = 1;
    }

    ds = dv7[0];
  }

  cosds = muDoubleScalarCos(ds);
  sinds = muDoubleScalarSin(ds);
  loop_ub = is1_size_idx_0 * is1_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    tmp_data[0] = 1;
    i0 = 1;
  }

  loop_ub = is1_size_idx_0 * is1_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    d_tmp_data[0] = -sinds;
    i0 = 1;
  }

  i0 = is1_size_idx_0 * is1_size_idx_1;
  loop_ub = is1_size_idx_0 * is1_size_idx_1;
  if (i0 != loop_ub) {
    emlrtSizeEqCheck1DR2012b(i0, loop_ub, &emlrtECI, sp);
  }

  dv8[0] = sinds;
  loop_ub = is1_size_idx_0 * is1_size_idx_1;
  i0 = 0;
  while (i0 <= loop_ub - 1) {
    dv8[tmp_data[0] - 1] = d_tmp_data[0];
    i0 = 1;
  }

  al21 = N * cosds - sinth1 * dv8[0];
  if (iv3[0] == 1) {
    phi2 = muDoubleScalarAtan(muDoubleScalarTan((1.5707963267948966 + s1) - ds) /
      onef);
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
    al21 = muDoubleScalarAtan(M / al21);
    if (al21 > 0.0) {
      al21 += 3.1415926535897931;
    }

    if (al12 < 0.0) {
      al21 -= 3.1415926535897931;
    }

    st.site = &g_emlrtRSI;
    cadjlon(&st, &al21);
    if (ellipse == 1) {
      b_ellipse = onef * M;
    } else {
      b_ellipse = M;
    }

    phi2 = muDoubleScalarAtan(-(sinth1 * cosds + N * dv8[0]) * muDoubleScalarSin
      (al21) / b_ellipse);
    de = muDoubleScalarAtan2(dv8[0] * sina12, costh1 * cosds - sinth1 * dv8[0] *
      cosa12);
    if (ellipse == 1) {
      if (iv1[0] == 1) {
        de += c1 * ((1.0 - c2) * ds + c2 * dv8[0] * muDoubleScalarCos(ss));
      } else {
        de -= c1 * ((1.0 - c2) * ds - c2 * dv8[0] * muDoubleScalarCos(ss));
      }
    }
  }

  lam2 = blon * 0.017453292519943295 + de;
  st.site = &h_emlrtRSI;
  cadjlon(&st, &lam2);
  *plon = lam2 * 57.295779513082323;
  *plat = phi2 * 57.295779513082323;
}

/* End of code generation (cget_latlon.c) */
