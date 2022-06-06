/* tess_stars2px.c
 * Christopher J. Burke
 * June 3, 2022
 * refactored for all cameras at once
 * Original C version
 * Alan M. Levine
 * July 27, 2018
 * Heritage: starspx9.c
 *
 * Compute predicted star positions in the focal plane and pixel coordinates.
 */

// August 3, 2017
// pixel numbers in fitpx[] start at zero (flight S/W convention)
//   -- Add one to each pixel number to obey ground FITS convention
// This addition of one is done below for consistency with ground FITS files.

// March 6, 2018 - Add columns in the output for sigx and sigy.

// March 19, 2018 - Add command line parameter to choose zero or one base for
//                  pixel numbering.
//                  zero base - numbering I use for pixels in a FITS image
//                  one base - numbering used for pixels in a ds9 display

// June 16, 2018 - Widen field definition in star_in_fov().

// June 23, 2018 - Read and process one star at a time, so there is no limit
//                 on the number of stars to be processed.

// July 27, 2018 - Add option to print out CCD pixel numbers as used by NASA Ames.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "vec.h"
#include "mat_ra3.h"

// #define NSTAR 10000

#define CCDWD_T 2048
#define CCDHT_T 2058
#define ROWA 44
#define ROWB 44
#define COLDK_T 20

#define NCAM 4  // no. of cameras
#define NCCD 4  // no. of CCDs per camera
#define NOPTCON 6

double eulcam[4][3], optcon[4][6], ccdxy0[4][4][2];
double pixsz[4][4][2], ccdang[4][4], ccdtilt[4][4][2];
double asymang[NCAM], asymfac[NCAM];

double ra_ab, dec_ab;
int nst, icam;
double rmat1[3][3], rmat2[3][3], rmat3[3][3], rmat4[3][3];
double ra_sc, dec_sc, roll_sc;

double fpxlo[4], fpxhi[4], fpylo[4], fpyhi[4];
double cpxlo[4], cpxhi[4], cpylo[4], cpyhi[4];
double dtor;

#define NSEC 69
// Hard coded spacecraft ra, dec, and rolls
double ras[NSEC] = {352.6844, 16.5571, 36.3138, 55.0070, 73.5382,\
                    92.0096,110.2559,128.1156,145.9071,\
                    165.0475,189.1247,229.5885,298.6671,\
                    276.7169,280.3985,282.4427,351.2381,\
                    16.1103, 60.2026,129.3867,171.7951,197.1008,\
                    217.2879,261.4516,265.6098,270.1381,\
                    326.8525,357.2944, 18.9190, 38.3564,\
                    57.6357, 77.1891, 96.5996,115.2951,\
                    133.2035,150.9497,170.2540,195.7176,\
                    242.1981,\
                    273.0766,277.6209, 13.0140, 49.5260,\
                    89.6066,130.2960,157.6997,143.3807,\
                    179.4254,202.6424,221.8575,239.4257,\
                    266.3618,270.8126,290.1210,307.8655,\
                    324.2778,344.2275,  9.3118, 52.9755,\
                    125.6742,118.0446,135.2412,153.0613,\
                    173.2653,201.6239,259.1702,326.7691,\
                    359.2829, 20.0449};
double decs[NSEC] = {-64.8531,-54.0160,-44.2590,-36.6420,-31.9349,\
                     -30.5839,-32.6344,-37.7370,-45.3044,\
                     -54.8165,-65.5369,-75.1256,-76.3281,\
                     62.4756, 64.0671, 66.1422, 57.8456,\
                     67.9575, 76.2343, 75.2520, 65.1924, 53.7434,\
                     43.8074, 63.1181, 61.9383, 61.5637,\
                     -72.4265,-63.0056,-52.8296,-43.3178,\
                     -35.7835,-31.3957,-30.7848,-33.7790,\
                     -39.6871,-47.7512,-57.3725,-67.8307,\
                     -76.3969,\
                     61.7450, 62.7640,  6.3337, 18.9737,\
                     24.1343, 19.0181, 10.0922, 73.1125,\
                     62.1038, 50.9532, 41.7577, 35.2333,\
                     61.8190, 61.5761, 32.6073, 37.6464,\
                     46.3448, 56.4121, 67.6524, 77.1746,\
                     77.3113,-36.0902,-42.2415,-50.6996,\
                     -60.8650,-71.5724,-78.7974,-74.2796,\
                     -64.2357,-54.2315};
double rolls[NSEC] = {222.1532,220.4335,213.0384,202.8302,191.0517,\
                      178.6367,166.4476,155.3091,145.9163,\
                      139.1724,138.0761,153.9773,198.9378,\
                      32.2329, 55.4277, 79.4699, 41.9686,\
                      40.5453, 19.6463,334.5689,317.9495,319.6992,\
                      327.4246,317.2624,339.5293,  0.6038,\
                      214.5061,222.5216,219.7970,212.0441,\
                      201.2334,188.6263,175.5369,163.1916,\
                      152.4006,143.7306,138.1685,139.3519,\
                      161.5986,\
                      14.1539, 37.2224,292.8009,284.9617,\
                      270.1557,255.0927,248.4063,327.1020,\
                      317.4166,321.3516,329.7340,339.8650,\
                      343.1429,  3.6838, 13.4565, 24.5369,\
                      36.2524, 44.0100, 45.3615, 26.5121,\
                      337.3244,162.2198,151.5884,142.7405,\
                      137.2810,140.7443,173.9147,217.4678,\
                      226.0975,222.7721};

/******************************************************************************/
// All input angles are in degrees.
// Focal plane geometry parameters are hardcoded
// using rfpg5_c#kb.txt files for parameters
// awk -v d=";" '{print $3 d}' rfpg5_c1kb.txt > temp.txt
// paste -d' ' param_load_template.txt temp.txt
// search replace icam with # of camera starting with 0
void set_fpg_parameters () {
// CAMERA 1
eulcam[0][0] =  0.101588;
eulcam[0][1] =  -36.022035;
eulcam[0][2] =  90.048315;
optcon[0][0] =  145.948116;
optcon[0][1] =  1.00000140;
optcon[0][2] =  0.24779006;
optcon[0][3] =  -0.22681254;
optcon[0][4] = 10.78243356;
optcon[0][5] =  -34.97817276;
asymang[0] =  0.00000000;
asymfac[0] =  1.00000000;
ccdxy0[0][0][0] =  31.573417;
ccdxy0[0][0][1] =  31.551637;
pixsz[0][0][0] =  0.015000;
pixsz[0][0][1] =  0.015000;
ccdang[0][0] =  179.980833;
ccdtilt[0][0][0] =  0.000000;
ccdtilt[0][0][1] = 0.000000;
ccdxy0[0][1][0] =  -0.906060;
ccdxy0[0][1][1] =  31.536148;
pixsz[0][1][0] =  0.015000;
pixsz[0][1][1] =  0.015000;
ccdang[0][1] =  180.000000;
ccdtilt[0][1][0] =  0.000000;
ccdtilt[0][1][1] = 0.000000;
ccdxy0[0][2][0] =  -31.652818;
ccdxy0[0][2][1] =  -31.438350;
pixsz[0][2][0] =  0.015000;
pixsz[0][2][1] =  0.015000;
ccdang[0][2] =  -0.024851;
ccdtilt[0][2][0] =  0.000000;
ccdtilt[0][2][1] = 0.000000;
ccdxy0[0][3][0] =  0.833161;
ccdxy0[0][3][1] =  -31.458180;
pixsz[0][3][0] =  0.015000;
pixsz[0][3][1] =  0.015000;
ccdang[0][3] =  0.001488;
ccdtilt[0][3][0] =  0.000000;
ccdtilt[0][3][1] = 0.000000;
// Camera 2
eulcam[1][0] =  -0.179412;
eulcam[1][1] =  -12.017260;
eulcam[1][2] =  90.046500;
optcon[1][0] =  145.989933;
optcon[1][1] =  1.00000140;
optcon[1][2] =  0.24069345;
optcon[1][3] =  0.15391120;
optcon[1][4] = 4.05433503;
optcon[1][5] =  3.43136895;
asymang[1] =  0.00000000;
asymfac[1] =  1.00000000;
ccdxy0[1][0][0] =  31.653635;
ccdxy0[1][0][1] =  31.470291;
pixsz[1][0][0] =  0.015000;
pixsz[1][0][1] =  0.015000;
ccdang[1][0] =  180.010890;
ccdtilt[1][0][0] =  0.000000;
ccdtilt[1][0][1] = 0.000000;
ccdxy0[1][1][0] =  -0.827405;
ccdxy0[1][1][1] =  31.491388;
pixsz[1][1][0] =  0.015000;
pixsz[1][1][1] =  0.015000;
ccdang[1][1] =  180.000000;
ccdtilt[1][1][0] =  0.000000;
ccdtilt[1][1][1] = 0.000000;
ccdxy0[1][2][0] =  -31.543794;
ccdxy0[1][2][1] =  -31.550699;
pixsz[1][2][0] =  0.015000;
pixsz[1][2][1] =  0.015000;
ccdang[1][2] =  -0.006624;
ccdtilt[1][2][0] =  0.000000;
ccdtilt[1][2][1] = 0.000000;
ccdxy0[1][3][0] =  0.922834;
ccdxy0[1][3][1] =  -31.557268;
pixsz[1][3][0] =  0.015000;
pixsz[1][3][1] =  0.015000;
ccdang[1][3] =  -0.015464;
ccdtilt[1][3][0] =  0.000000;
ccdtilt[1][3][1] = 0.000000;
// CAMERA 3
eulcam[2][0] =  0.066596;
eulcam[2][1] =  12.007750;
eulcam[2][2] =  -89.889085;
optcon[2][0] =  146.006602;
optcon[2][1] =  1.00000140;
optcon[2][2] =  0.23452229;
optcon[2][3] =  0.33552009;
optcon[2][4] = 1.92009863;
optcon[2][5] =  12.48880182;
asymang[2] =  0.00000000;
asymfac[2] =  1.00000000;
ccdxy0[2][0][0] =  31.615486;
ccdxy0[2][0][1] =  31.413644;
pixsz[2][0][0] =  0.015000;
pixsz[2][0][1] =  0.015000;
ccdang[2][0] =  179.993948;
ccdtilt[2][0][0] =  0.000000;
ccdtilt[2][0][1] = 0.000000;
ccdxy0[2][1][0] =  -0.832993;
ccdxy0[2][1][1] =  31.426621;
pixsz[2][1][0] =  0.015000;
pixsz[2][1][1] =  0.015000;
ccdang[2][1] =  180.000000;
ccdtilt[2][1][0] =  0.000000;
ccdtilt[2][1][1] = 0.000000;
ccdxy0[2][2][0] =  -31.548296;
ccdxy0[2][2][1] =  -31.606976;
pixsz[2][2][0] =  0.015000;
pixsz[2][2][1] =  0.015000;
ccdang[2][2] =  0.000298;
ccdtilt[2][2][0] =  0.000000;
ccdtilt[2][2][1] = 0.000000;
ccdxy0[2][3][0] =  0.896018;
ccdxy0[2][3][1] =  -31.569542;
pixsz[2][3][0] =  0.015000;
pixsz[2][3][1] =  0.015000;
ccdang[2][3] =  -0.006464;
ccdtilt[2][3][0] =  0.000000;
ccdtilt[2][3][1] = 0.000000;
// Camera 4
eulcam[3][0] =  0.030756;
eulcam[3][1] =  35.978116;
eulcam[3][2] =  -89.976802;
optcon[3][0] =  146.039793;
optcon[3][1] =  1.00000140;
optcon[3][2] =  0.23920416;
optcon[3][3] =  0.13349450;
optcon[3][4] = 4.77768896;
optcon[3][5] =  -1.75114744;
asymang[3] =  0.00000000;
asymfac[3] =  1.00000000;
ccdxy0[3][0][0] =  31.575820;
ccdxy0[3][0][1] =  31.316510;
pixsz[3][0][0] =  0.015000;
pixsz[3][0][1] =  0.015000;
ccdang[3][0] =  179.968217;
ccdtilt[3][0][0] =  0.000000;
ccdtilt[3][0][1] = 0.000000;
ccdxy0[3][1][0] =  -0.890877;
ccdxy0[3][1][1] =  31.363511;
pixsz[3][1][0] =  0.015000;
pixsz[3][1][1] =  0.015000;
ccdang[3][1] =  180.000000;
ccdtilt[3][1][0] =  0.000000;
ccdtilt[3][1][1] = 0.000000;
ccdxy0[3][2][0] =  -31.630470;
ccdxy0[3][2][1] =  -31.716942;
pixsz[3][2][0] =  0.015000;
pixsz[3][2][1] =  0.015000;
ccdang[3][2] =  -0.024359;
ccdtilt[3][2][0] =  0.000000;
ccdtilt[3][2][1] = 0.000000;
ccdxy0[3][3][0] =  0.824159;
ccdxy0[3][3][1] =  -31.728751;
pixsz[3][3][0] =  0.015000;
pixsz[3][3][1] =  0.015000;
ccdang[3][3] =  -0.024280;
ccdtilt[3][3][0] =  0.000000;
ccdtilt[3][3][1] = 0.000000;
}
/******************************************************************************/
// Print a 3x3 matrix

void prmat(double rm[3][3], FILE *fp)
{
  int i, j;

  fprintf(fp,"\n");
  for(i=0;i<3;++i) {
    fprintf(fp,"%d",i);
    for(j=0;j<3;++j) {
      fprintf(fp,"  %f",rm[i][j]);
    }
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
}

/******************************************************************************/
// Multiply a 3-vector in celestial coordinates by rmat1 to get a 3-vector
// in S/C coordinates.

void sky_to_sc_mat(double radecroll[3])
{
  double xeul[3];
  int i;

  xeul[0] = dtor*radecroll[0];
  xeul[1] = (M_PI/2.0) - (dtor*radecroll[1]);
  xeul[2] = dtor*radecroll[2];
  xeul[2] += M_PI;

  eulerm323(xeul,rmat1);
}

/******************************************************************************/
// Multiply a 3-vector in S/C coordinates by rmat2 to get a 3-vector
// in camera no. icam coordinates.

void sc_to_cam_mat(double euler[3])
{
  double xeul[3], angle, rm[3][3];
  int i;

  for(i=0;i<3;++i) {
    xeul[i] = dtor*euler[i];
  }

  eulerm323(xeul,rmat2);
  //prmat(rmat2,stderr);
}

/******************************************************************************/
/* return 1 if star is in the field of view of camera no. icam
 * return 0 otherwise
 */

int star_in_fov(double lngdeg, double latdeg)
{
  int ifov;
  double latr, lngr, v[3];

  lngr = dtor*lngdeg;
  latr = dtor*latdeg;

  ifov = 0;
  if(latdeg > 70.0) {
    dsphcr(lngr,latr,v);
    dnorm(v);
    // 12.0 -> 12.2 in the following 2 lines (June 16, 2018)
    if( (fabs(atan(v[0]/v[2])) <= (12.5*dtor)) &&
	(fabs(atan(v[1]/v[2])) <= (12.5*dtor)) ) {
      ifov = 1;
    }
  }

  return(ifov);
}

/******************************************************************************/
// Find new coordinates after rotating coordinate system by angle_deg

void xyrotate(double angle_deg, double xyin[2], double xyout[2])
{
  double ca, sa;

  ca = cos(dtor*angle_deg);
  sa = sin(dtor*angle_deg);
  xyout[0] = (ca*xyin[0]) + (sa*xyin[1]);
  xyout[1] = (-sa*xyin[0]) + (ca*xyin[1]);
}

/******************************************************************************/
// stretch predicted focal plane position by 'asymfac' parallel to
// azimuthal angle 'asymang'

void make_az_asym(double xyin[2], double xyout[2])
{
  double xyp[2], xypa[2];

  xyrotate(asymang[icam],xyin,xyp);
  xypa[0] = asymfac[icam]*xyp[0];
  xypa[1] = xyp[1];
  xyrotate(-asymang[icam],xypa,xyout);
}

/******************************************************************************/
/*  TSIG/AML optics model
 *
 * lng_deg, lat_deg are the spherical coordinates of the direction to the
 * star (degrees)
 *
 * xyfp[2] contains the x and y focal plane coordinates in mm
 */

void optics_fp(double lng_deg, double lat_deg, double xyfp[2])
{
  double thetar, tanth, cphi, sphi, ttsq, rfp0, rfp, xytmp[2];
  int i;

  thetar = (M_PI/2.0) - (lat_deg*dtor);
  tanth = tan(thetar);
  cphi = cos(dtor*lng_deg);
  sphi = sin(dtor*lng_deg);
  rfp0 = optcon[icam][0]*tanth;
  rfp = 0.0;
  for(i=1;i<NOPTCON;++i) {
    rfp += optcon[icam][i]*pow(tanth,2.0*(i-1));
  }
  xytmp[0] = -cphi*rfp0*rfp;
  xytmp[1] = -sphi*rfp0*rfp;
  make_az_asym(xytmp,xyfp);
}

/******************************************************************************/

void get_ra_dec_roll(double rm[3][3], double angs[3])
{
  double eul323[3];

  mateuler323(rm,eul323);  // angles in eul323 will be in radians
  angs[0] = eul323[0]/dtor;
  if (angs[0] < 0.0)
    angs[0] += 360.0;
  angs[1] = ((M_PI/2.0) - eul323[1])/dtor;
  angs[2] = eul323[2]/dtor;
}

/******************************************************************************/
/* Convert focal plane coordinates in mm to pseudo-equivalent TESS
 * camera image CCD pixel numbers and FITS file pixel numbers.
 */

// CCD tilt is ignored
// No checking is done of whether the star is on a CCD active area

// Camera number = icam + 1
// CCD number = iccd + 1

// pixel numbers in fitpx[] start at zero (flight S/W convention)
//   -- Add one to each pixel number to obey ground FITS convention

int mm_to_pix(double xy[2], double ccdpx[2], double fitpx[2])
{
  int iccd;
  double xya[2], xyb[2], xyccd[2];

  xya[0] = xy[0];
  xya[1] = xy[1];

  if(xya[0] >= 0.0) {
    if(xya[1] >= 0.0) {
      iccd = 0;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = (xyccd[0]/pixsz[icam][iccd][0]) - 0.5;
      ccdpx[1] = (xyccd[1]/pixsz[icam][iccd][1]) - 0.5;
      fitpx[0] = (CCDWD_T - ccdpx[0]) + CCDWD_T + 2*ROWA + ROWB - 1.0;
      fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0;
    }
    else {
      iccd = 3;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = (xyccd[0]/pixsz[icam][iccd][0]) - 0.5;
      ccdpx[1] = (xyccd[1]/pixsz[icam][iccd][1]) - 0.5;
      fitpx[0] = (ccdpx[0]) + CCDWD_T + 2*ROWA + ROWB;
      fitpx[1] = (ccdpx[1]);
    }
  }
  else {
    if(xya[1] >= 0.0) {
      iccd = 1;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = (xyccd[0]/pixsz[icam][iccd][0]) - 0.5;
      ccdpx[1] = (xyccd[1]/pixsz[icam][iccd][1]) - 0.5;
      fitpx[0] = (CCDWD_T - ccdpx[0]) + ROWA - 1.0;
      fitpx[1] = (CCDHT_T - ccdpx[1]) + CCDHT_T + 2*COLDK_T - 1.0;
    }
    else {
      iccd = 2;
      xyb[0] = xya[0] - ccdxy0[icam][iccd][0];
      xyb[1] = xya[1] - ccdxy0[icam][iccd][1];
      xyrotate(ccdang[icam][iccd],xyb,xyccd);
      ccdpx[0] = (xyccd[0]/pixsz[icam][iccd][0]) - 0.5;
      ccdpx[1] = (xyccd[1]/pixsz[icam][iccd][1]) - 0.5;
      fitpx[0] = (ccdpx[0]) + ROWA;
      fitpx[1] = (ccdpx[1]);
    }
  }
  /* fprintf(stderr,"mm_to_pix(): %d %f %f      %f %f      %f %f      %f %f\n",
	  iccd,xy[0],xy[1],xyccd[0],xyccd[1],ccdpx[0],ccdpx[1],fitpx[0],fitpx[1]);
  */
  return(iccd);
}

int ccdpx_to_outpx(double ccdpx[2], double napx[2], int *edgeWarn)
{
  // onSil flag indicates whether star is actually on silicon
  // edgePix - stars within edgePix set edgeWarn flag == 1
  int onSil=0, edgePix=6;
  double xUse, yUse, ymaxCoord, xmaxCoord, xMin;
  *edgeWarn = 0;

  napx[0] = -1.0;
  napx[1] = -1.0;
  xUse = ccdpx[0] + 45.0;
  yUse = ccdpx[1] + 1.0;
  xMin = 44.0;
  ymaxCoord = 2049.0;
  xmaxCoord = 2093.0;
  if (xUse>xMin && yUse>0.0 && xUse<xmaxCoord && yUse<ymaxCoord) {
    onSil = 1;
    if (xUse<=(xMin+edgePix)) {
      *edgeWarn=1;
    }
    if (xUse>=(xmaxCoord-edgePix)) {
      *edgeWarn=1;
    }
    if (yUse<=edgePix) {
      *edgeWarn=1;
    }
    if (yUse>=(ymaxCoord-edgePix)) {
      *edgeWarn=1;
    }
    napx[0] = xUse;
    napx[1] = yUse;
  }
  return(onSil);
}
/******************************************************************************/

void usage(const char *prog, int argc)
{
  fprintf(stderr,"Usage: %s RA_star Dec_star in [degrees]\n",prog);
  exit(-1);
}

/******************************************************************************/

void tess_stars2px(double ra_ab, double dec_ab)
{
  double radecroll[3], eulc[3], vstar[3], vcam[3], lng, lat, lngd, latd;
  double rdr[3];
  double xyfp[2], ccdpx[2], fitpx[2], napx[2];
  int iccd, i, jccd, kccd, j, edgeWarn, fndAny, iSec;


  dtor = M_PI/180.0;
  fndAny = 0;
  set_fpg_parameters();

  // Start loop over sectors
  for (iSec=0;iSec<NSEC;iSec++) {
    ra_sc = ras[iSec];
    dec_sc = decs[iSec];
    roll_sc = rolls[iSec];
    radecroll[0] = ra_sc;
    radecroll[1] = dec_sc;
    radecroll[2] = roll_sc;

    sky_to_sc_mat(radecroll);

    // Start loop over cameras
    for (icam=0;icam<4;icam++) {
      for(j=0;j<3;++j) {
        eulc[j] = eulcam[icam][j];
      }
      sc_to_cam_mat(eulc);
      matmat(rmat2,rmat1,rmat4);
      get_ra_dec_roll(rmat4,rdr);


      dsphcr(dtor*ra_ab,dtor*dec_ab,vstar);
      matvec(rmat4,vstar,vcam);
      dcrsph(vcam,&lng,&lat);
      lngd = lng/dtor;
      latd = lat/dtor;
      if(star_in_fov(lngd,latd) == 1) {
        optics_fp(lngd,latd,xyfp);
        iccd = mm_to_pix(xyfp,ccdpx,fitpx);
  	    if (ccdpx_to_outpx(ccdpx,napx,&edgeWarn) == 1) {
          if (fndAny == 0) {
            fprintf(stdout,"Sector Camera CCD ColPix RowPix EdgeWarn\n");
          }
          fndAny = 1;
  	      fprintf(stdout,"%3d %1d %1d %8.3f %8.3f %1d\n",
  		             iSec+1,icam+1,iccd+1,napx[0],napx[1],edgeWarn);
        }
      }

    }
  }

  if (fndAny == 0) {
    fprintf(stdout,"Star Not Observed By TESS\n");
  }

}
/******************************************************************************/
//
int main(int argc, char const *argv[]) {
  //if(argc != 3) {
  //  fprintf(stderr,"ERROR: argc = %d;\n",argc);
  //  usage(argv[0],argc);
  //}

  //ra_ab = atof(argv[1]);
  //dec_ab = atof(argv[2]);
  //fprintf(stdout, "Input ra: %11.5f dec: %11.5f\n", ra_ab, dec_ab);
  //tess_stars2px(ra_ab, dec_ab);
  //fprintf(stdout, "tess_stars2px Webassembly Ready");
  return 0;
}
