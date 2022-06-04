/* mat_ra3.c
 * Alan M. Levine
 * March 28, 2017 (from mat_ra2.c)
 *
 * 3x3 and other matrix routines in C
 *
 * Declarations are in mat_ra.h
 *
 * June 23, 1998 - Matrix allocation and deallocation and Gauss-Jordan 
 *                 solution of linear equations added
 * June 24, 1998 - ANSI-style declarations
 * July 31, 1998 - Add error return to gaussj()
 * March 5, 2008 - Add integer matrix allocation/deallocation
 *
 * May 13, 2015  - Make small change to mateuler313()
 *
 * March 28, 2017 - Fix description of mateuler321().
 *                  Add euler313(), eulerm123(), mateuler123().
 *                  Add euler323(), mateuler323().
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "vec.h"
#include "mat_ra.h"

#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=(temp);}

/******************************************************************************/
/* Matrix times a vector */

void matvec(double mat[3][3], double vin[3], double vout[3])
{
  int i,j;

  for(i=0;i<3;++i)
    {
      vout[i]=0.0;
      for(j=0;j<3;++j)
	vout[i] += mat[i][j]*vin[j];
    }
  return;
}

/******************************************************************************/
/* Matrix times a matrix */

void matmat(double mata[3][3], double matb[3][3], double matc[3][3])
{
  int i,j,k;

  for(i=0;i<3;++i)
    {
      for(j=0;j<3;++j)
	{
	  matc[i][j]=0.0;
	  for(k=0;k<3;++k)
	    matc[i][j] += mata[i][k]*matb[k][j];
	}
    }
  return;
}

/******************************************************************************/
/* Construct a rotation matrix that will rotate the coordinate system about
 * one of the coordinates axes by an angle (radians).
 * 
 * nax = 0, 1, or 2
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 *
 * Alternatively, the matrix rotates the vector by -angle around the given axis.
 */

void rotm1(int nax, double angle, double mat[3][3])
{
  double s;
  int i,j,n1,n2,n3;

  n1 = nax;
  n2 = (n1+1) % 3;
  n3 = (n2+1) % 3;
  s = sin(angle);

  for(i=0;i<3;++i)
    for(j=0;j<3;++j)
      mat[i][j] = 0.0;

  mat[n1][n1] = 1.0;
  mat[n3][n3] = mat[n2][n2] = cos(angle);
  mat[n2][n3] = s;
  mat[n3][n2] = -s;

  return;
}

/******************************************************************************/
/* Construct transpose of a matrix */

void trans(double mat[3][3], double trmat[3][3])
{
  int i,j;
  for(i=0;i<3;++i)
    for(j=0;j<3;++j)
      trmat[i][j] = mat[j][i];
  return;
}

/******************************************************************************/
/* Construct rotation matrix from 3-2-1 Euler angles.
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 */

void eulerm321(double euler[3], double rmat[3][3])
{
  double mat1[3][3], mat2[3][3], mata[3][3];

  rotm1(2,euler[0],mat1);
  rotm1(1,euler[1],mat2);
  matmat(mat2,mat1,mata);
  rotm1(0,euler[2],mat1);
  matmat(mat1,mata,rmat);

  return;
}

/******************************************************************************/
/* Construct rotation matrix from 3-1-3 Euler angles.
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 */

void eulerm313(double euler[3], double rmat[3][3])
{
  double mat1[3][3], mat2[3][3], mata[3][3];

  rotm1(2,euler[0],mat1);
  rotm1(0,euler[1],mat2);
  matmat(mat2,mat1,mata);
  rotm1(2,euler[2],mat1);
  matmat(mat1,mata,rmat);

  return;
}

/******************************************************************************/
/* Construct rotation matrix from 3-2-3 Euler angles.
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 */

void eulerm323(double euler[3], double rmat[3][3])
{
  double mat1[3][3], mat2[3][3], mata[3][3];

  rotm1(2,euler[0],mat1);
  rotm1(1,euler[1],mat2);
  matmat(mat2,mat1,mata);
  rotm1(2,euler[2],mat1);
  matmat(mat1,mata,rmat);

  return;
}

/******************************************************************************/
/* Construct rotation matrix from 1-2-3 Euler angles.
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 */

void eulerm123(double euler[3], double rmat[3][3])
{
  double mat1[3][3], mat2[3][3], mata[3][3];

  rotm1(0,euler[0],mat1);
  rotm1(1,euler[1],mat2);
  matmat(mat2,mat1,mata);
  rotm1(2,euler[2],mat1);
  matmat(mat1,mata,rmat);

  return;
}

/******************************************************************************/
/* Construct rotation matrix from 2-1-3 Euler angles.
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 */

void eulerm213(double euler[3], double rmat[3][3])
{
  double mat1[3][3], mat2[3][3], mata[3][3];

  rotm1(1,euler[0],mat1);
  rotm1(0,euler[1],mat2);
  matmat(mat2,mat1,mata);
  rotm1(2,euler[2],mat1);
  matmat(mat1,mata,rmat);

  return;
}

/******************************************************************************/
/* Construct 3-2-1 Euler angles from a rotation matrix
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 *
 * Matrix form is from "Spacecraft Attitude Determination and Control",
 * ed. J. Wertz
 */

void mateuler321(double rmat[3][3], double euler[3])
{
  double xlong, xlat, zp[3], zi[3];
  int i;

  for(i=0;i<3;++i)
      {
	  zp[i] = rmat[0][i];
	  zi[i] = rmat[2-i][2];
      }
  dcrsph(zp,&xlong,&xlat);
  euler[0] = xlong;
  if (euler[0] > M_PI)
      euler[0] -= 2.0*M_PI;
  if (euler[0] < -M_PI)
      euler[0] += 2.0*M_PI;
  euler[1] = -xlat;
  dcrsph(zi,&xlong,&xlat);
  euler[2] = xlong;
  if (euler[2] > M_PI)
      euler[2] -= 2.0*M_PI;
  if (euler[2] < -M_PI)
      euler[2] += 2.0*M_PI;

  return;
}

/******************************************************************************/
/* Construct 3-1-3 Euler angles from a rotation matrix
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 */

void mateuler313(double rmat[3][3], double euler[3])
{
  double xlong, xlat, zp[3], zi[3];
  int i;

  for(i=0;i<3;++i)
      {
	  zp[i] = rmat[2][i];
	  zi[i] = rmat[i][2];
      }
  dcrsph(zp,&xlong,&xlat);
  euler[0] = xlong + M_PI_2;
  if (euler[0] > M_PI)
      euler[0] -= 2.0*M_PI;
  if (euler[0] < -M_PI)
      euler[0] += 2.0*M_PI;
  euler[1] = M_PI_2 - xlat;
  dcrsph(zi,&xlong,&xlat);
  euler[2] = -xlong+M_PI_2;
  if (euler[2] > M_PI)
      euler[2] -= 2.0*M_PI;
  if (euler[2] < -M_PI)
      euler[2] += 2.0*M_PI;

  return;
}

/******************************************************************************/
/* Construct 3-2-3 Euler angles from a rotation matrix
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 *
 * Matrix form is from "Spacecraft Attitude Determination and Control",
 * ed. J. Wertz
 */

void mateuler323(double rmat[3][3], double euler[3])
{
  double xlong, xlat, zp[3], zi[3];
  int i;

  for(i=0;i<3;++i)
      {
	  zp[i] = rmat[2][i];
	  zi[i] = rmat[i][2];
      }
  dcrsph(zp,&xlong,&xlat);
  euler[0] = xlong;
  if (euler[0] > M_PI)
      euler[0] -= 2.0*M_PI;
  if (euler[0] < -M_PI)
      euler[0] += 2.0*M_PI;
  euler[1] = M_PI_2 - xlat;
  dcrsph(zi,&xlong,&xlat);
  euler[2] = -xlong+M_PI;
  if (euler[2] > M_PI)
      euler[2] -= 2.0*M_PI;
  if (euler[2] < -M_PI)
      euler[2] += 2.0*M_PI;

  return;
}

/******************************************************************************/
/* Construct 1-2-3 Euler angles from a rotation matrix
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 *
 * Matrix form is from "Spacecraft Attitude Determination and Control",
 * ed. J. Wertz
 */

void mateuler123(double rmat[3][3], double euler[3])
{
  double xlong, xlat, zp[3], zi[3];
  int i;

  for(i=0;i<3;++i)
      {
	  zp[i] = rmat[2][2-i];
	  zi[i] = rmat[i][0];
      }
  dcrsph(zp,&xlong,&xlat);
  euler[0] = -xlong;
  if (euler[0] > M_PI) 
      euler[0] -= 2.0*M_PI;
  if (euler[0] < -M_PI)
      euler[0] += 2.0*M_PI;
  euler[1] = xlat;
  dcrsph(zi,&xlong,&xlat);
  euler[2] = -xlong;
  if (euler[2] > M_PI)
      euler[2] -= 2.0*M_PI;
  if (euler[2] < -M_PI)
      euler[2] += 2.0*M_PI;

  return;
}

/******************************************************************************/
/* Construct 2-1-3 Euler angles from a rotation matrix
 *
 * The Euler angles specify the orientation of a new coordinate system
 * relative to an old coordinate system.
 *
 * The rotation matrix multiplies vectors expressed in the original coordinate
 * system to give components in the new coordinate system.
 *
 * Matrix form is from "Spacecraft Attitude Determination and Control",
 * ed. J. Wertz
 */

void mateuler213(double rmat[3][3], double euler[3])
{
  double xlong, xlat, zp[3], zi[3];
  int i;

  zp[0] = rmat[2][2];
  zp[1] = rmat[2][0];
  zp[2] = -rmat[2][1];
  zi[0] = rmat[1][1];
  zi[1] = rmat[0][1];
  zi[2] = -rmat[2][1];

  dcrsph(zp,&xlong,&xlat);
  euler[0] = xlong;
  if (euler[0] > M_PI) 
      euler[0] -= 2.0*M_PI;
  if (euler[0] < -M_PI)
      euler[0] += 2.0*M_PI;
  euler[1] = xlat;
  if (euler[1] > M_PI) 
      euler[1] -= 2.0*M_PI;
  if (euler[1] < -M_PI)
      euler[1] += 2.0*M_PI;
  dcrsph(zi,&xlong,&xlat);
  euler[2] = xlong;
  if (euler[2] > M_PI)
      euler[2] -= 2.0*M_PI;
  if (euler[2] < -M_PI)
      euler[2] += 2.0*M_PI;

  return;
}

/**************************************************************************/
/* matrix allocation/deallocation based on NRinC */

double **matrix(int r, int c)
{
    int i;
    double **m;
    
    if (!(m=(double **) malloc((unsigned)(r)*sizeof(double*))))
	fprintf(stderr, "allocation failure in matrix\n");
    
    for(i=0;i<r;i++) {
      if(!(m[i]=(double *) malloc((unsigned)(c)*sizeof(double))))
	fprintf(stderr, "allocation failure in matrix\n");
    }
    return(m);
}

/**************************************************************************/
  
void free_matrix(double **m, int r)
{ int i;
  for(i=r-1;i>=0;i--)
    free((char*)(m[i]));
  free((char*)(m));
}

/**************************************************************************/

void print_matrix(double **m, int r, int c, FILE *fp)
{
  int i, j;

  for (i=0;i<r;++i)
    {
      fprintf(fp,"%d",i);
      for(j=0;j<c;++j)
	{
	  fprintf(fp,"  %f",m[i][j]);
	}
      fprintf(fp,"\n");
    }
}

/**************************************************************************/
/* matrix allocation/deallocation based on NRinC */

int **imatrix(int r, int c)
{
    int i;
    int **m;
    
     if (!(m=(int **) malloc((unsigned)(r)*sizeof(int*))))
	fprintf(stderr, "allocation failure in imatrix\n");
    
    for(i=0;i<r;i++)
	{
	    if(!(m[i]=(int *) malloc((unsigned)(c)*sizeof(int))))
		fprintf(stderr, "allocation failure in imatrix\n");
	}
    return(m);
}

/**************************************************************************/
  
void free_imatrix(int **m, int r)
{ int i;
  for(i=r-1;i>=0;i--)
    free((char*)(m[i]));
  free((char*)(m));
}

/**************************************************************************/

void print_imatrix(int **m, int r, int c, FILE *fp)
{
  int i, j;

  for (i=0;i<r;++i)
    {
      fprintf(fp,"%d",i);
      for(j=0;j<c;++j)
	{
	  fprintf(fp,"  %d",m[i][j]);
	}
      fprintf(fp,"\n");
    }
}

/**************************************************************************/
/* gauss-jordan inversion subroutine, also from NR, somewhat awkward */

int gaussj(double **a, int n, double **b, int m)
{
  int indxc[50],indxr[50],ipiv[50];
  int i,icol,irow,j,k,l,ll;
  double big,dum,pivinv,fa;

  for(j=0;j<n;j++) ipiv[j]=0;
  for(i=0;i<n;i++) {
    big=0.0;
    for(j=0;j<n;j++)
      if (ipiv[j] !=1)
	for (k=0;k<n;k++) {
	  if (ipiv[k] == 0) {
	    if((fa=fabs(a[j][k]))>=big) {
	      big=fa;
	      irow=j;
	      icol=k;
	    }
	  } else if (ipiv[k] >1) {fprintf(stderr,"GAUSSJ1: singular matrix\n");
				  return(1);}
	}
    ++(ipiv[icol]);
    if (irow != icol) {
      for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]);
      for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l]);
    }
    indxr[i]=irow;
    indxc[i]=icol;
    if (a[icol][icol] == 0.0) {
      fprintf(stderr,"GAUSSJ2: singular matrix\n");
      return(1);
    }
    pivinv=1.0/a[icol][icol];
    a[icol][icol]=1.0;
    for (l=0;l<n;l++) a[icol][l] *= pivinv;
    for (l=0;l<m;l++) b[icol][l] *= pivinv;
      for (ll=0;ll<n;ll++)
	if(ll!=icol){
	  dum=a[ll][icol];
	  a[ll][icol]=0.0;
	  for(l=0;l<n;l++) a[ll][l] -=a[icol][l]*dum;
	  for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
	}
  }
   for (l=(n-1);l>=0;l--)
    {
      if (indxr[l] != indxc [l])
	for (k=0; k<n; k++)
	  SWAP (a[k][indxr[l]], a[k][indxc[l]]);
    }
  return(0);
}

/******************************************************************************/
