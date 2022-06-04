/* vec.c
 * Alan M. Levine
 * March 13, 1996
 *
 * Three-vector library
 *
 * Declarations are in "vec.h".
 *
 * April 16, 1998 - declarations converted to ANSI style
 */

#include <math.h>
#include "vec.h"
#define PI 3.14159265358979323846
#define DR (PI/180.0)

/*********************************************************************/
/* Cross product of 2 3-vectors - c[3] = a[3] cross b[3] */

void dcross(double a[3], double b[3], double c[3])
{
  c[0] = a[1]*b[2]-a[2]*b[1];
  c[1] = a[2]*b[0]-a[0]*b[2];
  c[2] = a[0]*b[1]-a[1]*b[0];

  return;
}

/*********************************************************************/

/* Cartesian vector to spherical coordinates conversion */
/* ra, dec are in radians */
void dcrsph(double vec[3], double *dra, double *ddec)
{
  double len, ra, dec;

  ra=0.0;
  dec=0.0;
  len = sqrt(dot(vec,vec));
  if (len > 0.0)
    {
      dec = asin(vec[2]/len);
      if ( (vec[0] != 0.0) || (vec[1] != 0.0) )
	{
	  ra = atan2(vec[1],vec[0]);
	  if (ra < 0.0) 
	    ra += 2.0*PI;
	}
    }
  *dra = ra;
  *ddec = dec;
  return;
}

/*********************************************************************/
/* Normalize a 3-vector to be of unit length - if length is non-zero */
void dnorm(double vec[3])
{
    int i;
    double len;

    len = vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
    len = sqrt(len);
    if (len > 0.0)
	{
	    for(i=0;i<3;i++)
		vec[i] = vec[i]/len;
	}
    return;
}

/*********************************************************************/
double dot(double a[3], double b[3])
{
    double dprod;
    dprod = a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
    return(dprod);
}

/*********************************************************************/
/* Rotate a vector around an axis by a given angle (in radians) */
void drotate(double vec[3], double axis[3], double angle, double rotvec[3])
{
  double vecx,vecz;
  double cang,sang,unaxis[3],y[3],x[3];
  int i;

  for(i=0;i<3;++i)
    unaxis[i] = axis[i];
  dnorm(unaxis);

  cang = cos(angle);
  sang = sin(angle);
  dcross(unaxis,vec,y);
  dnorm(y);
  dcross(y,unaxis,x);
  dnorm(x);
  vecx = dot(vec,x);
  vecz = dot(vec,unaxis);
  for(i=0;i<3;++i)
    rotvec[i] = vecz*unaxis[i]+vecx*(cang*x[i]+sang*y[i]);

  return;
}

/*********************************************************************/
void dsphcr(double ra, double dec, double vec[3])
{
    double sinra,cosra,sindec,cosdec;

    sinra = sin(ra);
    cosra = cos(ra);
    sindec = sin(dec);
    cosdec = cos(dec);

    vec[0] = cosra*cosdec;
    vec[1] = sinra*cosdec;
    vec[2] = sindec;

    return;
}

/*********************************************************************/
/* vec2 = factor * vec1 */

void scale_vec(double factor, double vec1[3], double vec2[3])
{
  vec2[0] = factor*vec1[0];
  vec2[1] = factor*vec1[1];
  vec2[2] = factor*vec1[2];
}

/*********************************************************************/

