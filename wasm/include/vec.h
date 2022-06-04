/* vec.h
 * Alan M. Levine
 * March 13, 1996
 *
 * Declarations for vec.c 
 *
 * ANSI-style declarations - April 16, 1998
 */

extern void dcross(double a[3], double b[3], double c[3]);
extern void dcrsph(double vec[3], double *dra, double *ddec);
extern void dnorm(double vec[3]);
extern double dot(double a[3], double b[3]);
extern void drotate(double vec[3], double axis[3], double angle,
		    double rotvec[3]);
extern void dsphcr(double ra, double dec, double vec[3]);
extern void scale_vec(double factor, double vec1[3], double vec2[3]);
