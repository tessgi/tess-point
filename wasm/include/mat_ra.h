/* mat_ra2.h
 * Alan M. Levine
 * June 24, 1998 (from mat.h)
 * May 13, 2015
 *
 * Declarations for mat_ra2.c
 */
extern void matvec(double mat[3][3], double vin[3], double vout[3]);
extern void matmat(double mata[3][3], double matb[3][3], double matc[3][3]);
extern void rotm1(int nax, double angle, double mat[3][3]);
extern void trans(double mat[3][3], double trmat[3][3]);
extern void eulerm321(double euler[3], double rmat[3][3]);
extern void mateuler313(double rmat[3][3], double euler[3]);
extern void mateuler321(double rmat[3][3], double euler[3]);
extern double **matrix(int r, int c);
extern void free_matrix(double **m, int r);
extern void print_matrix(double **m, int r, int c, FILE *fp);
extern int **imatrix(int r, int c);
extern void free_imatrix(int **m, int r);
extern void print_imatrix(int **m, int r, int c, FILE *fp);
extern int gaussj(double **a, int n, double **b, int m);
