#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

extern void nrerror (char error_text[]);
/* vector prototypes */

float *vector2 (float *, long nl, long nh);
void ez_vector(float **, long nh);
double *dvector2(double *,long nl, long nh);
void ez_dvector(double **, long nh);
int *ivector2 (int *,long nl, long nh);
void ez_ivector(int **, long nh);
long *lvector2 (long *,long nl, long nh);
void ez_lvector(long **, long nh);
unsigned int *uivector2(unsigned int *,long nl, long nh);
void ez_uivector(unsigned int **, long nh);
char *cvector2(char *,long nl, long nh);
void ez_cvector(char **, long nh);
unsigned char *ucvector2(unsigned char *,long nl, long nh);
void ez_ucvector(unsigned char **, long nh);
unsigned long *ulvector2(unsigned long *, long nl, long nh);
void ez_ulvector(unsigned long **, long nh);


/* matrix prototypes */

void ez_matrix(float ***v, long nh1, long nh2);
void ez_dmatrix(double ***v, long nh1, long nh2);
void ez_imatrix(int ***v, long nh1, long nh2);
void ez_lmatrix(long ***v, long nh1, long nh2);
void ez_cmatrix(char ***v, long nh1, long nh2);
float **matrix2(float **,long nrl, long nrh, long ncl, long nch);
double **dmatrix2(double **,long nrl, long nrh, long ncl, long nch);
int **imatrix2(int **,long nrl, long nrh, long ncl, long nch);
long **lmatrix2(long **,long nrl, long nrh, long ncl, long nch);
char **cmatrix2(char **,long nrl, long nrh, long ncl, long nch);

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch, long newrl, long newcl);
void free_submatrix(float **a, long nrl, long nrh, long ncl, long nch);
float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch);


/* tensor prototypes */

void ez_f3tensor(float ****v, long nh1, long nh2, long nh3);
void ez_d3tensor(double ****v, long nh1, long nh2, long nh3);
void ez_i3tensor(int ****v, long nh1, long nh2, long nh3);
void ez_l3tensor(long ****v, long nh1, long nh2, long nh3);
int ***i3tensor2(int ***,long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
long ***l3tensor2(long ***,long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
float ***f3tensor2(float ***,long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double ***d3tensor2(double ***,long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void deal_d3tensor(double ***t,long nrl, long ncl, long ndl);

double ***sd3dim( long d3rd, long **dims );
void free_sd3dim( double ***tofree, long d3rd, long **dims, int flag );
float ***sf3dim( long d3rd, long **dims );
void free_sf3dim( float ***tofree, long d3rd, long **dims, int flag );
int ***si3dim( long d3rd, long **dims );
void free_si3dim( int ***tofree, long d3rd, long **dims, int flag );
long ***sl3dim( long d3rd, long **dims );
void free_sl3dim( long ***tofree, long d3rd, long **dims, int flag );


/* 4-dimensional stuff prototypes */

void ez_f4dim(float *****v, long nh1, long nh2, long nh3, long nh4);
void ez_d4dim(double *****v, long nh1, long nh2, long nh3, long nh4);
void ez_i4dim(int *****v, long nh1, long nh2, long nh3, long nh4);
void ez_l4dim(long *****v, long nh1, long nh2, long nh3, long nh4);
int ****i4dim2(int ****,long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h);
long ****l4dim2(long ****,long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h);
float ****f4dim2(float ****,long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h);
double ****d4dim2(double ****,long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h);
double ****nd4dim( long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h );
void deal_nd4dim(double ****t, long d1l, long d1h, long d2l, long d3l, long d4l);
void ez_nd4dim(double *****t, long d1h_last, long d1h, long d2h, long d3h, long dh4);

double ****sd4dim( long d4th, long ***dims );
void free_sd4dim( double ****tofree, long d4th, long ***dims, int flag );
float ****sf4dim( long d4th, long ***dims );
void free_sf4dim( float ****tofree, long d4th, long ***dims, int flag );
int ****si4dim( long d4th, long ***dims );
void free_si4dim( int ****tofree, long d4th, long ***dims, int flag );
long ****sl4dim( long d4th, long ***dims );
void free_sl4dim( long ****tofree, long d4th, long ***dims, int flag );
