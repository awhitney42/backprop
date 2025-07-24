/* Public domain code as seen in Numerical Recipes in C */
#include <time.h>
#include <math.h>
#include "ez_alloc.h"

#define NR_END 0
#define FREE_ARG char *

/* vector allocation & deallocation */

/*****************************************************************************/
float *vector2(float *v, long nl, long nh)
/* allocate a float vector with range v[nl..nh] */
{
    v = (float *)realloc(v, (size_t)((nh - nl + 1 + NR_END) * sizeof(float)));
    if (!v)
        nrerror("allocation failure in vector2()");
    return v - nl + NR_END;
}
void ez_vector(float **v, long nh)
{
    *v = vector2(*v, 0, nh);
}

/*****************************************************************************/

double *dvector2(double *v, long nl, long nh)
/* allocate a double vector with range v[nl..nh] */
{
    v = (double *)realloc(v, (size_t)((nh - nl + 1 + NR_END) * sizeof(double)));
    if (!v)
        nrerror("allocation failure in dvector2()");
    return v - nl + NR_END;
}
void ez_dvector(double **v, long nh)
{
    *v = dvector2(*v, 0, nh);
}

/*****************************************************************************/

int *ivector2(int *v, long nl, long nh)
/* allocate an integer vector with range v[nl..nh] */
{
    v = (int *)realloc(v, (size_t)((nh - nl + 1 + NR_END) * sizeof(int)));
    if (!v)
        nrerror("allocation failure in ivector2()");
    return v - nl + NR_END;
}
void ez_ivector(int **v, long nh)
{
    *v = ivector2(*v, 0, nh);
}

/*****************************************************************************/

long *lvector2(long *v, long nl, long nh)
/* allocate a long integer vector with range v[nl..nh] */
{
    v = (long *)realloc(v, (size_t)((nh - nl + 1 + NR_END) * sizeof(long)));
    if (!v)
        nrerror("allocation failure in lvector2()");
    return v - nl + NR_END;
}
void ez_lvector(long **v, long nh)
{
    *v = lvector2(*v, 0, nh);
}

/*****************************************************************************/

unsigned int *uivector2(unsigned int *v, long nl, long nh)
/* allocate an unsigned int vector with range v[nl..nh] */
{
    v = (unsigned int *)realloc(v, (size_t)((nh - nl + 1 + NR_END) * sizeof(unsigned int)));
    if (!v)
        nrerror("allocation failure in uivector2()");
    return v - nl + NR_END;
}
void ez_uivector(unsigned int **v, long nh)
{
    *v = uivector2(*v, 0, nh);
}

/*****************************************************************************/

unsigned char *ucvector2(unsigned char *v, long nl, long nh)
/* allocates an unsigned char vector with range v[nl..nh] */
{
    v = (unsigned char *)realloc(v, (size_t)((nh - nl + 1 + NR_END) * sizeof(unsigned char)));
    if (!v)
        nrerror("allocation failure in ucvector2()");
    return v - nl + NR_END;
}
void ez_ucvector(unsigned char **v, long nh)
{
    *v = ucvector2(*v, 0, nh);
}

/*****************************************************************************/

char *cvector2(char *v, long nl, long nh)
/* allocates an char vector with range v[nl..nh] */
{
    v = (char *)realloc(v, (size_t)((nh - nl + 1 + NR_END) * sizeof(char)));
    if (!v)
        nrerror("allocation failure in cvector2()");
    return v - nl + NR_END;
}
void ez_cvector(char **v, long nh)
{
    *v = cvector2(*v, 0, nh);
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/* matrix allocation & deallocation */

void ez_imatrix(int ***v, long nh1, long nh2)
{
    *v = imatrix2(*v, 0, nh1, 0, nh2);
}

void ez_lmatrix(long ***v, long nh1, long nh2)
{
    *v = lmatrix2(*v, 0, nh1, 0, nh2);
}

void ez_matrix(float ***v, long nh1, long nh2)
{
    *v = matrix2(*v, 0, nh1, 0, nh2);
}

void ez_dmatrix(double ***v, long nh1, long nh2)
{
    *v = dmatrix2(*v, 0, nh1, 0, nh2);
}

void ez_cmatrix(char ***v, long nh1, long nh2)
{
    *v = cmatrix2(*v, 0, nh1, 0, nh2);
}

/*****************************************************************************/

int **imatrix2(int **m, long nrl, long nrh, long ncl, long nch)
/* allocate a matrix with range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;

    if (m != NULL)
    {
        free((FREE_ARG)(m[nrl] + ncl - NR_END));
        free((FREE_ARG)(m + nrl - NR_END));
    }

    /* allocate pointers to rows */
    m = (int **)malloc((size_t)((nrow + NR_END) * sizeof(int *)));
    if (!m)
    {
        fprintf(stderr, "d1 = %ld\n", nrow);
        nrerror("allocation failure 1 in imatrix2()");
    }
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (int *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(int)));
    if (!m[nrl])
    {
        fprintf(stderr, "d1 = %ld\td1 = %ld\n", nrow, ncol);
        nrerror("allocation failure 2 in imatrix ()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return a pointer to array of pointers to rows */
    return m;
}

long **lmatrix2(long **m, long nrl, long nrh, long ncl, long nch)
/* allocate a matrix with range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    if (m != NULL)
    {
        free((FREE_ARG)(m[nrl] + ncl - NR_END));
        free((FREE_ARG)(m + nrl - NR_END));
    }
    /* allocate pointers to rows */
    m = (long **)malloc((size_t)((nrow + NR_END) * sizeof(long *)));
    if (!m)
    {
        fprintf(stderr, "d1 = %ld\n", nrow);
        nrerror("allocation failure 1 in lmatrix2()");
    }
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (long *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(long)));
    if (!m[nrl])
    {
        fprintf(stderr, "d1 = %ld\td1 = %ld\n", nrow, ncol);
        nrerror("allocation failure 2 in lmatrix ()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return a pointer to array of pointers to rows */
    return m;
}

float **matrix2(float **m, long nrl, long nrh, long ncl, long nch)
/* allocate a matrix with range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
    if (m != NULL)
    {
        free((FREE_ARG)(m[nrl] + ncl - NR_END));
        free((FREE_ARG)(m + nrl - NR_END));
    }

    /* allocate pointers to rows */
    m = (float **)malloc((size_t)((nrow + NR_END) * sizeof(float *)));
    if (!m)
    {
        fprintf(stderr, "d1 = %ld\n", nrow);
        nrerror("allocation failure 1 in matrix2()");
    }
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (float *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float)));
    if (!m[nrl])
    {
        fprintf(stderr, "d1 = %ld\td1 = %ld\n", nrow, ncol);
        nrerror("allocation failure 2 in matrix ()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return a pointer to array of pointers to rows */
    return m;
}

double **dmatrix2(double **m, long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;

    if (m != NULL)
    {
        free((FREE_ARG)(m[nrl] + ncl - NR_END));
        free((FREE_ARG)(m + nrl - NR_END));
    }

    /* allocate pointers to rows */
    m = (double **)malloc((size_t)((nrow + NR_END) * sizeof(double *)));
    if (!m)
    {
        fprintf(stderr, "d1 = %ld\n", nrow);
        nrerror("allocation failure 1 in dmatrix2()");
    }
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (double *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double)));
    if (!m[nrl])
    {
        fprintf(stderr, "d1 = %ld\td1 = %ld\n", nrow, ncol);
        nrerror("allocation failure 2 in dmatrix ()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return a pointer to array of pointers to rows */
    return m;
}

char **cmatrix2(char **m, long nrl, long nrh, long ncl, long nch)
/* allocate a character matrix with range m[nrl..nrh][ncl..nch] */
{
    long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;

    if (m != NULL)
    {
        free((FREE_ARG)(m[nrl] + ncl - NR_END));
        free((FREE_ARG)(m + nrl - NR_END));
    }

    /* allocate pointers to rows */
    m = (char **)malloc((size_t)((nrow + NR_END) * sizeof(char *)));
    if (!m)
    {
        fprintf(stderr, "d1 = %ld\n", nrow);
        nrerror("allocation failure 1 in cmatrix2()");
    }
    m += NR_END;
    m -= nrl;

    /* allocate rows and set pointers to them */
    m[nrl] = (char *)malloc((size_t)((nrow * ncol + NR_END) * sizeof(char)));
    if (!m[nrl])
    {
        fprintf(stderr, "d1 = %ld\td1 = %ld\n", nrow, ncol);
        nrerror("allocation failure 2 in cmatrix ()");
    }
    m[nrl] += NR_END;
    m[nrl] -= ncl;

    for (i = nrl + 1; i <= nrh; i++)
        m[i] = m[i - 1] + ncol;

    /* return a pointer to array of pointers to rows */
    return m;
}
/* Public domain code as seen in Numerical Recipes in C */

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/* tensor (3-dimensional object) prototypes */

void ez_i3tensor(int ****v, long nh1, long nh2, long nh3)
{
    *v = i3tensor2(*v, 0, nh1, 0, nh2, 0, nh3);
}

void ez_l3tensor(long ****v, long nh1, long nh2, long nh3)
{
    *v = l3tensor2(*v, 0, nh1, 0, nh2, 0, nh3);
}

void ez_f3tensor(float ****v, long nh1, long nh2, long nh3)
{
    *v = f3tensor2(*v, 0, nh1, 0, nh2, 0, nh3);
}

void ez_d3tensor(double ****v, long nh1, long nh2, long nh3)
{
    *v = d3tensor2(*v, 0, nh1, 0, nh2, 0, nh3);
}

/*****************************************************************************/

float ***f3tensor2(float ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;

    if (t != NULL)
    {
        free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
        free((FREE_ARG)(t[nrl] + ncl - NR_END));
        free((FREE_ARG)(t + nrl - NR_END));
    }

    /* allocate pointers to pointers to rows */
    t = (float ***)malloc((size_t)((nrow + NR_END) * sizeof(float **)));
    if (!t)
        nrerror("allocation failure 1 in f3tensor2()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (float **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(float *)));
    if (!t[nrl])
        nrerror("allocation failure 2 in f3tensor2()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (float *)malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(float)));
    if (!t[nrl][ncl])
        nrerror("allocation failure 3 in f3tensor2()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ncl;

    for (j = ncl + 1; j < nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i < nrh; i++)
    {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j < nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

int ***i3tensor2(int ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a int 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;

    if (t != NULL)
    {
        free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
        free((FREE_ARG)(t[nrl] + ncl - NR_END));
        free((FREE_ARG)(t + nrl - NR_END));
    }

    /* allocate pointers to pointers to rows */
    t = (int ***)malloc((size_t)((nrow + NR_END) * sizeof(int **)));
    if (!t)
        nrerror("allocation failure 1 in f3tensor2()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (int **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(int *)));
    if (!t[nrl])
        nrerror("allocation failure 2 in f3tensor2()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (int *)malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(int)));
    if (!t[nrl][ncl])
        nrerror("allocation failure 3 in f3tensor2()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ncl;

    for (j = ncl + 1; j < nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i < nrh; i++)
    {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j < nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

long ***l3tensor2(long ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a long 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;

    if (t != NULL)
    {
        free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
        free((FREE_ARG)(t[nrl] + ncl - NR_END));
        free((FREE_ARG)(t + nrl - NR_END));
    }

    /* allocate pointers to pointers to rows */
    t = (long ***)malloc((size_t)((nrow + NR_END) * sizeof(long **)));
    if (!t)
        nrerror("allocation failure 1 in l3tensor2()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (long **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(long *)));
    if (!t[nrl])
        nrerror("allocation failure 2 in l3tensor2()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (long *)malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(long)));
    if (!t[nrl][ncl])
        nrerror("allocation failure 3 in l3tensor2()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ncl;

    for (j = ncl + 1; j < nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i < nrh; i++)
    {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j < nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

double ***d3tensor2(double ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;

    if (t != NULL)
    {
        free((FREE_ARG)(t[nrl][ncl] + ndl - NR_END));
        free((FREE_ARG)(t[nrl] + ncl - NR_END));
        free((FREE_ARG)(t + nrl - NR_END));
    }

    /* allocate pointers to pointers to rows */
    t = (double ***)malloc((size_t)((nrow + NR_END) * sizeof(double **)));
    if (!t)
        nrerror("allocation failure 1 in d3tensor2()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (double **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double *)));
    if (!t[nrl])
        nrerror("allocation failure 2 in d3tensor2()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (double *)malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(double)));
    if (!t[nrl][ncl])
        nrerror("allocation failure 3 in d3tensor2()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ncl;

    for (j = ncl + 1; j < nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i < nrh; i++)
    {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j < nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

double dsqrarg;

/* 4-dimensional object allocation routines, compliments of Rodney Schuler */

void ez_i4dim(int *****v, long nh1, long nh2, long nh3, long nh4)
{
    *v = i4dim2(*v, 0, nh1, 0, nh2, 0, nh3, 0, nh4);
}

void ez_l4dim(long *****v, long nh1, long nh2, long nh3, long nh4)
{
    *v = l4dim2(*v, 0, nh1, 0, nh2, 0, nh3, 0, nh4);
}

void ez_f4dim(float *****v, long nh1, long nh2, long nh3, long nh4)
{
    *v = f4dim2(*v, 0, nh1, 0, nh2, 0, nh3, 0, nh4);
}

void ez_d4dim(double *****v, long nh1, long nh2, long nh3, long nh4)
{
    *v = d4dim2(*v, 0L, nh1, 0L, nh2, 0L, nh3, 0L, nh4);
}

/*****************************************************************************/

float ****f4dim2(float ****t, long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h)
/* allocate a 4 dimensional array of floats with
range t[d1l..d1h][d2l..d2h][d3l..d3h][d4l..d4h] */
{
    long i, j, k;
    long gap1 = d1h - d1l + 1, gap2 = d2h - d2l + 1, gap3 = d3h - d3l + 1, gap4 = d4h - d4l + 1;
    float ***prev1, **prev2, *prev3;

    if (t != NULL)
    {
        free((FREE_ARG)(t[d1l][d2l][d3l] + d4l - NR_END));
        free((FREE_ARG)(t[d1l][d2l] + d3l - NR_END));
        free((FREE_ARG)(t[d1l] + d2l - NR_END));
        free((FREE_ARG)(t + d1l - NR_END));
    }

    /* allocate pointers to pointers to pointers to floats */
    t = (float ****)malloc((size_t)((gap1 + NR_END) * sizeof(float ***)));
    if (!t)
        nrerror("allocation failure 1 in f4tensor()");
    t += NR_END;
    t -= d1l;

    /* allocate pointers to pointers to floats */
    t[d1l] = (float ***)malloc((size_t)((gap1 * gap2 + NR_END) * sizeof(float **)));
    if (!t[d1l])
        nrerror("allocation failure 2 in f4tensor()");
    t[d1l] += NR_END;
    t[d1l] -= d2l;

    /* allocate pointers to floats */
    t[d1l][d2l] = (float **)malloc((size_t)((gap1 * gap2 * gap3 + NR_END) * sizeof(float *)));
    if (!t[d1l][d2l])
        nrerror("allocation failure 3 in f4tensor()");
    t[d1l][d2l] += NR_END;
    t[d1l][d2l] -= d3l;

    /* allocate floats */
    t[d1l][d2l][d3l] = (float *)malloc((size_t)((gap1 * gap2 * gap3 * gap4 + NR_END) * sizeof(float)));
    if (!t[d1l][d2l][d3l])
        nrerror("allocation failure 4 in f4tensor()");
    t[d1l][d2l][d3l] += NR_END;
    t[d1l][d2l][d3l] -= d4l;

    /* set first dim pointers */
    prev1 = t[d1l] - gap2;
    for (i = d1l; i <= d1h; i++)
    {
        t[i] = prev1 + gap2;
        prev1 = t[i];
    }

    /* set 2nd dim pointers */
    prev2 = t[d1l][d2l] - gap3;
    for (i = d1l; i <= d1h; i++)
        for (j = d2l; j <= d2h; j++)
        {
            t[i][j] = prev2 + gap3;
            prev2 = t[i][j];
        }

    /* set 3rd dim pointers */
    prev3 = t[d1l][d2l][d3l];
    for (i = d1l; i <= d1h; i++)
        for (j = d2l; j <= d2h; j++)
            for (k = d3l; k <= d3h; k++)
            {
                t[i][j][k] = prev3 + gap4;
                prev3 = t[i][j][k];
            }

    return t;
}

double ****d4dim2(double ****t, long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h)
{
    long i, j, k;
    long gap1 = d1h - d1l + 1, gap2 = d2h - d2l + 1, gap3 = d3h - d3l + 1, gap4 = d4h - d4l + 1;

    if (t != NULL)
    {
        free((FREE_ARG)(t[d1l][d2l][d3l] + d4l - NR_END));
        free((FREE_ARG)(t[d1l][d2l] + d3l - NR_END));
        free((FREE_ARG)(t[d1l] + d2l - NR_END));
        free((FREE_ARG)(t + d1l - NR_END));
    }

    t = (double ****)malloc((size_t)((gap1 + NR_END) * sizeof(double ***)));
    if (!t)
        nrerror("allocation failure 1 in d4tensor()");
    t += NR_END;
    t -= d1l;

    t[d1l] = (double ***)malloc((size_t)((gap1 * gap2 + NR_END) * sizeof(double **)));
    if (!t[d1l])
        nrerror("allocation failure 2 in d4tensor()");
    t[d1l] += NR_END;
    t[d1l] -= d2l;

    t[d1l][d2l] = (double **)malloc((size_t)((gap1 * gap2 * gap3 + NR_END) * sizeof(double *)));
    if (!t[d1l][d2l])
        nrerror("allocation failure 3 in d4tensor()");
    t[d1l][d2l] += NR_END;
    t[d1l][d2l] -= d3l;

    t[d1l][d2l][d3l] = (double *)malloc((size_t)((gap1 * gap2 * gap3 * gap4 + NR_END) * sizeof(double)));
    if (!t[d1l][d2l][d3l])
        nrerror("allocation failure 4 in d4tensor()");
    t[d1l][d2l][d3l] += NR_END;
    t[d1l][d2l][d3l] -= d4l;

    for (k = d3l + 1; k < d3h; k++)
        t[d1l][d2l][k] = t[d1l][d2l][k - 1] + gap4;

    for (j = d2l + 1; j < d2h; j++)
        t[d1l][j] = t[d1l][j - 1] + gap3 * gap4;

    for (i = d1l + 1; i < d1h; i++)
    {
        t[i] = t[i - 1] + gap2;
        t[i][d2l] = t[i - 1][d2l] + gap2 * gap3;
        for (j = d2l + 1; j < d2h; j++)
        {
            t[i][j] = t[i][j - 1] + gap3;
            t[i][j][d3l] = t[i][j - 1][d3l] + gap3 * gap4;
            for (k = d3l + 1; k < d3h; k++)
                t[i][j][k] = t[i][j][k - 1] + gap4;
        }
    }

    return t;
}

int ****i4dim2(int ****t, long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h)
/* allocate a 4 dimensional array of ints with
range t[d1l..d1h][d2l..d2h][d3l..d3h][d4l..d4h] */
{
    long i, j, k;
    long gap1 = d1h - d1l + 1, gap2 = d2h - d2l + 1, gap3 = d3h - d3l + 1, gap4 = d4h - d4l + 1;
    int ***prev1, **prev2, *prev3;

    if (t != NULL)
    {
        free((FREE_ARG)(t[d1l][d2l][d3l] + d4l - NR_END));
        free((FREE_ARG)(t[d1l][d2l] + d3l - NR_END));
        free((FREE_ARG)(t[d1l] + d2l - NR_END));
        free((FREE_ARG)(t + d1l - NR_END));
    }

    /* allocate pointers to pointers to pointers to ints */
    t = (int ****)malloc((size_t)((gap1 + NR_END) * sizeof(int ***)));
    if (!t)
        nrerror("allocation failure 1 in i4dim2()");
    t += NR_END;
    t -= d1l;

    /* allocate pointers to pointers to ints */
    t[d1l] = (int ***)malloc((size_t)((gap1 * gap2 + NR_END) * sizeof(int **)));
    if (!t[d1l])
        nrerror("allocation failure 2 in i4dim2()");
    t[d1l] += NR_END;
    t[d1l] -= d2l;

    /* allocate pointers to ints */
    t[d1l][d2l] = (int **)malloc((size_t)((gap1 * gap2 * gap3 + NR_END) * sizeof(int *)));
    if (!t[d1l][d2l])
        nrerror("allocation failure 3 in i4dim2()");
    t[d1l][d2l] += NR_END;
    t[d1l][d2l] -= d3l;

    /* allocate ints */
    t[d1l][d2l][d3l] = (int *)malloc((size_t)((gap1 * gap2 * gap3 * gap4 + NR_END) * sizeof(int)));
    if (!t[d1l][d2l][d3l])
        nrerror("allocation failure 4 in i4dim2()");
    t[d1l][d2l][d3l] += NR_END;
    t[d1l][d2l][d3l] -= d4l;

    /* set first dim pointers */
    prev1 = t[d1l] - gap2;
    for (i = d1l; i <= d1h; i++)
    {
        t[i] = prev1 + gap2;
        prev1 = t[i];
    }

    /* set 2nd dim pointers */
    prev2 = t[d1l][d2l] - gap3;
    for (i = d1l; i <= d1h; i++)
        for (j = d2l; j <= d2h; j++)
        {
            t[i][j] = prev2 + gap3;
            prev2 = t[i][j];
        }

    /* set 3rd dim pointers */
    prev3 = t[d1l][d2l][d3l];
    for (i = d1l; i <= d1h; i++)
        for (j = d2l; j <= d2h; j++)
            for (k = d3l; k <= d3h; k++)
            {
                t[i][j][k] = prev3 + gap4;
                prev3 = t[i][j][k];
            }

    return t;
}

long ****l4dim2(long ****t, long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h)
/* allocate a 4 dimensional array of longs with
range t[d1l..d1h][d2l..d2h][d3l..d3h][d4l..d4h] */
{
    long i, j, k;
    long gap1 = d1h - d1l + 1, gap2 = d2h - d2l + 1, gap3 = d3h - d3l + 1, gap4 = d4h - d4l + 1;
    long ***prev1, **prev2, *prev3;

    if (t != NULL)
    {
        free((FREE_ARG)(t[d1l][d2l][d3l] + d4l - NR_END));
        free((FREE_ARG)(t[d1l][d2l] + d3l - NR_END));
        free((FREE_ARG)(t[d1l] + d2l - NR_END));
        free((FREE_ARG)(t + d1l - NR_END));
    }

    /* allocate pointers to pointers to pointers to longs */
    t = (long ****)malloc((size_t)((gap1 + NR_END) * sizeof(long ***)));
    if (!t)
        nrerror("allocation failure 1 in l4dim2()");
    t += NR_END;
    t -= d1l;

    /* allocate pointers to pointers to longs */
    t[d1l] = (long ***)malloc((size_t)((gap1 * gap2 + NR_END) * sizeof(long **)));
    if (!t[d1l])
        nrerror("allocation failure 2 in l4dim2()");
    t[d1l] += NR_END;
    t[d1l] -= d2l;

    /* allocate pointers to longs */
    t[d1l][d2l] = (long **)malloc((size_t)((gap1 * gap2 * gap3 + NR_END) * sizeof(long *)));
    if (!t[d1l][d2l])
        nrerror("allocation failure 3 in l4dim2()");
    t[d1l][d2l] += NR_END;
    t[d1l][d2l] -= d3l;

    /* allocate longs */
    t[d1l][d2l][d3l] = (long *)malloc((size_t)((gap1 * gap2 * gap3 * gap4 + NR_END) * sizeof(long)));
    if (!t[d1l][d2l][d3l])
        nrerror("allocation failure 4 in l4dim2()");
    t[d1l][d2l][d3l] += NR_END;
    t[d1l][d2l][d3l] -= d3l;

    /* set first dim pointers */
    prev1 = t[d1l] - gap2;
    for (i = d1l; i <= d1h; i++)
    {
        t[i] = prev1 + gap2;
        prev1 = t[i];
    }

    /* set 2nd dim pointers */
    prev2 = t[d1l][d2l] - gap3;
    for (i = d1l; i <= d1h; i++)
        for (j = d2l; j <= d2h; j++)
        {
            t[i][j] = prev2 + gap3;
            prev2 = t[i][j];
        }

    /* set 3rd dim pointers */
    prev3 = t[d1l][d2l][d3l];
    for (i = d1l; i <= d1h; i++)
        for (j = d2l; j <= d2h; j++)
            for (k = d3l; k <= d3h; k++)
            {
                t[i][j][k] = prev3 + gap4;
                prev3 = t[i][j][k];
            }

    return t;
}
/******************************************************************************/
double ****nd4dim(long d1l, long d1h, long d2l, long d2h, long d3l, long d3h, long d4l, long d4h)
{
    long i, d1st;
    double ****t;

    d1st = d1h - d1l + 1;

    t = (double ****)malloc((size_t)((d1st + NR_END) * sizeof(double ***)));
    if (!t)
    {
        fprintf(stderr, "d1 = %ld\n", d1st);
        nrerror("allocation failure 1 in nd4dim()");
    }
    t += NR_END;

    for (i = 0; i < d1st; i++)
        t[i] = d3tensor(d2l, d2h, d3l, d3h, d4l, d4h);

    return t;
}

void deal_nd4dim(double ****t, long d1l, long d1h, long d2l, long d3l, long d4l)
{
    long i, d1st;

    d1st = d1h - d1l + 1;

    for (i = 0; i < d1st; i++)
        deal_d3tensor(t[i], d2l, d3l, d4l);

    free((FREE_ARG)(t + d1l - NR_END));
}

void ez_nd4dim(double *****t, long d1h_last, long d1h, long d2h, long d3h, long d4h)
{
    if ((*t) != NULL)
        deal_nd4dim(*t, 0L, d1h_last, 0L, 0L, 0L);
    *t = nd4dim(0L, d1h, 0L, d2h, 0L, d3h, 0L, d4h);
}
/******************************************************************************/
double ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
    double ***t;

    /* allocate pointers to pointers to rows */
    t = (double ***)malloc((size_t)((nrow + NR_END) * sizeof(double **)));
    if (!t)
        nrerror("allocation failure 1 in d3tensor()");
    t += NR_END;
    t -= nrl;

    /* allocate pointers to rows and set pointers to them */
    t[nrl] = (double **)malloc((size_t)((nrow * ncol + NR_END) * sizeof(double *)));
    if (!t[nrl])
        nrerror("allocation failure 2 in d3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;

    /* allocate rows and set pointers to them */
    t[nrl][ncl] = (double *)malloc((size_t)((nrow * ncol * ndep + NR_END) * sizeof(double)));
    if (!t[nrl][ncl])
        nrerror("allocation failure 3 in d3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ncl;

    for (j = ncl + 1; j < nch; j++)
        t[nrl][j] = t[nrl][j - 1] + ndep;
    for (i = nrl + 1; i < nrh; i++)
    {
        t[i] = t[i - 1] + ncol;
        t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
        for (j = ncl + 1; j < nch; j++)
            t[i][j] = t[i][j - 1] + ndep;
    }

    /* return pointer to array of pointers to rows */
    return t;
}
void deal_d3tensor(double ***t, long d2l, long d3l, long d4l)
{
    free((FREE_ARG)(t[d2l][d3l] + d4l - NR_END));
    free((FREE_ARG)(t[d2l] + d3l - NR_END));
    free((FREE_ARG)(t + d2l - NR_END));
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

void nrerror(char error_text[])
/* Numerical Recipies standard error handler */
{
    fprintf(stderr, "Numerical Recipies run-time error...\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}
