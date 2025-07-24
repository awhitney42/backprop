/* ***************************************************************** */
/* datafile.c                                                        */
/*                                                                   */
/* Adam P. Whitney                                                   */
/* ***************************************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "datafile.h"
#include "ez_alloc.h"

void dfile_info(char *filename, long *npatterns, long *ncolumns, long *linelength)
{

    FILE *fp = NULL;
    char *buffer = NULL;
    long id = 0, i, llTemp = 0;
    unsigned char c;

    ez_cvector(&buffer, MAX_BUFFER);

    (*npatterns) = 0;
    (*ncolumns) = 0;
    fp = fopen(filename, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Error opening file %s\n", filename);
        fflush(stderr);
        exit(1);
    }

    /* Determine the line length of the first data line. */
    id = 1;
    fgets(buffer, (MAX_BUFFER - 1), fp);
    while ((!feof(fp)) && (strchr(buffer, '#') != NULL))
    {
        fgets(buffer, (MAX_BUFFER - 1), fp);
        id++;
    }

    rewind(fp);

    for (i = 1; i < id; i++)
    {
        fgets(buffer, (MAX_BUFFER - 1), fp);
    }

    c = '\0';
    llTemp = 0;
    while (c != '\n')
    {
        c = fgetc(fp);
        llTemp++;
    }

    if (llTemp > MAX_BUFFER)
    {

        (*linelength) = llTemp + (long)((double)llTemp * 0.10);
        ez_cvector(&buffer, (*linelength));
    }
    else
    {

        (*linelength) = MAX_BUFFER;
    }

    rewind(fp);

    /* *** */

    while (!feof(fp))
    {
        fgets(buffer, ((*linelength) - 1), fp);
        if ((strchr(buffer, '#') == NULL) && (numerical(buffer)))
        {
            (*npatterns)++;
            if ((*ncolumns) == 0)
            {
                if (strtok(buffer, " \t\n") != NULL)
                {
                    (*ncolumns)++;
                    while (strtok(NULL, " \t\n") != NULL)
                    {
                        (*ncolumns)++;
                    }
                }
            }
        }
        buffer[0] = '\0';
    }
    fclose(fp);

    ez_cvector(&buffer, 0);

    return;
}

void dfile_read(char *filename, double ***d, long *npatterns, long *ncolumns, int printStatus)
{

    FILE *fp = NULL;
    int i, j;
    int reportFlag[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, reportCount = 1;
    char *buffer = NULL;
    int p;
    long linelength;

    dfile_info(filename, npatterns, ncolumns, &linelength);

    ez_cvector(&buffer, linelength);

    if (printStatus)
    {

        sprintf(buffer, "Reading %ld patterns in %s\r", (*npatterns), filename);
        p = strlen(buffer);
        fprintf(stdout, "%s", buffer);
        fflush(stdout);

        ez_dmatrix(d, (long)(*npatterns), (long)(*ncolumns));

        fp = fopen(filename, "r");
        if (fp == NULL)
        {
            fprintf(stderr, "Error opening file %s\n", filename);
            fflush(stderr);
            exit(1);
        }

        for (i = 0; i < (*npatterns); i++)
        {

            dfile_nextpattern(fp, ((*d)[i]), buffer, (*ncolumns), linelength);

            if (((((double)(i + 1) / (double)(*npatterns)) * 100) >= (10 * reportCount)) && (reportFlag[reportCount - 1] == 0))
            {
                sprintf(buffer, "Read %d%%", (reportCount * 10));
                fprintf(stdout, "%s", buffer);
                for (j = 0; j < (p - (int)strlen(buffer)); j++)
                    fprintf(stdout, " ");
                fprintf(stdout, "\r");
                fflush(stdout);
                reportFlag[reportCount - 1] = 1;
                reportCount++;
            }
        }

        fprintf(stdout, "Read %ld patterns in %s    \n", (*npatterns), filename);
    }
    else
    {

        ez_dmatrix(d, (long)(*npatterns), (long)(*ncolumns));

        fp = fopen(filename, "r");
        if (fp == NULL)
        {
            fprintf(stderr, "Error opening file %s\n", filename);
            fflush(stderr);
            exit(1);
        }

        for (i = 0; i < (*npatterns); i++)
        {

            dfile_nextpattern(fp, ((*d)[i]), buffer, (*ncolumns), linelength);
        }
    }

    ez_cvector(&buffer, 0);

    return;
}

int dfile_nextpattern(FILE *fp, double *vars, char *spattern, long ncolumns, long linelength)
{

    char *buffer = NULL;
    int i, j;

    ez_cvector(&buffer, linelength);

    if (feof(fp))
    {
        return 1;
    }

    fgets(buffer, (linelength - 1), fp);

    while (strchr(buffer, '#') != NULL)
    {
        fgets(buffer, (linelength - 1), fp);
    }

    strcpy(spattern, buffer);
    i = 0;
    vars[i] = atof(strtok(buffer, " \t\n"));
    i++;

    for (j = 1; j < ncolumns; j++)
    {

        vars[i] = atof(strtok(NULL, " \t\n"));
        i++;
    }

    ez_cvector(&buffer, 0);

    return 0;
}

void dfile_firstline(FILE *fp, char **titles, long ncolumns, long linelength)
{

    char *buffer = NULL;
    int i, j;

    ez_cvector(&buffer, linelength);

    if (feof(fp))
    {
        return;
    }

    fgets(buffer, (linelength - 1), fp);

    i = 0;

    strcpy(titles[i], strtok(buffer, " \t\n"));

    if (!(strlen(titles[i]) == 1 && titles[i][0] == '#'))
        i++;

    for (j = 1; j < ncolumns; j++)
    {

        strcpy(titles[i], strtok(NULL, " \t\n"));
        if (!(strlen(titles[i]) == 1 && titles[i][0] == '#'))
            i++;
    }

    ez_cvector(&buffer, 0);
    return;
}

void iParsePattern(char *s, int **d, int *n)
{

    int i, j;
    char buffer[MAX_BUFFER];

    (*n) = 0;

    /* Pass One: Determine number of tokens. */
    strcpy(buffer, s);
    if ((strchr(buffer, '#') == NULL) && (numerical(buffer)))
    {

        if ((*n) == 0)
        {

            if (strtok(buffer, " \t\n") != NULL)
            {

                (*n)++;

                while (strtok(NULL, " \t\n") != NULL)
                {
                    (*n)++;
                }
            }
        }
    }

    ez_ivector(d, ((long)(*n) - 1));

    /* Pass Two: Store tokens in d as integer values. */
    strcpy(buffer, s);

    i = 0;
    (*d)[i] = atoi(strtok(buffer, " \t\n"));
    i++;

    for (j = 1; j < (*n); j++)
    {

        (*d)[i] = atoi(strtok(NULL, " \t\n"));
        i++;
    }

    return;
}

void iParsePattern2(char *s, int **d)
{

    int j, n;
    char buffer[MAX_BUFFER];

    n = 0;

    /* Pass One: Determine number of tokens. */
    strcpy(buffer, s);
    if ((strchr(buffer, '#') == NULL) && (numerical(buffer)))
    {

        if (n == 0)
        {

            if (strtok(buffer, " \t\n") != NULL)
            {

                n++;

                while (strtok(NULL, " \t\n") != NULL)
                {
                    n++;
                }
            }
        }
    }

    ez_ivector(d, (long)(n - 1));

    /* Pass Two: Store tokens in d as integer values. */
    strcpy(buffer, s);

    (*d)[0] = atoi(strtok(buffer, " \t\n"));

    for (j = 1; j < n; j++)
    {

        (*d)[j] = atoi(strtok(NULL, " \t\n"));
    }

    return;
}

void dParsePattern(char *s, double **d, int *n)
{

    int i, j;
    char buffer[MAX_BUFFER];

    (*n) = 0;

    /* Pass One: Determine number of tokens. */
    strcpy(buffer, s);
    if ((strchr(buffer, '#') == NULL) && (numerical(buffer)))
    {

        if ((*n) == 0)
        {

            if (strtok(buffer, " \t\n") != NULL)
            {

                (*n)++;

                while (strtok(NULL, " \t\n") != NULL)
                {
                    (*n)++;
                }
            }
        }
    }

    ez_dvector(d, (long)((*n) - 1));

    /* Pass Two: Store tokens in d as double values. */
    strcpy(buffer, s);

    i = 0;
    (*d)[i] = atof(strtok(buffer, " \t\n"));
    i++;

    for (j = 1; j < (*n); j++)
    {

        (*d)[i] = atof(strtok(NULL, " \t\n"));
        i++;
    }

    return;
}

int numerical(char *s)
{
    if ((strchr(s, '0') == NULL) && (strchr(s, '1') == NULL) &&
        (strchr(s, '2') == NULL) && (strchr(s, '3') == NULL) &&
        (strchr(s, '4') == NULL) && (strchr(s, '5') == NULL) &&
        (strchr(s, '6') == NULL) && (strchr(s, '7') == NULL) &&
        (strchr(s, '8') == NULL) && (strchr(s, '9') == NULL))
        return 0;
    else
        return 1;
}

void parseFileName(char *fullFileName, char *path, char *fname, char *froot, char *fextension, char *OS)
{

    char c_dirDelimiter, s_dirDelimiter[2], *temp0 = NULL, *temp1 = NULL, *temp2 = NULL, *sPtr = NULL;

    if (strcmp(fullFileName, "") == 0)
    {
        strcpy(path, "");
        strcpy(fname, "");
        return;
    }

    temp0 = (char *)malloc(sizeof(char) * MAX_STRING);
    temp1 = (char *)malloc(sizeof(char) * MAX_STRING);
    temp2 = (char *)malloc(sizeof(char) * (strlen(fullFileName) + 1));

    sPtr = NULL;

    if (strcmp(OS, "DOS") == 0)
    {
        c_dirDelimiter = '\\';
        strcpy(s_dirDelimiter, "\\");
    }
    else
    {
        c_dirDelimiter = '/';
        strcpy(s_dirDelimiter, "/");
    }

    if (strchr(fullFileName, c_dirDelimiter) != NULL)
    {
        strcpy(temp2, fullFileName);
        strcpy(path, strtok(temp2, s_dirDelimiter));
        sPtr = strtok(NULL, s_dirDelimiter);
        while (sPtr != NULL)
        {
            strcpy(temp0, sPtr);
            strcpy(temp1, temp0);
            sPtr = strtok(NULL, s_dirDelimiter);
            if (sPtr != NULL)
            {
                strcpy(temp0, sPtr);
                strcat(path, s_dirDelimiter);
                strcat(path, temp1);
            }
            else
            {
                strcpy(fname, temp1);
            }
        }
        strcat(path, s_dirDelimiter);
    }
    else
    {
        strcpy(path, "");
        strcpy(fname, fullFileName);
    }

    if (strchr(fname, '.') != NULL)
    {
        strcpy(fextension, (strchr(fname, '.') + 1));
    }
    else
    {
        strcpy(fextension, "");
    }

    if (strchr(fname, '.') != NULL)
    {
        strncpy(froot, fname, (strchr(fname, '.') - fname));
        froot[(strchr(fname, '.') - fname)] = '\0';
    }
    else
    {
        strcpy(froot, fname);
    }

    free(temp0);
    free(temp1);
    free(temp2);

    return;
}

FILE *ftrap(char *filename, char *kind, char *how)
{
    FILE *op;
    if ((op = fopen(filename, how)) == NULL)
    {
        fprintf(stderr, "Could not open %s file '%s'\nin mode '%s'\n",
                kind, filename, how);
        fprintf(stderr, "Press <return> to continue.\n");
        getchar();
        exit(1);
    }
    return op;
}

FILE *ftrap2(const char *filename, const char *mode)
{
    FILE *fp;
    if ((fp = fopen(filename, mode)) == NULL)
    {
        fprintf(stderr, "Could not open file '%s'\nin mode '%s'\n", filename, mode);
        fprintf(stderr, "Press <return> to continue.\n");
        fflush(stderr);
        getchar();
        exit(1);
    }
    return fp;
}
