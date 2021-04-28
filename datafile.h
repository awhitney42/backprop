/* ************************************************** */
/* datafile.h                                         */
/*                                                    */
/* Adam P. Whitney                                    */
/* ************************************************** */

#ifndef MAX_STRING
#define MAX_STRING 2048
#endif
#ifndef MAX_BUFFER
#define MAX_BUFFER 4096
#endif


typedef struct {
	int noutputs;
	int npatterns;
	int ncolumns;
	double **inputs;
	double **outputs;
} dfileInfo;


/* dfile_read                                                        */
/*                                                                   */
/* Reads filename placing the (double) values into d.                */
/*                                                                   */
/* PRE:      d == NULL                                               */
/* MODIFIES: d, npatterns, ncolumns                                  */
/* ON ERROR: Writes to stderr then exit(1);                          */

void dfile_read(char *filename, double ***d, long *npatterns, long *ncolumns, int printStatus);



/* dfile_info                                                        */
/*                                                                   */
/* MODIFIES: npatterns, ncolumns, linelength                         */
/* ON ERROR: Writes to stderr then exit(1);                          */

void dfile_info(char *filename, long *npatterns, long *ncolumns, long *linelength);



/* dfile_nextpattern                                                 */
/*                                                                   */
/* Reads the next pattern in fp.                                     */
/* spattern contains the pattern as a character string.              */
/* vars contains a (double) vector of the pattern.                   */
/*                                                                   */
/* PRE:      vars is a double vector with at least ncolumns elements */
/* MODIFIES: vars, spattern                                          */
/* RETURN VAL: 1 if fp is at the end of file, otherwise 0            */

int dfile_nextpattern(FILE *fp, double *vars, char *spattern, long ncolumns, long linelength);



void dfile_firstline(FILE *fp, char **titles, long ncolumns, long linelength);




/* iParsePattern                                                     */
/*                                                                   */
/* Parses the character string s, placing the space/tab character    */
/* delimited integer values in the integer vector d.                 */
/* Stores the integer value count (size of d) in n.                  */
/*                                                                   */
/* PRE:      d == NULL                                               */
/* MODIFIES: d, n                                                    */

void iParsePattern(char *s, int **d, int *n);



/* iParsePattern2                                                    */
/*                                                                   */
/* Parses the character string s, placing the space/tab character    */
/* delimited integer values in the integer vector d.                 */
/* Stores the integer value count in d[0].                           */
/*                                                                   */
/* PRE:      d == NULL                                               */
/* MODIFIES: d                                                       */

void iParsePattern2(char *s, int **d);



/* dParsePattern                                                     */
/*                                                                   */
/* Parses the character string s, placing the space/tab character    */
/* delimited double values in the double vector d.                   */
/* Stores the double value count (size of d) in n.                   */
/*                                                                   */
/* PRE:      d == NULL                                               */
/* MODIFIES: d, n                                                    */

void dParsePattern(char *s, double **d, int *n);



/* numerical                                                         */
/*                                                                   */
/* Returns 1 if the character string s contains at least one numeric */
/* character, 0 otherwise.                                           */

int numerical(char *s);



/* parseFileName                                                     */
/*                                                                     */
/* Parses fullFileName to determine various useful file information.     */
/*                                                                         */
/* PRE:      OS == "DOS" to use directory delimiter character '\'            */
/*           Otherwise, directory delimiter character is '/'                   */
/* MODIFIES: path, fname, froot, fextension                                      */
/*                                                                                 */
/* Example: After calling the function like so, the following is true...           */
/*                                                                                 */
/* parseFileName("d:\data\working\abc.txt", path, fname, froot, fextension, "DOS");*/
/*                                                                                 */
/* path  == "d:\data\working"                                                      */
/* fname == "abc.txt"                                                              */
/* froot == "abc"                                                                  */
/* fextension == "txt"                                                             */

void parseFileName(char *fullFileName, char *path, char *fname, char *froot, char *fextension, char *OS);


/* ftrap                                                             */
/*                                                                   */
/* File opening trap, you give it the file name, the kind of file it */
/* is (for the error message, just in case it can't be opened) and   */
/* the mode (the same as with fopen)                                 */

FILE *ftrap( char *filename, char *kind, char *how );

FILE *ftrap2( const char *filename, const char *mode );

