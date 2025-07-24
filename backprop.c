/* Adam P. Whitney */
/* 13 May 2000     */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "ez_alloc.h"
#include "datafile.h"

#define NUM_LAYERS 3

#define TRAINING_EPOCHS 10000

extern float rand0(long *idum);
void readTrainFile(char *trainFilename, double **x, long npatterns, long ncolumns);
void normalizeData(double **x, int npatterns, int ncolumns, int activation_function);
double neuronResponse(double v, int activation_function);
double neuronResponsePrime(double v, int activation_function);
void inducedLocalField(int layer, int n);
void specifyArchitecture(int numInputs, int numHiddenNodes, int numOutputs);
void initializeWeights();
void trainNetwork(double **x, int npatterns, int activation_function);
void recallNetwork(double **x, int npatterns, int activation_function);
void forwardPass(int pat_num, double **x, int activation_function);
void backwardPass(int pat_num, double **x, int activation_function);
void writeWeights(char *filename);
void readWeights(char *filename);

/* Globals */

typedef struct
{
   double bias;
   double v_j;
   double y_j;
   double *w_ji;
   double d_j;
   double *delta_w_ji;

} node;

node **neuron = NULL;
int neuronsPerLayer[NUM_LAYERS];

/* Main Function */

int main(int argc, char **argv)
{
   long npatterns, ncolumns;
   int numInputs, numHiddenNodes, numOutputs;
   int i;
   /* int j; */
   long linelength = 0;
   char dataFilename[MAX_STRING], wtsFilename[MAX_STRING];
   double **x = NULL;
   int recall_flag = 0;
   int activation_function = 0;

   if (argc < 3)
   {
      fprintf(stderr, "usage: %s <data file> <wts file> <0=train/1=recall> [<0=logistic/1=hyperbolic tangent>] [<number of outputs>]\n", argv[0]);
      exit(1);
   }

   /* Process Arguments */

   strcpy(dataFilename, argv[1]);
   strcpy(wtsFilename, argv[2]);
   recall_flag = atoi(argv[3]);
   if (argc >= 5)
   {
      activation_function = atoi(argv[4]);
   }
   if (argc >= 6)
   {
      numOutputs = atoi(argv[5]);
   }
   else
   {
      numOutputs = 1;
   }

   /* Process Data */

   dfile_info(dataFilename, &npatterns, &ncolumns, &linelength);

   x = (double **)malloc(npatterns * sizeof(double *));
   x[0] = (double *)malloc(npatterns * ncolumns * sizeof(double));
   for (i = 1; i < npatterns; i++)
      x[i] = x[0] + i * ncolumns;

   readTrainFile(dataFilename, x, npatterns, ncolumns);

   //   normalizeData(x, npatterns, ncolumns, activation_function);

   /*
   for (i = 0; i < npatterns; i++) {
      for (j = 0; j < ncolumns; j++) {
    if (j == 0) {
       printf("%f", x[i][j]);
    } else {
       printf("\t%f", x[i][j]);
         }
      }
      printf("\n");
   }
   printf("\n");

   fflush(stdout);
   */

   /****************/

   numInputs = ncolumns - numOutputs;
   numHiddenNodes = 2 * numInputs;

   printf("Data Cols: %ld, Inputs: %d, Outputs: %d, Hidden Nodes: %d\n", ncolumns, numInputs, numOutputs, numHiddenNodes);

   specifyArchitecture(numInputs, numHiddenNodes, numOutputs);

   printf("Neural Net Created\n");

   initializeWeights();

   printf("Weights Initialized\n");

   readWeights(wtsFilename);

   if (recall_flag)
   {

      /* Recall the network. */

      free(x[0]);
      free(x);

      dfile_info(dataFilename, &npatterns, &ncolumns, &linelength);

      x = (double **)malloc(npatterns * sizeof(double *));
      x[0] = (double *)malloc(npatterns * ncolumns * sizeof(double));
      for (i = 1; i < npatterns; i++)
         x[i] = x[0] + i * ncolumns;

      readTrainFile(dataFilename, x, npatterns, ncolumns);

      //      normalizeData(x, npatterns, ncolumns, activation_function);
   }
   else
   {

      /* Train the network. */

      if (activation_function == 0)
      {
         printf("Activation Function: Logistic Function\n");
      }
      else
      {
         printf("Activation Function: Hyperbolic Tangent Function\n");
      }

      trainNetwork(x, npatterns, activation_function);

      writeWeights(wtsFilename);
   }

   recallNetwork(x, npatterns, activation_function);

   free(x[0]);
   free(x);

   /* getchar(); */

   exit(0);
}

void writeWeights(char *filename)
{

   FILE *fp;
   int i, j, n;

   fp = fopen(filename, "w");

   for (j = 1; j < NUM_LAYERS; j++)
   {
      for (n = 0; n < neuronsPerLayer[j]; n++)
      {

         // Here's the bias
         fprintf(fp, "%f ", neuron[j][n].bias);

         for (i = 0; i < neuronsPerLayer[j - 1]; i++)
         {

            fprintf(fp, "%f ", neuron[j][n].w_ji[i]);
         }
      }
      fprintf(fp, "\n");
   }

   fclose(fp);
}

void readWeights(char *filename)
{

   FILE *fp = NULL;
   int i, j, n;
   char buffer[MAX_STRING];

   double *wts = NULL;
   int wcount = 0, w_ix = 0;

   fp = fopen(filename, "r");

   if (fp != NULL)
   {

      for (j = 1; j < NUM_LAYERS; j++)
      {

         memset(buffer, '\0', MAX_STRING);

         fgets(buffer, MAX_STRING - 1, fp);

         dParsePattern(buffer, &wts, &wcount);

         w_ix = 0;
         for (n = 0; n < neuronsPerLayer[j]; n++)
         {

            // Here's the bias
            neuron[j][n].bias = wts[w_ix++];

            for (i = 0; i < neuronsPerLayer[j - 1]; i++)
            {

               neuron[j][n].w_ji[i] = wts[w_ix++];
            }
         }
      }

      fclose(fp);
   }
}

void forwardPass(int pat_num, double **x, int activation_function)
{

   int layer, n;

   //   printf ("%d ******\n", pat_num);

   for (layer = 0; layer < NUM_LAYERS; layer++)
   {

      for (n = 0; n < neuronsPerLayer[layer]; n++)
      {

         if (layer == 0)
         {
            neuron[layer][n].y_j = x[pat_num][n];
         }
         else
         {

            /* Calculate v_j */
            inducedLocalField(layer, n);
            /* Calculate y_j */
            neuron[layer][n].y_j = neuronResponse(neuron[layer][n].v_j, activation_function);
         }
      }
   }
}

void backwardPass(int pat_num, double **x, int activation_function)
{

   int layer, n, k, i;
   double e_j = 0.0, step_size = 0.01;

   for (layer = NUM_LAYERS - 1; layer > 0; layer--)
   {

      for (n = 0; n < neuronsPerLayer[layer]; n++)
      {

         if (layer == (NUM_LAYERS - 1))
         {

            /* Calculate d_j */

            e_j = x[pat_num][neuronsPerLayer[0]+n] - neuron[layer][n].y_j;
            // printf("%d Error: %f\n", pat_num, e_j);
            neuron[layer][n].d_j = e_j * neuronResponsePrime(neuron[layer][n].v_j, activation_function);
            // printf("%d v):       %f\n", pat_num, neuron[j][n].v_j);
         }
         else
         {

            /* Calculate d_j */
            neuron[layer][n].d_j = 0.0;
            for (k = 0; k < neuronsPerLayer[layer + 1]; k++)
            {
               neuron[layer][n].d_j += neuronResponsePrime(neuron[layer][n].v_j, activation_function) *
                                   (neuron[layer + 1][k].d_j * neuron[layer + 1][k].w_ji[n]);
            }
         }

         /* Calculate delta_w_ji and adjust w_ji */

         for (i = 0; i < neuronsPerLayer[layer - 1]; i++)
         {

            neuron[layer][n].delta_w_ji[i] = neuron[layer][n].d_j * neuron[layer - 1][i].y_j;
            neuron[layer][n].w_ji[i] += step_size * neuron[layer][n].delta_w_ji[i];

            // printf("%f ", neuron[layer][n].w_ji[i]);
         }
         // printf ("\n");

      } /* for (n...) */

      // printf ("\n");

   } /* for (layer...) */
}

void recallNetwork(double **x, int npatterns, int activation_function)
{

   int i, j;

   double yy, yd, a1, o1;
   double *yi, *yip, *yisq, *yipsq, *yiyip;
   double sigyisq, sigyipsq, sigyiyip;
   double *rms;
   double rtemp, rsquared;

   yi = (double *)malloc(sizeof(double) * neuronsPerLayer[NUM_LAYERS - 1]);
   yip = (double *)malloc(sizeof(double) * neuronsPerLayer[NUM_LAYERS - 1]);
   yisq = (double *)malloc(sizeof(double) * neuronsPerLayer[NUM_LAYERS - 1]);
   yipsq = (double *)malloc(sizeof(double) * neuronsPerLayer[NUM_LAYERS - 1]);
   yiyip = (double *)malloc(sizeof(double) * neuronsPerLayer[NUM_LAYERS - 1]);
   rms = (double *)malloc(sizeof(double) * neuronsPerLayer[NUM_LAYERS - 1]);

   for (j = 0; j < neuronsPerLayer[NUM_LAYERS - 1]; j++)
   {
      yi[j] = 0.0;
      yip[j] = 0.0;
      yisq[j] = 0.0;
      yipsq[j] = 0.0;
      yiyip[j] = 0.0;
      rms[j] = 0.0;
   }

   for (j = 0; j < neuronsPerLayer[0]; j++)
   {
      if (j == 0)
      {
         printf("#input%d", j + 1);
      }
      else
      {
         printf("\t\tinput%d", j + 1);
      }
   }
   for (j = 0; j < neuronsPerLayer[NUM_LAYERS - 1]; j++)
   {
      printf("\t\tdesired%d\tmodel%d", j + 1, j + 1);
   }
   printf("\n");

   for (i = 0; i < npatterns; i++)
   {
      forwardPass(i, x, activation_function);
      for (j = 0; j < neuronsPerLayer[0]; j++)
      {
         if (j == 0)
         {
            printf("%f", x[i][j]);
         }
         else
         {
            printf("\t%f", x[i][j]);
         }
      }

      for (j = 0; j < neuronsPerLayer[NUM_LAYERS - 1]; j++)
      {

         printf("\t%f\t%f", x[i][neuronsPerLayer[0]+j], neuron[NUM_LAYERS - 1][j].y_j);

         yy = neuron[NUM_LAYERS - 1][j].y_j; /* model output */
         yd = x[i][neuronsPerLayer[0]+j];      /* desired output */

         rtemp = yy - yd;
         rms[j] += rtemp * rtemp;

         a1 = yy;
         o1 = yd;
         yi[j] += a1;
         yip[j] += o1;
         yisq[j] += a1 * a1;
         yipsq[j] += o1 * o1;
         yiyip[j] += a1 * o1;
      }
      printf("\n");
   } /* for i = 0; i < npatters; i++ */

   for (j = 0; j < neuronsPerLayer[NUM_LAYERS - 1]; j++)
   {
      /*******************/
      /** calculate RMS **/
      /*******************/
      rms[j] = sqrt(rms[j] / (double)(npatterns));
      printf("#\n#\n# unnormalized rms for output %d: %lf\n", j+1, rms[j]);

      /*******************/
      /** calculate R^2 **/
      /*******************/
      sigyisq = yisq[j] - (yi[j] * yi[j]) / (double)npatterns;
      sigyipsq = yipsq[j] - (yip[j] * yip[j]) / (double)npatterns;
      sigyiyip = yiyip[j] - (yi[j] * yip[j]) / (double)npatterns;
      if ((sigyisq * sigyipsq) > 1e-5)
         rsquared = (sigyiyip / (sqrt(sigyisq * sigyipsq))) *
                  (sigyiyip / (sqrt(sigyisq * sigyipsq)));
      else
         rsquared = 0.0;
      printf("# R^2 for output %d:              %lf\n", j+1, rsquared);
   }

   free(yi);
   free(yip);
   free(yisq);
   free(yipsq);
   free(yiyip);
   free(rms);

}

void trainNetwork(double **x, int npatterns, int activation_function)
{

   int i, epochCount = 0;

   while (epochCount < TRAINING_EPOCHS)
   {

      for (i = 0; i < npatterns; i++)
      {
         forwardPass(i, x, activation_function);
         backwardPass(i, x, activation_function);
      }
      // printf ("\n");

      epochCount++;
      // if (epochCount % 100 == 0) {
      if (epochCount == TRAINING_EPOCHS)
      {
         printf("Finished epoch %d of %d.\n", epochCount, TRAINING_EPOCHS);
         fflush(stdout);
      }
   }

   printf("Finished training.\n");
}

void initializeWeights()
{

   int i, j, k;
   long q = 5;

   /*   srand( time(NULL) + clock() ); */

   for (i = 0; i < NUM_LAYERS; i++)
   {
      for (j = 0; j < neuronsPerLayer[i]; j++)
      {
         if (i == 0)
         {

            neuron[i][j].w_ji = NULL;
            neuron[i][j].delta_w_ji = NULL;
         }
         else
         {

            neuron[i][j].w_ji = (double *)malloc(sizeof(double) * neuronsPerLayer[i - 1]);
            neuron[i][j].delta_w_ji = (double *)malloc(sizeof(double) * neuronsPerLayer[i - 1]);
            for (k = 0; k < neuronsPerLayer[i - 1]; k++)
            {
               neuron[i][j].w_ji[k] = 2 * rand0(&q) - 1;
               neuron[i][j].delta_w_ji[k] = 0;
            }
            neuron[i][j].bias = 2 * rand0(&q) - 1;
         }
      }
   }
}

void specifyArchitecture(int numInputs, int numHiddenNodes, int numOutputs)
{

   int i;

   neuron = (node **)malloc(sizeof(node *) * NUM_LAYERS);

   for (i = 0; i < NUM_LAYERS; i++)
   {
      if (i == 0)
      {
         neuron[i] = (node *)malloc(sizeof(node) * numInputs);
         neuronsPerLayer[i] = numInputs;
      }
      else if (i == 1)
      {
         neuron[i] = (node *)malloc(sizeof(node) * numHiddenNodes);
         neuronsPerLayer[i] = numHiddenNodes;
      }
      else
      {
         neuron[i] = (node *)malloc(sizeof(node) * numOutputs);
         neuronsPerLayer[i] = numOutputs;
      }
   }
}

void inducedLocalField(int layer, int n)
{

   int m;
   neuron[layer][n].v_j = neuron[layer][n].bias;
   for (m = 0; m < neuronsPerLayer[layer - 1]; m++)
   {
      neuron[layer][n].v_j += neuron[layer][n].w_ji[m] * neuron[layer - 1][m].y_j;
   }
}

double neuronResponse(double v, int activation_function)
{

   double a, b;

   if (activation_function == 0)
   {
      a = 1.1;
      if (v == 0)
      {
         v = 0.000000001;
      }
      return (1 / (1 + exp(-1 * a * v)));
   }
   else
   {
      a = b = 1.;
      //      a = 1.7159, b = 0.666667;
      if (v == 0)
      {
         v = 0.000000001;
      }
      return (a * tanh(b * v));
   }
}

double neuronResponsePrime(double v, int activation_function)
{

   double a, b;

   if (activation_function == 0)
   {
      a = 1.1;
      if (v == 0)
      {
         v = 0.000000001;
      }
      return ((a * exp(-1 * a * v)) / ((1 + exp(-1 * a * v)) * (1 + exp(-1 * a * v))));
   }
   else
   {
      a = b = 1.;
      //      a = 1.7159, b = 0.666667;
      if (v == 0)
      {
         v = 0.000000001;
      }
      return (a * b * (1 - (tanh(b * v) * tanh(b * v))));
   }
}

void normalizeData(double **x, int npatterns, int ncolumns, int activation_function)
{

   int i, j;
   double *maxValues, *minValues;

   maxValues = (double *)malloc(sizeof(double) * ncolumns);
   minValues = (double *)malloc(sizeof(double) * ncolumns);

   for (i = 0; i < ncolumns; i++)
   {
      maxValues[i] = 0;
      minValues[i] = x[0][i];
   }

   /* Find the maximum and minimum values for each of the columns. */

   for (i = 0; i < npatterns; i++)
   {
      for (j = 0; j < ncolumns; j++)
      {
         if (x[i][j] > maxValues[j])
         {
            maxValues[j] = x[i][j];
         }
         if (x[i][j] < minValues[j])
         {
            minValues[j] = x[i][j];
         }
      }
   }

   for (i = 0; i < ncolumns; i++)
   {
      if (maxValues[i] == 0 && minValues[i] == 0)
      {
         fprintf(stderr, "maxValues[%d] == 0 and minValues[%d] == 0!\n", i, i);
         fprintf(stderr, "Hmmmm.  Aborting execution...\n");
         exit(1);
      }
   }

   /* Now we normalize the columns. */

   for (i = 0; i < npatterns; i++)
   {
      for (j = 0; j < ncolumns; j++)
      {
         if (j == ncolumns - 1)
         {
            if (activation_function == 1)
               x[i][j] = (x[i][j] - minValues[j]) / (maxValues[j] - minValues[j]) * 2 - 1;
            else
               /* Normalize from 0 to +1 */
               x[i][j] = (x[i][j] - minValues[j]) / (maxValues[j] - minValues[j]);
         }
         else
         {
            /* Normalize from -1 to +1 */
            x[i][j] = (x[i][j] - minValues[j]) / (maxValues[j] - minValues[j]) * 2 - 1;
         }
      }
   }

   free(maxValues);
   free(minValues);
}

void readTrainFile(char *trainFilename, double **x, long npatterns, long ncolumns)
{

   FILE *trainFile;
   int i;
   char spattern[MAX_STRING];

   trainFile = fopen(trainFilename, "r");
   if (trainFile == NULL)
   {
      fprintf(stderr, "Error opening the file %s\n", trainFilename);
      exit(1);
   }

   for (i = 0; i < npatterns; i++)
   {
      dfile_nextpattern(trainFile, x[i], spattern, ncolumns, MAX_BUFFER);
   }

   fclose(trainFile);
}
