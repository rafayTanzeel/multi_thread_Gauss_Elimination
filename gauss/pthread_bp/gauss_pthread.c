/* 
 * Original author:  UNKNOWN
 *
 * Modified:         Kai Shen (January 2010)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>

#ifndef _REENTRANT
#define _REENTRANT      /* basic 3-lines for threads */
#endif
#include <pthread.h>

// #define DEBUG

#define SWAP(a,b)       {double tmp; tmp = a; a = b; b = tmp;}

/* Solve the equation:
 *   matrix * X = R
 */

double **matrix, *X, *R;

/* Pre-set solution. */

double *X__;

int task_num;
int nsize;
struct timeval start, finish;

pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

void barrier (int expect)
{
    static int arrived = 0;

    pthread_mutex_lock (&mut);  //lock

    arrived++;
    if (arrived < expect)
        pthread_cond_wait (&cond, &mut);
    else {
        arrived = 0;        // reset the barrier before broadcast is important
        pthread_cond_broadcast (&cond);
    }

    pthread_mutex_unlock (&mut);    //unlock
}


/* Initialize the matirx. */

int initMatrix(const char *fname)
{
    FILE *file;
    int l1, l2, l3;
    double d;
    
    int i, j;
    double *tmp;
    char buffer[1024];

    if ((file = fopen(fname, "r")) == NULL) {
	fprintf(stderr, "The matrix file open error\n");
        exit(-1);
    }
    
    /* Parse the first line to get the matrix size. */
    fgets(buffer, 1024, file);
    sscanf(buffer, "%d %d %d", &l1, &l2, &l3);
    nsize = l1;
#ifdef DEBUG
    fprintf(stdout, "matrix size is %d\n", nsize);
#endif

    /* Initialize the space and set all elements to zero. */
    matrix = (double**)malloc(nsize*sizeof(double*));
    assert(matrix != NULL);
    tmp = (double*)malloc(nsize*nsize*sizeof(double));
    assert(tmp != NULL);    
    for (i = 0; i < nsize; i++) {
        matrix[i] = tmp;
        tmp = tmp + nsize;
    }
    for (i = 0; i < nsize; i++) {
        for (j = 0; j < nsize; j++) {
            matrix[i][j] = 0.0;
        }
    }

    /* Parse the rest of the input file to fill the matrix. */
    for (;;) {
	fgets(buffer, 1024, file);
	sscanf(buffer, "%d %d %lf", &l1, &l2, &d);
	if (l1 == 0) break;

	matrix[l1-1][l2-1] = d;
#ifdef DEBUG
	fprintf(stdout, "row %d column %d of matrix is %e\n", l1-1, l2-1, matrix[l1-1][l2-1]);
#endif
    }

    fclose(file);
    return nsize;
}

/* Initialize the right-hand-side following the pre-set solution. */

void initRHS()
{
    int i, j;

    X__ = (double*)malloc(nsize * sizeof(double));
    assert(X__ != NULL);
    for (i = 0; i < nsize; i++) {
	X__[i] = i+1;
    }

    R = (double*)malloc(nsize * sizeof(double));
    assert(R != NULL);
    for (i = 0; i < nsize; i++) {
	R[i] = 0.0;
	for (j = 0; j < nsize; j++) {
	    R[i] += matrix[i][j] * X__[j];
	}
    }
}

/* Initialize the results. */

void initResult()
{
    int i;

    X = (double*)malloc(nsize * sizeof(double));
    assert(X != NULL);
    for (i = 0; i < nsize; i++) {
	X[i] = 0.0;
    }
}

/* Get the pivot - the element on column with largest absolute value. */

void getPivot(int currow)
{
    int i, pivotrow;

    pivotrow = currow;
    for (i = currow+1; i < nsize; i++) {
        if (fabs(matrix[i][currow]) > fabs(matrix[pivotrow][currow])) {
            pivotrow = i;
        }
    }

    if (fabs(matrix[pivotrow][currow]) == 0.0) {
        fprintf(stderr, "The matrix is singular\n");
        exit(-1);
    }
    
    if (pivotrow != currow) {
#ifdef DEBUG
	fprintf(stdout, "pivot row at step %5d is %5d\n", currow, pivotrow);
#endif
        for (i = currow; i < nsize; i++) {
            SWAP(matrix[pivotrow][i],matrix[currow][i]);
        }
        SWAP(R[pivotrow],R[currow]);
    }
}

/* For all the rows, get the pivot and eliminate all rows and columns
 * for that particular pivot row. */

void computeGauss(int i, int begin, int end, int task_id)
{
    int j, k;
    double pivotval;

    // printf(" id: %i\n", task_id);

    if (task_id==0)
        getPivot(i);

    pivotval = matrix[i][i];
    barrier (task_num);

    /* Scale the main row. */
	if (pivotval != 1.0) {
	    if (task_id==0)
            matrix[i][i] = 1.0;
        barrier (task_num);
	    
        for (j = begin; j < end; j++) {
            if (j == i) continue;
            //printf("!!! id %i j: %i begin: %i end: %i mat: %f piv: %f \n", task_id, j, begin, end, matrix[i][j], pivotval);

            matrix[i][j] /= pivotval;
	    }

        if (task_id==0)
            R[i] /= pivotval;
	}

	/* Factorize the rest of the matrix. */
    for (j = i + 1; j < nsize; j++) {
        pivotval = matrix[j][i];
        barrier (task_num);

        if (task_id==0)
            matrix[j][i] = 0.0;
        barrier (task_num);

        for (k = begin; k < end; k++) {
            if (k == i) continue;
            //printf("!!! jk %f i: %i k: %i piv: %f ik: %f \n", matrix[j][k], i, k, pivotval, matrix[i][k]);

            matrix[j][k] -= pivotval * matrix[i][k];
        }

        if (task_id==0)
            R[j] -= pivotval * R[i];
    }

    barrier (task_num);

    // int x,y;
    // for (x=0;x<nsize;x++)
    // {
    //     for (y=0;y<nsize;y++)
    //     {
    //         printf("%f\n", matrix[x][y]);
    //     }
    // }

    // barrier (task_num);
}

/* Solve the equation. */

void solveGauss()
{
    int i, j;

    X[nsize-1] = R[nsize-1];
    for (i = nsize - 2; i >= 0; i --) {
        X[i] = R[i];
        for (j = nsize - 1; j > i; j--) {
            X[i] -= matrix[i][j] * X[j];
        }
    }

#ifdef DEBUG
    fprintf(stdout, "X = [");
    for (i = 0; i < nsize; i++) {
        fprintf(stdout, "%.6f ", X[i]);
    }
    fprintf(stdout, "];\n");
#endif
}




void *work_thread (void *lp)
{
    int task_id = *((int *) lp);
    int begin, end;

    int i;
    for (i = 0; i < nsize; i++) { //nsize

        /*get the divided task*/
        begin = ((nsize - i) * task_id) / task_num + i;
        end = ((nsize - i) * (task_id + 1)) / task_num + i;

        //fprintf (stderr, "thread %d: begin %d, end %d\n", task_id, begin, end);

        barrier (task_num);

        if (task_id==0 && i==0) {
            gettimeofday (&start, NULL);       
        }

        computeGauss(i, begin, end, task_id);

    }
    gettimeofday (&finish, NULL);
}

void
errexit (const char *err_str)
{
    fprintf (stderr, "%s", err_str);
    exit (1);
}


int main(int argc, char *argv[])
{
    int i;

    double error;

    pthread_attr_t attr;
    pthread_t *tid;
    int *id;
    
    if (argc != 3) {
	fprintf(stderr, "usage: %s <matrixfile> <numberofprocessors>\n", argv[0]);
	exit(-1);
    }

    task_num = atoi(argv[2]);
    nsize = initMatrix(argv[1]);
    initRHS();
    initResult();


    // create threads
    id = (int *) malloc (sizeof (int) * task_num);
    tid = (pthread_t *) malloc (sizeof (pthread_t) * task_num);
    if (!id || !tid)
        errexit ("out of shared memory");
    pthread_attr_init (&attr);
    pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);

    for (i = 1; i < task_num; i++) {
        id[i] = i;
        pthread_create (&tid[i], &attr, work_thread, &id[i]);
    }

    id[0]=0;
    work_thread(&id[0]);

    // wait for all threads to finish
    for (i = 1; i < task_num; i++)
    pthread_join (tid[i], NULL);


    solveGauss();
    
    fprintf(stdout, "Time:  %f seconds\n", (finish.tv_sec - start.tv_sec) + (finish.tv_usec - start.tv_usec)*0.000001);

    error = 0.0;
    for (i = 0; i < nsize; i++) {
	double error__ = (X__[i]==0.0) ? 1.0 : fabs((X[i]-X__[i])/X__[i]);
	if (error < error__) {
	    error = error__;
	}
    }
    fprintf(stdout, "Error: %e\n", error);

    return 0;
}
