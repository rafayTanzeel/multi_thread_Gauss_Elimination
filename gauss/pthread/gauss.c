/*
 * Original author:  UNKNOWN
 *
 * Modified:         Kai Shen (January 2010)
 */
#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>
#include <sched.h>
#include <pthread.h>

//#define DEBUG

#define SWAP(a,b)       {double tmp; tmp = a; a = b; b = tmp;}

/* Solve the equation:
 *   matrix * X = R
 */
struct timeval start, finish;
int task_num = 1;
int nsize = 0;
int nthreads = 1;
double **matrix, *X, *R;

/* Pre-set solution. */
double *X__;
double *PivotRow;
double *BufferRow;
int emptier=1;
/* Initialize the matirx. */

pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutPass = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t mutEmpty = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t cond = PTHREAD_COND_INITIALIZER;

void barrier (int expect);
void getPivot(int nsize, int currow);
double getPivotRowElement(int index);
void computeGauss(int nsize, int task_id);
void* work_thread (void *lp);
void solveGauss(int nsize);
int initMatrix(const char *fname);
void initRHS(int nsize);
void initResult(int nsize);
void emptyMatrix();
int evenCore(int task_id, int no_cores);
int oddCore(int task_id, int no_cores);
int amoebaCore(int task_id, int no_cores);

int main(int argc, char *argv[])
{
    int i;
    double error;
    pthread_attr_t attr;
    pthread_t *tid;
    int *id;

    if (argc > 3) {
    fprintf(stderr, "usage: %s <matrixfile>\n", argv[1]);
    exit(-1);
    }

    nsize = initMatrix(argv[1]);

    if(argc==3)
        task_num = atoi(argv[2]);

    initRHS(nsize);
    initResult(nsize);

    id = (int *) malloc (sizeof (int) * task_num);
    tid = (pthread_t *) malloc (sizeof (pthread_t) * task_num);

    PivotRow = (double *) malloc(sizeof(double) * nsize);
    BufferRow = (double *) malloc(sizeof(double) * nsize);

    for (i = 0; i < nsize; i++) {
        PivotRow[i]=0.0;
        BufferRow[i]=0.0;
     }

    if (!id || !tid){
        fprintf(stderr, "out of shared memory");
        exit(-1);
    }

    pthread_attr_init (&attr);
    pthread_attr_setscope (&attr, PTHREAD_SCOPE_SYSTEM);

    for (i = 1; i < task_num; i++) {
        id[i] = i;
        pthread_create (&tid[i], &attr, work_thread, &id[i]);
    }
    id[0]=0;
    work_thread(&id[0]);

    for (i = 1; i < task_num; i++)
        pthread_join (tid[i], NULL);

    emptyMatrix();

    solveGauss(nsize);

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


double getPivotRowElement(int index){
	if(PivotRow[index]==0.0){
		pthread_mutex_lock (&mutPass);
		if(PivotRow[index]==0.0){
			getPivot(nsize,index);
			PivotRow[index]=matrix[index][index];
		}
		pthread_mutex_unlock (&mutPass);
	}
	return PivotRow[index];
}

int evenCore(int task_id, int no_cores){
	if(task_id%2!=0){
		return (task_id+1)%no_cores;
	}
	return task_id%no_cores;
}


int oddCore(int task_id, int no_cores){
	if(task_id%2==0){
		return (task_id+1)%no_cores;
	}
	return task_id%no_cores;
}

int amoebaCore(int task_id, int no_cores){
    int even[12] = {0,2,4,6,8,10,12,14,16,18,20,22};
    int odd[12] = {1,3,5,7,9,11,13,15,17,19,21,23};
    int tid = task_id%no_cores;
    if (tid < 12) {
        return even[tid];
    }
    else if (tid < 24) {
        return odd[tid-12];
    }
    else {
        return -1;
    }
}


void* work_thread (void *lp)
{

    int task_id = *((int *) lp);

	int no_cores = sysconf(_SC_NPROCESSORS_ONLN);

	cpu_set_t cpuset;

	int cpu=task_id%no_cores;

	printf("CPU=%d\n",cpu);

	CPU_ZERO(&cpuset);
	CPU_SET( cpu , &cpuset);

	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);

    barrier (task_num);

    if(task_id==0)
        gettimeofday (&start, NULL);

    computeGauss(nsize, task_id);

//    barrier (task_num);
    gettimeofday (&finish, NULL);

    return NULL;

}


int initMatrix(const char *fname)
{

    FILE *file;
    int l1, l2, l3;
    double d;
    int nsize;
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

void initRHS(int nsize)
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

void initResult(int nsize)
{
    int i;

    X = (double*)malloc(nsize * sizeof(double));
    assert(X != NULL);
    for (i = 0; i < nsize; i++) {
    X[i] = 0.0;
    }
}

/* Get the pivot - the element on column with largest absolute value. */

void getPivot(int nsize, int currow)
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

void computeGauss(int nsize, int task_id)
{
    int i, j, k;
    double pivotval;

    for (i = 0; i < nsize; i++) { //i is rows
    /* Scale the main row. */

        pivotval = getPivotRowElement(i);

        if (pivotval != 1.0) {
            matrix[i][i] = 1.0;
            for (j = task_id + i + 1; j < nsize; j+=task_num) { //j is column
                matrix[i][j] /= pivotval;
                BufferRow[j]=matrix[j][i];
            }
            if(task_id==0){
                R[i] /= pivotval;
            }
        }

    /* Factorize the rest of the matrix. */
        for (j = i + 1; j < nsize; j++) {
            for (k = task_id + i + 1; k < nsize; k+=task_num) {
                matrix[j][k] -= matrix[j][i] * matrix[i][k];
            }
            if (task_id==0) {
                R[j] -= matrix[j][i] * R[i];
            }
        }
        barrier (task_num);
    }
}

void emptyMatrix() {
    int i, j;
    for (i = 0; i < nsize; i++) {
        for (j = i + 1; j < nsize; j++) {
            matrix[j][i] = 0.0;
        }
    }
}

/* Solve the equation. */

void solveGauss(int nsize)
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


