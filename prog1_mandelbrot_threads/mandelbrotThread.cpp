#include <stdio.h>
#include <pthread.h>
#include <algorithm>

#include "CycleTimer.h"

typedef struct {
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int* output;
    int threadId;
    int numThreads;
} WorkerArgs;


extern void mandelbrotSerial(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int numRows,
    int maxIterations,
    int output[]);


//
// workerThreadStart --
//
// Thread entrypoint.
void* workerThreadStart(void* threadArgs) {

    WorkerArgs* args = static_cast<WorkerArgs*>(threadArgs);
    double minThread = 1e30;
    int i;
    // TODO FOR CS348V STUDENTS: Implement worker thread here. Each
    // thread should make a call to mandelbrotSerial to compute a part
    // of the output image.  For example, in a program that uses two
    // threads, thread 0 could compute the top half of the image and
    // thread 1 could compute the bottom half.
    double startTime = CycleTimer::currentSeconds();
    /* 
    mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, 
            args->width, args->height, 
            args->threadId*args->height/args->numThreads, 
            args->height/args->numThreads, 
            args->maxIterations, args->output);
    */
    for(i=0; i<int(args->height/args->numThreads); ++i)
    {
        int pos = args->numThreads*i + args->threadId;
        mandelbrotSerial(args->x0, args->y0, args->x1, args->y1, 
            args->width, args->height, 
            pos, 
            1, 
            args->maxIterations, args->output);
    }
    double endTime = CycleTimer::currentSeconds();
    minThread = std::min(minThread, endTime - startTime);
    printf("Thread %d time: %.3f\n", args->threadId, minThread);
    return NULL;
}

//
// MandelbrotThread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Multi-threading performed via pthreads.
void mandelbrotThread(
    int numThreads,
    float x0, float y0, float x1, float y1,
    int width, int height,
    int maxIterations, int output[])
{
    const static int MAX_THREADS = 32;

    if (numThreads > MAX_THREADS)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
        exit(1);
    }

    pthread_t workers[MAX_THREADS];
    WorkerArgs args[MAX_THREADS];

    for (int i=0; i<numThreads; i++) {
      
        // TODO FOR CS348V STUDENTS: You may or may not wish to modify
        // the per-thread arguments here.  The code below copies the
        // same arguments for each thread
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = maxIterations;
        args[i].numThreads = numThreads;
        args[i].output = output;
      
        args[i].threadId = i;
    }

    // Spawn the worker threads.  Note that numThreads-1 pthreads
    // are created and the main app thread is used as a worker as
    // well.

    for (int i=1; i<numThreads; i++)
        pthread_create(&workers[i], NULL, workerThreadStart, &args[i]);
    workerThreadStart(&args[0]);

    // wait for worker threads to complete
    for (int i=1; i<numThreads; i++)
        pthread_join(workers[i], NULL);
}
