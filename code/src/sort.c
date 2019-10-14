#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <omp.h>
#include "sort.h"
#include "edgelist.h"

// Order edges by id of a source vertex,
// using the Counting Sort
// Complexity: O(E + V)
#define NUM_THREADS 16

void countSortEdgesBySource(struct Edge *edges_sorted, struct Edge *edges, int numVertices, int numEdges) {
    int i;
    int key;
    int pos;

    // auxiliary arrays, allocated at the start up of the program
    int *vertex_cnt = (int*)malloc(numVertices*sizeof(int)); // needed for Counting Sort

    for(i = 0; i < numVertices; ++i) {
        vertex_cnt[i] = 0;
    }

    // count occurrence of key: id of a source vertex
    for(i = 0; i < numEdges; ++i) {
        key = edges[i].src;
        vertex_cnt[key]++;
    }

    // transform to cumulative sum
    for(i = 1; i < numVertices; ++i) {
        vertex_cnt[i] += vertex_cnt[i - 1];
    }

    // fill-in the sorted array of edges
    for(i = numEdges - 1; i >= 0; --i) {
        key = edges[i].src;
        pos = vertex_cnt[key] - 1;
        edges_sorted[pos] = edges[i];
        vertex_cnt[key]--;
    }


    free(vertex_cnt);

}

void radixSortEdgesBySource(struct Edge *edges_sorted, struct Edge *edges, int numVertices, int numEdges) 
{
   int pos,i,key;
int *vertex_cnt = (int*)malloc(10*sizeof(int));
int max = edges[0].src;
for (i = 1; i<numEdges;i++)
{
 if (edges[i].src > max)
{
max = edges[i].src;
}
}
int x = 1;
while(max!=0)
{
   for(i = 0; i < 10; ++i) {
        vertex_cnt[i] = 0;
    }
    for(i = 0; i < numEdges; ++i) {
        key = (edges[i].src/x) % 10;
        vertex_cnt[key]++;
    }
    for(i = 1; i < 10; ++i) {
        vertex_cnt[i] += vertex_cnt[i - 1];
    }
    for(i = numEdges - 1; i >= 0; --i) {
        key = (edges[i].src/x) % 10;
        pos = vertex_cnt[key] - 1;
        edges_sorted[pos] = edges[i];
        vertex_cnt[key]--;
    }
max = max/10;
x = x * 10;
}
free(vertex_cnt);
}

void openmp_radixSortEdgesBySource(struct Edge *edges_sorted, struct Edge *edges, int numVertices, int numEdges) 
{
int max = edges[0].src;
//printf("%d", omp_get_max_threads());
for (int k = 1; k<numEdges;k++)
{
 if (edges[k].src > max)
{
max = edges[k].src;
}
}
int x = 1;
while(max!=0)
{
openmp_countSortEdgesBySource(edges_sorted,edges, numVertices, numEdges, x);
max = max/10;
x = x * 10;
}
}

void openmp_countSortEdgesBySource(struct Edge *edges_sorted, struct Edge *edges, int numVertices, int numEdges, int x)
{
   int base=0,i;
    int **vertex_cnt;
    vertex_cnt = (int **)malloc(NUM_THREADS*sizeof(int *));
    for(i=0;i<NUM_THREADS;i++)
    vertex_cnt[i] = (int *)malloc(10*sizeof(int));

omp_set_num_threads(NUM_THREADS);
int nthread;
#pragma omp parallel num_threads(NUM_THREADS)
{
int tid,nthrds,key,i,j,tid_start,tid_end;
tid = omp_get_thread_num();
nthrds = omp_get_num_threads();
tid_start = tid * (numEdges/nthrds);
tid_end = (tid + 1) * (numEdges/nthrds);
if(tid ==0) nthread = nthrds;
if(tid == (nthrds - 1)) tid_end = numEdges;
    for(i=0;i<10;i++){vertex_cnt[tid][i] = 0;}
    for(j = tid_start; j < tid_end; ++j)
    {
        key = (edges[j].src/x) % 10;
        vertex_cnt[tid][key]++;
    }
 #pragma omp barrier
}
   for(i=0;i<10;i++)
    {
        for(int j=0;j<NUM_THREADS;j++)
        {
            base = base + vertex_cnt[j][i];
            vertex_cnt[j][i] = base;
        }
    }

omp_set_num_threads(NUM_THREADS);

#pragma omp parallel num_threads(NUM_THREADS)
{
int pos,tid,nthrds,tid_start,tid_end,key;
tid = omp_get_thread_num();
nthrds = omp_get_num_threads();
tid_start = tid * (numEdges / nthrds );
tid_end = ( tid + 1 ) *( numEdges / nthrds);
if (tid == 0) nthread = nthrds;
if((nthrds - 1) == tid) tid_end = numEdges;

for(int j=(tid_end - 1);j>= tid_start;j--)
    {
        key = (edges[j].src/x) % 10;
        pos = ((vertex_cnt[tid][key]) - 1);
        edges_sorted[pos] = edges[j];
        vertex_cnt[tid][key] = (vertex_cnt[tid][key]) - 1;
    }
}

    free(vertex_cnt);
}
