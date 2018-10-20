#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <assert.h>
#include "mpi.h"
#include "wire_route.h"

static int g_num_rows = 0;
static int g_num_cols = 0;
static int g_delta = 0;
static int g_num_wires = 0;
wire_t *wires = NULL;
cost_t *costs = NULL;
cost_t *path_costs = NULL;
cost_t *temp_path_costs = NULL;

cost_t max(cost_t x, cost_t y) {
    return x > y ? x : y;
}

cost_t min(cost_t x, cost_t y) {
    return x < y ? x : y;
}

// Initialize problem
static inline void init(int numRows, int numCols, int delta, int numWires)
{
    g_num_rows = numRows;
    g_num_cols = numCols;
    g_delta = delta;
    g_num_wires = numWires;
    wires = (wire_t*)calloc(g_num_wires, sizeof(wire_t));
    costs = (cost_t*)calloc(g_num_rows * g_num_cols, sizeof(cost_t));
    path_costs = (cost_t*)malloc( (g_num_rows + g_num_cols)* sizeof(cost_t));
    temp_path_costs = (cost_t*)malloc( (g_num_rows + g_num_cols)* sizeof(cost_t));
}


// Initialize a given wire
static inline void initWire(int wireIndex, int x1, int y1, int x2, int y2)
{
    wire_t wire;
    wire.x1 = x1;
    wire.x2 = x2;
    wire.y1 = y1;
    wire.y2 = y2;
    wire.bend_x1 = -1;
    wire.bend_y1 = -1;
    wire.bend_x2 = -1;
    wire.bend_y2 = -1;
    wire.cost = -1;
    wires[wireIndex] = wire;    
}

// Return number of rows
static inline int getNumRows()
{
    //TODO Implement code here
    return g_num_rows;
}

// Return number of cols
static inline int getNumCols()
{
    //TODO Implement code here
    return g_num_cols;
}

// Return delta
static inline int getDelta()
{
	//TODO Implement code here
    return g_delta;
}

// Return number of wires
static inline int getNumWires()
{
    //TODO Implement code here
    return 0;
}

// Get cost array entry
static inline int getCost(int row, int col)
{
    //TODO Implement code here
    return 0;
}

// Get a wire placement. Returns number of points (should be 2-4 for 0-2 bends).
static inline int getWire(int wireIndex, int* x1, int* y1, int* x2, int* y2, int* x3, int* y3, int* x4, int* y4)
{
    //TODO Implement code here
    return 2;
}


static inline void find_min_path(wire_t wire, double anneal_prob, int *horizontal,
                int *vertical, int *max_overlap_horiz, int *max_overlap_vert) {
    int x1, x2, y1, y2;
    if (wire.x1 <= wire.x2) {
        x1 = wire.x1;
        y1 = wire.y1;
        x2 = wire.x2;
        y2 = wire.y2;
    } else {
        x1 = wire.x2;
        y1 = wire.y2;
        x2 = wire.x1;
        y2 = wire.y1;
    }
    int x_max = min(g_num_cols-1, (int)(g_delta/2) + x2);
    int x_min = max(0, x1 - (int)(g_delta/2));
    int y_max = min(g_num_rows-1, (int)(g_delta/2) + max(y1, y2));
    int y_min = max(0, min(y1, y2) - (int)(g_delta/2));
    if (wire.y1 == wire.y2) {
        x_min = x1;
        x_max = x2;
    }   
    if (wire.x1 == wire.x2) {
        y_min = y1;
        y_max = y2;
    }
    if (wire.cost != -1) { 
        wire.cost = find_wire_cost(wire);
        change_wire_route(wire, -1);
    }
    memset(horizontal, 0, sizeof(cost_t) * g_num_rows);
    memset(vertical, 0, sizeof(cost_t) * g_num_cols);
    memset(max_overlap_horiz, 0, sizeof(cost_t) * g_num_rows);
    memset(max_overlap_vert, 0, sizeof(cost_t) * g_num_cols);
    int i;
    for (i = y_min; i <= y_max; i++) {
        create_vert_horiz(i, wire, horizontal, vertical, 
                 max_overlap_horiz, max_overlap_vert);
    }
}

static inline void wire_routing(double anneal_prob, int *horizontal, int *vertical,
                                int *max_overlap_horiz, int *max_overlap_vert) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    for (int i = 0; i < g_num_wires; i++) {
        find_min_path(wires[i], anneal_prob, int *horizontal, int *vertical,
                                int *max_overlap_horiz, int *max_overlap_vert);    
    }
}

// Perform computation, including reading/writing output files
void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations)
{
    readInput(inputFilename);
    //TODO Implement code here
    //TODO Decide which processors should be reading/writing files
    writeCost(inputFilename, nproc);
    writeOutput(inputFilename, nproc);
}

// Read input file
void readInput(char* inputFilename)
{
    FILE* fp = fopen(inputFilename, "r");
    int numRows;
    int numCols;
    int delta;
    int numWires;
    int wireIndex;
    int x1;
    int y1;
    int x2;
    int y2;
    if (fp == NULL) {
        fprintf(stderr, "Failed to read input file %s\n", inputFilename);
        exit(-1);
    }
    if (fscanf(fp, "%d %d", &numCols, &numRows) != 2) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if ((numRows <= 0) || (numCols <= 0)) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (fscanf(fp, "%d", &delta) != 1) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (fscanf(fp, "%d", &numWires) != 1) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    if (numWires <= 0) {
        fprintf(stderr, "Invalid input file format\n");
        exit(-1);
    }
    init(numRows, numCols, delta, numWires);
    for (wireIndex = 0; wireIndex < numWires; wireIndex++) {
        if (fscanf(fp, "%d %d %d %d", &x1, &y1, &x2, &y2) != 4) {
            fprintf(stderr, "Invalid input file format\n");
            exit(-1);
        }
        initWire(wireIndex, x1, y1, x2, y2);
    }
    fclose(fp);
}

// Write cost array file based on input filename
void writeCost(char* inputFilename, int nproc)
{
    char* dname = strdup(inputFilename);
    char* bname = strdup(inputFilename);
    char* costsFilename = malloc(strlen(inputFilename) + 100);
    FILE* fp;
    int numRows = getNumRows();
    int numCols = getNumCols();
    int row;
    int col;
    assert(costsFilename != NULL);
    sprintf(costsFilename, "%s/costs_%s_%d.txt", dirname(dname), basename(bname), nproc);
    fp = fopen(costsFilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Failed to write costs file %s\n", costsFilename);
        exit(-1);
    }
    fprintf(fp, "%d %d\n", numCols, numRows);
    for (row = 0; row < numRows; row++) {
        for (col = 0; col < numCols; col++) {
            fprintf(fp, "%d ", getCost(row, col));
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    free(costsFilename);
    free(bname);
    free(dname);
}

// Write wire output file based on input filename
void writeOutput(char* inputFilename, int nproc)
{
    char* dname = strdup(inputFilename);
    char* bname = strdup(inputFilename);
    char* outputFilename = malloc(strlen(inputFilename) + 100);
    FILE* fp;
    int numRows = getNumRows();
    int numCols = getNumCols();
    int delta = getDelta();
    int numWires = getNumWires();
    int wireIndex;
    int numPoints;
    int x1;
    int y1;
    int x2;
    int y2;
    int x3;
    int y3;
    int x4;
    int y4;
    assert(outputFilename != NULL);
    sprintf(outputFilename, "%s/output_%s_%d.txt", dirname(dname), basename(bname), nproc);
    fp = fopen(outputFilename, "w");
    if (fp == NULL) {
        fprintf(stderr, "Failed to write output file %s\n", outputFilename);
        exit(-1);
    }
    fprintf(fp, "%d %d\n", numCols, numRows);
    fprintf(fp, "%d\n", delta);
    fprintf(fp, "%d\n", numWires);
    for (wireIndex = 0; wireIndex < numWires; wireIndex++) {
        numPoints = getWire(wireIndex, &x1, &y1, &x2, &y2, &x3, &y3, &x4, &y4);
        switch (numPoints) {
            case 2:
                fprintf(fp, "%d %d %d %d\n", x1, y1, x2, y2);
                break;

            case 3:
                fprintf(fp, "%d %d %d %d %d %d\n", x1, y1, x2, y2, x3, y3);
                break;

            case 4:
                fprintf(fp, "%d %d %d %d %d %d %d %d\n", x1, y1, x2, y2, x3, y3, x4, y4);
                break;

            default:
                assert(0); // invalid number of points
                break;
        }
    }
    fclose(fp);
    free(outputFilename);
    free(bname);
    free(dname);
}
