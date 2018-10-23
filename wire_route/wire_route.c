#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <libgen.h>
#include <assert.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include "mpi.h"
#include "wire_route.h"

//global variables
static int g_num_rows = 0;
static int g_num_cols = 0;
static int g_delta = 0;
static int g_num_wires = 0;
wire_t *wires = NULL;
cost_t *costs = NULL;

typedef struct
{
    int first;
    int second;
} pair_t;

// returns whether the pair x is better than the pair y, assuming the pairs 
// reprsent a wire path cost of the form (max_cost, agg_cost)
bool better_than(pair_t x, pair_t y) {
    return x.first < y.first || (x.first == y.first && x.second < y.second);
}

bool worse_than_equal_to(pair_t x, pair_t y) {
    return x.first > y.first || (x.first == y.first && x.second >= y.second);
}

int max(int x, int y) {
    return x > y ? x : y;
}

int min(int x, int y) {
    return x < y ? x : y;
}

void print(cost_t *array, int rows, int cols) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++)  
            printf("%d ", array[i*cols + j]);
        printf("\n");
    }
    printf("\n");
}

cost_t costs_max_overlap() {
    cost_t max_cost = 0;
    for (int i = 0; i < g_num_cols * g_num_rows; i++) {
        max_cost = max(max_cost, costs[i]);
    }
    return max_cost;
}

cost_t costs_agg_overlap() {
    cost_t agg_cost = 0;
    for (int i = 0; i < g_num_cols * g_num_rows; i++) {
        agg_cost += costs[i] * costs[i];
    }
    return agg_cost;
}

pair_t find_score_h(int i, int x1, int y1, int x2, int y2, pair_t best) {
    pair_t result;
    result.first = 0;
    result.second = 0;

    int dir = x1 < x2 ? 1 : -1;
    for (int x = x1; x != x2; x += dir) {
        cost_t cost = 1 + costs[g_num_cols * i + x];
        result.first = max(result.first, cost);
        result.second += cost;
    }
    if (worse_than_equal_to(result, best)) {
        return result;
    }



    dir = y1 < i ? 1 : -1;
    for (int y = y1; y != i; y += dir) {
        cost_t cost = 1 + costs[g_num_cols * y + x1];
        result.first = max(result.first, cost);
        result.second += cost;

        if (worse_than_equal_to(result, best)) {
            return result;
        }
    }
    
    dir = i < y2 ? 1 : -1;
    for (int y = i; y != y2; y += dir) {
        cost_t cost = 1 + costs[g_num_cols * y + x2];
        result.first = max(result.first, cost);
        result.second += cost;

        if (worse_than_equal_to(result, best)) {
            return result;
        }
    }
    cost_t cost = 1 + costs[g_num_cols * y2 + x2];
    result.first = max(result.first, cost);
    result.second += cost;

    return result;
}

pair_t find_score_v(int i, int x1, int y1, int x2, int y2, pair_t best) {
    pair_t result;
    result.first = 0;
    result.second = 0;

    int dir = x1 < i ? 1 : -1;
    for (int x = x1; x != i; x += dir) {
        cost_t cost = 1 + costs[g_num_cols * y1 + x];
        result.first = max(result.first, cost);
        result.second += cost;
    }

    if (worse_than_equal_to(result, best)) {
        return result;
    }

    dir = i < x2 ? 1 : -1;
    for (int x = i; x != x2; x += dir) {
        cost_t cost = 1 + costs[g_num_cols * y2 + x];
        result.first = max(result.first, cost);
        result.second += cost;
    }
    cost_t cost = 1 + costs[g_num_cols * y2 + x2];
    result.first = max(result.first, cost);
    result.second += cost;

    dir = y1 < y2 ? 1 : -1;
    for (int y = y1; y != y2; y += dir) {
        cost_t cost = 1 + costs[g_num_cols * y + i];
        result.first = max(result.first, cost);
        result.second += cost;

        if (worse_than_equal_to(result, best)) {
            return result;
        }
    }

    return result;
}

void change_wire_route_helper(int x1, int y1, int x2, int y2,
                              int increment) {
    if (x1 == x2 && y1 == y2)
        return; 
    if (y1 == y2) {
        if (x1 <= x2) {
            for (int i = x1; i < x2; i++) {
                costs[ y1 * g_num_cols + i] += increment;
            }
        } else {
            for (int i = x1; i > x2; i-=1) {
                costs[y1 * g_num_cols + i] += increment;
            }
        }
    } else {
        if (y1 <= y2) {
            for (int i = y1; i < y2; i++) {
                costs[i * g_num_cols + x1] += increment;
            }
        } else {
            for (int i = y1; i > y2; i--) {
                costs[i * g_num_cols + x1] += increment;
            }

        }
    }
}

void change_wire_route(wire_t wire, int increment) {
    if (wire.bend_x1 == -1 && wire.bend_y1 == -1) {
        change_wire_route_helper(wire.x1, wire.y1, wire.x2, 
            wire.y2, increment);
        costs[wire.y2 * g_num_cols + wire.x2] += increment; 
        return;
    }
    change_wire_route_helper(wire.x1, wire.y1, wire.bend_x1, wire.bend_y1,
                             increment);
    if (wire.bend_y2 >= 0 && wire.bend_x2 >= 0) {
        change_wire_route_helper(wire.bend_x1, wire.bend_y1,
                        wire.bend_x2, wire.bend_y2, increment);
        change_wire_route_helper(wire.bend_x2, wire.bend_y2,
                        wire.x2, wire.y2, increment);
    } else {
        change_wire_route_helper(wire.bend_x1, wire.bend_y1,
                        wire.x2, wire.y2, increment);
    }
    costs[wire.y2 * g_num_cols + wire.x2] += increment; 
}

void anneal(wire_t *wire) {
    // if (wire->cost != -1)
    int dx = abs(wire->x2 - wire->x1) + 1;
    int dy = abs(wire->y2 - wire->y1) + 1;
    int random_index = rand() % (dx + dy);
    if (random_index < dx) {
        // this will be a 'vertical' path
        if (random_index == 0) {
            wire->bend_x1 = wire->x1;
            wire->bend_y1 = wire->y2;
            wire->bend_x2 = -1;
            wire->bend_y2 = -1;
        }
        else {
            int dir = ceil((double)(wire->x2 - wire->x1) / dx);
            wire->bend_x1 = wire->x1 + dir * (random_index);
            wire->bend_y1 = wire->y1;
            wire->bend_x2 = wire->bend_x1 == wire->x2 ? -1 : wire->bend_x1;
            wire->bend_y2 = wire->bend_x1 == wire->x2 ? -1 : wire->y2;
        }
    }
    else {
        // this will be a 'horizontal' path
        random_index -= dx;
        if (random_index == 0) {
            wire->bend_x1 = wire->x2;
            wire->bend_y1 = wire->y1;
            wire->bend_x2 = -1;
            wire->bend_y2 = -1;
        }
        else {
            int dir = ceil((double)(wire->y2 - wire->y1) / dy);
            wire->bend_x1 = wire->x1;
            wire->bend_y1 = wire->y1 + dir * random_index;
            wire->bend_x2 = wire->bend_y1 == wire->y2 ? -1 : wire->x2;
            wire->bend_y2 = wire->bend_y1 == wire->y2 ? -1 : wire->bend_y1;
        }
    }
}

// find_mind_path_cost: takes in (wire.x1, wire.y1) to (wire.x2, wire.y2) and 
// adds the wire route to costs
void find_min_path(wire_t *wire, double anneal_prob) {
    /*double prob_sample = ((double)rand()) / ((double)RAND_MAX);
    if (prob_sample < anneal_prob) {
        change_wire_route(*wire, -1);
        anneal(wire);
        change_wire_route(*wire, 1);
        return;
    }*/
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    int x1, x2, y1, y2;
    
    if (wire->x1 <= wire->x2) {
        x1 = wire->x1;
        y1 = wire->y1;
        x2 = wire->x2;
        y2 = wire->y2;
    } else {
        x1 = wire->x2;
        y1 = wire->y2;
        x2 = wire->x1;
        y2 = wire->y1;
    }
    // printf("wire: %d %d %d %d\n", x1, y1, x2, y2);
    int x_max = min(g_num_cols-1, (int)(g_delta/2) + x2);
    int x_min = max(0, x1 - (int)(g_delta/2));
    int y_max = min(g_num_rows-1, (int)(g_delta/2) + max(y1, y2));
    int y_min = max(0, min(y1, y2) - (int)(g_delta/2));
    if (wire->y1 == wire->y2) {
        x_min = x1;
        x_max = x2;
    }   
    if (wire->x1 == wire->x2) {
        y_min = y1;
        y_max = y2;
    }

    change_wire_route(*wire, -1);
    //printf("%d %d %d %d\n", wire->bend_x1, wire->bend_y1, wire->bend_x2,
    //        wire->bend_y2);
    // print(costs, g_num_rows, g_num_cols);
   
    // MAIN IDEA: 
    // horizontal keeps track of the path cost of a wire route that has a 
    // bend at row x
    // vertical keeps track of the path cost of a wire route that has a bend 
    // at column y
   
    int new_bendx1 = -1;
    int new_bendy1 = -1;
    int new_bendx2 = -1;
    int new_bendy2 = -1;

    // find optimal wire route
    pair_t best_score;
    best_score.first = INT_MAX;
    best_score.second = INT_MAX;
    
    int flag = 0;
    int index = -1;

    double frac1 = (double)world_rank / (double)world_size;
    double frac2 = (double)(world_rank+1) / (double)world_size;

    int y_start = y_min + ceil(frac1 * (double)(y_max - y_min));
    int y_end = y_min + ceil(frac2 * (double)(y_max - y_min));
    // horizontal
    for (int i = y_start; i <= y_end; i++) {
        pair_t current_score = find_score_h(i, x1, y1, x2, y2, best_score);
        if (better_than(current_score, best_score)) {
            index = i; 
            best_score = current_score;
        }
    }
    int x_start = x_min + ceil(frac1 * (double)(x_max - x_min ));
    int x_end = x_min + ceil(frac2 * (double)(x_max - x_min ));

   for (int i = x_start; i < x_end; i++) {
        pair_t current_score = find_score_v(i, x1, y1, x2, y2, best_score);
        if (better_than(current_score, best_score)) {
            best_score = current_score;
            flag = 1;
            index = i;
        }
    }
    // printf("%d %d %d %d\n", wire->x1, wire->y1, wire->x2, wire->y2);
    // printf("y: %d %d x: %d %d, score %d, flag %d, index %d\n", y_start, y_end,
            // x_start, x_end, best_score.first, flag, index);

    // send data to master process
    int num = 4;
    int temp[num];
    temp[0] = best_score.first;
    temp[1] = best_score.second;
    temp[2] = flag;
    temp[3] = index;
    MPI_Send(&temp, num, MPI_INT, 0, 0, MPI_COMM_WORLD);
    // gather data to master process
    int *buffer = NULL;
    if (world_rank == 0) {
        buffer = malloc(sizeof(int) * num * world_size);
    }   
    MPI_Gather(&temp, num, MPI_INT, buffer, num, MPI_INT, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        /*for (int i = 0; i < world_size; i++) {
            for (int j = 0; j < num; j++) {
                printf("%d ", buffer[num*i + j]);
            }
            printf("\n");
        }*/
        int best = INT_MAX;
        int agg = INT_MAX;
        int final_flag = -1;
        int final_ind = -1;
        for (int i = 0; i < num*world_size; i+=4) {
            if (best > buffer[i]) {
                best = buffer[i];
                agg = buffer[i+1];
                final_flag = buffer[i+2];
                final_ind = buffer[i+3];
            } else if (best == buffer[i] && agg > buffer[i+1]) {
                best = buffer[i];
                agg = buffer[i+1];
                final_flag = buffer[i+2];
                final_ind = buffer[i+3];
            }
        }
        printf("%d %d\n", final_ind, final_flag);

        if (final_flag) { // vertical
            if ( !(wire->x1 == wire->x2 && final_ind == wire->x1) ) {
                new_bendx1 = final_ind;
                if (final_ind == wire->x1) {
                    new_bendy1 = wire->y2;
                } else {
                    // if horizontal segment aligns w/ y2 or != y1
                    new_bendy1 = wire->y1;
                } if (final_ind != wire->x1 && final_ind != wire->x2) {
                    new_bendy2 = wire->y2;
                    new_bendx2 = final_ind;
                }
            }
        } else {
            if ( !(wire->y1 == wire->y2 && final_ind == wire->y1) ) {
                // if horizontal segment is align y1
                new_bendy1 = final_ind;
                if (final_ind == wire->y1) {
                    new_bendx1 = wire->x2;
                } else {
                    // if horizontal segment aligns w/ y2 or != y1
                    new_bendx1 = wire->x1;
                } if (final_ind != wire->y1 && final_ind != wire->y2) {
                    new_bendy2 = final_ind;
                    new_bendx2 = wire->x2;
                }
            }
        }
        wire->bend_x1 = new_bendx1;
        wire->bend_y1 = new_bendy1;
        wire->bend_x2 = new_bendx2;
        wire->bend_y2 = new_bendy2;
        // printf("%d %d %d %d\n", wire->x1, wire->y1, wire->x2, wire->y2);
        // printf("%d %d %d %d\n", new_bendx1, new_bendy1, new_bendx2, new_bendy2);
        change_wire_route(*wire, 1);
    }
}

// Initialize problem
static inline void init(int numRows, int numCols, int delta, int numWires)
{
    // only initialize random seed once
    srand(time(NULL));
    g_num_rows = numRows;
    g_num_cols = numCols;
    g_delta = delta;
    g_num_wires = numWires;
    wires = (wire_t*)calloc(g_num_wires, sizeof(wire_t));
    costs = (cost_t*)calloc(g_num_rows * g_num_cols, sizeof(cost_t));
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
    wires[wireIndex] = wire;    
}

// Return number of rows
static inline int getNumRows()
{
    return g_num_rows;
}

// Return number of cols
static inline int getNumCols()
{
    return g_num_cols;
}

// Return delta
static inline int getDelta()
{
    return g_delta;
}

// Return number of wires
static inline int getNumWires()
{
    return g_num_wires;
}

// Get cost array entry
static inline int getCost(int row, int col)
{
    return costs[g_num_cols * row + col];
}

// Get a wire placement. Returns number of points (should be 2-4 for 0-2 bends).
static inline int getWire(int wireIndex, int* x1, int* y1, int* x2, int* y2, 
                          int* x3, int* y3, int* x4, int* y4)
{
    wire_t wire = wires[wireIndex];
    *x1 = wire.x1;
    *y1 = wire.y1;
    if(wire.bend_x1 == -1 || wire.bend_y1 == -1) {
        *x2 = wire.x2;
        *y2 = wire.y2;
        return 2;
    } else if (wire.bend_x2 == -1 || wire.bend_y2 == -1) {
        *x2 = wire.bend_x1;
        *y2 = wire.bend_y1;
        *x3 = wire.x2;
        *y3 = wire.y2;
        return 3;
    } else {
        *x2 = wire.bend_x1;
        *y2 = wire.bend_y1;
        *x3 = wire.bend_x2;
        *y3 = wire.bend_y2;
        *x4 = wire.x2;
        *y4 = wire.y2;
        return 4;
    }
}

static inline void wire_routing(double anneal_prob) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    for (int i = 0; i < g_num_wires; i++) {
        find_min_path(&wires[i], anneal_prob);
    }
}

// Perform computation, including reading/writing output files
void compute(int procID, int nproc, char* inputFilename, double prob, 
             int numIterations)
{
    readInput(inputFilename);

    for (int i = 0; i < g_num_wires; i++) {
        anneal(&wires[i]);
        change_wire_route(wires[i], 1);
    }

    for (int i = 0; i < numIterations; i++) {
        wire_routing(prob); 
    }

    //TODO Implement code here
    //TODO Decide which processors should be reading/writing files
    writeCost(inputFilename, nproc);
    writeOutput(inputFilename, nproc);

    cost_t max_cost = costs_max_overlap();
    printf("Max cost: %d\n", max_cost);
    

    cost_t agg_cost = costs_agg_overlap();
    printf("Agg cost: %d\n", agg_cost);
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
