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
cost_t *costs_T = NULL; // transposed cost matrix

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

static inline void swap(int i, int j) {
    wire_t temp = wires[i];
    wires[i] = wires[j];
    wires[j] = temp;
}

static void sort_wires() {
    return;

    for (int i = 0; i < g_num_wires; i++) {
        int random_index = rand() % (g_num_wires);
        swap(i, random_index);
    }
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
        // cost_t cost = 1 + costs[g_num_cols * y + x1];
        cost_t cost = 1 + costs_T[x1 * g_num_rows + y];
        result.first = max(result.first, cost);
        result.second += cost;

    }
    if (worse_than_equal_to(result, best)) {
        return result;
    }

    
    dir = i < y2 ? 1 : -1;
    for (int y = i; y != y2; y += dir) {
        // cost_t cost = 1 + costs[g_num_cols * y + x2];
        cost_t cost = 1 + costs_T[x2 * g_num_rows + y];
        result.first = max(result.first, cost);
        result.second += cost;
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

    if (worse_than_equal_to(result, best)) {
        return result;
    }


    dir = y1 < y2 ? 1 : -1;
    for (int y = y1; y != y2; y += dir) {
        // cost_t cost = 1 + costs[g_num_cols * y + i];
        cost_t cost = 1 + costs_T[i * g_num_rows + y];
        result.first = max(result.first, cost);
        result.second += cost;
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
                costs[y1 * g_num_cols + i] += increment;
                costs_T[i * g_num_rows + y1] += increment;
            }
        } else {
            for (int i = x1; i > x2; i-=1) {
                costs[y1 * g_num_cols + i] += increment;
                costs_T[i * g_num_rows + y1] += increment;
            }
        }
    } else {
        if (y1 <= y2) {
            for (int i = y1; i < y2; i++) {
                costs[i * g_num_cols + x1] += increment;
                costs_T[x1 * g_num_rows + i] += increment;
            }
        } else {
            for (int i = y1; i > y2; i--) {
                costs[i * g_num_cols + x1] += increment;
                costs_T[x1 * g_num_rows + i] += increment;
            }

        }
    }
}

void change_wire_route(wire_t wire, int increment) {
    if (wire.bend_x1 == -1 && wire.bend_y1 == -1) {
        change_wire_route_helper(wire.x1, wire.y1, wire.x2, 
            wire.y2, increment);
        costs[wire.y2 * g_num_cols + wire.x2] += increment; 
        costs_T[wire.x2 * g_num_rows + wire.y2] += increment;
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
    costs_T[wire.x2 * g_num_rows + wire.y2] += increment;
}

pair_t anneal(wire_t wire) {

    int x1 = wire.x1;
    int y1 = wire.y1;
    int x2 = wire.x2;
    int y2 = wire.y2;

    int x_max = min(g_num_cols-1, (int)(g_delta/2) + max(x1, x2));
    int x_min = max(0, min(x1, x2) - (int)(g_delta/2));
    int y_max = min(g_num_rows-1, (int)(g_delta/2) + max(y1, y2));
    int y_min = max(0, min(y1, y2) - (int)(g_delta/2));

    int dx = x_max - x_min + 1;
    int dy = y_max - y_min + 1;
    int random_index = rand() % (dx + dy);
    pair_t new_path;
    if (random_index < dx) {
        // this will be a 'vertical' path        
        new_path.first = 1;
        new_path.second = x_min + random_index;
    }
    else {
        // this will be a 'horizontal' path
        random_index -= dx;
        new_path.first = 0;
        new_path.second = y_min + random_index;
    }

    return new_path;
}

void set_wire(wire_t *wire, int vertical, int location) {
    int new_bendx1 = -1;
    int new_bendy1 = -1;
    int new_bendx2 = -1;
    int new_bendy2 = -1;

    if (vertical == 1) {
        if ( !(wire->x1 == wire->x2 && location == wire->x1) ) {
            new_bendx1 = location;
            if (location == wire->x1) {
                new_bendy1 = wire->y2;
            } else {
                new_bendy1 = wire->y1;
            } if (location != wire->x1 && location != wire->x2) {
                new_bendy2 = wire->y2;
                new_bendx2 = location;
            }
        }
    } else { // horizontal
        if ( !(wire->y1 == wire->y2 && location == wire->y1) ) {
            new_bendy1 = location;
            if (location == wire->y1) {
                new_bendx1 = wire->x2;
            } else {
                new_bendx1 = wire->x1;
            } if (location != wire->y1 && location != wire->y2) {
                new_bendy2 = location;
                new_bendx2 = wire->x2;
            }
        }
    }

    wire->bend_x1 = new_bendx1;
    wire->bend_y1 = new_bendy1;
    wire->bend_x2 = new_bendx2;
    wire->bend_y2 = new_bendy2;
}

// find_mind_path_cost: takes in (wire.x1, wire.y1) to (wire.x2, wire.y2) and 
// adds the wire route to costs
pair_t find_min_path(wire_t wire, double anneal_prob, 
                     int base_rank, int top_rank) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    MPI_Status status;

    pair_t new_path;

    int annealed_size = 2;
    int annealed_path[annealed_size];

    if (world_rank == base_rank) {
        double prob_sample = ((double)rand()) / ((double)RAND_MAX);
        if (prob_sample < anneal_prob) {
            new_path = anneal(wire);

            annealed_path[0] = new_path.first;
            annealed_path[1] = new_path.second;

            for (int i = base_rank + 1; i < top_rank; i++){
                MPI_Send(&annealed_path, annealed_size, MPI_INT, i, 0, 
                         MPI_COMM_WORLD);
            } 
            return new_path;
        } else {
            annealed_path[0] = -1;
            annealed_path[1] = -1;

            for (int i = base_rank + 1; i < top_rank; i++){
                MPI_Send(&annealed_path, annealed_size, MPI_INT, i, 0, 
                         MPI_COMM_WORLD);
            } 
        }
    } else {
        MPI_Recv(&annealed_path, annealed_size, MPI_INT, base_rank, 0,
                      MPI_COMM_WORLD, &status);

        if (annealed_path[0] != -1 && annealed_path[1] != -1) {
            new_path.first = annealed_path[0];
            new_path.second = annealed_path[1];
            return new_path;
        }
    }

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
 
    // find optimal wire route
    pair_t best_score;
    best_score.first = INT_MAX;
    best_score.second = INT_MAX;    
    int vertical = 0;
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
            vertical = 1;
            index = i;
        }
    }

    int message_size = 4;
    int temp[message_size];

    if (world_rank == base_rank) {
        // gather data to master process
        pair_t temp_score;
        for (int i = base_rank + 1; i < top_rank; i++) {
            MPI_Recv(&temp, message_size, MPI_INT, i, 1, MPI_COMM_WORLD, 
                     &status);
            temp_score.first = temp[0];
            temp_score.second = temp[1];
            if (better_than(temp_score, best_score)) {
                best_score = temp_score;
                vertical = temp[2];
                index = temp[3];
            }
        }

        temp[0] = best_score.first;
        temp[1] = best_score.second;
        temp[2] = vertical;
        temp[3] = index;

        for (int i = base_rank + 1; i < top_rank; i++){
            MPI_Send(&temp, message_size, MPI_INT, i, 2,
                     MPI_COMM_WORLD);
        }
    } else {
        // send data to master process
        temp[0] = best_score.first;
        temp[1] = best_score.second;
        temp[2] = vertical;
        temp[3] = index;

        MPI_Send(&temp, message_size, MPI_INT, base_rank, 1, MPI_COMM_WORLD);

        MPI_Recv(&temp, message_size, MPI_INT, base_rank, 2,
                 MPI_COMM_WORLD, &status);
    }

    new_path.first = temp[2];
    new_path.second = temp[3];
    return new_path;
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
    costs_T = (cost_t*)calloc(g_num_rows * g_num_cols, sizeof(cost_t));
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
    MPI_Status status;

    int parallelized_wires = 2;
    parallelized_wires = min(parallelized_wires, world_size);
    int processes_per_wire = (world_size + parallelized_wires - 1) / 
                              parallelized_wires;

    // stores all the new paths
    int new_paths_size = 2 * parallelized_wires;
    int new_paths[new_paths_size];

    // stores the new path computed by this process
    int new_path[2] = {-1, -1};

    for (int i = 0; i < g_num_wires; i += parallelized_wires) {

        // wire_index is the index of the wire in this batch
        int wire_index = world_rank / processes_per_wire;
        int base_rank = wire_index * processes_per_wire;
        int top_rank = min(base_rank + processes_per_wire, world_size);

        wire_index += i; // wire_index is now the index of the wire in 
                             // the wires array

        // compute the wire's new path
        if (wire_index < g_num_wires) {
            change_wire_route(wires[wire_index], -1);

            pair_t new_path_pair = find_min_path(wires[wire_index],
                                    anneal_prob, base_rank, top_rank);


            new_path[0] = new_path_pair.first;
            new_path[1] = new_path_pair.second;
    
            set_wire(&wires[wire_index], new_path[0], new_path[1]);
            change_wire_route(wires[wire_index], 1);  
        }

        // report new path to the master
        if (world_rank > 0 && world_rank % processes_per_wire == 0) {
            MPI_Send(&new_path, 2, MPI_INT, 0, 42,
                     MPI_COMM_WORLD); 
        }


        if (world_rank == 0) {
            // collect all the new paths from the workers
            new_paths[0] = new_path[0];
            new_paths[1] = new_path[1];

            for (int j = processes_per_wire; j < world_size; 
                        j += processes_per_wire) {
                MPI_Recv(&new_path, 2, MPI_INT, j, 42, MPI_COMM_WORLD, &status);
                int path_j = j / processes_per_wire;
                new_paths[2 * path_j] = new_path[0];
                new_paths[2 * path_j + 1] = new_path[1];
            }
        }

        MPI_Bcast(&new_paths, new_paths_size, MPI_INT, 0, MPI_COMM_WORLD);

        for (int j = 0; j < new_paths_size; j += 2) {
            wire_index = i + j / 2;
            if (new_paths[j] != -1 && new_paths[j + 1] != -1) {
                change_wire_route(wires[wire_index], -1);
                set_wire(&wires[wire_index], new_paths[j], new_paths[j+1]);
                change_wire_route(wires[wire_index], 1);
            }
        }
    }
}

// Perform computation, including reading/writing output files
void compute(int procID, int nproc, char* inputFilename, double prob, 
             int numIterations) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    readInput(inputFilename);

    srand(0);
    sort_wires();

    // have all wires initially have the wire path that is just a
    // horizontal line and then a vertical line
    for (int i = 0; i < g_num_wires; i++) {
        set_wire(&wires[i], 0, wires[i].y1);
        change_wire_route(wires[i], 1);
    }

    for (int i = 0; i < numIterations; i++) {
        wire_routing(prob); 
    }

    if (world_rank == 0) {
        writeCost(inputFilename, nproc);
        writeOutput(inputFilename, nproc);

        cost_t max_cost = costs_max_overlap();
        printf("Max cost: %d\n", max_cost);
    
        cost_t agg_cost = costs_agg_overlap();
        printf("Agg cost: %d\n", agg_cost);
    }
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
