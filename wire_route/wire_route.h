#ifndef _WIRE_ROUTE_H
#define _WIRE_ROUTE_H

#include <stdint.h>

typedef int cost_t;
typedef struct
{
    int16_t x1;
    int16_t x2;
    int16_t y1;
    int16_t y2;
    int16_t bend_x1;
    int16_t bend_x2;
    int16_t bend_y1;
    int16_t bend_y2;
} wire_t;

// Perform computation, including reading/writing output files
void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations);

// Read input file
void readInput(char* inputFilename);

// Write cost array file based on input filename
void writeCost(char* inputFilename, int nproc);

// Write wire output file based on input filename
void writeOutput(char* inputFilename, int nproc);

#endif // _WIRE_ROUTE_H
