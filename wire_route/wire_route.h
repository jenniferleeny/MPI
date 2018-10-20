#ifndef _WIRE_ROUTE_H
#define _WIRE_ROUTE_H

// Perform computation, including reading/writing output files
void compute(int procID, int nproc, char* inputFilename, double prob, int numIterations);

// Read input file
void readInput(char* inputFilename);

// Write cost array file based on input filename
void writeCost(char* inputFilename, int nproc);

// Write wire output file based on input filename
void writeOutput(char* inputFilename, int nproc);

#endif // _WIRE_ROUTE_H
