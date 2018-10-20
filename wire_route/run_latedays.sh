#!/bin/bash

# Max allowed nodes.
MAX_ALLOWED_NODES=4

# Each node has 24 cores.
TOTAL_PROCESSORS_PER_NODE=24

# 5 minute time limit.
WALLTIME=5

# Ensure 4 arguments for nodes, processors per node, input file, and stdout/err file.
if [ $# -ne 4 ]; then
  echo "Usage: $(basename $0) nodes processors_per_node input_file stdout_file"
  exit $E_BADARGS
fi

# Get command line arguments.
NODES=$1
PROCESSORS_PER_NODE=$2
INPUT_FILE=$3
STDOUT_FILE=$4

# Validate arguments.
if [ $NODES -le 0 ] || [ $NODES -gt $MAX_ALLOWED_NODES ]; then
    echo "ERROR: Only $MAX_ALLOWED_NODES nodes allowed."
    exit $E_BADARGS
fi
if [ $PROCESSORS_PER_NODE -le 0 ] || [ $PROCESSORS_PER_NODE -gt $TOTAL_PROCESSORS_PER_NODE ]; then
    echo "ERROR: Each node only has $TOTAL_PROCESSORS_PER_NODE cores."
    exit $E_BADARGS
fi
if [ ! -f $INPUT_FILE ]; then
    echo "ERROR: Input file $INPUT_FILE does not exist."
    exit $E_BADARGS
fi
if [ ! -f "./wire_route" ]; then
    echo "ERROR: ./wire_route program does not exist."
    exit $E_BADARGS
fi

# Submit the job.  No need to modify this.
qsub -l walltime=0:$WALLTIME:00,nodes=$NODES:ppn=$TOTAL_PROCESSORS_PER_NODE -F "$NODES $PROCESSORS_PER_NODE $INPUT_FILE" -j oe -o $STDOUT_FILE latedays.qsub
