#!/bin/bash

# Folder containing the .inp files
WORKDIR="/home/duminy/Documents/01_Converter_Simulation/jobs"

# Abaqus scratch folder
SCRATCHDIR="/home/duminy/Documents/01_Converter_Simulation/scratch"

cd "$WORKDIR" || exit 1

mkdir -p "$SCRATCHDIR"

echo "Running all Abaqus input files in:"
echo "$WORKDIR"
echo

for inpfile in *.inp
do
    # If no .inp file exists, skip
    [ -e "$inpfile" ] || continue

    jobname="${inpfile%.inp}"

    echo "========================================"
    echo "Running job: $jobname"
    echo "Input file: $inpfile"
    echo "========================================"

    abaqus job="$jobname" input="$inpfile" scratch="$SCRATCHDIR" cpus=4 interactive

    echo
    echo "Finished job: $jobname"
    echo
done

echo "All jobs finished."