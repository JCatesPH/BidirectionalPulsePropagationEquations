#!/bin/bash

make clean; make

if [ "$?" -ne 0 ]; then
    printf "\nMake returned errors..\n";
    # your termination logic here
    exit 1;
fi

printf "\n\n\n"
printf "========================================\n"


for i in 1 2 3 4 5 6 7 8 9 10
do
    echo "Using paramfile$i.in as input."
    # Make directories for output if they do not exist
    [ ! -d "DATA$i" ] && mkdir "DATA$i"
    [ ! -d "DATA$i/figs" ] && mkdir "DATA$i/figs"

    # Run the simulation with given input file.
    ./test.out "paramfiles/paramfile$i.in" > "./DATA$i/stdout.txt"
    echo "Simulation complete."

    python ./plots_Andrew.py "DATA$i"
    cp "DATA$i/figs/Ew_both_amp.png" "AndrewFigs/ampSpectrum$i.png"
    #python ./plots.py "DATA$i"
    #cp "DATA$i/figs/Ew_both.png" "AndrewFigs/ampSpectrum$i.png"
done

