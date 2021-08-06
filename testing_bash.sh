#!/bin/bash

make clean; make

if [ "$?" -ne 0 ]; then
    printf "\nMake returned errors..\n";
    # your termination logic here
    exit 1;
fi

printf "\n\n\n"
printf "========================================\n"

./test.out | tee "./DATA/stdout.txt"


