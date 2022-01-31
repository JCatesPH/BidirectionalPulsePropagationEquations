#!/bin/bash

make cleanall; make

PARAMFILE_PREFIX="paramfiles/paramfile_Drude"
DATA_PREFIX="DATA/Drude"

if [ "$?" -ne 0 ]; then
    printf "\nMake returned errors..\n";
    # your termination logic here
    exit 1;
fi

printf "\n\n\n"
printf "========================================\n"

print_paramfile() {
printf "# This is a test\n" > "${PARAMFILE_PREFIX}_${1}umL_${2}nE.in"
printf "Verbosity=4\n" >> "${PARAMFILE_PREFIX}_${1}umL_${2}nE.in"
printf "outputPath=${DATA_PREFIX}_${1}umL_${2}nE/\n" >> "${PARAMFILE_PREFIX}_${1}umL_${2}nE.in"
printf "numThreads=36\n" >> "${PARAMFILE_PREFIX}_${1}umL_${2}nE.in"

cat << EOF >> "${PARAMFILE_PREFIX}_${1}umL_${2}nE.in"
###########################
# Pulse parameters
meanPumpIntensity         = 10.0e9
twoColorRelativeIntensity = 0.0
twoColorPhase             = 0.0
pulseDuration             = 50.0e-15
fundamentalWavelength     = 10.0e-6

###########################
# Time domain parameters
numTimePoints             = 8192
timeDomainSize            = 6.6e-12
omegLowerCutoff           = 5
omegUpperCutoff           = 1500
tukeyWindowAlpha          = 0.0

###########################
# Spatial domain parameters
sampleLayerThickness      = ${1}.0e-6
LHSbufferThickness        = 10.0e-6
RHSbufferThickness        = 10.0e-6
initialZStep              = 5e-9
initialEDensity           = 1e${2}
EOF
}
printf "DATA_PREFIX = $DATA_PREFIX\n"
printf "PARAMFILE_PREFIX = $PARAMFILE_PREFIX\n"

for i in 50 #30 40 #10 20
do
    for j in 23 24 25 26
    do
        echo "Using ${PARAMFILE_PREFIX}_${i}umL_${j}nE.in as input."
        # Make directories for output if they do not exist
        [ ! -d "${DATA_PREFIX}_${i}umL_${j}nE" ] && mkdir "${DATA_PREFIX}_${i}umL_${j}nE"
        [ ! -d "${DATA_PREFIX}_${i}umL_${j}nE/figs" ] && mkdir "${DATA_PREFIX}_${i}umL_${j}nE/figs"
        [ ! -d "${PARAMFILE_PREFIX}_${i}umL_${j}nE.in" ] && print_paramfile $i $j

        #wait 9410
        #wait 9801
        # Run the simulation with given input file.
        #wait ppid; 
	./test.out "${PARAMFILE_PREFIX}_${i}umL_${j}nE.in" &> "${DATA_PREFIX}_${i}umL_${j}nE/stdout.txt"
        #ppid=$!
        
        #disown
        echo "Submitted job for L=$i um, nE=1e$j m^-3."

        #python ./plots_Berge.py "DATA$i"
        #cp "DATA$i/figs/Ew_both_amp.png" "AndrewFigs/ampSpectrum$i.png"
        #python ./plots.py "DATA$i"
        #cp "DATA$i/figs/Ew_both.png" "AndrewFigs/ampSpectrum$i.png"
        printf "\n========================================\n"
        
    done
done
