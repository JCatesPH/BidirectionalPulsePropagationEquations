#!/bin/bash

make cleanall; make

PARAMFILE_PREFIX="paramfiles/paramfile_Drude"
DATA_PREFIX="DATA/GridSearch_032622/Drude"

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
printf "numThreads=6\n" >> "${PARAMFILE_PREFIX}_${1}umL_${2}nE.in"

cat << EOF >> "${PARAMFILE_PREFIX}_${1}umL_${2}nE.in"
###########################
# Pulse parameters
meanPumpIntensity         = 10.0e9
pulseDuration             = 40.0e-15
fundamentalWavelength     = 1.9e-6
twoColorRelativeIntensity = 0.0
twoColorPhase             = 0.0
pulseWaist                = 0.0

###########################
# Time domain parameters
numTimePoints             = 16384
timeDomainSize            = 2.4e-12
omegLowerCutoff           = 107
omegUpperCutoff           = 1514
tukeyWindowAlpha          = 0.5

###########################
# Spatial domain parameters
sampleLayerThickness      = ${1}.0e-6
LHSbufferThickness        = 10.0e-6
RHSbufferThickness        = 10.0e-6
initialZStep              = 20e-9
initialEDensity           = 1e${2}
numTransversePoints         = 1
transverseDomainSize        = 0.0
EOF
}
printf "DATA_PREFIX = $DATA_PREFIX\n"
printf "PARAMFILE_PREFIX = $PARAMFILE_PREFIX\n"

for i in 10 20 30 40 50
do
    for j in 23 24 25
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
        echo "Finished job for L=$i um, nE=1e$j m^-3."

        #python ./plots_Berge.py "DATA$i"
        #cp "DATA$i/figs/Ew_both_amp.png" "AndrewFigs/ampSpectrum$i.png"
        #python ./plots.py "DATA$i"
        #cp "DATA$i/figs/Ew_both.png" "AndrewFigs/ampSpectrum$i.png"
        printf "\n========================================\n"
        
    done
done
