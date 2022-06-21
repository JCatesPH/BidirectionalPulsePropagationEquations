#!/bin/bash
DATADIR='/home/jmcates/BPPE-Plasma/DATA/silica_fdtdComp_LJet_061322'

for j in {0..9}
do
    echo "$j"
    cp "${DATADIR}/Spectrum_iteration_${j}_Reflected.dat" "${DATADIR}/Spectrum_iteration_0${j}_Reflected.dat"
    cp "${DATADIR}/Spectrum_iteration_${j}_Transmitted.dat" "${DATADIR}/Spectrum_iteration_0${j}_Transmitted.dat"


done
