#!/bin/bash

macroName="createCalibFile.C"
nCores=10
logdir="calibLogs"


mkdir -p ${logdir}

echo "Starting batch MoNA calibration, ${nCores} cores"
echo "Logs output to ${logdir}/"
echo ""

function wait_for_slot() {
    while [ "$(jobs -r | wc -l)" -ge "${nCores}" ]; do
        sleep 1
    done
}

for file in "$@"; do

    wait_for_slot

    base=$(basename "$file" .root) 
    logfile="${logdir}/${base}.log"

    echo "Launching calibration for ${file}"

    root -l -b -q "${macroName}(${file})" > "${logfile}" 2>&1 &

done

wait

echo ""
echo "Calibrations complete"
