#!/bin/bash

if ! command -v openmc &>/dev/null; then
  echo "Error: 'openmc' command not found. Exiting." >&2
  exit 1
fi

# initialsteptime=100
timeend=$(echo "663 * 110 * 3600 * 24 / 4" | bc)
timestep=$(echo "110 * 3600 * 24" | bc)
startstep=229
shufflemult=4
shuffledays=$(echo "$timestep / $shufflemult" | bc)
steps=$(echo "$timeend / $shuffledays" | bc)

mkdir -p depletion

# Outer loop to handle restarts
while true; do
    # Inner loop to process steps
    for (( step=startstep; step<=steps; step++ )); do
        echo "##===================================== Step: ${step} =====================================##"

        # Run the Python script and capture its exit status
        python deplete.py $step $timestep |& tee "./depletion/depletion_${step}.txt"
        py_exit=${PIPESTATUS[0]}  # Get the exit code of the Python script
        echo "ERROR" $py_exit
        
        # Handle Python exit codes
        if [ $py_exit -eq 134 ]; then
            # Restart from previous step (step - 1)
            startstep=$((step - 1))
            # Ensure startstep doesn't go below 0
            [ $startstep -lt 0 ] && startstep=0
            rm -r ./depletion/step_${startstep}
            echo "Restarting from step $startstep..."
            break   # Exit for loops to restart
        elif [ $py_exit -eq 137 ]; then
            # Critical error: Exit immediately without cleanup
            echo "Keyboard Exit (exit code 137). Terminating."
            exit 1
        fi
        
        # Remove statepoint files (existing cleanup)
        rm ./depletion/step_${step}/statepoint.*
        
        # Check for errors in the Python script or cleanup
        if [ $? -ne 0 ]; then
            if grep -q "Microsoft" /proc/version; then
                echo "Critical error detected. Shutting down..."
                # shutdown.exe /s  # Uncomment for Windows shutdown
                exit 1
            else
                echo "Critical error detected. Terminating."
                exit 1
            fi
        fi
    done

    # If all steps completed without restarting, exit the outer loop
    if (( step >= steps )); then
        break
    fi
done

if grep -q "Microsoft" /proc/version; then
    echo "Done running in WSL, Shutting down the system..."
    # Call shutdown.exe for Windows (shut down immediately)
    # shutdown.exe /s
fi
