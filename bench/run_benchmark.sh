#!/bin/bash

# Run all examples in examples.txt with the choosen CAS system

# Supported CAS systems:
# 1. gamba
# 2. magma
# 3. msolve
# 4. julia
# 5. maple
# 6. singular

# Disable Hyperthreading
echo forceoff | sudo tee /sys/devices/system/cpu/smt/control > /dev/null 2>&1
# Disable Turbo Boost (Intel)
echo "1" | sudo tee /sys/devices/system/cpu/intel_pstate/no_turbo > /dev/null 2>&1
# Disable Turbo Boost
echo "0" | sudo tee /sys/devices/system/cpu/cpufreq/boost > /dev/null 2>&1
# Set CPU scaling governor to performance
sudo cpupower frequency-set -g performance > /dev/null 2>&1

# Path to the current script location relative to pwd
SCRIPT_PATH=$(dirname "${BASH_SOURCE[0]}")
# Path to the examples files relative to the script path
EXAMPLES_FILE="$SCRIPT_PATH""/examples.txt"
# Examples folder relarive to the script path
EXAMPLES_DIR="$SCRIPT_PATH""/../examples/"

# Check that at least one arguments is supplied
if [[ $# -eq 0 ]]; then
    echo "No arguments supplied"
    exit -1
fi

# Path to the script calling the select CAS system
RUN_SCRIPT="$SCRIPT_PATH""/""$1""_run.sh"

# Check is CAS system is supported in the benchmark
if [[ ! -f "$RUN_SCRIPT" ]]; then
    echo "$1 not supported"
    exit -1
fi

# Loop through the lines in examples.txt
while read line; do
    # Ignore empty or commented lines
    [[ "$line" =~ ^#.*$ || ! -n "$line" ]] && continue

    # Print example name
    echo "$line"
    # Run the example through script
    timeout 24h "$RUN_SCRIPT" "$EXAMPLES_DIR""$line"".txt"
    # Return value
    RETVAL=$?

    # Check if it ended with a timeout
    if [[ $RETVAL == 124 ]]; then
        echo "Timeout"
    # Killed by memory exahustion
    elif [[ $RETVAL != 0 ]]; then
        echo "Memory exhausted"
    # Otherwise, run it x2 times for a best of 3
    else
        "$RUN_SCRIPT" "$EXAMPLES_DIR""$line"".txt"
        "$RUN_SCRIPT" "$EXAMPLES_DIR""$line"".txt"
    fi

    # Delimiter
    echo "======================="
done <$EXAMPLES_FILE
