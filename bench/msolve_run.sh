#!/bin/bash

# Run an example using msolve and get the time & memory usage

# check that at least one arguments is supplied
if [[ $# -eq 0 ]]; then
    echo "No arguments supplied"
    exit -1
fi

# set binary path correctly according to env. variable
if [[ -z "${MSOLVE_BINARY}" ]]; then
  BINARY="msolve"
else
  BINARY="${MSOLVE_BINARY}"
fi

# check that binary exists
if [[ ! $(command -v "$BINARY") ]]; then
    echo "Binary $BINARY not found"
    exit -1
fi

# run the example
OUTPUT=$(taskset -c 0 "$BINARY" -g 2 -l 44 -v 1 -f "$1" -o /dev/null 2>&1)
# store return value
RETVAL=$?

# output time
echo "$OUTPUT" | sed -n 's/.*overall(elapsed)//p' | tr -d '[:space:]'
echo ""

exit $RETVAL
