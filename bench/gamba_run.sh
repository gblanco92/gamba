#!/bin/bash

# Run an example and get the time & memory usage

# check that at least one arguments is supplied
if [[ $# -eq 0 ]]; then
    echo "No arguments supplied"
    exit -1
fi

# set binary path correctly according to env. variable
if [[ -z "${GAMBA_BINARY}" ]]; then
  BINARY="gamba"
else
  BINARY="${GAMBA_BINARY}"
fi

# check that binary exists
if [[ ! $(command -v "$BINARY") ]]; then
    echo "Binary $BINARY not found"
    exit -1
fi

# run the example
OUTPUT=$(taskset -c 0 "$BINARY" -i "$1")
# store return value
RETVAL=$?

# output time
echo "$OUTPUT" | sed -n 's/.*overall (wall)//p' | tr -d '[:space:]'
echo ""
# output memory
echo "$OUTPUT" | grep "reduce basis" | head -n 1 | tr -d '[:space:]' | sed -e 's/.*sec\(.*\)iB.*/\1/'
echo "iB"

exit $RETVAL
