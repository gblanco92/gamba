#!/bin/bash

# Run an example and get the time & memory usage

# check that at least one arguments is supplied
if [[ $# -eq 0 ]]; then
    echo "No arguments supplied"
    exit -1
fi

# set binary path correctly according to env. variable
if [[ -z "${AXF4_BINARY}" ]]; then
  BINARY="axf4"
else
  BINARY="${AXF4_BINARY}"
fi

# check that binary exists
if [[ ! $(command -v "$BINARY") ]]; then
    echo "Binary $BINARY not found"
    exit -1
fi

# get variable names without spaces and possible last comma
VARS=$(cat $1 | head -n 1 | sed 's/,$//' | tr -d '[:space:]')
# count number of variables
NUMVARS=$(echo $VARS | tr ',' '\n' | wc -l)
# get field characteristic
CHAR=$(cat $1 | head -n 2 | tail -n 1 | sed 's/,$//')


TMPFILE=$(mktemp)
# tmp file (remove comas & spaces)
cat "$1" | tail -n +3 | sed 's/,$//' | sed 's/ //g' > "$TMPFILE"
# run the example
taskset -c 0 "$BINARY" -p "$CHAR" -v ["$VARS"] $TMPFILE
#OUTPUT=$(taskset -c 0 "$BINARY" -p "$CHAR" -v ["$VARS"] $TMPFILE)
echo "$OUTPUT" | tail -n 2 | head -n 1 | cut -d ',' -f 2
# store return value
RETVAL=$?

# remove tmp file
rm -f "$TMPFILE"
rm -f "$TMPFILE.out"

exit $RETVAL
