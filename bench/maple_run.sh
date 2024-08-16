#!/bin/bash

# Run an example using msolve and get the time & memory usage

# check that at least one arguments is supplied
if [[ $# -eq 0 ]]; then
    echo "No arguments supplied"
    exit -1
fi

# set binary path correctly according to env. variable
if [[ -z "${MAPLE_BINARY}" ]]; then
  BINARY="maple"
else
  BINARY="${MAPLE_BINARY}"
fi

# check that binary exists
if [[ ! $(command -v "$BINARY") ]]; then
    echo "Binary $BINARY not found"
    exit -1
fi

# get variable names without spaces and possible last comma
VARS=$(cat $1 | head -n 1 | sed 's/,$//' | tr -d '[:space:]')
# add quotes to variables names
VARS2=$(echo "$VARS" | sed -e 's/,/","/g' -e 's/^/"/g' -e 's/$/"/g')
# get field characteristic
CHAR=$(cat $1 | head -n 2 | tail -n 1 | sed 's/,$//')
# print ideal generators, remove possible last comma
GENS=$(tail -n +3 $1 | sed '$ s/,$//')

# import necessary Maple libraries
IMPORT="with(Groebner):
infolevel[GroebnerBasis] := 5:"
# print polynomial ring
# form polynomial basis
IDEAL="G := [$GENS]:";
# call groebner basis routine
BASIS="Basis(G, tdeg($VARS), characteristic = $CHAR, method = fgb):"

# input commands
INPUT=$(echo "$IMPORT"$'\n'"$IDEAL"$'\n'"$BASIS"$'\n')
# run the example
OUTPUT=$(echo "$INPUT" | taskset -c 0 "$BINARY" /dev/stdin 2>&1)
# store return value
RETVAL=$?

# print timing
echo "$OUTPUT" | grep -e "total time:" | sed 's/^.*:         //'

exit $RETVAL
