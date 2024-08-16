#!/bin/bash

# Run an example using msolve and get the time & memory usage

# check that at least one arguments is supplied
if [[ $# -eq 0 ]]; then
    echo "No arguments supplied"
    exit -1
fi

# set binary path correctly according to env. variable
if [[ -z "${JULIA_BINARY}" ]]; then
  BINARY="julia"
else
  BINARY="${JULIA_BINARY}"
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

# import necessary Julia libraries
IMPORT="using BenchmarkTools; using AbstractAlgebra; using Groebner;"
# print polynomial ring
RING="R, ($VARS) = polynomial_ring(GF($CHAR), [$VARS2]);";
# form polynomial basis
IDEAL="system = [$GENS];";
# call groebner basis routine
BASIS="timing = @timed result = groebner(system, ordering=DegRevLex(), linalg=:randomized, threaded=:no);"

# input commands (call groebner twice to compile, how stupid...)
INPUT=$(echo "$IMPORT"$'\n'"$RING"$'\n'"$IDEAL"$'\n'"$BASIS"$'\n'"$BASIS"$'\n'"println(timing.time, \" sec\");"$'\n'"println(timing.bytes, \" bytes\");")

# run the example
OUTPUT=$(echo "$INPUT" | taskset -c 0 "$BINARY" -L /dev/stdin 2>&1)
# store return value
RETVAL=$?

# print timing & memory usage
echo "$OUTPUT"

exit $RETVAL
