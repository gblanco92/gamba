#!/bin/bash

# Run an example in GamBa's format using Magma and get the time & memory usage

# check that at least one arguments is supplied
if [[ $# -eq 0 ]]; then
    echo "No arguments supplied"
    exit -1
fi

# set binary path correctly according to env. variable
if [[ -z "${MAGMA_BINARY}" ]]; then
  BINARY="magma"
else
  BINARY="${MAGMA_BINARY}"
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

# print polynomial ring
RING="P<$VARS> := PolynomialRing(FiniteField($CHAR), $NUMVARS, \"grevlex\");"
# print ideal generators, remove possible last comma
GENS=$(tail -n +3 $1 | sed '$ s/,$//')
# set verbosity level
VERBOSE="SetVerbose(\"Faugere\", 1);"
# compute Groebner basis
BASIS="Groebner(I : Al := \"Direct\", Faugere := true);"

# input commands
INPUT=$(echo "$RING"$'\n'"I:=Ideal([$GENS]);"$'\n'"$VERBOSE"$'\n'"$BASIS")

# call magma
OUTPUT=$(echo "$INPUT" | taskset -c 0 $BINARY /dev/stdin | grep -e "real time" -e "memory usage")
# store return value
RETVAL=$?

# output time
echo "$OUTPUT" | sed -n 's/.*real time: //p' /dev/stdin
# output memory
echo "$OUTPUT" | sed -n 's/.*memory usage: //p' /dev/stdin

exit $RETVAL
