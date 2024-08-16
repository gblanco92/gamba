#!/bin/bash

# gamba binary directory & example filename path
BINARY_DIR="$1"
EXAMPLE_NAME="$2"
NUM_ELIM="0"

# set lexicographical order if requested
if [[ $# -eq 3 ]]; then
    NUM_ELIM=$3
    POSTFIX="-elim$3"
fi

# run with default max. spairs
OUT1=$(${BINARY_DIR}/gamba -i ../examples/${EXAMPLE_NAME}.txt -e $NUM_ELIM -o /dev/stderr 2>&1 > /dev/null)

# diff output with stored result
diff <(echo "$OUT1") <(tar -xOzf ./results/${EXAMPLE_NAME}${POSTFIX}.out.txt.tar.gz)

# store return value
RETVAL1=$?

# run with max. spairs equal to 0
OUT2=$(${BINARY_DIR}/gamba -i ../examples/${EXAMPLE_NAME}.txt -e $NUM_ELIM -o /dev/stderr --max-spairs=0 2>&1 > /dev/null)

# diff output with stored result
diff <(echo "$OUT2") <(tar -xOzf ./results/${EXAMPLE_NAME}${POSTFIX}.out.txt.tar.gz)

# store return value
RETVAL2=$?

# both runs must be succesfull
! (($RETVAL1 || $RETVAL2))
