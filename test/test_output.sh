#!/bin/bash

# gamba binary directory & example filename path
BINARY_DIR="$1"
EXAMPLE_NAME="$2"

MON_ORDER="grevlex"
NUM_ELIM="0"
WEIGHTS=""

# set monomial order if requested
if [[ -n $3 ]]; then
    MON_ORDER=$3
    POSTFIX="-$3"
fi

if [[ $MON_ORDER == "blockelim" && -n $4 ]]; then
    NUM_ELIM=$4
    POSTFIX="-$3$4"
elif [[ $MON_ORDER == "grevlexw" && -n $4 ]]; then
    WEIGHTS="-w $4"
    POSTFIX="-$3$4"
fi

# run with default max. spairs
OUT1=$(${BINARY_DIR}/gamba -i ../examples/${EXAMPLE_NAME}.txt -d $MON_ORDER -e $NUM_ELIM $WEIGHTS -v -10 -o /dev/stderr 2>&1 > /dev/null)

# diff output with stored result
diff <(echo "$OUT1") <(tar -xOzf ./results/${EXAMPLE_NAME}${POSTFIX}.out.txt.tar.gz)

# store return value
RETVAL1=$?

# run with max. spairs equal to 0
OUT2=$(${BINARY_DIR}/gamba -i ../examples/${EXAMPLE_NAME}.txt -d $MON_ORDER -e $NUM_ELIM $WEIGHTS -v -10 -o /dev/stderr --max-spairs=0 2>&1 > /dev/null)

# diff output with stored result
diff <(echo "$OUT2") <(tar -xOzf ./results/${EXAMPLE_NAME}${POSTFIX}.out.txt.tar.gz)

# store return value
RETVAL2=$?

# both runs must be succesfull
! (($RETVAL1 || $RETVAL2))
