#!/bin/bash

# This script takes as input a GamBa's example file, runs the example using
# magma, converts the output back to GamBa's format and compress it to a file

# by default the ordering is grevlex
if [[ $# -eq 1 ]]; then
    ORDER="\"grevlex\""
    POSTFIX=""
fi

# set lexicographical order if requested
if [[ $# -eq 2 && ($2 -eq "lex" || $2 -eq "lexic") ]]; then
    ORDER="\"lex\""
    POSTFIX="-lex"
fi

# set block elim order if requested
if [[ $# -eq 2 ]]; then
    ORDER="\"elim\", $2"
    POSTFIX="-elim$2"
fi

# get filename without extension
FILENAME=$(basename -- "$1")
# get file extension
EXTENSION="${FILENAME##*.}"
# remove extension from filename
FILENAME="${FILENAME%.*}"

# create temporary file
TMPFILE="$FILENAME$POSTFIX.out.txt"

# get variable names without spaces and possible last comma
VARS=$(cat $1 | head -n 1 | sed 's/,$//' | tr -d '[:space:]')
# count number of variables
NUMVARS=$(echo $VARS | tr ',' '\n' | wc -l)
# get field characteristic
CHAR=$(cat $1 | head -n 2 | tail -n 1 | sed 's/,$//')

# print polynomial ring
RING="P<$VARS> := PolynomialRing(FiniteField($CHAR), $NUMVARS, ${ORDER});"
# print ideal generators, remove possible last comma
GENS=$(tail -n +3 $1 | sed '$ s/,$//')
# compute Groebner basis
BASIS="I := GroebnerBasis(I); print I;"

# input commands
INPUT=$(echo "$RING"$'\n'"I:=[$GENS];"$'\n'"$BASIS")

# print variables
echo "$VARS," > "$TMPFILE"
# print field characteristic
echo "$CHAR" >> "$TMPFILE"

# call magma
echo "$INPUT" | magma /dev/stdin |
# remove extra metada before and after generators
sed -n -e '/^\[$/,/^\]$/{ /^\[$/d; /^\]$/d; p; }' |
# remove all spaces
sed 's/[[:space:]]*//g' |
# append comma to last generator
sed '$ s/$/,/' |
# remove all new lines
sed -z 's/\n//g' |
# add new lines only after commas
sed 's/,/,\n/g' |
# reverse lines
tac |
# remove comma from last generator
sed '$ s/.$//' >> "$TMPFILE"

# save and compress the results
tar -czvf "${FILENAME}${POSTFIX}.out.txt.tar.gz" "$TMPFILE"
# remove tmpfile
rm -f "$TMPFILE"
