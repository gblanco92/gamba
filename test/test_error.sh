#!/bin/bash

BINARY_DIR="$1"
EXAMPLE_NAME="$2"

${BINARY_DIR}/gamba -i ../examples/${EXAMPLE_NAME}.txt 2>&1 > /dev/null \
    | diff - ./results/${EXAMPLE_NAME}.error.txt

exit ${PIPESTATUS[1]}
