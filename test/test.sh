#!/usr/bin/env bash

FILE="$1"

if [[ -f "${FILE}" ]]; then
    if [[ -s "${FILE}" ]]; then
        echo "TEST OK"
    else
        echo "TEST FAILED: file ${FILE} is empty"
        exit 1
    fi
else
    echo "TEST FAILED: file ${FILE} is missing"
    exit 1
fi
