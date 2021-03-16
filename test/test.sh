#!/usr/bin/env bash

HIC_FILE="$1"

if [[ -f "${HIC_FILE}" ]]; then
    if [[ -s "${HIC_FILE}" ]]; then
        echo "TEST OK"
    else
        echo "TEST FAILED: file ${HIC_FILE} is empty"
        exit 1
    fi
else
    echo "TEST FAILED: file ${HIC_FILE} is missing"
    exit 1
fi
