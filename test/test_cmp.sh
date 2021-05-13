#!/usr/bin/env bash

if cmp -s "$1" "$2"; then
    echo "TEST OK: File $1 is the same as expected $2"
else
    echo "TEST FAILED: File $1 is not the same as expected  $2"
    exit 1
fi