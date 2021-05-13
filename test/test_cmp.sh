#!/usr/bin/env bash

diff -w "$1" "$2";>/dev/null;CMP=$?
# if cmp -s "$1" "$2"; then
if [ ${CMP} -eq 0 ]; then
    echo "TEST OK: File $1 is the same as expected $2"
else
    echo "TEST FAILED: File $1 is not the same as expected  $2"
    exit 1
fi